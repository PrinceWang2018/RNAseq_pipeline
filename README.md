# RNAseq_pipeline

A lightweight, single-script pipeline for **paired-end bulk RNA-seq** on HPC clusters.
It turns raw paired FASTQs into a gene count matrix and a QC report:

```
FASTQ pairs ──► Trim Galore (+FastQC) ──► STAR ──► samtools index ──► featureCounts ──► MultiQC
                 01_trimmed/               02_alignment/                03_counts/          04_multiqc/
```

- **Any common FASTQ layout.** FASTQs can sit in one folder or in one folder per sample. Multiple lanes are merged automatically.
- **No genome version to choose.** You give the path to any STAR index. The GTF is read from the index unless you pass one.
- **One conda environment** holds every tool.
- **Works with any scheduler.** It is plain Bash. Threads are read from SLURM, PBS or LSF, and a SLURM template is included.
- **Resumable.** Finished steps are recorded, so resubmitting the same command carries on where it stopped.
- **Useful extras.** It writes a clean `gene_counts.tsv` matrix, a strandedness check, software versions and per-step logs.

---

## 1. Installation

Requirements: Linux, Bash ≥ 4.2, and [conda](https://github.com/conda-forge/miniforge) (Miniforge/mamba recommended).

```bash
git clone https://github.com/<your-account>/RNAseq_pipeline.git
cd RNAseq_pipeline

# Option A: helper script (uses mamba if available)
./install.sh                 # creates env "rnaseq"
# ./install.sh myenv         # custom name
# ./install.sh -p /shared/envs/rnaseq   # install to a prefix (e.g. project space)

# Option B: manually
mamba env create -f environment.yml    # or: conda env create -f environment.yml

# Verify
conda activate rnaseq
./rnaseq_pipeline.sh --check
```

Software in the environment (`environment.yml`):

| Tool | Version | Purpose |
|---|---|---|
| STAR | 2.7.11b (pinned) | alignment |
| Trim Galore | ≥ 0.6.7 (0.6.x and 2.x both supported) | adapter & quality trimming |
| FastQC | ≥ 0.12 | read QC |
| samtools | ≥ 1.17 | BAM indexing |
| Subread (featureCounts) | ≥ 2.0.3 | gene counting |
| MultiQC | ≥ 1.19 | summary report |
| pigz | — | parallel (de)compression |

> **Compute nodes without internet:** create the environment on a login node. It lives in your conda directory and can be activated from any node.

---

## 2. Prepare a STAR index (once per genome)

If your group already has a STAR index built with **STAR ≥ 2.7.4a**, you can use it as it is. Otherwise, build one:

```bash
conda activate rnaseq
STAR --runMode genomeGenerate \
     --runThreadN 16 \
     --genomeDir /path/to/STAR_index_GRCh38 \
     --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
     --sjdbGTFfile gencode.v44.primary_assembly.annotation.gtf \
     --sjdbOverhang 149        # read length - 1
```

- A human or mouse index needs about **32 GB of RAM** to build and to align.
- If the index was built with `--sjdbGTFfile`, you can leave out `--gtf` and the pipeline uses the same GTF. Pass `--gtf` if that file has moved, or if you want a different annotation.
- The FASTA and the GTF must use the same chromosome names (`chr1` vs `1`).

---

## 3. Input

### A. All FASTQs in one folder

```
fastq/
├── Ctrl1_R1_001.fastq.gz   Ctrl1_R2_001.fastq.gz
├── Ctrl2_1.fq.gz           Ctrl2_2.fq.gz
├── Treat1_L001_R1.fastq.gz Treat1_L001_R2.fastq.gz   ┐ merged into
└── Treat1_L002_R1.fastq.gz Treat1_L002_R2.fastq.gz   ┘ sample "Treat1"
```
The sample name is the file name with the read tag and extension removed. A trailing lane tag (`_L001`, `_L1`) is also removed, so a sample's lanes are merged.

### B. One folder per sample (e.g. `01.RawData` from a sequencing provider)

```
01.RawData/
├── Ctrl1/  Ctrl1_1.fq.gz  Ctrl1_2.fq.gz
└── KO1/    KO1_L1_1.fq.gz KO1_L1_2.fq.gz KO1_L2_1.fq.gz KO1_L2_2.fq.gz
```
The sample name is the **folder name**. All read pairs inside a folder are merged into that sample.

Both layouts are passed with `--input` and can even be mixed.

**Recognised read-pair names** (extensions `.fastq.gz`, `.fq.gz`, `.fastq`, `.fq`):

| R1 | R2 |
|---|---|
| `*_R1_001.fastq.gz` | `*_R2_001.fastq.gz` |
| `*_R1.fq.gz` | `*_R2.fq.gz` |
| `*_1.fq.gz` | `*_2.fq.gz` |
| `*.R1.fq.gz` | `*.R2.fq.gz` |
| `*.1.fq.gz` | `*.2.fq.gz` |

### C. Sample sheet (any other naming, or to rename samples)

CSV or TSV with a header. Repeat a sample on several rows to merge lanes. Relative paths are resolved from the sheet's folder.

```csv
sample,fastq_1,fastq_2
Ctrl_rep1,raw/A01_S1_R1.fq.gz,raw/A01_S1_R2.fq.gz
KO_rep1,raw/B07_lane1_R1.fq.gz,raw/B07_lane1_R2.fq.gz
KO_rep1,raw/B07_lane2_R1.fq.gz,raw/B07_lane2_R2.fq.gz
```

Tip: each run writes the samples it found to `pipeline_info/samples.tsv`. You can edit that file and pass it back with `--samplesheet`.

**Always check sample detection first with a dry run:**

```bash
./rnaseq_pipeline.sh -i fastq/ -o results/ -x /path/to/STAR_index --dry-run
```

---

## 4. Usage

```bash
conda activate rnaseq
./rnaseq_pipeline.sh \
    --input      /path/to/fastq \
    --outdir     /path/to/results \
    --star-index /path/to/STAR_index \
    --threads    16
```

| Option | Default | Description |
|---|---|---|
| `-i, --input DIR` | | FASTQ folder (layout A or B) |
| `--samplesheet FILE` | | `sample,fastq_1,fastq_2` sheet (instead of `--input`) |
| `-o, --outdir DIR` | *required* | output directory |
| `-x, --star-index DIR` | *required* | STAR index directory |
| `-a, --gtf FILE` | from index | GTF annotation (plain or gzipped) |
| `-s, --strandedness N` | `0` | featureCounts `-s`: 0 unstranded, 1 stranded, 2 reverse |
| `--feature-type STR` | `exon` | featureCounts `-t` |
| `--attribute STR` | `gene_id` | featureCounts `-g` (e.g. `gene_name`) |
| `--featurecounts-extra "ARGS"` | | e.g. `"-M --fraction"` |
| `--trim-quality N` | `20` | Trim Galore `-q` |
| `--trim-length N` | `20` | Trim Galore `--length` |
| `--trim-stringency N` | `3` | Trim Galore `--stringency` |
| `--trim-error-rate F` | `0.1` | Trim Galore `-e` |
| `--skip-trim` | | align input FASTQs directly |
| `--skip-fastqc` | | no FastQC on trimmed reads |
| `--delete-trimmed` | | delete trimmed FASTQs once a sample is aligned |
| `--star-extra "ARGS"` | | e.g. `"--limitBAMsortRAM 30000000000"` |
| `-t, --threads N` | `$SLURM_CPUS_PER_TASK` / `nproc` | threads |
| `--skip-multiqc` | | no MultiQC report |
| `-n, --dry-run` | | detect samples, validate inputs, check software, exit |
| `--force` | | ignore finished-step markers and rerun all |
| `--check` | | check software and exit |

### SLURM

Edit the paths at the top of `slurm_submit.sh`, then:

```bash
sbatch slurm_submit.sh
```

### PBS / Torque

```bash
#!/bin/bash
#PBS -N rnaseq
#PBS -l nodes=1:ppn=16,mem=64gb,walltime=72:00:00
cd $PBS_O_WORKDIR
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate rnaseq
bash /path/to/rnaseq_pipeline.sh -i fastq/ -o results/ -x /path/to/STAR_index -t 16
```

**Resources:** use at least 32–40 GB of memory for human or mouse, plus about 4 GB per sorting thread. Samples are processed one after another, so the run time grows linearly with the number of samples.

---

## 5. Output

```
results/
├── 01_trimmed/                  trimmed FASTQs, Trim Galore reports, FastQC
├── 02_alignment/
│   ├── <sample>.bam / .bam.bai  coordinate-sorted, indexed BAM (read group = sample)
│   └── <sample>/                STAR Log.final.out, SJ.out.tab, ReadsPerGene.out.tab
├── 03_counts/
│   ├── gene_counts.tsv          ★ count matrix: gene_id + one column per sample
│   ├── featureCounts.txt        raw featureCounts output (with gene lengths)
│   ├── featureCounts.txt.summary
│   └── strandedness_check.tsv   library strandedness inferred from STAR
├── 04_multiqc/multiqc_report.html
├── logs/                        per-sample, per-tool logs
├── pipeline_info/               params.txt, software_versions.txt, samples.tsv, pipeline.done
└── pipeline.log
```

`gene_counts.tsv` can be loaded straight into DESeq2 or edgeR:

```r
counts <- read.delim("results/03_counts/gene_counts.tsv", row.names = 1, check.names = FALSE)
```

---

## 6. Notes

### Strandedness
The wrong `-s` value can throw away most of your counts. When the STAR index contains annotation, `03_counts/strandedness_check.tsv` compares the counts for each setting and suggests a value. The log also warns if that suggestion differs from `--strandedness`. Most dUTP kits, such as Illumina TruSeq Stranded, are reverse stranded (`-s 2`). If the value was wrong, rerun with the correct `-s`. Only featureCounts runs again, because trimming and alignment are already marked as finished.

### Resuming and rerunning
A marker is written to `<outdir>/.done/` after each sample finishes a step. Resubmitting the same command skips finished work. The counting step reruns by itself when the sample list, GTF, or counting options change. Use `--force` to redo everything, for example after changing trimming options.

### Counting mode
The pipeline runs featureCounts in paired-end mode (`-p --countReadPairs`), which counts **fragments**, not reads.

---

## 7. Troubleshooting

| Symptom | Fix |
|---|---|
| `Genome version ... is INCOMPATIBLE` | The index was built with STAR < 2.7.4a. Rebuild it with this environment. |
| `not enough memory for BAM sorting` | `--star-extra "--limitBAMsortRAM 30000000000"`, or request more memory |
| STAR is killed or `std::bad_alloc` | Request more memory: about 32 GB for human or mouse |
| Nearly all reads are `Unassigned_NoFeatures` | Chromosome names differ between the GTF and the index, or `--strandedness` is wrong |
| A FASTQ is "ignored" | Its name does not match the recognised patterns. Rename it or use `--samplesheet`. |
| `R1 ... but its mate ... is missing` | The R2 file is missing or named differently |

---

## Citation

If you use this pipeline, please cite the underlying tools: STAR (Dobin *et al.*, 2013), Cutadapt (Martin, 2011), Trim Galore, FastQC, SAMtools (Danecek *et al.*, 2021), featureCounts (Liao *et al.*, 2014) and MultiQC (Ewels *et al.*, 2016).

## License

GPL-3.0. See [LICENSE](LICENSE).
