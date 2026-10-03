#!/usr/bin/env bash
#
# rnaseq_pipeline.sh - paired-end bulk RNA-seq pipeline
#   Trim Galore (+FastQC) -> STAR -> samtools index -> featureCounts -> MultiQC
#
# Run `rnaseq_pipeline.sh --help` for usage.

set -Eeuo pipefail

VERSION="1.0.0"
SCRIPT_NAME=$(basename "$0")

# ----------------------------------------------------------------------------
# Defaults
# ----------------------------------------------------------------------------
INPUT_DIR=""
SAMPLESHEET=""
OUTDIR=""
STAR_INDEX=""
GTF=""
THREADS="${SLURM_CPUS_PER_TASK:-${NCPUS:-${PBS_NUM_PPN:-${LSB_DJOB_NUMPROC:-}}}}"
[[ -n $THREADS ]] || THREADS=$(nproc 2>/dev/null || echo 8)
STRAND=0
FEATURE_TYPE="exon"
ATTR_TYPE="gene_id"
TRIM_QUALITY=20
TRIM_LENGTH=20
TRIM_STRINGENCY=3
TRIM_ERROR=0.1
STAR_EXTRA=""
FC_EXTRA=""
SKIP_TRIM=0
SKIP_FASTQC=0
SKIP_MULTIQC=0
DELETE_TRIMMED=0
DRY_RUN=0
FORCE=0
CHECK_ONLY=0

# Recognised read-pair naming: "<R1 tag>:<R2 tag>", matched before the extension.
PAIR_TAGS=("_R1_001:_R2_001" "_R1:_R2" "_1:_2" ".R1:.R2" ".1:.2")
FASTQ_EXTS=(".fastq.gz" ".fq.gz" ".fastq" ".fq")

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
ts()   { date '+%Y-%m-%d %H:%M:%S'; }
log()  { printf '[%s] %s\n' "$(ts)" "$*"; }
warn() { printf '[%s] WARNING: %s\n' "$(ts)" "$*" >&2; }
die()  { printf '[%s] ERROR: %s\n' "$(ts)" "$*" >&2; exit 1; }

usage() {
    cat <<EOF
$SCRIPT_NAME v$VERSION - paired-end RNA-seq: Trim Galore -> STAR -> featureCounts -> MultiQC

Usage:
  $SCRIPT_NAME -i <fastq_dir> -o <outdir> -x <star_index> [-a <gtf>] [options]
  $SCRIPT_NAME --samplesheet <samples.csv> -o <outdir> -x <star_index> [options]
  $SCRIPT_NAME --check

Input (choose one):
  -i, --input DIR          Folder with paired FASTQs. Both layouts are accepted:
                             DIR/S1_R1.fq.gz, DIR/S1_R2.fq.gz ...       (flat)
                             DIR/S1/S1_1.fq.gz, DIR/S1/S1_2.fq.gz ...   (one folder per sample)
      --samplesheet FILE   CSV/TSV with columns: sample,fastq_1,fastq_2
                           (relative paths are relative to the sheet's folder)

Required:
  -o, --outdir DIR         Output directory
  -x, --star-index DIR     STAR genome index directory

Annotation / counting:
  -a, --gtf FILE           GTF annotation (plain or .gz). If omitted, the GTF recorded
                           in the STAR index (genomeParameters.txt) is used.
  -s, --strandedness N     featureCounts -s: 0 unstranded, 1 stranded, 2 reversely
                           stranded [default: $STRAND]
      --feature-type STR   featureCounts -t [default: $FEATURE_TYPE]
      --attribute STR      featureCounts -g [default: $ATTR_TYPE]
      --featurecounts-extra "ARGS"  Extra arguments passed to featureCounts

Trimming:
      --trim-quality N     Trim Galore -q [default: $TRIM_QUALITY]
      --trim-length N      Trim Galore --length [default: $TRIM_LENGTH]
      --trim-stringency N  Trim Galore --stringency [default: $TRIM_STRINGENCY]
      --trim-error-rate F  Trim Galore -e [default: $TRIM_ERROR]
      --skip-trim          Align the input FASTQs directly
      --skip-fastqc        Do not run FastQC on trimmed reads
      --delete-trimmed     Delete trimmed FASTQs after a sample is aligned

Alignment:
      --star-extra "ARGS"  Extra arguments passed to STAR

General:
  -t, --threads N          Threads [default: \$SLURM_CPUS_PER_TASK or nproc; now $THREADS]
      --skip-multiqc       Do not run MultiQC
  -n, --dry-run            Detect samples, validate inputs and exit
      --force              Ignore finished-step markers and rerun everything
      --check              Check that all required software is available and exit
  -h, --help               Show this help
  -v, --version            Show version

Finished steps are recorded in <outdir>/.done/, so re-running the same command
resumes after the last completed step.
EOF
}

# Expand a possibly-empty array safely under `set -u` on bash 4.2.
# Usage: "${arr[@]+"${arr[@]}"}"

abs_dir() { (cd "$1" 2>/dev/null && pwd) || die "Directory not found: $1"; }

abs_path() {
    local p=$1
    [[ $p == /* ]] || p="$PWD/$p"
    # Normalise trailing slashes without resolving symlinks.
    while [[ $p == */ && $p != / ]]; do p=${p%/}; done
    printf '%s\n' "$p"
}

run_logged() {
    # run_logged <logfile> <command...>
    local logfile=$1
    shift
    printf '[%s] CMD: %s\n' "$(ts)" "$*" >>"$logfile"
    if ! "$@" >>"$logfile" 2>&1; then
        printf '[%s] ERROR: command failed: %s\n' "$(ts)" "$1" >&2
        printf '---- last 30 lines of %s ----\n' "$logfile" >&2
        tail -n 30 "$logfile" >&2
        exit 1
    fi
}

done_marker() { printf '%s/.done/%s.%s\n' "$OUTDIR" "$1" "$2"; }

is_done() {
    # is_done <step> <id> [content]: true if the marker exists (and matches content)
    local marker
    marker=$(done_marker "$1" "$2")
    [[ -f $marker ]] || return 1
    [[ $# -lt 3 ]] || [[ "$(<"$marker")" == "$3" ]]
}

mark_done() {
    local marker
    marker=$(done_marker "$1" "$2")
    printf '%s\n' "${3:-}" >"$marker"
}

tool_version() {
    local out
    case $1 in
        STAR)          out=$(STAR --version 2>&1) ;;
        trim_galore)   out=$(trim_galore --version 2>&1 | grep -Eo '[0-9]+\.[0-9]+\.[0-9]+' | head -n1) ;;
        cutadapt)      out=$(cutadapt --version 2>&1) ;;
        fastqc)        out=$(fastqc --version 2>&1) ;;
        samtools)      out=$(samtools --version 2>&1 | head -n1) ;;
        featureCounts) out=$(featureCounts -v 2>&1 | grep -Eo 'v[0-9][0-9.]*' | head -n1) ;;
        multiqc)       out=$(multiqc --version 2>&1) ;;
        pigz)          out=$(pigz --version 2>&1) ;;
    esac || true
    printf '%s\n' "${out:-unknown}"
}

check_tools() {
    # cutadapt is only needed by Trim Galore < 2.0 (its conda package pulls it in).
    local required=(STAR samtools featureCounts) optional=() t missing=0
    if ((SKIP_TRIM)); then
        optional+=(trim_galore cutadapt pigz)
    else
        required+=(trim_galore)
        optional+=(cutadapt pigz)
    fi
    if ((SKIP_FASTQC || SKIP_TRIM)); then optional+=(fastqc); else required+=(fastqc); fi
    if ((SKIP_MULTIQC)); then optional+=(multiqc); else required+=(multiqc); fi

    for t in "${required[@]}"; do
        if command -v "$t" >/dev/null 2>&1; then
            printf '  %-14s %-10s %s\n' "$t" "ok" "$(tool_version "$t")"
        else
            printf '  %-14s %s\n' "$t" "MISSING (required)"
            missing=1
        fi
    done
    for t in "${optional[@]}"; do
        if command -v "$t" >/dev/null 2>&1; then
            printf '  %-14s %-10s %s\n' "$t" "ok" "$(tool_version "$t")"
        else
            printf '  %-14s %s\n' "$t" "not found (optional)"
        fi
    done
    return "$missing"
}

# ----------------------------------------------------------------------------
# Sample discovery
# ----------------------------------------------------------------------------
declare -A SAMPLE_R1=() SAMPLE_R2=()
SAMPLE_ORDER=()

add_pair() {
    local s=$1 r1=$2 r2=$3
    [[ $s =~ ^[A-Za-z0-9._-]+$ ]] ||
        die "Invalid sample name '$s' (allowed characters: letters, digits, '.', '_', '-')"
    [[ ! $r1$r2 =~ [[:space:],] ]] || die "FASTQ paths must not contain spaces or commas: $r1 $r2"
    [[ -f $r1 ]] || die "FASTQ not found: $r1"
    [[ -f $r2 ]] || die "FASTQ not found: $r2"
    [[ $r1 != "$r2" ]] || die "R1 and R2 are the same file for sample $s: $r1"
    if [[ -z ${SAMPLE_R1[$s]+x} ]]; then
        SAMPLE_ORDER+=("$s")
        SAMPLE_R1[$s]=$r1
        SAMPLE_R2[$s]=$r2
    else
        [[ ",${SAMPLE_R1[$s]}," != *",$r1,"* ]] || die "Duplicate FASTQ for sample $s: $r1"
        SAMPLE_R1[$s]+=",$r1"
        SAMPLE_R2[$s]+=",$r2"
    fi
}

# match_mate <basename> <1|2>: if the name is a read-<n> file, set M_PREFIX and M_MATE
match_mate() {
    local b=$1 which=$2 ext tag mine mate
    for ext in "${FASTQ_EXTS[@]}"; do
        [[ $b == *"$ext" ]] || continue
        for tag in "${PAIR_TAGS[@]}"; do
            if [[ $which == 1 ]]; then mine=${tag%%:*}; mate=${tag##*:}; else mine=${tag##*:}; mate=${tag%%:*}; fi
            if [[ $b == ?*"$mine$ext" ]]; then
                M_PREFIX=${b%"$mine$ext"}
                M_MATE="$M_PREFIX$mate$ext"
                return 0
            fi
        done
        return 1
    done
    return 1
}

discover_from_dir() {
    local dir=$1 f b parent s
    local -a files=()
    local -A used=()

    while IFS= read -r -d '' f; do
        [[ $f == "$OUTDIR"/* ]] && continue
        files+=("$f")
    done < <(find -L "$dir" -mindepth 1 -maxdepth 2 -type f \
        \( -name '*.fastq.gz' -o -name '*.fq.gz' -o -name '*.fastq' -o -name '*.fq' \) \
        -print0 2>/dev/null | sort -z)

    ((${#files[@]} > 0)) || die "No FASTQ files (*.fastq.gz, *.fq.gz, *.fastq, *.fq) found in $dir or its subfolders"

    for f in "${files[@]}"; do
        b=$(basename "$f")
        parent=$(dirname "$f")
        match_mate "$b" 1 || continue
        [[ -f "$parent/$M_MATE" ]] || die "Found R1 '$f' but its mate '$parent/$M_MATE' is missing"
        if [[ $parent == "$dir" ]]; then
            # Flat layout: sample = file prefix, minus a lane tag (_L001 / _L1) so lanes merge.
            s=$(sed -E 's/_L[0-9]{1,3}$//' <<<"$M_PREFIX")
        else
            # Per-sample folders: sample = folder name; all pairs inside are merged.
            s=$(basename "$parent")
        fi
        add_pair "$s" "$f" "$parent/$M_MATE"
        used[$f]=1
        used[$parent/$M_MATE]=1
    done

    for f in "${files[@]}"; do
        [[ -n ${used[$f]+x} ]] && continue
        b=$(basename "$f")
        if match_mate "$b" 2; then
            die "Found R2 '$f' but its mate '$(dirname "$f")/$M_MATE' is missing"
        fi
        warn "Ignoring FASTQ with unrecognised pair naming: $f"
    done
    ((${#SAMPLE_ORDER[@]} > 0)) || die "No read pairs recognised in $dir. Supported R1/R2 tags: ${PAIR_TAGS[*]} (or use --samplesheet)"
}

discover_from_samplesheet() {
    local sheet=$1 base line s r1 r2 extra lineno=0
    [[ -f $sheet ]] || die "Sample sheet not found: $sheet"
    base=$(abs_dir "$(dirname "$sheet")")
    while IFS= read -r line || [[ -n $line ]]; do
        lineno=$((lineno + 1))
        line=${line%$'\r'}
        [[ -z ${line//[[:space:]]/} || $line == \#* ]] && continue
        IFS=$',\t' read -r s r1 r2 extra <<<"$line"
        if ((lineno == 1)) && [[ $s == sample ]]; then continue; fi
        [[ -n $s && -n $r1 && -n $r2 ]] || die "$sheet line $lineno: expected 'sample,fastq_1,fastq_2'"
        [[ $r1 == /* ]] || r1="$base/$r1"
        [[ $r2 == /* ]] || r2="$base/$r2"
        add_pair "$s" "$r1" "$r2"
    done <"$sheet"
    ((${#SAMPLE_ORDER[@]} > 0)) || die "No samples found in $sheet"
}

print_samples() {
    local s n
    printf '  %-28s %-6s %s\n' "SAMPLE" "PAIRS" "FIRST R1"
    for s in "${SAMPLE_ORDER[@]}"; do
        n=$(tr ',' '\n' <<<"${SAMPLE_R1[$s]}" | wc -l)
        printf '  %-28s %-6s %s\n' "$s" "$n" "${SAMPLE_R1[$s]%%,*}"
    done
}

write_samplesheet() {
    local out=$1 s i
    local -a r1s r2s
    printf 'sample\tfastq_1\tfastq_2\n' >"$out"
    for s in "${SAMPLE_ORDER[@]}"; do
        IFS=',' read -r -a r1s <<<"${SAMPLE_R1[$s]}"
        IFS=',' read -r -a r2s <<<"${SAMPLE_R2[$s]}"
        for i in "${!r1s[@]}"; do
            printf '%s\t%s\t%s\n' "$s" "${r1s[$i]}" "${r2s[$i]}" >>"$out"
        done
    done
}

# ----------------------------------------------------------------------------
# Pipeline steps
# ----------------------------------------------------------------------------

# prepare_inputs <sample>: set IN_R1/IN_R2, concatenating multiple lanes if needed
prepare_inputs() {
    local s=$1 ext f
    local -a r1s r2s
    IFS=',' read -r -a r1s <<<"${SAMPLE_R1[$s]}"
    IFS=',' read -r -a r2s <<<"${SAMPLE_R2[$s]}"
    IS_MERGED=0
    if ((${#r1s[@]} == 1)); then
        IN_R1=${r1s[0]}
        IN_R2=${r2s[0]}
        return
    fi

    if [[ ${r1s[0]} == *.gz ]]; then ext=".fastq.gz"; else ext=".fastq"; fi
    for f in "${r1s[@]}" "${r2s[@]}"; do
        if [[ ($ext == .fastq.gz && $f != *.gz) || ($ext == .fastq && $f == *.gz) ]]; then
            die "Sample $s mixes gzipped and uncompressed FASTQs; cannot merge lanes"
        fi
    done

    IN_R1="$MERGE_DIR/${s}_R1$ext"
    IN_R2="$MERGE_DIR/${s}_R2$ext"
    IS_MERGED=1
    if [[ -s $IN_R1 && -s $IN_R2 ]] && is_done merge "$s"; then return; fi
    log "  merging ${#r1s[@]} FASTQ pairs for $s"
    mkdir -p "$MERGE_DIR"
    cat "${r1s[@]}" >"$IN_R1.tmp" && mv "$IN_R1.tmp" "$IN_R1"
    cat "${r2s[@]}" >"$IN_R2.tmp" && mv "$IN_R2.tmp" "$IN_R2"
    mark_done merge "$s"
}

step_trim() {
    local s=$1
    local -a qc=() gz=()
    if is_done trim "$s"; then log "  [skip] $s already trimmed"; return; fi
    prepare_inputs "$s"
    ((SKIP_FASTQC)) || qc=(--fastqc)
    # Trim Galore >= 2.0 gzips by default and deprecates --gzip.
    [[ $(tool_version trim_galore) =~ ^[01]\. ]] && gz=(--gzip)
    log "  trimming $s"
    run_logged "$LOG_DIR/$s.trim_galore.log" \
        trim_galore --paired \
        --quality "$TRIM_QUALITY" --phred33 \
        --stringency "$TRIM_STRINGENCY" \
        --length "$TRIM_LENGTH" \
        -e "$TRIM_ERROR" \
        ${gz[@]+"${gz[@]}"} --cores "$TG_CORES" \
        ${qc[@]+"${qc[@]}"} \
        --basename "$s" \
        -o "$TRIM_DIR" \
        "$IN_R1" "$IN_R2"
    [[ -s $TRIM_DIR/${s}_val_1.fq.gz && -s $TRIM_DIR/${s}_val_2.fq.gz ]] ||
        die "Trim Galore finished but trimmed files for $s are missing (see $LOG_DIR/$s.trim_galore.log)"
    if ((IS_MERGED)); then rm -f "$IN_R1" "$IN_R2" "$(done_marker merge "$s")"; fi
    mark_done trim "$s"
}

step_align() {
    local s=$1 r1 r2 prefix
    local -a rfc=() qm=() extra=()
    if is_done align "$s" && [[ -s $ALIGN_DIR/$s.bam ]]; then log "  [skip] $s already aligned"; return; fi

    IS_MERGED=0
    if ((SKIP_TRIM)); then
        prepare_inputs "$s"
        r1=$IN_R1
        r2=$IN_R2
    else
        r1="$TRIM_DIR/${s}_val_1.fq.gz"
        r2="$TRIM_DIR/${s}_val_2.fq.gz"
        [[ -s $r1 && -s $r2 ]] || die "Trimmed reads for $s not found ($r1). Rerun with --force."
    fi
    [[ $r1 == *.gz ]] && rfc=(--readFilesCommand zcat)
    [[ -f $STAR_INDEX/exonGeTrInfo.tab ]] && qm=(--quantMode GeneCounts)
    [[ -n $STAR_EXTRA ]] && read -r -a extra <<<"$STAR_EXTRA"

    mkdir -p "$ALIGN_DIR/$s"
    prefix="$ALIGN_DIR/$s/$s."
    rm -rf "${prefix}_STARtmp"
    log "  aligning $s"
    run_logged "$LOG_DIR/$s.STAR.log" \
        STAR --runThreadN "$THREADS" \
        --genomeDir "$STAR_INDEX" \
        --readFilesIn "$r1" "$r2" \
        ${rfc[@]+"${rfc[@]}"} \
        --outFileNamePrefix "$prefix" \
        --outSAMtype BAM SortedByCoordinate \
        --outBAMsortingThreadN "$SORT_THREADS" \
        --outSAMattrRGline "ID:$s" "SM:$s" \
        --outSAMstrandField intronMotif \
        --outSAMattributes NH HI AS NM MD XS \
        ${qm[@]+"${qm[@]}"} \
        ${extra[@]+"${extra[@]}"}
    mv "${prefix}Aligned.sortedByCoord.out.bam" "$ALIGN_DIR/$s.bam"
    rm -rf "${prefix}_STARtmp"
    run_logged "$LOG_DIR/$s.samtools_index.log" samtools index -@ "$THREADS" "$ALIGN_DIR/$s.bam"

    if ((IS_MERGED)); then rm -f "$IN_R1" "$IN_R2" "$(done_marker merge "$s")"; fi
    if ((DELETE_TRIMMED && !SKIP_TRIM)); then rm -f "$r1" "$r2"; fi
    mark_done align "$s"
}

strandedness_report() {
    # STAR ReadsPerGene columns: 2 unstranded, 3 = featureCounts -s 1, 4 = featureCounts -s 2
    local out="$COUNT_DIR/strandedness_check.tsv" s f guess mismatch=0
    printf 'sample\tunstranded\tforward_s1\treverse_s2\tsuggested_s\n' >"$out"
    for s in "${SAMPLE_ORDER[@]}"; do
        f="$ALIGN_DIR/$s/$s.ReadsPerGene.out.tab"
        [[ -f $f ]] || return 0
        guess=$(awk -v s="$s" 'NR>4{u+=$2; f+=$3; r+=$4}
            END{g=0; if(f+r>0){ if(f/(f+r)>0.8) g=1; else if(r/(f+r)>0.8) g=2 }
                printf "%s\t%d\t%d\t%d\t%d\n", s, u, f, r, g}' "$f")
        printf '%s\n' "$guess" >>"$out"
        [[ ${guess##*$'\t'} == "$STRAND" ]] || mismatch=1
    done
    log "  strandedness check written to $out"
    if ((mismatch)); then
        warn "Library strandedness suggested by STAR differs from --strandedness $STRAND for some samples:"
        column -t "$out" >&2 2>/dev/null || cat "$out" >&2
    fi
}

step_count() {
    local s fc_help sig out="$COUNT_DIR/featureCounts.txt"
    local -a bams=() pe=(-p) extra=()
    for s in "${SAMPLE_ORDER[@]}"; do bams+=("$ALIGN_DIR/$s.bam"); done
    sig="samples=${SAMPLE_ORDER[*]} gtf=$GTF s=$STRAND t=$FEATURE_TYPE g=$ATTR_TYPE extra=$FC_EXTRA"
    if is_done count all "$sig" && [[ -s $COUNT_DIR/gene_counts.tsv ]]; then
        log "  [skip] counts are up to date"
        return
    fi

    strandedness_report

    # subread >= 2.0.2 needs --countReadPairs to count fragments instead of reads.
    fc_help=$(featureCounts 2>&1 || true)
    [[ $fc_help == *--countReadPairs* ]] && pe+=(--countReadPairs)
    [[ -n $FC_EXTRA ]] && read -r -a extra <<<"$FC_EXTRA"

    log "  featureCounts on ${#bams[@]} BAM files"
    run_logged "$LOG_DIR/featureCounts.log" \
        featureCounts -T "$THREADS" "${pe[@]}" \
        -t "$FEATURE_TYPE" -g "$ATTR_TYPE" -s "$STRAND" \
        -a "$GTF" -o "$out" \
        ${extra[@]+"${extra[@]}"} \
        "${bams[@]}"

    # Clean matrix: gene_id + one column per sample, named by sample.
    awk -v names="${SAMPLE_ORDER[*]}" 'BEGIN{FS=OFS="\t"; n=split(names, a, " ")}
        /^#/ {next}
        !hdr {printf "gene_id"; for(i=1;i<=n;i++) printf "%s%s", OFS, a[i]; print ""; hdr=1; next}
        {printf "%s", $1; for(i=7;i<=NF;i++) printf "%s%s", OFS, $i; print ""}' \
        "$out" >"$COUNT_DIR/gene_counts.tsv"
    log "  count matrix: $COUNT_DIR/gene_counts.tsv"
    mark_done count all "$sig"
}

step_multiqc() {
    local d
    local -a dirs=()
    for d in "$TRIM_DIR" "$ALIGN_DIR" "$COUNT_DIR"; do [[ -d $d ]] && dirs+=("$d"); done
    log "  running MultiQC"
    mkdir -p "$MULTIQC_DIR"
    run_logged "$LOG_DIR/multiqc.log" \
        multiqc --force --outdir "$MULTIQC_DIR" --filename multiqc_report.html "${dirs[@]}"
    log "  report: $MULTIQC_DIR/multiqc_report.html"
}

# ----------------------------------------------------------------------------
# Argument parsing
# ----------------------------------------------------------------------------
LONG_OPTS="input:,samplesheet:,outdir:,star-index:,gtf:,threads:,strandedness:,feature-type:,attribute:,"
LONG_OPTS+="featurecounts-extra:,trim-quality:,trim-length:,trim-stringency:,trim-error-rate:,skip-trim,"
LONG_OPTS+="skip-fastqc,delete-trimmed,star-extra:,skip-multiqc,dry-run,force,check,help,version"
PARSED=$(getopt -o i:o:x:a:t:s:nhv --long "$LONG_OPTS" -n "$SCRIPT_NAME" -- "$@") || { usage >&2; exit 1; }
eval set -- "$PARSED"
while true; do
    case $1 in
        -i | --input)          INPUT_DIR=$2; shift 2 ;;
        --samplesheet)         SAMPLESHEET=$2; shift 2 ;;
        -o | --outdir)         OUTDIR=$2; shift 2 ;;
        -x | --star-index)     STAR_INDEX=$2; shift 2 ;;
        -a | --gtf)            GTF=$2; shift 2 ;;
        -t | --threads)        THREADS=$2; shift 2 ;;
        -s | --strandedness)   STRAND=$2; shift 2 ;;
        --feature-type)        FEATURE_TYPE=$2; shift 2 ;;
        --attribute)           ATTR_TYPE=$2; shift 2 ;;
        --featurecounts-extra) FC_EXTRA=$2; shift 2 ;;
        --trim-quality)        TRIM_QUALITY=$2; shift 2 ;;
        --trim-length)         TRIM_LENGTH=$2; shift 2 ;;
        --trim-stringency)     TRIM_STRINGENCY=$2; shift 2 ;;
        --trim-error-rate)     TRIM_ERROR=$2; shift 2 ;;
        --skip-trim)           SKIP_TRIM=1; shift ;;
        --skip-fastqc)         SKIP_FASTQC=1; shift ;;
        --delete-trimmed)      DELETE_TRIMMED=1; shift ;;
        --star-extra)          STAR_EXTRA=$2; shift 2 ;;
        --skip-multiqc)        SKIP_MULTIQC=1; shift ;;
        -n | --dry-run)        DRY_RUN=1; shift ;;
        --force)               FORCE=1; shift ;;
        --check)               CHECK_ONLY=1; shift ;;
        -h | --help)           usage; exit 0 ;;
        -v | --version)        echo "$SCRIPT_NAME $VERSION"; exit 0 ;;
        --)                    shift; break ;;
        *)                     die "Unexpected option: $1" ;;
    esac
done
(($# == 0)) || die "Unexpected positional arguments: $* (see --help)"

if ((CHECK_ONLY)); then
    echo "Software check:"
    if check_tools; then echo "All required software found."; exit 0; fi
    echo "Some required software is missing. Did you 'conda activate rnaseq'?" >&2
    exit 1
fi

# ----------------------------------------------------------------------------
# Validation
# ----------------------------------------------------------------------------
[[ -n $INPUT_DIR || -n $SAMPLESHEET ]] || { usage >&2; die "Provide --input or --samplesheet"; }
[[ -z $INPUT_DIR || -z $SAMPLESHEET ]] || die "Use either --input or --samplesheet, not both"
[[ -n $OUTDIR ]] || die "Missing required option --outdir"
[[ -n $STAR_INDEX ]] || die "Missing required option --star-index"
[[ $THREADS =~ ^[1-9][0-9]*$ ]] || die "--threads must be a positive integer (got '$THREADS')"
[[ $STRAND =~ ^[012]$ ]] || die "--strandedness must be 0, 1 or 2 (got '$STRAND')"

OUTDIR=$(abs_path "$OUTDIR")
STAR_INDEX=$(abs_dir "$STAR_INDEX")
for f in SA Genome genomeParameters.txt; do
    [[ -f $STAR_INDEX/$f ]] || die "$STAR_INDEX does not look like a STAR index (missing $f)"
done

if [[ -z $GTF ]]; then
    GTF=$(awk '$1 == "sjdbGTFfile" {print $2}' "$STAR_INDEX/genomeParameters.txt")
    [[ -n $GTF && $GTF != "-" ]] ||
        die "No --gtf given and the STAR index was built without a GTF. Please pass --gtf."
    [[ -f $GTF ]] || die "GTF recorded in the STAR index ($GTF) is not accessible. Please pass --gtf."
    GTF_SOURCE="auto-detected from STAR index"
else
    GTF_SOURCE="user supplied"
fi
[[ -f $GTF ]] || die "GTF not found: $GTF"
GTF=$(abs_path "$GTF")

if [[ -n $SAMPLESHEET ]]; then
    discover_from_samplesheet "$SAMPLESHEET"
    INPUT_DESC="samplesheet $(abs_path "$SAMPLESHEET")"
else
    INPUT_DIR=$(abs_dir "$INPUT_DIR")
    discover_from_dir "$INPUT_DIR"
    INPUT_DESC="$INPUT_DIR"
fi

TG_CORES=$((THREADS / 4))
((TG_CORES >= 1)) || TG_CORES=1
((TG_CORES <= 8)) || TG_CORES=8
SORT_THREADS=$((THREADS < 8 ? THREADS : 8))

TRIM_DIR="$OUTDIR/01_trimmed"
ALIGN_DIR="$OUTDIR/02_alignment"
COUNT_DIR="$OUTDIR/03_counts"
MULTIQC_DIR="$OUTDIR/04_multiqc"
MERGE_DIR="$OUTDIR/00_merged_fastq"
LOG_DIR="$OUTDIR/logs"
INFO_DIR="$OUTDIR/pipeline_info"

print_header() {
    echo "============================================================"
    echo " $SCRIPT_NAME v$VERSION"
    echo "============================================================"
    echo " Input        : $INPUT_DESC"
    echo " Output       : $OUTDIR"
    echo " STAR index   : $STAR_INDEX"
    echo " GTF          : $GTF ($GTF_SOURCE)"
    echo " Threads      : $THREADS"
    echo " Strandedness : $STRAND"
    echo " Trimming     : $( ((SKIP_TRIM)) && echo skipped || echo "q=$TRIM_QUALITY length=$TRIM_LENGTH stringency=$TRIM_STRINGENCY e=$TRIM_ERROR")"
    echo " Samples (${#SAMPLE_ORDER[@]}):"
    print_samples
    echo "============================================================"
}

if ((DRY_RUN)); then
    print_header
    echo "Software check:"
    check_tools || warn "Some required software is missing"
    echo "Dry run finished; nothing was executed."
    exit 0
fi

check_tools >/dev/null || { check_tools; die "Required software is missing. Did you activate the conda environment?"; }

# ----------------------------------------------------------------------------
# Run
# ----------------------------------------------------------------------------
mkdir -p "$OUTDIR" "$LOG_DIR" "$INFO_DIR" "$OUTDIR/.done" "$ALIGN_DIR" "$COUNT_DIR"
((SKIP_TRIM)) || mkdir -p "$TRIM_DIR"
((FORCE)) && rm -f "$OUTDIR"/.done/*

exec > >(tee -a "$OUTDIR/pipeline.log") 2>&1
trap 'die "Pipeline failed at line $LINENO. See logs in $LOG_DIR"' ERR

print_header
log "Pipeline started"
write_samplesheet "$INFO_DIR/samples.tsv"
{
    echo "command: $0 $*"
    echo "date: $(date)"
    echo "host: $(hostname)"
    for v in INPUT_DESC OUTDIR STAR_INDEX GTF THREADS STRAND FEATURE_TYPE ATTR_TYPE FC_EXTRA \
        TRIM_QUALITY TRIM_LENGTH TRIM_STRINGENCY TRIM_ERROR SKIP_TRIM SKIP_FASTQC STAR_EXTRA; do
        printf '%s: %s\n' "$v" "${!v}"
    done
    grep -E '^versionGenome' "$STAR_INDEX/genomeParameters.txt" | sed 's/^/star_index_/' || true
} >"$INFO_DIR/params.txt"
{
    echo "rnaseq_pipeline $VERSION"
    for t in STAR trim_galore cutadapt fastqc samtools featureCounts multiqc pigz; do
        command -v "$t" >/dev/null 2>&1 && printf '%s\t%s\n' "$t" "$(tool_version "$t")"
    done
} >"$INFO_DIR/software_versions.txt"

if ((SKIP_TRIM)); then
    log "Step 1/4: trimming skipped"
else
    log "Step 1/4: adapter/quality trimming (Trim Galore)"
    for s in "${SAMPLE_ORDER[@]}"; do step_trim "$s"; done
fi

log "Step 2/4: alignment (STAR) and BAM indexing"
for s in "${SAMPLE_ORDER[@]}"; do step_align "$s"; done

rmdir "$MERGE_DIR" 2>/dev/null || true

log "Step 3/4: gene counting (featureCounts)"
step_count

if ((SKIP_MULTIQC)); then
    log "Step 4/4: MultiQC skipped"
else
    log "Step 4/4: QC report (MultiQC)"
    step_multiqc
fi

touch "$INFO_DIR/pipeline.done"
log "Pipeline finished successfully"
