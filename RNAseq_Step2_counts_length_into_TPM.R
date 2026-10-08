##counts_length_into_TPM
###############counts_rpkm##################
mycounts<-read.table("/home/wzx/project16t_new/Personal_cohort_PB_SOC_FT/RNA-seq/Merge_featurecounts/merge_add_ref_48_samples_featurecounts.txt", header = T)
n = grep(pattern = "PAR_Y",mycounts$Geneid)
mycounts <- mycounts[-n,]
name1 <- mycounts[,1]
name2 <- gsub("\\..*$","",name1)
rownames(mycounts) <- name2
mycounts <- mycounts[,-1:-5]
kb <- mycounts$Length /1000
countdata <- mycounts[,2:ncol(mycounts)]
rpk<-countdata / kb
tpm <- t(t(rpk)/colSums(rpk)*1000000)
# fpkm <- t(t(rpk)/colSums(countdata)*1000000)
#ID_NAME
ENSG_NAME <- row.names(tpm)
library("clusterProfiler")
library(org.Hs.eg.db)
gene.df <- bitr(ENSG_NAME, fromType = "ENSEMBL", #fromType是指你的数据ID类型是属于哪一类的
                toType = c("SYMBOL"), #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                OrgDb = org.Hs.eg.db)#Orgdb是指对应的注释包是哪个
head(gene.df)
gene.df <- data.frame (gene.df)
mycounts <- cbind(rownames(mycounts),mycounts)
merge<- merge(mycounts,gene.df, by.y="ENSEMBL", by.x="rownames(mycounts)",all.x= TRUE)
# fpkm <- data.frame(fpkm)
# fpkm <- cbind(rownames(fpkm),fpkm)
tpm <- data.frame(tpm)
tpm <- cbind(rownames(tpm),tpm)
merge<- merge(merge,tpm, by.y="rownames(tpm)", by.x="rownames(mycounts)",all.x= TRUE)
##diff
degseq2<- read.table("./DEG_RESULT_new")
degseq2 <- cbind(rownames(degseq2),degseq2)
# n = grep(pattern = "PAR_Y",degseq2$`rownames(degseq2)`)
# degseq2 <- degseq2[-n,]
name1 <- degseq2[,1]
name2 <- gsub("\\..*$","",name1)
rownames(degseq2) <- name2
degseq2 <- degseq2[,-1]
degseq2 <- cbind(rownames(degseq2),degseq2)
merge<- merge(merge,degseq2, by.y="rownames(degseq2)", by.x="rownames(mycounts)",all.x= TRUE)
write.table(merge,"./merge_tpm_SOC_FT_ALL_20240111.txt", sep = "\t", row.names = FALSE)












#simple_method
setwd("~/Project/ESR1_COUNT_FPKM/")
mycounts <- read.table("GSE115481_Feature_Counts_All.txt", header = T )
rownames(mycounts) <- mycounts[,1]
mycounts <- mycounts[,-1]

kb <- mycounts$Length /1000
kb
countdata <- mycounts[,2:4]
rpk<-countdata / kb
tpm <- t(t(rpk)/colSums(rpk)*1000000)
head(tpm)

fpkm <- t(t(rpk)/colSums(countdata)*1000000)
write.table(fpkm,file="esr1_fpkm.txt",sep="\t",quote = F)


#calculate_whole_exon
setwd("~/home/wzx/Project/FPKM_COUNTS_GEO_SOC_20200425")
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("/home/wzx/reference_total_genome_rna_gtf/gencode.v19.annotation.gtf",format = "gtf")
exons_gene <- exonsBy(txdb, by = "gene")
exons_gene_lens <- lapply(exons_gene,function(x){sum(width(reduce(x)))})
lengthdata<- as.data.frame(exons_gene_lens)
t1 <- t(data.frame(lengthdata,row.names = "length"))
t1 <- cbind(rownames(t1),t1)
length_exons <- data.frame(t1)

t1$newname = row.names(t1)
write.table(t1,"/home/wzx/Project/FPKM_COUNTS_GEO_SOC_20200425/exonlenth_hg19.txt", sep = "\t")
exonlength=read.table("/home/wzx/Project/FPKM_COUNTS_GEO_SOC_20200425/exonlenth_hg19.txt")
#mycounts<-read.table("/home/wzx/Project/FPKM_COUNTS_GEO_SOC_20200425/counts/GEO_FT_SOC_33.txt", header = T)
#n = grep(pattern = "PAR_Y",mycounts$Geneid)
#mycounts <- mycounts[-n,]
#name1 <- mycounts[,1]
#name2 <- gsub("_.*$","",name1)
#rownames(mycounts) <- name2
#mycounts <- mycounts[,-1:-5]
#mycounts <- cbind(rownames(mycounts),mycounts)
#merge<- merge(length_exons, mycounts, by.x="V1", by.y="rownames(mycounts)")
########FUCK!!!###exons in the counts files is the same as recalcualte of the sum of exons###########

