##jimmy_rna_featurecounts
rm(list = ls())
options(stringsAsFactors= F)
setwd('/home/wzx/project4t/QJC_SMARTSEQ_20210106/featurecounts/')
a=read.table('QJC_SMARTSEQ_ko_wt_22.txt', head=T)
n = grep(pattern = "PAR_Y",a$Geneid)
a <- a[-n,]
meta=a[,1:6]
exprSet=a[,7:ncol(a)]
colnames(exprSet)
colnames(exprSet)<-c("KO1","KO2","WT1","WT2")
rowname1<-gsub("\\..*$","",meta[,1])
rownames(exprSet)<-rowname1
write.table(exprSet, file = "exprSet_counts.txt", sep = "\t")
group_list=c('KO','KO','WT','WT')  #OR colname of specific table

library(corrplot)
library(pheatmap)
pheatmap(scale(cor(log2(exprSet+1))))

## hclust
#colnames(exprSet)=paste(group_list,1:ncol(exprSet),sep='_')
# Define nodePar
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19),
                cex = 0.7, col = "blue")
hc=hclust(dist(t(log2(exprSet+1))))
par(mar=c(5,5,5,10))
#png('hclust.png',res=120)
plot(as.dendrogram(hc), nodePar = nodePar, horiz = TRUE)
dev.off()

#analysis
##################################################################################
####################################DESEQ2########################################
##################################################################################
library("DESeq2","edgeR","limma")
suppressMessages(library(DESeq2))
colData <- data.frame(row.names = colnames(exprSet),
                       group_list=group_list)
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~ group_list)           #CALCULATE
dds <- DESeq(dds)

##diff_gene
res <-results(dds,
              contrast = c("group_list","KO","WT"))  ##TREAT / CONTROL  #READ RESULT
resOrdered <- res[order(res$padj),]
head(resOrdered)
DEG_result=as.data.frame(resOrdered)

#DEG_result = na.omit(DEG_result)  #TO DELETE NA ROW
nrDEG=DEG_result
##HEATMAP
library(pheatmap)
choose_gene=head(rownames(nrDEG),100)  ###50 maybe better
choose_matrix=exprSet[choose_gene,]
choose_matrix=t(scale(t(choose_matrix)))
pheatmap(choose_matrix,filename='DEG_TOP100_HEATMAP.png',show_rownames = F)

logFC_cutoff <- with(DEG_result, mean(abs(log2FoldChange))+ 2*sd(abs(log2FoldChange)))
logFC_cutoff=0.67
DEG_result$change = as.factor(ifelse(DEG_result$pvalue < 0.05 &abs(DEG_result$log2FoldChange) > logFC_cutoff,
                              ifelse(DEG_result$log2FoldChange > logFC_cutoff, "UP","DOWN"),'NOT')
)
write.table(DEG_result,"DEG_RESULT_new",sep = "\t")
this_tile <- paste0('Cutoff for log FC is ', logFC_cutoff)
##vocalno
library(ggplot2)
g = ggplot(data=DEG_result,aes(x=log2FoldChange,y=-log10(pvalue),
                               color=change))+
  geom_point(alpha=0.4, size=1.75)+ xlim(-15,15)+
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 fold change")+ylab("-log10 p-value")+
  ggtitle(this_tile)+theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red')) 
print(g)
ggsave(g,filename = 'volcano.pdf',device = "pdf",width= 5.5,height = 4)
#
pdf("qc_dispersions.pdf",height=4,width=4,pointsize = 0.1)
plotDispEsts(dds,main="Dispersion plot")
dev.off()
##normaliation
rld <- rlogTransformation(dds)
exprMatrix_rlog=assay(rld)
exprMatrix_rlog_frame<- as.data.frame(exprMatrix_rlog)
write.table(exprMatrix_rlog_frame,'exprMatrix.rlog.txt', sep = "\t")
#check_normalization
pdf("DEseq_RAWvsNORM.pdf",height=8,width=8)
par(cex=0.7)
n.sample=ncol(exprSet)
if(n.sample>40) par(cex=0.5)
cols <- rainbow(n.sample*1.2)
par(mfrow=c(2,2))
boxplot(exprSet,col=cols, main="expression value",las=2)
boxplot(exprMatrix_rlog,col=cols,main="expression value",las=2)
hist(as.matrix(exprSet))
hist(exprMatrix_rlog)
dev.off()
##########################################################################
#############################edgeR########################################
##########################################################################
library(edgeR)
d <- DGEList(counts=exprSet,group=factor(group_list))
keep <- rowSums(cpm(d)>1) >= 2
table(keep)
d <- d[keep, , keep.lib.sizes=FALSE]
d$samples$lib.size <- colSums(d$counts)
d <- calcNormFactors(d)
d$samples

design <- model.matrix(~0+factor(group_list))
dge=d
rownames(design)<-colnames(dge)
colnames(design)<-levels(factor(group_list))
dge <- estimateGLMCommonDisp(dge,design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

fit <- glmFit(dge, design)
# https://www.biostars.org/p/110861/
lrt <- glmLRT(fit,  contrast=c(0,1))    ##belong to design identification
nrDEG=topTags(lrt, n=nrow(dge))
nrDEG=as.data.frame(nrDEG)
head(nrDEG)
edgeR_nrDEG= nrDEG
write.table(edgeR_nrDEG,"edgeR_nrDEG.txt", sep = "\t")

##################################################################3
###################### then for limma/voom #######################
##################################################################
suppressMessages(library(limma))
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)
design
dge <- DGEList(counts=exprSet)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log=TRUE, prior.count=3)

v <- voom(dge,design,plot=TRUE, normalize="quantile")
fit <- lmFit(v, design)

group_list
cont.matrix=makeContrasts(contrasts=c('SOC-FT'),levels = design)
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)

tempOutput = topTable(fit2, coef='SOC-FT', n=Inf)
DEG_limma_voom = na.omit(tempOutput)
head(DEG_limma_voom)
write.table(DEG_limma_voom,"DEG_limma_voom.txt", sep = "\t")