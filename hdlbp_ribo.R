##hdlbp riboseq

setwd("/Volumes/landthaler/pcp/projects/miha/HDLBP/ribo/")
# setwd("E:/work/hdlbp/ribo/")

files<-list.files(getwd(),pattern=".*genes.*")
files<-files[!grepl("Mock", files)&!grepl("si", files)]
dfs<-lapply(files,read.delim,header=T)
names(dfs)<-gsub("\\..*","",files)

tab<-data.frame(sapply(dfs, function(x) cbind(x[,5])))
row.names(tab)<-dfs[[1]][,1]
colnames(tab)[1:2]<-c("293_1","293_2")

library(corrplot)
corrplot(cor(tab, method="spearman"),method="number",type="upper")
corrplot(cor(tab, method="pearson"),method="number",type="upper")


setwd("/Volumes/landthaler/pcp/projects/miha/HDLBP/frac/")
# setwd("E:/work/hdlbp/frac/")

files<-list.files(getwd(),pattern="T_.*genes.*")
dfs<-lapply(files,read.delim,header=T)
names(dfs)<-gsub("\\..*","",files)
rna<-data.frame(sapply(dfs, function(x) cbind(x[,5])))
row.names(rna)<-dfs[[1]][,1]

corrplot(cor(rna, method="spearman"),method="number",type="upper")
corrplot(cor(rna, method="pearson"),method="number",type="upper")


colnames(rna)<-paste0("rna_",colnames(tab)[1:6])
colnames(tab)<-paste0("ribo_",colnames(tab))

cnt<-merge(tab,rna,by="row.names")
row.names(cnt)<-cnt$Row.names
cnt<-cnt[,-1]

corrplot(cor(cnt))

library(ggplot2)
library(DESeq2)
cnt<-cnt[,colnames(cnt)[!grepl("no_CHX", colnames(cnt))]]
coldata<-data.frame(SampleID=colnames(cnt),
                    Condition=gsub("_.*","",gsub(".*;","",sub("_",";", colnames(cnt)))),
                    SeqType=c(rep("ribo",6),rep("rna",6)),
                    Batch=rep(c(1,2),6))
coldata$ConditionMerge<-ifelse(grepl("293", colnames(cnt)),"WT","KO")

coldata<-apply(coldata,2,as.factor)

dds<-DESeqDataSetFromMatrix(countData = round(cnt),
                            colData = coldata, 
                            design =~ SeqType + Condition +SeqType:Condition)
dds<-estimateSizeFactors(dds)
#sizeFactors(dds)<-c(estimateSizeFactorsForMatrix(cnt[,1:6]), estimateSizeFactorsForMatrix(rna) )

vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Condition", "SeqType"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=SeqType)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 
ggplot(pcaData[grepl("ribo",pcaData$group),], aes(PC1, PC2, color=Condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 
ggplot(pcaData[grepl("rna",pcaData$group),], aes(PC1, PC2, color=Condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 

vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Condition", "SeqType","Batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData[grepl("ribo",pcaData$group),], aes(PC1, PC2, color=Condition, shape=Batch)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 

dds$Condition = relevel(dds$Condition,"293") #everything relative to Mock
dds$SeqType = relevel(dds$SeqType,"rna") #everything relative to RNAseq
dds<-DESeq(dds, test="LRT", reduced=~SeqType+Condition)

resultsNames(dds)
# Condition2.SeqTypeRibo.seq means Changes in Ribo-seq levels in Condition2 vs Condition1 accounting for changes in RNA-seq levels in Condition2 vs Condition1
res <- results(dds, contrast=list("SeqTyperibo.Conditionguide1")) 
summary(res)
DESeq2::plotMA(res)
lfc.ribo.rna.guide1.293<-data.frame(res[,c(1,2,6)])
colnames(lfc.ribo.rna.guide1.293)<-paste0(colnames(lfc.ribo.rna.guide1.293),".ribo.rna.guide1.293")

res <- results(dds, contrast=list("SeqTyperibo.Conditionguide2")) 
summary(res)
DESeq2::plotMA(res)
DESeq2::plotDispEsts(dds)
lfc.ribo.rna.guide2.293<-data.frame(res[,c(1,2,6)])
colnames(lfc.ribo.rna.guide2.293)<-paste0(colnames(lfc.ribo.rna.guide2.293),".ribo.rna.guide2.293")


dds<-DESeqDataSetFromMatrix(countData = round(cnt),
                            colData = coldata, 
                            design =~  SeqType+ConditionMerge+ SeqType:ConditionMerge)

dds<-estimateSizeFactors(dds)
#sizeFactors(dds)<-c(estimateSizeFactorsForMatrix(rpf), estimateSizeFactorsForMatrix(rna) )


dds$ConditionMerge = relevel(dds$ConditionMerge,"WT") #everything relative to wt_as
dds$SeqType = relevel(dds$SeqType,"rna") #everything relative to RNAseq
dds<-DESeq(dds)
resultsNames(dds)
# Condition2.SeqTypeRibo.seq means Changes in Ribo-seq levels in Condition2 vs Condition1 accounting for changes in RNA-seq levels in Condition2 vs Condition1
res <- results(dds, contrast=list("SeqTyperibo.ConditionMergeKO")) 
summary(res)
DESeq2::plotMA(res)
lfc.ribo.rna.KO.WT<-data.frame(res[,c(1,2,6)])
colnames(lfc.ribo.rna.KO.WT)<-paste0(colnames(lfc.ribo.rna.KO.WT),".ribo.rna.KO.WT")



cnt<-cnt[,colnames(cnt)[grepl("ribo", colnames(cnt))]]
coldata<-data.frame(SampleID=colnames(cnt),
                    Condition=gsub("_.*","",gsub(".*;","",sub("_",";", colnames(cnt)))),
                    SeqType=rep("ribo",6),
                    Batch=rep(c(1,2),3))
coldata$ConditionMerge<-ifelse(grepl("293", colnames(cnt)),"WT","KO")

coldata<-apply(coldata,2,as.factor)

dds<-DESeqDataSetFromMatrix(countData = round(cnt),
                            colData = coldata, 
                            design =~ Condition)
dds<-estimateSizeFactors(dds)

dds$Condition = relevel(dds$Condition,"293") #everything relative to Mock
dds<-DESeq(dds)

resultsNames(dds)
res <- results(dds, contrast=c("Condition","guide1","293")) 
summary(res)
DESeq2::plotMA(res)
lfc.rpf.guide1.293<-data.frame(res[,c(1,2,6)])
colnames(lfc.rpf.guide1.293)<-paste0(colnames(lfc.rpf.guide1.293),".rpf.guide1.293")

res <- results(dds, contrast=c("Condition","guide2","293")) 
summary(res)
DESeq2::plotMA(res)
lfc.rpf.guide2.293<-data.frame(res[,c(1,2,6)])
colnames(lfc.rpf.guide2.293)<-paste0(colnames(lfc.rpf.guide2.293),".rpf.guide2.293")

dds<-DESeqDataSetFromMatrix(countData = round(cnt),
                            colData = coldata, 
                            design =~ ConditionMerge)
dds<-estimateSizeFactors(dds)

dds$ConditionMerge = relevel(dds$ConditionMerge,"WT") #everything relative to Mock
dds<-DESeq(dds)

resultsNames(dds)
res <- results(dds, contrast=c("ConditionMerge","KO","WT")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.rpf.KO.WT<-data.frame(res[,c(1,2,6)])
colnames(lfc.rpf.KO.WT)<-paste0(colnames(lfc.rpf.KO.WT),".rpf.KO.WT")


###tpm
setwd("/Volumes/landthaler/pcp/projects/miha/HDLBP/ribo/")
# setwd("E:/work/hdlbp/ribo/")

files<-list.files(getwd(),pattern=".*genes.*")[1:6]
dfs<-lapply(files,read.delim,header=T)
names(dfs)<-gsub("\\..*","",files)

tab<-data.frame(sapply(dfs, function(x) cbind(x[,6])))
row.names(tab)<-dfs[[1]][,1]
colnames(tab)[1:2]<-c("293_1","293_2")

library(corrplot)
corrplot(cor(tab, method="spearman"),method="number",type="upper")
corrplot(cor(tab, method="pearson"),method="number",type="upper")


setwd("/Volumes/landthaler/pcp/projects/miha/HDLBP/frac/")
# setwd("E:/work/hdlbp/frac/")

files<-list.files(getwd(),pattern="T_.*genes.*")
dfs<-lapply(files,read.delim,header=T)
names(dfs)<-gsub("\\..*","",files)
rna<-data.frame(sapply(dfs, function(x) cbind(x[,6])))
row.names(rna)<-dfs[[1]][,1]

corrplot(cor(rna, method="spearman"),method="number",type="upper")
corrplot(cor(rna, method="pearson"),method="number",type="upper")


colnames(rna)<-paste0("rna_",colnames(tab))
colnames(tab)<-paste0("ribo_",colnames(tab))

rna[rna==0]<-NA
te<-tab/rna
colnames(te)<-gsub("ribo_","te_",colnames(te))

tpm<-cbind(tab, rna, te)




###knockdown
setwd("/Volumes/landthaler/pcp/projects/miha/HDLBP/ribo/")
# setwd("E:/work/hdlbp/ribo/")

files<-list.files(getwd(),pattern=".*genes.*")
files<-files[grepl("Mock", files)|grepl("si", files)]
dfs<-lapply(files,read.delim,header=T)
names(dfs)<-gsub("\\..*","",files)

tab<-data.frame(sapply(dfs, function(x) cbind(x[,5])))
row.names(tab)<-dfs[[1]][,1]

library(corrplot)
corrplot(cor(tab, method="spearman"),method="number",type="upper")
corrplot(cor(tab, method="pearson"),method="number",type="upper")

setwd("/Volumes/landthaler/pcp/projects/miha/HDLBP/frac/")
# setwd("E:/work/hdlbp/frac/")

files<-list.files(getwd(),pattern="^Mock|^si")
dfs<-lapply(files,read.delim,header=T)
names(dfs)<-gsub("\\..*","",files)
rna<-data.frame(sapply(dfs, function(x) cbind(x[,5])))
row.names(rna)<-dfs[[1]][,1]

corrplot(cor(rna, method="spearman"),method="number",type="upper")
corrplot(cor(rna, method="pearson"),method="number",type="upper")


colnames(rna)<-paste0("rna_",colnames(rna))
colnames(tab)<-paste0("ribo_",colnames(tab))

cnt<-merge(tab,rna,by="row.names")
row.names(cnt)<-cnt$Row.names
cnt<-cnt[,-1]

corrplot(cor(cnt), type="upper", method="square")

library(ggplot2)
library(DESeq2)

coldata<-data.frame(SampleID=colnames(cnt),
                    Condition=ifelse(grepl("Mock", colnames(cnt)),"Mock",gsub("_.*","",gsub(".*;","",sub("_",";", colnames(cnt))))),
                    SeqType=c(rep("ribo",6),rep("rna",6)),
                    Batch=rep(c(1,2),6))
coldata$ConditionMerge<-ifelse(grepl("Mock", colnames(cnt)),"mock","si")

coldata<-apply(coldata,2,as.factor)

dds<-DESeqDataSetFromMatrix(countData = round(cnt),
                            colData = coldata, 
                            design =~ SeqType + Condition +SeqType:Condition)
# dds<-estimateSizeFactors(dds)
sizeFactors(dds)<-c(estimateSizeFactorsForMatrix(cnt[,1:6]), estimateSizeFactorsForMatrix(cnt[,7:12]) )

vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Condition", "SeqType"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=SeqType)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 
ggplot(pcaData[grepl("ribo",pcaData$group),], aes(PC1, PC2, color=Condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 
ggplot(pcaData[grepl("rna",pcaData$group),], aes(PC1, PC2, color=Condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 




dds$Condition = relevel(dds$Condition,"Mock") #everything relative to Mock
dds$SeqType = relevel(dds$SeqType,"rna") #everything relative to RNAseq
dds<-DESeq(dds, test="LRT", reduced=~SeqType+Condition)

resultsNames(dds)
# Condition2.SeqTypeRibo.seq means Changes in Ribo-seq levels in Condition2 vs Condition1 accounting for changes in RNA-seq levels in Condition2 vs Condition1
res <- results(dds, contrast=list("SeqTyperibo.Conditionsi1")) 
summary(res)
DESeq2::plotMA(res)
lfc.ribo.rna.si1.mock<-data.frame(res[,c(1,2,6)])
colnames(lfc.ribo.rna.si1.mock)<-paste0(colnames(lfc.ribo.rna.si1.mock),".ribo.rna.si1.mock")

res <- results(dds, contrast=list("SeqTyperibo.Conditionsi2")) 
summary(res)
DESeq2::plotMA(res)
lfc.ribo.rna.si2.mock<-data.frame(res[,c(1,2,6)])
colnames(lfc.ribo.rna.si2.mock)<-paste0(colnames(lfc.ribo.rna.si2.mock),".ribo.rna.si2.mock")


dds<-DESeqDataSetFromMatrix(countData = round(cnt),
                            colData = coldata, 
                            design =~ SeqType + ConditionMerge +SeqType:ConditionMerge)
# dds<-estimateSizeFactors(dds)
sizeFactors(dds)<-c(estimateSizeFactorsForMatrix(cnt[,1:6]), estimateSizeFactorsForMatrix(cnt[,7:12]) )

dds$ConditionMerge = relevel(dds$ConditionMerge,"mock") #everything relative to Mock
dds$SeqType = relevel(dds$SeqType,"rna") #everything relative to RNAseq
dds<-DESeq(dds, test="LRT", reduced=~SeqType+ConditionMerge)

resultsNames(dds)
# Condition2.SeqTypeRibo.seq means Changes in Ribo-seq levels in Condition2 vs Condition1 accounting for changes in RNA-seq levels in Condition2 vs Condition1
res <- results(dds, contrast=list("SeqTyperibo.ConditionMergesi")) 
summary(res)
DESeq2::plotMA(res)
lfc.ribo.rna.si.Mock<-data.frame(res[,c(1,2,6)])
colnames(lfc.ribo.rna.si.Mock)<-paste0(colnames(lfc.ribo.rna.si.Mock),".ribo.rna.si.Mock")


cnt<-cnt[,colnames(cnt)[grepl("ribo_", colnames(cnt))]]
coldata<-data.frame(SampleID=colnames(cnt),
                    Condition=ifelse(grepl("Mock", colnames(cnt)),"Mock",gsub("_.*","",gsub(".*;","",sub("_",";", colnames(cnt))))),
                    SeqType=rep("ribo",6),
                    Batch=rep(c(1,2),3))
coldata$ConditionMerge<-ifelse(grepl("Mock", colnames(cnt)),"mock","si")

coldata<-apply(coldata,2,as.factor)

dds<-DESeqDataSetFromMatrix(countData = round(cnt),
                            colData = coldata, 
                            design =~ Condition)

vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Condition", "Batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Batch)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 

dds$Condition = relevel(dds$Condition,"Mock") #everything relative to Mock
dds<-DESeq(dds)

resultsNames(dds)
# Condition2.SeqTypeRibo.seq means Changes in Ribo-seq levels in Condition2 vs Condition1 accounting for changes in RNA-seq levels in Condition2 vs Condition1
res <- results(dds, contrast=c("Condition", "si1", "Mock")) 
summary(res)
DESeq2::plotMA(res)
lfc.rpf.si1.mock<-data.frame(res[,c(1,2,6)])
colnames(lfc.rpf.si1.mock)<-paste0(colnames(lfc.rpf.si1.mock),".rpf.si1.mock")

res <- results(dds, contrast=c("Condition", "si2", "Mock")) 
summary(res)
DESeq2::plotMA(res)
lfc.rpf.si2.mock<-data.frame(res[,c(1,2,6)])
colnames(lfc.rpf.si2.mock)<-paste0(colnames(lfc.rpf.si2.mock),".rpf.si2.mock")



dds<-DESeqDataSetFromMatrix(countData = round(cnt),
                            colData = coldata, 
                            design =~ ConditionMerge)

dds<-estimateSizeFactors(dds)
dds$ConditionMerge = relevel(dds$ConditionMerge,"mock") #everything relative to Mock
dds<-DESeq(dds)


resultsNames(dds)
# Condition2.SeqTypeRibo.seq means Changes in Ribo-seq levels in Condition2 vs Condition1 accounting for changes in RNA-seq levels in Condition2 vs Condition1
res <- results(dds, contrast=c("ConditionMerge", "si", "mock")) 
summary(res)
DESeq2::plotMA(res)
lfc.rpf.si.mock<-data.frame(res[,c(1,2,6)])
colnames(lfc.rpf.si.mock)<-paste0(colnames(lfc.rpf.si.mock),".rpf.si.mock")

# riboResults<-cbind(tpm, lfc.ribo.rna.guide1.293, lfc.ribo.rna.guide2.293, 
#                    lfc.ribo.rna.KO.WT, lfc.rpf.guide1.293, lfc.rpf.guide2.293,
#                    lfc.rpf.KO.WT,
#                    lfc.ribo.rna.si1.mock, lfc.ribo.rna.si2.mock,
#                    lfc.ribo.rna.si.Mock, 
#                    lfc.rpf.si1.mock,lfc.rpf.si2.mock,lfc.rpf.si.mock)

riboResults<-cbind(lfc.ribo.rna.si1.mock, lfc.ribo.rna.si2.mock,
                   lfc.ribo.rna.si.Mock, 
                   lfc.rpf.si1.mock,lfc.rpf.si2.mock,lfc.rpf.si.mock)

rel<-subset(riboResults, rowMeans(riboResults[,c("rna_293_1","rna_293_2")])>=10)
corrplot(cor(rel[,13:18], use = "pairwise.complete.obs"))
corrplot(cor(rel[,colnames(rel)[grepl("log2Fold",colnames(rel))]], use = "pairwise.complete.obs"))






#fin variable from fractionation script
fin<-read.delim("~/Google Drive/hdlbp/master_table_hdlbp.txt", header=T)
fin<-merge(fin, riboResults, by.x="gene_id", by.y="row.names") 
library(ggplot2)

rel<-subset(fin, tpm_cutoff>=10)
ggplot(rel, aes(log2FoldChange.rpf.KO.WT, colour=changeMem_enrichment.KO))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))
ggplot(rel, aes(log2FoldChange.rpf.si.mock, colour=changeMem_enrichment.KO))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))
ggplot(rel, aes(log2FoldChange.ribo.rna.si.Mock, colour=changeMem_enrichment.KO))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))

ggplot(rel, aes(log2FoldChange.ribo.rna.KO.WT, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))

ggplot(rel, aes(log2FoldChange.rpf.KO.WT, colour=changeMem_enrichment.KO))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))

ggplot(rel, aes(log2FoldChange.rpf.KO.WT, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))

wilcox.test(rel[rel$changeMem_enrichment.KO=="downMem_enrichment.KO","log2FoldChange.rpf.KO.WT"],
            rel[rel$changeMem_enrichment.KO=="all_membrane_localized","log2FoldChange.rpf.KO.WT"])

ggplot(rel, aes(log2FoldChange.tot.cyt.KO.293, colour=changeMem_enrichment.KO))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))

#write.table(fin, "master_table_hdlbp.txt", quote=F, sep="\t", row.names=F, col.names = T)
