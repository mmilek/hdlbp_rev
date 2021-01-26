library(ggplot2)
library(reshape2)
library(corrplot)
library(DESeq2)
setwd("/Volumes/landthaler/pcp/projects/miha/HDLBP/frac/")

setwd("D:/landthaler/HDLBP/frac/")
# setwd("E:/work/hdlbp/frac/")


files<-list.files(getwd(), pattern="genes")

dat<-lapply(files, read.delim)
names(dat)<-gsub("\\..*","",files)

tab<-data.frame(sapply(dat, function(x) cbind(x[,5])))
row.names(tab)<-dat[[1]][,1]

corrplot(cor(tab), type="upper", method="number")


tpm<-data.frame(sapply(dat, function(x) cbind(x[,6])))
row.names(tpm)<-dat[[1]][,1]

corrplot(cor(tpm), type="upper", method="square")

##PCA
subs<-subset(tab, select=colnames(tab)[grepl("[CMT]_",colnames(tab))])
# subs<-tab[,c(1:14,17:22,27:ncol(tab)) ]
coldata<-data.frame(SampleID=colnames(subs),
                    Condition=sub(";","_",gsub("_.*","",sub("_",";",colnames(subs)))),
                    Batch=substr(colnames(subs), nchar(colnames(subs)), nchar(colnames(subs))),
                    Fraction=gsub("_.*","",colnames(subs)))

coldata<-apply(coldata,2,as.factor)
# coldata[7,2]<-"T_crispr"
# coldata[8,2]<-"T_crispr"

corrplot(cor(subs, use="pairwise.complete.obs"), type="upper", method="color",tl.col = "black")




dds<-DESeqDataSetFromMatrix(countData = round(subs),
                            colData = coldata, 
                            design = ~Batch+Condition)


vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Fraction", "Batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Fraction, shape=Batch)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 

subs<-subset(tab, select=colnames(tab)[grepl("[CMT]_293",colnames(tab))])
# subs<-tab[,c(1:14,17:22,27:ncol(tab)) ]
coldata<-data.frame(SampleID=colnames(subs),
                    Condition=sub(";","_",gsub("_.*","",sub("_",";",colnames(subs)))),
                    Batch=substr(colnames(subs), nchar(colnames(subs)), nchar(colnames(subs))),
                    Fraction=gsub("_.*","",colnames(subs)))

coldata<-apply(coldata,2,as.factor)
# coldata[7,2]<-"T_crispr"
# coldata[8,2]<-"T_crispr"

corrplot(cor(subs, use="pairwise.complete.obs"), type="upper", method="color",tl.col = "black", addCoef.col = "white")




dds<-DESeqDataSetFromMatrix(countData = round(subs),
                            colData = coldata, 
                            design = ~Batch+Condition)


vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Fraction", "Batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Fraction, shape=Batch)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 






membr<-tab[,colnames(tab)[grepl("M_",colnames(tab))]]
coldata<-data.frame(SampleID=colnames(membr),
                    Condition=sub(";","_",gsub("_.*","",sub("_",";",colnames(membr)))),
                    Batch=substr(colnames(membr), nchar(colnames(membr)), nchar(colnames(membr))))

coldata<-apply(coldata,2,as.factor)

dds<-DESeqDataSetFromMatrix(countData = round(membr),
                            colData = coldata, 
                            design = ~Batch+Condition)


vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Condition", "Batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Batch)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 


cytos<-tab[,colnames(tab)[grepl("C_",colnames(tab))]]
coldata<-data.frame(SampleID=colnames(cytos),
                    Condition=sub(";","_",gsub("_.*","",sub("_",";",colnames(cytos)))),
                    Batch=substr(colnames(cytos), nchar(colnames(cytos)), nchar(colnames(cytos))))

coldata<-apply(coldata,2,as.factor)

dds<-DESeqDataSetFromMatrix(countData = round(cytos),
                            colData = coldata, 
                            design = ~Batch+Condition)


vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Condition", "Batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Batch)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 





nucs<-tab[,colnames(tab)[grepl("N_",colnames(tab))]]
coldata<-data.frame(SampleID=colnames(nucs),
                    Condition=sub(";","_",gsub("_.*","",sub("_",";",colnames(nucs)))),
                    Batch=substr(colnames(nucs), nchar(colnames(nucs)), nchar(colnames(nucs))))

coldata<-apply(coldata,2,as.factor)

dds<-DESeqDataSetFromMatrix(countData = round(nucs),
                            colData = coldata, 
                            design = ~Batch+Condition)


vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Condition", "Batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Batch)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 



tots<-tab[,colnames(tab)[grepl("T_",colnames(tab))]]
coldata<-data.frame(SampleID=colnames(tots),
                    Condition=sub(";","_",gsub("_.*","",sub("_",";",colnames(tots)))),
                    Batch=substr(colnames(tots), nchar(colnames(tots)), nchar(colnames(tots))))

coldata<-apply(coldata,2,as.factor)

dds<-DESeqDataSetFromMatrix(countData = round(tots),
                            colData = coldata, 
                            design = ~Batch+Condition)


vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Condition", "Batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Batch)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 


##add other info

genes<-read.table("/Volumes/landthaler/pcp/projects/miha/DDR/nascent/RSEM/geneInfo.txt", header=T)[,c(1,4:5,10,12)]
# genes<-read.table("geneInfo.txt", header=T)[,c(1,4:5,10,12)]

tar<-read.delim("neel_HDLBP2_targets.txt", header=T)
colnames(tar)[1]<-"gene_name"
tar$target<-"target"
sig<-read.csv("mart_export.txt", header=T)


###attempt to make final table with fold changes etc
fin<-merge(tpm, genes, by.x="row.names", by.y= "Gene", all.x=T)
fin<-merge(fin, sig, by.x="Symbol", by.y= "Gene.name", all.x=T)[,-c(43,44)]
fin<-merge(fin, tar, by.x="Symbol", by.y= "gene_name", all.x=T)

fin<-subset(fin, !duplicated(Row.names))

fin$target<-ifelse(is.na(fin$target), "nontarget", "target")

fin$tpm_cutoff<-rowMeans(fin[,c("T_293_1" ,"T_293_2")])
# fin<-subset(fin, tpm_cutoff>=10)


fin$Transmembrane.helices<-ifelse(is.na(fin$Transmembrane.helices), "noTM", 
                                  ifelse(fin$Transmembrane.helices=="", "noTM", "TMhelix"))

fin$Cleavage.site..Signalp.<-ifelse(is.na(fin$Cleavage.site..Signalp.), "noSignalP", 
                                    ifelse(fin$Cleavage.site..Signalp.=="", "noSignalP", "SignalP"))

fin$Low.complexity..Seg.<-ifelse(is.na(fin$Low.complexity..Seg.), "noLoComp", 
                                 ifelse(fin$Low.complexity..Seg.=="", "noLoComp", "LoComp"))


rel<-subset(fin, tpm_cutoff>1)
mel<-melt(rel, measure.vars = colnames(fin[,3:34]), id.vars = c("Symbol", "Annotation","target","Transmembrane.helices", 
                                                                "Cleavage.site..Signalp.","Low.complexity..Seg."))
ggplot(mel, aes(variable, log2(value), fill=Transmembrane.helices))+geom_boxplot()
ggplot(mel, aes(variable, log2(value), fill=Cleavage.site..Signalp.))+geom_boxplot()


####deseq simple comparisons within fractions
cnt<-tab[,colnames(tab)[grepl("M_", colnames(tab))]]

coldata<-data.frame(SampleID=colnames(cnt),
                    Condition=ifelse(grepl("293", colnames(cnt)),"WT","KO"),
                    Fraction=gsub("_.*","",colnames(cnt)),
                    Batch=gsub(".*_","",colnames(cnt)))
coldata<-apply(coldata,2,as.factor)

dds<-DESeqDataSetFromMatrix(countData = round(cnt),
                            colData = coldata, 
                            design =~ Condition )

dds<-estimateSizeFactors(dds)

dds$Condition = relevel(dds$Condition,"WT")
dds<-DESeq(dds)

resultsNames(dds)
res <- results(dds, contrast=c("Condition","KO","WT")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.mem.KO.WT<-data.frame(res[,c(1,2,6)])
colnames(lfc.mem.KO.WT)<-paste0(colnames(lfc.mem.KO.WT),".mem.KO.WT")


cnt<-tab[,colnames(tab)[grepl("C_", colnames(tab))]]

coldata<-data.frame(SampleID=colnames(cnt),
                    Condition=ifelse(grepl("293", colnames(cnt)),"WT","KO"),
                    Fraction=gsub("_.*","",colnames(cnt)),
                    Batch=gsub(".*_","",colnames(cnt)))
coldata<-apply(coldata,2,as.factor)

dds<-DESeqDataSetFromMatrix(countData = round(cnt),
                            colData = coldata, 
                            design =~ Condition )

dds<-estimateSizeFactors(dds)

dds$Condition = relevel(dds$Condition,"WT")
dds<-DESeq(dds)

resultsNames(dds)
res <- results(dds, contrast=c("Condition","KO","WT")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.cyt.KO.WT<-data.frame(res[,c(1,2,6)])
colnames(lfc.cyt.KO.WT)<-paste0(colnames(lfc.cyt.KO.WT),".cyt.KO.WT")


cnt<-tab[,colnames(tab)[grepl("N_", colnames(tab))]]

coldata<-data.frame(SampleID=colnames(cnt),
                    Condition=ifelse(grepl("293", colnames(cnt)),"WT","KO"),
                    Fraction=gsub("_.*","",colnames(cnt)),
                    Batch=gsub(".*_","",colnames(cnt)))
coldata<-apply(coldata,2,as.factor)

dds<-DESeqDataSetFromMatrix(countData = round(cnt),
                            colData = coldata, 
                            design =~ Condition )

dds<-estimateSizeFactors(dds)

dds$Condition = relevel(dds$Condition,"WT")
dds<-DESeq(dds)

resultsNames(dds)
res <- results(dds, contrast=c("Condition","KO","WT")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.nuc.KO.WT<-data.frame(res[,c(1,2,6)])
colnames(lfc.nuc.KO.WT)<-paste0(colnames(lfc.nuc.KO.WT),".nuc.KO.WT")

cnt<-tab[,colnames(tab)[grepl("T_", colnames(tab))]]

coldata<-data.frame(SampleID=colnames(cnt),
                    Condition=ifelse(grepl("293", colnames(cnt)),"WT","KO"),
                    Fraction=gsub("_.*","",colnames(cnt)),
                    Batch=gsub(".*_","",colnames(cnt)))
coldata<-apply(coldata,2,as.factor)

dds<-DESeqDataSetFromMatrix(countData = round(cnt),
                            colData = coldata, 
                            design =~ Condition )

dds<-estimateSizeFactors(dds)

dds$Condition = relevel(dds$Condition,"WT")
dds<-DESeq(dds)

resultsNames(dds)
res <- results(dds, contrast=c("Condition","KO","WT")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.tot.KO.WT<-data.frame(res[,c(1,2,6)])
colnames(lfc.tot.KO.WT)<-paste0(colnames(lfc.tot.KO.WT),".tot.KO.WT")



dese<-cbind(lfc.mem.KO.WT, lfc.cyt.KO.WT, lfc.nuc.KO.WT, lfc.tot.KO.WT)

fin<-merge(fin, dese, by.x="Row.names", by.y="row.names")


###deseq pairwise between fractions
cnt<-tab[,colnames(tab)[grepl("C_", colnames(tab))|grepl("M_", colnames(tab))]]


coldata<-data.frame(SampleID=colnames(cnt),
                    Condition=sub(";","_",gsub("_.*","",sub("_",";",colnames(cnt)))),
                    Fraction=gsub("_.*","",colnames(cnt)),
                    Batch=gsub(".*_","",colnames(cnt)))
coldata<-apply(coldata,2,as.factor)

dds<-DESeqDataSetFromMatrix(countData = round(cnt),
                            colData = coldata, 
                            design =~ Condition )

dds<-estimateSizeFactors(dds)
# sizeFactors(dds)<-c(estimateSizeFactorsForMatrix(cnt[,1:6]), estimateSizeFactorsForMatrix(cnt[,7:12]) )


dds$Condition = relevel(dds$Condition,"C_293") #everything relative to Mock
dds<-DESeq(dds)

resultsNames(dds)
res <- results(dds, contrast=c("Condition","M_293","C_293")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.mem.cyt.293<-data.frame(res[,c(1,2,6)])
colnames(lfc.mem.cyt.293)<-paste0(colnames(lfc.mem.cyt.293),".mem.cyt.293")

res <- results(dds, contrast=c("Condition","M_2A15","C_2A15")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.mem.cyt.2A15<-data.frame(res[,c(1,2,6)])
colnames(lfc.mem.cyt.2A15)<-paste0(colnames(lfc.mem.cyt.2A15),".mem.cyt.2A15")

res <- results(dds, contrast=c("Condition","M_3C2","C_3C2")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.mem.cyt.3C2<-data.frame(res[,c(1,2,6)])
colnames(lfc.mem.cyt.3C2)<-paste0(colnames(lfc.mem.cyt.3C2),".mem.cyt.3C2")

dese<-cbind(lfc.mem.cyt.293, lfc.mem.cyt.2A15, lfc.mem.cyt.3C2)

fin<-merge(fin, dese, by.x="Row.names", by.y="row.names")

rel<-subset(fin, tpm_cutoff>=10)

ggplot()+geom_point(aes(rel$log2FoldChange.mem.cyt.2A15, rel$log2FoldChange.mem.cyt.293 ), colour="grey")+
  geom_point(aes(rel[rel$target=="target","log2FoldChange.mem.cyt.2A15"], rel[rel$target=="target","log2FoldChange.mem.cyt.293"] ),colour="orange")+
  geom_text(label=row.names(rel))+geom_abline(slope=1)

ggplot()+geom_point(aes(rel$log2FoldChange.mem.cyt.3C2, rel$log2FoldChange.mem.cyt.293 ), colour="grey")+
  geom_point(aes(rel[rel$target=="target","log2FoldChange.mem.cyt.2A15"], rel[rel$target=="target","log2FoldChange.mem.cyt.293"] ), colour="orange")+
  geom_abline(slope=1)

subset(rel, log2FoldChange.mem.cyt.2A15-log2FoldChange.mem.cyt.293<(-0.3) & target=="target" & tpm_cutoff>10, select=c("Symbol", "conv_CDS","reads_CDS"))

ggplot()+ geom_point(aes(rel[rel$target=="target","log2FoldChange.mem.cyt.2A15"], rel[rel$target=="target","conv_CDS"] ), colour="orange")+
  geom_abline(slope=1)


cnt<-tab[,colnames(tab)[grepl("C_", colnames(tab))|grepl("N_", colnames(tab))]]
coldata<-data.frame(SampleID=colnames(cnt),
                    Condition=sub(";","_",gsub("_.*","",sub("_",";",colnames(cnt)))),
                    Fraction=gsub("_.*","",colnames(cnt)),
                    Batch=gsub(".*_","",colnames(cnt)))
coldata<-apply(coldata,2,as.factor)

dds<-DESeqDataSetFromMatrix(countData = round(cnt),
                            colData = coldata, 
                            design =~ Condition )
dds<-estimateSizeFactors(dds)


dds$Condition = relevel(dds$Condition,"C_293")
dds<-DESeq(dds)

resultsNames(dds)
res <- results(dds, contrast=c("Condition","N_293","C_293")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.nuc.cyt.293<-data.frame(res[,c(1,2,6)])
colnames(lfc.nuc.cyt.293)<-paste0(colnames(lfc.nuc.cyt.293),".nuc.cyt.293")

res <- results(dds, contrast=c("Condition","N_2A15","C_2A15")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.nuc.cyt.2A15<-data.frame(res[,c(1,2,6)])
colnames(lfc.nuc.cyt.2A15)<-paste0(colnames(lfc.nuc.cyt.2A15),".nuc.cyt.2A15")

res <- results(dds, contrast=c("Condition","N_3C2","C_3C2")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.nuc.cyt.3C2<-data.frame(res[,c(1,2,6)])
colnames(lfc.nuc.cyt.3C2)<-paste0(colnames(lfc.nuc.cyt.3C2),".nuc.cyt.3C2")

dese<-cbind(lfc.nuc.cyt.293, lfc.nuc.cyt.2A15, lfc.nuc.cyt.3C2)

fin<-merge(fin, dese, by.x="Row.names", by.y="row.names")



cnt<-tab[,colnames(tab)[grepl("M_", colnames(tab))|grepl("T_", colnames(tab))]]
coldata<-data.frame(SampleID=colnames(cnt),
                    Condition=sub(";","_",gsub("_.*","",sub("_",";",colnames(cnt)))),
                    Fraction=gsub("_.*","",colnames(cnt)),
                    Batch=gsub(".*_","",colnames(cnt)))
coldata<-apply(coldata,2,as.factor)

dds<-DESeqDataSetFromMatrix(countData = round(cnt),
                            colData = coldata, 
                            design =~ Condition )
dds<-estimateSizeFactors(dds)


dds$Condition = relevel(dds$Condition,"T_293")
dds<-DESeq(dds)

resultsNames(dds)
res <- results(dds, contrast=c("Condition","M_293","T_293")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.mem.tot.293<-data.frame(res[,c(1,2,6)])
colnames(lfc.mem.tot.293)<-paste0(colnames(lfc.mem.tot.293),".mem.tot.293")

res <- results(dds, contrast=c("Condition","M_2A15","T_2A15")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.mem.tot.2A15<-data.frame(res[,c(1,2,6)])
colnames(lfc.mem.tot.2A15)<-paste0(colnames(lfc.mem.tot.2A15),".mem.tot.2A15")

res <- results(dds, contrast=c("Condition","M_3C2","T_3C2")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.mem.tot.3C2<-data.frame(res[,c(1,2,6)])
colnames(lfc.mem.tot.3C2)<-paste0(colnames(lfc.mem.tot.3C2),".mem.tot.3C2")

dese<-cbind(lfc.mem.tot.293, lfc.mem.tot.2A15, lfc.mem.tot.3C2)

fin<-merge(fin, dese, by.x="Row.names", by.y="row.names")


cnt<-tab[,colnames(tab)[grepl("N_", colnames(tab))|grepl("T_", colnames(tab))]]
coldata<-data.frame(SampleID=colnames(cnt),
                    Condition=sub(";","_",gsub("_.*","",sub("_",";",colnames(cnt)))),
                    Fraction=gsub("_.*","",colnames(cnt)),
                    Batch=gsub(".*_","",colnames(cnt)))
coldata<-apply(coldata,2,as.factor)

dds<-DESeqDataSetFromMatrix(countData = round(cnt),
                            colData = coldata, 
                            design =~ Condition )
dds<-estimateSizeFactors(dds)


dds$Condition = relevel(dds$Condition,"T_293")
dds<-DESeq(dds)

resultsNames(dds)
res <- results(dds, contrast=c("Condition","N_293","T_293")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.nuc.tot.293<-data.frame(res[,c(1,2,6)])
colnames(lfc.nuc.tot.293)<-paste0(colnames(lfc.nuc.tot.293),".nuc.tot.293")

res <- results(dds, contrast=c("Condition","N_2A15","T_2A15")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.nuc.tot.2A15<-data.frame(res[,c(1,2,6)])
colnames(lfc.nuc.tot.2A15)<-paste0(colnames(lfc.nuc.tot.2A15),".nuc.tot.2A15")

res <- results(dds, contrast=c("Condition","N_3C2","T_3C2")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.nuc.tot.3C2<-data.frame(res[,c(1,2,6)])
colnames(lfc.nuc.tot.3C2)<-paste0(colnames(lfc.nuc.tot.3C2),".nuc.tot.3C2")

dese<-cbind(lfc.nuc.tot.293, lfc.nuc.tot.2A15, lfc.nuc.tot.3C2)

fin<-merge(fin, dese, by.x="Row.names", by.y="row.names")



cnt<-tab[,colnames(tab)[grepl("C_", colnames(tab))|grepl("T_", colnames(tab))]]
coldata<-data.frame(SampleID=colnames(cnt),
                    Condition=sub(";","_",gsub("_.*","",sub("_",";",colnames(cnt)))),
                    Fraction=gsub("_.*","",colnames(cnt)),
                    Batch=gsub(".*_","",colnames(cnt)))
coldata<-apply(coldata,2,as.factor)

dds<-DESeqDataSetFromMatrix(countData = round(cnt),
                            colData = coldata, 
                            design =~ Condition )
dds<-estimateSizeFactors(dds)


dds$Condition = relevel(dds$Condition,"T_293")
dds<-DESeq(dds)

resultsNames(dds)
res <- results(dds, contrast=c("Condition","C_293","T_293")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.cyt.tot.293<-data.frame(res[,c(1,2,6)])
colnames(lfc.cyt.tot.293)<-paste0(colnames(lfc.cyt.tot.293),".cyt.tot.293")

res <- results(dds, contrast=c("Condition","C_2A15","T_2A15")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.cyt.tot.2A15<-data.frame(res[,c(1,2,6)])
colnames(lfc.cyt.tot.2A15)<-paste0(colnames(lfc.cyt.tot.2A15),".cyt.tot.2A15")

res <- results(dds, contrast=c("Condition","C_3C2","T_3C2")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.cyt.tot.3C2<-data.frame(res[,c(1,2,6)])
colnames(lfc.cyt.tot.3C2)<-paste0(colnames(lfc.cyt.tot.3C2),".cyt.tot.3C2")

dese<-cbind(lfc.cyt.tot.293, lfc.cyt.tot.2A15, lfc.cyt.tot.3C2)

fin<-merge(fin, dese, by.x="Row.names", by.y="row.names")




###multifactor design for differences between KO and WT
cnt<-tab[,colnames(tab)[grepl("M_", colnames(tab))|grepl("T_", colnames(tab))|grepl("C_", colnames(tab))|grepl("N_", colnames(tab))]]

coldata<-data.frame(SampleID=colnames(cnt),
                    Condition=gsub("_.*","",gsub(".*;","",sub("_",";",colnames(cnt)))),
                    Fraction=gsub("_.*","",colnames(cnt)),
                    Batch=gsub(".*_","",colnames(cnt)))
coldata<-apply(coldata,2,as.factor)



dds<-DESeqDataSetFromMatrix(countData = round(cnt),
                            colData = coldata, 
                            design =~ Fraction + Condition +Fraction:Condition)
dds<-estimateSizeFactors(dds)

dds$Condition = relevel(dds$Condition,"293") #everything relative to 293
dds$Fraction = relevel(dds$Fraction,"C") #everything relative to C
dds<-DESeq(dds, test="LRT", reduced=~Fraction+Condition)

resultsNames(dds)

res <- results(dds, contrast=list("FractionM.Condition2A15")) 
summary(res)
DESeq2::plotMA(res)
lfc.mem.cyt.2A15.293<-data.frame(res[,c(1,2,6)])
colnames(lfc.mem.cyt.2A15.293)<-paste0(colnames(lfc.mem.cyt.2A15.293),".mem.cyt.2A15.293")

res <- results(dds, contrast=list("FractionM.Condition3C2")) 
summary(res)
DESeq2::plotMA(res)
lfc.mem.cyt.3C2.293<-data.frame(res[,c(1,2,6)])
colnames(lfc.mem.cyt.3C2.293)<-paste0(colnames(lfc.mem.cyt.3C2.293),".mem.cyt.3C2.293")

res <- results(dds, contrast=list("FractionN.Condition2A15")) 
summary(res)
DESeq2::plotMA(res)
lfc.nuc.cyt.2A15.293<-data.frame(res[,c(1,2,6)])
colnames(lfc.nuc.cyt.2A15.293)<-paste0(colnames(lfc.nuc.cyt.2A15.293),".nuc.cyt.2A15.293")

res <- results(dds, contrast=list("FractionN.Condition3C2")) 
summary(res)
DESeq2::plotMA(res)
lfc.nuc.cyt.3C2.293<-data.frame(res[,c(1,2,6)])
colnames(lfc.nuc.cyt.3C2.293)<-paste0(colnames(lfc.nuc.cyt.3C2.293),".nuc.cyt.3C2.293")





dds$Condition = relevel(dds$Condition,"293") #everything relative to 293
dds$Fraction = relevel(dds$Fraction,"T") #everything relative to T
dds<-DESeq(dds, test="LRT", reduced=~Fraction+Condition)

resultsNames(dds)

res <- results(dds, contrast=list("FractionM.Condition2A15")) 
summary(res)
DESeq2::plotMA(res)
lfc.mem.tot.2A15.293<-data.frame(res[,c(1,2,6)])
colnames(lfc.mem.tot.2A15.293)<-paste0(colnames(lfc.mem.tot.2A15.293),".mem.tot.2A15.293")

res <- results(dds, contrast=list("FractionM.Condition3C2")) 
summary(res)
DESeq2::plotMA(res)
lfc.mem.tot.3C2.293<-data.frame(res[,c(1,2,6)])
colnames(lfc.mem.tot.3C2.293)<-paste0(colnames(lfc.mem.tot.3C2.293),".mem.tot.3C2.293")

res <- results(dds, contrast=list("FractionN.Condition2A15")) 
summary(res)
DESeq2::plotMA(res)
lfc.nuc.tot.2A15.293<-data.frame(res[,c(1,2,6)])
colnames(lfc.nuc.tot.2A15.293)<-paste0(colnames(lfc.nuc.tot.2A15.293),".nuc.tot.2A15.293")

res <- results(dds, contrast=list("FractionN.Condition3C2")) 
summary(res)
DESeq2::plotMA(res)
lfc.nuc.tot.3C2.293<-data.frame(res[,c(1,2,6)])
colnames(lfc.nuc.tot.3C2.293)<-paste0(colnames(lfc.nuc.tot.3C2.293),".nuc.tot.3C2.293")

res <- results(dds, contrast=list("FractionC.Condition2A15")) 
summary(res)
DESeq2::plotMA(res)
lfc.cyt.tot.2A15.293<-data.frame(res[,c(1,2,6)])
colnames(lfc.cyt.tot.2A15.293)<-paste0(colnames(lfc.cyt.tot.2A15.293),".cyt.tot.2A15.293")

res <- results(dds, contrast=list("FractionC.Condition3C2")) 
summary(res)
DESeq2::plotMA(res)
lfc.cyt.tot.3C2.293<-data.frame(res[,c(1,2,6)])
colnames(lfc.cyt.tot.3C2.293)<-paste0(colnames(lfc.cyt.tot.3C2.293),".cyt.tot.3C2.293")


dese<-cbind(lfc.mem.cyt.2A15.293, lfc.mem.cyt.3C2.293, 
            lfc.nuc.cyt.2A15.293, lfc.nuc.cyt.3C2.293, 
            lfc.mem.tot.2A15.293, lfc.mem.tot.3C2.293, 
            lfc.nuc.tot.2A15.293, lfc.nuc.tot.3C2.293, 
            lfc.cyt.tot.2A15.293, lfc.cyt.tot.3C2.293)

fin<-merge(fin, dese, by.x="Row.names", by.y="row.names")




###merge KO
cnt<-tab[,colnames(tab)[grepl("C_", colnames(tab))|grepl("M_", colnames(tab))]]


coldata<-data.frame(SampleID=colnames(cnt),
                    Condition=ifelse(grepl("293", colnames(cnt)),paste0(gsub("_.*","",colnames(cnt)),"_WT"), paste0(gsub("_.*","",colnames(cnt)),"_KO")),
                    Fraction=gsub("_.*","",colnames(cnt)),
                    Batch=gsub(".*_","",colnames(cnt)))
coldata<-apply(coldata,2,as.factor)

dds<-DESeqDataSetFromMatrix(countData = round(cnt),
                            colData = coldata, 
                            design =~ Condition )

dds<-estimateSizeFactors(dds)
# sizeFactors(dds)<-c(estimateSizeFactorsForMatrix(cnt[,1:6]), estimateSizeFactorsForMatrix(cnt[,7:12]) )


dds$Condition = relevel(dds$Condition,"C_WT") 
dds<-DESeq(dds)

resultsNames(dds)
res <- results(dds, contrast=c("Condition","M_WT","C_WT")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.mem.cyt.WT<-data.frame(res[,c(1,2,6)])
colnames(lfc.mem.cyt.WT)<-paste0(colnames(lfc.mem.cyt.WT),".mem.cyt.WT")

resultsNames(dds)
res <- results(dds, contrast=c("Condition","M_KO","C_KO")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.mem.cyt.KO<-data.frame(res[,c(1,2,6)])
colnames(lfc.mem.cyt.KO)<-paste0(colnames(lfc.mem.cyt.KO),".mem.cyt.KO")

dese<-cbind(lfc.mem.cyt.WT, lfc.mem.cyt.KO)

fin<-merge(fin, dese, by.x="Row.names", by.y="row.names")



###merge KO and multifactor

cnt<-tab[,colnames(tab)[grepl("M_", colnames(tab))|grepl("T_", colnames(tab))|grepl("C_", colnames(tab))|grepl("N_", colnames(tab))]]

coldata<-data.frame(SampleID=colnames(cnt),
                    Condition=ifelse(grepl("293", colnames(cnt)),"WT", "KO"),
                    Fraction=gsub("_.*","",colnames(cnt)),
                    Batch=gsub(".*_","",colnames(cnt)))
coldata<-apply(coldata,2,as.factor)



dds<-DESeqDataSetFromMatrix(countData = round(cnt),
                            colData = coldata, 
                            design =~ Fraction + Condition +Fraction:Condition)
dds<-estimateSizeFactors(dds)

dds$Condition = relevel(dds$Condition,"WT") #everything relative to 293
dds$Fraction = relevel(dds$Fraction,"C") #everything relative to C
dds<-DESeq(dds, test="LRT", reduced=~Fraction+Condition)

resultsNames(dds)

res <- results(dds, contrast=list("FractionM.ConditionKO")) 
summary(res)
DESeq2::plotMA(res)
lfc.mem.cyt.KO.293<-data.frame(res[,c(1,2,6)])
colnames(lfc.mem.cyt.KO.293)<-paste0(colnames(lfc.mem.cyt.KO.293),".mem.cyt.KO.293")

res <- results(dds, contrast=list("FractionN.ConditionKO")) 
summary(res)
DESeq2::plotMA(res)
lfc.nuc.cyt.KO.293<-data.frame(res[,c(1,2,6)])
colnames(lfc.nuc.cyt.KO.293)<-paste0(colnames(lfc.nuc.cyt.KO.293),".nuc.cyt.KO.293")

res <- results(dds, contrast=list("FractionT.ConditionKO")) 
summary(res)
DESeq2::plotMA(res)
lfc.tot.cyt.KO.293<-data.frame(res[,c(1,2,6)])
colnames(lfc.tot.cyt.KO.293)<-paste0(colnames(lfc.tot.cyt.KO.293),".tot.cyt.KO.293")

dese<-cbind(lfc.mem.cyt.KO.293, lfc.nuc.cyt.KO.293, lfc.tot.cyt.KO.293)

fin<-merge(fin, dese, by.x="Row.names", by.y="row.names")

rel<-subset(fin, tpm_cutoff>=10)

#major plots

ggplot()+geom_point(aes(rel$log2FoldChange.cyt.tot.2A15, rel$log2FoldChange.cyt.tot.293 ), colour="grey")+
  geom_point(aes(rel[rel$target=="target","log2FoldChange.cyt.tot.2A15"], rel[rel$target=="target","log2FoldChange.cyt.tot.293"] ),colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(rel$log2FoldChange.cyt.tot.3C2, rel$log2FoldChange.cyt.tot.293 ), colour="grey")+
  geom_point(aes(rel[rel$target=="target","log2FoldChange.cyt.tot.2A15"], rel[rel$target=="target","log2FoldChange.cyt.tot.293"] ), colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(rel$log2FoldChange.nuc.tot.2A15, rel$log2FoldChange.nuc.tot.293 ), colour="grey")+
  geom_point(aes(rel[rel$target=="target","log2FoldChange.nuc.tot.2A15"], rel[rel$target=="target","log2FoldChange.nuc.tot.293"] ),colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(rel$log2FoldChange.nuc.tot.3C2, rel$log2FoldChange.nuc.tot.293 ), colour="grey")+
  geom_point(aes(rel[rel$target=="target","log2FoldChange.nuc.tot.2A15"], rel[rel$target=="target","log2FoldChange.nuc.tot.293"] ), colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(rel$log2FoldChange.mem.tot.2A15, rel$log2FoldChange.mem.tot.293 ), colour="grey")+
  geom_point(aes(rel[rel$target=="target","log2FoldChange.mem.tot.2A15"], rel[rel$target=="target","log2FoldChange.mem.tot.293"] ),colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(rel$log2FoldChange.mem.tot.3C2, rel$log2FoldChange.mem.tot.293 ), colour="grey")+
  geom_point(aes(rel[rel$target=="target","log2FoldChange.mem.tot.2A15"], rel[rel$target=="target","log2FoldChange.mem.tot.293"] ), colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(rel$log2FoldChange.mem.cyt.2A15, rel$log2FoldChange.mem.cyt.293 ), colour="grey")+
  geom_point(aes(rel[rel$target=="target","log2FoldChange.mem.cyt.2A15"], rel[rel$target=="target","log2FoldChange.mem.cyt.293"] ),colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(rel$log2FoldChange.mem.cyt.3C2, rel$log2FoldChange.mem.cyt.293 ), colour="grey")+
  geom_point(aes(rel[rel$target=="target","log2FoldChange.mem.cyt.2A15"], rel[rel$target=="target","log2FoldChange.mem.cyt.293"] ), colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(rel$log2FoldChange.nuc.cyt.2A15, rel$log2FoldChange.nuc.cyt.293 ), colour="grey")+
  geom_point(aes(rel[rel$target=="target","log2FoldChange.nuc.cyt.2A15"], rel[rel$target=="target","log2FoldChange.nuc.cyt.293"] ),colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(rel$log2FoldChange.nuc.cyt.3C2, rel$log2FoldChange.nuc.cyt.293 ), colour="grey")+
  geom_point(aes(rel[rel$target=="target","log2FoldChange.nuc.cyt.2A15"], rel[rel$target=="target","log2FoldChange.nuc.cyt.293"] ), colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(rel$log2FoldChange.mem.cyt.KO, rel$log2FoldChange.mem.cyt.WT ), colour="grey")+
  geom_point(aes(rel[rel$target=="target","log2FoldChange.mem.cyt.KO"], rel[rel$target=="target","log2FoldChange.mem.cyt.WT"] ),colour="orange")+
  #geom_text(aes(rel[rel$target=="target","log2FoldChange.mem.cyt.KO"], rel[rel$target=="target","log2FoldChange.mem.cyt.WT"], label=rel[rel$target=="target","Symbol"]), size=3)+
  geom_abline(slope=1)+ylim(-2,10)

subset(rel, log2FoldChange.mem.cyt.2A15-log2FoldChange.mem.cyt.293<(-0.3) & target=="target" & tpm_cutoff>10, select=c("Symbol", "conv_CDS","reads_CDS"))

ggplot()+ geom_point(aes(rel[rel$target=="target","log2FoldChange.mem.cyt.2A15"], rel[rel$target=="target","conv_CDS"] ), colour="orange")+
  geom_abline(slope=1)

rel<-subset(fin, tpm_cutoff>=10)
ggplot()+geom_point(aes(rel$log2FoldChange.nuc.cyt.2A15.293, rel$log2FoldChange.nuc.cyt.3C2.293 ), colour="grey")+
  geom_point(aes(rel[rel$target=="target","log2FoldChange.nuc.cyt.2A15.293"], rel[rel$target=="target","log2FoldChange.nuc.cyt.3C2.293"] ), colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(log2(rel$baseMean.mem.cyt.2A15.293), rel$log2FoldChange.mem.cyt.2A15.293 ), colour="grey")+
  geom_point(aes(log2(rel[rel$target=="target","baseMean.mem.cyt.2A15.293"]), rel[rel$target=="target","log2FoldChange.mem.cyt.2A15.293"] ), colour="orange")+
  geom_text(aes(log2(rel[rel$target=="target","baseMean.mem.cyt.2A15.293"]), rel[rel$target=="target","log2FoldChange.mem.cyt.2A15.293"], label=rel[rel$target=="target","Symbol"]), size=3)+
  geom_abline(slope=0)+ylim(-6,6)

ggplot()+geom_point(aes(log2(rel$baseMean.mem.cyt.3C2.293), rel$log2FoldChange.mem.cyt.3C2.293 ), colour="grey")+
  geom_point(aes(log2(rel[rel$target=="target","baseMean.mem.cyt.3C2.293"]), rel[rel$target=="target","log2FoldChange.mem.cyt.3C2.293"] ), colour="orange")+
  geom_abline(slope=0)+ylim(-6,6)

ggplot()+geom_point(aes(log2(rel$baseMean.mem.cyt.KO.293), rel$log2FoldChange.mem.cyt.KO.293 ), colour="grey")+
  geom_point(aes(log2(rel[rel$target=="target","baseMean.mem.cyt.KO.293"]), rel[rel$target=="target","log2FoldChange.mem.cyt.KO.293"] ), colour="orange")+
  #geom_text(aes(log2(rel[rel$target=="target","baseMean.mem.cyt.KO.293"]), rel[rel$target=="target","log2FoldChange.mem.cyt.KO.293"], label=rel[rel$target=="target","Symbol"]), size=3)+
  geom_abline(slope=0)+ylim(-6,6)





fin$changeMem_sig.2A15<-ifelse(fin$log2FoldChange.mem.cyt.2A15.293<(-0.1) & fin$padj.mem.cyt.2A15.293<0.1 & 
                                 fin$log2FoldChange.mem.cyt.2A15>2 & fin$log2FoldChange.mem.cyt.WT>2 &
                                       fin$tpm_cutoff>=10 ,
                                     "downMem_sig.2A15", NA)

nrow(subset(fin, changeMem_sig.2A15=="downMem_sig.2A15"))
nrow(subset(fin, changeMem_sig.2A15=="downMem_sig.2A15"&target=="target"))


fin$changeMem_sig.3C2<-ifelse(fin$log2FoldChange.mem.cyt.3C2.293<(-0.1) & fin$padj.mem.cyt.3C2.293<0.1 & 
                                fin$log2FoldChange.mem.cyt.3C2>2 & fin$log2FoldChange.mem.cyt.WT>2 &
                                        fin$tpm_cutoff>=10 ,
                                      "downMem_sig.3C2", NA)
nrow(subset(fin, changeMem_sig.3C2=="downMem_sig.3C2"))
nrow(subset(fin, changeMem_sig.3C2=="downMem_sig.3C2" & target=="target"))

nrow(subset(fin, changeMem_sig.2A15=="downMem_sig.2A15" & changeMem_sig.3C2=="downMem_sig.3C2" ))



fin$changeMem_enrichment.2A15<-ifelse(fin$log2FoldChange.mem.cyt.2A15.293<(-.1) & 
                                        fin$log2FoldChange.mem.cyt.2A15>2 & fin$log2FoldChange.mem.cyt.WT>2 &
                                 fin$tpm_cutoff>=10 ,
                               "downMem_enrichment.2A15", NA)
nrow(subset(fin, changeMem_enrichment.2A15=="downMem_enrichment.2A15"))
nrow(subset(fin, changeMem_enrichment.2A15=="downMem_enrichment.2A15" & target=="target"))


fin$changeMem_enrichment.3C2<-ifelse(fin$log2FoldChange.mem.cyt.3C2.293<(-.1) & 
                                       fin$log2FoldChange.mem.cyt.3C2>2 & fin$log2FoldChange.mem.cyt.WT>2 &
                                fin$tpm_cutoff>=10 ,
                              "downMem_enrichment.3C2", NA)
nrow(subset(fin, changeMem_enrichment.3C2=="downMem_enrichment.3C2"))
nrow(subset(fin, changeMem_enrichment.3C2=="downMem_enrichment.3C2" & target=="target"))

nrow(subset(fin, changeMem_enrichment.2A15=="downMem_enrichment.2A15" & changeMem_enrichment.3C2=="downMem_enrichment.3C2" ))


fin$changeMem_sig.KO<-ifelse(fin$log2FoldChange.mem.cyt.KO.293<(-0.1) & fin$padj.mem.cyt.KO.293<0.1 & 
                               fin$log2FoldChange.mem.cyt.KO>2 & fin$log2FoldChange.mem.cyt.WT>2 &
                                 fin$tpm_cutoff>=10 ,
                               "downMem_sig.KO", NA)
nrow(subset(fin, changeMem_sig.KO=="downMem_sig.KO"))
nrow(subset(fin, changeMem_sig.KO=="downMem_sig.KO" & target=="target"))

fin$changeMem_enrichment.KO<-ifelse(fin$log2FoldChange.mem.cyt.KO.293<(-0.1) &  
                                      fin$log2FoldChange.mem.cyt.KO>2 & fin$log2FoldChange.mem.cyt.WT>2 &
                               fin$tpm_cutoff>=10 ,
                             "downMem_enrichment.KO", NA)

fin$changeMem_enrichment.KO<-ifelse(fin$log2FoldChange.mem.cyt.KO.293<(-0.1) &  
                                      fin$log2FoldChange.mem.cyt.KO>2 & fin$log2FoldChange.mem.cyt.WT>2 &
                                      fin$tpm_cutoff>=10 ,
                                    "downMem_enrichment.KO", 
                                    ifelse(fin$log2FoldChange.mem.cyt.KO>2 & fin$log2FoldChange.mem.cyt.WT>2 &
                                             fin$tpm_cutoff>=10, "all_membrane_localized",NA))

nrow(subset(fin, changeMem_enrichment.KO=="downMem_enrichment.KO"))
nrow(subset(fin, changeMem_enrichment.KO=="downMem_enrichment.KO"& target=="target"))
subset(fin, changeMem_enrichment.KO=="downMem_enrichment.KO"& target=="target", 
       select=c("Symbol","conv_CDS","log2FoldChange.mem.cyt.KO.293"))
subset(fin, changeMem_enrichment.KO=="downMem_enrichment.KO"& target=="target", 
       select=c("Symbol"))

nrow(subset(fin, changeMem_sig.KO=="downMem_sig.KO" & changeMem_enrichment.KO=="downMem_enrichment.KO"))

rel<-subset(fin, tpm_cutoff>=10)

ggplot()+geom_point(aes(rel$log2FoldChange.mem.cyt.2A15, rel$log2FoldChange.mem.cyt.293 ), colour="grey")+
  geom_point(aes(rel[rel$changeMem_sig.2A15=="downMem_sig.2A15","log2FoldChange.mem.cyt.2A15"], rel[rel$changeMem_sig.2A15=="downMem_sig.2A15","log2FoldChange.mem.cyt.293"] ),colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(log2(rel$baseMean.mem.KO.WT), rel$log2FoldChange.mem.KO.WT ), colour="grey")+
  geom_point(aes(log2(rel[rel$target=="target","baseMean.mem.KO.WT"]), rel[rel$target=="target","log2FoldChange.mem.KO.WT"] ),colour="orange")+
  geom_abline(slope=0)

ggplot()+geom_point(aes(log2(rel$baseMean.tot.KO.WT), rel$log2FoldChange.tot.KO.WT ), colour="grey")+
  geom_point(aes(log2(rel[rel$target=="target","baseMean.tot.KO.WT"]), rel[rel$target=="target","log2FoldChange.tot.KO.WT"] ),colour="orange")+
  geom_abline(slope=0)


ggplot()+geom_point(aes(rel$log2FoldChange.mem.cyt.KO, rel$log2FoldChange.mem.cyt.WT ), colour="grey")+
  geom_point(aes(rel[rel$changeMem_sig.KO=="downMem_sig.KO","log2FoldChange.mem.cyt.KO"], rel[rel$changeMem_sig.KO=="downMem_sig.KO","log2FoldChange.mem.cyt.WT"] ),colour="orange")+
  geom_abline(slope=1)


ggplot()+geom_point(aes(rel$log2FoldChange.mem.cyt.2A15, rel$log2FoldChange.mem.cyt.293 ), colour="grey")+
  geom_point(aes(rel[rel$changeMem_enrichment.2A15=="downMem_enrichment.2A15","log2FoldChange.mem.cyt.2A15"], rel[rel$changeMem_enrichment.2A15=="downMem_enrichment.2A15","log2FoldChange.mem.cyt.293"] ),colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(rel$log2FoldChange.mem.cyt.KO, rel$log2FoldChange.mem.cyt.WT ), colour="grey")+
  geom_point(aes(rel[rel$changeMem_enrichment.KO=="downMem_enrichment.KO","log2FoldChange.mem.cyt.KO"], rel[rel$changeMem_enrichment.KO=="downMem_enrichment.KO","log2FoldChange.mem.cyt.WT"] ),colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(rel$log2FoldChange.cyt.tot.2A15, rel$log2FoldChange.cyt.tot.293 ), colour="grey")+
  geom_point(aes(rel[rel$changeMem_enrichment.2A15=="downMem_enrichment.2A15","log2FoldChange.cyt.tot.2A15"], rel[rel$changeMem_enrichment.2A15=="downMem_enrichment.2A15","log2FoldChange.cyt.tot.293"] ),colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(rel$log2FoldChange.cyt.tot.2A15, rel$log2FoldChange.cyt.tot.293 ), colour="grey")+
  geom_point(aes(rel[rel$changeMem_enrichment.KO=="downMem_enrichment.KO","log2FoldChange.cyt.tot.2A15"], rel[rel$changeMem_enrichment.KO=="downMem_enrichment.KO","log2FoldChange.cyt.tot.293"] ),colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(rel$log2FoldChange.cyt.tot.3C2, rel$log2FoldChange.cyt.tot.293 ), colour="grey")+
  geom_point(aes(rel[rel$changeMem_enrichment.KO=="downMem_enrichment.KO","log2FoldChange.cyt.tot.3C2"], rel[rel$changeMem_enrichment.KO=="downMem_enrichment.KO","log2FoldChange.cyt.tot.293"] ),colour="orange")+
  geom_abline(slope=1)


ggplot()+geom_point(aes(rel$log2FoldChange.nuc.cyt.2A15, rel$log2FoldChange.nuc.cyt.293 ), colour="grey")+
  geom_point(aes(rel[rel$changeMem_enrichment.KO=="downMem_enrichment.KO","log2FoldChange.nuc.cyt.2A15"], rel[rel$changeMem_enrichment.KO=="downMem_enrichment.KO","log2FoldChange.nuc.cyt.293"] ),colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(rel$log2FoldChange.nuc.cyt.3C2, rel$log2FoldChange.nuc.cyt.293 ), colour="grey")+
  geom_point(aes(rel[rel$changeMem_enrichment.KO=="downMem_enrichment.KO","log2FoldChange.nuc.cyt.3C2"], rel[rel$changeMem_enrichment.KO=="downMem_enrichment.KO","log2FoldChange.nuc.cyt.293"] ),colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(rel$lfc.tot.2A15_1, rel$lfc.tot.3C2_1 ), colour="grey")+
  geom_point(aes(rel[rel$changeMem_enrichment.2A15=="downMem_enrichment.2A15","lfc.tot.2A15_1"], rel[rel$changeMem_enrichment.2A15=="downMem_enrichment.2A15","lfc.tot.3C2_1"] ),colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(rel$log2FoldChange.mem.tot.3C2, rel$log2FoldChange.mem.tot.293 ), colour="grey")+
  geom_point(aes(rel[rel$changeMem_enrichment.3C2=="downMem_enrichment.3C2","log2FoldChange.mem.tot.3C2"], rel[rel$changeMem_enrichment.3C2=="downMem_enrichment.3C2","log2FoldChange.mem.tot.293"] ),colour="orange")+
  geom_abline(slope=1)



ggplot(rel, aes(log2FoldChange.mem.cyt.KO.293, colour=changeMem_enrichment.KO))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))
ggplot(rel, aes(log2FoldChange.nuc.cyt.KO.293, colour=changeMem_enrichment.KO))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))
ggplot(rel, aes(log2FoldChange.tot.cyt.KO.293, colour=changeMem_enrichment.KO))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))

ggplot(rel, aes(log2FoldChange.cyt.tot.2A15.293, colour=changeMem_enrichment.KO))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(rel, aes(log2FoldChange.cyt.tot.3C2.293, colour=changeMem_enrichment.KO))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(rel, aes(log2FoldChange.tot.cyt.KO.293, colour=changeMem_enrichment.KO))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))

ggplot(rel, aes(log2FoldChange.mem.KO.WT, colour=changeMem_enrichment.KO))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(rel, aes(log2FoldChange.cyt.KO.WT, colour=changeMem_enrichment.KO))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(rel, aes(log2FoldChange.nuc.KO.WT, colour=changeMem_enrichment.KO))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(rel, aes(log2FoldChange.tot.KO.WT, colour=changeMem_enrichment.KO))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))

ggplot(rel, aes(log2FoldChange.tot.KO.WT, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(rel, aes(log2FoldChange.mem.KO.WT, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(rel, aes(log2FoldChange.cyt.KO.WT, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(rel, aes(log2FoldChange.nuc.KO.WT, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))


ggplot(rel, aes(log2FoldChange.nuc.cyt.KO.293, colour=changeMem_enrichment.KO))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))
ggplot(rel, aes(log2FoldChange.mem.tot.2A15.293, colour=changeMem_enrichment.KO))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))


ggplot(rel, aes(lfc.cyt.2A15_2, colour=changeMem_enrichment.KO))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))
ggplot(rel, aes(lfc.cyt.3C2_1, colour=changeMem_enrichment.KO))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))
ggplot(rel, aes(lfc.cyt.3C2_2, colour=changeMem_enrichment.KO))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))


subset(rel, log2FoldChange.mem.cyt.2A15-log2FoldChange.mem.cyt.293<(-0.1) & 
         log2FoldChange.mem.cyt.3C2-log2FoldChange.mem.cyt.293<(-0.1) & 
         target=="target" & tpm_cutoff>10, select=c("Symbol", "conv_CDS","reads_CDS","hits_CDS","log2FoldChange.mem.cyt.2A15.293"))

subset(rel, log2FoldChange.mem.cyt.2A15.293<0 & padj.mem.cyt.2A15.293<0.1 & 
         target=="target" & tpm_cutoff>10, select=c("Symbol", "conv_CDS","reads_CDS","hits_CDS","log2FoldChange.mem.cyt.2A15.293"))
subset(fin, log2FoldChange.mem.cyt.2A15.293<(-0.1) & 
         target=="target" & tpm_cutoff>10, select=c("Symbol", "conv_CDS","reads_CDS","hits_CDS","log2FoldChange.mem.cyt.2A15.293"))


subset(rel, log2FoldChange.mem.cyt.3C2.293<0 & padj.mem.cyt.3C2.293<0.1 & 
         target=="target" & tpm_cutoff>10, select=c("Symbol", "conv_CDS","reads_CDS","hits_CDS","log2FoldChange.mem.cyt.3C2.293"))

subset(rel,padj.mem.cyt.KO.293<0.1 & log2FoldChange.mem.cyt.KO.293<0 , select=c("Symbol","log2FoldChange.mem.cyt.KO.293", "log2FoldChange.mem.cyt.293"))





###get per replicate fold changes
cnt<-tab[,colnames(tab)[grepl("C_", colnames(tab))|grepl("M_", colnames(tab))]]


coldata<-data.frame(SampleID=colnames(cnt),
                    Condition=sub(";","_",gsub("_.*","",sub("_",";",colnames(cnt)))),
                    Fraction=gsub("_.*","",colnames(cnt)),
                    Batch=gsub(".*_","",colnames(cnt)))
coldata<-apply(coldata,2,as.factor)

dds<-DESeqDataSetFromMatrix(countData = round(cnt),
                            colData = coldata, 
                            design =~ Condition )

dds<-estimateSizeFactors(dds)
# sizeFactors(dds)<-c(estimateSizeFactorsForMatrix(cnt[,1:6]), estimateSizeFactorsForMatrix(cnt[,7:12]) )


dds$Condition = relevel(dds$Condition,"C_293") #everything relative to Mock
dds<-DESeq(dds)

resultsNames(dds)
res <- results(dds, contrast=c("Condition","M_293","C_293")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
coun<-data.frame(counts(dds, normalize=T))[,c(1,2,7,8)]
coun$log2FoldChange.mem.cyt.293_1<-ifelse(coun$C_293_1==0 | coun$M_293_1==0, NA, 
                                          log2(coun$M_293_1/coun$C_293_1))
coun$log2FoldChange.mem.cyt.293_2<-ifelse(coun$C_293_2==0 | coun$M_293_2==0, NA, 
                                          log2(coun$M_293_2/coun$C_293_2))
ggplot(coun, aes(log2FoldChange.mem.cyt.293_1, log2FoldChange.mem.cyt.293_2))+geom_point()
# write.table(coun[,5:6], "~/Google Drive/hdlbp/mem.cyt.293.lfc.replicates.txt", quote=F, sep="\t", row.names=T )

fin<-merge(fin, coun, by.x="Row.names", by.y="row.names")
