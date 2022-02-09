library(ggplot2)
library(reshape2)
library(corrplot)
library(DESeq2)
setwd("/Volumes/landthaler/pcp/projects/miha/HDLBP/frac/")
# setwd("E:/work/hdlbp/frac/")


files<-list.files(getwd(), pattern="genes")

dat<-lapply(files, read.delim)
names(dat)<-gsub("\\..*","",files)

tab<-data.frame(sapply(dat, function(x) cbind(x[,5])))
row.names(tab)<-dat[[1]][,1]

corrplot(cor(tab), type="upper")


tpm<-data.frame(sapply(dat, function(x) cbind(x[,6])))
row.names(tpm)<-dat[[1]][,1]

corrplot(cor(tpm), type="upper")


coldata<-data.frame(SampleID=colnames(tab),
                    Condition=sub(";","_",gsub("_.*","",sub("_",";",colnames(tab)))),
                    Batch=substr(colnames(tab), nchar(colnames(tab)), nchar(colnames(tab))))

coldata<-apply(coldata,2,as.factor)
coldata[7,2]<-"crispr"
coldata[8,2]<-"crispr"

dds<-DESeqDataSetFromMatrix(countData = round(tab),
                            colData = coldata, 
                            design = ~Batch+Condition)


vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Condition", "Batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Batch)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 



subs<-tab[,c(1:14,17:22,27:ncol(tab)) ]
coldata<-data.frame(SampleID=colnames(subs),
                    Condition=sub(";","_",gsub("_.*","",sub("_",";",colnames(subs)))),
                    Batch=substr(colnames(subs), nchar(colnames(subs)), nchar(colnames(subs))))

coldata<-apply(coldata,2,as.factor)
coldata[7,2]<-"T_crispr"
coldata[8,2]<-"T_crispr"


dds<-DESeqDataSetFromMatrix(countData = round(subs),
                            colData = coldata, 
                            design = ~Batch+Condition)


vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Condition", "Batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Batch)) +
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



si<-tab[,c(15:16,23:26) ]
coldata<-data.frame(SampleID=colnames(si),
                    Condition=c("Mock", "Mock", "si1","si1","si2","si2"),
                    Batch=substr(colnames(si), nchar(colnames(si)), nchar(colnames(si))))

coldata<-apply(coldata,2,as.factor)


dds<-DESeqDataSetFromMatrix(countData = round(si),
                            colData = coldata, 
                            design = ~Batch+Condition)


vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Condition", "Batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Batch)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 


genes<-read.table("/Volumes/landthaler/pcp/projects/miha/DDR/nascent/RSEM/geneInfo.txt", header=T)[,c(1,4:5,10,12)]
# genes<-read.table("geneInfo.txt", header=T)[,c(1,4:5,10,12)]

tar<-read.delim("neel_HDLBP2_targets.txt", header=T)
colnames(tar)[1]<-"gene_name"
tar$target<-"target"
sig<-read.csv("mart_export.txt", header=T)

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


mel<-melt(fin, measure.vars = colnames(fin[,3:34]), id.vars = c("Symbol", "Annotation","target","Transmembrane.helices", 
                                                                "Cleavage.site..Signalp.","Low.complexity..Seg."))

ggplot(mel, aes(variable,log2(value), fill=target))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(mel, aes(variable,log2(value), fill=Transmembrane.helices))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(mel, aes(variable,log2(value), fill=Cleavage.site..Signalp.))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(mel, aes(variable,log2(value), fill=Low.complexity..Seg.))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
ggplot(mel[mel$Annotation=="protein_coding",], aes(variable,log2(value), fill=target))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(mel[mel$Annotation=="protein_coding",], aes(variable,log2(value), fill=Transmembrane.helices))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(mel[mel$Annotation=="protein_coding",], aes(variable,log2(value), fill=Cleavage.site..Signalp.))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))

cyto<-fin[,colnames(fin)[grepl("C_", colnames(fin))]]
mem<-fin[,colnames(fin)[grepl("M_", colnames(fin))]]
nuc<-fin[,colnames(fin)[grepl("N_", colnames(fin))]]
tot<-fin[,colnames(fin)[grepl("T_", colnames(fin))]]
tot[tot==0]<-NA

lfc.cyt<-cyto/tot
diff.cyt<-lfc.cyt
diff.cyt[diff.cyt<1]<-NA
diff.cyt$diff.C.2A15_1<-diff.cyt$C_2A15_1-diff.cyt$C_293_1
diff.cyt$diff.C.2A15_2<-diff.cyt$C_2A15_2-diff.cyt$C_293_2
diff.cyt$diff.C.3C2_1<-diff.cyt$C_3C2_1-diff.cyt$C_293_1
diff.cyt$diff.C.3C2_2<-diff.cyt$C_3C2_2-diff.cyt$C_293_2
lfc.cyt[lfc.cyt==0]<-NA
lfc.cyt<-log2(lfc.cyt)
colnames(lfc.cyt)<-paste0("lfc.cyt.tot.",colnames(lfc.cyt))

lfc.mem<-mem/tot
diff.mem<-lfc.mem
diff.mem[diff.mem<1]<-NA
diff.mem$diff.M.2A15_1<-diff.mem$M_2A15_1-diff.mem$M_293_1
diff.mem$diff.M.2A15_2<-diff.mem$M_2A15_2-diff.mem$M_293_2
diff.mem$diff.M.3C2_1<-diff.mem$M_3C2_1-diff.mem$M_293_1
diff.mem$diff.M.3C2_2<-diff.mem$M_3C2_2-diff.mem$M_293_2
lfc.mem[lfc.mem==0]<-NA
lfc.mem<-log2(lfc.mem)
colnames(lfc.mem)<-paste0("lfc.mem.tot.",colnames(lfc.mem))


lfc.nuc<-nuc/tot
diff.nuc<-lfc.nuc
diff.nuc[diff.nuc<1]<-NA
diff.nuc$diff.N.2A15_1<-diff.nuc$N_2A15_1-diff.nuc$N_293_1
diff.nuc$diff.N.2A15_2<-diff.nuc$N_2A15_2-diff.nuc$N_293_2
diff.nuc$diff.N.3C2_1<-diff.nuc$N_3C2_1-diff.nuc$N_293_1
diff.nuc$diff.N.3C2_2<-diff.nuc$N_3C2_2-diff.nuc$N_293_2
lfc.nuc[lfc.nuc==0]<-NA
lfc.nuc<-log2(lfc.nuc)
colnames(lfc.nuc)<-paste0("lfc.nuc.tot.",colnames(lfc.nuc))

nuc[nuc==0]<-NA

lfc.cyt.nuc<-cyto/nuc
lfc.cyt.nuc[lfc.cyt.nuc==0]<-NA
lfc.cyt.nuc<-log2(lfc.cyt.nuc)
colnames(lfc.cyt.nuc)<-paste0("lfc.cyt.nuc.",colnames(lfc.cyt.nuc))

lfc.mem.nuc<-mem/nuc
lfc.mem.nuc[lfc.mem.nuc==0]<-NA
lfc.mem.nuc<-log2(lfc.mem.nuc)
colnames(lfc.mem.nuc)<-paste0("lfc.mem.nuc.",colnames(lfc.mem.nuc))

cyto<-fin[,colnames(fin)[grepl("C_", colnames(fin))]]
mem<-fin[,colnames(fin)[grepl("M_", colnames(fin))]]
nuc<-fin[,colnames(fin)[grepl("N_", colnames(fin))]]
tot<-fin[,colnames(fin)[grepl("T_", colnames(fin))]]

cyto[cyto==0]<-NA

lfc.nuc.cyt<-nuc/cyto
lfc.nuc.cyt[lfc.nuc.cyt==0]<-NA
lfc.nuc.cyt<-log2(lfc.nuc.cyt)
colnames(lfc.nuc.cyt)<-paste0("lfc.nuc.cyt.",colnames(lfc.nuc.cyt))

lfc.mem.cyt<-mem/cyto
lfc.mem.cyt[lfc.mem.cyt==0]<-NA
lfc.mem.cyt<-log2(lfc.mem.cyt)
colnames(lfc.mem.cyt)<-paste0("lfc.mem.cyt.",colnames(lfc.mem.cyt))


tot$lfc.tot.2A15_1<-log2(tot$T_2A15_1/tot$T_293_1)
tot$lfc.tot.2A15_2<-log2(tot$T_2A15_2/tot$T_293_2)
tot$lfc.tot.3C2_1<-log2(tot$T_3C2_1/tot$T_293_1)
tot$lfc.tot.3C2_2<-log2(tot$T_3C2_2/tot$T_293_2)

cyto$lfc.dcyto.2A15_1<-log2(cyto$C_2A15_1/cyto$C_293_1)
cyto$lfc.dcyto.2A15_2<-log2(cyto$C_2A15_2/cyto$C_293_2)
cyto$lfc.dcyto.3C2_1<-log2(cyto$C_3C2_1/cyto$C_293_1)
cyto$lfc.dcyto.3C2_2<-log2(cyto$C_3C2_2/cyto$C_293_2)

mem$lfc.dmem.2A15_1<-log2(mem$M_2A15_1/mem$M_293_1)
mem$lfc.dmem.2A15_2<-log2(mem$M_2A15_2/mem$M_293_2)
mem$lfc.dmem.3C2_1<-log2(mem$M_3C2_1/mem$M_293_1)
mem$lfc.dmem.3C2_2<-log2(mem$M_3C2_2/mem$M_293_2)

nuc$lfc.dnuc.2A15_1<-log2(nuc$N_2A15_1/nuc$N_293_1)
nuc$lfc.dnuc.2A15_2<-log2(nuc$N_2A15_2/nuc$N_293_2)
nuc$lfc.dnuc.3C2_1<-log2(nuc$N_3C2_1/nuc$N_293_1)
nuc$lfc.dnuc.3C2_2<-log2(nuc$N_3C2_2/nuc$N_293_2)

fin<-cbind(fin, lfc.cyt, lfc.mem, lfc.nuc, lfc.cyt.nuc, lfc.mem.nuc, lfc.nuc.cyt, lfc.mem.cyt,
           tot[,colnames(tot)[grepl("lfc", colnames(tot))]],
           cyto[,colnames(cyto)[grepl("dcyto", colnames(cyto))]],
           mem[,colnames(mem)[grepl("dmem", colnames(mem))]],
           nuc[,colnames(nuc)[grepl("dnuc", colnames(nuc))]],
           diff.cyt[,colnames(diff.cyt)[grepl("diff", colnames(diff.cyt))]],
           diff.mem[,colnames(diff.mem)[grepl("diff", colnames(diff.mem))]],
           diff.nuc[,colnames(diff.nuc)[grepl("diff", colnames(diff.nuc))]])

mel<-melt(fin, measure.vars = colnames(fin)[grepl("lfc.",colnames(fin))], id.vars = c("Symbol", "Annotation","target","Transmembrane.helices", 
                                                                "Cleavage.site..Signalp.","Low.complexity..Seg."))

ggplot(mel, aes(variable,value, fill=target))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(mel, aes(variable,value, fill=Transmembrane.helices))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(mel, aes(variable,value, fill=Cleavage.site..Signalp.))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot()+geom_point(aes(fin$lfc.mem.tot.M_293_1, fin$lfc.mem.tot.M_2A15_1 ), colour="grey")+
  geom_point(aes(fin[fin$target=="target","lfc.mem.tot.M_293_1"], fin[fin$target=="target","lfc.mem.tot.M_2A15_1"] ), colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(fin$lfc.mem.nuc.M_293_1, fin$lfc.mem.nuc.M_2A15_1 ), colour="grey")+
  geom_point(aes(fin[fin$target=="target","lfc.mem.nuc.M_293_1"], fin[fin$target=="target","lfc.mem.nuc.M_2A15_1"] ), colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(fin$lfc.mem.cyt.M_293_1, fin$lfc.mem.cyt.M_2A15_1 ), colour="grey")+
  geom_point(aes(fin[fin$target=="target","lfc.mem.cyt.M_293_1"], fin[fin$target=="target","lfc.mem.cyt.M_2A15_1"] ), colour="orange")+
  geom_abline(slope=1)


ggplot()+geom_point(aes(fin$lfc.mem.cyt.M_293_2, fin$lfc.mem.cyt.M_2A15_2 ), colour="grey")+
  geom_point(aes(fin[fin$target=="target","lfc.mem.cyt.M_293_2"], fin[fin$target=="target","lfc.mem.cyt.M_2A15_2"] ), colour="orange")+
  geom_abline(slope=1)


ggplot()+geom_point(aes(fin$lfc.nuc.cyt.N_293_1, fin$lfc.nuc.cyt.N_2A15_1 ), colour="grey")+
  geom_point(aes(fin[fin$target=="target","lfc.nuc.cyt.N_293_1"], fin[fin$target=="target","lfc.nuc.cyt.N_2A15_1"] ), colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(fin$lfc.dnuc.2A15_1, fin$lfc.dmem.2A15_1 ), colour="grey")+
  geom_point(aes(fin[fin$target=="target","lfc.dnuc.2A15_1"], fin[fin$target=="target","lfc.dmem.2A15_1"] ), colour="orange")+
  geom_abline(slope=1)



ggplot()+geom_point(aes(fin$lfc.mem.cyt.M_293_1, fin$lfc.mem.cyt.M_3C2_1 ), colour="grey")+
  geom_point(aes(fin[fin$target=="target","lfc.mem.cyt.M_293_1"], fin[fin$target=="target","lfc.mem.cyt.M_3C2_1"] ), colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(fin$lfc.mem.cyt.M_293_2, fin$lfc.mem.cyt.M_3C2_2 ), colour="grey")+
  geom_point(aes(fin[fin$target=="target","lfc.mem.cyt.M_293_2"], fin[fin$target=="target","lfc.mem.cyt.M_3C2_2"] ), colour="orange")+
  geom_abline(slope=1)



ggplot()+geom_point(aes(fin$lfc.mem.nuc.M_293_1, fin$lfc.mem.nuc.M_2A15_1 ), colour="grey")+
  geom_point(aes(fin[fin$target=="target","lfc.mem.nuc.M_293_1"], fin[fin$target=="target","lfc.mem.nuc.M_2A15_1"] ), colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(fin$lfc.mem.nucM_293_2, fin$lfc.mem.nucM_2A15_2 ), colour="grey")+
  geom_point(aes(fin[fin$target=="target","lfc.mem.nucM_293_2"], fin[fin$target=="target","lfc.mem.nucM_2A15_2"] ), colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(fin$lfc.mem.nucM_293_1, fin$lfc.mem.nucM_3C2_1 ), colour="grey")+
  geom_point(aes(fin[fin$target=="target","lfc.mem.nucM_293_1"], fin[fin$target=="target","lfc.mem.nucM_3C2_1"] ), colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(fin$lfc.cyt.nucC_293_1, fin$lfc.cyt.nucC_3C2_1 ), colour="grey")+
  geom_point(aes(fin[fin$target=="target","lfc.cyt.nucC_293_1"], fin[fin$target=="target","lfc.cyt.nucC_3C2_1"] ), colour="orange")+
  geom_abline(slope=1)


ggplot()+geom_point(aes(fin$lfc.mem.M_293_2, fin$lfc.mem.M_2A15_2 ), colour="grey")+
  geom_point(aes(fin[fin$target=="target","lfc.mem.M_293_2"], fin[fin$target=="target","lfc.mem.M_2A15_2"] ), colour="orange")+
  geom_abline(slope=1)


ggplot()+geom_point(aes(fin$lfc.tot.2A15_1, fin$lfc.mem.M_2A15_1 ), colour="grey")+
  geom_point(aes(fin[fin$target=="target","lfc.tot.2A15_1"], fin[fin$target=="target","lfc.mem.M_2A15_1"] ), colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(fin$lfc.tot.2A15_1, fin$lfc.tot.3C2_1 ), colour="grey")+
  geom_point(aes(fin[fin$target=="target","lfc.tot.2A15_1"], fin[fin$target=="target","lfc.tot.3C2_1"] ), colour="orange")+
  geom_abline(slope=1)


ggplot()+geom_point(aes(fin$lfc.tot.2A15_2, fin$lfc.mem.M_2A15_2 ), colour="grey")+
  geom_point(aes(fin[fin$target=="target","lfc.tot.2A15_2"], fin[fin$target=="target","lfc.mem.M_2A15_2"] ), colour="orange")+
  geom_abline(slope=1)


ggplot()+geom_point(aes(fin$lfc.tot.2A15_1, fin$lfc.mem.M_293_1 ), colour="grey")+
  geom_point(aes(fin[fin$target=="target","lfc.tot.2A15_1"], fin[fin$target=="target","lfc.mem.M_293_1"] ), colour="orange")+
  geom_abline(slope=1)


ggplot()+geom_point(aes(fin$lfc.mem.M_293_1, fin$lfc.mem.M_293_2 ), colour="grey")+
  geom_point(aes(fin[fin$target=="target","lfc.mem.M_293_1"], fin[fin$target=="target","lfc.mem.M_293_2"] ), colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(fin$lfc.mem.M_293_2, fin$lfc.mem.M_2A15_2 ), colour="grey")+
  geom_point(aes(fin[fin$target=="target","lfc.mem.M_293_2"], fin[fin$target=="target","lfc.mem.M_2A15_2"] ), colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(fin$lfc.mem.M_293_1, fin$lfc.mem.M_3C2_1 ), colour="grey")+
  geom_point(aes(fin[fin$target=="target","lfc.mem.M_293_1"], fin[fin$target=="target","lfc.mem.M_3C2_1"] ), colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(fin$lfc.mem.M_293_2, fin$lfc.mem.M_3C2_2 ), colour="grey")+
  geom_point(aes(fin[fin$target=="target","lfc.mem.M_293_2"], fin[fin$target=="target","lfc.mem.M_3C2_2"] ), colour="orange")+
  geom_abline(slope=1)



ggplot()+geom_point(aes(fin[fin$Annotation=="protein_coding","lfc.mem.tot.M_293_1"], fin[fin$Annotation=="protein_coding","lfc.mem.tot.M_2A15_1"] ), colour="grey")+
  geom_point(aes(fin[fin$target=="target" & fin$Annotation=="protein_coding","lfc.mem.tot.M_293_1"], fin[fin$target=="target"& fin$Annotation=="protein_coding","lfc.mem.tot.M_2A15_1"] ), colour="orange")+
  geom_abline(slope=1)


ggplot()+geom_point(aes(fin[fin$Annotation=="protein_coding","lfc.cyt.C_293_1"], fin[fin$Annotation=="protein_coding","lfc.cyt.C_2A15_1"] ), colour="grey")+
  geom_point(aes(fin[fin$target=="target" & fin$Annotation=="protein_coding","lfc.cyt.C_293_1"], fin[fin$target=="target"& fin$Annotation=="protein_coding","lfc.cyt.C_2A15_1"] ), colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(fin[fin$Annotation=="protein_coding","lfc.nuc.N_293_1"], fin[fin$Annotation=="protein_coding","lfc.nuc.N_2A15_1"] ), colour="grey")+
  geom_point(aes(fin[fin$target=="target" & fin$Annotation=="protein_coding","lfc.nuc.N_293_1"], fin[fin$target=="target"& fin$Annotation=="protein_coding","lfc.nuc.N_2A15_1"] ), colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(fin[fin$Annotation=="protein_coding","lfc.mem.M_293_2"], fin[fin$Annotation=="protein_coding","lfc.mem.M_2A15_2"] ), colour="grey")+
  geom_point(aes(fin[fin$target=="target" & fin$Annotation=="protein_coding","lfc.mem.M_293_2"], fin[fin$target=="target"& fin$Annotation=="protein_coding","lfc.mem.M_2A15_2"] ), colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(fin[fin$Annotation=="protein_coding","lfc.cyt.C_293_2"], fin[fin$Annotation=="protein_coding","lfc.cyt.C_2A15_2"] ), colour="grey")+
  geom_point(aes(fin[fin$target=="target" & fin$Annotation=="protein_coding","lfc.cyt.C_293_2"], fin[fin$target=="target"& fin$Annotation=="protein_coding","lfc.cyt.C_2A15_2"] ), colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(fin[fin$Annotation=="protein_coding","lfc.nuc.N_293_2"], fin[fin$Annotation=="protein_coding","lfc.nuc.N_2A15_2"] ), colour="grey")+
  geom_point(aes(fin[fin$target=="target" & fin$Annotation=="protein_coding","lfc.nuc.N_293_2"], fin[fin$target=="target"& fin$Annotation=="protein_coding","lfc.nuc.N_2A15_2"] ), colour="orange")+
  geom_abline(slope=1)


ggplot()+geom_point(aes(fin[fin$Annotation=="protein_coding","lfc.mem.M_293_1"], fin[fin$Annotation=="protein_coding","lfc.mem.M_2A15_1"] ), colour="grey")+
  geom_point(aes(fin[fin$Transmembrane.helices=="TMhelix" & fin$Annotation=="protein_coding","lfc.mem.M_293_1"], fin[fin$Transmembrane.helices=="TMhelix"& fin$Annotation=="protein_coding","lfc.mem.M_2A15_1"] ), colour="orange")+
  geom_abline(slope=1)

ggplot()+geom_point(aes(fin[fin$Annotation=="protein_coding","lfc.mem.M_293_1"], fin[fin$Annotation=="protein_coding","lfc.mem.M_2A15_1"] ), colour="grey")+
  geom_point(aes(fin[fin$Cleavage.site..Signalp.=="SignalP" & fin$Annotation=="protein_coding","lfc.mem.M_293_1"], fin[fin$Cleavage.site..Signalp.=="SignalP"& fin$Annotation=="protein_coding","lfc.mem.M_2A15_1"] ), colour="orange")+
  geom_abline(slope=1)


ggplot(fin, aes(log2(diff.M.2A15_1), colour=target))+stat_ecdf()
ggplot(fin, aes(log2(diff.M.2A15_2), colour=target))+stat_ecdf()
ggplot(fin, aes(log2(diff.M.3C2_1), colour=target))+stat_ecdf()
ggplot(fin, aes(log2(diff.M.3C2_2), colour=target))+stat_ecdf()



ggplot(fin, aes(diff.C.2A15_1, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(fin, aes(diff.C.3C2_1, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(fin, aes(diff.C.2A15_2, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(fin, aes(diff.C.3C2_2, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))

ggplot(fin, aes(diff.N.2A15_1, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(fin, aes(diff.N.3C2_1, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(fin, aes(diff.N.2A15_2, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(fin, aes(diff.N.3C2_2, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))

wilcox.test(fin[fin$target=="target", "diff.M.2A15_1"], fin[fin$target=="nontarget", "diff.M.2A15_1"])
wilcox.test(fin[fin$target=="target", "diff.M.2A15_2"], fin[fin$target=="nontarget", "diff.M.2A15_2"])
wilcox.test(fin[fin$target=="target", "diff.M.3C2_1"], fin[fin$target=="nontarget", "diff.M.3C2_1"])
wilcox.test(fin[fin$target=="target", "diff.M.3C2_2"], fin[fin$target=="nontarget", "diff.M.3C2_2"])

wilcox.test(fin[fin$target=="target", "diff.C.2A15_1"], fin[fin$target=="nontarget", "diff.C.2A15_1"])
wilcox.test(fin[fin$target=="target", "diff.C.2A15_2"], fin[fin$target=="nontarget", "diff.C.2A15_2"])
wilcox.test(fin[fin$target=="target", "diff.C.3C2_1"], fin[fin$target=="nontarget", "diff.C.3C2_1"])
wilcox.test(fin[fin$target=="target", "diff.C.3C2_2"], fin[fin$target=="nontarget", "diff.C.3C2_2"])

ggplot(fin, aes(lfc.tot.3C2_2, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-4,4))
ggplot(fin, aes(lfc.tot.3C2_1, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-4,4))
ggplot(fin, aes(lfc.tot.2A15_2, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-4,4))
ggplot(fin, aes(lfc.tot.2A15_1, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-4,4))
ggplot(fin, aes(lfc.dmem.3C2_1, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-4,4))
ggplot(fin, aes(lfc.dmem.2A15_1, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-4,4))
ggplot(fin, aes(lfc.dnuc.3C2_1, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-4,4))
ggplot(fin, aes(lfc.dnuc.2A15_1, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-4,4))
ggplot(fin, aes(lfc.dnuc.2A15_2, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-4,4))
ggplot(fin, aes(lfc.dnuc.3C2_2, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-4,4))
ggplot(fin, aes(lfc.dcyto.3C2_1, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-4,4))
ggplot(fin, aes(lfc.dcyto.2A15_1, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-4,4))
ggplot(fin, aes(lfc.dcyto.3C2_1, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-4,4))
ggplot(fin, aes(lfc.dcyto.2A15_1, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-4,4))
ggplot(fin, aes(lfc.mem.nucM_293_1, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-4,4))
ggplot(fin, aes(lfc.mem.nucM_2A15_1, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-4,4))
ggplot(fin, aes(lfc.mem.nucM_3C2_1, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-4,4))
ggplot(fin, aes(lfc.cyt.nucC_293_1, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-4,4))
ggplot(fin, aes(lfc.cyt.nucC_2A15_1, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-4,4))
ggplot(fin, aes(lfc.cyt.nucC_3C2_1, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-4,4))

gplot(fin, aes(lfc.mem.cyt.M_293_1, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-3,10))
ggplot(fin, aes(lfc.mem.cyt.M_2A15_1, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-3,10))


ggplot(fin, aes(lfc.mem.M_293_1, colour=target))+stat_ecdf()
ggplot(fin, aes(lfc.mem.M_293_2, colour=target))+stat_ecdf()
ggplot(fin, aes(lfc.cyt.C_293_1, colour=target))+stat_ecdf()
ggplot(fin, aes(lfc.cyt.C_293_2, colour=target))+stat_ecdf()
ggplot(fin, aes(lfc.nuc.N_293_2, colour=target))+stat_ecdf()
ggplot(fin, aes(lfc.nuc.N_293_2, colour=target))+stat_ecdf()
ggplot(fin, aes(lfc.mem.M_2A15_1, colour=target))+stat_ecdf()
ggplot(fin, aes(lfc.mem.M_2A15_2, colour=target))+stat_ecdf()
ggplot(fin, aes(lfc.cyt.C_2A15_1, colour=target))+stat_ecdf()
ggplot(fin, aes(lfc.cyt.C_2A15_2, colour=target))+stat_ecdf()
ggplot(fin, aes(lfc.nuc.N_2A15_1, colour=target))+stat_ecdf()
ggplot(fin, aes(lfc.nuc.N_2A15_2, colour=target))+stat_ecdf()


# library(GGally)
# ggpairs(fin[,c(55:60,47)], aes(colour = factor(target), alpha = 0.4))




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
lfc.mem.cyt.293<-data.frame(res[,c(2,6)])
colnames(lfc.mem.cyt.293)<-paste0(colnames(lfc.mem.cyt.293),".mem.cyt.293")

res <- results(dds, contrast=c("Condition","M_2A15","C_2A15")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.mem.cyt.2A15<-data.frame(res[,c(2,6)])
colnames(lfc.mem.cyt.2A15)<-paste0(colnames(lfc.mem.cyt.2A15),".mem.cyt.2A15")

res <- results(dds, contrast=c("Condition","M_3C2","C_3C2")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.mem.cyt.3C2<-data.frame(res[,c(2,6)])
colnames(lfc.mem.cyt.3C2)<-paste0(colnames(lfc.mem.cyt.3C2),".mem.cyt.3C2")

dese<-cbind(lfc.mem.cyt.293, lfc.mem.cyt.2A15, lfc.mem.cyt.3C2)

fin<-merge(fin, dese, by.x="Row.names", by.y="row.names")

fin<-subset(fin, tpm_cutoff>10)

ggplot()+geom_point(aes(fin$log2FoldChange.mem.cyt.2A15, fin$log2FoldChange.mem.cyt.293 ), colour="grey")+
  geom_point(aes(fin[fin$target=="target","log2FoldChange.mem.cyt.2A15"], fin[fin$target=="target","log2FoldChange.mem.cyt.293"] ),colour="orange")+
  geom_text(label=row.names(fin))+geom_abline(slope=1)

ggplot()+geom_point(aes(fin$log2FoldChange.mem.cyt.3C2, fin$log2FoldChange.mem.cyt.293 ), colour="grey")+
  geom_point(aes(fin[fin$target=="target","log2FoldChange.mem.cyt.2A15"], fin[fin$target=="target","log2FoldChange.mem.cyt.293"] ), colour="orange")+
  geom_abline(slope=1)

subset(fin, log2FoldChange.mem.cyt.2A15-log2FoldChange.mem.cyt.293<(-0.3) & target=="target" & tpm_cutoff>10, select=c("Symbol", "conv_CDS","reads_CDS"))

ggplot()+ geom_point(aes(fin[fin$target=="target","log2FoldChange.mem.cyt.2A15"], fin[fin$target=="target","conv_CDS"] ), colour="orange")+
  geom_abline(slope=1)


cnt<-tab[,colnames(tab)[grepl("M_", colnames(tab))]]

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

dds$Condition = relevel(dds$Condition,"M_293") #everything relative to Mock
dds<-DESeq(dds)

resultsNames(dds)
res <- results(dds, contrast=c("Condition","M_2A15","M_293")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.mem.2A15.293<-data.frame(res[,c(2,6)])
colnames(lfc.mem.2A15.293)<-paste0(colnames(lfc.mem.2A15.293),".mem.2A15.293")

res <- results(dds, contrast=c("Condition","M_3C2","M_293")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.mem.3C2.293<-data.frame(res[,c(2,6)])
colnames(lfc.mem.3C2.293)<-paste0(colnames(lfc.mem.3C2.293),".mem.3C2.293")

dese<-cbind(lfc.mem.2A15.293, lfc.mem.3C2.293)

fin<-merge(fin, dese, by.x="Row.names", by.y="row.names")

ggplot()+geom_point(aes(fin$log2FoldChange.mem.2A15.293, fin$log2FoldChange.mem.cyt.293 ), colour="grey")+
  geom_point(aes(fin[fin$target=="target","log2FoldChange.mem.2A15.293"], fin[fin$target=="target","log2FoldChange.mem.cyt.293"] ), colour="orange")+
  geom_abline(slope=1)+xlim(-3,3)+ylim(-5,10)

ggplot()+geom_point(aes(fin$log2FoldChange.mem.2A15.293, fin$log2FoldChange.mem.cyt.2A15 ), colour="grey")+
  geom_point(aes(fin[fin$target=="target","log2FoldChange.mem.2A15.293"], fin[fin$target=="target","log2FoldChange.mem.cyt.2A15"] ), colour="orange")+
  geom_abline(slope=1)+xlim(-3,3)+ylim(-5,10)



cnt<-tab[,colnames(tab)[grepl("M_", colnames(tab))|grepl("T_", colnames(tab))]]

coldata<-data.frame(SampleID=colnames(cnt),
                    Condition=gsub("_.*","",gsub(".*;","",sub("_",";",colnames(cnt)))),
                    Fraction=gsub("_.*","",colnames(cnt)),
                    Batch=gsub(".*_","",colnames(cnt)))
coldata<-apply(coldata,2,as.factor)

dds<-DESeqDataSetFromMatrix(countData = round(cnt),
                            colData = coldata, 
                            design =~ Fraction + Condition +Fraction:Condition)
dds<-estimateSizeFactors(dds)
# sizeFactors(dds)<-c(estimateSizeFactorsForMatrix(cnt[,1:6]), estimateSizeFactorsForMatrix(cnt[,7:12]) )

vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Condition", "Fraction","Batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Batch)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 

ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Fraction)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 


dds$Condition = relevel(dds$Condition,"293") #everything relative to Mock
dds$Fraction = relevel(dds$Fraction,"T") #everything relative to RNAseq
dds<-DESeq(dds, test="LRT", reduced=~Fraction+Condition)

resultsNames(dds)
res <- results(dds, contrast=list("FractionM.Condition2A15")) 
summary(res)
res<-res[order(res$padj),]
DESeq2::plotMA(res)
DESeq2::plotDispEsts(dds)
res <- results(dds, contrast=list("FractionM.Condition3C2")) 
summary(res)
res<-res[order(res$padj),]
DESeq2::plotMA(res)






#systematic deseq2

cnt<-tab[,colnames(tab)[grepl("C_", colnames(tab))|grepl("N_", colnames(tab))|grepl("M_", colnames(tab))|grepl("T_", colnames(tab))]]

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


dds$Condition = relevel(dds$Condition,"C_293") #everything relative to C_293
dds<-DESeq(dds)

resultsNames(dds)
res <- results(dds, contrast=c("Condition","M_293","C_293")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-6,6))
lfc.mem.cyt.293<-data.frame(res[,c(2,6)])
colnames(lfc.mem.cyt.293)<-paste0(colnames(lfc.mem.cyt.293),".mem.cyt.293")

dds$Condition = relevel(dds$Condition,"C_2A15") #everything relative to C_2A15
dds<-DESeq(dds)
res <- results(dds, contrast=c("Condition","M_2A15","C_2A15")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.mem.cyt.2A15<-data.frame(res[,c(2,6)])
colnames(lfc.mem.cyt.2A15)<-paste0(colnames(lfc.mem.cyt.2A15),".mem.cyt.2A15")

dds$Condition = relevel(dds$Condition,"C_3C2") #everything relative to C_3C2
dds<-DESeq(dds)
res <- results(dds, contrast=c("Condition","M_3C2","C_3C2")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.mem.cyt.3C2<-data.frame(res[,c(2,6)])
colnames(lfc.mem.cyt.3C2)<-paste0(colnames(lfc.mem.cyt.3C2),".mem.cyt.3C2")




dds$Condition = relevel(dds$Condition,"C_293") #everything relative to C_293
dds<-DESeq(dds)
res <- results(dds, contrast=c("Condition","N_293","C_293")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-6,6))
lfc.nuc.cyt.293<-data.frame(res[,c(2,6)])
colnames(lfc.nuc.cyt.293)<-paste0(colnames(lfc.nuc.cyt.293),".nuc.cyt.293")

dds$Condition = relevel(dds$Condition,"C_2A15") #everything relative to C_293
dds<-DESeq(dds)
res <- results(dds, contrast=c("Condition","N_2A15","C_2A15")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.nuc.cyt.2A15<-data.frame(res[,c(2,6)])
colnames(lfc.nuc.cyt.2A15)<-paste0(colnames(lfc.nuc.cyt.2A15),".nuc.cyt.2A15")

dds$Condition = relevel(dds$Condition,"C_3C2") #everything relative to C_293
dds<-DESeq(dds)
res <- results(dds, contrast=c("Condition","N_3C2","C_3C2")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.nuc.cyt.3C2<-data.frame(res[,c(2,6)])
colnames(lfc.nuc.cyt.3C2)<-paste0(colnames(lfc.nuc.cyt.3C2),".nuc.cyt.3C2")


dds$Condition = relevel(dds$Condition,"T_293") #everything relative to T_293
dds<-DESeq(dds)
resultsNames(dds)
res <- results(dds, contrast=c("Condition","M_293","T_293")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-6,6))
lfc.mem.tot.293<-data.frame(res[,c(2,6)])
colnames(lfc.mem.tot.293)<-paste0(colnames(lfc.mem.tot.293),".mem.tot.293")

dds$Condition = relevel(dds$Condition,"T_2A15") #everything relative to C_293
dds<-DESeq(dds)
res <- results(dds, contrast=c("Condition","M_2A15","T_2A15")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.mem.tot.2A15<-data.frame(res[,c(2,6)])
colnames(lfc.mem.tot.2A15)<-paste0(colnames(lfc.mem.tot.2A15),".mem.tot.2A15")

dds$Condition = relevel(dds$Condition,"T_3C2") #everything relative to C_293
dds<-DESeq(dds)
res <- results(dds, contrast=c("Condition","M_3C2","T_3C2")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.mem.tot.3C2<-data.frame(res[,c(2,6)])
colnames(lfc.mem.tot.3C2)<-paste0(colnames(lfc.mem.tot.3C2),".mem.tot.3C2")

dds$Condition = relevel(dds$Condition,"T_293") #everything relative to C_293
dds<-DESeq(dds)
res <- results(dds, contrast=c("Condition","N_293","T_293")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-6,6))
lfc.nuc.tot.293<-data.frame(res[,c(2,6)])
colnames(lfc.nuc.tot.293)<-paste0(colnames(lfc.nuc.tot.293),".nuc.tot.293")

dds$Condition = relevel(dds$Condition,"T_2A15") #everything relative to C_293
dds<-DESeq(dds)
res <- results(dds, contrast=c("Condition","N_2A15","T_2A15")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.nuc.tot.2A15<-data.frame(res[,c(2,6)])
colnames(lfc.nuc.tot.2A15)<-paste0(colnames(lfc.nuc.tot.2A15),".nuc.tot.2A15")

dds$Condition = relevel(dds$Condition,"T_3C2") #everything relative to C_293
dds<-DESeq(dds)
res <- results(dds, contrast=c("Condition","N_3C2","T_3C2")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.nuc.tot.3C2<-data.frame(res[,c(2,6)])
colnames(lfc.nuc.tot.3C2)<-paste0(colnames(lfc.nuc.tot.3C2),".nuc.tot.3C2")

dds$Condition = relevel(dds$Condition,"T_293") #everything relative to C_293
dds<-DESeq(dds)
res <- results(dds, contrast=c("Condition","C_293","T_293")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-6,6))
lfc.cyt.tot.293<-data.frame(res[,c(2,6)])
colnames(lfc.cyt.tot.293)<-paste0(colnames(lfc.cyt.tot.293),".cyt.tot.293")

dds$Condition = relevel(dds$Condition,"T_2A15") #everything relative to C_293
dds<-DESeq(dds)
res <- results(dds, contrast=c("Condition","C_2A15","T_2A15")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.cyt.tot.2A15<-data.frame(res[,c(2,6)])
colnames(lfc.cyt.tot.2A15)<-paste0(colnames(lfc.cyt.tot.2A15),".cyt.tot.2A15")

dds$Condition = relevel(dds$Condition,"T_3C2") #everything relative to C_293
dds<-DESeq(dds)
res <- results(dds, contrast=c("Condition","C_3C2","T_3C2")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-4,4))
lfc.cyt.tot.3C2<-data.frame(res[,c(2,6)])
colnames(lfc.cyt.tot.3C2)<-paste0(colnames(lfc.cyt.tot.3C2),".cyt.tot.3C2")



###multifactor design for differences between 2A14 and 3C2
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
lfc.mem.cyt.2A15.293<-data.frame(res[,c(2,6)])
colnames(lfc.mem.cyt.2A15.293)<-paste0(colnames(lfc.mem.cyt.2A15.293),".mem.cyt.2A15.293")

res <- results(dds, contrast=list("FractionM.Condition3C2")) 
summary(res)
DESeq2::plotMA(res)
lfc.mem.cyt.3C2.293<-data.frame(res[,c(2,6)])
colnames(lfc.mem.cyt.3C2.293)<-paste0(colnames(lfc.mem.cyt.3C2.293),".mem.cyt.3C2.293")

res <- results(dds, contrast=list("FractionN.Condition2A15")) 
summary(res)
DESeq2::plotMA(res)
lfc.nuc.cyt.2A15.293<-data.frame(res[,c(2,6)])
colnames(lfc.nuc.cyt.2A15.293)<-paste0(colnames(lfc.nuc.cyt.2A15.293),".nuc.cyt.2A15.293")

res <- results(dds, contrast=list("FractionN.Condition3C2")) 
summary(res)
DESeq2::plotMA(res)
lfc.nuc.cyt.3C2.293<-data.frame(res[,c(2,6)])
colnames(lfc.nuc.cyt.3C2.293)<-paste0(colnames(lfc.nuc.cyt.3C2.293),".nuc.cyt.3C2.293")





dds$Condition = relevel(dds$Condition,"293") #everything relative to 293
dds$Fraction = relevel(dds$Fraction,"T") #everything relative to T
dds<-DESeq(dds, test="LRT", reduced=~Fraction+Condition)

resultsNames(dds)

res <- results(dds, contrast=list("FractionM.Condition2A15")) 
summary(res)
DESeq2::plotMA(res)
lfc.mem.tot.2A15.293<-data.frame(res[,c(2,6)])
colnames(lfc.mem.tot.2A15.293)<-paste0(colnames(lfc.mem.tot.2A15.293),".mem.tot.2A15.293")

res <- results(dds, contrast=list("FractionM.Condition3C2")) 
summary(res)
DESeq2::plotMA(res)
lfc.mem.tot.3C2.293<-data.frame(res[,c(2,6)])
colnames(lfc.mem.tot.3C2.293)<-paste0(colnames(lfc.mem.tot.3C2.293),".mem.tot.3C2.293")

res <- results(dds, contrast=list("FractionN.Condition2A15")) 
summary(res)
DESeq2::plotMA(res)
lfc.nuc.tot.2A15.293<-data.frame(res[,c(2,6)])
colnames(lfc.nuc.tot.2A15.293)<-paste0(colnames(lfc.nuc.tot.2A15.293),".nuc.tot.2A15.293")

res <- results(dds, contrast=list("FractionN.Condition3C2")) 
summary(res)
DESeq2::plotMA(res)
lfc.nuc.tot.3C2.293<-data.frame(res[,c(2,6)])
colnames(lfc.nuc.tot.3C2.293)<-paste0(colnames(lfc.nuc.tot.3C2.293),".nuc.tot.3C2.293")

res <- results(dds, contrast=list("FractionC.Condition2A15")) 
summary(res)
DESeq2::plotMA(res)
lfc.cyt.tot.2A15.293<-data.frame(res[,c(2,6)])
colnames(lfc.cyt.tot.2A15.293)<-paste0(colnames(lfc.cyt.tot.2A15.293),".cyt.tot.2A15.293")

res <- results(dds, contrast=list("FractionC.Condition3C2")) 
summary(res)
DESeq2::plotMA(res)
lfc.cyt.tot.3C2.293<-data.frame(res[,c(2,6)])
colnames(lfc.cyt.tot.3C2.293)<-paste0(colnames(lfc.cyt.tot.3C2.293),".cyt.tot.3C2.293")






dese<-cbind(lfc.mem.cyt.293, lfc.mem.cyt.2A15, lfc.mem.cyt.3C2, 
            lfc.nuc.cyt.293, lfc.nuc.cyt.2A15, lfc.nuc.cyt.3C2, 
            lfc.mem.tot.293, lfc.mem.tot.2A15, lfc.mem.tot.3C2, 
            lfc.nuc.tot.293, lfc.nuc.tot.2A15, lfc.nuc.tot.3C2,
            lfc.cyt.tot.293, lfc.cyt.tot.2A15, lfc.cyt.tot.3C2,
            lfc.mem.cyt.2A15.293, lfc.mem.cyt.3C2.293, 
            lfc.nuc.cyt.2A15.293, lfc.nuc.cyt.3C2.293, 
            lfc.mem.tot.2A15.293, lfc.mem.tot.3C2.293, 
            lfc.nuc.tot.2A15.293, lfc.nuc.tot.3C2.293, 
            lfc.cyt.tot.2A15.293, lfc.cyt.tot.3C2.293)
  
  


fin<-merge(fin, dese, by.x="Row.names", by.y="row.names")

rel<-subset(fin, tpm_cutoff>=20 & Annotation=="protein_coding")
ggplot()+geom_point(aes(rel$log2FoldChange.mem.cyt.293, rel$log2FoldChange.mem.cyt.2A15 ), colour="grey")+
  geom_point(aes(rel[rel$target=="target","log2FoldChange.mem.cyt.293"], rel[rel$target=="target","log2FoldChange.mem.cyt.2A15"] ), colour="orange")+
  geom_abline(slope=1)+xlim(-3,15)+ylim(-3,15)
ggplot()+geom_point(aes(rel$log2FoldChange.mem.cyt.293, rel$log2FoldChange.mem.cyt.3C2 ), colour="grey")+
  geom_point(aes(rel[rel$target=="target","log2FoldChange.mem.cyt.293"], rel[rel$target=="target","log2FoldChange.mem.cyt.3C2"] ), colour="orange")+
  geom_abline(slope=1)+xlim(-3,15)+ylim(-3,15)
