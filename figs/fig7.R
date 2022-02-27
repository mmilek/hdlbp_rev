#fig7
library(ggplot2)

dat<-read.delim("data/xenograft_all.txt", stringsAsFactors = F, dec = ",")

dat$status<-factor(dat$status, levels=c("WT", "KO"))
dat$time<-factor(dat$time, levels=c("7", "14", "21"))

#fig7f left panel
ggplot(dat, aes(time, log2(Volume), colour=status))+geom_boxplot()+
  geom_point( position = position_jitterdodge(jitter.width = 0.5))+xlab("time (days)")+ylab("log2 tumor volume")+
  scale_colour_manual(values=c("dodgerblue3", "orange3"))

ggplot(dat, aes(time, (Volume), colour=status))+geom_boxplot()+
  geom_point( position = position_jitterdodge(jitter.width = 0.5))+xlab("time (days)")+ylab("tumor volume")+
  scale_colour_manual(values=c("dodgerblue3", "orange3"))
ggplot(dat, aes(time, log2(Height), colour=status))+geom_boxplot()+
  geom_point( position = position_jitterdodge(jitter.width = 0.5))+xlab("time (days)")+ylab("log2 tumor height")+
  scale_colour_manual(values=c("dodgerblue3", "orange3"))
ggplot(dat, aes(time, (Width), colour=status))+geom_boxplot()+
  geom_point( position = position_jitterdodge(jitter.width = 0.5))+xlab("time (days)")+ylab("tumor width")+
  scale_colour_manual(values=c("dodgerblue3", "orange3"))

ggplot(dat, aes(time, (Length), colour=status))+geom_boxplot()+
  geom_point( position = position_jitterdodge(jitter.width = 0.5))+xlab("time (days)")+ylab("log2 tumor length")+
  scale_colour_manual(values=c("dodgerblue3", "orange3"))

#fig7f right panel
ggplot(dat, aes(status, (tumor_weight_mg), colour=status))+geom_boxplot()+
  geom_point( position = position_jitterdodge(jitter.width = 0.8))+xlab("")+ylab("tumor weight (mg)")+
  scale_colour_manual(values=c("dodgerblue3", "orange3"))

# write.table(dat, "source_data/fig7f.txt", quote=F, sep="\t", row.names = F)


t.test(dat[dat$status=="WT" & dat$time=="7", "Volume"], dat[dat$status=="KO" & dat$time=="7", "Volume"])
t.test(dat[dat$status=="WT" & dat$time=="14", "Volume"], dat[dat$status=="KO" & dat$time=="14", "Volume"])
t.test(dat[dat$status=="WT" & dat$time=="21", "Volume"], dat[dat$status=="KO" & dat$time=="21", "Volume"])

t.test(dat[dat$status=="WT" & dat$time=="7", "Length"], dat[dat$status=="KO" & dat$time=="7", "Length"])
t.test(dat[dat$status=="WT" & dat$time=="14", "Length"], dat[dat$status=="KO" & dat$time=="14", "Length"])
t.test(dat[dat$status=="WT" & dat$time=="21", "Length"], dat[dat$status=="KO" & dat$time=="21", "Length"])

t.test(dat[dat$status=="WT" & dat$time=="7", "Width"], dat[dat$status=="KO" & dat$time=="7", "Width"])
t.test(dat[dat$status=="WT" & dat$time=="14", "Width"], dat[dat$status=="KO" & dat$time=="14", "Width"])
t.test(dat[dat$status=="WT" & dat$time=="21", "Width"], dat[dat$status=="KO" & dat$time=="21", "Width"])

t.test(dat[dat$status=="WT" & dat$time=="7", "Height"], dat[dat$status=="KO" & dat$time=="7", "Height"])
t.test(dat[dat$status=="WT" & dat$time=="14", "Height"], dat[dat$status=="KO" & dat$time=="14", "Height"])
t.test(dat[dat$status=="WT" & dat$time=="21", "Height"], dat[dat$status=="KO" & dat$time=="21", "Height"])


t.test(dat[dat$status=="WT" & dat$time=="21", "tumor_weight_mg"], dat[dat$status=="KO" & dat$time=="21", "tumor_weight_mg"])


length(na.omit(dat[dat$status=="WT" & dat$time=="7", "Volume"]))
length(na.omit(dat[dat$status=="KO" & dat$time=="7", "Volume"]))
length(na.omit(dat[dat$status=="WT" & dat$time=="14", "Volume"]))
length(na.omit(dat[dat$status=="KO" & dat$time=="14", "Volume"]))
length(na.omit(dat[dat$status=="WT" & dat$time=="21", "Volume"]))
length(na.omit(dat[dat$status=="KO" & dat$time=="21", "Volume"]))
length(na.omit(dat[dat$status=="WT" & dat$time=="21", "tumor_weight_mg"]))
length(na.omit(dat[dat$status=="KO" & dat$time=="21", "tumor_weight_mg"]))


#figS7
tab<-read.delim("data/tumor_rnaseq_read_counts.txt", header=T)
tpm<-read.delim("data/tumor_rnaseq_tpm.txt", header=T)

##PCA

subs<-tab
row.names(subs)<-subs$gene_id
subs<-subs[,-1]
coldata<-data.frame(SampleID=colnames(subs),
                    Condition=sub("_[0-9]","",gsub(".*;","",sub("_",";",colnames(subs)))),
                    Batch=substr(colnames(subs), nchar(colnames(subs)), nchar(colnames(subs))),
                    Source=gsub("_.*","",colnames(subs)))

coldata<-apply(coldata,2,as.factor)

corrplot(cor(subs, use="pairwise.complete.obs"), type="upper", method="color",tl.col = "black")


dds<-DESeqDataSetFromMatrix(countData = round(subs),
                            colData = coldata, 
                            design = ~Batch+Condition+Source)


vsd <- vst(dds, blind=FALSE)

#figS7c lower panel
pcaData <- plotPCA(vsd, intgroup=c("Batch", "Source"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Batch, shape=Source)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 

# write.table(pcaData, "source_data/fig7c_lower.txt", quote=F, sep="\t", row.names =F )

#figS7c upper panel
pcaData <- plotPCA(vsd, intgroup=c("Batch", "Condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Batch, shape=Condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 

# write.table(pcaData, "source_data/fig7c_upper.txt", quote=F, sep="\t", row.names =F )


####deseq simple comparisons within fractions
cnt<-subs[,colnames(subs)[grepl("Tumor_", colnames(subs))]]

coldata<-data.frame(SampleID=colnames(cnt),
                    Condition=sub("_[0-9]","",gsub(".*;","",sub("_",";",colnames(cnt)))),
                    Batch=substr(colnames(cnt), nchar(colnames(cnt)), nchar(colnames(cnt))),
                    Source=gsub("_.*","",colnames(cnt)))

coldata<-apply(coldata,2,as.factor)


dds<-DESeqDataSetFromMatrix(countData = round(cnt),
                            colData = coldata, 
                            design =~ Condition )

dds<-estimateSizeFactors(dds)

dds$Condition = relevel(dds$Condition,"WT")
dds<-DESeq(dds)

resultsNames(dds)
res <- results(dds, contrast=c("Condition","HDLBP_KO","WT")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-8,8))
lfc.tumor.KO.WT<-data.frame(res[,c(1,2,6)])

colnames(lfc.tumor.KO.WT)<-paste0(colnames(lfc.tumor.KO.WT),".tumor.KO.WT")
counts.tumor.KO.WT<-cbind(counts(dds, normalize=F), counts(dds, normalize=T))
colnames(counts.tumor.KO.WT)[13:ncol(counts.tumor.KO.WT)]<-paste0("norm.", colnames(counts.tumor.KO.WT)[13:ncol(counts.tumor.KO.WT)])


cnt<-subs[,colnames(subs)[grepl("A_", colnames(subs))]]

coldata<-data.frame(SampleID=colnames(cnt),
                    Condition=sub("_[0-9]","",gsub(".*;","",sub("_",";",colnames(cnt)))),
                    Batch=substr(colnames(cnt), nchar(colnames(cnt)), nchar(colnames(cnt))),
                    Source=gsub("_.*","",colnames(cnt)))

coldata<-apply(coldata,2,as.factor)


dds<-DESeqDataSetFromMatrix(countData = round(cnt),
                            colData = coldata, 
                            design =~ Condition )

dds<-estimateSizeFactors(dds)

dds$Condition = relevel(dds$Condition,"WT")
dds<-DESeq(dds)

resultsNames(dds)
res <- results(dds, contrast=c("Condition","HDLBP_KO","WT")) 
summary(res)
DESeq2::plotMA(res, ylim=c(-8,8))
lfc.a549.KO.WT<-data.frame(res[,c(1,2,6)])

colnames(lfc.a549.KO.WT)<-paste0(colnames(lfc.a549.KO.WT),".a549.KO.WT")
counts.a549.KO.WT<-cbind(counts(dds, normalize=F), counts(dds, normalize=T))
colnames(counts.a549.KO.WT)[7:ncol(counts.a549.KO.WT)]<-paste0("norm.", colnames(counts.a549.KO.WT)[7:ncol(counts.a549.KO.WT)])

row.names(tpm)<-tpm$gene_id
tpm<-tpm[,-1]
colnames(tpm)<-paste0("tpm.", colnames(tpm))

fin<-merge(tpm, counts.a549.KO.WT, by="row.names")
fin<-merge(fin, lfc.a549.KO.WT, by.x="Row.names", by.y="row.names")
fin<-merge(fin, counts.tumor.KO.WT, by.x="Row.names", by.y="row.names")
fin<-merge(fin, lfc.tumor.KO.WT, by.x="Row.names", by.y="row.names")
colnames(fin)[1]<-"gene_id"

thr<-10
fin<-subset(fin, tpm.A_WT_1 >=thr |
              tpm.A_WT_2 >=thr |
              tpm.A_WT_3 >=thr |
              tpm.A_HDLBP_KO_1 >=thr |
              tpm.A_HDLBP_KO_2 >=thr |
              tpm.A_HDLBP_KO_3 >=thr |
              tpm.Tumor_WT_1 >=thr |
              tpm.Tumor_WT_2 >=thr |
              tpm.Tumor_WT_3 >=thr |
              tpm.Tumor_WT_4 >=thr |
              tpm.Tumor_WT_5 >=thr |
              tpm.Tumor_WT_6 >=thr |
              tpm.Tumor_WT_7 >=thr |
              tpm.Tumor_WT_8 >=thr |
              tpm.Tumor_HDLBP_KO_2 >=thr |
              tpm.Tumor_HDLBP_KO_3 >=thr |
              tpm.Tumor_HDLBP_KO_4 >=thr |
              tpm.Tumor_HDLBP_KO_5 >=thr )


mas<-read.delim("data/hdlbp_master_table_with_classes_uniq_tsig.txt", header=T)
inf<-subset(mas, select=c("gene_id","Symbol","gene_biotype", 
                          "tsig" , "loc_tar_CDS", "localization_cat" ))

exp<-merge(fin, inf, by="gene_id", all.x=T)

exp$reg_a549<-ifelse(exp$log2FoldChange.a549.KO.WT>=1 &
                       exp$padj.a549.KO.WT<0.1 , "up_a549",
                     ifelse(exp$log2FoldChange.a549.KO.WT<=(-1) &
                              exp$padj.a549.KO.WT<0.1, "down_a549",
                            ifelse(exp$log2FoldChange.a549.KO.WT<1 & 
                                     exp$log2FoldChange.a549.KO.WT>(-1) &
                                     exp$padj.a549.KO.WT>=0.1, "unchanged_a549", NA)))
exp$reg_tumor<-ifelse(exp$log2FoldChange.tumor.KO.WT>=1 &
                        exp$padj.tumor.KO.WT<0.1 , "up_tumor",
                      ifelse(exp$log2FoldChange.tumor.KO.WT<=(-1) &
                               exp$padj.tumor.KO.WT<0.1, "down_tumor",
                             ifelse(exp$log2FoldChange.tumor.KO.WT<1 & 
                                      exp$log2FoldChange.tumor.KO.WT>(-1) &
                                      exp$padj.tumor.KO.WT>=0.1 , "unchanged_tumor", "not_quantified")))

exp$reg_both<-paste0(exp$reg_a549,";", exp$reg_tumor)
ggplot(subset(exp, !is.na(reg_a549)), aes(reg_a549))+geom_bar()+coord_flip()
ggplot(subset(exp, !is.na(reg_tumor)), aes(reg_tumor))+geom_bar()+coord_flip()
ggplot(subset(exp, !is.na(reg_tumor) | !is.na(reg_a549)), aes(reg_both))+geom_bar()+coord_flip()



#supplemental tables
reg_tumor_up<-subset(exp, reg_tumor=="up_tumor" & Annotation=="protein_coding",
                     select = colnames(exp)[c(1,68,59:61,62:67,69:71,8:19)])
reg_tumor_up<-reg_tumor_up[order(reg_tumor_up$log2FoldChange.tumor.KO.WT, decreasing = T),]
reg_tumor_down<-subset(exp, reg_tumor=="down_tumor" & Annotation=="protein_coding",
                       select = colnames(exp)[c(1,68,59:61,62:67,69:71,8:19)])
reg_tumor_down<-reg_tumor_down[order(reg_tumor_down$log2FoldChange.tumor.KO.WT, decreasing = F),]

#write.table(reg_tumor_up, "data_tables/tumor_up.txt", quote=F, sep="\t", row.names = F, col.names = T)
#write.table(reg_tumor_down, "data_tables/tumor_down.txt", quote=F, sep="\t", row.names = F, col.names = T)



#fig7h
ggplot(subset(exp, gene_biotype=="protein_coding" & !is.na(reg_tumor) & !is.na(localization_cat)), aes(log2FoldChange.tumor.KO.WT, colour=localization_cat))+stat_ecdf()+
  scale_color_manual(values=c("dodgerblue3", "orange3"))+coord_cartesian(xlim=c(-2.5,2.5))+
  geom_hline(yintercept=0.5, lty=2, colour="grey")

sd<-subset(exp, gene_biotype=="protein_coding" & !is.na(reg_tumor) & !is.na(localization_cat) & !is.na(log2FoldChange.tumor.KO.WT),
           select=c("gene_id", "Symbol", "log2FoldChange.tumor.KO.WT", "reg_tumor", "localization_cat"))
# write.table(sd, "source_data/fig7h.txt", quote=F, sep="\t", row.names = F)

nrow(subset(exp, gene_biotype=="protein_coding" & 
              !is.na(reg_tumor) & localization_cat=="cytosolic" & 
              !is.na(log2FoldChange.tumor.KO.WT)))
nrow(subset(exp, gene_biotype=="protein_coding" & 
              !is.na(reg_tumor) & localization_cat=="membrane" & 
              !is.na(log2FoldChange.tumor.KO.WT)))

g1<-subset(exp, gene_biotype=="protein_coding" & 
             !is.na(reg_tumor) & localization_cat=="cytosolic" & 
             !is.na(log2FoldChange.tumor.KO.WT), select="log2FoldChange.tumor.KO.WT")[,1]
g2<-subset(exp, gene_biotype=="protein_coding" & 
             !is.na(reg_tumor) & localization_cat=="membrane" & 
             !is.na(log2FoldChange.tumor.KO.WT), select="log2FoldChange.tumor.KO.WT")[,1]

#does not run
# wilcox.test(g1[!duplicated(g1)],g2[!duplicated(g2)], exact = T)

wilcox.test(subset(exp, gene_biotype=="protein_coding" & 
                     !is.na(reg_tumor) & localization_cat=="cytosolic" & 
                     !is.na(log2FoldChange.tumor.KO.WT), select="log2FoldChange.tumor.KO.WT")[,1],
            subset(exp, gene_biotype=="protein_coding" & 
                     !is.na(reg_tumor) & localization_cat=="membrane" & 
                     !is.na(log2FoldChange.tumor.KO.WT), select="log2FoldChange.tumor.KO.WT")[,1], 
            exact=F)

nrow(subset(exp, gene_biotype=="protein_coding" & 
              reg_tumor=="down_tumor" & 
              !is.na(log2FoldChange.tumor.KO.WT)))
nrow(subset(exp, gene_biotype=="protein_coding" & 
              reg_tumor=="up_tumor" & 
              !is.na(log2FoldChange.tumor.KO.WT)))
nrow(subset(exp, gene_biotype=="protein_coding" & 
              reg_tumor=="unchanged_tumor" & 
              !is.na(log2FoldChange.tumor.KO.WT)))

nrow(subset(exp, gene_biotype=="protein_coding" & 
              reg_tumor=="not_quantified" & 
              !is.na(log2FoldChange.tumor.KO.WT)))





