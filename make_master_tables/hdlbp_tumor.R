library(ggplot2)
library(reshape2)
library(corrplot)
library(DESeq2)
# setwd("D:/Documents/hdlbp_git/tumor/")
setwd("../../tumor/")

files<-list.files(getwd(), pattern="genes")

dat<-lapply(files, read.delim)
names(dat)<-gsub("\\..*","",files)

tab<-data.frame(sapply(dat, function(x) cbind(x[,5])))
tab$gene_id<-dat[[1]][,1]
tab<-tab[,c(19,1:18)]

corrplot(cor(tab), type="upper", method="color")
# write.table(tab,"data/tumor_rnaseq_read_counts.txt", quote=F, sep="\t", row.names = F)

tpm<-data.frame(sapply(dat, function(x) cbind(x[,6])))
tpm$gene_id<-dat[[1]][,1]
tpm<-tpm[,c(19,1:18)]

corrplot(cor(tpm), type="upper", method="color")
# write.table(tpm,"data/tumor_rnaseq_tpm.txt", quote=F, sep="\t", row.names = F)

#figS7
tab<-read.delim("data/tumor_rnaseq_read_counts.txt", header=T)
tpm<-read.delim("data/tumor_rnaseq_tpm.txt", header=T)

##PCA
# subs<-subset(tab, select=colnames(tab)[grepl("[CMT]_",colnames(tab))])
subs<-tab
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

# genes<-read.table("geneInfo.txt", header=T)[,c(1,4:5,10)]




####deseq simple comparisons within fractions
cnt<-tab[,colnames(tab)[grepl("Tumor_", colnames(tab))]]

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


cnt<-tab[,colnames(tab)[grepl("A_", colnames(tab))]]

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


setwd("D:/google_drive//hdlbp/")
mas<-read.delim("hdlbp_master_table_with_classes_uniq_tsig.txt", header=T)
inf<-subset(mas, select=c("gene_id","gene_biotype", 
                          "tsig" , "loc_tar_CDS", "localization_cat" ))
setwd("E:/work/hdlbp/tumor")
# genes<-read.table("geneInfo.txt", header=T)[,c(1,4:5,10)]

exp<-merge(fin, inf, by="gene_id", all.x=T)
# exp<-merge(exp, genes, by.x="gene_id", by.y="Gene", all.x=T)

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



# reg_tumor_up<-subset(exp, reg_tumor=="up_tumor" & Annotation=="protein_coding", 
#                      select = colnames(exp)[c(1,68,59:61,62:67,69:71,8:19)])
# reg_tumor_up<-reg_tumor_up[order(reg_tumor_up$log2FoldChange.tumor.KO.WT, decreasing = T),]
# reg_tumor_down<-subset(exp, reg_tumor=="down_tumor" & Annotation=="protein_coding", 
#                      select = colnames(exp)[c(1,68,59:61,62:67,69:71,8:19)])
# reg_tumor_down<-reg_tumor_down[order(reg_tumor_down$log2FoldChange.tumor.KO.WT, decreasing = F),]

# reg_tumor_unchanged<-subset(exp, reg_tumor=="unchanged_tumor", 
#                        select = colnames(exp)[c(1,62,59:61,63:69,8:19)])
# 
# 
# reg_a549_up<-subset(exp, reg_a549=="up_a549", 
#                      select = colnames(exp)[c(1,68,59:61,63:67,8:19)])
# reg_a549_down<-subset(exp, reg_a549=="down_a549", 
#                        select = colnames(exp)[c(1,62,59:61,63:69,8:19)])
# reg_a549_unchanged<-subset(exp, reg_a549=="unchanged_a549", 
#                             select = colnames(exp)[c(1,62,59:61,63:69,8:19)])
write.table(reg_tumor_up, "data_tables/tumor_up.txt", quote=F, sep="\t", row.names = F, col.names = T)
write.table(reg_tumor_down, "data_tables/tumor_down.txt", quote=F, sep="\t", row.names = F, col.names = T)
# write.table(reg_tumor_unchanged, "data_tables/tumor_unchanged.txt", quote=F, sep="\t", row.names = F, col.names = T)
# write.table(reg_a549_up, "data_tables/a549_up.txt", quote=F, sep="\t", row.names = F, col.names = T)
# write.table(reg_a549_down, "data_tables/a549_down.txt", quote=F, sep="\t", row.names = F, col.names = T)
# write.table(reg_a549_unchanged, "data_tables/a549_unchanged.txt", quote=F, sep="\t", row.names = F, col.names = T)

ggplot(subset(exp, Annotation=="protein_coding" & !is.na(reg_tumor) & !is.na(localization_cat)), aes(log2FoldChange.tumor.KO.WT, colour=localization_cat))+stat_ecdf()+coord_cartesian(xlim=c(-5,5))+
  scale_color_manual(values=c("dodgerblue3", "orange3"))+coord_cartesian(xlim=c(-2.5,2.5))+
  geom_hline(yintercept=0.5, lty=2, colour="grey")
ggplot(subset(exp, Annotation=="protein_coding" & !is.na(localization_cat)), aes(log2FoldChange.tumor.KO.WT, colour=localization_cat))+stat_ecdf()+coord_cartesian(xlim=c(-5,5))
ggplot(subset(exp, !is.na(reg_tumor) & !is.na(localization_cat)), aes(log2FoldChange.tumor.KO.WT, colour=localization_cat))+stat_ecdf()+coord_cartesian(xlim=c(-5,5))
ggplot(subset(exp, !is.na(localization_cat)), aes(log2FoldChange.tumor.KO.WT, colour=localization_cat))+stat_ecdf()+coord_cartesian(xlim=c(-5,5))

nrow(subset(exp, Annotation=="protein_coding" & 
              !is.na(reg_tumor) & localization_cat=="cytosolic" & 
              !is.na(log2FoldChange.tumor.KO.WT)))
nrow(subset(exp, Annotation=="protein_coding" & 
              !is.na(reg_tumor) & localization_cat=="membrane" & 
              !is.na(log2FoldChange.tumor.KO.WT)))

wilcox.test(subset(exp, gene_biotype=="protein_coding" & 
                     !is.na(reg_tumor) & localization_cat=="cytosolic" & 
                     !is.na(log2FoldChange.tumor.KO.WT), select="log2FoldChange.tumor.KO.WT")[,1],
            subset(exp, gene_biotype=="protein_coding" & 
                     !is.na(reg_tumor) & localization_cat=="membrane" & 
                     !is.na(log2FoldChange.tumor.KO.WT), select="log2FoldChange.tumor.KO.WT")[,1], 
            exact=F)

nrow(subset(exp, Annotation=="protein_coding" & 
              reg_tumor=="down_tumor" & 
              !is.na(log2FoldChange.tumor.KO.WT)))
nrow(subset(exp, Annotation=="protein_coding" & 
              reg_tumor=="up_tumor" & 
              !is.na(log2FoldChange.tumor.KO.WT)))
nrow(subset(exp, Annotation=="protein_coding" & 
              reg_tumor=="unchanged_tumor" & 
              !is.na(log2FoldChange.tumor.KO.WT)))

nrow(subset(exp, Annotation=="protein_coding" & 
              reg_tumor=="not_quantified" & 
              !is.na(log2FoldChange.tumor.KO.WT)))

ggplot(subset(exp, localization_cat=="membrane" &Annotation=="protein_coding" & !is.na(reg_tumor) & !is.na(loc_tar_CDS)), aes(log2FoldChange.tumor.KO.WT, colour=loc_tar_CDS))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))+
  scale_color_manual(values=c("darkgreen","dodgerblue3","orange3","black"))+coord_cartesian(xlim=c(-2,1.5))+
  geom_hline(yintercept=0.5, lty=2, colour="grey")
ggplot(subset(exp, Annotation=="protein_coding" & !is.na(loc_tar_CDS)), aes(log2FoldChange.tumor.KO.WT, colour=loc_tar_CDS))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(subset(exp, !is.na(reg_tumor) & !is.na(loc_tar_CDS)), aes(log2FoldChange.tumor.KO.WT, colour=loc_tar_CDS))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(subset(exp, !is.na(loc_tar_CDS)), aes(log2FoldChange.tumor.KO.WT, colour=loc_tar_CDS))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))

nrow(subset(exp, Annotation=="protein_coding" & 
              loc_tar_CDS=="membrane_tc<1.66"&
              !is.na(reg_tumor) & localization_cat=="membrane" & 
              !is.na(log2FoldChange.tumor.KO.WT)))
nrow(subset(exp, Annotation=="protein_coding" & 
              loc_tar_CDS=="membrane_tc>1.66 & tc<5.56"&
              !is.na(reg_tumor) & localization_cat=="membrane" & 
              !is.na(log2FoldChange.tumor.KO.WT)))
nrow(subset(exp, Annotation=="protein_coding" & 
              loc_tar_CDS=="membrane_tc>5.56 & tc<65.26"&
              !is.na(reg_tumor) & localization_cat=="membrane" & 
              !is.na(log2FoldChange.tumor.KO.WT)))
nrow(subset(exp, Annotation=="protein_coding" & 
              loc_tar_CDS=="nontarget_membrane"&
              !is.na(reg_tumor) & localization_cat=="membrane" & 
              !is.na(log2FoldChange.tumor.KO.WT)))

wilcox.test(subset(exp, Annotation=="protein_coding" & 
                     loc_tar_CDS=="nontarget_membrane"&
                     !is.na(reg_tumor) & localization_cat=="membrane" & 
                     !is.na(log2FoldChange.tumor.KO.WT), select="log2FoldChange.tumor.KO.WT")[,1],
            subset(exp, Annotation=="protein_coding" & 
                     loc_tar_CDS=="membrane_tc>5.56 & tc<65.26"&
                     !is.na(reg_tumor) & localization_cat=="membrane" & 
                     !is.na(log2FoldChange.tumor.KO.WT), select="log2FoldChange.tumor.KO.WT")[,1])            

wilcox.test(subset(exp, Annotation=="protein_coding" & 
                     loc_tar_CDS=="membrane_tc<1.66"&
                     !is.na(reg_tumor) & localization_cat=="membrane" & 
                     !is.na(log2FoldChange.tumor.KO.WT), select="log2FoldChange.tumor.KO.WT")[,1],
            subset(exp, Annotation=="protein_coding" & 
                     loc_tar_CDS=="membrane_tc>5.56 & tc<65.26"&
                     !is.na(reg_tumor) & localization_cat=="membrane" & 
                     !is.na(log2FoldChange.tumor.KO.WT), select="log2FoldChange.tumor.KO.WT")[,1])            

wilcox.test(subset(exp, Annotation=="protein_coding" & 
                     loc_tar_CDS=="membrane_tc>1.66 & tc<5.56"&
                     !is.na(reg_tumor) & localization_cat=="membrane" & 
                     !is.na(log2FoldChange.tumor.KO.WT), select="log2FoldChange.tumor.KO.WT")[,1],
            subset(exp, Annotation=="protein_coding" & 
                     loc_tar_CDS=="membrane_tc>5.56 & tc<65.26"&
                     !is.na(reg_tumor) & localization_cat=="membrane" & 
                     !is.na(log2FoldChange.tumor.KO.WT), select="log2FoldChange.tumor.KO.WT")[,1])            
wilcox.test(subset(exp, Annotation=="protein_coding" & 
                     loc_tar_CDS=="membrane_tc<1.66"&
                     !is.na(reg_tumor) & localization_cat=="membrane" & 
                     !is.na(log2FoldChange.tumor.KO.WT), select="log2FoldChange.tumor.KO.WT")[,1],
            subset(exp, Annotation=="protein_coding" & 
                     loc_tar_CDS=="nontarget_membrane"&
                     !is.na(reg_tumor) & localization_cat=="membrane" & 
                     !is.na(log2FoldChange.tumor.KO.WT), select="log2FoldChange.tumor.KO.WT")[,1])            

wilcox.test(subset(exp, Annotation=="protein_coding" & 
                     loc_tar_CDS=="membrane_tc>1.66 & tc<5.56"&
                     !is.na(reg_tumor) & localization_cat=="membrane" & 
                     !is.na(log2FoldChange.tumor.KO.WT), select="log2FoldChange.tumor.KO.WT")[,1],
            subset(exp, Annotation=="protein_coding" & 
                     loc_tar_CDS=="nontarget_membrane"&
                     !is.na(reg_tumor) & localization_cat=="membrane" & 
                     !is.na(log2FoldChange.tumor.KO.WT), select="log2FoldChange.tumor.KO.WT")[,1])            


ggplot(subset(exp, Annotation=="protein_coding" ), 
       aes(log10(baseMean.tumor.KO.WT), log2FoldChange.tumor.KO.WT, color=reg_tumor))+
  geom_point(shape=1, size=1)+coord_cartesian(ylim=c(-10,10), xlim=c(0,6))+scale_color_manual(values=c("orange3", "grey",  "black", "dodgerblue3"))+geom_hline(yintercept=0, lty=2, colour="grey")
sub<-subset(exp, Annotation=="protein_coding")
plot(log10(sub$baseMean.tumor.KO.WT), sub$log2FoldChange.tumor.KO.WT, cex=0.6)
identify(log10(sub$baseMean.tumor.KO.WT), sub$log2FoldChange.tumor.KO.WT, labels=sub$Symbol)


ggplot(subset(exp), 
       aes(log10(baseMean.tumor.KO.WT), log2FoldChange.tumor.KO.WT, color=reg_tumor))+
  geom_point(shape=1, size=1)+coord_cartesian(ylim=c(-10,10), xlim=c(-1,6))




exp$reg_tumor_p<-ifelse(exp$padj.tumor.KO.WT<0.1 , "regulated_tumor",
                             "unchanged_tumor")
ggplot(subset(exp), 
       aes(log10(baseMean.tumor.KO.WT), log2FoldChange.tumor.KO.WT, color=reg_tumor_p))+
  geom_point(shape=1, size=1)+coord_cartesian(ylim=c(-10,10), xlim=c())

ggplot(subset(exp, Annotation=="protein_coding"), 
       aes(log10(baseMean.tumor.KO.WT), log2FoldChange.tumor.KO.WT, color=reg_tumor_p))+
  geom_point(shape=1, size=1)+coord_cartesian(ylim=c(-10,10), xlim=c())

exp$reg_tumor_pc<-ifelse(exp$log2FoldChange.tumor.KO.WT>=1 &
                        exp$padj.tumor.KO.WT<0.1 & exp$Annotation=="protein_coding", "up_tumor",
                      ifelse(exp$log2FoldChange.tumor.KO.WT<=(-1) &
                               exp$padj.tumor.KO.WT<0.1& exp$Annotation=="protein_coding", "down_tumor",
                             ifelse(exp$log2FoldChange.tumor.KO.WT<1 & 
                                      exp$log2FoldChange.tumor.KO.WT>(-1) &
                                      exp$padj.tumor.KO.WT>=0.1& exp$Annotation=="protein_coding", "unchanged_tumor", "not_mRNA")))
exp$reg_tumor_pc<-factor(exp$reg_tumor_pc, levels=(c("up_tumor", "down_tumor", "unchanged_tumor", "not_mRNA")))
ggplot(subset(exp), 
       aes(log10(baseMean.tumor.KO.WT), log2FoldChange.tumor.KO.WT, color=reg_tumor_pc))+
  geom_point(shape=1, size=1, alpha=0.3)+coord_cartesian(ylim=c(-10,10), xlim=c())+
  scale_color_manual(values=c("orange3", "dodgerblue3",  "black", "grey"))
ggplot(subset(exp, Annotation=="protein_coding"), 
       aes(log10(baseMean.tumor.KO.WT), log2FoldChange.tumor.KO.WT, color=reg_tumor_pc))+
  geom_point(shape=1, size=1, alpha=0.3)+coord_cartesian(ylim=c(-10,10), xlim=c())+
  scale_color_manual(values=c("orange3", "dodgerblue3",  "black", "grey"))

#go plot
gos<-read.delim("data_tables/go_table.txt", header=T)[1:12,]
gos$Functional.Category<-factor(gos$Functional.Category, 
                                levels=rev(gos$Functional.Category))
ggplot(gos, aes(Functional.Category,-log10(Enrichment.FDR), fill=regulation))+
  geom_bar(stat="identity")+coord_flip()+
  scale_fill_manual(values=c("orange3", "dodgerblue4"))
