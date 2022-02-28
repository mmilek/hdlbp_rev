library(ggplot2)
library(reshape2)
library(corrplot)
library(DESeq2)

#figS1b figS1c

dat<-read.delim("geo_processed_data/processed_data_rnaseq_fractionation.txt", header=T)

tab<-subset(dat, select=c("gene_id", colnames(dat)[grepl("^[CMT]_293", colnames(dat))]))

row.names(tab)<-tab$gene_id
tab<-tab[,-1]

#figS3d
corrplot(cor(tab), type="upper", method="color", addCoef.col = "white")

# write.table(cor(tab), "source_data/figS1c.txt", quote=F, sep="\t", row.names=F)

coldata<-data.frame(SampleID=colnames(tab),
                    Fraction=sub(";","_",gsub("_.*","",sub("_",";",colnames(tab)))),
                    Batch=substr(colnames(tab), nchar(colnames(tab)), nchar(colnames(tab))))

coldata<-apply(coldata,2,as.factor)

dds<-DESeqDataSetFromMatrix(countData = round(tab),
                            colData = coldata, 
                            design = ~Batch+Fraction)

vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Fraction", "Batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

#figS3b
ggplot(pcaData, aes(PC1, PC2, color=Fraction, shape=Batch)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 

# write.table(pcaData, "source_data/figS1b.txt", quote=F, sep="\t", row.names=F)

#figS1f

dat<-read.delim("data/hdlbp_master_table_with_classes_uniq.txt", header=T)

dat$localization_cat<-as.factor(dat$localization_cat)

ggplot(subset(dat, tpm_cutoff>=1 & !is.na(localization_cat)& Annotation=="protein_coding"),
       aes(log2(tc_CDS), log2(tpm_cutoff), colour=localization_cat))+geom_point(shape=1, alpha=0.7)+
  scale_colour_manual(values=c("dodgerblue4", "orange2"))+xlim(0,15)+ylim(0,15)

ggplot(subset(dat, tpm_cutoff>=1 & localization_cat=="cytosolic" & Annotation=="protein_coding"),
       aes(log2(tc_CDS), log2(tpm_cutoff), colour=localization_cat))+geom_point(shape=1, alpha=0.7)+
  scale_colour_manual(values=c("dodgerblue4"))+xlim(0,15)+ylim(0,15)

ggplot(subset(dat, tpm_cutoff>=1 & localization_cat=="membrane" & Annotation=="protein_coding"),
       aes(log2(tc_CDS), log2(tpm_cutoff), colour=localization_cat))+geom_point(shape=1, alpha=0.7)+
  scale_colour_manual(values=c("orange2"))+xlim(0,15)+ylim(0,15)

#figS1g

ggplot(subset(dat, tpm_cutoff>=1 & !is.na(localization_cat)& Annotation=="protein_coding"),
       aes(log2(tc_CDS_norm), log2(tpm_cutoff), colour=localization_cat))+geom_point(shape=1, alpha=0.7)+
  scale_colour_manual(values=c("dodgerblue4", "orange2"))+xlim(-10.5,8)+ylim(0,15)

ggplot(subset(dat, tpm_cutoff>=1 & localization_cat=="cytosolic" & Annotation=="protein_coding"),
       aes(log2(tc_CDS_norm), log2(tpm_cutoff), colour=localization_cat))+geom_point(shape=1, alpha=0.7)+
  scale_colour_manual(values=c("dodgerblue4"))+xlim(-10.5,8)+ylim(0,15)

ggplot(subset(dat, tpm_cutoff>=1 & localization_cat=="membrane" & Annotation=="protein_coding"),
       aes(log2(tc_CDS_norm), log2(tpm_cutoff), colour=localization_cat))+geom_point(shape=1, alpha=0.7)+
  scale_colour_manual(values=c("orange2"))+xlim(-10.5,8)+ylim(0,15)

sd<-subset(dat, tpm_cutoff>=1 & !is.na(localization_cat)& Annotation=="protein_coding" & !is.na(tc_CDS_norm) & !is.na(tc_CDS), 
           select=c("gene_id", "Symbol","tc_CDS", "tc_CDS_norm", "tpm_cutoff", "localization_cat", "Annotation"))
# write.table(sd, "source_data/figS1fg.txt", quote=F, sep="\t", row.names=F)


#figS1h

dat<-read.delim("data/hdlbp_master_table_with_classes_uniq.txt", header=T)

dat$localization_cat<-as.factor(dat$localization_cat)

mel<-subset(dat, !is.na(localization_cat) & Annotation=="protein_coding" & tpm_cutoff>=10 )
mel<-reshape2::melt(mel, measure.vars=c("tc_transcript_norm", "tc_CDS_norm",  "tc_UTR3_norm"),
                    id.vars=c("gene_id","Symbol", "localization_cat"),
                              variable.name="region",
                              value.name="par_clip_enrichment")
mel$region<-sub("_norm","",sub("tc_","",mel$region))
mel$region<-factor(mel$region, levels=c("transcript", "CDS", "UTR3"))
mel<-subset(mel, !is.na(par_clip_enrichment))

ggplot(mel, aes(localization_cat, log2(par_clip_enrichment), fill=localization_cat))+
  geom_violin(scale="count",na.rm=T,position=position_dodge())+
  geom_boxplot(width=0.1,na.rm=T, position=position_dodge(width=0.9), outlier.shape = NA)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+scale_fill_manual(values=c("dodgerblue2","orange3"))+
  facet_wrap(~region)

summary(subset(mel, region=="transcript" & !is.na(par_clip_enrichment))$localization_cat)
summary(subset(mel, region=="CDS" & !is.na(par_clip_enrichment))$localization_cat)
summary(subset(mel, region=="UTR3" & !is.na(par_clip_enrichment))$localization_cat)

# write.table(mel, "source_data/figS1h.txt", quote=F, sep="\t", row.names=F, col.names=T)

#figS1i

ggplot(subset(dat, tpm_cutoff>=10 & Annotation=="protein_coding"),
       aes(log2FoldChange.mem.cyt.293_1,log2FoldChange.mem.cyt.293_2,colour=log2(tc_CDS_norm)))+
       geom_point(shape=1, alpha=0.7)+
       scale_colour_gradient(low="dodgerblue2", high="orange3", limits=c(-5,5))

sd<-subset(dat, tpm_cutoff>=10 & Annotation=="protein_coding" & !is.na(log2FoldChange.mem.cyt.293_1) & !is.na(log2FoldChange.mem.cyt.293_2),
           select=c("gene_id", "Symbol", "tpm_cutoff", "log2FoldChange.mem.cyt.293_1", "log2FoldChange.mem.cyt.293_2", "tc_CDS_norm"))

# write.table(sd, "source_data/figS1i.txt", quote=F, sep="\t", row.names=F, col.names=T)
