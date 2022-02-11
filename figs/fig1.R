library(ggplot2)
library(reshape2)


#fig1a and fig1b

#prepare data, correct for MitoCarta mRNAs

dat<-read.delim("data/hdlbp_master_table_with_classes_uniq_tsig.txt", header=T)
mito<-read.delim("data/mitocarta2.txt", header=T)
colnames(mito)<-c("gene_id", "mito_new")
dat$ens_gene_id<-gsub("\\..*", "", dat$gene_id)
dat<-merge(dat, mito, by.x="ens_gene_id", by.y="gene_id", all.x=T)
dat$mito_new<-ifelse(is.na(dat$mito_new), "other", dat$mito_new)
dat$tsig_new<-ifelse(  (dat$mito_new=="MitoCarta" &
                          dat$tsig=="mem_notsig") |
                         (dat$mito_new=="MitoCarta" &
                            dat$tsig=="MitoCarta") |
                         (dat$mito_new=="MitoCarta" &
                            dat$tsig=="cytANDmem_notsig")
                       , "MitoCarta",
                       ifelse(dat$tsig=="cyt_notsig", "cyt_notsig",
                              ifelse(dat$tsig=="mem_notsig", "mem_notsig",
                                     ifelse(dat$tsig=="TMhelix-only", "TMhelix-only",
                                            ifelse(dat$tsig=="SignalP-TM", "SignalP-TM",
                                                   ifelse(dat$tsig=="TailAnchored", "TailAnchored",
                                                                 ifelse(dat$tsig=="MitoEncoded", "MitoEncoded",
                                                                        ifelse(dat$tsig=="SignalP-noTM-only",
                                                                               "SignalP-noTM-only", "other"))))))))

dat$tsig_new<-factor(dat$tsig_new, levels=c("cyt_notsig", "MitoCarta", "TailAnchored", "MitoEncoded", "mem_notsig",
                                            "SignalP-noTM-only" ,"SignalP-TM", "TMhelix-only", "other"))

dat$localization_cat<-as.factor(dat$localization_cat)

#fig1a

ggplot(subset(dat, tpm_cutoff>=10 & gene_biotype=="protein_coding"), aes(log2FoldChange.mem.cyt.293, fill=localization_cat))+geom_histogram(bins=250)+
  scale_fill_manual(values=c("dodgerblue4", "orange2"))

summary(subset(dat, tpm_cutoff>=10 & gene_biotype=="protein_coding" & !is.na(log2FoldChange.mem.cyt.293))$localization_cat)

sd<-subset(dat, tpm_cutoff>=10 & gene_biotype=="protein_coding" & !is.na(log2FoldChange.mem.cyt.293), select=c("gene_id", "Symbol", "log2FoldChange.mem.cyt.293", "localization_cat"))
# write.table(sd, "source_data/fig1a.txt", quote=F, sep="\t", row.names=F, col.names=T)
ggplot(sd, aes(log2FoldChange.mem.cyt.293, fill=localization_cat))+geom_histogram(bins=250)+
  scale_fill_manual(values=c("dodgerblue4", "orange2"))


#fig1b
ggplot(subset(dat, tpm_cutoff>=10 &  gene_biotype=="protein_coding" & tsig_new!="other"),
       aes(tsig_new,log2FoldChange.mem.cyt.293, fill=tsig_new))+geom_violin()+
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
  scale_fill_manual(values=c("dodgerblue4", "dodgerblue4","dodgerblue4","darkgreen","orange2", "darkred", "darkred","darkred"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  

summary(subset(dat, tpm_cutoff>=10 &  gene_biotype=="protein_coding" & !is.na(log2FoldChange.mem.cyt.293))$tsig_new)

sd<-subset(dat, tpm_cutoff>=10 &  gene_biotype=="protein_coding" & !is.na(log2FoldChange.mem.cyt.293) & tsig_new!="other", select=c("gene_id", "Symbol", "log2FoldChange.mem.cyt.293", "localization_cat", "tsig_new"))
# write.table(sd, "source_data/fig1b.txt", quote=F, sep="\t", row.names=F, col.names=T)
ggplot(sd,
       aes(tsig_new,log2FoldChange.mem.cyt.293, fill=tsig_new))+geom_violin()+
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
  scale_fill_manual(values=c("dodgerblue4", "dodgerblue4","dodgerblue4","darkgreen","orange2", "darkred", "darkred","darkred"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


#fig1d
ggplot(subset(dat, !is.na(localization_cat)& gene_biotype=="protein_coding"),
       aes(log2(tc_CDS_norm), log2(tpm_cutoff), colour=localization_cat))+geom_point(shape=1, alpha=0.7)+
       scale_colour_manual(values=c("dodgerblue4", "orange2"))


ggplot(subset(dat, !is.na(localization_cat)& gene_biotype=="protein_coding"),
       aes(log2(tc_CDS_norm), colour=localization_cat))+geom_density()+
       scale_colour_manual(values=c("dodgerblue4", "orange2"))


sd<-subset(dat, !is.na(localization_cat)& gene_biotype=="protein_coding" & !is.na(tc_CDS_norm) & !is.na(tpm_cutoff) , select=c("gene_id", "Symbol", "tc_CDS_norm", "tpm_cutoff", "localization_cat", "tsig_new"))
# write.table(sd, "source_data/fig1d.txt", quote=F, sep="\t", row.names=F, col.names=T)
ggplot(sd,
       aes(log2(tc_CDS_norm), log2(tpm_cutoff), colour=localization_cat))+geom_point(shape=1, alpha=0.7)+
  scale_colour_manual(values=c("dodgerblue4", "orange2"))

#fig1e

nrow(subset(dat, tpm_cutoff>=10 & gene_biotype=="protein_coding" & !is.na(tc_transcript_norm) & localization_cat=="membrane"))
nrow(subset(dat, tpm_cutoff>=10 & gene_biotype=="protein_coding" & !is.na(tc_transcript_norm) & localization_cat=="cytosolic"))

nrow(subset(dat, tpm_cutoff>=10  & gene_biotype=="protein_coding"  & is.na(tc_transcript_norm) & localization_cat=="membrane"))
nrow(subset(dat, tpm_cutoff>=10  & gene_biotype=="protein_coding"  & is.na(tc_transcript_norm) & localization_cat=="cytosolic"))


#fig1f
ggplot(subset(dat, tpm_cutoff>=10 & tsig_new!="other"& tsig_new!="MitoEncoded"),
       aes(tsig_new,log2(tc_CDS_norm), fill=tsig_new))+geom_violin()+
  geom_boxplot(width=0.2, fill="white", outlier.shape = NA)+
  scale_fill_manual(values=c("dodgerblue4", "dodgerblue4","dodgerblue4","orange2", "darkred", "darkred","darkred"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

summary(subset(dat, tpm_cutoff>=10 & tsig_new!="other"& tsig_new!="MitoEncoded" & !is.na(tc_CDS_norm))$tsig_new)




