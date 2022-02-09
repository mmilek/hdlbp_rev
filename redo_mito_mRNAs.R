library(ggplot2)


dat<-read.delim("data/hdlbp_master_table_with_classes_uniq_tsig.txt", header=T)
library(ggplot2)
ggplot(dat, aes(tsig, log2FoldChange.mem.cyt.293))+geom_violin()+coord_flip()

ggplot(subset(dat, tpm_cutoff>=10),
       aes(log2FoldChange.mem.cyt.293, colour=tc_CDS_norm_cat))+stat_ecdf()+facet_wrap(~tsig)

mito<-read.delim("data/mitocarta2.txt", header=T)
colnames(mito)<-c("gene_id", "mito_new")
dat$ens_gene_id<-gsub("\\..*", "", dat$gene_id)
dat<-merge(dat, mito, by.x="ens_gene_id", by.y="gene_id", all.x=T)
dat$mito_new<-ifelse(is.na(dat$mito_new), "other", dat$mito_new)
ggplot(subset(dat, tpm_cutoff>=10),
       aes(mito_new,log2FoldChange.mem.cyt.293, fill=mito_new))+geom_violin()
dat$tsig_new<-ifelse((dat$mito=="MitoCarta" & 
                        dat$tsig=="mem_notsig") |
                       (dat$mito=="MitoCarta" & 
                          dat$tsig=="MitoCarta") |
                       (dat$mito=="MitoCarta" & 
                          dat$tsig=="cytANDmem_notsig")
                     , "MitoCarta", dat$tsig)
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
                       ifelse(dat$tsig=="TailAnchored", "TailAnchroed",
                       ifelse(dat$tsig=="cytANDmem_notsig", "cytANDmem_notsig",
                       ifelse(dat$tsig=="MitoEncoded", "MitoEncoded",
                      ifelse(dat$tsig=="SignalP-noTM-only",
                             "SignalP-noTM-only", "other")))))))))

ggplot(subset(dat, tpm_cutoff>=10),
       aes(tsig_new,log2FoldChange.mem.cyt.293, fill=tsig_new))+geom_violin()

nrow(subset(dat, mito=="MitoCarta"))
nrow(subset(dat, tsig_new=="MitoCarta"))
nrow(subset(dat, mito_new=="MitoCarta" & tsig=="cytANDmem_notsig"))
nrow(subset(dat, mito_new=="MitoCarta" & tsig=="MitoCarta"))
nrow(subset(dat, mito_new=="MitoCarta" & tsig=="mem_notsig"))

subset(dat, tsig_new=="MitoCarta" & log2FoldChange.mem.cyt.293>1.5, select="gene_name")
colnames(dat)

mito<-read.delim("data/mitocarta2.txt", header=T)
dat$ens_gene_id<-gsub("\\..*", "", dat$gene_id)
dat<-merge(dat, mito, by.x="ens_gene_id", by.y="gene_id", all.x=T)
ggplot(subset(dat, tpm_cutoff>=10),
       aes(tsig.y,log2FoldChange.mem.cyt.293, fill=tsig.y))+geom_violin()
nrow(subset(dat, tpm_cutoff>=10 & log2FoldChange.mem.cyt.293>1.5 & tsig.y=="MitoCarta"))
nrow(subset(dat, tpm_cutoff>=10 & log2FoldChange.mem.cyt.293<0.5 & tsig.y=="MitoCarta"))
nrow(subset(dat, tpm_cutoff>=10 & tsig.y=="MitoCarta"))

