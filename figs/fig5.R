library(ggplot2)
library(reshape2)


#fig5b

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



dat$loc_tar_CDS<-as.factor(dat$loc_tar_CDS)
ggplot(subset(dat, localization_cat=="membrane" & !is.na(loc_tar_CDS)),
       aes(log2FoldChange.ribo.rna.KO.WT, colour=loc_tar_CDS))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))

summary(subset(dat, !is.na(loc_tar_CDS) & tpm_cutoff>=10 &localization_cat=="membrane" & !is.na(log2FoldChange.ribo.rna.KO.WT))$loc_tar_CDS)

sd<-subset(dat, localization_cat=="membrane" & !is.na(log2FoldChange.ribo.rna.KO.WT)& !is.na(loc_tar_CDS),
           select=c("gene_id", "Symbol", "tpm_cutoff", "localization_cat",
           "loc_tar_CDS", "log2FoldChange.ribo.rna.KO.WT"))
# write.table(sd, "source_data/fig5b.txt", quote=F, sep="\t", row.names=F)

#start fig5c,fig5d

PG_summ<-read.delim("data/protein_groups_psilac.txt", header=T)
man<-read.delim("data/hdlbp_master_table_with_classes_uniq_tsig.txt", header=T)
pas<-merge(man, PG_summ, by.x="Symbol", by.y="Gene.names", all.x=T)

mel<-reshape2::melt(pas, measure.vars = colnames(pas)[grepl("iBAQ.L", colnames(pas))], id.vars = c("Symbol", "localization_cat","tsig"))
plot<-subset(mel, !is.na(localization_cat))
ggplot(plot, aes(variable, value, fill=localization_cat))+geom_boxplot()+coord_flip()
plot<-subset(mel, !is.na(tsig))
ggplot(plot, aes(variable, value, fill=tsig))+geom_boxplot()+coord_flip()


pas$mean_l_mem<-rowMeans(pas[,colnames(pas)[grepl("iBAQ.L.Memb", colnames(pas))]])

mel<-reshape2::melt(pas, measure.vars = colnames(pas)[grepl("Ratio.H.M.normalized.Cyto", colnames(pas))], id.vars = c("Symbol", "localization_cat","mean_l_mem", "tsig"))
plot<-subset(mel, !is.na(localization_cat))
ggplot(plot, aes(variable, value, fill=localization_cat))+geom_boxplot()+coord_flip()

pas$Ratio.H.M.normalized.Memb_4h_Reverse_1<--pas$Ratio.H.M.normalized.Memb_4h_Reverse_1
pas$Ratio.H.M.normalized.Memb_4h_Reverse_2<--pas$Ratio.H.M.normalized.Memb_4h_Reverse_2
# pas$Ratio.H.M.normalized.Memb_2h_Reverse_1<--pas$Ratio.H.M.normalized.Memb_2h_Reverse_1
# pas$Ratio.H.M.normalized.Memb_2h_Reverse_2<--pas$Ratio.H.M.normalized.Memb_2h_Reverse_2

pas$mean_memb_4h<-rowMeans(pas[,colnames(pas)[grepl("Ratio.H.M.normalized.Memb_4h", colnames(pas))]])

# pas$mean_memb_2h<-rowMeans(pas[,colnames(pas)[grepl("Ratio.H.M.normalized.Memb_2h", colnames(pas))]])

ggplot(subset(pas, !is.na(localization_cat)), aes(mean_memb_4h, colour=localization_cat))+stat_ecdf()+coord_cartesian(xlim=c(-.7,.7))+xlab("H/M KO/WT")
# ggplot(subset(pas, !is.na(localization_cat)), aes(mean_memb_2h, colour=localization_cat))+stat_ecdf()

ks.test(subset(pas, localization_cat=="membrane" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, localization_cat=="cytosolic" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
length(subset(pas,localization_cat=="membrane" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
length(subset(pas,localization_cat=="cytosolic" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])

ggplot(subset(pas, !is.na(tc_transcript_norm_cat)), aes(mean_memb_4h, colour=tc_transcript_norm_cat))+stat_ecdf()+coord_cartesian(xlim=c(-.7,.7))+xlab("H/M KO/WT")

#fig5d
ggplot(subset(pas, !is.na(tc_CDS_norm_cat)), aes(mean_memb_4h, colour=tc_CDS_norm_cat))+stat_ecdf()+coord_cartesian(xlim=c(-.7,.7))+xlab("H/M KO/WT")



ks.test(subset(pas, tc_transcript_norm_cat=="tc>1.12 & tc<132.32" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, tc_transcript_norm_cat=="nontarget" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
ks.test(subset(pas, tc_transcript_norm_cat=="tc>1.12 & tc<132.32" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, tc_transcript_norm_cat=="tc>0.3 & tc<1.12" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
ks.test(subset(pas, tc_transcript_norm_cat=="tc>1.12 & tc<132.32" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, tc_transcript_norm_cat=="tc<0.3" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
ks.test(subset(pas, tc_transcript_norm_cat=="nontarget" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, tc_transcript_norm_cat=="tc>0.3 & tc<1.12" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
ks.test(subset(pas, tc_transcript_norm_cat=="nontarget" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, tc_transcript_norm_cat=="tc<0.3" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
ks.test(subset(pas, tc_transcript_norm_cat=="tc>0.3 & tc<1.12" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, tc_transcript_norm_cat=="tc<0.3" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])

ks.test(subset(pas, tc_CDS_norm_cat=="tc>1.04 & tc<65.26" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, tc_CDS_norm_cat=="nontarget" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
ks.test(subset(pas, tc_CDS_norm_cat=="tc>1.04 & tc<65.26" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, tc_CDS_norm_cat=="tc>0.27 & tc<1.04" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
ks.test(subset(pas, tc_CDS_norm_cat=="tc>1.04 & tc<65.26" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, tc_CDS_norm_cat=="tc<0.27" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
ks.test(subset(pas, tc_CDS_norm_cat=="nontarget" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, tc_CDS_norm_cat=="tc>0.27 & tc<1.04" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
ks.test(subset(pas, tc_CDS_norm_cat=="nontarget" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, tc_CDS_norm_cat=="tc<0.27" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
ks.test(subset(pas, tc_CDS_norm_cat=="tc>0.27 & tc<1.04" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, tc_CDS_norm_cat=="tc<0.27" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])



ks.test(subset(pas, tc_transcript_norm_cat=="tc>1.12 & tc<132.32" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1), select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1],
        subset(pas, tc_transcript_norm_cat=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1), select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1])
ks.test(subset(pas, tc_transcript_norm_cat=="tc>1.12 & tc<132.32" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_2), select="Ratio.H.M.normalized.Memb_4h_Forward_2")[,1],
        subset(pas, tc_transcript_norm_cat=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_2), select="Ratio.H.M.normalized.Memb_4h_Forward_2")[,1])
ks.test(subset(pas, tc_transcript_norm_cat=="tc>1.12 & tc<132.32" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_2), select="Ratio.H.M.normalized.Memb_4h_Reverse_2")[,1],
        subset(pas, tc_transcript_norm_cat=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_2), select="Ratio.H.M.normalized.Memb_4h_Reverse_2")[,1])
ks.test(subset(pas, tc_transcript_norm_cat=="tc>1.12 & tc<132.32" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_1), select="Ratio.H.M.normalized.Memb_4h_Reverse_1")[,1],
        subset(pas, tc_transcript_norm_cat=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_1), select="Ratio.H.M.normalized.Memb_4h_Reverse_1")[,1])

length(subset(pas, tc_transcript_norm_cat=="tc>1.12 & tc<132.32" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
length(subset(pas, tc_transcript_norm_cat=="tc>0.3 & tc<1.12" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
length(subset(pas, tc_transcript_norm_cat=="tc<0.3" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
length(subset(pas, tc_transcript_norm_cat=="nontarget" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])

length(subset(pas, tc_CDS_norm_cat=="tc>1.04 & tc<65.26" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
length(subset(pas, tc_CDS_norm_cat=="tc>0.27 & tc<1.04" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
length(subset(pas, tc_CDS_norm_cat=="tc<0.27" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
length(subset(pas, tc_CDS_norm_cat=="nontarget" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])

ggplot(subset(pas, !is.na(tc_transcript_norm_cat)), aes(Ratio.H.M.normalized.Memb_4h_Forward_1, colour=tc_transcript_norm_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1.5,1.5))
ggplot(subset(pas, !is.na(tc_transcript_norm_cat)), aes(Ratio.H.M.normalized.Memb_4h_Forward_2, colour=tc_transcript_norm_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1.5,1.5))
ggplot(subset(pas, !is.na(tc_transcript_norm_cat)), aes(-Ratio.H.M.normalized.Memb_4h_Reverse_1, colour=tc_transcript_norm_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1.5,1.5))
ggplot(subset(pas, !is.na(tc_transcript_norm_cat)), aes(-Ratio.H.M.normalized.Memb_4h_Reverse_2, colour=tc_transcript_norm_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1.5,1.5))

pas$target<-ifelse(is.na(pas$tc_transcript_norm_cat), NA,
                   ifelse(pas$tc_transcript_norm_cat=="nontarget", "nontarget", "target"))
pas$target_loc<-ifelse(is.na(pas$target), NA,
                       ifelse(pas$target=="nontarget", "nontarget", 
                              ifelse(pas$target=="target" & pas$localization_cat=="cytosolic", "cyto_target",
                                     ifelse(pas$target=="target" & pas$localization_cat=="membrane", "mem_target", NA))))

ggplot(subset(pas, gene_biotype=="protein_coding" & !is.na(target)), aes(mean_memb_4h, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))
ggplot(subset(pas, gene_biotype=="protein_coding" & !is.na(target)), aes(Ratio.H.M.normalized.Memb_4h_Forward_1, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))
ggplot(subset(pas, gene_biotype=="protein_coding" & !is.na(target)), aes(Ratio.H.M.normalized.Memb_4h_Forward_2, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))
ggplot(subset(pas, gene_biotype=="protein_coding" & !is.na(target)), aes(-Ratio.H.M.normalized.Memb_4h_Reverse_1, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))
ggplot(subset(pas, gene_biotype=="protein_coding" & !is.na(target)), aes(-Ratio.H.M.normalized.Memb_4h_Reverse_2, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))

#fig5c
ggplot(subset(pas, gene_biotype=="protein_coding" & !is.na(target_loc)), aes(mean_memb_4h, colour=target_loc))+stat_ecdf()+coord_cartesian(xlim=c(-.7,.7))+xlab("mean H/M (KO/WT)")

sd<-subset(pas, !is.na(tc_CDS_norm_cat) & !is.na(mean_memb_4h),
           select=c("gene_id", "Symbol", "gene_biotype", "tc_CDS_norm_cat", "target_loc", "mean_memb_4h"))
# write.table(sd, "source_data/fig5bc.txt", quote=F, sep="\t", row.names=F)

ggplot(subset(sd, gene_biotype=="protein_coding" & !is.na(target_loc)), aes(mean_memb_4h, colour=target_loc))+stat_ecdf()+coord_cartesian(xlim=c(-.7,.7))+xlab("mean H/M (KO/WT)")
ggplot(subset(sd, !is.na(tc_CDS_norm_cat)), aes(mean_memb_4h, colour=tc_CDS_norm_cat))+stat_ecdf()+coord_cartesian(xlim=c(-.7,.7))+xlab("H/M KO/WT")

ggplot(subset(pas, gene_biotype=="protein_coding" & !is.na(target_loc)), aes(Ratio.H.M.normalized.Memb_4h_Forward_1, colour=target_loc))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))
ggplot(subset(pas, gene_biotype=="protein_coding" & !is.na(target_loc)), aes(Ratio.H.M.normalized.Memb_4h_Forward_2, colour=target_loc))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))
ggplot(subset(pas, gene_biotype=="protein_coding" & !is.na(target_loc)), aes(Ratio.H.M.normalized.Memb_4h_Reverse_1, colour=target_loc))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))
ggplot(subset(pas, gene_biotype=="protein_coding" & !is.na(target_loc)), aes(Ratio.H.M.normalized.Memb_4h_Reverse_2, colour=target_loc))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))


length(subset(pas, gene_biotype=="protein_coding" & !is.na(mean_memb_4h) & target=="target", select="mean_memb_4h")[,1])
length(subset(pas, gene_biotype=="protein_coding" & !is.na(mean_memb_4h) & target=="nontarget", select="mean_memb_4h")[,1])

length(subset(pas, gene_biotype=="protein_coding" & !is.na(mean_memb_4h) & target_loc=="mem_target", select="mean_memb_4h")[,1])
length(subset(pas, gene_biotype=="protein_coding" & !is.na(mean_memb_4h) & target_loc=="cyto_target", select="mean_memb_4h")[,1])
length(subset(pas, gene_biotype=="protein_coding" & !is.na(mean_memb_4h) & target_loc=="nontarget", select="mean_memb_4h")[,1])


ks.test(subset(pas, target=="target" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1), select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1],
        subset(pas, target=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1), select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1])
ks.test(subset(pas, target=="target" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_2), select="Ratio.H.M.normalized.Memb_4h_Forward_2")[,1],
        subset(pas, target=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_2), select="Ratio.H.M.normalized.Memb_4h_Forward_2")[,1])
ks.test(subset(pas, target=="target" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_1), select="Ratio.H.M.normalized.Memb_4h_Reverse_1")[,1],
        subset(pas, target=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_1), select="Ratio.H.M.normalized.Memb_4h_Reverse_1")[,1])
ks.test(subset(pas, target=="target" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_2), select="Ratio.H.M.normalized.Memb_4h_Reverse_2")[,1],
        subset(pas, target=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_2), select="Ratio.H.M.normalized.Memb_4h_Reverse_2")[,1])

ks.test(subset(pas, target_loc=="mem_target" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, target_loc=="cyto_target" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
ks.test(subset(pas, target_loc=="mem_target" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, target_loc=="nontarget" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
ks.test(subset(pas, target_loc=="cyto_target" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, target_loc=="nontarget" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])

ks.test(subset(pas, target_loc=="mem_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1), select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1],
        subset(pas, target_loc=="cyto_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1), select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1])
ks.test(subset(pas, target_loc=="mem_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_2), select="Ratio.H.M.normalized.Memb_4h_Forward_2")[,1],
        subset(pas, target_loc=="cyto_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_2), select="Ratio.H.M.normalized.Memb_4h_Forward_2")[,1])
ks.test(subset(pas, target_loc=="mem_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_1), select="Ratio.H.M.normalized.Memb_4h_Reverse_1")[,1],
        subset(pas, target_loc=="cyto_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_1), select="Ratio.H.M.normalized.Memb_4h_Reverse_1")[,1])
ks.test(subset(pas, target_loc=="mem_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_2), select="Ratio.H.M.normalized.Memb_4h_Reverse_2")[,1],
        subset(pas, target_loc=="cyto_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_2), select="Ratio.H.M.normalized.Memb_4h_Reverse_2")[,1])

ks.test(subset(pas, target_loc=="mem_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1), select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1],
        subset(pas, target_loc=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1), select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1])
ks.test(subset(pas, target_loc=="mem_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_2), select="Ratio.H.M.normalized.Memb_4h_Forward_2")[,1],
        subset(pas, target_loc=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_2), select="Ratio.H.M.normalized.Memb_4h_Forward_2")[,1])
ks.test(subset(pas, target_loc=="mem_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_1), select="Ratio.H.M.normalized.Memb_4h_Reverse_1")[,1],
        subset(pas, target_loc=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_1), select="Ratio.H.M.normalized.Memb_4h_Reverse_1")[,1])
ks.test(subset(pas, target_loc=="mem_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_2), select="Ratio.H.M.normalized.Memb_4h_Reverse_2")[,1],
        subset(pas, target_loc=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_2), select="Ratio.H.M.normalized.Memb_4h_Reverse_2")[,1])

