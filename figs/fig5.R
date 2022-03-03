library(ggplot2)
library(reshape2)
library(zoo)

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

#fig5i begin

mas<-read.delim("data/hdlbp_master_table_with_classes_uniq.txt", header=T)
inf<-subset(mas, select=c("gene_id","Symbol",
                          "tpm_cutoff","tc_CDS_norm","localization_cat",
                          "mean_te_293","loc_tar_CDS","tc_CDS_norm_cat",
                          "tc_transcript_norm",
                          "log2FoldChange.mem.cyt.293",
                          "log2FoldChange.mem.cyt.KO.293",
                          "log2FoldChange.ribo.rna.KO.WT",
                          "Transmembrane.helices" ,"Cleavage.site..Signalp." ))

trans<-read.delim("data/considered_gene_names.txt", header=T)[,c(1:10,16,18)]
colnames(trans)[10]<-"gene_id"

inf<-merge(inf, trans, by="gene_id")

fin<-read.delim("data/psite_position_counts_chx.txt", header=T)
colnames(fin)[2:3]<-c("293_1", "293_2")

tsig<-read.delim("data/signalp_tm_positions.txt", header=T) # download this from http://grch37.ensembl.org/biomart/martview careful with id versions
genc<-read.delim("data/gencode_v19_gene_id_to_gene_name_all.txt", header=F)
nrow(merge(tsig, genc, by.x="Gene.stable.ID.version", by.y="V1"))

tsig$Transmembrane.helices<-as.character(tsig$Transmembrane.helices)
tsig$Cleavage.site..Signalp.<-as.character(tsig$Cleavage.site..Signalp.)

tsig$tsig<-ifelse(tsig$Transmembrane.helices=="TMhelix" & tsig$Cleavage.site..Signalp.=="SignalP-TM", "SignalP-TM-TMhelix",
                  ifelse(tsig$Transmembrane.helices!="TMhelix" & tsig$Cleavage.site..Signalp.=="SignalP-TM", "SignalP-TM-only",
                         ifelse(tsig$Transmembrane.helices!="TMhelix" & tsig$Cleavage.site..Signalp.=="SignalP-noTM", "SignalP-noTM-only",
                                ifelse(tsig$Transmembrane.helices=="TMhelix" & tsig$Cleavage.site..Signalp.=="SignalP-noTM", "SignalP-noTM-TM",
                                       ifelse(tsig$Transmembrane.helices=="TMhelix" & tsig$Cleavage.site..Signalp.!="SignalP-noTM" & tsig$Cleavage.site..Signalp.!="SignalP-TM", "TMhelix-only","none")))))

nrow(subset(tsig, Transmembrane.helices=="TMhelix" & Cleavage.site..Signalp.=="SignalP-TM"))
nrow(subset(tsig, Transmembrane.helices!="TMhelix" & Cleavage.site..Signalp.=="SignalP-TM"))
nrow(subset(tsig, Transmembrane.helices!="TMhelix" & Cleavage.site..Signalp.=="SignalP-noTM"))
nrow(subset(tsig, Transmembrane.helices=="TMhelix" & Cleavage.site..Signalp.=="SignalP-noTM"))

ggplot(tsig, aes(tsig))+geom_bar()
ggplot(tsig, aes(Cleavage.site..Signalp..start))+geom_histogram()
ggplot(tsig, aes(Cleavage.site..Signalp..end))+geom_histogram()
ggplot(tsig, aes(Transmembrane.helices.start))+geom_histogram(bins=250)
ggplot(tsig, aes(Transmembrane.helices.end))+geom_histogram(bins=250)

tsig<-subset(tsig, tsig!="none")

tsig$start<-ifelse(grepl("SignalP", tsig$tsig), tsig$Cleavage.site..Signalp..start, tsig$Transmembrane.helices.start)

dms<-aggregate(start~ Transcript.stable.ID.version+Gene.stable.ID.version+tsig, data=tsig, min)

ggplot(dms, aes(start))+geom_histogram(bins=1000)+facet_wrap(~tsig)+coord_cartesian(xlim=c(-10,100))

glock<-which(dms$Transcript.stable.ID.version %in% inf$transcript)

tarsig<-dms[glock,]

ggplot(tarsig, aes(tsig))+geom_bar()

tsig_annot<-merge(inf, tarsig, by.x="gene_id", by.y="Gene.stable.ID.version", all.x=T)
tsig_annot$tsig<-ifelse(is.na(tsig_annot$tsig) & tsig_annot$localization_cat=="membrane", "mem_notsig", 
                        ifelse(is.na(tsig_annot$tsig) & tsig_annot$localization_cat=="cytosolic", "cyt_notsig",tsig_annot$tsig))
ggplot(tsig_annot, aes(tsig))+geom_bar()

nrow(subset(tsig_annot, as.character(transcript)==as.character(Transcript.stable.ID.version)))

man<-subset(tsig_annot, !is.na(tsig))

man$nucstart<-(man$start-1)*3
man$minusstart<-ifelse(man$nucstart==0,0, -man$nucstart)
man$plusstart<-ifelse(man$l_cds>=1502, 1502+man$minusstart,man$l_cds+man$minusstart )

man$dist_stop<-man$l_cds-man$nucstart

man$tsig<-ifelse(man$dist_stop<=150 & man$tsig=="TMhelix-only" , "TailAnchored", man$tsig)
man$tsig<-ifelse(grepl("MT-", man$Symbol), "MitoEncoded", man$tsig)
man$tsig<-ifelse(grepl("SignalP", man$tsig), "SignalP", man$tsig)

mito<-read.delim("data/mitocarta2.txt", header=T)
colnames(mito)<-c("gene_id", "mito")
man$ens_gene_id<-gsub("\\..*","",man$gene_id)
stud<-merge(man, mito, by.x="ens_gene_id", by.y="gene_id", all.x=T)
stud$tsig<-ifelse(is.na(stud$mito), stud$tsig,
                  ifelse(stud$tsig=="cyt_notsig" & stud$mito=="MitoCarta", "MitoCarta", stud$tsig))
man<-stud

ggplot(man, aes(tsig))+geom_bar()
ggplot(man[man$tpm_cutoff>=10,],aes(tsig,log2(tc_CDS_norm),fill=tsig))+geom_boxplot()+geom_point(aes(fill = tsig), size = 1, shape = 1, position = position_jitterdodge())
ggplot(man[man$tpm_cutoff>=10,],aes(tsig,log2(tc_CDS_norm),fill=tsig))+geom_violin(scale="area",na.rm=T)+geom_boxplot(fill="white", na.rm=T, width=0.3)

ggplot(man[man$tpm_cutoff>=10,],aes(tsig,log2(tc_CDS_norm),fill=tsig))+geom_boxplot()+geom_point(aes(fill = tsig), size = 1, shape = 1, position = position_jitterdodge())

ggplot(man[man$tpm_cutoff>=10,],aes(log2FoldChange.mem.cyt.293,log2(mean_te_293),colour=tsig))+geom_point()


man<-subset(man, tsig!="TailAnchored")
man<-subset(man, !is.na(minusstart) & !is.na(plusstart))
seqs<-man[,c("transcript", "minusstart","plusstart")]
seqs<-subset(seqs, !is.na(minusstart) & !is.na(plusstart))
lseq<-apply(seqs, 1, function(x) seq(x[2], x[3]))

dseq<-data.frame(transcript=rep(unique(man$transcript), man$plusstart-man$minusstart+1), pos_from_tsig=unlist(lseq))
dseq<-subset(dseq, pos_from_tsig>=(-300))

dseq$psite_id_tsig<-paste0(dseq$transcript,"_",dseq$pos_from_tsig)

df<-merge(fin, tarsig, by.x="transcript", by.y="Transcript.stable.ID.version")

df$pos_from_tsig<-df$pos_from_start-(df$start-1)*3
df$psite_id_tsig<-paste0(df$transcript,"_",df$pos_from_tsig)
df$frame_start<-ifelse(df$pos_from_tsig%%3==0, df$pos_from_tsig/3,
                       ifelse(df$pos_from_tsig%%3==1, ((df$pos_from_tsig)-1)/3, ((df$pos_from_tsig)-2)/3))

tust<-merge(dseq, df, by="psite_id_tsig", all.x=T)
tust$pos_from_tsig<-as.numeric(gsub(".*_","",tust$psite_id_tsig))
tust$frame_start<-ifelse(tust$pos_from_tsig%%3==0, tust$pos_from_tsig/3,
                         ifelse(tust$pos_from_tsig%%3==1, ((tust$pos_from_tsig)-1)/3, ((tust$pos_from_tsig)-2)/3))

tust<-subset(tust, select=c("psite_id_tsig","transcript.x","pos_from_tsig","frame_start","293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2"))
values<-tust[,c("293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2")]
values[is.na(values)]<-0
tust<-cbind(tust[,1:4],values)
tust$id<-paste0(tust$transcript.x,"_",tust$frame_start)

codons<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~id+transcript.x, data=tust, sum)
codons$frame_start<-as.numeric(gsub(".*_","",codons$id))

excluded<-which(codons$frame_start<=40 & codons$frame_start>=20)
excluded<-codons[excluded,]
normi<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~transcript.x, data=excluded, mean, na.rm=T)

normi[normi==0]<-NA

tes<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~id, data=tust, sum)

tes$transcript<-gsub("_.*","",tes$id)
tes$codon<-as.numeric(gsub(".*_","",tes$id))

mask<-merge(tes, normi, by.x="transcript", by.y="transcript.x")
first<-as.matrix(mask[,3:8])
second<-as.matrix(mask[10:15])
nor<-first/second

nor<-cbind(mask[,c(1,2,9)], nor)
colnames(nor)<-gsub("\\.x","",colnames(nor))

avg<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~codon, data=nor, mean, na.rm=T)
mel<-melt(avg, measure.vars = c("293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2"), id.vars = "codon")

ggplot(mel, aes(codon, value, colour=variable))+geom_line()+coord_cartesian(xlim=c(-0,500) )
ggplot(mel[mel$codon>2,], aes(codon, value, colour=variable))+geom_line()+coord_cartesian(xlim=c(0,500 ))


##tm or sp
tm<-subset(tsig_annot, tsig=="TMhelix-only" & tpm_cutoff>=10 & localization_cat=="membrane")
sp<-subset(tsig_annot, grepl("SignalP", tsig_annot$tsig) & tpm_cutoff>=10 & localization_cat=="membrane" & !is.na(loc_tar_CDS) & loc_tar_CDS!="nontarget_membrane")
# tm<-subset(tsig_annot, tsig=="TMhelix-only" & tpm_cutoff>=10 )
# sp<-subset(tsig_annot, grepl("SignalP", tsig_annot$tsig) & tpm_cutoff>=10)

tm_list<-subset(tm, select=c("gene_id", "Symbol", "Transcript.stable.ID.version", 
                             "l_utr5", "l_cds", "l_utr3", "start"))
sp_list<-subset(sp, select=c("gene_id", "Symbol", "Transcript.stable.ID.version", 
                             "l_utr5", "l_cds", "l_utr3", "start"))

# make lists of genes for browsing
# write.table(tm_list, "data/tm_list.txt", quote=F, sep="\t", row.names=F)
# write.table(sp_list, "data/sp_list.txt", quote=F, sep="\t", row.names=F)

tm<-which(nor$transcript %in% tm$transcript)
sp<-which(nor$transcript %in% sp$transcript)

tm<-nor[tm,]
sp<-nor[sp,]

length(unique(tm$transcript))
length(unique(sp$transcript))

list1<-data.frame(transcript=unique(sp$transcript), tsig="sp")
list2<-data.frame(transcript=unique(tm$transcript), tsig="tm")
list<-rbind(list1, list2)

# write.table(list, "tm_sp_ribo_tsig_list.txt", quote=F, sep="\t", row.names=F)

tm$localization<-"tm"
sp$localization<-"sp"

locn<-rbind(tm, sp)

avg<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~codon+localization, data=locn, mean, na.rm=T)
mel<-melt(avg, measure.vars = c("293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2"), id.vars = c("codon","localization"))

ggplot(mel, aes(codon, value, colour=localization))+geom_line(aes(y=rollmean(value, 5, na.pad=T)))+facet_wrap(~variable)


ggplot(mel[ mel$localization=="sp",], aes(codon, value, colour=localization))+geom_line(aes(y=rollmean(value, 5, na.pad=T)))+facet_wrap(~variable)+
        geom_hline(yintercept = 2.5, lty=2)+geom_vline(xintercept = 40, lty=2)+coord_cartesian(ylim=c(0,5),xlim=c(0,250))
ggplot(mel[ mel$localization=="tm",], aes(codon, value, colour=localization))+geom_line(aes(y=rollmean(value, 5, na.pad=T)))+facet_wrap(~variable)+
        geom_hline(yintercept = 2, lty=2)+geom_vline(xintercept = 40, lty=2)+coord_cartesian(ylim=c(0,4),xlim=c(-100,250))

# make source data
# write.table(mel, "../../hdlbp_rev/hdlbp_rev/source_data/fig5i.txt", quote=F, sep="\t", row.names = F)




