library(ggplot2)


tc<-read.delim("data/reproducible.hdlbp.TCseq.bed", header=F)
colnames(tc)<-c("transcript_id","tc_start","tc_stop","tc_num1","all_reads1","tc_num2","all_reads2","seq")

tsig_annot<-read.delim("data/tsig_annot.txt", header=T)

dat<-merge(tc, tsig_annot, by="transcript_id", by.y="transcript")
length(unique(dat$transcript_id))
length(unique(dat$gene_id))


norm<-colSums(dat[,c(4,6)], na.rm=T)

dat<-cbind(dat, dat[,c(4,6)]/norm*1e6)
colnames(dat)[36:ncol(dat)]<-paste0("norm.",colnames(dat)[36:ncol(dat)])


dat$cat_tc<-ifelse(dat$tc_stop<=dat$l_utr5, "utr5", 
                   ifelse(dat$tc_stop<=(dat$l_utr5+dat$l_cds), "cds", 
                          ifelse(dat$tc_stop<=(dat$l_utr5+dat$l_cds+dat$l_utr3), "utr3", NA)))
dat<-subset(dat, !is.na(cat_tc))

dat$tc_from_start<-ifelse(dat$cat_tc=="utr5", dat$tc_stop,
                          ifelse(dat$cat_tc=="cds", dat$tc_stop-dat$cds_start,
                                 ifelse(dat$cat_tc=="utr3", dat$tc_stop-dat$cds_stop, NA)))

dat$tc_from_stop<-ifelse(dat$cat_tc=="utr5", dat$tc_stop-dat$l_utr5-1,
                         ifelse(dat$cat_tc=="cds", dat$tc_stop-dat$cds_stop-2,
                                ifelse(dat$cat_tc=="utr3", dat$tc_stop-dat$l_tr-2, NA)))

dat$tc_id<-paste0(dat$transcript_id, "_", dat$tc_from_start,"_",dat$tc_from_stop,"_", dat$cat_tc)

#get counts per region
regionCounts<-aggregate(cbind(norm.tc_num1, norm.tc_num2)~cat_tc+gene_id, data=dat, sum, na.rm=T)
regionCounts<-split(regionCounts, f=regionCounts$cat_tc)
regionCounts<- Reduce(function(...) merge(..., by="gene_id",all=T), regionCounts)
regionCounts<-regionCounts[,-c(2,5,8)]
colnames(regionCounts)[2:ncol(regionCounts)]<-paste0(rep(c("cds.", "utr3.", "utr5."),each=2), c("norm.tc_num1", "norm.tc_num2"))
regionCounts$trans.norm.tc_num1<-rowSums(regionCounts[,c("cds.norm.tc_num1", "utr3.norm.tc_num1", "utr5.norm.tc_num1")], na.rm=T)
regionCounts$trans.norm.tc_num2<-rowSums(regionCounts[,c("cds.norm.tc_num2", "utr3.norm.tc_num2", "utr5.norm.tc_num2")], na.rm=T)
library(corrplot)
regionCounts$utr3_vs_cds1<-regionCounts$utr3.norm.tc_num1/regionCounts$cds.norm.tc_num1
regionCounts$utr3_vs_cds2<-regionCounts$utr3.norm.tc_num2/regionCounts$cds.norm.tc_num2



# #with imputation
# regionCounts$utr3_vs_cds1<-ifelse(!is.na(regionCounts$utr3.norm.tc_num1) & is.na(regionCounts$cds.norm.tc_num1), 
#                                   regionCounts$utr3.norm.tc_num1/min(regionCounts[!is.na(regionCounts$cds.norm.tc_num1), "cds.norm.tc_num1"]),
#                                   ifelse(is.na(regionCounts$utr3.norm.tc_num1) & !is.na(regionCounts$cds.norm.tc_num1), 
#                                          min(regionCounts[!is.na(regionCounts$utr3.norm.tc_num1), "utr3.norm.tc_num1"])/regionCounts$cds.norm.tc_num1,
#                                          ifelse(!is.na(regionCounts$utr3.norm.tc_num1) & !is.na(regionCounts$cds.norm.tc_num1),
#                                                 regionCounts$utr3.norm.tc_num1/regionCounts$cds.norm.tc_num1, NA)))
# regionCounts$utr3_vs_cds2<-ifelse(!is.na(regionCounts$utr3.norm.tc_num2) & is.na(regionCounts$cds.norm.tc_num2), 
#                                   regionCounts$utr3.norm.tc_num2/min(regionCounts[!is.na(regionCounts$cds.norm.tc_num2), "cds.norm.tc_num2"]),
#                                   ifelse(is.na(regionCounts$utr3.norm.tc_num2) & !is.na(regionCounts$cds.norm.tc_num2), 
#                                          min(regionCounts[!is.na(regionCounts$utr3.norm.tc_num2), "utr3.norm.tc_num2"])/regionCounts$cds.norm.tc_num2,
#                                          ifelse(!is.na(regionCounts$utr3.norm.tc_num2) & !is.na(regionCounts$cds.norm.tc_num2),
#                                                 regionCounts$utr3.norm.tc_num2/regionCounts$cds.norm.tc_num2, NA)))
# 
# corrplot(cor(regionCounts[2:ncol(regionCounts)], use = "pairwise.complete.obs"), type="upper", method="color")
# 
# ggplot(regionCounts, aes(log2(utr3_vs_cds1), log2(utr3_vs_cds2)))+geom_point()


regionCounts<-merge(regionCounts, tsig_annot, by="gene_id")
regionCounts$cds.norm.tc_num1.len<-regionCounts$cds.norm.tc_num1/regionCounts$l_cds*1e3
regionCounts$cds.norm.tc_num2.len<-regionCounts$cds.norm.tc_num2/regionCounts$l_cds*1e3
regionCounts$utr3.norm.tc_num1.len<-regionCounts$utr3.norm.tc_num1/regionCounts$l_utr3*1e3
regionCounts$utr3.norm.tc_num2.len<-regionCounts$utr3.norm.tc_num2/regionCounts$l_utr3*1e3
regionCounts$utr5.norm.tc_num1.len<-regionCounts$utr5.norm.tc_num1/regionCounts$l_utr5*1e3
regionCounts$utr5.norm.tc_num2.len<-regionCounts$utr5.norm.tc_num2/regionCounts$l_utr5*1e3
regionCounts$trans.norm.tc_num1.len<-regionCounts$trans.norm.tc_num1/regionCounts$l_tr*1e3
regionCounts$trans.norm.tc_num2.len<-regionCounts$trans.norm.tc_num2/regionCounts$l_tr*1e3

regionCounts$utr3_vs_cds1.norm<-regionCounts$utr3.norm.tc_num1.len/regionCounts$cds.norm.tc_num1.len
regionCounts$utr3_vs_cds2.norm<-regionCounts$utr3.norm.tc_num2.len/regionCounts$cds.norm.tc_num2.len


regs<-subset(regionCounts, select=c("gene_id", "trans.norm.tc_num1","trans.norm.tc_num2","cds.norm.tc_num1","cds.norm.tc_num2","utr3.norm.tc_num1","utr3.norm.tc_num2",
                                    "utr3_vs_cds1", "utr3_vs_cds2","utr3_vs_cds1.norm", "utr3_vs_cds2.norm" ))

regs$trans.norm.tc.mean<-rowMeans(regs[,c("trans.norm.tc_num1", "trans.norm.tc_num2")])
regs$cds.norm.tc.mean<-rowMeans(regs[,c("cds.norm.tc_num1", "cds.norm.tc_num2")])
regs$utr3.norm.tc.mean<-rowMeans(regs[,c("utr3.norm.tc_num1", "utr3.norm.tc_num2")])
regs$utr3_vs_cds.mean<-rowMeans(regs[,c("utr3_vs_cds1", "utr3_vs_cds2")])
regs$utr3_vs_cds.norm.mean<-rowMeans(regs[,c("utr3_vs_cds1.norm", "utr3_vs_cds2.norm")])

ggplot(regs, aes(log2(utr3_vs_cds.mean),log2(trans.norm.tc.mean) ))+geom_point()
ggplot(regs, aes(log2(utr3_vs_cds.norm.mean),log2(trans.norm.tc.mean) ))+geom_point()



###hdlbp targets according to tsignals
library(riboWaltz)
library(GenomicRanges)
library(ggplot2)
library(reshape2)
# data(mm81cdna)
# data(reads_list)
# setwd("/Volumes/landthaler/pcp/projects/miha/koshi_riboseq/ribo/mapping_rerun/bed/")
# setwd("E:/Google Drive/koshi_revision/riboprof/")
# setwd("E:/koshi_codon/agami/")
# check annotation gtf and transcript fa
# anH<-create_annotation("gencode.v19.annotation.gtf", dataSource="gencode.v19", organism="Homo sapiens")
# sequences_biost <- Biostrings::readDNAStringSet("hg19bt1.transcripts.fa",
#                                                 format = "fasta", use.names = TRUE)
# length(names(sequences_biost) %in% anH$transcript)
# nrow(an[names(sequences_biost) %in% anH$transcript,])
# length(sequences_biost[as.character(anH$transcript)])

# save.image("ribowaltz_human.RData")

# setwd("~/Google Drive/koshi_revision/riboprof/")
# # setwd("E:/Google Drive/koshi_revision/riboprof/")
# load("ribowaltz_human.RData")


mas<-read.delim("data/hdlbp_master_table_with_classes_uniq.txt", header=T)
inf<-subset(mas, select=c("gene_id","Symbol",
                          "tpm_cutoff","tc_transcript_norm","tc_CDS_norm","tc_UTR3_norm","localization_cat","log2FoldChange.tot.KO.WT","log2FoldChange.mem.KO.WT","log2FoldChange.cyt.KO.WT","log2FoldChange.tot.cyt.KO.293","log2FoldChange.nuc.cyt.KO.293","log2FoldChange.nuc.tot.2A15.293","log2FoldChange.nuc.tot.3C2.293","log2FoldChange.mem.tot.2A15.293","log2FoldChange.mem.tot.3C2.293", 
                          "mean_te_293","loc_tar_transcript","tc_transcript_norm_cat","loc_tar_CDS","tc_CDS_norm_cat","loc_tar_UTR3",
                          "tc_transcript_norm", "tc_UTR3_norm",
                          "log2FoldChange.mem.cyt.293_1",
                          "log2FoldChange.mem.cyt.293_2",
                          "log2FoldChange.mem.cyt.293",
                          "log2FoldChange.mem.cyt.KO.293",
                          "log2FoldChange.ribo.rna.KO.WT",
                          "Transmembrane.helices" ,"Cleavage.site..Signalp." ))



trans<-read.delim("data/considered_gene_names.txt", header=T)[,c(1:10,16,18)]
colnames(trans)[10]<-"gene_id"

inf<-merge(inf, trans, by="gene_id")



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

##make bed file for TM,SP
# bedi<-merge(anH, tsig, by.x="transcript", by.y="Transcript.stable.ID.version")
# sps<-subset(bedi, !is.na(bedi$Cleavage.site..Signalp..start))
# tms<-subset(bedi, !is.na(bedi$Transmembrane.helices.start))
# 
# sps$start<-sps$l_utr5+3*(sps$Cleavage.site..Signalp..start-1)
# sps$stop<-sps$l_utr5+3*(sps$Cleavage.site..Signalp..end-1)
# sps$num<-1
# 
# tms$start<-tms$l_utr5+3*(tms$Transmembrane.helices.start-1)
# tms$stop<-tms$l_utr5+3*(tms$Transmembrane.helices.end-1)
# tms$num<-1
# 
# bedi<-rbind(sps,tms)
#write.table(bedi[,c("transcript", "start","stop","tsig")], "tsig.bed", quote=F, sep="\t", row.names=F, col.names = F)

#continue
tsig$tmcount<-1
dms<-aggregate(start~ Transcript.stable.ID.version+Gene.stable.ID.version+tsig, data=tsig, min)
tmc<-subset(tsig, tsig!="SignalP-noTM-only")
tmc<-aggregate(tmcount~ Transcript.stable.ID.version+Gene.stable.ID.version, data=tmc, sum)

# dms<-subset(dms, !duplicated(Gene.stable.ID.version))

ggplot(dms, aes(start))+geom_histogram(bins=1000)+facet_wrap(~tsig)+coord_cartesian(xlim=c(-10,100))
ggplot(tmc, aes(tmcount))+geom_histogram()

glock<-which(dms$Transcript.stable.ID.version %in% inf$transcript)
glock2<-which(tmc$Transcript.stable.ID.version %in% inf$transcript)

tarsig<-dms[glock,]
tarsig2<-tmc[glock2,]
tarsig<-merge(tarsig, tarsig2, by="Gene.stable.ID.version", all.x=T)
tarsig<-tarsig[,-5]
colnames(tarsig)[2]<-"Transcript.stable.ID.version"

ggplot(tarsig, aes(tsig))+geom_bar()

tsig_annot<-merge(inf, tarsig, by.x="gene_id", by.y="Gene.stable.ID.version", all.x=T)
tsig_annot$tsig<-ifelse(is.na(tsig_annot$tsig) & tsig_annot$localization_cat=="membrane", "mem_notsig", 
                        ifelse(is.na(tsig_annot$tsig) & tsig_annot$localization_cat=="cytosolic", "cyt_notsig",
                               ifelse(is.na(tsig_annot$localization_cat) & tsig_annot$tpm_cutoff>=10 , "cytANDmem_notsig", tsig_annot$tsig)))
ggplot(tsig_annot, aes(tsig))+geom_bar()
nrow(subset(tsig_annot, as.character(transcript)==as.character(Transcript.stable.ID.version)))

man<-subset(tsig_annot, !is.na(tsig))

man$nucstart<-(man$start-1)*3
man$minusstart<-ifelse(man$nucstart==0,0, -man$nucstart)
man$plusstart<-ifelse(man$l_cds>=1502, 1502+man$minusstart,man$l_cds+man$minusstart )

man$dist_stop<-man$l_cds-man$nucstart

man$tsig<-ifelse(man$dist_stop<=150 & man$tsig=="TMhelix-only" , "TailAnchored", man$tsig)
man$tsig<-ifelse(grepl("MT-", man$Symbol), "MitoEncoded", man$tsig)

man$tsig<-ifelse(man$tsig=="SignalP-noTM-TM" |  man$tsig=="SignalP-TM-TMhelix" |  man$tsig=="SignalP-TM-only" , "SignalP-TM", man$tsig)
mito<-read.delim("data/mitocarta2.txt", header=T)
man$ens_gene_id<-gsub("\\..*","",man$gene_id)
stud<-merge(man, mito, by.x="ens_gene_id", by.y="gene_id", all.x=T)
stud$tsig.x<-ifelse(is.na(stud$mito), stud$tsig.x,
                  ifelse(stud$tsig.x=="cyt_notsig" & stud$mito=="MitoCarta", "MitoCarta", stud$tsig))

man<-stud

man<-merge(man, regs, by="gene_id", all.x=T)

ggplot(man[man$tpm_cutoff>=10 &  man$gene_biotype=="protein_coding" ,],
       aes(loc_tar_transcript, log2(utr3_vs_cds1.norm),fill=loc_tar_transcript))+geom_violin(scale="area",na.rm=T)+
  geom_boxplot(fill="white", na.rm=T, width=0.15, outlier.shape = NA)+theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" ,],
       aes(loc_tar_transcript, log2(utr3_vs_cds2.norm),fill=loc_tar_transcript))+geom_violin(scale="area",na.rm=T)+
  geom_boxplot(fill="white", na.rm=T, width=0.15, outlier.shape = NA)+theme(axis.text.x = element_text(angle = 45, hjust = 1))




ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat),],
       aes(tc_transcript_norm_cat, log2(utr3_vs_cds1.norm),fill=localization_cat))+geom_violin(scale="area",na.rm=T,position=position_dodge())+
  geom_boxplot(width=0.1,na.rm=T, position=position_dodge(width=0.9), outlier.shape = NA)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+scale_fill_manual(values=c("dodgerblue2","orange3"))

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat),],
       aes(tc_transcript_norm_cat, log2(tpm_cutoff),fill=localization_cat))+geom_violin(scale="area",na.rm=T,position=position_dodge())+
  geom_boxplot(width=0.1,na.rm=T, position=position_dodge(width=0.9), outlier.shape = NA)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+scale_fill_manual(values=c("dodgerblue2","orange3"))


ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" ,],
       aes(tsig, log2(utr3_vs_cds1.norm),fill=tsig))+geom_violin(scale="area",na.rm=T,position=position_dodge())+
  geom_boxplot(width=0.1,na.rm=T, position=position_dodge(width=0.9), outlier.shape = NA)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values=c("dodgerblue2","dodgerblue2","dodgerblue2","green4","orange3","red3","red3","red3"))

ggplot(man, aes(log2(utr3_vs_cds.mean),log2(trans.norm.tc.mean) ))+geom_point()
ggplot(man, aes(log2(utr3_vs_cds.norm.mean),log2(trans.norm.tc.mean) ))+geom_point()

man$trans.norm.tc.mean.exp<- ifelse(man$tpm_cutoff==0, NA, man$trans.norm.tc.mean/man$tpm_cutoff)
man$cds.norm.tc.mean.exp<- ifelse(man$tpm_cutoff==0, NA, man$cds.norm.tc.mean/man$tpm_cutoff)
man$utr3.norm.tc.mean.exp<- ifelse(man$tpm_cutoff==0, NA, man$utr3.norm.tc.mean/man$tpm_cutoff)

ggplot(man, aes(log2(utr3_vs_cds.norm.mean),log2(trans.norm.tc.mean/tpm_cutoff) ))+geom_point()
#fig2b
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat) ,], aes(-log2(utr3_vs_cds.norm.mean),log2(trans.norm.tc.mean.exp), colour=localization_cat ))+geom_point(shape=1, alpha=0.7)+scale_colour_manual(values=c("dodgerblue2","orange3"))
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat),], aes(log2(utr3_vs_cds.norm.mean), colour=localization_cat ))+geom_density()

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat) ,], aes(-log2(utr3_vs_cds.norm.mean),log2(cds.norm.tc.mean.exp), colour=localization_cat ))+geom_point(shape=1, alpha=0.7)+scale_colour_manual(values=c("dodgerblue2","orange3"))
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat) ,], aes(-log2(utr3_vs_cds.norm.mean),log2(utr3.norm.tc.mean.exp), colour=localization_cat ))+geom_point(shape=1, alpha=0.7)+scale_colour_manual(values=c("dodgerblue2","orange3"))


thr<-10
quants<-round(quantile(man[man$tpm_cutoff>=thr,"trans.norm.tc.mean.exp"], probs=c(0,1/3,2/3,1), na.rm=T),digits=2)
man$trans.norm.tc.mean.exp_cat<-ifelse(is.na(man$trans.norm.tc.mean.exp) & man$tpm_cutoff<thr, NA,
                                ifelse(is.na(man$trans.norm.tc.mean.exp) & man$tpm_cutoff>=thr, "nontarget",
                                       ifelse(man$tpm_cutoff>=thr & man$trans.norm.tc.mean.exp<=quants[2] & man$trans.norm.tc.mean.exp>quants[1], paste0("tc<",quants[2]),
                                              ifelse(man$tpm_cutoff>=thr & man$trans.norm.tc.mean.exp<=quants[3] & man$trans.norm.tc.mean.exp>quants[2], paste0("tc>",quants[2]," & tc<",quants[3]),
                                                     ifelse(man$tpm_cutoff>=thr & man$trans.norm.tc.mean.exp<=quants[4] & man$trans.norm.tc.mean.exp>quants[3], paste0("tc>",quants[3]," & tc<",quants[4]), NA)))))

quants_mem<-round(quantile(man[man$tpm_cutoff>=thr &man$localization_cat=="membrane","utr3_vs_cds.norm.mean"], probs=c(0,1/3,2/3,1), na.rm=T),digits=2)
quants_cyt<-round(quantile(man[man$tpm_cutoff>=thr &man$localization_cat=="cytosolic","utr3_vs_cds.norm.mean"], probs=c(0,1/3,2/3,1), na.rm=T),digits=2)
man$loc_utr3_vs_cds<-ifelse(is.na(man$utr3_vs_cds.norm.mean) & man$tpm_cutoff<thr , NA,
                             ifelse(is.na(man$utr3_vs_cds.norm.mean) & man$tpm_cutoff>=thr &man$localization_cat=="membrane", NA,
                                    ifelse(man$tpm_cutoff>=thr & man$utr3_vs_cds.norm.mean<=quants_mem[2] & man$utr3_vs_cds.norm.mean>quants_mem[1] &man$localization_cat=="membrane", paste0("membrane_utr3_cds<",quants_mem[2]),
                                           ifelse(man$tpm_cutoff>=thr & man$utr3_vs_cds.norm.mean<=quants_mem[3] & man$utr3_vs_cds.norm.mean>quants_mem[2] &man$localization_cat=="membrane", paste0("membrane_utr3_cds>",quants_mem[2]," & tc<",quants_mem[3]),
                                                  ifelse(man$tpm_cutoff>=thr & man$utr3_vs_cds.norm.mean<=quants_mem[4] & man$utr3_vs_cds.norm.mean>quants_mem[3] &man$localization_cat=="membrane", paste0("membrane_utr3_cds>",quants_mem[3]," & tc<",quants_mem[4]),
                                                         ifelse(is.na(man$utr3_vs_cds.norm.mean) & man$tpm_cutoff>=thr &man$localization_cat=="cytosolic", NA,
                                                                ifelse(man$tpm_cutoff>=thr & man$utr3_vs_cds.norm.mean<=quants_cyt[2] & man$utr3_vs_cds.norm.mean>quants_cyt[1] &man$localization_cat=="cytosolic", paste0("cytosolic_utr3_cds<",quants_cyt[2]),
                                                                       ifelse(man$tpm_cutoff>=thr & man$utr3_vs_cds.norm.mean<=quants_cyt[3] & man$utr3_vs_cds.norm.mean>quants_cyt[2] &man$localization_cat=="cytosolic", paste0("cytosolic_utr3_cds>",quants_cyt[2]," & tc<",quants_cyt[3]),
                                                                              ifelse(man$tpm_cutoff>=thr & man$utr3_vs_cds.norm.mean<=quants_cyt[4] & man$utr3_vs_cds.norm.mean>quants_cyt[3] &man$localization_cat=="cytosolic", paste0("cytosolic_utr3_cds>",quants_cyt[3]," & tc<",quants_cyt[4]),NA)))))))))


quants<-round(quantile(man[man$tpm_cutoff>=thr & man$localization_cat=="cytosolic","trans.norm.tc.mean.exp"], probs=c(0,1/3,2/3,1), na.rm=T),digits=2)
man$cyt.trans.norm.tc.mean.exp_cat<-ifelse(is.na(man$trans.norm.tc.mean.exp) & man$tpm_cutoff<thr, NA,
                                       ifelse(is.na(man$trans.norm.tc.mean.exp) & man$tpm_cutoff>=thr & man$localization_cat=="cytosolic", "nontarget",
                                              ifelse(man$tpm_cutoff>=thr & man$trans.norm.tc.mean.exp<=quants[2] & man$trans.norm.tc.mean.exp>quants[1] & man$localization_cat=="cytosolic", paste0("tc<",quants[2]),
                                                     ifelse(man$tpm_cutoff>=thr & man$trans.norm.tc.mean.exp<=quants[3] & man$trans.norm.tc.mean.exp>quants[2]& man$localization_cat=="cytosolic", paste0("tc>",quants[2]," & tc<",quants[3]),
                                                            ifelse(man$tpm_cutoff>=thr & man$trans.norm.tc.mean.exp<=quants[4] & man$trans.norm.tc.mean.exp>quants[3]& man$localization_cat=="cytosolic", paste0("tc>",quants[3]," & tc<",quants[4]), NA)))))

quants<-round(quantile(man[man$tpm_cutoff>=thr & man$localization_cat=="cytosolic","cds.norm.tc.mean.exp"], probs=c(0,1/3,2/3,1), na.rm=T),digits=2)
man$cyt.cds.norm.tc.mean.exp_cat<-ifelse(is.na(man$cds.norm.tc.mean.exp) & man$tpm_cutoff<thr, NA,
                                           ifelse(is.na(man$cds.norm.tc.mean.exp) & man$tpm_cutoff>=thr & man$localization_cat=="cytosolic", "nontarget",
                                                  ifelse(man$tpm_cutoff>=thr & man$cds.norm.tc.mean.exp<=quants[2] & man$cds.norm.tc.mean.exp>quants[1] & man$localization_cat=="cytosolic", paste0("tc<",quants[2]),
                                                         ifelse(man$tpm_cutoff>=thr & man$cds.norm.tc.mean.exp<=quants[3] & man$cds.norm.tc.mean.exp>quants[2]& man$localization_cat=="cytosolic", paste0("tc>",quants[2]," & tc<",quants[3]),
                                                                ifelse(man$tpm_cutoff>=thr & man$cds.norm.tc.mean.exp<=quants[4] & man$cds.norm.tc.mean.exp>quants[3]& man$localization_cat=="cytosolic", paste0("tc>",quants[3]," & tc<",quants[4]), NA)))))
quants<-round(quantile(man[man$tpm_cutoff>=thr & man$localization_cat=="cytosolic","utr3.norm.tc.mean.exp"], probs=c(0,1/3,2/3,1), na.rm=T),digits=2)
man$cyt.utr3.norm.tc.mean.exp_cat<-ifelse(is.na(man$utr3.norm.tc.mean.exp) & man$tpm_cutoff<thr, NA,
                                         ifelse(is.na(man$utr3.norm.tc.mean.exp) & man$tpm_cutoff>=thr & man$localization_cat=="cytosolic", "nontarget",
                                                ifelse(man$tpm_cutoff>=thr & man$utr3.norm.tc.mean.exp<=quants[2] & man$utr3.norm.tc.mean.exp>quants[1] & man$localization_cat=="cytosolic", paste0("tc<",quants[2]),
                                                       ifelse(man$tpm_cutoff>=thr & man$utr3.norm.tc.mean.exp<=quants[3] & man$utr3.norm.tc.mean.exp>quants[2]& man$localization_cat=="cytosolic", paste0("tc>",quants[2]," & tc<",quants[3]),
                                                              ifelse(man$tpm_cutoff>=thr & man$utr3.norm.tc.mean.exp<=quants[4] & man$utr3.norm.tc.mean.exp>quants[3]& man$localization_cat=="cytosolic", paste0("tc>",quants[3]," & tc<",quants[4]), NA)))))

quants<-round(quantile(man[man$tpm_cutoff>=thr & man$localization_cat=="membrane","cds.norm.tc.mean.exp"], probs=c(0,1/3,2/3,1), na.rm=T),digits=2)
man$mem.cds.norm.tc.mean.exp_cat<-ifelse(is.na(man$cds.norm.tc.mean.exp) & man$tpm_cutoff<thr, NA,
                                         ifelse(is.na(man$cds.norm.tc.mean.exp) & man$tpm_cutoff>=thr & man$localization_cat=="membrane", "nontarget",
                                                ifelse(man$tpm_cutoff>=thr & man$cds.norm.tc.mean.exp<=quants[2] & man$cds.norm.tc.mean.exp>quants[1] & man$localization_cat=="membrane", paste0("tc<",quants[2]),
                                                       ifelse(man$tpm_cutoff>=thr & man$cds.norm.tc.mean.exp<=quants[3] & man$cds.norm.tc.mean.exp>quants[2]& man$localization_cat=="membrane", paste0("tc>",quants[2]," & tc<",quants[3]),
                                                              ifelse(man$tpm_cutoff>=thr & man$cds.norm.tc.mean.exp<=quants[4] & man$cds.norm.tc.mean.exp>quants[3]& man$localization_cat=="membrane", paste0("tc>",quants[3]," & tc<",quants[4]), NA)))))



ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)  & !is.na(man$cyt.trans.norm.tc.mean.exp_cat)& man$localization_cat!="membrane",],
       aes(-log2FoldChange.tot.cyt.KO.293, colour=cyt.trans.norm.tc.mean.exp_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("black","red3","dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)  & !is.na(man$cyt.trans.norm.tc.mean.exp_cat)& man$localization_cat!="membrane",],
       aes(log2FoldChange.nuc.tot.2A15.293, colour=cyt.trans.norm.tc.mean.exp_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("black","red3","dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)  & !is.na(man$cyt.trans.norm.tc.mean.exp_cat)& man$localization_cat!="membrane",],
       aes(log2FoldChange.nuc.tot.3C2.293, colour=cyt.trans.norm.tc.mean.exp_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("black","red3","dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))


ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane",],
       aes(log2FoldChange.cyt.KO.WT, colour=cyt.trans.norm.tc.mean.exp_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("black","red3","dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))


ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane",],
       aes(-log2FoldChange.tot.cyt.KO.293, colour=cyt.cds.norm.tc.mean.exp_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("black","red3","dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane",],
       aes(-log2FoldChange.tot.cyt.KO.293, colour=cyt.utr3.norm.tc.mean.exp_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("black","red3","dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane",],
       aes(log2FoldChange.mem.cyt.KO.293, colour=cyt.trans.norm.tc.mean.exp_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("black","red3","dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane",],
       aes(log2FoldChange.mem.cyt.KO.293, colour=cyt.cds.norm.tc.mean.exp_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("black","red3","dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane",],
       aes(log2FoldChange.mem.cyt.KO.293, colour=cyt.utr3.norm.tc.mean.exp_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("black","red3","dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat=="membrane",],
       aes(log2FoldChange.mem.cyt.KO.293, colour=mem.cds.norm.tc.mean.exp_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("black","red3","dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)  ,],
       aes(-log2FoldChange.tot.cyt.KO.293, colour=localization_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)  ,],
       aes(-log2FoldChange.tot.cyt.KO.293, colour=localization_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)  ,],
       aes(log2FoldChange.mem.tot.3C2.293, colour=localization_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("black","red3","dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-1,1))

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$trans.norm.tc.mean.exp)  ,],
       aes(log2FoldChange.mem.tot.2A15.293, -log2FoldChange.tot.cyt.KO.293, colour=log2(cds.norm.tc.mean.exp)))+geom_point(shape=1, alpha=0.6)+geom_abline()+coord_cartesian(xlim=c(-1,1),ylim=c(-1,1))+facet_wrap(~localization_cat)
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$utr3_vs_cds.norm.mean)  ,],
       aes(log2FoldChange.mem.tot.2A15.293, -log2FoldChange.tot.cyt.KO.293, colour=-log2(utr3_vs_cds.norm.mean)))+geom_point(shape=1, alpha=0.6)+geom_abline()+coord_cartesian(xlim=c(-1,1),ylim=c(-1,1))+facet_wrap(~localization_cat)
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$utr3_vs_cds.norm.mean)  ,],
       aes(log2FoldChange.mem.cyt.KO.293, -log2FoldChange.tot.cyt.KO.293, colour=-log2(utr3_vs_cds.norm.mean)))+geom_point(shape=1, alpha=0.6)+geom_abline()+coord_cartesian(xlim=c(-2,2),ylim=c(-2,2))+facet_wrap(~localization_cat)
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$utr3_vs_cds.norm.mean)  ,],
       aes(log2FoldChange.mem.cyt.KO.293, log2FoldChange.tot.KO.WT, colour=-log2(utr3_vs_cds.norm.mean)))+geom_point(shape=1, alpha=0.6)+geom_abline()+coord_cartesian(xlim=c(-1,1),ylim=c(-1,1))+facet_wrap(~localization_cat)
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$utr3_vs_cds.norm.mean)  ,],
       aes(log2FoldChange.mem.cyt.KO.293, log2FoldChange.cyt.KO.WT, colour=-log2(utr3_vs_cds.norm.mean)))+geom_point(shape=1, alpha=0.6)+geom_abline()+coord_cartesian(xlim=c(-1,1),ylim=c(-1,1))+facet_wrap(~localization_cat)
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$utr3_vs_cds.norm.mean)  ,],
       aes(log2FoldChange.mem.cyt.KO.293, log2FoldChange.ribo.rna.KO.WT, colour=-log2(utr3_vs_cds.norm.mean)))+geom_point(shape=1, alpha=0.6)+geom_abline()+coord_cartesian(xlim=c(-1,1),ylim=c(-1,1))+facet_wrap(~localization_cat)
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$trans.norm.tc.mean.exp)  ,],
       aes(log2FoldChange.mem.cyt.KO.293, log2FoldChange.ribo.rna.KO.WT, colour=log2(trans.norm.tc.mean.exp)))+geom_point(shape=1, alpha=0.6)+geom_abline()+coord_cartesian(xlim=c(-1,1),ylim=c(-1,1))+facet_wrap(~localization_cat)
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$cds.norm.tc.mean.exp)  ,],
       aes(log2FoldChange.mem.cyt.KO.293, log2FoldChange.ribo.rna.KO.WT, colour=log2(cds.norm.tc.mean.exp)))+geom_point(shape=1, alpha=0.6)+geom_abline()+coord_cartesian(xlim=c(-1,1),ylim=c(-1,1))+facet_wrap(~localization_cat)

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)  ,],
       aes(log2FoldChange.cyt.KO.WT, colour=localization_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("black","red3","dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat) & man$localization_cat=="membrane"  ,],
       aes(-log2FoldChange.tot.cyt.KO.293, colour=loc_tar_CDS))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("black","red3","dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat) & man$localization_cat=="cytosolic"  ,],
       aes(-log2FoldChange.tot.cyt.KO.293, colour=loc_tar_CDS))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("black","red3","dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat) & man$localization_cat=="membrane"  ,],
       aes(log2FoldChange.mem.KO.WT, colour=loc_tar_CDS))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("black","red3","dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat) & man$localization_cat=="cytosolic"  ,],
       aes(log2FoldChange.mem.KO.WT, colour=loc_tar_CDS))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("black","red3","dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))


ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane",],
       aes(log2FoldChange.ribo.rna.KO.WT, colour=cyt.trans.norm.tc.mean.exp_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("black","red3","dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane",],
       aes(log2FoldChange.ribo.rna.KO.WT, colour=cyt.utr3.norm.tc.mean.exp_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("black","red3","dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane",],
       aes(log2FoldChange.ribo.rna.KO.WT, colour=cyt.cds.norm.tc.mean.exp_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("black","red3","dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat=="membrane",],
       aes(log2FoldChange.ribo.rna.KO.WT, colour=mem.cds.norm.tc.mean.exp_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("black","red3","dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))


count(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane" & !is.na(man$log2FoldChange.tot.cyt.KO.293),], vars="cyt.trans.norm.tc.mean.exp_cat")
count(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane" & !is.na(man$log2FoldChange.tot.cyt.KO.293),], vars="cyt.cds.norm.tc.mean.exp_cat")
count(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane" & !is.na(man$log2FoldChange.tot.cyt.KO.293),], vars="cyt.utr3.norm.tc.mean.exp_cat")
count(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat=="membrane" & !is.na(man$log2FoldChange.tot.cyt.KO.293),], vars="mem.cds.norm.tc.mean.exp_cat")
count(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat=="membrane" & !is.na(man$log2FoldChange.tot.cyt.KO.293),], vars="loc_tar_CDS")

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane",],
       aes(log2FoldChange.ribo.rna.KO.WT, colour=cyt.trans.norm.tc.mean.exp_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("black","red3","dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane",],
       aes(log2FoldChange.ribo.rna.KO.WT, colour=cyt.cds.norm.tc.mean.exp_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("black","red3","dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane",],
       aes(log2FoldChange.ribo.rna.KO.WT, colour=trans.norm.tc.mean.exp_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("black","red3","dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane",],
       aes(log2FoldChange.ribo.rna.KO.WT, colour=loc_tar_UTR3))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("black","red3","dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane",],
       aes(-log2FoldChange.tot.cyt.KO.293, colour=trans.norm.tc.mean.exp_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("black","red3","dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))
count(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane" & !is.na(man$log2FoldChange.tot.cyt.KO.293),], vars="trans.norm.tc.mean.exp_cat")

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane",], aes(trans.norm.tc.mean.exp_cat))+geom_bar()

ggplot(man[man$tpm_cutoff>=1 & man$localization_cat!="cytosolic" & !is.na(man$loc_tar_CDS) & man$gene_biotype=="protein_coding" & man$tsig!="MitoEncoded",],
       aes(log2FoldChange.ribo.rna.KO.WT,colour=loc_tar_CDS))+stat_ecdf()+
  scale_colour_manual(values=c("red3","dodgerblue2","orange3","black"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))
count(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat=="membrane" & !is.na(man$log2FoldChange.ribo.rna.KO.WT),], vars="loc_tar_CDS")


ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane",],
       aes(log2FoldChange.cyt.KO.WT, colour=trans.norm.tc.mean.exp_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("black","red3","dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))
count(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane" & !is.na(man$log2FoldChange.cyt.KO.WT),], vars="trans.norm.tc.mean.exp_cat")
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane",],
       aes(log2FoldChange.tot.KO.WT, colour=trans.norm.tc.mean.exp_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("black","red3","dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane",],
       aes(log2FoldChange.nuc.cyt.KO.293, colour=trans.norm.tc.mean.exp_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("black","red3","dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))


ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane",],
       aes(-log2FoldChange.tot.cyt.KO.293, colour=loc_utr3_vs_cds))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+
  scale_colour_manual(values=c("black","red3","dodgerblue2","orange3"))+coord_cartesian(xlim=c(-1,1))+coord_cartesian(xlim=c(-.5,1))


ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat) &!is.na(man$loc_tar_CDS) & man$localization_cat!="cytosolic",],
       aes(log2FoldChange.ribo.rna.KO.WT, colour=loc_tar_CDS))+stat_ecdf()+coord_cartesian(xlim=c(-.5,.5))
  
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)  &!is.na(man$loc_tar_UTR3) & man$localization_cat!="membrane",],
       aes(log2FoldChange.cyt.KO.WT, colour=loc_tar_UTR3))+stat_ecdf()+coord_cartesian(xlim=c(-.5,.5))

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)  &!is.na(man$loc_tar_UTR3) & man$localization_cat!="membrane",],
       aes(log2FoldChange.cyt.KO.WT, colour=trans.norm.tc.mean.exp_cat))+stat_ecdf()+coord_cartesian(xlim=c(-.5,.5))
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)  &!is.na(man$loc_tar_UTR3) & man$localization_cat!="membrane",],
       aes(log2FoldChange.cyt.KO.WT, colour=loc_tar_CDS))+stat_ecdf()+coord_cartesian(xlim=c(-.5,.5))
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)  &!is.na(man$loc_tar_UTR3) & man$localization_cat!="membrane",],
       aes(log2FoldChange.cyt.KO.WT, colour=loc_tar_UTR3))+stat_ecdf()+coord_cartesian(xlim=c(-.5,.5))
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)  &!is.na(man$loc_tar_UTR3) & man$localization_cat!="membrane",],
       aes(log2FoldChange.cyt.KO.WT, colour=loc_utr3_vs_cds))+stat_ecdf()+coord_cartesian(xlim=c(-.5,.5))

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)  &!is.na(man$loc_tar_UTR3) & man$localization_cat!="cytosolic",],
       aes(log2FoldChange.ribo.rna.KO.WT, colour=loc_utr3_vs_cds))+stat_ecdf()+coord_cartesian(xlim=c(-.5,.5))
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)  &!is.na(man$loc_tar_UTR3) & man$localization_cat!="membrane",],
       aes(log2FoldChange.ribo.rna.KO.WT, colour=loc_utr3_vs_cds))+stat_ecdf()+coord_cartesian(xlim=c(-.5,.5))
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane",],
       aes(log2FoldChange.cyt.KO.WT, colour=trans.norm.tc.mean.exp_cat))+stat_ecdf()+coord_cartesian(xlim=c(-.5,.5))
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane",],
       aes(-log2FoldChange.tot.cyt.KO.293, colour=trans.norm.tc.mean.exp_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane",],
       aes(-log2FoldChange.tot.cyt.KO.293, colour=loc_tar_UTR3))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat=="membrane",],
       aes(log2FoldChange.ribo.rna.KO.WT, colour=trans.norm.tc.mean.exp_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))


ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)  & man$localization_cat!="membrane",],
       aes(-log2FoldChange.tot.cyt.KO.293, colour=trans.norm.tc.mean.exp_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)  ,],
       aes(-log2FoldChange.tot.cyt.KO.293, colour=loc_utr3_vs_cds))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)  ,],
       aes(log2FoldChange.ribo.rna.KO.WT, colour=loc_utr3_vs_cds))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))


ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat) & man$trans.norm.tc.mean.exp_cat!="nontarget",],
         aes(trans.norm.tc.mean.exp_cat, -log2(utr3_vs_cds.norm.mean),fill=localization_cat))+geom_violin(scale="area",na.rm=T,position=position_dodge())+
  geom_boxplot(width=0.1,na.rm=T, position=position_dodge(width=0.9), outlier.shape = NA)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values=c("dodgerblue2","orange3"))

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat) & man$trans.norm.tc.mean.exp_cat!="nontarget",],
       aes(cyt.trans.norm.tc.mean.exp_cat, -log2(utr3_vs_cds.norm.mean),fill=localization_cat))+geom_violin(scale="area",na.rm=T,position=position_dodge())+
  geom_boxplot(width=0.1,na.rm=T, position=position_dodge(width=0.9), outlier.shape = NA)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values=c("dodgerblue2","orange3"))

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat) & man$trans.norm.tc.mean.exp_cat!="nontarget",],
       aes(tsig, log2(utr3_vs_cds.norm.mean),fill=tsig))+geom_violin(scale="area",na.rm=T,position=position_dodge())+
  geom_boxplot(width=0.1,na.rm=T, position=position_dodge(width=0.9), outlier.shape = NA)+theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat) ,],
       aes(log2FoldChange.ribo.rna.KO.WT, colour=loc_utr3_vs_cds))+stat_ecdf()+coord_cartesian(xlim=c(-0.5,0.5))
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat) & man$trans.norm.tc.mean.exp_cat!="nontarget",],
       aes(loc_utr3_vs_cds, log2FoldChange.ribo.rna.KO.WT,fill=loc_utr3_vs_cds))+geom_violin(scale="area",na.rm=T,position=position_dodge())+
  geom_boxplot(width=0.1,na.rm=T, position=position_dodge(width=0.9), outlier.shape = NA)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+coord_cartesian(ylim=c(-.5,.5))



ggplot(man, aes(tsig))+geom_bar()
ggplot(man[man$tpm_cutoff>=10,],aes(tsig,log2(tc_CDS_norm),fill=tsig))+geom_boxplot()+geom_point(aes(fill = tsig), size = 1, shape = 1, position = position_jitterdodge())
man$tsig<-factor(man$tsig, levels=c("cyt_notsig", "MitoCarta", "TailAnchored", "cytANDmem_notsig","MitoEncoded", "mem_notsig","SignalP-noTM-only","SignalP-TM", "TMhelix-only"))
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & man$tsig!="cytANDmem_notsig"& man$tsig!="MitoEncoded",],
       aes(tsig,log2(tc_transcript_norm),fill=tsig))+geom_violin(scale="area",na.rm=T)+
  geom_boxplot(fill="white", na.rm=T, width=0.3, outlier.shape = NA)+
  scale_fill_manual(values=c("dodgerblue2","dodgerblue2","dodgerblue2","orange3","red3","red3","red3"))+theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & man$tsig!="cytANDmem_notsig"& man$tsig!="MitoEncoded",],
       aes(tsig,log2(tc_UTR3_norm),fill=tsig))+geom_violin(scale="area",na.rm=T)+
  geom_boxplot(fill="white", na.rm=T, width=0.3, outlier.shape = NA)+
  scale_fill_manual(values=c("dodgerblue2","dodgerblue2","dodgerblue2","orange3","red3","red3","red3"))+theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & man$tsig!="cytANDmem_notsig"& man$tsig!="MitoEncoded",],
       aes(tsig,log2(tc_CDS_norm),fill=tsig))+geom_violin(scale="area",na.rm=T)+
  geom_boxplot(fill="white", na.rm=T, width=0.3, outlier.shape = NA)+
  scale_fill_manual(values=c("dodgerblue2","dodgerblue2","dodgerblue2","orange3","red3","red3","red3"))+theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & man$tsig!="cytANDmem_notsig",],
       aes(tsig,(log2FoldChange.mem.cyt.293),fill=tsig))+geom_violin(scale="area",na.rm=T)+
  geom_boxplot(fill="white", na.rm=T, width=0.15, outlier.shape = NA)+
  scale_fill_manual(values=c("dodgerblue2","dodgerblue2","dodgerblue2","green4","orange3","red3","red3","red3"))+theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & man$tsig!="cytANDmem_notsig",],
       aes(tsig,log2(mean_te_293),fill=tsig))+geom_violin(scale="area",na.rm=T)+
  geom_boxplot(fill="white", na.rm=T, width=0.15, outlier.shape = NA)+
  scale_fill_manual(values=c("dodgerblue2","dodgerblue2","dodgerblue2","green4","orange3","red3","red3","red3"))+theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & man$tsig!="MitoEncoded",],
       aes(tsig,log2(utr3_vs_cds.norm.mean),fill=tsig))+geom_violin(scale="area",na.rm=T)+
  geom_boxplot(fill="white", na.rm=T, width=0.15, outlier.shape = NA)+
  scale_fill_manual(values=c("dodgerblue2","dodgerblue2","dodgerblue2","green4","orange3","red3","red3","red3"))+theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & man$tsig!="MitoEncoded",],
       aes(log2FoldChange.ribo.rna.KO.WT,colour=tsig))+stat_ecdf()+coord_cartesian(xlim=c(-.5,.5))


ggplot(man[man$tpm_cutoff>=1 & man$localization_cat!="cytosolic" & !is.na(man$loc_tar_CDS) & man$gene_biotype=="protein_coding" & man$tsig!="MitoEncoded",],
       aes(log2FoldChange.ribo.rna.KO.WT,colour=loc_tar_CDS))+stat_ecdf()+
  scale_colour_manual(values=c("red3","dodgerblue2","orange3","black"))+coord_cartesian(xlim=c(-1,1))

ggplot(man[man$tpm_cutoff>=1 & man$localization_cat!="membrane" & !is.na(man$loc_tar_CDS) & man$gene_biotype=="protein_coding" & man$tsig!="MitoEncoded",],
       aes(log2FoldChange.tot.cyt.KO.293,colour=loc_tar_UTR3))+stat_ecdf()+
  scale_colour_manual(values=c("red3","dodgerblue2","orange3","black"))+coord_cartesian(xlim=c(-1,1))


ggplot(man[man$tpm_cutoff>=1 & man$localization_cat!="cytosolic" & !is.na(man$loc_tar_CDS) & man$gene_biotype=="protein_coding" & man$tsig!="MitoEncoded",],
       aes(log2FoldChange.mem.cyt.KO.293,colour=loc_tar_CDS))+stat_ecdf()+
  scale_colour_manual(values=c("red3","dodgerblue2","orange3","black"))


ggplot(man[man$tpm_cutoff>=1 & man$localization_cat!="cytosolic" & !is.na(man$loc_tar_CDS) & man$gene_biotype=="protein_coding" & man$tsig!="MitoEncoded",],
       aes(log2FoldChange.mem.cyt.KO.293,colour=loc_tar_transcript))+stat_ecdf()+
  scale_colour_manual(values=c("red3","dodgerblue2","orange3","black"))


ggplot(man[man$tpm_cutoff>=1 & man$localization_cat!="cytosolic" & !is.na(man$loc_tar_CDS) & man$gene_biotype=="protein_coding" & man$tsig!="MitoEncoded",],
       aes(log2FoldChange.ribo.rna.KO.WT,colour=loc_tar_transcript))+stat_ecdf()+
  scale_colour_manual(values=c("red3","dodgerblue2","orange3","black"))

  
  
ggplot(man[!is.na(man$utr3_vs_cds.norm.mean),], aes(log2FoldChange.mem.cyt.293_1, log2FoldChange.mem.cyt.293_2, colour=-log2(utr3_vs_cds.norm.mean) ))+geom_point(shape=1, alpha=0.6)
ggplot(man[!is.na(man$utr3_vs_cds.norm.mean),], aes(log2FoldChange.mem.cyt.293, -log2FoldChange.tot.cyt.KO.293, colour=-log2(utr3_vs_cds.norm.mean) ))+geom_point(shape=1, alpha=0.6)+coord_cartesian(ylim=c(-1,1))

ggplot(man[!is.na(man$utr3_vs_cds.norm.mean),], aes(log2(cds.norm.tc.mean.exp),log2(utr3.norm.tc.mean.exp), colour=-log2(utr3_vs_cds.norm.mean) ))+geom_point(shape=1, alpha=0.6)
ggplot(man[!is.na(man$utr3_vs_cds.norm.mean),], aes(log2(cds.norm.tc.mean.exp),log2(utr3.norm.tc.mean.exp), colour=log2FoldChange.mem.cyt.293 ))+geom_point(shape=1, alpha=0.6)+geom_abline(slope=1)
ggplot(man[!is.na(man$utr3_vs_cds.norm.mean),], aes(log2(cds.norm.tc.mean.exp),log2(utr3.norm.tc.mean.exp), colour=-log2FoldChange.tot.cyt.KO.293 ))+geom_point(shape=1, alpha=0.6)+geom_abline(slope=1)+scale_color_gradient(limits = c(-.5,.5))
ggplot(man[!is.na(man$utr3_vs_cds.norm.mean),], aes(log2(cds.norm.tc.mean.exp),log2(utr3.norm.tc.mean.exp), colour=localization_cat ))+geom_point(shape=1, alpha=0.6)+geom_abline(slope=1)
ggplot(man[!is.na(man$utr3_vs_cds.norm.mean),], aes(log2(cds.norm.tc.mean.exp),log2(utr3.norm.tc.mean.exp), colour=log2FoldChange.ribo.rna.KO.WT ))+geom_point(shape=1, alpha=0.6)+geom_abline(slope=1)+scale_color_gradient(limits = c(-.5,.5))


ggplot(man[!is.na(man$utr3_vs_cds.norm.mean),], aes(log2(cds.norm.tc.mean.exp),log2(utr3.norm.tc.mean.exp), colour=cyt.cds.norm.tc.mean.exp_cat ))+geom_point(shape=1, alpha=0.6)+geom_abline(slope=1)
ggplot(man[!is.na(man$utr3_vs_cds.norm.mean),], aes(log2(cds.norm.tc.mean.exp),log2(utr3.norm.tc.mean.exp), colour=cyt.utr3.norm.tc.mean.exp_cat ))+geom_point(shape=1, alpha=0.6)+geom_abline(slope=1)


ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane",],
       aes(-log2(utr3_vs_cds.norm.mean), colour=cyt.cds.norm.tc.mean.exp_cat))+stat_ecdf()

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane",],
       aes(-log2(utr3_vs_cds.norm.mean), colour=cyt.utr3.norm.tc.mean.exp_cat))+stat_ecdf()
ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane",],
       aes(log2(cds.norm.tc.mean.exp), colour=cyt.utr3.norm.tc.mean.exp_cat))+stat_ecdf()

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$localization_cat)   & man$localization_cat!="membrane",],
       aes(log2FoldChange.ribo.rna.KO.WT, colour=cyt.utr3.norm.tc.mean.exp_cat))+stat_ecdf()+facet_wrap(~cyt.cds.norm.tc.mean.exp_cat)+coord_cartesian(xlim=c(-.5,.5))

ggplot(man[man$tpm_cutoff>=10 & man$gene_biotype=="protein_coding" & !is.na(man$log2FoldChange.mem.cyt.293_1) & !is.na(man$log2FoldChange.mem.cyt.293_2) & !is.na(man$cds.norm.tc.mean.exp) ,],
       aes(log2FoldChange.mem.cyt.293_1,log2FoldChange.mem.cyt.293_2, colour=log2(cds.norm.tc.mean.exp)))+geom_point(shape=1, alpha=0.5)+scale_colour_gradient(low="dodgerblue3",high="orange3" ,limits=c(-5,5))+geom_abline(slope=1, lty=2, colour="grey")

ggplot(man[ man$gene_biotype=="protein_coding" & !is.na(man$log2FoldChange.mem.cyt.293_1) & !is.na(man$log2FoldChange.mem.cyt.293_2) & !is.na(man$cds.norm.tc.mean.exp) ,],
       aes(log2FoldChange.mem.cyt.293, log2(tpm_cutoff),colour=log2(cds.norm.tc.mean.exp)))+geom_point(shape=1, alpha=0.5)+scale_colour_gradient(low="dodgerblue3",high="orange3" ,limits=c(-5,5))


# write.table(man, "hdlbp_master_table_with_classes_uniq_tsig.txt", quote=F, sep="\t", row.names=F)
