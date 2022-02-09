library(ggplot2)
library(reshape2)
library(corrplot)
library(Biostrings)
library(dplyr)
library(ggseqlogo)

len<-read.delim("data/transcript_cds_utr_lenght.txt", header=T)
rownames(len)<-len$transcript
len$cds_start<-len$l_utr5+1
len$cds_stop<-len$l_utr5+len$l_cds
len$frame_start<-len$l_utr5
len$frame_stop<-len$cds_start


#selecting highest expressed isoform
rsem1<-read.delim("data/T_293_1.isoforms.results", header=T)
colnames(rsem1)[2:length(colnames(rsem1))]<-paste0(colnames(rsem1)[2:length(colnames(rsem1))],"_","T_293_1")
rsem2<-read.delim("data/T_293_2.isoforms.results", header=T)
colnames(rsem2)[2:length(colnames(rsem2))]<-paste0(colnames(rsem2)[2:length(colnames(rsem2))],"_","T_293_2")
rsem<-merge(rsem1,rsem2, by="transcript_id")
rsem$TPM<-rowMeans(rsem[,c("TPM_T_293_1","TPM_T_293_2")])
rsem<-rsem[order(rsem$gene_id_T_293_1, rsem$TPM, decreasing = T),]
rsem<-subset(rsem, select=c("transcript_id","gene_id_T_293_1", "TPM"))
rsem<-subset(rsem,!duplicated(gene_id_T_293_1))
rsem<-rsem[order(rsem$TPM, decreasing = T),]
length(unique(rsem$gene_id_T_293_1))
colnames(rsem)[2:3]<-c("gene_id","TPM_transcript")

ham<-merge(len, rsem, by.x="transcript", by.y="transcript_id")
nrow(subset(ham,  l_utr5 > 0 & l_cds > 0 & l_cds%%3 == 0 & l_utr3 > 0))

ham<-subset(ham,  l_utr5 > 0 & l_cds > 0 & l_cds%%3 == 0 & l_utr3 > 0)

relevantBed<-ham
relevantBed$start<-1
relevantBed$num<-1
relevantBed$strand<-"+"
relevantBed<-subset(relevantBed, TPM_transcript>=1)
relevantBed<-relevantBed[order(relevantBed$transcript, relevantBed$start),]

#write.table(relevantBed[,c(1, 12, 2, 10, 13,14)], "considered.transcripts.bed", quote=F, sep="\t", row.names = F, col.names = F)


tc<-read.delim("data/reproducible.hdlbp.TCseq.bed", header=F)
colnames(tc)<-c("transcript_id","tc_start","tc_stop","tc_num1","all_reads1","tc_num2","all_reads2","seq")
library(DESeq2)
fac<-estimateSizeFactorsForMatrix(tc[,c("all_reads1","all_reads2")])
tc$norm_tc_num1<-tc$tc_num1/fac[1]
tc$norm_tc_num2<-tc$tc_num2/fac[2]
tc<-subset(tc, tc_num1/all_reads1<0.95 & tc_num2/all_reads2<0.95) # here can be even more stringent

dat<-merge(tc, ham, by.x="transcript_id", by.y="transcript")

fastapath<-"~/hdlbp_git/hg19bt1.transcripts.fa" # here provide local transcript fasta file - e.g. from RSEM index, cannot include in github repo due to file size
seqTrans <- Biostrings::readDNAStringSet(fastapath, format = "fasta", use.names = TRUE)

sub<- seqTrans[ham$transcript]

seqCds <- Biostrings::subseq(sub, start = ham$cds_start, end = ham$cds_stop)
seqUtr3<- Biostrings::subseq(sub, start = ham$cds_stop+1, end = ham$l_tr)
seqUtr5<- Biostrings::subseq(sub, start = 1, end = ham$l_utr5)

dat$tc_region<-ifelse(dat$tc_stop>=dat$cds_start & dat$tc_stop<=dat$cds_stop, "cds",
                      ifelse(dat$tc_stop<dat$cds_start, "utr5", 
                             ifelse(dat$tc_stop>dat$cds_stop, "utr3", NA)))
dat$tc_region<-factor(dat$tc_region, levels=c("utr5", "cds", "utr3"))

mas<-read.delim("data/hdlbp_master_table_with_classes_uniq.txt", header=T)
inf<-subset(mas, select=c("gene_id","Symbol",
                          "tpm_cutoff","tc_CDS_norm","localization_cat",
                          "mean_te_293","loc_tar_CDS","tc_CDS_norm_cat",
                          "tc_transcript_norm",
                          "log2FoldChange.mem.cyt.293",
                          "log2FoldChange.mem.cyt.KO.293",
                          "log2FoldChange.ribo.rna.KO.WT"))

dat<-merge(dat, inf, by="gene_id")

#figS2a
ggplot(subset(dat, tpm_cutoff>=10 & !is.na(localization_cat)),
       aes(log2(norm_tc_num1), log2(norm_tc_num2)))+geom_point(shape=1, size=0.7)+facet_wrap(~tc_region)+geom_abline(lty=2, colour="grey")
cor(subset(dat, tpm_cutoff>=10 & !is.na(localization_cat) & 
             tc_region=="utr5" , 
           select=c("norm_tc_num1", "norm_tc_num2")))
cor(subset(dat, tpm_cutoff>=10 & !is.na(localization_cat) & 
             tc_region=="cds" , 
           select=c("norm_tc_num1", "norm_tc_num2")))
cor(subset(dat, tpm_cutoff>=10 & !is.na(localization_cat) & 
             tc_region=="utr3" , 
           select=c("norm_tc_num1", "norm_tc_num2")))