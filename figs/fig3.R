# this script requires a transcriptome fasta file (line 74). 
# due to size limitations, we could not provide it in the 
# github repository. we used a transcriptome fasta file generated 
# by RSEM during index build from gencode v19 gtf annotation file.

library(ggplot2)
library(reshape2)
library(corrplot)
library(Biostrings)
library(dplyr)
library(ggseqlogo)
library(ggpubr)

#fig3i
dat<-read.delim("data/kd.txt", header=T)

mel<-melt(dat[,1:4])

mel$RNA_oligo<-factor(mel$RNA_oligo, levels=c("H40", "H41", "H44", "H42", "H43"))

ggbarplot(mel, x = "RNA_oligo", y = "value",
          ylab= "Kd (nM)", xlab = "RNA_oligonucleotide", color = "RNA_oligo", fill = "RNA_oligo",
          add = c("mean_sd"), palette = "jco",
          position = position_dodge(.8)) +
  geom_jitter(aes(RNA_oligo, value, fill = RNA_oligo), shape = 21, color = "black",  position = position_jitterdodge(jitter.height = -2, jitter.width = 1.5))

# the rest of fig3
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

##7mer crosslinked
mer<-dat

mer$count<-ifelse(is.na(mer$localization_cat), NA,
                  ifelse(mer$norm_tc_num1>2 & mer$norm_tc_num2>2 & mer$TPM_transcript>=10 , 1, NA))
mer<-subset(mer, !is.na(count))
mer$tc_seq_pos1_start<-ifelse(mer$tc_start-5<1, NA, mer$tc_start-5)
mer$tc_seq_pos1_stop<-ifelse(mer$tc_start+1>mer$l_tr, NA, mer$tc_start+1)
mer$tc_seq_pos2_start<-ifelse(mer$tc_start-4<1, NA, mer$tc_start-4)
mer$tc_seq_pos2_stop<-ifelse(mer$tc_start+2>mer$l_tr, NA, mer$tc_start+2)
mer$tc_seq_pos3_start<-ifelse(mer$tc_start-3<1, NA, mer$tc_start-3)
mer$tc_seq_pos3_stop<-ifelse(mer$tc_start+3>mer$l_tr, NA, mer$tc_start+3)
mer$tc_seq_pos4_start<-ifelse(mer$tc_start-2<1, NA, mer$tc_start-2)
mer$tc_seq_pos4_stop<-ifelse(mer$tc_start+4>mer$l_tr, NA, mer$tc_start+4)
mer$tc_seq_pos5_start<-ifelse(mer$tc_start-1<1, NA, mer$tc_start-1)
mer$tc_seq_pos5_stop<-ifelse(mer$tc_start+5>mer$l_tr, NA, mer$tc_start+5)
mer$tc_seq_pos6_start<-ifelse(mer$tc_start<1, NA, mer$tc_start)
mer$tc_seq_pos6_stop<-ifelse(mer$tc_start+6>mer$l_tr, NA, mer$tc_start+6)
mer$tc_seq_pos7_start<-ifelse(mer$tc_start+1<1, NA, mer$tc_start+1)
mer$tc_seq_pos7_stop<-ifelse(mer$tc_start+7>mer$l_tr, NA, mer$tc_start+7)



mer$tc_seq_pos1<-toupper(as.character(Biostrings::subseq(sub[as.character(mer$transcript_id)],
                                                         start=mer$tc_seq_pos1_start,end=mer$tc_seq_pos1_stop)))
mer$tc_seq_pos2<-toupper(as.character(Biostrings::subseq(sub[as.character(mer$transcript_id)],
                                                         start=mer$tc_seq_pos2_start,end=mer$tc_seq_pos2_stop)))
mer$tc_seq_pos3<-toupper(as.character(Biostrings::subseq(sub[as.character(mer$transcript_id)],
                                                         start=mer$tc_seq_pos3_start,end=mer$tc_seq_pos3_stop)))
mer$tc_seq_pos4<-toupper(as.character(Biostrings::subseq(sub[as.character(mer$transcript_id)],
                                                         start=mer$tc_seq_pos4_start,end=mer$tc_seq_pos4_stop)))
mer$tc_seq_pos5<-toupper(as.character(Biostrings::subseq(sub[as.character(mer$transcript_id)],
                                                         start=mer$tc_seq_pos5_start,end=mer$tc_seq_pos5_stop)))
mer$tc_seq_pos6<-toupper(as.character(Biostrings::subseq(sub[as.character(mer$transcript_id)],
                                                         start=mer$tc_seq_pos6_start,end=mer$tc_seq_pos6_stop)))
mer$tc_seq_pos7<-toupper(as.character(Biostrings::subseq(sub[as.character(mer$transcript_id)],
                                                         start=mer$tc_seq_pos7_start,end=mer$tc_seq_pos7_stop)))

mer<-melt(mer, measure.vars = colnames(mer)[48:54], id.vars = colnames(mer)[1:33], value.name = "seqs", variable.name = "tc_seq_pos")

countmers<-aggregate(count~seqs+tc_region+localization_cat, sum, data=mer)
# countmers$seqs<-gsub("T", "U", countmers$seqs)
countmers$frac_tot<-countmers$count/sum(countmers$count)
agr<-aggregate(frac_tot~seqs, sum, data=countmers)
agr<-agr[order(agr$frac_tot, decreasing = T),]
agr<-agr[1:10, "seqs"]
rows<-which(countmers$seqs %in% agr)
countmers_plot<-countmers[rows,]
countmers_plot$seqs<-factor(countmers_plot$seqs, levels=rev(agr))
countmers_plot$localization_cat<-factor(countmers_plot$localization_cat, levels=c("membrane", "cytosolic"))

#fig3a
ggplot(subset(countmers_plot, tc_region!="utr5"), aes(gsub("T","U",seqs), frac_tot, fill=tc_region))+geom_bar(stat="identity", position="dodge",)+coord_flip()+facet_wrap(~localization_cat)+
  theme(text = element_text(size=8))+scale_fill_manual(values = c("dodgerblue3", "orange2"))

sd<-subset(countmers_plot, tc_region!="utr5")
sd$seqs<-gsub("T", "U", sd$seqs)
#write.table(sd, "source_data/fig3a.txt", quote=F, sep="\t", row.names=F)


#z-scores. due to size limitations we only provide files up to 8-mer.

files<-list.files("data/", pattern="whole", full.names = T)
zs<-lapply(files, read.delim, header=T)
names(zs)<-gsub(".*;","",gsub("-.*", "", sub("-",";", files)))

zs<-lapply(zs, function(x) x[order(x$labs, decreasing=T),])
zs<-lapply(zs, subset, !duplicated(kmer))

diffs_CdsMem_Utr3Cyt<-sapply(zs, function(x) x$freqCds_mem - x$freqUtr3_cyt)
zscores_CdsMem_Utr3Cyt<-sapply(diffs_CdsMem_Utr3Cyt, function(x) scale(x, center=T,scale=T)[,1])
diffs_CdsMem_CdsCyt<-sapply(zs, function(x) x$freqCds_mem - x$freqCds_cyt)
zscores_CdsMem_CdsCyt<-sapply(diffs_CdsMem_CdsCyt, function(x) scale(x, center=T,scale=T)[,1])
diffs_CdsUtr3Mem_CdsUtr3Cyt<-sapply(zs, function(x) x$freqCdsUtr3_mem - x$freqCdsUtr3_cyt)
zscores_CdsUtr3Mem_CdsUtr3Cyt<-sapply(diffs_CdsUtr3Mem_CdsUtr3Cyt, function(x) scale(x, center=T,scale=T)[,1])
diffs_Utr3Mem_Utr3Cyt<-sapply(zs, function(x) x$freqUtr3_mem - x$freqUtr3_cyt)
zscores_Utr3Mem_Utr3Cyt<-sapply(diffs_Utr3Mem_Utr3Cyt, function(x) scale(x, center=T,scale=T)[,1])
diffs_CdsCyt_Utr3Cyt<-sapply(zs, function(x) x$freqCds_cyt - x$freqUtr3_cyt)
zscores_CdsCyt_Utr3Cyt<-sapply(diffs_CdsCyt_Utr3Cyt, function(x) scale(x, center=T,scale=T)[,1])
diffs_CdsMem_Utr3Mem<-sapply(zs, function(x) x$freqCds_mem - x$freqUtr3_mem)
zscores_CdsMem_Utr3Mem<-sapply(diffs_CdsMem_Utr3Mem, function(x) scale(x, center=T,scale=T)[,1])
diffs_CdsCyt_Utr3Mem<-sapply(zs, function(x) x$freqCds_cyt - x$freqUtr3_mem)
zscores_CdsCyt_Utr3Mem<-sapply(diffs_CdsCyt_Utr3Mem, function(x) scale(x, center=T,scale=T)[,1])


zs<-Map(cbind, zs, diffs_CdsMem_Utr3Cyt=diffs_CdsMem_Utr3Cyt, zscores_CdsMem_Utr3Cyt=zscores_CdsMem_Utr3Cyt,
        diffs_CdsMem_CdsCyt=diffs_CdsMem_CdsCyt, zscores_CdsMem_CdsCyt=zscores_CdsMem_CdsCyt,
        diffs_CdsUtr3Mem_CdsUtr3Cyt=diffs_CdsUtr3Mem_CdsUtr3Cyt, zscores_CdsUtr3Mem_CdsUtr3Cyt=zscores_CdsUtr3Mem_CdsUtr3Cyt,
        diffs_Utr3Mem_Utr3Cyt=diffs_Utr3Mem_Utr3Cyt, zscores_Utr3Mem_Utr3Cyt=zscores_Utr3Mem_Utr3Cyt,
        diffs_CdsCyt_Utr3Cyt=diffs_CdsCyt_Utr3Cyt, zscores_CdsCyt_Utr3Cyt=zscores_CdsCyt_Utr3Cyt,
        diffs_CdsMem_Utr3Mem=diffs_CdsMem_Utr3Mem, zscores_CdsMem_Utr3Mem=zscores_CdsMem_Utr3Mem,
        diffs_CdsCyt_Utr3Mem=diffs_CdsCyt_Utr3Mem, zscores_CdsCyt_Utr3Mem=zscores_CdsCyt_Utr3Mem)

zs<-do.call("rbind", zs)

zs$length<-gsub("\\..*", "", row.names(zs))
zs$type=ifelse(is.na(zs$labs), "other", "top40")
zs$type<-factor(zs$type, levels=c("top40", "other"))
ggplot(zs, aes(zscores_CdsMem_Utr3Cyt, colour=length))+geom_density()+facet_wrap(~type, scales="free")+coord_cartesian(xlim=c(-5,12))
ggplot(zs, aes(zscores_CdsMem_CdsCyt, colour=length))+geom_density()+facet_wrap(~type, scales="free")+coord_cartesian(xlim=c(-2,10))
ggplot(zs, aes(zscores_CdsUtr3Mem_CdsUtr3Cyt, colour=length))+geom_density()+facet_wrap(~type, scales="free")+coord_cartesian(xlim=c(-5,10))
ggplot(zs, aes(zscores_Utr3Mem_Utr3Cyt, colour=length))+geom_density()+facet_wrap(~type, scales="free")+coord_cartesian(xlim=c(-5,4))
ggplot(zs, aes(zscores_CdsCyt_Utr3Cyt, colour=length))+geom_density()+facet_wrap(~type, scales="free")+coord_cartesian(xlim=c(-5,1.25))

ggplot(zs, aes(zscores_CdsMem_CdsCyt, fill=length))+geom_histogram(bins=500)+facet_wrap(~type+length, scales="free_y", ncol=6)+
  geom_vline(xintercept=0, lty=2, colour="grey")+coord_cartesian(xlim=c(-3,10))

#fig3c
p<-ggplot(zs, aes(zscores_CdsUtr3Mem_CdsUtr3Cyt, colour=length))+stat_ecdf()+facet_wrap(~type, scales="free")+coord_cartesian(xlim=c(-3,15))
print(p)

sd<-subset(zs, type=="top40",select=c("kmer", "type", "length", colnames(zs)[grepl("zscores", colnames(zs))]))
ggplot(sd, aes(zscores_CdsUtr3Mem_CdsUtr3Cyt, colour=length))+stat_ecdf()+facet_wrap(~type, scales="free")+coord_cartesian(xlim=c(-3,15))

# write.table(sd, "source_data/fig3c.txt", quote=F, sep="\t", row.names=F)

#figS3c top
ggplot(subset(zs, type=="top40"), aes(zscores_CdsMem_CdsCyt, colour=length))+stat_ecdf()+facet_wrap(~type, scales="free")+
  coord_cartesian(xlim=c(-3,20))
#figS3c bottom
ggplot(subset(zs, type=="top40"), aes(zscores_CdsCyt_Utr3Cyt, colour=length))+stat_ecdf()+facet_wrap(~type, scales="free")+
  coord_cartesian(xlim=c(-3,15))
#figS3d
ggplot(subset(zs, type=="top40"), aes(zscores_CdsMem_Utr3Cyt, colour=length))+stat_ecdf()+facet_wrap(~type, scales="free")+
  coord_cartesian(xlim=c(-3,15))

wilcox.test(subset(zs, type=="top40" & length=="4mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="5mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="4mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="6mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="4mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="7mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="4mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="8mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
# wilcox.test(subset(zs, type=="top40" & length=="4mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
#             subset(zs, type=="top40" & length=="9mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
# wilcox.test(subset(zs, type=="top40" & length=="4mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
#             subset(zs, type=="top40" & length=="10mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
# wilcox.test(subset(zs, type=="top40" & length=="4mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
#             subset(zs, type=="top40" & length=="11mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
# wilcox.test(subset(zs, type=="top40" & length=="4mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
#             subset(zs, type=="top40" & length=="12mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="5mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="6mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="5mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="7mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="5mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="8mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
# wilcox.test(subset(zs, type=="top40" & length=="5mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
#             subset(zs, type=="top40" & length=="9mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
# wilcox.test(subset(zs, type=="top40" & length=="5mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
#             subset(zs, type=="top40" & length=="10mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
# wilcox.test(subset(zs, type=="top40" & length=="5mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
#             subset(zs, type=="top40" & length=="11mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
# wilcox.test(subset(zs, type=="top40" & length=="5mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
#             subset(zs, type=="top40" & length=="12mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="6mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="7mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="6mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="8mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
# wilcox.test(subset(zs, type=="top40" & length=="6mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
#             subset(zs, type=="top40" & length=="9mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
# wilcox.test(subset(zs, type=="top40" & length=="6mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
#             subset(zs, type=="top40" & length=="10mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
# wilcox.test(subset(zs, type=="top40" & length=="6mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
#             subset(zs, type=="top40" & length=="11mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
# wilcox.test(subset(zs, type=="top40" & length=="6mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
#             subset(zs, type=="top40" & length=="12mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value

wilcox.test(subset(zs, type=="top40" & length=="7mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="8mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
# wilcox.test(subset(zs, type=="top40" & length=="7mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
#             subset(zs, type=="top40" & length=="9mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
# wilcox.test(subset(zs, type=="top40" & length=="7mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
#             subset(zs, type=="top40" & length=="10mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
# wilcox.test(subset(zs, type=="top40" & length=="7mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
#             subset(zs, type=="top40" & length=="11mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
# wilcox.test(subset(zs, type=="top40" & length=="7mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
#             subset(zs, type=="top40" & length=="12mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
# 
# wilcox.test(subset(zs, type=="top40" & length=="8mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
#             subset(zs, type=="top40" & length=="9mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
# wilcox.test(subset(zs, type=="top40" & length=="8mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
#             subset(zs, type=="top40" & length=="10mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
# wilcox.test(subset(zs, type=="top40" & length=="8mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
#             subset(zs, type=="top40" & length=="11mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
# wilcox.test(subset(zs, type=="top40" & length=="8mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
#             subset(zs, type=="top40" & length=="12mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
# 
# wilcox.test(subset(zs, type=="top40" & length=="9mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
#             subset(zs, type=="top40" & length=="10mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
# wilcox.test(subset(zs, type=="top40" & length=="9mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
#             subset(zs, type=="top40" & length=="11mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
# wilcox.test(subset(zs, type=="top40" & length=="9mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
#             subset(zs, type=="top40" & length=="12mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
# 
# wilcox.test(subset(zs, type=="top40" & length=="10mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
#             subset(zs, type=="top40" & length=="11mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
# wilcox.test(subset(zs, type=="top40" & length=="10mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
#             subset(zs, type=="top40" & length=="12mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
# 
# wilcox.test(subset(zs, type=="top40" & length=="11mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
#             subset(zs, type=="top40" & length=="12mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value

#fig3d: p values not done yet

setwd("E:/work//hdlbp/kmers_mem_cyt/")
files<-list.files(getwd(), pattern="pval")
ps<-lapply(files, read.delim, header=F)
names(ps)<-gsub("\\..*","", gsub("pval_", "", files))
kmer<-c(rep(4, 8), rep(5, 7), rep(6, 6), rep(7, 5), rep(8, 4), rep(9, 3), rep(10, 2), rep(11, 1))
kmer2<-c(5:12, 6:12, 7:12, 8:12, 9:12, 10:12, 11:12,12)
ps<-Map(cbind, ps, kmer=rep(list(kmer),length(ps)))
ps<-Map(cbind, ps, kmer2=rep(list(kmer2),length(ps)))
ones<-rep(1, length(unique(ps[[1]]$kmer))+1)
ones<-data.frame(V1=ones, kmer=4:12, kmer2=4:12)
ps<-Map(rbind, ps, kmer2=rep(list(ones),length(ps)))
ps<-do.call("rbind", ps)
ps$comparison<-gsub("\\..*", "", row.names(ps))
colnames(ps)[1]<-"p"
ps$comparison<-factor(ps$comparison, levels=c("CdsUtr3Mem_CdsUtr3Cyt", "CdsMem_CdsCyt","CdsMem_Utr3Cyt","CdsCyt_Utr3Cyt", "Utr3Mem_Utr3Cyt" ))
ggplot(ps, aes(factor(kmer), factor(kmer2)))+geom_tile(aes(fill=-log10(p)), na.rm = T)+facet_wrap(~comparison, nrow=1)
ggplot(ps, aes(factor(kmer2), factor(kmer)))+geom_tile(aes(fill=-log10(p)), na.rm = T)+facet_wrap(~comparison, nrow=1)

p_CdsMem_CdsCyt<-sapply(zs, function(x) wilcox.test(subset(x, !is.na(labs), select="freqCds_mem")[,1],
                                                    subset(x, !is.na(labs), select="freqCds_cyt")[,1] )$p.value)
p_CdsMem_CdsCyt_other<-sapply(zs, function(x) wilcox.test(subset(x, is.na(labs), select="freqCds_mem")[,1],
                                                          subset(x, is.na(labs), select="freqCds_cyt")[,1] )$p.value)
p_CdsUtr3Mem_CdsUtr3Cyt<-sapply(zs, function(x) wilcox.test(subset(x, !is.na(labs), select="freqCdsUtr3_mem")[,1],
                                                            subset(x, !is.na(labs), select="freqCdsUtr3_cyt")[,1] )$p.value)
p_CdsUtr3Mem_CdsUtr3Cyt_other<-sapply(zs, function(x) wilcox.test(subset(x, is.na(labs), select="freqCdsUtr3_mem")[,1],
                                                                  subset(x, is.na(labs), select="freqCdsUtr3_cyt")[,1] )$p.value)

p_CdsMem_Utr3Cyt<-sapply(zs, function(x) wilcox.test(subset(x, !is.na(labs), select="freqCds_mem")[,1],
                                                     subset(x, is.na(labs), select="freqUtr3_cyt")[,1] )$p.value)
p_CdsMem_Utr3Cyt<-sapply(zs, function(x) wilcox.test(subset(x, !is.na(labs), select="diffs_CdsMem_Utr3Cyt")[,1],
                                                     subset(x, is.na(labs), select="diffs_CdsMem_Utr3Cyt")[,1] )$p.value)
p_CdsUtr3Mem_CdsUtr3Cyt<-sapply(zs, function(x) wilcox.test(subset(x, !is.na(labs), select="diffs_CdsUtr3Mem_CdsUtr3Cyt")[,1],
                                                            subset(x, is.na(labs), select="diffs_CdsUtr3Mem_CdsUtr3Cyt")[,1] )$p.value)

dats<-zs[[7]]
ggplot(dats, aes(diffs_CdsMem_CdsCyt, freqCds_mem+freqCds_cyt,colour=labs))+geom_point(shape=1)+geom_text(aes(label=labs), size=3)+theme(legend.position = "none")
ggplot(dats, aes(diffs_CdsMem_Utr3Cyt, freqCds_mem+freqUtr3_cyt,colour=labs))+geom_point(shape=1)+geom_text(aes(label=labs), size=3)+theme(legend.position = "none")+coord_cartesian(xlim=c(-0.001,0.001), ylim=c(0,0.001))


ggplot(dats, aes(freqCds_mem, freqCds_cyt,colour=labs))+geom_point(shape=1, size=1)+geom_text(aes(label=labs), size=3)+geom_abline()+theme(legend.position = "none")
ggplot(dats, aes(freqCds_mem, freqUtr3_cyt,colour=labs))+geom_point(shape=1, size=1)+geom_text(aes(label=labs), size=3)+geom_abline()
ggplot(dats, aes(freqCdsUtr3_mem, freqCdsUtr3_cyt,colour=labs))+geom_point(shape=1, size=1)+geom_text(aes(label=labs), size=3)+geom_abline()+theme(legend.position = "none")+coord_cartesian(xlim=c(0,0.0006), ylim=c(0,0.0006))

dats<-zs[[9]]
ggplot(dats, aes(diffs_CdsMem_CdsCyt, freqCds_mem+freqCds_cyt,colour=labs))+geom_point(shape=1)+geom_text(aes(label=labs), size=3)+theme(legend.position = "none")
ggplot(dats, aes(diffs_CdsMem_Utr3Cyt, freqCds_mem+freqUtr3_cyt,colour=labs))+geom_point(shape=1)+geom_text(aes(label=labs), size=3)+theme(legend.position = "none")+coord_cartesian(xlim=c(-0.001,0.001), ylim=c(0,0.001))


ggplot(dats, aes(freqCds_mem, freqCds_cyt,colour=labs))+geom_point(shape=1, size=1)+geom_text(aes(label=labs), size=3)+geom_abline()+theme(legend.position = "none")
ggplot(dats, aes(freqCds_mem, freqUtr3_cyt,colour=labs))+geom_point(shape=1, size=1)+geom_text(aes(label=labs), size=3)+geom_abline()

#fig3b
##7mer crosslinked with motif
mer<-dat

mer$count<-ifelse(is.na(mer$localization_cat), NA,
                  ifelse(mer$norm_tc_num1>2 & mer$norm_tc_num2>2 & mer$TPM_transcript>=10 , 1, NA))
mer<-subset(mer, !is.na(count))
mer$seq<-toupper(mer$seq)
countmers<-aggregate(count~seq+tc_region+localization_cat, sum, data=mer)

countmers$frac_tot<-countmers$count/sum(countmers$count)
agr<-aggregate(cbind(count,frac_tot)~seq+localization_cat+tc_region, sum, data=countmers)
agr<-subset(agr, localization_cat=="membrane" & tc_region=="cds")
agr<-agr[order(agr$frac_tot, decreasing = T),]
agr<-agr[1:5, "seq"]
rows<-which(countmers$seq %in% agr)
countmers_plot<-countmers[rows,]
countmers_plot$seq<-factor(countmers_plot$seq, levels=rev(agr))
# ggplot(subset(countmers_plot, tc_region!="utr5"), aes(seq, frac_tot, fill=tc_region))+geom_bar(stat="identity", position="dodge",)+coord_flip()+facet_wrap(~localization_cat)+
#   theme(text = element_text(size=8))

library(ggseqlogo)
p<-ggplot(subset(countmers_plot, localization_cat=="membrane" & tc_region=="cds"))+geom_logo(gsub("T","U",as.character(countmers_plot$seq)), method="bits", seq_type = "rna")+theme_logo()
print(p)

sd<-unique(countmers_plot$seq)
sd



#begin multivalency fig3e fig3f
pfgr<-read.delim("data/h_values_transcriptome_with_annotation_kmer.txt", header=T)
pfgr$hval<-1
pfgr$dup_id<-paste0(pfgr$second.X.seqnames,"_", pfgr$second.mcols.norm_tc_num1,
                    "_", pfgr$second.mcols.norm_tc_num2,"_", pfgr$second.mcols.localization_cat,
                    "_", pfgr$second.mcols.tc_region,"_", pfgr$distance,"_",pfgr$kmer_pos,"_",pfgr$kmer,"_",pfgr$kmer)
pfgr<-subset(pfgr, !duplicated(dup_id))
# hs<-aggregate(hval~second.X.seqnames+kmer+second.mcols.localization_cat+second.mcols.tc_region, data=pfgr, sum)
hs<-aggregate(hval~second.X.seqnames+kmer+second.mcols.localization_cat+second.mcols.tc_region, 
              data=subset(pfgr, abs(distance)>4), sum)

mean_hs<-aggregate(hval~kmer+second.mcols.localization_cat+second.mcols.tc_region, data=hs, mean)
max_hs<-aggregate(hval~kmer+second.mcols.localization_cat+second.mcols.tc_region, data=hs, max)
sum_hs<-aggregate(hval~kmer+second.mcols.localization_cat+second.mcols.tc_region, data=hs, sum)


mean_hs_all<-aggregate(hval~kmer, data=hs, mean)
mean_hs_all<-mean_hs_all[order(mean_hs_all$hval, decreasing = T),]
mean_hs_all<-mean_hs_all[1:10,]
filter<-which(mean_hs$kmer %in% mean_hs_all$kmer)
mean_hs_plot<-mean_hs[filter,]
ggplot(mean_hs_plot, aes(kmer,hval, fill=second.mcols.tc_region))+
  geom_bar(stat="identity", position="dodge")+facet_wrap(~second.mcols.localization_cat)+coord_flip()

max_hs_all<-aggregate(hval~kmer, data=hs, max)
max_hs_all<-max_hs_all[order(max_hs_all$hval, decreasing = T),]
max_hs_all<-max_hs_all[1:10,]
filter<-which(max_hs$kmer %in% max_hs_all$kmer)
max_hs_plot<-max_hs[filter,]
ggplot(max_hs_plot, aes(kmer,hval, fill=second.mcols.tc_region))+
  geom_bar(stat="identity", position="dodge")+facet_wrap(~second.mcols.localization_cat)+coord_flip()


sum_hs_all<-aggregate(hval~kmer, data=hs, sum)
sum_hs_all<-sum_hs_all[order(sum_hs_all$hval, decreasing = T),]
sum_hs_all<-sum_hs_all[1:10,]
filter<-which(sum_hs$kmer %in% sum_hs_all$kmer)
sum_hs_plot<-sum_hs[filter,]
sum_hs_plot$kmer<-factor(sum_hs_plot$kmer, levels=rev(sum_hs_all$kmer))
ggplot(subset(sum_hs_plot, second.mcols.tc_region!="utr5"), aes(gsub("T","U",kmer),hval, fill=second.mcols.tc_region))+
  geom_bar(stat="identity", position="dodge")+facet_wrap(~second.mcols.localization_cat)+coord_flip()+scale_fill_manual(values=c("dodgerblue4", "orange3"))


freq_hs<-aggregate(hval~kmer+second.mcols.localization_cat+second.mcols.tc_region, data=hs, sum)
freq_hs$freq<-freq_hs$hval/sum(freq_hs$hval)
freq_hs_all<-aggregate(freq~kmer, data=freq_hs, sum)
freq_hs_all<-freq_hs_all[order(freq_hs_all$freq, decreasing = T),]
freq_hs_all<-freq_hs_all[1:10,]
filter<-which(freq_hs$kmer %in% freq_hs_all$kmer)
freq_hs_plot<-freq_hs[filter,]
freq_hs_plot$kmer<-factor(freq_hs_plot$kmer, levels=rev(unique(freq_hs_all$kmer)))
freq_hs_plot$second.mcols.localization_cat<-factor(freq_hs_plot$second.mcols.localization_cat, levels=c("membrane", "cytosolic"))

#figS3f right panel
ggplot(subset(freq_hs_plot, second.mcols.tc_region!="utr5"), aes(kmer,freq, fill=second.mcols.tc_region))+
  geom_bar(stat="identity", position="dodge")+facet_wrap(~second.mcols.localization_cat)+coord_flip()+scale_fill_manual(values=c("dodgerblue4", "orange3"))
sd<-subset(freq_hs_plot, second.mcols.tc_region!="utr5")
# write.table(sd, "source_data/figs3f_right.txt", quote=F, sep="\t", row.names=F)


hs$id<-paste0(hs$second.X.seqnames, "_", hs$kmer)
hs$log2hval<-ifelse(hs$hval<=0, NA, log2(hs$hval))


hs_plot<-which(hs$kmer %in% freq_hs_plot$kmer)
hs_plot<-hs[hs_plot,]
hs_plot$kmer<-factor(hs_plot$kmer, levels=rev(freq_hs_all$kmer))

ggplot(subset(hs_plot, second.mcols.tc_region!="utr5"), aes(kmer,log2(hval), fill=second.mcols.tc_region)) +
  geom_boxplot()+facet_wrap(~second.mcols.localization_cat)+coord_flip()
ggplot(subset(hs_plot, second.mcols.tc_region!="utr5"), aes(kmer,log2(hval), fill=second.mcols.localization_cat)) +
  geom_boxplot()+facet_wrap(~second.mcols.tc_region)+coord_flip()+geom_hline(yintercept=3, lty=2, colour="grey")+scale_fill_manual(values=c("dodgerblue3", "orange3"))


hs_plot<-subset(hs_plot, select=c("id", "hval"))
tab<-aggregate(hval~second.X.seqnames, data=hs_plot,sum)
colnames(tab)[1]<-"transcript"
#write.table(tab, "hval_perGene.txt", quote=F, sep="\t", row.names=F)
pfgr$id<-paste0(pfgr$second.X.seqnames, "_", pfgr$kmer)
pfgr<-pfgr[,-9]
pfgr<-merge(pfgr, hs_plot, by="id", all.x=T)
pfgr<-subset(pfgr, !duplicated(id))

sub_pfgr<-which(pfgr$kmer %in% freq_hs_plot$kmer)
sub_pfgr<-pfgr[sub_pfgr,]
sub_pfgr$hval<-ifelse(is.na(sub_pfgr$hval), 0, sub_pfgr$hval)

min_tc<-aggregate(cbind( second.mcols.norm_tc_num1, 
                         second.mcols.norm_tc_num2)~second.X.seqnames, data=sub_pfgr, sum)
min_tc<-subset(min_tc, second.mcols.norm_tc_num1>=2 &
                 second.mcols.norm_tc_num2>=2)
min_tc<-which(sub_pfgr$second.X.seqnames %in% min_tc$second.X.seqnames)
sub_pfgr<-sub_pfgr[min_tc,]

allhs<-aggregate(hval~second.X.seqnames+kmer, sum, data=sub_pfgr)

quants<-quantile(allhs$hval, probs=c(0,1/5,2/5,3/5,4/5,1))
quants<-quantile(subset(sub_pfgr,second.mcols.norm_tc_num1>0 & !is.na(second.mcols.norm_tc_num1)&
                          second.mcols.norm_tc_num2>0 & !is.na(second.mcols.norm_tc_num2), select="hval")[,1], probs=c(0,1/5,2/5,3/5,4/5,1))
sub_pfgr$hval_cat<-ifelse(sub_pfgr$hval>quants[5], "hval_hi", 
                          ifelse(sub_pfgr$hval>quants[4], "hval_midhi", 
                                 ifelse(sub_pfgr$hval>quants[3], "hval_midmid",
                                        ifelse(sub_pfgr$hval>quants[2], "hval_midlo",
                                               ifelse(sub_pfgr$hval>=0, "hval_lo", NA)))))
sub_pfgr$hval_cat<-factor(sub_pfgr$hval_cat, levels=c("hval_hi", "hval_midhi", "hval_midmid", "hval_midlo", "hval_lo"))

ggplot(sub_pfgr, aes(hval_cat))+geom_bar()
ggplot(sub_pfgr , aes(hval_cat, log2(second.mcols.norm_tc_num1), fill=second.mcols.localization_cat))+geom_boxplot()+facet_wrap(~second.mcols.tc_region)
ggplot(sub_pfgr , aes(hval_cat, log2(second.mcols.norm_tc_num2), fill=second.mcols.localization_cat))+geom_boxplot()+facet_wrap(~second.mcols.tc_region)
ggplot(sub_pfgr , aes(hval_cat, log2(second.mcols.norm_tc_num1), fill=hval_cat))+geom_violin(scale = "count")+geom_boxplot(width=0.1, outlier.shape = NA) +theme(axis.text.x = element_text(angle = 45, hjust=1))+
  scale_fill_brewer(palette=4)
ggplot(sub_pfgr , aes(hval_cat, log2(second.mcols.norm_tc_num2), fill=hval_cat))+geom_violin(scale = "count")+geom_boxplot(width=0.1, outlier.shape = NA) +theme(axis.text.x = element_text(angle = 45, hjust=1))+
  scale_fill_brewer(palette=4)

#fig3e
ggplot(sub_pfgr , aes(hval_cat, log2(rowMeans(sub_pfgr[,c("second.mcols.norm_tc_num1","second.mcols.norm_tc_num2")])), fill=hval_cat))+geom_violin(scale = "count")+geom_boxplot(width=0.1, outlier.shape = NA) +theme(axis.text.x = element_text(angle = 45, hjust=1))+
  scale_fill_brewer(palette=4)


summary(subset(sub_pfgr, !is.na(log2(rowMeans(sub_pfgr[,c("second.mcols.norm_tc_num1","second.mcols.norm_tc_num2")]))))$hval_cat)
sd<-subset(sub_pfgr, !is.na(log2(rowMeans(sub_pfgr[,c("second.mcols.norm_tc_num1","second.mcols.norm_tc_num2")]))),
           select=c("second.X.seqnames", "second.mcols.norm_tc_num1", "second.mcols.norm_tc_num2",
                    "second.mcols.localization_cat", "second.mcols.tc_region", "kmer_pos","kmer","hval", "hval_cat"))
# write.table(sd,"source_data/fig3e.txt", quote=F, sep="\t", row.names=F)

##most crosslinked 4mers


# pfgr<-melt(dfgr, measure.vars = colnames(dfgr)[grepl("pos[0-9]$",colnames(dfgr))],
# id.vars = c("second.X.seqnames","second.mcols.norm_tc_num1", "second.mcols.norm_tc_num2", "second.mcols.localization_cat", "second.mcols.tc_region", "distance"),
# value.name = "kmer", variable.name = "kmer_pos")
pfgr<-read.delim("data/h_values_transcriptome_with_annotation_kmer.txt", header=T)
pfgr$hval<-1
pfgr$dup_id<-paste0(pfgr$second.X.seqnames,"_", pfgr$second.mcols.norm_tc_num1,
                    "_", pfgr$second.mcols.norm_tc_num2,"_", pfgr$second.mcols.localization_cat,
                    "_", pfgr$second.mcols.tc_region,"_", pfgr$distance,"_",pfgr$kmer_pos,"_",pfgr$kmer,"_",pfgr$kmer)
pfgr<-subset(pfgr, !duplicated(dup_id))
# hs<-aggregate(hval~second.X.seqnames+kmer+second.mcols.localization_cat+second.mcols.tc_region, data=pfgr, sum)
hs<-aggregate(hval~second.X.seqnames+kmer+second.mcols.localization_cat+second.mcols.tc_region, 
              data=subset(pfgr, abs(distance)>4), sum)

mean_hs<-aggregate(hval~kmer+second.mcols.localization_cat+second.mcols.tc_region, data=hs, mean)
max_hs<-aggregate(hval~kmer+second.mcols.localization_cat+second.mcols.tc_region, data=hs, max)

xlinked<-c("TTCT","CTTC","TCTT","TTCC","TCCT","CTCT","TTTC","CTTT","TCTC","TTTT") #most crosslinked kmers
# xlinked<-unique(freq_hs_plot$kmer)
# xlinked<-"CTTC"

filter<-which(mean_hs$kmer %in% xlinked)
mean_hs_plot<-mean_hs[filter,]
ggplot(mean_hs_plot, aes(kmer,hval, fill=second.mcols.tc_region))+
  geom_bar(stat="identity", position="dodge")+facet_wrap(~second.mcols.localization_cat)+coord_flip()

filter<-which(max_hs$kmer %in% xlinked)
max_hs_plot<-max_hs[filter,]
ggplot(max_hs_plot, aes(kmer,hval, fill=second.mcols.tc_region))+
  geom_bar(stat="identity", position="dodge")+facet_wrap(~second.mcols.localization_cat)+coord_flip()

hs$id<-paste0(hs$second.X.seqnames, "_", hs$kmer)
hs$log2hval<-ifelse(hs$hval<=0, NA, log2(hs$hval))

hs_plot<-which(hs$kmer %in% xlinked)
hs_plot<-hs[hs_plot,]
ggplot(subset(hs_plot, second.mcols.tc_region!="utr5"), aes(kmer,log2hval, fill=second.mcols.tc_region)) +
  geom_boxplot()+facet_wrap(~second.mcols.localization_cat)+coord_flip()
ggplot(subset(hs_plot, second.mcols.tc_region!="utr5"), aes(kmer,log2hval, fill=second.mcols.localization_cat)) + geom_boxplot()+facet_wrap(~second.mcols.tc_region)+coord_flip()

hs_plot<-subset(hs_plot, select=c("id", "hval"))
pfgr$id<-paste0(pfgr$second.X.seqnames, "_", pfgr$kmer)
pfgr<-pfgr[,-9]
pfgr<-merge(pfgr, hs_plot, by="id", all.x=T)
pfgr<-subset(pfgr, !duplicated(id))

sub_pfgr<-which(pfgr$kmer %in% xlinked)
sub_pfgr<-pfgr[sub_pfgr,]
sub_pfgr$hval<-ifelse(is.na(sub_pfgr$hval), 0, sub_pfgr$hval)

min_tc<-aggregate(cbind( second.mcols.norm_tc_num1, 
                         second.mcols.norm_tc_num2)~second.X.seqnames, data=sub_pfgr, sum)
min_tc<-subset(min_tc, second.mcols.norm_tc_num1>=2 &
                 second.mcols.norm_tc_num2>=2)
min_tc<-which(sub_pfgr$second.X.seqnames %in% min_tc$second.X.seqnames)
sub_pfgr<-sub_pfgr[min_tc,]

allhs<-aggregate(hval~second.X.seqnames+kmer, sum, data=sub_pfgr)

quants<-quantile(allhs$hval, probs=c(0,1/5,2/5,3/5,4/5,1))
quants<-quantile(sub_pfgr$hval, probs=c(0,1/5,2/5,3/5,4/5,1))

sub_pfgr$hval_cat<-ifelse(sub_pfgr$hval>quants[5], "hval_hi", 
                          ifelse(sub_pfgr$hval>quants[4], "hval_midhi", 
                                 ifelse(sub_pfgr$hval>quants[3], "hval_midmid",
                                        ifelse(sub_pfgr$hval>quants[2], "hval_midlo",
                                               ifelse(sub_pfgr$hval>=0, "hval_lo", NA)))))
sub_pfgr$hval_cat<-factor(sub_pfgr$hval_cat, levels=c("hval_hi", "hval_midhi", "hval_midmid", "hval_midlo", "hval_lo"))

ggplot(sub_pfgr, aes(hval_cat))+geom_bar()
ggplot(sub_pfgr , aes(hval_cat, log2(second.mcols.norm_tc_num1), fill=second.mcols.localization_cat))+geom_boxplot()+facet_wrap(~second.mcols.tc_region)
ggplot(sub_pfgr , aes(hval_cat, log2(second.mcols.norm_tc_num2), fill=second.mcols.localization_cat))+geom_boxplot()+facet_wrap(~second.mcols.tc_region)
ggplot(sub_pfgr , aes(hval_cat, log2(second.mcols.norm_tc_num1), fill=hval_cat))+geom_boxplot()
ggplot(sub_pfgr , aes(hval_cat, log2(second.mcols.norm_tc_num2), fill=hval_cat))+geom_boxplot()
ggplot(sub_pfgr , aes(hval_cat, log2(second.mcols.norm_tc_num2), fill=hval_cat))+geom_violin()
ggplot(sub_pfgr, aes(distance))+geom_bar()


plot<-aggregate(second.mcols.norm_tc_num1~distance+
                  second.mcols.localization_cat+second.mcols.tc_region+hval_cat, data=sub_pfgr, function(x) sum(x)/sum(sub_pfgr[,"second.mcols.norm_tc_num1"])*100)


ggplot(plot, aes(hval_cat))+geom_bar()

plot$id<-paste0(plot$second.mcols.localization_cat,";",
                plot$second.mcols.tc_region, ";",
                plot$hval_cat)
library(reshape2)
fill<-dcast(plot, distance~id, value.var = "second.mcols.norm_tc_num1")
fill[is.na(fill)]<-0
plot<-melt(fill, measure.vars = colnames(fill)[2:ncol(fill)],
           id.vars="distance", value.name = "tc", variable.name = "id")
plot$localization_cat<-gsub(";.*","",plot$id)
plot$tc_region<-gsub(".*;","",gsub(";hval.*", "", plot$id))
plot$hval_cat<-gsub(".*;","",plot$id)
plot$hval_cat<-factor(plot$hval_cat, levels=c("hval_hi", "hval_midhi", "hval_midmid", "hval_midlo", "hval_lo"))
ggplot(plot , aes(hval_cat, log2(tc), fill=hval_cat))+geom_violin(scale = "count")+facet_wrap(~localization_cat+tc_region)


#fig3f replicate 1
ggplot(subset(plot, ((tc_region=="cds" & localization_cat=="membrane") | (tc_region=="utr3" & localization_cat=="cytosolic") )),
       aes(distance, tc, colour=hval_cat))+geom_line()+facet_wrap(~localization_cat+tc_region+hval_cat, scales="free_y", ncol=5)

sd<-subset(plot, ((tc_region=="cds" & localization_cat=="membrane") | (tc_region=="utr3" & localization_cat=="cytosolic") ))
# write.table(sd, "source_data/fig3f.txt", quote=F, sep="\t", row.names=F)

ggplot(subset(plot, tc_region!="utr5"),aes(distance, tc, colour=hval_cat))+geom_line()+facet_wrap(~localization_cat+tc_region+hval_cat, scales="free_y", ncol=5)



plot<-aggregate(second.mcols.norm_tc_num2~distance+
                  second.mcols.localization_cat+second.mcols.tc_region+hval_cat, data=sub_pfgr, function(x) sum(x)/sum(sub_pfgr[,"second.mcols.norm_tc_num2"])*100)


ggplot(plot, aes(hval_cat))+geom_bar()

plot$id<-paste0(plot$second.mcols.localization_cat,";",
                plot$second.mcols.tc_region, ";",
                plot$hval_cat)
library(reshape2)
fill<-dcast(plot, distance~id, value.var = "second.mcols.norm_tc_num2")
fill[is.na(fill)]<-0
plot<-melt(fill, measure.vars = colnames(fill)[2:ncol(fill)],
           id.vars="distance", value.name = "tc", variable.name = "id")
plot$localization_cat<-gsub(";.*","",plot$id)
plot$tc_region<-gsub(".*;","",gsub(";hval.*", "", plot$id))
plot$hval_cat<-gsub(".*;","",plot$id)
plot$hval_cat<-factor(plot$hval_cat, levels=c("hval_hi", "hval_midhi", "hval_midmid", "hval_midlo", "hval_lo"))
ggplot(plot , aes(hval_cat, log2(tc), fill=hval_cat))+geom_violin(scale = "count")+facet_wrap(~localization_cat+tc_region)

#fig3f replicate 2
ggplot(subset(plot, ((tc_region=="cds" & localization_cat=="membrane") | (tc_region=="utr3" & localization_cat=="cytosolic") )),
       aes(distance, tc, colour=hval_cat))+geom_line()+facet_wrap(~localization_cat+tc_region+hval_cat, scales="free_y", ncol=5)

ggplot(subset(plot, tc_region!="utr5"),aes(distance, tc, colour=hval_cat))+geom_line()+facet_wrap(~localization_cat+tc_region+hval_cat, scales="free_y", ncol=5)


