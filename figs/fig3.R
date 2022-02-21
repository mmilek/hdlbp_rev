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

sd<-subset(dat, tpm_cutoff>=10 & !is.na(localization_cat) & 
             !is.na(log2(norm_tc_num1)) & !is.na(log2(norm_tc_num2)), 
           select=c("gene_id", "transcript_id", "Symbol", "tc_region","tc_start", "tc_stop",
                    "norm_tc_num1","norm_tc_num2" ))
#write.table(sd, "source_data/s2a.txt", quote=F, sep="\t", row.names=F)

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


##redo 7mer enrichment in CDS/UTR3 mem/cyt
#run until tcCds


mas<-read.delim("data/hdlbp_master_table_with_classes_uniq_tsig.txt", header=T)

mem<-subset(mas, localization_cat=="membrane" & gene_biotype=="protein_coding" & tpm_cutoff>=10, select="transcript")[,1]
cyt<-subset(mas, localization_cat=="cytosolic" & gene_biotype=="protein_coding" & tpm_cutoff>=10, select="transcript")[,1]


seqCds_mem<-which(names(seqCds) %in% mem)
seqCds_mem<-seqCds[seqCds_mem]
seqCds_cyt<-which(names(seqCds) %in% cyt)
seqCds_cyt<-seqCds[seqCds_cyt]

seqUtr3_mem<-which(names(seqUtr3) %in% mem)
seqUtr3_mem<-seqUtr3[seqUtr3_mem]
seqUtr3_cyt<-which(names(seqUtr3) %in% cyt)
seqUtr3_cyt<-seqUtr3[seqUtr3_cyt]
kmer<-7
randsCds_mem <- (Biostrings::oligonucleotideFrequency(seqCds_mem, width=kmer, step = 1, as.prob = F, with.labels = TRUE, simplify.as="collapsed"))
freqCds_mem<-randsCds_mem/sum(randsCds_mem)

randsCds_cyt <- (Biostrings::oligonucleotideFrequency(seqCds_cyt, width=kmer, step = 1, as.prob = F, with.labels = TRUE, simplify.as="collapsed"))
freqCds_cyt<-randsCds_cyt/sum(randsCds_cyt)

randsUtr3_mem <- (Biostrings::oligonucleotideFrequency(seqUtr3_mem, width=kmer, step = 1, as.prob = F, with.labels = TRUE, simplify.as="collapsed"))
freqUtr3_mem<-randsUtr3_mem/sum(randsUtr3_mem)

randsUtr3_cyt <- (Biostrings::oligonucleotideFrequency(seqUtr3_cyt, width=kmer, step = 1, as.prob = F, with.labels = TRUE, simplify.as="collapsed"))
freqUtr3_cyt<-randsUtr3_cyt/sum(randsUtr3_cyt)

freqCdsUtr3_mem<-(randsCds_mem+randsUtr3_mem)/sum(randsCds_mem+randsUtr3_mem)
freqCdsUtr3_cyt<-(randsCds_cyt+randsUtr3_cyt)/sum(randsCds_cyt+randsUtr3_cyt)


freqs<-data.frame(kmer=names(randsCds_mem),
                  freqCds_mem=freqCds_mem,
                  freqCds_cyt=freqCds_cyt,
                  freqUtr3_mem=freqUtr3_mem,
                  freqUtr3_cyt=freqUtr3_cyt,
                  freqCdsUtr3_mem=freqCdsUtr3_mem,
                  freqCdsUtr3_cyt=freqCdsUtr3_cyt)
labs<-which(freqs$kmer %in% agr)
labs<-freqs[labs,]
labs$labs<-labs$kmer
freqs$labs<-NA
freqs<-rbind(freqs, labs)

mer_freq<-melt(freqs, measure.vars = colnames(freqs)[grepl("freq", colnames(freqs))], id.vars="kmer", value.name = "freq_whole", variable.name = "loc_reg")
mer_freq$id<-paste0(mer_freq$kmer, "_", mer_freq$loc_reg)

countmers$id<-ifelse(countmers$tc_region=="cds" & countmers$localization_cat=="cytosolic", paste0(countmers$seq, "_freqCds_cyt"),
                     ifelse(countmers$tc_region=="cds" & countmers$localization_cat=="membrane", paste0(countmers$seq, "_freqCds_mem"),
                            ifelse(countmers$tc_region=="utr3" & countmers$localization_cat=="cytosolic", paste0(countmers$seq, "_freqUtr3_cyt"),
                                   ifelse(countmers$tc_region=="utr3" & countmers$localization_cat=="membrane", paste0(countmers$seq, "_freqUtr3_mem"),NA))))

mer_freq<-merge(mer_freq, countmers, by="id")

labs<-which(mer_freq$kmer %in% agr)
labs<-mer_freq[labs,]
labs$labs<-labs$kmer
mer_freq$labs<-NA
mer_freq<-rbind(mer_freq, labs)

#figS3B
ggplot(mer_freq,aes(frac_tot, freq_whole))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_hline(yintercept=0.00015, lty=2, colour="grey")+theme(legend.position = "none")+
  facet_wrap(~localization_cat+tc_region)+coord_cartesian(ylim=c(0,0.001))

sd<-subset(mer_freq, !is.na(frac_tot) & !is.na(freq_whole), select=colnames(mer_freq)[2:length(colnames(mer_freq))])
sd$kmer<-gsub("T", "U", sd$kmer)

#write.table(sd, "source_data/figs3b.txt", quote=F, sep="\t", row.names=F)



# ggplot(mer_freq,aes(frac_tot, freq_whole))+geom_point()+geom_abline()+theme(legend.position = "none")+
#   facet_wrap(~localization_cat+tc_region)+coord_cartesian(ylim=c(0,0.001))

# mer_freq$top40<-ifelse(is.na(mer_freq$labs), "other", "top40")
# ggplot(mer_freq,aes( freq_whole, colour=top40))+geom_density()+
#   facet_wrap(~localization_cat+tc_region)+coord_cartesian(xlim=c(0,0.0005))
# 
# memCds<-subset(mer_freq, loc_reg=="freqCds_mem")
# cytCds<-subset(mer_freq, loc_reg=="freqCds_cyt")
# colnames(memCds)<-paste0(colnames(memCds), "_memCds")
# colnames(cytCds)<-paste0(colnames(cytCds), "_cytCds")
# 
# smer<-merge(memCds, cytCds, by.x="kmer_memCds", by.y="kmer_cytCds")
# ggplot(smer, aes(freq_whole_memCds-freq_whole_cytCds, frac_tot_memCds))+
#   geom_text(aes(label=kmer_memCds, colour=labs_memCds), size=2)+
#   geom_vline(xintercept=0,lty=2)+theme(legend.position = "none")
# 
# ggplot(smer, aes(freq_whole_memCds-freq_whole_cytCds, frac_tot_cytCds))+
#   geom_text(aes(label=kmer_memCds, colour=labs_memCds), size=2)+
#   geom_vline(xintercept=0,lty=2)+theme(legend.position = "none")
# 
# 
# memCds_bound<-subset(countmers, localization_cat=="membrane" & tc_region=="cds")
# cytUtr3_bound<-subset(countmers, localization_cat=="cytosolic" & tc_region=="utr3")
# colnames(memCds_bound)<-paste0(colnames(memCds_bound), "_memCds")
# colnames(cytUtr3_bound)<-paste0(colnames(cytUtr3_bound), "_cytUtr3")
# 
# smer<-merge(memCds_bound, cytUtr3_bound, by.x="seqs_memCds", by.y="seqs_cytUtr3")
# ggplot(smer, aes(frac_tot_memCds, frac_tot_cytUtr3))+
#   geom_text(aes(label=seqs_memCds), size=2)+
#   geom_abline(slope=1)
# ggplot(smer, aes(frac_tot_memCds, frac_tot_cytUtr3))+
#   geom_point()+
#   geom_abline(slope=1)
# 
# memCds_bound<-subset(countmers, localization_cat=="membrane" & tc_region=="cds")
# cytCds_bound<-subset(countmers, localization_cat=="cytosolic" & tc_region=="cds")
# colnames(memCds_bound)<-paste0(colnames(memCds_bound), "_memCds")
# colnames(cytCds_bound)<-paste0(colnames(cytCds_bound), "_cytCds")
# 
# smer<-merge(memCds_bound, cytCds_bound, by.x="seqs_memCds", by.y="seqs_cytCds")
# ggplot(smer, aes(frac_tot_memCds, frac_tot_cytCds))+
#   geom_text(aes(label=seqs_memCds), size=2)+
#   geom_abline(slope=1)

#7-mer per gene

ag1_cds<-dcast(subset(mer, tc_region=="cds", select=c("transcript_id","localization_cat", "seqs", "norm_tc_num1")),transcript_id+localization_cat~seqs,sum,value.var="norm_tc_num1")
ag2_cds<-dcast(subset(mer, tc_region=="cds", select=c("transcript_id","localization_cat", "seqs", "norm_tc_num2")), transcript_id+localization_cat~seqs,sum,value.var="norm_tc_num2")

ag1_utr3<-dcast(subset(mer, tc_region=="utr3", select=c("transcript_id","localization_cat", "seqs", "norm_tc_num1")), transcript_id+localization_cat~seqs,sum,value.var="norm_tc_num1")
ag2_utr3<-dcast(subset(mer, tc_region=="utr3", select=c("transcript_id","localization_cat", "seqs", "norm_tc_num2")), transcript_id+localization_cat~seqs,sum,value.var="norm_tc_num2")

ag1_cds_len<-which(len$transcript %in% ag1_cds$transcript_id)
ag2_cds_len<-which(len$transcript %in% ag2_cds$transcript_id)
ag1_utr3_len<-which(len$transcript %in% ag1_utr3$transcript_id)
ag2_utr3_len<-which(len$transcript %in% ag2_utr3$transcript_id)

ag1_cds_len<-len[ag1_cds_len, c("transcript", "l_cds")]
ag2_cds_len<-len[ag2_cds_len, c("transcript", "l_cds")]
ag1_utr3_len<-len[ag1_utr3_len, c("transcript", "l_utr3")]
ag2_utr3_len<-len[ag2_utr3_len, c("transcript", "l_utr3")]

ag1_cds<-cbind(ag1_cds[,1:2],ag1_cds[,3:ncol(ag1_cds)]/ag1_cds_len$l_cds*1e3)
ag2_cds<-cbind(ag2_cds[,1:2],ag2_cds[,3:ncol(ag2_cds)]/ag2_cds_len$l_cds*1e3)
ag1_utr3<-cbind(ag1_utr3[,1:2],ag1_utr3[,3:ncol(ag1_utr3)]/ag1_utr3_len$l_utr3*1e3)
ag2_utr3<-cbind(ag2_utr3[,1:2],ag2_utr3[,3:ncol(ag2_utr3)]/ag2_utr3_len$l_utr3*1e3)

ag1_cds_tpm<-which(mas$transcript %in% ag1_cds$transcript_id)
ag2_cds_tpm<-which(mas$transcript %in% ag2_cds$transcript_id)
ag1_utr3_tpm<-which(mas$transcript %in% ag1_utr3$transcript_id)
ag2_utr3_tpm<-which(mas$transcript %in% ag2_utr3$transcript_id)

ag1_cds_tpm<-mas[ag1_cds_tpm, c("transcript", "tpm_cutoff")]
ag1_cds_tpm<-ag1_cds_tpm[order(match(ag1_cds_tpm$transcript, ag1_cds$transcript_id)),]
ag2_cds_tpm<-mas[ag2_cds_tpm, c("transcript", "tpm_cutoff")]
ag2_cds_tpm<-ag2_cds_tpm[order(match(ag2_cds_tpm$transcript, ag2_cds$transcript_id)),]
ag1_utr3_tpm<-mas[ag1_utr3_tpm, c("transcript", "tpm_cutoff")]
ag1_utr3_tpm<-ag1_utr3_tpm[order(match(ag1_utr3_tpm$transcript, ag1_utr3$transcript_id)),]
ag2_utr3_tpm<-mas[ag2_utr3_tpm, c("transcript", "tpm_cutoff")]
ag2_utr3_tpm<-ag2_utr3_tpm[order(match(ag2_utr3_tpm$transcript, ag2_utr3$transcript_id)),]

ag1_cds<-cbind(ag1_cds[,1:2],ag1_cds[,3:ncol(ag1_cds)]/ag1_cds_tpm$tpm_cutoff*1e3)
ag2_cds<-cbind(ag2_cds[,1:2],ag2_cds[,3:ncol(ag2_cds)]/ag2_cds_tpm$tpm_cutoff*1e3)
ag1_utr3<-cbind(ag1_utr3[,1:2],ag1_utr3[,3:ncol(ag1_utr3)]/ag1_utr3_tpm$tpm_cutoff*1e3)
ag2_utr3<-cbind(ag2_utr3[,1:2],ag2_utr3[,3:ncol(ag2_utr3)]/ag2_utr3_tpm$tpm_cutoff*1e3)

ag1_cds$region<-"cds1"
ag2_cds$region<-"cds2"
ag1_utr3$region<-"utr3_1"
ag2_utr3$region<-"utr3_2"

p1_cds<-melt(ag1_cds, measure.vars = colnames(ag1_cds)[3:(ncol(ag1_cds)-1)], 
             id.vars = c("transcript_id", "localization_cat", "region"), value.name = "xl", variable.name = "kmer")
p2_cds<-melt(ag2_cds, measure.vars = colnames(ag2_cds)[3:(ncol(ag2_cds)-1)], 
             id.vars = c("transcript_id", "localization_cat", "region"), value.name = "xl", variable.name = "kmer")
p1_utr3<-melt(ag1_utr3, measure.vars = colnames(ag1_utr3)[3:(ncol(ag1_utr3)-1)], 
              id.vars = c("transcript_id", "localization_cat", "region"), value.name = "xl", variable.name = "kmer")
p2_utr3<-melt(ag2_utr3, measure.vars = colnames(ag2_utr3)[3:(ncol(ag2_utr3)-1)], 
              id.vars = c("transcript_id", "localization_cat", "region"), value.name = "xl", variable.name = "kmer")

plot<-rbind(p1_cds, p2_cds, p1_utr3, p2_utr3)
plot<-subset(plot, xl>0)
plot$log2xl<-ifelse(plot$xl<=0, NA, log2(plot$xl))
plot$count<-ifelse(plot$xl>0, 1, NA)
meds<-aggregate(log2xl~kmer+localization_cat+region, data=plot, median, na.rm=T)
meds<-meds[order(meds$log2xl, decreasing=T),]
meds$id<-paste0(meds$kmer,"_", meds$localization_cat,"_", meds$region)
freq<-aggregate(count~kmer+localization_cat+region, data=plot, sum)
freq<-freq[order(freq$count, decreasing=T),]
freq$id<-paste0(freq$kmer,"_", freq$localization_cat,"_", freq$region)

ssr<-merge(meds, freq, by="id")
ggplot(ssr, aes(log2xl, count, fill=region.x ))+geom_text(aes(label=kmer.x), size=2)+facet_wrap(~localization_cat.x)
ssr<-subset(ssr, count>20)
ssr<-ssr[order(ssr$count, decreasing=T),]

top40<-subset(ssr, localization_cat.x=="membrane" & region.x=="cds1")$kmer.x[1:10]
plot_sub<-which(plot$kmer %in% top40)
plot_sub<-plot[plot_sub,]
plot_sub$kmer<-factor(plot_sub$kmer, levels=rev(top40))
plot_sub$kmer<-factor(plot_sub$kmer, levels=(top40))
plot_sub$localization_cat<-factor(plot_sub$localization_cat, levels=c("membrane", "cytosolic"))

#figS3a
ggplot(plot_sub, aes(kmer, log2xl, fill=region))+geom_hline(yintercept = 5, lty=2, colour="grey")+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~localization_cat, ncol=1)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_fill_manual(values=c("dodgerblue4", "dodgerblue3", "orange2", "orange3"))

sd<-subset(plot_sub, !is.na(log2xl) & !is.na(kmer) )
sd$kmer<-gsub("T", "U", sd$kmer)

#write.table(sd, "source_data/figs3a.txt", quote=F, sep="\t", row.names=F)



summary(subset(plot_sub, localization_cat =="membrane" &region=="cds1" & !is.na(log2xl))$kmer)
summary(subset(plot_sub, localization_cat =="membrane" &region=="cds2" & !is.na(log2xl))$kmer)
summary(subset(plot_sub, localization_cat =="membrane" &region=="utr3_1" & !is.na(log2xl))$kmer)
summary(subset(plot_sub, localization_cat =="membrane" &region=="utr3_2" & !is.na(log2xl))$kmer)
summary(subset(plot_sub, localization_cat =="cytosolic" &region=="cds1" & !is.na(log2xl))$kmer)
summary(subset(plot_sub, localization_cat =="cytosolic" &region=="cds2" & !is.na(log2xl))$kmer)
summary(subset(plot_sub, localization_cat =="cytosolic" &region=="utr3_1" & !is.na(log2xl))$kmer)
summary(subset(plot_sub, localization_cat =="cytosolic" &region=="utr3_2" & !is.na(log2xl))$kmer)

ssr<-merge(meds, freq, by="id")
ssr_sub<-which(ssr$kmer.x %in% top40)
ssr_sub<-ssr[ssr_sub,]
ssr_sub$kmer.x<-factor(ssr_sub$kmer.x, levels=(top40))
ssr_sub$localization_cat.x<-factor(ssr_sub$localization_cat.x, levels=c("membrane", "cytosolic"))

#figS3a
ggplot(ssr_sub, aes(kmer.x, count, fill=region.x))+geom_bar(stat = "identity", position = "dodge")+
  facet_wrap(~localization_cat.x, ncol=1)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_fill_manual(values=c("dodgerblue4", "dodgerblue3", "orange2", "orange3"))

sd<-subset(ssr_sub, !is.na(kmer.x) & !is.na(count) , select=c("kmer.x", "localization_cat.x", 
                                                              "region.x", "count"))
sd$kmer.x<-gsub("T", "U", sd$kmer)

#write.table(sd, "source_data/figs3a_right.txt", quote=F, sep="\t", row.names=F)

##4mer crosslinked
mer<-dat

mer$count<-ifelse(is.na(mer$localization_cat), NA,
                  ifelse(mer$norm_tc_num1>2 & mer$norm_tc_num2>2 & mer$TPM_transcript>=10 , 1, NA))
mer<-subset(mer, !is.na(count))


mer$tc_seq_pos1_start<-ifelse(mer$tc_start-2<1, NA, mer$tc_start-2)
mer$tc_seq_pos1_stop<-ifelse(mer$tc_start+1>mer$l_tr, NA, mer$tc_start+1)
mer$tc_seq_pos2_start<-ifelse(mer$tc_start-1<1, NA, mer$tc_start-1)
mer$tc_seq_pos2_stop<-ifelse(mer$tc_start+2>mer$l_tr, NA, mer$tc_start+2)
mer$tc_seq_pos3_start<-ifelse(mer$tc_start<1, NA, mer$tc_start)
mer$tc_seq_pos3_stop<-ifelse(mer$tc_start+3>mer$l_tr, NA, mer$tc_start+3)
mer$tc_seq_pos4_start<-ifelse(mer$tc_start+1<1, NA, mer$tc_start+1)
mer$tc_seq_pos4_stop<-ifelse(mer$tc_start+4>mer$l_tr, NA, mer$tc_start+4)



mer$tc_seq_pos1<-toupper(as.character(Biostrings::subseq(sub[as.character(mer$transcript_id)],
                                                         start=mer$tc_seq_pos1_start,end=mer$tc_seq_pos1_stop)))
mer$tc_seq_pos2<-toupper(as.character(Biostrings::subseq(sub[as.character(mer$transcript_id)],
                                                         start=mer$tc_seq_pos2_start,end=mer$tc_seq_pos2_stop)))
mer$tc_seq_pos3<-toupper(as.character(Biostrings::subseq(sub[as.character(mer$transcript_id)],
                                                         start=mer$tc_seq_pos3_start,end=mer$tc_seq_pos3_stop)))
mer$tc_seq_pos4<-toupper(as.character(Biostrings::subseq(sub[as.character(mer$transcript_id)],
                                                         start=mer$tc_seq_pos4_start,end=mer$tc_seq_pos4_stop)))

mer<-melt(mer, measure.vars = colnames(mer)[42:45], id.vars = colnames(mer)[1:33], value.name = "seqs", variable.name = "tc_seq_pos")


# countmers_gene<-dcast(subset(mer, !is.na(count), select=c("transcript_id", "seq", "count")),transcript_id~seq,sum,value.var="count")

countmers<-aggregate(count~seqs+tc_region+localization_cat, sum, data=mer)

countmers$frac_tot<-countmers$count/sum(countmers$count)
agr<-aggregate(frac_tot~seqs, sum, data=countmers)
agr<-agr[order(agr$frac_tot, decreasing = T),]
agr<-agr[1:10, "seqs"]
rows<-which(countmers$seqs %in% agr)
countmers_plot<-countmers[rows,]
countmers_plot$seqs<-factor(countmers_plot$seqs, levels=rev(agr))
countmers_plot$localization_cat<-factor(countmers_plot$localization_cat, levels=c("membrane", "cytosolic"))
ggplot(subset(countmers_plot, tc_region!="utr5"), aes(seqs, frac_tot, fill=tc_region))+geom_bar(stat="identity", position="dodge",)+coord_flip()+facet_wrap(~localization_cat)+
  theme(text = element_text(size=8))+scale_fill_manual(values = c("dodgerblue4", "orange3"))

library(ggseqlogo)
ggplot(subset(countmers_plot, localization_cat=="cytosolic"))+geom_logo(gsub("T","U",as.character(countmers_plot$seqs)), method="bits", seq_type = "rna")+theme_logo()

#write.table(countmers, "./kmers_mem_cyt/bound-4mers-freqs-01032020.txt", quote=F, sep="\t", row.names = F)

##redo4mer enrichment in CDS/UTR3 mem/cyt
#run until tcCds


setwd("E:/Google Drive/hdlbp/")
# mas<-read.delim("hdlbp_master_table_with_classes.txt", header=T)
mas<-read.delim("hdlbp_master_table_with_classes_uniq_tsig.txt", header=T)
inf<-subset(mas, select=c("gene_id","Symbol",
                          "tpm_cutoff","tc_CDS_norm","localization_cat",
                          "mean_te_293","loc_tar_CDS","tc_CDS_norm_cat",
                          "tc_transcript_norm",
                          "log2FoldChange.mem.cyt.293",
                          "log2FoldChange.mem.cyt.KO.293",
                          "log2FoldChange.ribo.rna.KO.WT"))


mem<-subset(mas, localization_cat=="membrane" & gene_biotype=="protein_coding" & tpm_cutoff>=10, select="transcript")[,1]
cyt<-subset(mas, localization_cat=="cytosolic" & gene_biotype=="protein_coding" & tpm_cutoff>=10, select="transcript")[,1]


seqCds_mem<-which(names(seqCds) %in% mem)
seqCds_mem<-seqCds[seqCds_mem]
seqCds_cyt<-which(names(seqCds) %in% cyt)
seqCds_cyt<-seqCds[seqCds_cyt]

seqUtr3_mem<-which(names(seqUtr3) %in% mem)
seqUtr3_mem<-seqUtr3[seqUtr3_mem]
seqUtr3_cyt<-which(names(seqUtr3) %in% cyt)
seqUtr3_cyt<-seqUtr3[seqUtr3_cyt]
kmer<-4
randsCds_mem <- (Biostrings::oligonucleotideFrequency(seqCds_mem, width=kmer, step = 1, as.prob = F, with.labels = TRUE))
freqCds_mem<-colSums(randsCds_mem)/sum(colSums(randsCds_mem))

randsCds_cyt <- (Biostrings::oligonucleotideFrequency(seqCds_cyt, width=kmer, step = 1, as.prob = F, with.labels = TRUE))
freqCds_cyt<-colSums(randsCds_cyt)/sum(colSums(randsCds_cyt))

randsUtr3_mem <- (Biostrings::oligonucleotideFrequency(seqUtr3_mem, width=kmer, step = 1, as.prob = F, with.labels = TRUE))
freqUtr3_mem<-colSums(randsUtr3_mem)/sum(colSums(randsUtr3_mem))

randsUtr3_cyt <- (Biostrings::oligonucleotideFrequency(seqUtr3_cyt, width=kmer, step = 1, as.prob = F, with.labels = TRUE))
freqUtr3_cyt<-colSums(randsUtr3_cyt)/sum(colSums(randsUtr3_cyt))

freqCdsUtr3_mem<-(colSums(randsCds_mem)+colSums(randsUtr3_mem))/sum(colSums(randsCds_mem)+colSums(randsUtr3_mem))
freqCdsUtr3_cyt<-(colSums(randsCds_cyt)+colSums(randsUtr3_cyt))/sum(colSums(randsCds_cyt)+colSums(randsUtr3_cyt))


freqs<-data.frame(kmer=names(colSums(randsCds_mem)),
                  freqCds_mem=freqCds_mem,
                  freqCds_cyt=freqCds_cyt,
                  freqUtr3_mem=freqUtr3_mem,
                  freqUtr3_cyt=freqUtr3_cyt,
                  freqCdsUtr3_mem=freqCdsUtr3_mem,
                  freqCdsUtr3_cyt=freqCdsUtr3_cyt)
labs<-which(freqs$kmer %in% agr)
labs<-freqs[labs,]
labs$labs<-labs$kmer
freqs$labs<-NA
freqs<-rbind(freqs, labs)
ggplot(freqs, aes(freqCds_mem, freqCds_cyt))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()
ggplot(freqs, aes(freqCds_mem, freqUtr3_mem))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()
ggplot(freqs, aes(freqCds_mem, freqUtr3_mem, colour=labs))+geom_point()+geom_abline()
ggplot(freqs, aes(freqCds_mem, freqUtr3_cyt))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()
ggplot(freqs, aes(freqCdsUtr3_mem, freqCdsUtr3_cyt))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()
ggplot(freqs, aes(freqCdsUtr3_mem, freqCdsUtr3_cyt, colour=labs))+geom_point()+geom_abline()

ggplot(freqs, aes(freqCds_cyt, freqUtr3_cyt))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()
ggplot(freqs, aes(freqCds_cyt, freqUtr3_cyt, colour=labs))+geom_point()+geom_abline()
ggplot(freqs, aes(freqUtr3_mem, freqUtr3_cyt, colour=labs))+geom_point()+geom_abline()
ggplot(freqs, aes(freqCds_cyt, freqUtr3_cyt))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()

#write.table(freqs, "./kmers_mem_cyt/whole-4mers-freqs-01032020.txt", quote=F, sep="\t", row.names=F)

mer_freq<-melt(freqs, measure.vars = colnames(freqs)[grepl("freq", colnames(freqs))], id.vars="kmer", value.name = "freq_whole", variable.name = "loc_reg")
mer_freq$id<-paste0(mer_freq$kmer, "_", mer_freq$loc_reg)

countmers$id<-ifelse(countmers$tc_region=="cds" & countmers$localization_cat=="cytosolic", paste0(countmers$seq, "_freqCds_cyt"),
                     ifelse(countmers$tc_region=="cds" & countmers$localization_cat=="membrane", paste0(countmers$seq, "_freqCds_mem"),
                            ifelse(countmers$tc_region=="utr3" & countmers$localization_cat=="cytosolic", paste0(countmers$seq, "_freqUtr3_cyt"),
                                   ifelse(countmers$tc_region=="utr3" & countmers$localization_cat=="membrane", paste0(countmers$seq, "_freqUtr3_mem"),NA))))

mer_freq<-merge(mer_freq, countmers, by="id")

labs<-which(mer_freq$kmer %in% agr)
labs<-mer_freq[labs,]
labs$labs<-labs$kmer
mer_freq$labs<-NA
mer_freq<-rbind(mer_freq, labs)
ggplot(mer_freq,aes(frac_tot, freq_whole))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()+theme(legend.position = "none")+
  facet_wrap(~localization_cat+tc_region)

ggplot(mer_freq,aes(frac_tot, freq_whole))+geom_point()+geom_abline()+theme(legend.position = "none")+
  facet_wrap(~localization_cat+tc_region)

##select also unbound kmers with similar frequency as the bound ones
ags<-(unique(subset(mer_freq, freq_whole<0.005 & freq_whole>0.003 & frac_tot<0.0025, select="kmer")[,1]))
ags<-ags[grepl("AG",ags)]
ags<-ags[!grepl("CT",ags)]
ags<-ags[!grepl("TC",ags)]

mer_freq$top40<-ifelse(is.na(mer_freq$labs), "other", "top40")
ggplot(mer_freq,aes( freq_whole, colour=top40))+geom_density()+
  facet_wrap(~localization_cat+tc_region)+coord_cartesian(xlim=c(0,0.02))

memCds<-subset(mer_freq, loc_reg=="freqCds_mem")
cytCds<-subset(mer_freq, loc_reg=="freqCds_cyt")
colnames(memCds)<-paste0(colnames(memCds), "_memCds")
colnames(cytCds)<-paste0(colnames(cytCds), "_cytCds")

smer<-merge(memCds, cytCds, by.x="kmer_memCds", by.y="kmer_cytCds")
ggplot(smer, aes(freq_whole_memCds-freq_whole_cytCds, frac_tot_memCds))+
  geom_text(aes(label=kmer_memCds, colour=labs_memCds), size=2)+
  geom_vline(xintercept=0,lty=2)+theme(legend.position = "none")

ggplot(smer, aes(freq_whole_memCds-freq_whole_cytCds, frac_tot_cytCds))+
  geom_text(aes(label=kmer_memCds, colour=labs_memCds), size=2)+
  geom_vline(xintercept=0,lty=2)+theme(legend.position = "none")


memCds_bound<-subset(countmers, localization_cat=="membrane" & tc_region=="cds")
cytUtr3_bound<-subset(countmers, localization_cat=="cytosolic" & tc_region=="utr3")
colnames(memCds_bound)<-paste0(colnames(memCds_bound), "_memCds")
colnames(cytUtr3_bound)<-paste0(colnames(cytUtr3_bound), "_cytUtr3")

smer<-merge(memCds_bound, cytUtr3_bound, by.x="seqs_memCds", by.y="seqs_cytUtr3")
ggplot(smer, aes(frac_tot_memCds, frac_tot_cytUtr3))+
  geom_text(aes(label=seqs_memCds), size=2)+
  geom_abline(slope=1)
ggplot(smer, aes(frac_tot_memCds, frac_tot_cytUtr3))+
  geom_point()+
  geom_abline(slope=1)

memCds_bound<-subset(countmers, localization_cat=="membrane" & tc_region=="cds")
cytCds_bound<-subset(countmers, localization_cat=="cytosolic" & tc_region=="cds")
colnames(memCds_bound)<-paste0(colnames(memCds_bound), "_memCds")
colnames(cytCds_bound)<-paste0(colnames(cytCds_bound), "_cytCds")

smer<-merge(memCds_bound, cytCds_bound, by.x="seqs_memCds", by.y="seqs_cytCds")
ggplot(smer, aes(frac_tot_memCds, frac_tot_cytCds))+
  geom_text(aes(label=seqs_memCds), size=2)+
  geom_abline(slope=1)

#4-mer per gene

ag1_cds<-dcast(subset(mer, tc_region=="cds", select=c("transcript_id","localization_cat", "seqs", "norm_tc_num1")),transcript_id+localization_cat~seqs,sum,value.var="norm_tc_num1")
ag2_cds<-dcast(subset(mer, tc_region=="cds", select=c("transcript_id","localization_cat", "seqs", "norm_tc_num2")), transcript_id+localization_cat~seqs,sum,value.var="norm_tc_num2")

ag1_utr3<-dcast(subset(mer, tc_region=="utr3", select=c("transcript_id","localization_cat", "seqs", "norm_tc_num1")), transcript_id+localization_cat~seqs,sum,value.var="norm_tc_num1")
ag2_utr3<-dcast(subset(mer, tc_region=="utr3", select=c("transcript_id","localization_cat", "seqs", "norm_tc_num2")), transcript_id+localization_cat~seqs,sum,value.var="norm_tc_num2")

ag1_cds_len<-which(len$transcript %in% ag1_cds$transcript_id)
ag2_cds_len<-which(len$transcript %in% ag2_cds$transcript_id)
ag1_utr3_len<-which(len$transcript %in% ag1_utr3$transcript_id)
ag2_utr3_len<-which(len$transcript %in% ag2_utr3$transcript_id)

ag1_cds_len<-len[ag1_cds_len, c("transcript", "l_cds")]
ag2_cds_len<-len[ag2_cds_len, c("transcript", "l_cds")]
ag1_utr3_len<-len[ag1_utr3_len, c("transcript", "l_utr3")]
ag2_utr3_len<-len[ag2_utr3_len, c("transcript", "l_utr3")]

ag1_cds<-cbind(ag1_cds[,1:2],ag1_cds[,3:ncol(ag1_cds)]/ag1_cds_len$l_cds*1e3)
ag2_cds<-cbind(ag2_cds[,1:2],ag2_cds[,3:ncol(ag2_cds)]/ag2_cds_len$l_cds*1e3)
ag1_utr3<-cbind(ag1_utr3[,1:2],ag1_utr3[,3:ncol(ag1_utr3)]/ag1_utr3_len$l_utr3*1e3)
ag2_utr3<-cbind(ag2_utr3[,1:2],ag2_utr3[,3:ncol(ag2_utr3)]/ag2_utr3_len$l_utr3*1e3)

ag1_cds_tpm<-which(mas$transcript %in% ag1_cds$transcript_id)
ag2_cds_tpm<-which(mas$transcript %in% ag2_cds$transcript_id)
ag1_utr3_tpm<-which(mas$transcript %in% ag1_utr3$transcript_id)
ag2_utr3_tpm<-which(mas$transcript %in% ag2_utr3$transcript_id)

ag1_cds_tpm<-mas[ag1_cds_tpm, c("transcript", "tpm_cutoff")]
ag1_cds_tpm<-ag1_cds_tpm[order(match(ag1_cds_tpm$transcript, ag1_cds$transcript_id)),]
ag2_cds_tpm<-mas[ag2_cds_tpm, c("transcript", "tpm_cutoff")]
ag2_cds_tpm<-ag2_cds_tpm[order(match(ag2_cds_tpm$transcript, ag2_cds$transcript_id)),]
ag1_utr3_tpm<-mas[ag1_utr3_tpm, c("transcript", "tpm_cutoff")]
ag1_utr3_tpm<-ag1_utr3_tpm[order(match(ag1_utr3_tpm$transcript, ag1_utr3$transcript_id)),]
ag2_utr3_tpm<-mas[ag2_utr3_tpm, c("transcript", "tpm_cutoff")]
ag2_utr3_tpm<-ag2_utr3_tpm[order(match(ag2_utr3_tpm$transcript, ag2_utr3$transcript_id)),]

ag1_cds<-cbind(ag1_cds[,1:2],ag1_cds[,3:ncol(ag1_cds)]/ag1_cds_tpm$tpm_cutoff)
ag2_cds<-cbind(ag2_cds[,1:2],ag2_cds[,3:ncol(ag2_cds)]/ag2_cds_tpm$tpm_cutoff)
ag1_utr3<-cbind(ag1_utr3[,1:2],ag1_utr3[,3:ncol(ag1_utr3)]/ag1_utr3_tpm$tpm_cutoff)
ag2_utr3<-cbind(ag2_utr3[,1:2],ag2_utr3[,3:ncol(ag2_utr3)]/ag2_utr3_tpm$tpm_cutoff)

ag1_cds$region<-"cds1"
ag2_cds$region<-"cds2"
ag1_utr3$region<-"utr3_1"
ag2_utr3$region<-"utr3_2"

p1_cds<-melt(ag1_cds, measure.vars = colnames(ag1_cds)[3:(ncol(ag1_cds)-1)], 
             id.vars = c("transcript_id", "localization_cat", "region"), value.name = "xl", variable.name = "kmer")
p2_cds<-melt(ag2_cds, measure.vars = colnames(ag2_cds)[3:(ncol(ag2_cds)-1)], 
             id.vars = c("transcript_id", "localization_cat", "region"), value.name = "xl", variable.name = "kmer")
p1_utr3<-melt(ag1_utr3, measure.vars = colnames(ag1_utr3)[3:(ncol(ag1_utr3)-1)], 
              id.vars = c("transcript_id", "localization_cat", "region"), value.name = "xl", variable.name = "kmer")
p2_utr3<-melt(ag2_utr3, measure.vars = colnames(ag2_utr3)[3:(ncol(ag2_utr3)-1)], 
              id.vars = c("transcript_id", "localization_cat", "region"), value.name = "xl", variable.name = "kmer")

plot<-rbind(p1_cds, p2_cds, p1_utr3, p2_utr3)
plot<-subset(plot, xl>0)
plot$log2xl<-ifelse(plot$xl<=0, NA, log2(plot$xl))
plot$count<-ifelse(plot$xl>0, 1, NA)
meds<-aggregate(log2xl~kmer+localization_cat+region, data=plot, median, na.rm=T)
meds<-meds[order(meds$log2xl, decreasing=T),]
meds$id<-paste0(meds$kmer,"_", meds$localization_cat,"_", meds$region)
freq<-aggregate(count~kmer+localization_cat+region, data=plot, sum)
freq<-freq[order(freq$count, decreasing=T),]
freq$id<-paste0(freq$kmer,"_", freq$localization_cat,"_", freq$region)

ssr<-merge(meds, freq, by="id")
ggplot(ssr, aes(log2xl, count, fill=region.x ))+geom_text(aes(label=kmer.x), size=2)+facet_wrap(~localization_cat.x)
ssr<-subset(ssr, count>20)
ssr<-ssr[order(ssr$log2xl, decreasing=T),]

top40<-subset(ssr, localization_cat.x=="membrane" & region.x=="cds1")$kmer.x[1:10]
plot_sub<-which(plot$kmer %in% top40)
plot_sub<-plot[plot_sub,]
plot_sub$kmer<-factor(plot_sub$kmer, levels=rev(top40))
ggplot(plot_sub, aes(kmer, log2xl, fill=region))+geom_boxplot()+facet_wrap(~localization_cat)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
ssr<-merge(meds, freq, by="id")
ssr_sub<-which(ssr$kmer.x %in% top40)
ssr_sub<-ssr[ssr_sub,]
ssr_sub$kmer.x<-factor(ssr_sub$kmer.x, levels=rev(top40))
ggplot(ssr_sub, aes(kmer.x, count, fill=region.x))+geom_bar(stat = "identity", position = "dodge")+facet_wrap(~localization_cat.x)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))


##9mer crosslinked
mer<-dat

mer$count<-ifelse(is.na(mer$localization_cat), NA,
                  ifelse(mer$norm_tc_num1>2 & mer$norm_tc_num2>2 & mer$TPM_transcript>=10 , 1, NA))
mer<-subset(mer, !is.na(count))
mer$tc_seq_pos1_start<-ifelse(mer$tc_start-7<1, NA, mer$tc_start-7)
mer$tc_seq_pos1_stop<-ifelse(mer$tc_start+1>mer$l_tr, NA, mer$tc_start+1)
mer$tc_seq_pos2_start<-ifelse(mer$tc_start-6<1, NA, mer$tc_start-6)
mer$tc_seq_pos2_stop<-ifelse(mer$tc_start+2>mer$l_tr, NA, mer$tc_start+2)
mer$tc_seq_pos3_start<-ifelse(mer$tc_start-5<1, NA, mer$tc_start-5)
mer$tc_seq_pos3_stop<-ifelse(mer$tc_start+3>mer$l_tr, NA, mer$tc_start+3)
mer$tc_seq_pos4_start<-ifelse(mer$tc_start-4<1, NA, mer$tc_start-4)
mer$tc_seq_pos4_stop<-ifelse(mer$tc_start+4>mer$l_tr, NA, mer$tc_start+4)
mer$tc_seq_pos5_start<-ifelse(mer$tc_start-3<1, NA, mer$tc_start-3)
mer$tc_seq_pos5_stop<-ifelse(mer$tc_start+5>mer$l_tr, NA, mer$tc_start+5)
mer$tc_seq_pos6_start<-ifelse(mer$tc_start-2<1, NA, mer$tc_start-2)
mer$tc_seq_pos6_stop<-ifelse(mer$tc_start+6>mer$l_tr, NA, mer$tc_start+6)
mer$tc_seq_pos7_start<-ifelse(mer$tc_start-1<1, NA, mer$tc_start-1)
mer$tc_seq_pos7_stop<-ifelse(mer$tc_start+7>mer$l_tr, NA, mer$tc_start+7)
mer$tc_seq_pos8_start<-ifelse(mer$tc_start<1, NA, mer$tc_start)
mer$tc_seq_pos8_stop<-ifelse(mer$tc_start+8>mer$l_tr, NA, mer$tc_start+8)
mer$tc_seq_pos9_start<-ifelse(mer$tc_start+1<1, NA, mer$tc_start+1)
mer$tc_seq_pos9_stop<-ifelse(mer$tc_start+9>mer$l_tr, NA, mer$tc_start+9)




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
mer$tc_seq_pos8<-toupper(as.character(Biostrings::subseq(sub[as.character(mer$transcript_id)],
                                                         start=mer$tc_seq_pos8_start,end=mer$tc_seq_pos8_stop)))
mer$tc_seq_pos9<-toupper(as.character(Biostrings::subseq(sub[as.character(mer$transcript_id)],
                                                         start=mer$tc_seq_pos9_start,end=mer$tc_seq_pos9_stop)))

mer<-melt(mer, measure.vars = colnames(mer)[52:60], id.vars = colnames(mer)[1:33], value.name = "seqs", variable.name = "tc_seq_pos")


# countmers_gene<-dcast(subset(mer, !is.na(count), select=c("transcript_id", "seq", "count")),transcript_id~seq,sum,value.var="count")

countmers<-aggregate(count~seqs+tc_region+localization_cat, sum, data=mer)

countmers$frac_tot<-countmers$count/sum(countmers$count)
agr<-aggregate(frac_tot~seqs, sum, data=countmers)
agr<-agr[order(agr$frac_tot, decreasing = T),]
agr<-agr[1:40, "seqs"]
rows<-which(countmers$seqs %in% agr)
countmers_plot<-countmers[rows,]
countmers_plot$seqs<-factor(countmers_plot$seqs, levels=rev(agr))
ggplot(subset(countmers_plot, tc_region!="utr5"), aes(seqs, frac_tot, fill=tc_region))+geom_bar(stat="identity", position="dodge",)+coord_flip()+facet_wrap(~localization_cat)+
  theme(text = element_text(size=8))

library(ggseqlogo)
ggplot(subset(countmers_plot, localization_cat=="cytosolic"))+geom_logo(gsub("T","U",as.character(countmers_plot$seqs)), method="bits", seq_type = "rna")+theme_logo()

#write.table(countmers, "./kmers_mem_cyt/bound-9mers-freqs-01032020.txt", quote=F, sep="\t", row.names = F)

##redo 9mer enrichment in CDS/UTR3 mem/cyt
#run until tcCds


setwd("E:/Google Drive/hdlbp/")
# mas<-read.delim("hdlbp_master_table_with_classes.txt", header=T)
mas<-read.delim("data/hdlbp_master_table_with_classes_uniq_tsig.txt", header=T)
inf<-subset(mas, select=c("gene_id","Symbol",
                          "tpm_cutoff","tc_CDS_norm","localization_cat",
                          "mean_te_293","loc_tar_CDS","tc_CDS_norm_cat",
                          "tc_transcript_norm",
                          "log2FoldChange.mem.cyt.293",
                          "log2FoldChange.mem.cyt.KO.293",
                          "log2FoldChange.ribo.rna.KO.WT"))


mem<-subset(mas, localization_cat=="membrane" & gene_biotype=="protein_coding" & tpm_cutoff>=10, select="transcript")[,1]
cyt<-subset(mas, localization_cat=="cytosolic" & gene_biotype=="protein_coding" & tpm_cutoff>=10, select="transcript")[,1]


seqCds_mem<-which(names(seqCds) %in% mem)
seqCds_mem<-seqCds[seqCds_mem]
seqCds_cyt<-which(names(seqCds) %in% cyt)
seqCds_cyt<-seqCds[seqCds_cyt]

seqUtr3_mem<-which(names(seqUtr3) %in% mem)
seqUtr3_mem<-seqUtr3[seqUtr3_mem]
seqUtr3_cyt<-which(names(seqUtr3) %in% cyt)
seqUtr3_cyt<-seqUtr3[seqUtr3_cyt]
kmer<-9
randsCds_mem <- (Biostrings::oligonucleotideFrequency(seqCds_mem, width=kmer, step = 1, as.prob = F, with.labels = TRUE, simplify.as="collapsed"))
freqCds_mem<-randsCds_mem/sum(randsCds_mem)

randsCds_cyt <- (Biostrings::oligonucleotideFrequency(seqCds_cyt, width=kmer, step = 1, as.prob = F, with.labels = TRUE, simplify.as="collapsed"))
freqCds_cyt<-randsCds_cyt/sum(randsCds_cyt)

randsUtr3_mem <- (Biostrings::oligonucleotideFrequency(seqUtr3_mem, width=kmer, step = 1, as.prob = F, with.labels = TRUE, simplify.as="collapsed"))
freqUtr3_mem<-randsUtr3_mem/sum(randsUtr3_mem)

randsUtr3_cyt <- (Biostrings::oligonucleotideFrequency(seqUtr3_cyt, width=kmer, step = 1, as.prob = F, with.labels = TRUE, simplify.as="collapsed"))
freqUtr3_cyt<-randsUtr3_cyt/sum(randsUtr3_cyt)

freqCdsUtr3_mem<-(randsCds_mem+randsUtr3_mem)/sum(randsCds_mem+randsUtr3_mem)
freqCdsUtr3_cyt<-(randsCds_cyt+randsUtr3_cyt)/sum(randsCds_cyt+randsUtr3_cyt)


freqs<-data.frame(kmer=names(randsCds_mem),
                  freqCds_mem=freqCds_mem,
                  freqCds_cyt=freqCds_cyt,
                  freqUtr3_mem=freqUtr3_mem,
                  freqUtr3_cyt=freqUtr3_cyt,
                  freqCdsUtr3_mem=freqCdsUtr3_mem,
                  freqCdsUtr3_cyt=freqCdsUtr3_cyt)
labs<-which(freqs$kmer %in% agr)
labs<-freqs[labs,]
labs$labs<-labs$kmer
freqs$labs<-NA
freqs<-rbind(freqs, labs)
ggplot(freqs, aes(freqCds_mem, freqCds_cyt))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()
ggplot(freqs, aes(freqCds_mem, freqUtr3_mem))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()
ggplot(freqs, aes(freqCds_mem, freqUtr3_mem, colour=labs))+geom_point()+geom_abline()
ggplot(freqs, aes(freqCds_mem, freqUtr3_cyt))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()
ggplot(freqs, aes(freqCdsUtr3_mem, freqCdsUtr3_cyt))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()+coord_cartesian(xlim=c(0,0.0002), ylim=c(0,0.0002))
ggplot(freqs, aes(freqCdsUtr3_mem, freqCdsUtr3_cyt, colour=labs))+geom_point()+geom_abline()

ggplot(freqs, aes(freqCds_cyt, freqUtr3_cyt))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()
ggplot(freqs, aes(freqCds_cyt, freqUtr3_cyt, colour=labs))+geom_point()+geom_abline()
ggplot(freqs, aes(freqUtr3_mem, freqUtr3_cyt, colour=labs))+geom_point()+geom_abline()
ggplot(freqs, aes(freqCds_cyt, freqUtr3_cyt))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()

#write.table(freqs, "./kmers_mem_cyt/whole-9mers-freqs-01032020.txt", quote=F, sep="\t", row.names=F)
mer_freq<-melt(freqs, measure.vars = colnames(freqs)[grepl("freq", colnames(freqs))], id.vars="kmer", value.name = "freq_whole", variable.name = "loc_reg")
mer_freq$id<-paste0(mer_freq$kmer, "_", mer_freq$loc_reg)

countmers$id<-ifelse(countmers$tc_region=="cds" & countmers$localization_cat=="cytosolic", paste0(countmers$seq, "_freqCds_cyt"),
                     ifelse(countmers$tc_region=="cds" & countmers$localization_cat=="membrane", paste0(countmers$seq, "_freqCds_mem"),
                            ifelse(countmers$tc_region=="utr3" & countmers$localization_cat=="cytosolic", paste0(countmers$seq, "_freqUtr3_cyt"),
                                   ifelse(countmers$tc_region=="utr3" & countmers$localization_cat=="membrane", paste0(countmers$seq, "_freqUtr3_mem"),NA))))

mer_freq<-merge(mer_freq, countmers, by="id")

labs<-which(mer_freq$kmer %in% agr)
labs<-mer_freq[labs,]
labs$labs<-labs$kmer
mer_freq$labs<-NA
mer_freq<-rbind(mer_freq, labs)
mer_freq$kmer<-gsub("T","U", mer_freq$kmer)
mer_freq$labs<-gsub("T","U", mer_freq$labs)
mer_freq$localization_cat<-factor(mer_freq$localization_cat, levels =c("membrane", "cytosolic"))
mer_freq$tc_region<-factor(mer_freq$tc_region, levels =c("cds", "utr3"))

ggplot(mer_freq,aes(frac_tot, freq_whole))+geom_hline(lty=2,yintercept=0.00002)+
  geom_text(aes(label=kmer, colour=labs), size=2)+
  theme(legend.position = "none")+
  facet_wrap(~localization_cat+tc_region)+coord_cartesian(ylim=c(0,0.0001))

ggplot(mer_freq,aes(frac_tot, freq_whole, colour=labs))+geom_hline(lty=2,yintercept=0.00002)+
  geom_point(shape=1, size=1)+geom_text(aes(label=labs), size=2)+
  theme(legend.position = "none")+
  facet_wrap(~localization_cat+tc_region)+coord_cartesian(ylim=c(0,0.0001))

mer_freq$top40<-ifelse(is.na(mer_freq$labs), "other", "top40")
ggplot(mer_freq,aes( freq_whole, colour=top40))+geom_density()+
  facet_wrap(~localization_cat+tc_region)+coord_cartesian(xlim=c(0,0.00005))
ggplot(mer_freq,aes( freq_whole, colour=loc_reg))+stat_ecdf()+
  coord_cartesian(xlim=c(0,0.00005))+facet_wrap(~top40)
ggplot(mer_freq,aes( frac_tot, colour=loc_reg))+stat_ecdf()+
  coord_cartesian(xlim=c(0,0.0005))+facet_wrap(~top40)

wilcox.test(subset(mer_freq, top40=="top40" & loc_reg=="freqCds_mem", select="freq_whole")[,1], 
            subset(mer_freq, top40=="top40"& loc_reg=="freqUtr3_cyt", select="freq_whole")[,1])
wilcox.test(subset(mer_freq, top40=="top40" & loc_reg=="freqCds_mem", select="freq_whole")[,1], 
            subset(mer_freq, top40=="top40"& loc_reg=="freqCds_cyt", select="freq_whole")[,1])
wilcox.test(subset(mer_freq, top40=="top40" & loc_reg=="freqUtr3_mem", select="freq_whole")[,1], 
            subset(mer_freq, top40=="top40"& loc_reg=="freqUtr3_cyt", select="freq_whole")[,1])
wilcox.test(subset(mer_freq, top40=="top40" & loc_reg=="freqCds_mem", select="freq_whole")[,1], 
            subset(mer_freq, top40=="top40"& loc_reg=="freqUtr3_cyt", select="freq_whole")[,1])

memCds<-subset(mer_freq, loc_reg=="freqCds_mem")
cytCds<-subset(mer_freq, loc_reg=="freqCds_cyt")
colnames(memCds)<-paste0(colnames(memCds), "_memCds")
colnames(cytCds)<-paste0(colnames(cytCds), "_cytCds")

smer<-merge(memCds, cytCds, by.x="kmer_memCds", by.y="kmer_cytCds")
ggplot(smer, aes(freq_whole_memCds-freq_whole_cytCds, frac_tot_memCds))+
  geom_text(aes(label=kmer_memCds, colour=labs_memCds), size=2)+
  geom_vline(xintercept=0,lty=2)+theme(legend.position = "none")

ggplot(smer, aes(freq_whole_memCds-freq_whole_cytCds, frac_tot_cytCds))+
  geom_text(aes(label=kmer_memCds, colour=labs_memCds), size=2)+
  geom_vline(xintercept=0,lty=2)+theme(legend.position = "none")

ggplot(smer, aes(freq_whole_memCds,freq_whole_cytCds))+
  geom_text(aes(label=kmer_memCds, colour=labs_memCds), size=2)+
  geom_abline(slope=1,lty=2)+theme(legend.position = "none")

memCds<-subset(mer_freq, loc_reg=="freqCds_mem")
cytUtr3<-subset(mer_freq, loc_reg=="freqUtr3_cyt")
colnames(memCds)<-paste0(colnames(memCds), "_memCds")
colnames(cytUtr3)<-paste0(colnames(cytUtr3), "_cytUtr3")

smer<-merge(memCds, cytUtr3, by.x="kmer_memCds", by.y="kmer_cytUtr3")
ggplot(smer, aes(freq_whole_memCds-freq_whole_cytUtr3, frac_tot_memCds))+
  geom_text(aes(label=kmer_memCds, colour=labs_memCds), size=2)+
  geom_vline(xintercept=0,lty=2)+theme(legend.position = "none")

ggplot(smer, aes(freq_whole_memCds-freq_whole_cytUtr3, frac_tot_cytUtr3))+
  geom_text(aes(label=kmer_memCds, colour=labs_memCds), size=2)+
  geom_vline(xintercept=0,lty=2)+theme(legend.position = "none")

ggplot(smer, aes(freq_whole_memCds,freq_whole_cytUtr3))+
  geom_text(aes(label=kmer_memCds, colour=labs_memCds), size=2)+
  geom_abline(slope=1,lty=2)+theme(legend.position = "none")


memCds_bound<-subset(countmers, localization_cat=="membrane" & tc_region=="cds")
cytUtr3_bound<-subset(countmers, localization_cat=="cytosolic" & tc_region=="utr3")
colnames(memCds_bound)<-paste0(colnames(memCds_bound), "_memCds")
colnames(cytUtr3_bound)<-paste0(colnames(cytUtr3_bound), "_cytUtr3")

smer<-merge(memCds_bound, cytUtr3_bound, by.x="seqs_memCds", by.y="seqs_cytUtr3")
ggplot(smer, aes(frac_tot_memCds, frac_tot_cytUtr3))+
  geom_text(aes(label=seqs_memCds), size=2)+
  geom_abline(slope=1)
ggplot(smer, aes(frac_tot_memCds, frac_tot_cytUtr3))+
  geom_point()+
  geom_abline(slope=1)

memCds_bound<-subset(countmers, localization_cat=="membrane" & tc_region=="cds")
cytCds_bound<-subset(countmers, localization_cat=="cytosolic" & tc_region=="cds")
colnames(memCds_bound)<-paste0(colnames(memCds_bound), "_memCds")
colnames(cytCds_bound)<-paste0(colnames(cytCds_bound), "_cytCds")

smer<-merge(memCds_bound, cytCds_bound, by.x="seqs_memCds", by.y="seqs_cytCds")
ggplot(smer, aes(frac_tot_memCds, frac_tot_cytCds))+
  geom_text(aes(label=seqs_memCds), size=2)+
  geom_abline(slope=1)

#9-mer per gene

ag1_cds<-dcast(subset(mer, tc_region=="cds", select=c("transcript_id","localization_cat", "seqs", "norm_tc_num1")),transcript_id+localization_cat~seqs,sum,value.var="norm_tc_num1")
ag2_cds<-dcast(subset(mer, tc_region=="cds", select=c("transcript_id","localization_cat", "seqs", "norm_tc_num2")), transcript_id+localization_cat~seqs,sum,value.var="norm_tc_num2")

ag1_utr3<-dcast(subset(mer, tc_region=="utr3", select=c("transcript_id","localization_cat", "seqs", "norm_tc_num1")), transcript_id+localization_cat~seqs,sum,value.var="norm_tc_num1")
ag2_utr3<-dcast(subset(mer, tc_region=="utr3", select=c("transcript_id","localization_cat", "seqs", "norm_tc_num2")), transcript_id+localization_cat~seqs,sum,value.var="norm_tc_num2")

ag1_cds_len<-which(len$transcript %in% ag1_cds$transcript_id)
ag2_cds_len<-which(len$transcript %in% ag2_cds$transcript_id)
ag1_utr3_len<-which(len$transcript %in% ag1_utr3$transcript_id)
ag2_utr3_len<-which(len$transcript %in% ag2_utr3$transcript_id)

ag1_cds_len<-len[ag1_cds_len, c("transcript", "l_cds")]
ag2_cds_len<-len[ag2_cds_len, c("transcript", "l_cds")]
ag1_utr3_len<-len[ag1_utr3_len, c("transcript", "l_utr3")]
ag2_utr3_len<-len[ag2_utr3_len, c("transcript", "l_utr3")]

ag1_cds<-cbind(ag1_cds[,1:2],ag1_cds[,3:ncol(ag1_cds)]/ag1_cds_len$l_cds*1e3)
ag2_cds<-cbind(ag2_cds[,1:2],ag2_cds[,3:ncol(ag2_cds)]/ag2_cds_len$l_cds*1e3)
ag1_utr3<-cbind(ag1_utr3[,1:2],ag1_utr3[,3:ncol(ag1_utr3)]/ag1_utr3_len$l_utr3*1e3)
ag2_utr3<-cbind(ag2_utr3[,1:2],ag2_utr3[,3:ncol(ag2_utr3)]/ag2_utr3_len$l_utr3*1e3)

ag1_cds_tpm<-which(mas$transcript %in% ag1_cds$transcript_id)
ag2_cds_tpm<-which(mas$transcript %in% ag2_cds$transcript_id)
ag1_utr3_tpm<-which(mas$transcript %in% ag1_utr3$transcript_id)
ag2_utr3_tpm<-which(mas$transcript %in% ag2_utr3$transcript_id)

ag1_cds_tpm<-mas[ag1_cds_tpm, c("transcript", "tpm_cutoff")]
ag1_cds_tpm<-ag1_cds_tpm[order(match(ag1_cds_tpm$transcript, ag1_cds$transcript_id)),]
ag2_cds_tpm<-mas[ag2_cds_tpm, c("transcript", "tpm_cutoff")]
ag2_cds_tpm<-ag2_cds_tpm[order(match(ag2_cds_tpm$transcript, ag2_cds$transcript_id)),]
ag1_utr3_tpm<-mas[ag1_utr3_tpm, c("transcript", "tpm_cutoff")]
ag1_utr3_tpm<-ag1_utr3_tpm[order(match(ag1_utr3_tpm$transcript, ag1_utr3$transcript_id)),]
ag2_utr3_tpm<-mas[ag2_utr3_tpm, c("transcript", "tpm_cutoff")]
ag2_utr3_tpm<-ag2_utr3_tpm[order(match(ag2_utr3_tpm$transcript, ag2_utr3$transcript_id)),]

ag1_cds<-cbind(ag1_cds[,1:2],ag1_cds[,3:ncol(ag1_cds)]/ag1_cds_tpm$tpm_cutoff)
ag2_cds<-cbind(ag2_cds[,1:2],ag2_cds[,3:ncol(ag2_cds)]/ag2_cds_tpm$tpm_cutoff)
ag1_utr3<-cbind(ag1_utr3[,1:2],ag1_utr3[,3:ncol(ag1_utr3)]/ag1_utr3_tpm$tpm_cutoff)
ag2_utr3<-cbind(ag2_utr3[,1:2],ag2_utr3[,3:ncol(ag2_utr3)]/ag2_utr3_tpm$tpm_cutoff)

ag1_cds$region<-"cds1"
ag2_cds$region<-"cds2"
ag1_utr3$region<-"utr3_1"
ag2_utr3$region<-"utr3_2"

p1_cds<-melt(ag1_cds, measure.vars = colnames(ag1_cds)[3:(ncol(ag1_cds)-1)], 
             id.vars = c("transcript_id", "localization_cat", "region"), value.name = "xl", variable.name = "kmer")
p2_cds<-melt(ag2_cds, measure.vars = colnames(ag2_cds)[3:(ncol(ag2_cds)-1)], 
             id.vars = c("transcript_id", "localization_cat", "region"), value.name = "xl", variable.name = "kmer")
p1_utr3<-melt(ag1_utr3, measure.vars = colnames(ag1_utr3)[3:(ncol(ag1_utr3)-1)], 
              id.vars = c("transcript_id", "localization_cat", "region"), value.name = "xl", variable.name = "kmer")
p2_utr3<-melt(ag2_utr3, measure.vars = colnames(ag2_utr3)[3:(ncol(ag2_utr3)-1)], 
              id.vars = c("transcript_id", "localization_cat", "region"), value.name = "xl", variable.name = "kmer")

plot<-rbind(p1_cds, p2_cds, p1_utr3, p2_utr3)
plot<-subset(plot, xl>0)
plot$log2xl<-ifelse(plot$xl<=0, NA, log2(plot$xl))
plot$count<-ifelse(plot$xl>0, 1, NA)
meds<-aggregate(log2xl~kmer+localization_cat+region, data=plot, median, na.rm=T)
meds<-meds[order(meds$log2xl, decreasing=T),]
meds$id<-paste0(meds$kmer,"_", meds$localization_cat,"_", meds$region)
freq<-aggregate(count~kmer+localization_cat+region, data=plot, sum)
freq<-freq[order(freq$count, decreasing=T),]
freq$id<-paste0(freq$kmer,"_", freq$localization_cat,"_", freq$region)

ssr<-merge(meds, freq, by="id")
ggplot(ssr, aes(log2xl, count, fill=region.x ))+geom_text(aes(label=kmer.x), size=2)+facet_wrap(~localization_cat.x)
ssr<-subset(ssr, count>20)
ssr<-ssr[order(ssr$count, decreasing=T),]

top40<-subset(ssr, localization_cat.x=="membrane" & region.x=="cds1")$kmer.x[1:10]
plot_sub<-which(plot$kmer %in% top40)
plot_sub<-plot[plot_sub,]
plot_sub$kmer<-factor(plot_sub$kmer, levels=rev(top40))
ggplot(plot_sub, aes(kmer, log2xl, fill=region))+geom_boxplot()+facet_wrap(~localization_cat)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
ssr<-merge(meds, freq, by="id")
ssr_sub<-which(ssr$kmer.x %in% top40)
ssr_sub<-ssr[ssr_sub,]
ssr_sub$kmer.x<-factor(ssr_sub$kmer.x, levels=rev(top40))
ggplot(ssr_sub, aes(kmer.x, count, fill=region.x))+geom_bar(stat = "identity", position = "dodge")+facet_wrap(~localization_cat.x)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))

#z-scores
setwd("D:/google_drive/hdlbp/kmers_mem_cyt/")
# setwd("E:/work/hdlbp/kmers_mem_cyt/")

files<-list.files(getwd(), pattern="whole")
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

ggplot(zs, aes(zscores_CdsMem_CdsCyt, colour=length))+stat_ecdf()+facet_wrap(~type, scales="free")+coord_cartesian(xlim=c(-3,13))
ggplot(zs, aes(zscores_CdsMem_Utr3Cyt, colour=length))+stat_ecdf()+facet_wrap(~type, scales="free")+coord_cartesian(xlim=c(-5,7))
ggplot(zs, aes(zscores_CdsCyt_Utr3Cyt, colour=length))+stat_ecdf()+facet_wrap(~type, scales="free")+coord_cartesian(xlim=c(-5,7))
ggplot(zs, aes(zscores_Utr3Mem_Utr3Cyt, colour=length))+stat_ecdf()+facet_wrap(~type, scales="free")+coord_cartesian(xlim=c(-5,4))

ggplot(subset(zs, type=="top40"), aes(zscores_CdsUtr3Mem_CdsUtr3Cyt, colour=length))+stat_ecdf()+facet_wrap(~type, scales="free")+
  coord_cartesian(xlim=c(-3,15))
ggplot(subset(zs, type=="top40"), aes(zscores_CdsMem_CdsCyt, colour=length))+stat_ecdf()+facet_wrap(~type, scales="free")+
  coord_cartesian(xlim=c(-3,20))
ggplot(subset(zs, type=="top40"), aes(zscores_CdsMem_Utr3Cyt, colour=length))+stat_ecdf()+facet_wrap(~type, scales="free")+
  coord_cartesian(xlim=c(-3,15))
ggplot(subset(zs, type=="top40"), aes(zscores_CdsCyt_Utr3Cyt, colour=length))+stat_ecdf()+facet_wrap(~type, scales="free")+
  coord_cartesian(xlim=c(-3,15))
ggplot(subset(zs, type=="top40"), aes(zscores_Utr3Mem_Utr3Cyt, colour=length))+stat_ecdf()+facet_wrap(~type, scales="free")+
  coord_cartesian(xlim=c(-5,5))
ggplot(subset(zs, type=="top40"), aes(zscores_Utr3Mem_Utr3Cyt, colour=length))+geom_density()+facet_wrap(~type, scales="free")+
  coord_cartesian(xlim=c(-5,5))

wilcox.test(subset(zs, type=="top40" & length=="4mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="5mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="4mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="6mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="4mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="7mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="4mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="8mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="4mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="9mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="4mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="10mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="4mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="11mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="4mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="12mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="5mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="6mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="5mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="7mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="5mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="8mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="5mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="9mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="5mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="10mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="5mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="11mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="5mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="12mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="6mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="7mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="6mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="8mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="6mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="9mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="6mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="10mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="6mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="11mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="6mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="12mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value

wilcox.test(subset(zs, type=="top40" & length=="7mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="8mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="7mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="9mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="7mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="10mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="7mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="11mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="7mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="12mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value

wilcox.test(subset(zs, type=="top40" & length=="8mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="9mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="8mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="10mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="8mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="11mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="8mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="12mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value

wilcox.test(subset(zs, type=="top40" & length=="9mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="10mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="9mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="11mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="9mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="12mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value

wilcox.test(subset(zs, type=="top40" & length=="10mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="11mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value
wilcox.test(subset(zs, type=="top40" & length=="10mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="12mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value

wilcox.test(subset(zs, type=="top40" & length=="11mers", select="zscores_CdsCyt_Utr3Cyt")[,1], 
            subset(zs, type=="top40" & length=="12mers", select="zscores_CdsCyt_Utr3Cyt")[,1])$p.value


wilcox.test(subset(zs, type=="top40" & length=="4mers", select="zscores_CdsCdsCyt_CdsUtr3Cyt")[,1], 
            subset(zs, type=="other" & length=="4mers", select="zscores_CdsCdsCyt_CdsUtr3Cyt")[,1])
wilcox.test(subset(zs, type=="top40" & length=="5mers", select="zscores_CdsCdsCyt_CdsUtr3Cyt")[,1], 
            subset(zs, type=="other" & length=="5mers", select="zscores_CdsCdsCyt_CdsUtr3Cyt")[,1])
wilcox.test(subset(zs, type=="top40" & length=="6mers", select="zscores_CdsCdsCyt_CdsUtr3Cyt")[,1], 
            subset(zs, type=="other" & length=="6mers", select="zscores_CdsCdsCyt_CdsUtr3Cyt")[,1])

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


##5mer crosslinked
mer<-dat

mer$count<-ifelse(is.na(mer$localization_cat), NA,
                  ifelse(mer$norm_tc_num1>2 & mer$norm_tc_num2>2 & mer$TPM_transcript>=10 , 1, NA))
mer<-subset(mer, !is.na(count))


mer$tc_seq_pos1_start<-ifelse(mer$tc_start-3<1, NA, mer$tc_start-3)
mer$tc_seq_pos1_stop<-ifelse(mer$tc_start+1>mer$l_tr, NA, mer$tc_start+1)
mer$tc_seq_pos2_start<-ifelse(mer$tc_start-2<1, NA, mer$tc_start-2)
mer$tc_seq_pos2_stop<-ifelse(mer$tc_start+2>mer$l_tr, NA, mer$tc_start+2)
mer$tc_seq_pos3_start<-ifelse(mer$tc_start-1<1, NA, mer$tc_start-1)
mer$tc_seq_pos3_stop<-ifelse(mer$tc_start+3>mer$l_tr, NA, mer$tc_start+3)
mer$tc_seq_pos4_start<-ifelse(mer$tc_start<1, NA, mer$tc_start)
mer$tc_seq_pos4_stop<-ifelse(mer$tc_start+4>mer$l_tr, NA, mer$tc_start+4)
mer$tc_seq_pos5_start<-ifelse(mer$tc_start+1<1, NA, mer$tc_start+1)
mer$tc_seq_pos5_stop<-ifelse(mer$tc_start+5>mer$l_tr, NA, mer$tc_start+5)



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

mer<-melt(mer, measure.vars = colnames(mer)[44:48], id.vars = colnames(mer)[1:33], value.name = "seqs", variable.name = "tc_seq_pos")


# countmers_gene<-dcast(subset(mer, !is.na(count), select=c("transcript_id", "seq", "count")),transcript_id~seq,sum,value.var="count")

countmers<-aggregate(count~seqs+tc_region+localization_cat, sum, data=mer)

countmers$frac_tot<-countmers$count/sum(countmers$count)
agr<-aggregate(frac_tot~seqs, sum, data=countmers)
agr<-agr[order(agr$frac_tot, decreasing = T),]
agr<-agr[1:40, "seqs"]
rows<-which(countmers$seqs %in% agr)
countmers_plot<-countmers[rows,]
countmers_plot$seqs<-factor(countmers_plot$seqs, levels=rev(agr))
ggplot(subset(countmers_plot, tc_region!="utr5"), aes(seqs, frac_tot, fill=tc_region))+geom_bar(stat="identity", position="dodge",)+coord_flip()+facet_wrap(~localization_cat)+
  theme(text = element_text(size=8))

library(ggseqlogo)
ggplot(subset(countmers_plot, localization_cat=="cytosolic"))+geom_logo(gsub("T","U",as.character(countmers_plot$seqs)), method="bits", seq_type = "rna")+theme_logo()

#write.table(countmers, "./kmers_mem_cyt/bound-4mers-freqs-01032020.txt", quote=F, sep="\t", row.names = F)

##redo5mer enrichment in CDS/UTR3 mem/cyt
#run until tcCds


setwd("E:/Google Drive/hdlbp/")
# mas<-read.delim("hdlbp_master_table_with_classes.txt", header=T)
mas<-read.delim("hdlbp_master_table_with_classes_uniq_tsig.txt", header=T)
inf<-subset(mas, select=c("gene_id","Symbol",
                          "tpm_cutoff","tc_CDS_norm","localization_cat",
                          "mean_te_293","loc_tar_CDS","tc_CDS_norm_cat",
                          "tc_transcript_norm",
                          "log2FoldChange.mem.cyt.293",
                          "log2FoldChange.mem.cyt.KO.293",
                          "log2FoldChange.ribo.rna.KO.WT"))


mem<-subset(mas, localization_cat=="membrane" & gene_biotype=="protein_coding" & tpm_cutoff>=10, select="transcript")[,1]
cyt<-subset(mas, localization_cat=="cytosolic" & gene_biotype=="protein_coding" & tpm_cutoff>=10, select="transcript")[,1]


seqCds_mem<-which(names(seqCds) %in% mem)
seqCds_mem<-seqCds[seqCds_mem]
seqCds_cyt<-which(names(seqCds) %in% cyt)
seqCds_cyt<-seqCds[seqCds_cyt]

seqUtr3_mem<-which(names(seqUtr3) %in% mem)
seqUtr3_mem<-seqUtr3[seqUtr3_mem]
seqUtr3_cyt<-which(names(seqUtr3) %in% cyt)
seqUtr3_cyt<-seqUtr3[seqUtr3_cyt]
kmer<-5
randsCds_mem <- (Biostrings::oligonucleotideFrequency(seqCds_mem, width=kmer, step = 1, as.prob = F, with.labels = TRUE))
freqCds_mem<-colSums(randsCds_mem)/sum(colSums(randsCds_mem))

randsCds_cyt <- (Biostrings::oligonucleotideFrequency(seqCds_cyt, width=kmer, step = 1, as.prob = F, with.labels = TRUE))
freqCds_cyt<-colSums(randsCds_cyt)/sum(colSums(randsCds_cyt))

randsUtr3_mem <- (Biostrings::oligonucleotideFrequency(seqUtr3_mem, width=kmer, step = 1, as.prob = F, with.labels = TRUE))
freqUtr3_mem<-colSums(randsUtr3_mem)/sum(colSums(randsUtr3_mem))

randsUtr3_cyt <- (Biostrings::oligonucleotideFrequency(seqUtr3_cyt, width=kmer, step = 1, as.prob = F, with.labels = TRUE))
freqUtr3_cyt<-colSums(randsUtr3_cyt)/sum(colSums(randsUtr3_cyt))

freqCdsUtr3_mem<-(colSums(randsCds_mem)+colSums(randsUtr3_mem))/sum(colSums(randsCds_mem)+colSums(randsUtr3_mem))
freqCdsUtr3_cyt<-(colSums(randsCds_cyt)+colSums(randsUtr3_cyt))/sum(colSums(randsCds_cyt)+colSums(randsUtr3_cyt))


freqs<-data.frame(kmer=names(colSums(randsCds_mem)),
                  freqCds_mem=freqCds_mem,
                  freqCds_cyt=freqCds_cyt,
                  freqUtr3_mem=freqUtr3_mem,
                  freqUtr3_cyt=freqUtr3_cyt,
                  freqCdsUtr3_mem=freqCdsUtr3_mem,
                  freqCdsUtr3_cyt=freqCdsUtr3_cyt)
labs<-which(freqs$kmer %in% agr)
labs<-freqs[labs,]
labs$labs<-labs$kmer
freqs$labs<-NA
freqs<-rbind(freqs, labs)
ggplot(freqs, aes(freqCds_mem, freqCds_cyt))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()
ggplot(freqs, aes(freqCds_mem, freqUtr3_mem))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()
ggplot(freqs, aes(freqCds_mem, freqUtr3_mem, colour=labs))+geom_point()+geom_abline()
ggplot(freqs, aes(freqCds_mem, freqUtr3_cyt))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()
ggplot(freqs, aes(freqCdsUtr3_mem, freqCdsUtr3_cyt))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()
ggplot(freqs, aes(freqCdsUtr3_mem, freqCdsUtr3_cyt, colour=labs))+geom_point()+geom_abline()

ggplot(freqs, aes(freqCds_cyt, freqUtr3_cyt))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()
ggplot(freqs, aes(freqCds_cyt, freqUtr3_cyt, colour=labs))+geom_point()+geom_abline()
ggplot(freqs, aes(freqUtr3_mem, freqUtr3_cyt, colour=labs))+geom_point()+geom_abline()
ggplot(freqs, aes(freqCds_cyt, freqUtr3_cyt))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()

#write.table(freqs, "./kmers_mem_cyt/whole-5mers-freqs-01032020.txt", quote=F, sep="\t", row.names=F)

mer_freq<-melt(freqs, measure.vars = colnames(freqs)[grepl("freq", colnames(freqs))], id.vars="kmer", value.name = "freq_whole", variable.name = "loc_reg")
mer_freq$id<-paste0(mer_freq$kmer, "_", mer_freq$loc_reg)

countmers$id<-ifelse(countmers$tc_region=="cds" & countmers$localization_cat=="cytosolic", paste0(countmers$seq, "_freqCds_cyt"),
                     ifelse(countmers$tc_region=="cds" & countmers$localization_cat=="membrane", paste0(countmers$seq, "_freqCds_mem"),
                            ifelse(countmers$tc_region=="utr3" & countmers$localization_cat=="cytosolic", paste0(countmers$seq, "_freqUtr3_cyt"),
                                   ifelse(countmers$tc_region=="utr3" & countmers$localization_cat=="membrane", paste0(countmers$seq, "_freqUtr3_mem"),NA))))

mer_freq<-merge(mer_freq, countmers, by="id")

labs<-which(mer_freq$kmer %in% agr)
labs<-mer_freq[labs,]
labs$labs<-labs$kmer
mer_freq$labs<-NA
mer_freq<-rbind(mer_freq, labs)
ggplot(mer_freq,aes(frac_tot, freq_whole))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()+theme(legend.position = "none")+
  facet_wrap(~localization_cat+tc_region)

ggplot(mer_freq,aes(frac_tot, freq_whole))+geom_point()+geom_abline()+theme(legend.position = "none")+
  facet_wrap(~localization_cat+tc_region)

mer_freq$top40<-ifelse(is.na(mer_freq$labs), "other", "top40")
ggplot(mer_freq,aes( freq_whole, colour=top40))+geom_density()+
  facet_wrap(~localization_cat+tc_region)+coord_cartesian(xlim=c(0,0.02))

memCds<-subset(mer_freq, loc_reg=="freqCds_mem")
cytCds<-subset(mer_freq, loc_reg=="freqCds_cyt")
colnames(memCds)<-paste0(colnames(memCds), "_memCds")
colnames(cytCds)<-paste0(colnames(cytCds), "_cytCds")

smer<-merge(memCds, cytCds, by.x="kmer_memCds", by.y="kmer_cytCds")
ggplot(smer, aes(freq_whole_memCds-freq_whole_cytCds, frac_tot_memCds))+
  geom_text(aes(label=kmer_memCds, colour=labs_memCds), size=2)+
  geom_vline(xintercept=0,lty=2)+theme(legend.position = "none")

ggplot(smer, aes(freq_whole_memCds-freq_whole_cytCds, frac_tot_cytCds))+
  geom_text(aes(label=kmer_memCds, colour=labs_memCds), size=2)+
  geom_vline(xintercept=0,lty=2)+theme(legend.position = "none")


memCds_bound<-subset(countmers, localization_cat=="membrane" & tc_region=="cds")
cytUtr3_bound<-subset(countmers, localization_cat=="cytosolic" & tc_region=="utr3")
colnames(memCds_bound)<-paste0(colnames(memCds_bound), "_memCds")
colnames(cytUtr3_bound)<-paste0(colnames(cytUtr3_bound), "_cytUtr3")

smer<-merge(memCds_bound, cytUtr3_bound, by.x="seqs_memCds", by.y="seqs_cytUtr3")
ggplot(smer, aes(frac_tot_memCds, frac_tot_cytUtr3))+
  geom_text(aes(label=seqs_memCds), size=2)+
  geom_abline(slope=1)
ggplot(smer, aes(frac_tot_memCds, frac_tot_cytUtr3))+
  geom_point()+
  geom_abline(slope=1)

memCds_bound<-subset(countmers, localization_cat=="membrane" & tc_region=="cds")
cytCds_bound<-subset(countmers, localization_cat=="cytosolic" & tc_region=="cds")
colnames(memCds_bound)<-paste0(colnames(memCds_bound), "_memCds")
colnames(cytCds_bound)<-paste0(colnames(cytCds_bound), "_cytCds")

smer<-merge(memCds_bound, cytCds_bound, by.x="seqs_memCds", by.y="seqs_cytCds")
ggplot(smer, aes(frac_tot_memCds, frac_tot_cytCds))+
  geom_text(aes(label=seqs_memCds), size=2)+
  geom_abline(slope=1)

#5-mer per gene

ag1_cds<-dcast(subset(mer, tc_region=="cds", select=c("transcript_id","localization_cat", "seqs", "norm_tc_num1")),transcript_id+localization_cat~seqs,sum,value.var="norm_tc_num1")
ag2_cds<-dcast(subset(mer, tc_region=="cds", select=c("transcript_id","localization_cat", "seqs", "norm_tc_num2")), transcript_id+localization_cat~seqs,sum,value.var="norm_tc_num2")

ag1_utr3<-dcast(subset(mer, tc_region=="utr3", select=c("transcript_id","localization_cat", "seqs", "norm_tc_num1")), transcript_id+localization_cat~seqs,sum,value.var="norm_tc_num1")
ag2_utr3<-dcast(subset(mer, tc_region=="utr3", select=c("transcript_id","localization_cat", "seqs", "norm_tc_num2")), transcript_id+localization_cat~seqs,sum,value.var="norm_tc_num2")

ag1_cds_len<-which(len$transcript %in% ag1_cds$transcript_id)
ag2_cds_len<-which(len$transcript %in% ag2_cds$transcript_id)
ag1_utr3_len<-which(len$transcript %in% ag1_utr3$transcript_id)
ag2_utr3_len<-which(len$transcript %in% ag2_utr3$transcript_id)

ag1_cds_len<-len[ag1_cds_len, c("transcript", "l_cds")]
ag2_cds_len<-len[ag2_cds_len, c("transcript", "l_cds")]
ag1_utr3_len<-len[ag1_utr3_len, c("transcript", "l_utr3")]
ag2_utr3_len<-len[ag2_utr3_len, c("transcript", "l_utr3")]

ag1_cds<-cbind(ag1_cds[,1:2],ag1_cds[,3:ncol(ag1_cds)]/ag1_cds_len$l_cds*1e3)
ag2_cds<-cbind(ag2_cds[,1:2],ag2_cds[,3:ncol(ag2_cds)]/ag2_cds_len$l_cds*1e3)
ag1_utr3<-cbind(ag1_utr3[,1:2],ag1_utr3[,3:ncol(ag1_utr3)]/ag1_utr3_len$l_utr3*1e3)
ag2_utr3<-cbind(ag2_utr3[,1:2],ag2_utr3[,3:ncol(ag2_utr3)]/ag2_utr3_len$l_utr3*1e3)

ag1_cds_tpm<-which(mas$transcript %in% ag1_cds$transcript_id)
ag2_cds_tpm<-which(mas$transcript %in% ag2_cds$transcript_id)
ag1_utr3_tpm<-which(mas$transcript %in% ag1_utr3$transcript_id)
ag2_utr3_tpm<-which(mas$transcript %in% ag2_utr3$transcript_id)

ag1_cds_tpm<-mas[ag1_cds_tpm, c("transcript", "tpm_cutoff")]
ag1_cds_tpm<-ag1_cds_tpm[order(match(ag1_cds_tpm$transcript, ag1_cds$transcript_id)),]
ag2_cds_tpm<-mas[ag2_cds_tpm, c("transcript", "tpm_cutoff")]
ag2_cds_tpm<-ag2_cds_tpm[order(match(ag2_cds_tpm$transcript, ag2_cds$transcript_id)),]
ag1_utr3_tpm<-mas[ag1_utr3_tpm, c("transcript", "tpm_cutoff")]
ag1_utr3_tpm<-ag1_utr3_tpm[order(match(ag1_utr3_tpm$transcript, ag1_utr3$transcript_id)),]
ag2_utr3_tpm<-mas[ag2_utr3_tpm, c("transcript", "tpm_cutoff")]
ag2_utr3_tpm<-ag2_utr3_tpm[order(match(ag2_utr3_tpm$transcript, ag2_utr3$transcript_id)),]

ag1_cds<-cbind(ag1_cds[,1:2],ag1_cds[,3:ncol(ag1_cds)]/ag1_cds_tpm$tpm_cutoff)
ag2_cds<-cbind(ag2_cds[,1:2],ag2_cds[,3:ncol(ag2_cds)]/ag2_cds_tpm$tpm_cutoff)
ag1_utr3<-cbind(ag1_utr3[,1:2],ag1_utr3[,3:ncol(ag1_utr3)]/ag1_utr3_tpm$tpm_cutoff)
ag2_utr3<-cbind(ag2_utr3[,1:2],ag2_utr3[,3:ncol(ag2_utr3)]/ag2_utr3_tpm$tpm_cutoff)

ag1_cds$region<-"cds1"
ag2_cds$region<-"cds2"
ag1_utr3$region<-"utr3_1"
ag2_utr3$region<-"utr3_2"

p1_cds<-melt(ag1_cds, measure.vars = colnames(ag1_cds)[3:(ncol(ag1_cds)-1)], 
             id.vars = c("transcript_id", "localization_cat", "region"), value.name = "xl", variable.name = "kmer")
p2_cds<-melt(ag2_cds, measure.vars = colnames(ag2_cds)[3:(ncol(ag2_cds)-1)], 
             id.vars = c("transcript_id", "localization_cat", "region"), value.name = "xl", variable.name = "kmer")
p1_utr3<-melt(ag1_utr3, measure.vars = colnames(ag1_utr3)[3:(ncol(ag1_utr3)-1)], 
              id.vars = c("transcript_id", "localization_cat", "region"), value.name = "xl", variable.name = "kmer")
p2_utr3<-melt(ag2_utr3, measure.vars = colnames(ag2_utr3)[3:(ncol(ag2_utr3)-1)], 
              id.vars = c("transcript_id", "localization_cat", "region"), value.name = "xl", variable.name = "kmer")

plot<-rbind(p1_cds, p2_cds, p1_utr3, p2_utr3)
plot<-subset(plot, xl>0)
plot$log2xl<-ifelse(plot$xl<=0, NA, log2(plot$xl))
plot$count<-ifelse(plot$xl>0, 1, NA)
meds<-aggregate(log2xl~kmer+localization_cat+region, data=plot, median, na.rm=T)
meds<-meds[order(meds$log2xl, decreasing=T),]
meds$id<-paste0(meds$kmer,"_", meds$localization_cat,"_", meds$region)
freq<-aggregate(count~kmer+localization_cat+region, data=plot, sum)
freq<-freq[order(freq$count, decreasing=T),]
freq$id<-paste0(freq$kmer,"_", freq$localization_cat,"_", freq$region)

ssr<-merge(meds, freq, by="id")
ggplot(ssr, aes(log2xl, count, fill=region.x ))+geom_text(aes(label=kmer.x), size=2)+facet_wrap(~localization_cat.x)
ssr<-subset(ssr, count>20)
ssr<-ssr[order(ssr$log2xl, decreasing=T),]

top40<-subset(ssr, localization_cat.x=="membrane" & region.x=="cds1")$kmer.x[1:10]
plot_sub<-which(plot$kmer %in% top40)
plot_sub<-plot[plot_sub,]
plot_sub$kmer<-factor(plot_sub$kmer, levels=rev(top40))
ggplot(plot_sub, aes(kmer, log2xl, fill=region))+geom_boxplot()+facet_wrap(~localization_cat)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
ssr<-merge(meds, freq, by="id")
ssr_sub<-which(ssr$kmer.x %in% top40)
ssr_sub<-ssr[ssr_sub,]
ssr_sub$kmer.x<-factor(ssr_sub$kmer.x, levels=rev(top40))
ggplot(ssr_sub, aes(kmer.x, count, fill=region.x))+geom_bar(stat = "identity", position = "dodge")+facet_wrap(~localization_cat.x)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))


##6mer crosslinked
mer<-dat

mer$count<-ifelse(is.na(mer$localization_cat), NA,
                  ifelse(mer$norm_tc_num1>2 & mer$norm_tc_num2>2 & mer$TPM_transcript>=10 , 1, NA))
mer<-subset(mer, !is.na(count))


mer$tc_seq_pos1_start<-ifelse(mer$tc_start-4<1, NA, mer$tc_start-4)
mer$tc_seq_pos1_stop<-ifelse(mer$tc_start+1>mer$l_tr, NA, mer$tc_start+1)
mer$tc_seq_pos2_start<-ifelse(mer$tc_start-3<1, NA, mer$tc_start-3)
mer$tc_seq_pos2_stop<-ifelse(mer$tc_start+2>mer$l_tr, NA, mer$tc_start+2)
mer$tc_seq_pos3_start<-ifelse(mer$tc_start-2<1, NA, mer$tc_start-2)
mer$tc_seq_pos3_stop<-ifelse(mer$tc_start+3>mer$l_tr, NA, mer$tc_start+3)
mer$tc_seq_pos4_start<-ifelse(mer$tc_start-1<1, NA, mer$tc_start-1)
mer$tc_seq_pos4_stop<-ifelse(mer$tc_start+4>mer$l_tr, NA, mer$tc_start+4)
mer$tc_seq_pos5_start<-ifelse(mer$tc_start<1, NA, mer$tc_start)
mer$tc_seq_pos5_stop<-ifelse(mer$tc_start+5>mer$l_tr, NA, mer$tc_start+5)
mer$tc_seq_pos6_start<-ifelse(mer$tc_start+1<1, NA, mer$tc_start+1)
mer$tc_seq_pos6_stop<-ifelse(mer$tc_start+6>mer$l_tr, NA, mer$tc_start+6)



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

mer<-melt(mer, measure.vars = colnames(mer)[46:51], id.vars = colnames(mer)[1:33], value.name = "seqs", variable.name = "tc_seq_pos")


# countmers_gene<-dcast(subset(mer, !is.na(count), select=c("transcript_id", "seq", "count")),transcript_id~seq,sum,value.var="count")

countmers<-aggregate(count~seqs+tc_region+localization_cat, sum, data=mer)

countmers$frac_tot<-countmers$count/sum(countmers$count)
agr<-aggregate(frac_tot~seqs, sum, data=countmers)
agr<-agr[order(agr$frac_tot, decreasing = T),]
agr<-agr[1:40, "seqs"]
rows<-which(countmers$seqs %in% agr)
countmers_plot<-countmers[rows,]
countmers_plot$seqs<-factor(countmers_plot$seqs, levels=rev(agr))
ggplot(subset(countmers_plot, tc_region!="utr5"), aes(seqs, frac_tot, fill=tc_region))+geom_bar(stat="identity", position="dodge",)+coord_flip()+facet_wrap(~localization_cat)+
  theme(text = element_text(size=8))

library(ggseqlogo)
ggplot(subset(countmers_plot, localization_cat=="cytosolic"))+geom_logo(gsub("T","U",as.character(countmers_plot$seqs)), method="bits", seq_type = "rna")+theme_logo()

#write.table(countmers, "./kmers_mem_cyt/bound-6mers-freqs-01032020.txt", quote=F, sep="\t", row.names = F)

##redo6mer enrichment in CDS/UTR3 mem/cyt
#run until tcCds


setwd("E:/Google Drive/hdlbp/")
# mas<-read.delim("hdlbp_master_table_with_classes.txt", header=T)
mas<-read.delim("hdlbp_master_table_with_classes_uniq_tsig.txt", header=T)
inf<-subset(mas, select=c("gene_id","Symbol",
                          "tpm_cutoff","tc_CDS_norm","localization_cat",
                          "mean_te_293","loc_tar_CDS","tc_CDS_norm_cat",
                          "tc_transcript_norm",
                          "log2FoldChange.mem.cyt.293",
                          "log2FoldChange.mem.cyt.KO.293",
                          "log2FoldChange.ribo.rna.KO.WT"))


mem<-subset(mas, localization_cat=="membrane" & gene_biotype=="protein_coding" & tpm_cutoff>=10, select="transcript")[,1]
cyt<-subset(mas, localization_cat=="cytosolic" & gene_biotype=="protein_coding" & tpm_cutoff>=10, select="transcript")[,1]


seqCds_mem<-which(names(seqCds) %in% mem)
seqCds_mem<-seqCds[seqCds_mem]
seqCds_cyt<-which(names(seqCds) %in% cyt)
seqCds_cyt<-seqCds[seqCds_cyt]

seqUtr3_mem<-which(names(seqUtr3) %in% mem)
seqUtr3_mem<-seqUtr3[seqUtr3_mem]
seqUtr3_cyt<-which(names(seqUtr3) %in% cyt)
seqUtr3_cyt<-seqUtr3[seqUtr3_cyt]
kmer<-6
randsCds_mem <- (Biostrings::oligonucleotideFrequency(seqCds_mem, width=kmer, step = 1, as.prob = F, with.labels = TRUE))
freqCds_mem<-colSums(randsCds_mem)/sum(colSums(randsCds_mem))

randsCds_cyt <- (Biostrings::oligonucleotideFrequency(seqCds_cyt, width=kmer, step = 1, as.prob = F, with.labels = TRUE))
freqCds_cyt<-colSums(randsCds_cyt)/sum(colSums(randsCds_cyt))

randsUtr3_mem <- (Biostrings::oligonucleotideFrequency(seqUtr3_mem, width=kmer, step = 1, as.prob = F, with.labels = TRUE))
freqUtr3_mem<-colSums(randsUtr3_mem)/sum(colSums(randsUtr3_mem))

randsUtr3_cyt <- (Biostrings::oligonucleotideFrequency(seqUtr3_cyt, width=kmer, step = 1, as.prob = F, with.labels = TRUE))
freqUtr3_cyt<-colSums(randsUtr3_cyt)/sum(colSums(randsUtr3_cyt))

freqCdsUtr3_mem<-(colSums(randsCds_mem)+colSums(randsUtr3_mem))/sum(colSums(randsCds_mem)+colSums(randsUtr3_mem))
freqCdsUtr3_cyt<-(colSums(randsCds_cyt)+colSums(randsUtr3_cyt))/sum(colSums(randsCds_cyt)+colSums(randsUtr3_cyt))


freqs<-data.frame(kmer=names(colSums(randsCds_mem)),
                  freqCds_mem=freqCds_mem,
                  freqCds_cyt=freqCds_cyt,
                  freqUtr3_mem=freqUtr3_mem,
                  freqUtr3_cyt=freqUtr3_cyt,
                  freqCdsUtr3_mem=freqCdsUtr3_mem,
                  freqCdsUtr3_cyt=freqCdsUtr3_cyt)
labs<-which(freqs$kmer %in% agr)
labs<-freqs[labs,]
labs$labs<-labs$kmer
freqs$labs<-NA
freqs<-rbind(freqs, labs)
ggplot(freqs, aes(freqCds_mem, freqCds_cyt))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()
ggplot(freqs, aes(freqCds_mem, freqUtr3_mem))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()
ggplot(freqs, aes(freqCds_mem, freqUtr3_mem, colour=labs))+geom_point()+geom_abline()
ggplot(freqs, aes(freqCds_mem, freqUtr3_cyt))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()
ggplot(freqs, aes(freqCdsUtr3_mem, freqCdsUtr3_cyt))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()
ggplot(freqs, aes(freqCdsUtr3_mem, freqCdsUtr3_cyt, colour=labs))+geom_point()+geom_abline()

ggplot(freqs, aes(freqCds_cyt, freqUtr3_cyt))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()
ggplot(freqs, aes(freqCds_cyt, freqUtr3_cyt, colour=labs))+geom_point()+geom_abline()
ggplot(freqs, aes(freqUtr3_mem, freqUtr3_cyt, colour=labs))+geom_point()+geom_abline()
ggplot(freqs, aes(freqCds_cyt, freqUtr3_cyt))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()

#write.table(freqs, "./kmers_mem_cyt/whole-6mers-freqs-01032020.txt", quote=F, sep="\t", row.names=F)

mer_freq<-melt(freqs, measure.vars = colnames(freqs)[grepl("freq", colnames(freqs))], id.vars="kmer", value.name = "freq_whole", variable.name = "loc_reg")
mer_freq$id<-paste0(mer_freq$kmer, "_", mer_freq$loc_reg)

countmers$id<-ifelse(countmers$tc_region=="cds" & countmers$localization_cat=="cytosolic", paste0(countmers$seq, "_freqCds_cyt"),
                     ifelse(countmers$tc_region=="cds" & countmers$localization_cat=="membrane", paste0(countmers$seq, "_freqCds_mem"),
                            ifelse(countmers$tc_region=="utr3" & countmers$localization_cat=="cytosolic", paste0(countmers$seq, "_freqUtr3_cyt"),
                                   ifelse(countmers$tc_region=="utr3" & countmers$localization_cat=="membrane", paste0(countmers$seq, "_freqUtr3_mem"),NA))))

mer_freq<-merge(mer_freq, countmers, by="id")

labs<-which(mer_freq$kmer %in% agr)
labs<-mer_freq[labs,]
labs$labs<-labs$kmer
mer_freq$labs<-NA
mer_freq<-rbind(mer_freq, labs)
ggplot(mer_freq,aes(frac_tot, freq_whole))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()+theme(legend.position = "none")+
  facet_wrap(~localization_cat+tc_region)

ggplot(mer_freq,aes(frac_tot, freq_whole))+geom_point()+geom_abline()+theme(legend.position = "none")+
  facet_wrap(~localization_cat+tc_region)

mer_freq$top40<-ifelse(is.na(mer_freq$labs), "other", "top40")
ggplot(mer_freq,aes( freq_whole, colour=top40))+geom_density()+
  facet_wrap(~localization_cat+tc_region)+coord_cartesian(xlim=c(0,0.02))

memCds<-subset(mer_freq, loc_reg=="freqCds_mem")
cytCds<-subset(mer_freq, loc_reg=="freqCds_cyt")
colnames(memCds)<-paste0(colnames(memCds), "_memCds")
colnames(cytCds)<-paste0(colnames(cytCds), "_cytCds")

smer<-merge(memCds, cytCds, by.x="kmer_memCds", by.y="kmer_cytCds")
ggplot(smer, aes(freq_whole_memCds-freq_whole_cytCds, frac_tot_memCds))+
  geom_text(aes(label=kmer_memCds, colour=labs_memCds), size=2)+
  geom_vline(xintercept=0,lty=2)+theme(legend.position = "none")

ggplot(smer, aes(freq_whole_memCds-freq_whole_cytCds, frac_tot_cytCds))+
  geom_text(aes(label=kmer_memCds, colour=labs_memCds), size=2)+
  geom_vline(xintercept=0,lty=2)+theme(legend.position = "none")


memCds_bound<-subset(countmers, localization_cat=="membrane" & tc_region=="cds")
cytUtr3_bound<-subset(countmers, localization_cat=="cytosolic" & tc_region=="utr3")
colnames(memCds_bound)<-paste0(colnames(memCds_bound), "_memCds")
colnames(cytUtr3_bound)<-paste0(colnames(cytUtr3_bound), "_cytUtr3")

smer<-merge(memCds_bound, cytUtr3_bound, by.x="seqs_memCds", by.y="seqs_cytUtr3")
ggplot(smer, aes(frac_tot_memCds, frac_tot_cytUtr3))+
  geom_text(aes(label=seqs_memCds), size=2)+
  geom_abline(slope=1)
ggplot(smer, aes(frac_tot_memCds, frac_tot_cytUtr3))+
  geom_point()+
  geom_abline(slope=1)

memCds_bound<-subset(countmers, localization_cat=="membrane" & tc_region=="cds")
cytCds_bound<-subset(countmers, localization_cat=="cytosolic" & tc_region=="cds")
colnames(memCds_bound)<-paste0(colnames(memCds_bound), "_memCds")
colnames(cytCds_bound)<-paste0(colnames(cytCds_bound), "_cytCds")

smer<-merge(memCds_bound, cytCds_bound, by.x="seqs_memCds", by.y="seqs_cytCds")
ggplot(smer, aes(frac_tot_memCds, frac_tot_cytCds))+
  geom_text(aes(label=seqs_memCds), size=2)+
  geom_abline(slope=1)

#6-mer per gene

ag1_cds<-dcast(subset(mer, tc_region=="cds", select=c("transcript_id","localization_cat", "seqs", "norm_tc_num1")),transcript_id+localization_cat~seqs,sum,value.var="norm_tc_num1")
ag2_cds<-dcast(subset(mer, tc_region=="cds", select=c("transcript_id","localization_cat", "seqs", "norm_tc_num2")), transcript_id+localization_cat~seqs,sum,value.var="norm_tc_num2")

ag1_utr3<-dcast(subset(mer, tc_region=="utr3", select=c("transcript_id","localization_cat", "seqs", "norm_tc_num1")), transcript_id+localization_cat~seqs,sum,value.var="norm_tc_num1")
ag2_utr3<-dcast(subset(mer, tc_region=="utr3", select=c("transcript_id","localization_cat", "seqs", "norm_tc_num2")), transcript_id+localization_cat~seqs,sum,value.var="norm_tc_num2")

ag1_cds_len<-which(len$transcript %in% ag1_cds$transcript_id)
ag2_cds_len<-which(len$transcript %in% ag2_cds$transcript_id)
ag1_utr3_len<-which(len$transcript %in% ag1_utr3$transcript_id)
ag2_utr3_len<-which(len$transcript %in% ag2_utr3$transcript_id)

ag1_cds_len<-len[ag1_cds_len, c("transcript", "l_cds")]
ag2_cds_len<-len[ag2_cds_len, c("transcript", "l_cds")]
ag1_utr3_len<-len[ag1_utr3_len, c("transcript", "l_utr3")]
ag2_utr3_len<-len[ag2_utr3_len, c("transcript", "l_utr3")]

ag1_cds<-cbind(ag1_cds[,1:2],ag1_cds[,3:ncol(ag1_cds)]/ag1_cds_len$l_cds*1e3)
ag2_cds<-cbind(ag2_cds[,1:2],ag2_cds[,3:ncol(ag2_cds)]/ag2_cds_len$l_cds*1e3)
ag1_utr3<-cbind(ag1_utr3[,1:2],ag1_utr3[,3:ncol(ag1_utr3)]/ag1_utr3_len$l_utr3*1e3)
ag2_utr3<-cbind(ag2_utr3[,1:2],ag2_utr3[,3:ncol(ag2_utr3)]/ag2_utr3_len$l_utr3*1e3)

ag1_cds_tpm<-which(mas$transcript %in% ag1_cds$transcript_id)
ag2_cds_tpm<-which(mas$transcript %in% ag2_cds$transcript_id)
ag1_utr3_tpm<-which(mas$transcript %in% ag1_utr3$transcript_id)
ag2_utr3_tpm<-which(mas$transcript %in% ag2_utr3$transcript_id)

ag1_cds_tpm<-mas[ag1_cds_tpm, c("transcript", "tpm_cutoff")]
ag1_cds_tpm<-ag1_cds_tpm[order(match(ag1_cds_tpm$transcript, ag1_cds$transcript_id)),]
ag2_cds_tpm<-mas[ag2_cds_tpm, c("transcript", "tpm_cutoff")]
ag2_cds_tpm<-ag2_cds_tpm[order(match(ag2_cds_tpm$transcript, ag2_cds$transcript_id)),]
ag1_utr3_tpm<-mas[ag1_utr3_tpm, c("transcript", "tpm_cutoff")]
ag1_utr3_tpm<-ag1_utr3_tpm[order(match(ag1_utr3_tpm$transcript, ag1_utr3$transcript_id)),]
ag2_utr3_tpm<-mas[ag2_utr3_tpm, c("transcript", "tpm_cutoff")]
ag2_utr3_tpm<-ag2_utr3_tpm[order(match(ag2_utr3_tpm$transcript, ag2_utr3$transcript_id)),]

ag1_cds<-cbind(ag1_cds[,1:2],ag1_cds[,3:ncol(ag1_cds)]/ag1_cds_tpm$tpm_cutoff)
ag2_cds<-cbind(ag2_cds[,1:2],ag2_cds[,3:ncol(ag2_cds)]/ag2_cds_tpm$tpm_cutoff)
ag1_utr3<-cbind(ag1_utr3[,1:2],ag1_utr3[,3:ncol(ag1_utr3)]/ag1_utr3_tpm$tpm_cutoff)
ag2_utr3<-cbind(ag2_utr3[,1:2],ag2_utr3[,3:ncol(ag2_utr3)]/ag2_utr3_tpm$tpm_cutoff)

ag1_cds$region<-"cds1"
ag2_cds$region<-"cds2"
ag1_utr3$region<-"utr3_1"
ag2_utr3$region<-"utr3_2"

p1_cds<-melt(ag1_cds, measure.vars = colnames(ag1_cds)[3:(ncol(ag1_cds)-1)], 
             id.vars = c("transcript_id", "localization_cat", "region"), value.name = "xl", variable.name = "kmer")
p2_cds<-melt(ag2_cds, measure.vars = colnames(ag2_cds)[3:(ncol(ag2_cds)-1)], 
             id.vars = c("transcript_id", "localization_cat", "region"), value.name = "xl", variable.name = "kmer")
p1_utr3<-melt(ag1_utr3, measure.vars = colnames(ag1_utr3)[3:(ncol(ag1_utr3)-1)], 
              id.vars = c("transcript_id", "localization_cat", "region"), value.name = "xl", variable.name = "kmer")
p2_utr3<-melt(ag2_utr3, measure.vars = colnames(ag2_utr3)[3:(ncol(ag2_utr3)-1)], 
              id.vars = c("transcript_id", "localization_cat", "region"), value.name = "xl", variable.name = "kmer")

plot<-rbind(p1_cds, p2_cds, p1_utr3, p2_utr3)
plot<-subset(plot, xl>0)
plot$log2xl<-ifelse(plot$xl<=0, NA, log2(plot$xl))
plot$count<-ifelse(plot$xl>0, 1, NA)
meds<-aggregate(log2xl~kmer+localization_cat+region, data=plot, median, na.rm=T)
meds<-meds[order(meds$log2xl, decreasing=T),]
meds$id<-paste0(meds$kmer,"_", meds$localization_cat,"_", meds$region)
freq<-aggregate(count~kmer+localization_cat+region, data=plot, sum)
freq<-freq[order(freq$count, decreasing=T),]
freq$id<-paste0(freq$kmer,"_", freq$localization_cat,"_", freq$region)

ssr<-merge(meds, freq, by="id")
ggplot(ssr, aes(log2xl, count, fill=region.x ))+geom_text(aes(label=kmer.x), size=2)+facet_wrap(~localization_cat.x)
ssr<-subset(ssr, count>20)
ssr<-ssr[order(ssr$log2xl, decreasing=T),]

top40<-subset(ssr, localization_cat.x=="membrane" & region.x=="cds1")$kmer.x[1:10]
plot_sub<-which(plot$kmer %in% top40)
plot_sub<-plot[plot_sub,]
plot_sub$kmer<-factor(plot_sub$kmer, levels=rev(top40))
ggplot(plot_sub, aes(kmer, log2xl, fill=region))+geom_boxplot()+facet_wrap(~localization_cat)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
ssr<-merge(meds, freq, by="id")
ssr_sub<-which(ssr$kmer.x %in% top40)
ssr_sub<-ssr[ssr_sub,]
ssr_sub$kmer.x<-factor(ssr_sub$kmer.x, levels=rev(top40))
ggplot(ssr_sub, aes(kmer.x, count, fill=region.x))+geom_bar(stat = "identity", position = "dodge")+facet_wrap(~localization_cat.x)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))




#7-mer per gene

ag1_cds<-dcast(subset(mer, tc_region=="cds", select=c("transcript_id","localization_cat", "seqs", "norm_tc_num1")),transcript_id+localization_cat~seqs,sum,value.var="norm_tc_num1")
ag2_cds<-dcast(subset(mer, tc_region=="cds", select=c("transcript_id","localization_cat", "seqs", "norm_tc_num2")), transcript_id+localization_cat~seqs,sum,value.var="norm_tc_num2")

ag1_utr3<-dcast(subset(mer, tc_region=="utr3", select=c("transcript_id","localization_cat", "seqs", "norm_tc_num1")), transcript_id+localization_cat~seqs,sum,value.var="norm_tc_num1")
ag2_utr3<-dcast(subset(mer, tc_region=="utr3", select=c("transcript_id","localization_cat", "seqs", "norm_tc_num2")), transcript_id+localization_cat~seqs,sum,value.var="norm_tc_num2")

ag1_cds_len<-which(len$transcript %in% ag1_cds$transcript_id)
ag2_cds_len<-which(len$transcript %in% ag2_cds$transcript_id)
ag1_utr3_len<-which(len$transcript %in% ag1_utr3$transcript_id)
ag2_utr3_len<-which(len$transcript %in% ag2_utr3$transcript_id)

ag1_cds_len<-len[ag1_cds_len, c("transcript", "l_cds")]
ag2_cds_len<-len[ag2_cds_len, c("transcript", "l_cds")]
ag1_utr3_len<-len[ag1_utr3_len, c("transcript", "l_utr3")]
ag2_utr3_len<-len[ag2_utr3_len, c("transcript", "l_utr3")]

ag1_cds<-cbind(ag1_cds[,1:2],ag1_cds[,3:ncol(ag1_cds)]/ag1_cds_len$l_cds*1e3)
ag2_cds<-cbind(ag2_cds[,1:2],ag2_cds[,3:ncol(ag2_cds)]/ag2_cds_len$l_cds*1e3)
ag1_utr3<-cbind(ag1_utr3[,1:2],ag1_utr3[,3:ncol(ag1_utr3)]/ag1_utr3_len$l_utr3*1e3)
ag2_utr3<-cbind(ag2_utr3[,1:2],ag2_utr3[,3:ncol(ag2_utr3)]/ag2_utr3_len$l_utr3*1e3)

ag1_cds_tpm<-which(mas$transcript %in% ag1_cds$transcript_id)
ag2_cds_tpm<-which(mas$transcript %in% ag2_cds$transcript_id)
ag1_utr3_tpm<-which(mas$transcript %in% ag1_utr3$transcript_id)
ag2_utr3_tpm<-which(mas$transcript %in% ag2_utr3$transcript_id)

ag1_cds_tpm<-mas[ag1_cds_tpm, c("transcript", "tpm_cutoff")]
ag1_cds_tpm<-ag1_cds_tpm[order(match(ag1_cds_tpm$transcript, ag1_cds$transcript_id)),]
ag2_cds_tpm<-mas[ag2_cds_tpm, c("transcript", "tpm_cutoff")]
ag2_cds_tpm<-ag2_cds_tpm[order(match(ag2_cds_tpm$transcript, ag2_cds$transcript_id)),]
ag1_utr3_tpm<-mas[ag1_utr3_tpm, c("transcript", "tpm_cutoff")]
ag1_utr3_tpm<-ag1_utr3_tpm[order(match(ag1_utr3_tpm$transcript, ag1_utr3$transcript_id)),]
ag2_utr3_tpm<-mas[ag2_utr3_tpm, c("transcript", "tpm_cutoff")]
ag2_utr3_tpm<-ag2_utr3_tpm[order(match(ag2_utr3_tpm$transcript, ag2_utr3$transcript_id)),]

ag1_cds<-cbind(ag1_cds[,1:2],ag1_cds[,3:ncol(ag1_cds)]/ag1_cds_tpm$tpm_cutoff)
ag2_cds<-cbind(ag2_cds[,1:2],ag2_cds[,3:ncol(ag2_cds)]/ag2_cds_tpm$tpm_cutoff)
ag1_utr3<-cbind(ag1_utr3[,1:2],ag1_utr3[,3:ncol(ag1_utr3)]/ag1_utr3_tpm$tpm_cutoff)
ag2_utr3<-cbind(ag2_utr3[,1:2],ag2_utr3[,3:ncol(ag2_utr3)]/ag2_utr3_tpm$tpm_cutoff)

ag1_cds$region<-"cds1"
ag2_cds$region<-"cds2"
ag1_utr3$region<-"utr3_1"
ag2_utr3$region<-"utr3_2"

p1_cds<-melt(ag1_cds, measure.vars = colnames(ag1_cds)[3:(ncol(ag1_cds)-1)], 
             id.vars = c("transcript_id", "localization_cat", "region"), value.name = "xl", variable.name = "kmer")
p2_cds<-melt(ag2_cds, measure.vars = colnames(ag2_cds)[3:(ncol(ag2_cds)-1)], 
             id.vars = c("transcript_id", "localization_cat", "region"), value.name = "xl", variable.name = "kmer")
p1_utr3<-melt(ag1_utr3, measure.vars = colnames(ag1_utr3)[3:(ncol(ag1_utr3)-1)], 
              id.vars = c("transcript_id", "localization_cat", "region"), value.name = "xl", variable.name = "kmer")
p2_utr3<-melt(ag2_utr3, measure.vars = colnames(ag2_utr3)[3:(ncol(ag2_utr3)-1)], 
              id.vars = c("transcript_id", "localization_cat", "region"), value.name = "xl", variable.name = "kmer")

plot<-rbind(p1_cds, p2_cds, p1_utr3, p2_utr3)
plot<-subset(plot, xl>0)
plot$log2xl<-ifelse(plot$xl<=0, NA, log2(plot$xl))
plot$count<-ifelse(plot$xl>0, 1, NA)
meds<-aggregate(log2xl~kmer+localization_cat+region, data=plot, median, na.rm=T)
meds<-meds[order(meds$log2xl, decreasing=T),]
meds$id<-paste0(meds$kmer,"_", meds$localization_cat,"_", meds$region)
freq<-aggregate(count~kmer+localization_cat+region, data=plot, sum)
freq<-freq[order(freq$count, decreasing=T),]
freq$id<-paste0(freq$kmer,"_", freq$localization_cat,"_", freq$region)

ssr<-merge(meds, freq, by="id")
ggplot(ssr, aes(log2xl, count, fill=region.x ))+geom_text(aes(label=kmer.x), size=2)+facet_wrap(~localization_cat.x)
ssr<-subset(ssr, count>20)
ssr<-ssr[order(ssr$count, decreasing=T),]

top40<-subset(ssr, localization_cat.x=="membrane" & region.x=="cds1")$kmer.x[1:10]
plot_sub<-which(plot$kmer %in% top40)
plot_sub<-plot[plot_sub,]
plot_sub$kmer<-factor(plot_sub$kmer, levels=rev(top40))
ggplot(plot_sub, aes(kmer, log2xl, fill=region))+geom_boxplot()+facet_wrap(~localization_cat)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
ssr<-merge(meds, freq, by="id")
ssr_sub<-which(ssr$kmer.x %in% top40)
ssr_sub<-ssr[ssr_sub,]
ssr_sub$kmer.x<-factor(ssr_sub$kmer.x, levels=rev(top40))
ggplot(ssr_sub, aes(kmer.x, count, fill=region.x))+geom_bar(stat = "identity", position = "dodge")+facet_wrap(~localization_cat.x)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))




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





