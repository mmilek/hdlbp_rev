# transcriptome fasta is needed for this figure
# figS6a

library(ggplot2)
library(reshape2)
library(corrplot)
library(Biostrings)
library(dplyr)
library(DESeq2)

len<-read.delim("data/transcript_cds_utr_lenght.txt", header=T)
rownames(len)<-len$transcript
len$cds_start<-len$l_utr5+1
len$cds_stop<-len$l_utr5+len$l_cds
len$frame_start<-len$l_utr5
len$frame_stop<-len$cds_start

# selecting highest expressed isoform
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

tc<-read.delim("data/reproducible.hdlbp.TCseq.bed", header=F)
colnames(tc)<-c("transcript_id","tc_start","tc_stop","tc_num1","all_reads1","tc_num2","all_reads2","seq")

fac<-estimateSizeFactorsForMatrix(tc[,c("all_reads1","all_reads2")])
tc$norm_tc_num1<-tc$tc_num1/fac[1]
tc$norm_tc_num2<-tc$tc_num2/fac[2]
ggplot(tc, aes(log2(norm_tc_num1),log2(norm_tc_num2)))+geom_point()+geom_abline(slope=1)
tc<-subset(tc, tc_num1/all_reads1!=1 & tc_num2/all_reads2!=1)

dat<-merge(tc, ham, by.x="transcript_id", by.y="transcript")
length(unique(dat$transcript_id))
length(unique(dat$gene_id))

dat$tc_from_start<-dat$tc_stop-dat$cds_start+1
dat<-subset(dat, tc_from_start<=l_cds & tc_from_start>0)
dat$codon_start<-ifelse((dat$tc_from_start-1)%%3==0, dat$tc_from_start, 
                        ifelse((dat$tc_from_start-2)%%3==0,dat$tc_from_start-1,
                               dat$tc_from_start-2))

# due to size limitations we cannot provide the transcriptome fasta file. please 
# provide this using local path. We used RSEM to produce 
# transcriptome fasta using gencode v19 annotation

fastapath<-"~/hdlbp_git/hg19bt1.transcripts.fa"
seqTrans <- Biostrings::readDNAStringSet(fastapath, format = "fasta", use.names = TRUE)

sub<- seqTrans[ham$transcript]

seqCds <- Biostrings::subseq(sub, start = ham$cds_start, end = ham$cds_stop)
seqUtr3<- Biostrings::subseq(sub, start = ham$cds_stop+1, end = ham$l_tr)
seqUtr5<- Biostrings::subseq(sub, start = 1, end = ham$l_utr5)

tcCds<-seqCds[unique(dat$transcript_id)]

dat$tc_codon<-as.character(Biostrings::subseq(tcCds[as.character(dat$transcript_id)],
                                              start=dat$codon_start,end=dat$codon_start+2))

dat$tc_codon_pos<-ifelse((dat$tc_from_start-dat$codon_start)==0,1,
                         ifelse((dat$tc_from_start-dat$codon_start)==1,2,3))


ag1<-dcast(dat[,c("gene_id","tc_codon","norm_tc_num1")], gene_id~tc_codon,sum,value.var="norm_tc_num1")
ag2<-dcast(dat[,c("gene_id","tc_codon","norm_tc_num2")], gene_id~tc_codon,sum,value.var="norm_tc_num2")

ag1$tc_codon_tot<-rowSums(ag1[2:ncol(ag1)],na.rm=T)
ag2$tc_codon_tot<-rowSums(ag2[2:ncol(ag2)],na.rm=T)
ag1$tc_codon_tot<-ifelse(ag1$tc_codon_tot==0,NA,ag1$tc_codon_tot)
ag2$tc_codon_tot<-ifelse(ag2$tc_codon_tot==0,NA,ag2$tc_codon_tot)

seq_freq <- (Biostrings::trinucleotideFrequency(tcCds, step = 3, as.prob = F, with.labels = TRUE))
rownames(seq_freq)<-names(tcCds)
seq_freq<-subset(seq_freq, select=colnames(ag1)[2:(ncol(ag1)-1)])
seq_freq[seq_freq==0]<-NA
fin1<-ag1[2:(ncol(ag1)-1)]/seq_freq*1e6

fin1<-cbind(gene_id=ag1$gene_id, fin1, tc_codon_tot=ag1$tc_codon_tot)
fin1<-subset(fin1, select=colnames(fin1)[!grepl("TGA",colnames(fin1)) & !grepl("TAA",colnames(fin1)) & !grepl("TAG",colnames(fin1))])
heat1<-melt(fin1, measure.vars = colnames(fin1)[2:35] , id.vars="gene_id")

fin2<-ag2[2:(ncol(ag2)-1)]/seq_freq*1e6
fin2<-cbind(gene_id=ag2$gene_id, fin2, tc_codon_tot=ag2$tc_codon_tot)
fin2<-subset(fin2, select=colnames(fin2)[!grepl("TGA",colnames(fin2)) & !grepl("TAA",colnames(fin2)) & !grepl("TAG",colnames(fin2))])
heat2<-melt(fin2, measure.vars = colnames(fin2)[2:35] , id.vars="gene_id")

colnames(fin1)[2:(ncol(fin1))]<-paste0("rep1_",colnames(fin1)[2:(ncol(fin1))])
colnames(fin2)[2:(ncol(fin2))]<-paste0("rep2_",colnames(fin2)[2:(ncol(fin2))])
fin<-merge(fin1, fin2, by="gene_id")

mas<-read.delim("data/hdlbp_master_table_with_classes_uniq.txt", header=T)
inf<-subset(mas, select=c("gene_id","Symbol",
                          "tpm_cutoff","tc_CDS_norm","localization_cat",
                          "mean_te_293","loc_tar_CDS","tc_CDS_norm_cat",
                          "tc_transcript_norm",
                          "log2FoldChange.mem.cyt.293",
                          "log2FoldChange.mem.cyt.KO.293",
                          "log2FoldChange.ribo.rna.KO.WT"))

fin<-merge(fin, inf, by="gene_id")

#correlation crosslinks per codon between replicates
ggplot(fin, aes(log2(rep1_tc_codon_tot), log2(rep2_tc_codon_tot)))+geom_point()+geom_abline()

ave<-(fin[,2:36]+fin[,37:71])/2
colnames(ave)<-gsub("rep1_","",colnames(ave))
ave<-cbind(gene_id=fin$gene_id, ave, fin[,72:ncol(fin)])
colnames(ave)[2:35]<-gsub("T","U",colnames(ave)[2:35])
cols<-data.frame(codon=colnames(ave)[2:35], seq=seq(1,length(colnames(ave)[2:35])))

cod<-read.delim("data/codon_table.txt")
codcols<-merge(cols, cod, by="codon")
codcols$cod_aa<-paste0(codcols$codon,";",codcols$aa2)
colnames(ave)[2:35]<-codcols$cod_aa
usg<-read.delim("trna/codonUsage.txt")
usg$codon<-gsub("T","U",usg$codon)
usg$codon<-paste0(usg$codon,";",usg$aa)
usg$dum<-1

ave<-subset(ave, tc_codon_tot>=10) ## thresholdat least 10 T-C per CDS
ggplot(ave, aes(localization_cat, log2(tc_codon_tot)))+geom_boxplot()
nas<-ave[,2:35]
nas[is.na(nas)]<-0
ave<-cbind(gene_id=ave[,1],nas,ave[36:ncol(ave)])
heat<-melt(ave, measure.vars = colnames(ave)[c(2:35)] , 
           id.vars=c("gene_id","tpm_cutoff","tc_CDS_norm",
                     "localization_cat","mean_te_293","loc_tar_CDS","tc_CDS_norm_cat",
                     "log2FoldChange.mem.cyt.293",
                     "log2FoldChange.mem.cyt.KO.293",
                     "log2FoldChange.ribo.rna.KO.WT",
                     "tc_transcript_norm",
                     "tc_codon_tot"))

# normalize codon crosslinks to expression level
heat$norm_cod_xl<-ifelse(heat$tpm_cutoff==0, NA, heat$value/heat$tpm_cutoff)


heat$gene_id<-factor(heat$gene_id, levels=unique(heat$gene_id[order(heat$log2FoldChange.mem.cyt.293, decreasing = T)]), ordered = T)

# figS6a
ggplot(heat, aes(gene_id, variable))+geom_tile(aes(fill=log2(norm_cod_xl)))+
  scale_fill_continuous(na.value = 'black')+theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
                                                  axis.text.x=element_blank(),
                                                  axis.ticks.x=element_blank())

sd<-ave[,c(1:39, 44, 45)]

# write.table(sd, "source_data/figS6a.txt", quote=F, sep="\t", row.names = F)

ave$dum<-1
ave$gene_id<-factor(ave$gene_id, levels=unique(ave$gene_id[order(ave$log2FoldChange.mem.cyt.293, decreasing = T)]), ordered = T)

# figS6a
ggplot(ave, aes(gene_id, dum))+geom_tile(aes(fill=log2FoldChange.mem.cyt.293))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
# figS6a
# total transcript
ggplot(ave, aes(gene_id, dum))+geom_tile(aes(fill=log2(tc_transcript_norm)))+scale_fill_continuous(na.value = 'black')+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
# CDS only
ggplot(ave, aes(gene_id, dum))+geom_tile(aes(fill=log2(tc_CDS_norm)))+scale_fill_continuous(na.value = 'black')+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
usg<-subset(usg, grepl("U", usg$codon))
usg<-subset(usg, !grepl("Stop", usg$codon))

# figS6a
ggplot(usg,aes(dum, codon))+geom_tile(aes(fill=codonUsage))

