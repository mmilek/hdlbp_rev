library(ggplot2)
library(reshape2)
library(corrplot)
library(Biostrings)
library(dplyr)
library(DESeq2)

# setwd("/Volumes/landthaler/pcp/projects/miha/HDLBP/all_clip_data/mapping_trans/")
# setwd("/Volumes/landthaler/pcp/projects/miha/HDLBP/all_clip_data/reclip/mapping_trans/")
# setwd("F:/landthaler/HDLBP/all_clip_data/reclip/mapping_trans/")
len<-read.delim("data/transcript_cds_utr_lenght.txt", header=T)
rownames(len)<-len$transcript
len$cds_start<-len$l_utr5+1
len$cds_stop<-len$l_utr5+len$l_cds
len$frame_start<-len$l_utr5
len$frame_stop<-len$cds_start

#make a bed file of all orfs for use in IGV
# bed<-len[rep(row.names(len), round(len$l_cds/3, digits=0)),]
# bed<-bed %>%
#   group_by(transcript) %>%
#   mutate(number = 1:n()) 
# bed[543,]
# bed$frame_start<-bed$frame_start+(bed$number-1)*3
# bed$frame_stop<-bed$frame_start+1
# bed$num<-1
# bed$num<-as.integer(bed$num)
# bed$strand<-"+"
# bed$orf<-"orf"
# bed<-bed[,c("transcript","frame_start","frame_stop","orf","num","strand")]
# write.table(bed, "transcript_orf.bed", quote=F, sep="\t", row.names=F, col.names = F)

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
# fin1<-ag1[2:(ncol(ag1)-1)]/ag1$tc_codon_tot*1e6
# fin1<-ag1[2:(ncol(ag1)-1)]/ag1$tc_codon_tot/seq_freq*1e6
fin1<-ag1[2:(ncol(ag1)-1)]/seq_freq*1e6

fin1<-cbind(gene_id=ag1$gene_id, fin1, tc_codon_tot=ag1$tc_codon_tot)
fin1<-subset(fin1, select=colnames(fin1)[!grepl("TGA",colnames(fin1)) & !grepl("TAA",colnames(fin1)) & !grepl("TAG",colnames(fin1))])
heat1<-melt(fin1, measure.vars = colnames(fin1)[2:35] , id.vars="gene_id")

# fin2<-ag2[2:(ncol(ag2)-1)]/ag2$tc_codon_tot*1e6
# fin2<-ag2[2:(ncol(ag2)-1)]/ag2$tc_codon_tot/seq_freq*1e6
fin2<-ag2[2:(ncol(ag2)-1)]/seq_freq*1e6
fin2<-cbind(gene_id=ag2$gene_id, fin2, tc_codon_tot=ag2$tc_codon_tot)
fin2<-subset(fin2, select=colnames(fin2)[!grepl("TGA",colnames(fin2)) & !grepl("TAA",colnames(fin2)) & !grepl("TAG",colnames(fin2))])
heat2<-melt(fin2, measure.vars = colnames(fin2)[2:35] , id.vars="gene_id")

# ggplot(heat1, aes(gene_id, variable))+geom_tile(aes(fill=log2(value)))+scale_fill_continuous(na.value = 'black')
# ggplot(heat2, aes(gene_id, variable))+geom_tile(aes(fill=log2(value)))+scale_fill_continuous(na.value = 'black')

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

ggplot(fin, aes(log2(rep1_tc_codon_tot), log2(rep2_tc_codon_tot)))+geom_point()+geom_abline()


ave<-(fin[,2:36]+fin[,37:71])/2
colnames(ave)<-gsub("rep1_","",colnames(ave))
ave<-cbind(gene_id=fin$gene_id, ave, fin[,72:ncol(fin)])
colnames(ave)[2:35]<-gsub("T","U",colnames(ave)[2:35])
cols<-data.frame(codon=colnames(ave)[2:35], seq=seq(1,length(colnames(ave)[2:35])))
# setwd("/Volumes/landthaler/pcp/projects/miha/HDLBP/all_clip_data/mapping_trans/")
# setwd("/Volumes/landthaler/pcp/projects/miha/HDLBP/all_clip_data/reclip/mapping_trans/")
setwd("F:/landthaler/HDLBP/all_clip_data/reclip//mapping_trans/")
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

#normalize codon crosslinks to expression level
heat$norm_cod_xl<-ifelse(heat$tpm_cutoff==0, NA, heat$value/heat$tpm_cutoff)


heat$gene_id<-factor(heat$gene_id, levels=unique(heat$gene_id[order(heat$log2FoldChange.mem.cyt.293, decreasing = T)]), ordered = T)
ggplot(heat, aes(gene_id, variable))+geom_tile(aes(fill=log2(value)))+
  scale_fill_continuous(limits=c(0,20),na.value = 'black')+theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
                                                                 axis.text.x=element_blank(),
                                                                 axis.ticks.x=element_blank())
ggplot(heat, aes(gene_id, variable))+geom_tile(aes(fill=log2(value/tpm_cutoff)))+
  scale_fill_continuous(na.value = 'black')+theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
                                                                 axis.text.x=element_blank(),
                                                                 axis.ticks.x=element_blank())
ggplot(heat, aes(gene_id, variable))+geom_tile(aes(fill=log2(value/tpm_cutoff)))+
  scale_fill_continuous(na.value = 'black')+theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
                                                  axis.text.x=element_blank(),
                                                  axis.ticks.x=element_blank())
ggplot(heat, aes(gene_id, variable))+geom_tile(aes(fill=log2(norm_cod_xl)))+
  scale_fill_continuous(na.value = "black")+theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
                                                  axis.text.x=element_blank(),
                                                  axis.ticks.x=element_blank())
ggplot(heat, aes(gene_id, variable))+geom_tile(aes(fill=log2(norm_cod_xl)))+
  scale_fill_continuous()+theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
                                                  axis.text.x=element_blank(),
                                                  axis.ticks.x=element_blank())



plot<-subset(heat, grepl("UUC", heat$variable))
ggplot(plot, aes(variable, log2(norm_cod_xl), fill=tc_CDS_norm_cat))+geom_boxplot()+coord_flip()
ggplot(subset(heat, !is.na(tc_CDS_norm_cat) & tc_CDS_norm_cat!="nontarget"), aes(variable, log2(value), fill=tc_CDS_norm_cat))+geom_boxplot()+coord_flip()
ggplot(subset(heat, !is.na(localization_cat)), aes(variable, log2(value), fill=localization_cat))+geom_boxplot()+coord_flip()

meds<-aggregate(norm_cod_xl~variable+localization_cat, data=heat, mean, na.rm=T)
meds<-meds[order(meds$norm_cod_xl, decreasing = F),]

meds<-aggregate(value~variable, data=heat, median, na.rm=T)
meds<-meds[order(meds$value, decreasing = F),]

perTrans<-aggregate(value~gene_id+localization_cat, data=heat, mean, na.rm=T)
ggplot(perTrans, aes(localization_cat, log2(value)))+geom_boxplot()

heat$variable<-factor(heat$variable, levels=meds$variable)
ggplot(heat, aes(variable, log2(value/tpm_cutoff)))+geom_boxplot()+coord_flip()+geom_hline(yintercept = 8)
ggplot(subset(heat, !is.na(localization_cat)), aes(variable, log2(value/tpm_cutoff), fill=localization_cat))+geom_boxplot()+coord_flip()+geom_hline(yintercept = 8)

ggplot(subset(heat, !is.na(localization_cat)), aes(x=reorder( variable, log2(norm_cod_xl), FUN=median), y=log2(norm_cod_xl), fill=localization_cat))+geom_boxplot()+coord_flip()
ggplot(subset(heat, !is.na(localization_cat)), aes(x=reorder( variable, log2(norm_cod_xl), FUN=median), y=log2(norm_cod_xl)))+geom_boxplot()+coord_flip()

ggplot(subset(heat, !is.na(tc_CDS_norm_cat)), aes(variable, log2(value), fill=tc_CDS_norm_cat))+geom_boxplot()+coord_flip()
ave$dum<-1
ave$gene_id<-factor(ave$gene_id, levels=unique(ave$gene_id[order(ave$log2FoldChange.mem.cyt.293, decreasing = T)]), ordered = T)

ggplot(ave, aes(gene_id, dum))+geom_tile(aes(fill=log2FoldChange.mem.cyt.293))+
  theme(axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank())
ggplot(ave, aes(gene_id, dum))+geom_tile(aes(fill=log2(tc_transcript_norm)))+scale_fill_continuous(na.value = 'black')+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggplot(ave, aes(gene_id, dum))+geom_tile(aes(fill=log2(tc_CDS_norm)))+scale_fill_continuous(na.value = 'black')+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
usg<-subset(usg, grepl("U", usg$codon))
usg<-subset(usg, !grepl("Stop", usg$codon))
ggplot(usg,aes(dum, codon))+geom_tile(aes(fill=codonUsage))


heat$gene_id<-factor(heat$gene_id, levels=unique(heat$gene_id[order(heat$tc_transcript_norm, decreasing = T)]), ordered = T)
ggplot(heat, aes(gene_id, variable))+geom_tile(aes(fill=log2(value)))+
  scale_fill_continuous(limits=c(0,20),na.value = 'black')+theme(axis.title.x=element_blank(),
                                                                 axis.text.x=element_blank(),
                                                                 axis.ticks.x=element_blank())

ggplot(heat, aes(gene_id, variable))+geom_tile(aes(fill=log2(value/tpm_cutoff)))+
  scale_fill_continuous(limits=c(0,25),na.value = 'grey')+theme(axis.title.x=element_blank(),
                                                                 axis.text.x=element_blank(),
                                                                 axis.ticks.x=element_blank())
ggplot(heat, aes(gene_id, variable))+geom_tile(aes(fill=log2(value/tpm_cutoff)))+
  scale_fill_continuous(limits=c(0,25),na.value = 'black')+facet_wrap(~localization_cat, scales="free_x")+theme(axis.title.x=element_blank(),
                                                                axis.text.x=element_blank(),
                                                                axis.ticks.x=element_blank())

# hi<-subset(heat, loc_tar_CDS=="membrane_tc>10.82 & tc<181.46")
# mid<-subset(heat, loc_tar_CDS=="membrane_tc>2.81 & tc<10.82")
# lo<-subset(heat, loc_tar_CDS=="membrane_tc<2.81")
hi<-subset(heat, loc_tar_CDS=="membrane_tc>5.56 & tc<65.26")
mid<-subset(heat, loc_tar_CDS=="membrane_tc>1.66 & tc<5.56")
lo<-subset(heat, loc_tar_CDS=="membrane_tc<1.66")
hi$gene_id<-factor(hi$gene_id, levels=unique(hi$gene_id[order(hi$tc_CDS_norm, decreasing = T)]), ordered = T)
mid$gene_id<-factor(mid$gene_id, levels=unique(mid$gene_id[order(mid$tc_CDS_norm, decreasing = T)]), ordered = T)
lo$gene_id<-factor(lo$gene_id, levels=unique(lo$gene_id[order(lo$tc_CDS_norm, decreasing = T)]), ordered = T)

ggplot(hi, aes(gene_id, variable))+geom_tile(aes(fill=log2(value)))+
  scale_fill_continuous(limits=c(0,20),na.value = 'black')+theme(axis.title.x=element_blank(),
                                                                 axis.text.x=element_blank(),
                                                                 axis.ticks.x=element_blank())
ggplot(mid, aes(gene_id, variable))+geom_tile(aes(fill=log2(value)))+
  scale_fill_continuous(limits=c(0,20),na.value = 'black')+theme(axis.title.x=element_blank(),
                                                                 axis.text.x=element_blank(),
                                                                 axis.ticks.x=element_blank())
ggplot(lo, aes(gene_id, variable))+geom_tile(aes(fill=log2(value)))+
  scale_fill_continuous(limits=c(0,20),na.value = 'black')+theme(axis.title.x=element_blank(),
                                                                 axis.text.x=element_blank(),
                                                                 axis.ticks.x=element_blank())

hisc<-subset(ave, loc_tar_CDS=="membrane_tc>5.56 & tc<65.26")
hisc$gene_id<-factor(hisc$gene_id, levels=unique(hisc$gene_id[order(hisc$tc_CDS_norm, decreasing = T)]), ordered = T)
ggplot(hisc, aes(gene_id, dum))+geom_tile(aes(fill=log2(tc_CDS_norm)))+scale_fill_continuous(limits=c(-4,8),na.value = 'black')+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

midsc<-subset(ave, loc_tar_CDS=="membrane_tc>1.66 & tc<5.56")
midsc$gene_id<-factor(midsc$gene_id, levels=unique(midsc$gene_id[order(midsc$tc_CDS_norm, decreasing = T)]), ordered = T)
ggplot(midsc, aes(gene_id, dum))+geom_tile(aes(fill=log2(tc_CDS_norm)))+scale_fill_continuous(limits=c(-4,8),na.value = 'black')+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

losc<-subset(ave, loc_tar_CDS=="membrane_tc<1.66")
losc$gene_id<-factor(losc$gene_id, levels=unique(losc$gene_id[order(losc$tc_CDS_norm, decreasing = T)]), ordered = T)
ggplot(losc, aes(gene_id, dum))+geom_tile(aes(fill=log2(tc_CDS_norm)))+scale_fill_continuous(limits=c(-4,8),na.value = 'black')+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())



###compare between random trimers and codons

ave<-merge(ave, rsem, by="gene_id")
relseq<-tcCds[ave$transcript_id]

bcods <- (Biostrings::trinucleotideFrequency(relseq, step = 3, as.prob = F, with.labels = TRUE))
rands <- (Biostrings::trinucleotideFrequency(relseq, step = 1, as.prob = F, with.labels = TRUE))

bcods<-data.frame(transcript_id=names(relseq), bcods)
rands<-data.frame(transcript_id=names(relseq), rands)

bcods<-merge(bcods, ave, by="transcript_id")
rands<-merge(rands, ave, by="transcript_id")

bcods<-aggregate(cbind(AAA,AAC,AAG,AAT,ACA,ACC,ACG,ACT,AGA,AGC,AGG,AGT,ATA,ATC,ATG,ATT,CAA,CAC,CAG,CAT,CCA,CCC,CCG,CCT,CGA,CGC,CGG,CGT,CTA,CTC,CTG,CTT,GAA,GAC,GAG,GAT,GCA,GCC,GCG,GCT,GGA,GGC,GGG,GGT,GTA,GTC,GTG,GTT,TAA,TAC,TAG,TAT,TCA,TCC,TCG,TCT,TGA,TGC,TGG,TGT,TTA,TTC,TTG,TTT)~loc_tar_CDS, data=bcods, sum)
rands<-aggregate(cbind(AAA,AAC,AAG,AAT,ACA,ACC,ACG,ACT,AGA,AGC,AGG,AGT,ATA,ATC,ATG,ATT,CAA,CAC,CAG,CAT,CCA,CCC,CCG,CCT,CGA,CGC,CGG,CGT,CTA,CTC,CTG,CTT,GAA,GAC,GAG,GAT,GCA,GCC,GCG,GCT,GGA,GGC,GGG,GGT,GTA,GTC,GTG,GTT,TAA,TAC,TAG,TAT,TCA,TCC,TCG,TCT,TGA,TGC,TGG,TGT,TTA,TTC,TTG,TTT)~loc_tar_CDS, data=rands, sum)

bcodsSum<-rowSums(bcods[,2:ncol(bcods)])
randsSum<-rowSums(rands[,2:ncol(rands)])

bcods<-cbind(loc_tar_CDS=bcods$loc_tar_CDS,bcods[,2:ncol(bcods)]/bcodsSum)
rands<-cbind(loc_tar_CDS=rands$loc_tar_CDS,rands[,2:ncol(rands)]/randsSum)

colnames(bcods)[2:ncol(bcods)]<-gsub("T","U", colnames(bcods)[2:ncol(bcods)])
colnames(rands)[2:ncol(rands)]<-gsub("T","U", colnames(rands)[2:ncol(rands)])

uid<-merge(cod, data.frame(codon=colnames(bcods)[2:ncol(bcods)], seq=seq(1,64)), by="codon")
uid$codon_aa<-paste0(uid$codon,";", uid$aa2)

colnames(bcods)[2:ncol(bcods)]<-uid$codon_aa
colnames(rands)[2:ncol(rands)]<-uid$codon_aa


mel<-cbind(melt(bcods, id.vars = "loc_tar_CDS", variable.name = "bcods", value.name = "freq_bcods"),
           melt(rands, id.vars = "loc_tar_CDS", variable.name = "rands", value.name = "freq_rands"))

ggplot(mel, aes(freq_bcods, freq_rands))+geom_text(aes(label=bcods), size=3)+geom_abline(slope=1)+facet_grid(~loc_tar_CDS)

mel$diff<-mel$freq_bcods-mel$freq_rands
mel$loc_tar_CDS<-factor(mel$loc_tar_CDS, levels=
                          c("membrane_tc>10.82 & tc<181.46", "membrane_tc>2.81 & tc<10.82", "membrane_tc<2.81", "nontarget_membrane",
                            "cytosolic_tc>0.83 & tc<166.4","cytosolic_tc>0.3 & tc<0.83","cytosolic_tc<0.3", "nontarget_cytosolic"))

ggplot(mel, aes(loc_tar_CDS, bcods))+geom_tile(aes(fill=diff))+theme(axis.text.x = element_text(angle = 90)) 

ggplot(mel, aes(loc_tar_CDS, bcods))+geom_tile(aes(fill=freq_bcods))+theme(axis.text.x = element_text(angle = 90))+
  scale_fill_continuous(limits=c(0,0.06),na.value = 'black')
ggplot(mel, aes(loc_tar_CDS, bcods))+geom_tile(aes(fill=freq_rands))+theme(axis.text.x = element_text(angle = 90))+
  scale_fill_continuous(limits=c(0,0.06),na.value = 'black')

#only look at U codons and no Stop codons
mel<-subset(mel, grepl("U",mel$bcods))
mel<-subset(mel, !grepl("Stop",mel$bcods))
ggplot(mel, aes(loc_tar_CDS, bcods))+geom_tile(aes(fill=diff))+theme(axis.text.x = element_text(angle = 90)) 

ggplot(mel, aes(loc_tar_CDS, bcods))+geom_tile(aes(fill=freq_bcods))+theme(axis.text.x = element_text(angle = 90))+
  scale_fill_continuous(limits=c(0,0.06),na.value = 'black')
ggplot(mel, aes(loc_tar_CDS, bcods))+geom_tile(aes(fill=freq_rands))+theme(axis.text.x = element_text(angle = 90))+
  scale_fill_continuous(limits=c(0,0.06),na.value = 'black')



#correlate codon crosslinks and fold changes
sel<-subset(ave, localization_cat=="membrane")
res<-apply(sel[2:35], 2, function(x) cor(x, sel$log2FoldChange.ribo.rna.KO.WT,use="pairwise.complete.obs", method="spearman") )
res<-data.frame(codon=names(res),spearman=res)
ggplot(res, aes(codon, spearman))+geom_bar(stat="identity")+coord_flip()

res<-apply(sel[2:35], 2, function(x) cor(x, sel$log2FoldChange.mem.cyt.293,use="pairwise.complete.obs", method="spearman") )
res<-data.frame(codon=names(res),spearman=res)
ggplot(res, aes(codon, spearman))+geom_bar(stat="identity")+coord_flip()

res<-apply(sel[2:35], 2, function(x) cor(x, sel$log2FoldChange.mem.cyt.KO.293,use="pairwise.complete.obs", method="spearman") )
res<-data.frame(codon=names(res),spearman=res)
ggplot(res, aes(codon, spearman))+geom_bar(stat="identity")+coord_flip()

res<-apply(sel[2:35], 2, function(x) cor(x, sel$mean_te_293,use="pairwise.complete.obs", method="spearman") )
res<-data.frame(codon=names(res),spearman=res)
ggplot(res, aes(codon, spearman))+geom_bar(stat="identity")+coord_flip()




###6mers
bcods <- (Biostrings::oligonucleotideFrequency(relseq, width=6, step = 3, as.prob = F, with.labels = TRUE))
rands <- (Biostrings::oligonucleotideFrequency(relseq, width=6, step = 1, as.prob = F, with.labels = TRUE))

bcods<-data.frame(transcript_id=names(relseq), bcods)
rands<-data.frame(transcript_id=names(relseq), rands)

bcods<-merge(bcods, ave, by="transcript_id")
rands<-merge(rands, ave, by="transcript_id")

bcods<-apply(bcods[2:(ncol(bcods)-50)], 2, function(x) aggregate(x~bcods$loc_tar_CDS, data=bcods, sum))
rands<-apply(rands[2:(ncol(rands)-50)], 2, function(x) aggregate(x~rands$loc_tar_CDS, data=rands, sum))
names_bcods<-bcods[[1]][1]
names_rands<-rands[[1]][1]

bcods<-bind_cols(lapply(bcods, function(x) x[,2]))
rands<-bind_cols(lapply(rands, function(x) x[,2]))

bcodsSum<-rowSums(bcods)
randsSum<-rowSums(rands)

bcods<-cbind(loc_tar_CDS=names_bcods,bcods/bcodsSum)
rands<-cbind(loc_tar_CDS=names_rands,rands/randsSum)

colnames(bcods)[2:ncol(bcods)]<-gsub("T","U", colnames(bcods)[2:ncol(bcods)])
colnames(rands)[2:ncol(rands)]<-gsub("T","U", colnames(rands)[2:ncol(rands)])




mel<-cbind(melt(bcods, id.vars = "bcods$loc_tar_CDS", variable.name = "bcods", value.name = "freq_bcods"),
           melt(rands, id.vars = "rands$loc_tar_CDS", variable.name = "rands", value.name = "freq_rands"))

ggplot(mel, aes(freq_bcods, freq_rands))+geom_text(aes(label=bcods), size=2)+geom_abline(slope=1)+facet_wrap(~`bcods$loc_tar_CDS`)

mel$diff<-mel$freq_bcods-mel$freq_rands
mel$loc_tar_CDS<-factor(mel$`bcods$loc_tar_CDS`, levels=
                          c("membrane_tc>10.82 & tc<181.46", "membrane_tc>2.81 & tc<10.82", "membrane_tc<2.81", "nontarget_membrane",
                            "cytosolic_tc>0.83 & tc<166.4","cytosolic_tc>0.3 & tc<0.83","cytosolic_tc<0.3", "nontarget_cytosolic"))

ggplot(mel, aes(loc_tar_CDS, bcods))+geom_tile(aes(fill=diff))+theme(axis.text.x = element_text(angle = 90)) 

ggplot(mel, aes(loc_tar_CDS, bcods))+geom_tile(aes(fill=freq_bcods))+theme(axis.text.x = element_text(angle = 90))+
  scale_fill_continuous(na.value = 'black')
ggplot(mel, aes(loc_tar_CDS, bcods))+geom_tile(aes(fill=freq_rands))+theme(axis.text.x = element_text(angle = 90))+
  scale_fill_continuous(na.value = 'black')


bcodsmem<-melt(bcods[grepl("membrane", bcods$`bcods$loc_tar_CDS`),], id.vars = "bcods$loc_tar_CDS", variable.name = "bcods", value.name = "freq_bcods")
colnames(bcodsmem)<-paste0(colnames(bcodsmem),"_membrane")
bcodscyt<-melt(bcods[grepl("cytosol", bcods$`bcods$loc_tar_CDS`),], id.vars = "bcods$loc_tar_CDS", variable.name = "bcods", value.name = "freq_bcods")
colnames(bcodscyt)<-paste0(colnames(bcodscyt),"_cytosolic")
bcodsmem<-aggregate(freq_bcods_membrane~bcods_membrane, data=bcodsmem,mean)
bcodscyt<-aggregate(freq_bcods_cytosolic~bcods_cytosolic, data=bcodscyt,mean)

loccods<-cbind(bcodsmem,bcodscyt)

ggplot(loccods, aes(freq_bcods_membrane, freq_bcods_cytosolic))+geom_text(aes(label=bcods_membrane), size=2)+geom_abline(slope=1)

randsmem<-melt(rands[grepl("membrane", rands$`rands$loc_tar_CDS`),], id.vars = "rands$loc_tar_CDS", variable.name = "rands", value.name = "freq_rands")
colnames(randsmem)<-paste0(colnames(randsmem),"_membrane")
randscyt<-melt(rands[grepl("cytosol", rands$`rands$loc_tar_CDS`),], id.vars = "rands$loc_tar_CDS", variable.name = "rands", value.name = "freq_rands")
colnames(randscyt)<-paste0(colnames(randscyt),"_cytosolic")
randsmem<-aggregate(freq_rands_membrane~rands_membrane, data=randsmem,mean)
randscyt<-aggregate(freq_rands_cytosolic~rands_cytosolic, data=randscyt,mean)

locrands<-cbind(randsmem,randscyt)
ggplot(locrands, aes(freq_rands_membrane, freq_rands_cytosolic))+geom_text(aes(label=rands_membrane), size=2)+geom_abline(slope=1)


###4mers
bcods <- (Biostrings::oligonucleotideFrequency(relseq, width=4, step = 3, as.prob = F, with.labels = TRUE))
rands <- (Biostrings::oligonucleotideFrequency(relseq, width=4, step = 1, as.prob = F, with.labels = TRUE))

bcods<-data.frame(transcript_id=names(relseq), bcods)
rands<-data.frame(transcript_id=names(relseq), rands)

bcods<-merge(bcods, ave, by="transcript_id")
rands<-merge(rands, ave, by="transcript_id")

bcods<-apply(bcods[2:(ncol(bcods)-50)], 2, function(x) aggregate(x~bcods$loc_tar_CDS, data=bcods, sum))
rands<-apply(rands[2:(ncol(rands)-50)], 2, function(x) aggregate(x~rands$loc_tar_CDS, data=rands, sum))
names_bcods<-bcods[[1]][1]
names_rands<-rands[[1]][1]

bcods<-bind_cols(lapply(bcods, function(x) x[,2]))
rands<-bind_cols(lapply(rands, function(x) x[,2]))

bcodsSum<-rowSums(bcods)
randsSum<-rowSums(rands)

bcods<-cbind(loc_tar_CDS=names_bcods,bcods/bcodsSum)
rands<-cbind(loc_tar_CDS=names_rands,rands/randsSum)

colnames(bcods)[2:ncol(bcods)]<-gsub("T","U", colnames(bcods)[2:ncol(bcods)])
colnames(rands)[2:ncol(rands)]<-gsub("T","U", colnames(rands)[2:ncol(rands)])




mel<-cbind(melt(bcods, id.vars = "bcods$loc_tar_CDS", variable.name = "bcods", value.name = "freq_bcods"),
           melt(rands, id.vars = "rands$loc_tar_CDS", variable.name = "rands", value.name = "freq_rands"))

ggplot(mel, aes(freq_bcods, freq_rands))+geom_text(aes(label=bcods), size=2)+geom_abline(slope=1)+facet_wrap(~`bcods$loc_tar_CDS`)

mel$diff<-mel$freq_bcods-mel$freq_rands
mel$loc_tar_CDS<-factor(mel$`bcods$loc_tar_CDS`, levels=
                          c("membrane_tc>5.56 & tc<65.26", "membrane_tc>1.66 & tc<5.56", "membrane_tc<1.66", "nontarget_membrane",
                            "cytosolic_tc>0.47 & tc<27.84","cytosolic_tc>0.18 & tc<0.47","cytosolic_tc<0.18", "nontarget_cytosolic"))

ggplot(mel, aes(loc_tar_CDS, bcods))+geom_tile(aes(fill=diff))+theme(axis.text.x = element_text(angle = 90)) 

ggplot(mel, aes(loc_tar_CDS, bcods))+geom_tile(aes(fill=freq_bcods))+theme(axis.text.x = element_text(angle = 90))+
  scale_fill_continuous(na.value = 'black')
ggplot(mel, aes(loc_tar_CDS, bcods))+geom_tile(aes(fill=freq_rands))+theme(axis.text.x = element_text(angle = 90))+
  scale_fill_continuous(na.value = 'black')


bcodsmem<-melt(bcods[grepl("membrane", bcods$`bcods$loc_tar_CDS`),], id.vars = "bcods$loc_tar_CDS", variable.name = "bcods", value.name = "freq_bcods")
colnames(bcodsmem)<-paste0(colnames(bcodsmem),"_membrane")
bcodscyt<-melt(bcods[grepl("cytosol", bcods$`bcods$loc_tar_CDS`),], id.vars = "bcods$loc_tar_CDS", variable.name = "bcods", value.name = "freq_bcods")
colnames(bcodscyt)<-paste0(colnames(bcodscyt),"_cytosolic")
bcodsmem<-aggregate(freq_bcods_membrane~bcods_membrane, data=bcodsmem,mean)
bcodscyt<-aggregate(freq_bcods_cytosolic~bcods_cytosolic, data=bcodscyt,mean)

loccods<-cbind(bcodsmem,bcodscyt)

ggplot(loccods, aes(freq_bcods_membrane, freq_bcods_cytosolic))+geom_text(aes(label=bcods_membrane), size=2)+geom_abline(slope=1)

randsmem<-melt(rands[grepl("membrane", rands$`rands$loc_tar_CDS`),], id.vars = "rands$loc_tar_CDS", variable.name = "rands", value.name = "freq_rands")
colnames(randsmem)<-paste0(colnames(randsmem),"_membrane")
randscyt<-melt(rands[grepl("cytosol", rands$`rands$loc_tar_CDS`),], id.vars = "rands$loc_tar_CDS", variable.name = "rands", value.name = "freq_rands")
colnames(randscyt)<-paste0(colnames(randscyt),"_cytosolic")
randsmem<-aggregate(freq_rands_membrane~rands_membrane, data=randsmem,mean)
randscyt<-aggregate(freq_rands_cytosolic~rands_cytosolic, data=randscyt,mean)

locrands<-cbind(randsmem,randscyt)
ggplot(locrands, aes(freq_rands_membrane, freq_rands_cytosolic))+geom_text(aes(label=rands_membrane), size=3)+geom_abline(slope=1, lty=2, colour="grey")


##compare 3'UTR 4mers
# setwd("~/Google Drive/hdlbp/")
setwd("E:/Google Drive/hdlbp/")

mas<-read.delim("hdlbp_master_table_with_classes_uniq_tsig.txt", header=T)
inf<-subset(mas, select=c("gene_id","Symbol","transcript" ,
                          "tpm_cutoff","tc_CDS_norm","localization_cat",
                          "mean_te_293","loc_tar_CDS","tc_CDS_norm_cat",
                          "tc_transcript_norm",
                          "log2FoldChange.mem.cyt.293",
                          "log2FoldChange.mem.cyt.KO.293",
                          "log2FoldChange.ribo.rna.KO.WT"))
# tcUtr3<-seqUtr3[unique(dat$transcript_id)]
# 
# relseq<-tcUtr3[ave$transcript_id]

# rands <- (Biostrings::oligonucleotideFrequency(relseq, width=4, step = 1, as.prob = F, with.labels = TRUE))
seqUtr3stop<- Biostrings::subseq(sub, start = ham$cds_stop+1-30, end = ham$l_tr)#here include last 30 nt from CDS

rands <- (Biostrings::oligonucleotideFrequency(seqUtr3stop, width=4, step = 1, as.prob = F, with.labels = TRUE))


# rands<-data.frame(transcript_id=names(relseq), rands)
rands<-data.frame(transcript_id=names(seqUtr3stop), rands)

rands<-merge(rands, inf, by.x="transcript_id", by.y="transcript")


rands<-apply(rands[,2:257], 2, function(x) aggregate(x~rands$localization_cat, data=rands, sum))
names_rands<-rands[[1]][1]

rands<-bind_cols(lapply(rands, function(x) x[,2]))

randsSum<-rowSums(rands)


rands<-cbind(localization_cat=names_rands,rands/randsSum)

colnames(rands)[2:ncol(rands)]<-gsub("T","U", colnames(rands)[2:ncol(rands)])




mel<-melt(rands, id.vars = "rands$localization_cat", variable.name = "rands", value.name = "freq_rands")


mel$localization_cat<-factor(mel$`rands$localization_cat`, levels=
                          c("membrane", "cytosolic"))
utr3rands<-mel
utr3rands$id<-paste0(utr3rands$rands,"_",utr3rands$localization_cat)

ggplot(mel, aes(localization_cat, rands))+geom_tile(aes(fill=freq_rands))+theme(axis.text.x = element_text(angle = 90))+
  scale_fill_continuous(na.value = 'black')



randsmem<-melt(rands[grepl("membrane", rands$`rands$localization_cat`),], id.vars = "rands$localization_cat", variable.name = "rands", value.name = "freq_rands")
colnames(randsmem)<-paste0(colnames(randsmem),"_membrane")
randscyt<-melt(rands[grepl("cytosol", rands$`rands$localization_cat`),], id.vars = "rands$localization_cat", variable.name = "rands", value.name = "freq_rands")
colnames(randscyt)<-paste0(colnames(randscyt),"_cytosolic")
randsmem<-aggregate(freq_rands_membrane~rands_membrane, data=randsmem,mean)
randscyt<-aggregate(freq_rands_cytosolic~rands_cytosolic, data=randscyt,mean)

locrands<-cbind(randsmem,randscyt)
ggplot(locrands, aes(freq_rands_membrane, freq_rands_cytosolic))+geom_text(aes(label=rands_membrane), size=3)+geom_abline(slope=1, lty=2, colour="grey")


##compare CDS 4mers
# tcCds<-seqCds[unique(dat$transcript_id)]
# 
# relseq<-tcCds[ave$transcript_id]
# rands <- (Biostrings::oligonucleotideFrequency(relseq, width=4, step = 1, as.prob = F, with.labels = TRUE))

rands <- (Biostrings::oligonucleotideFrequency(seqCds, width=4, step = 1, as.prob = F, with.labels = TRUE))

# rands<-data.frame(transcript_id=names(relseq), rands)

rands<-data.frame(transcript_id=names(seqCds), rands)

rands<-merge(rands, inf, by.x="transcript_id", by.y="transcript")


rands<-apply(rands[,2:257], 2, function(x) aggregate(x~rands$localization_cat, data=rands, sum))
names_rands<-rands[[1]][1]

rands<-bind_cols(lapply(rands, function(x) x[,2]))

randsSum<-rowSums(rands)


rands<-cbind(localization_cat=names_rands,rands/randsSum)

colnames(rands)[2:ncol(rands)]<-gsub("T","U", colnames(rands)[2:ncol(rands)])




mel<-melt(rands, id.vars = "rands$localization_cat", variable.name = "rands", value.name = "freq_rands")


mel$localization_cat<-factor(mel$`rands$localization_cat`, levels=
                               c("membrane", "cytosolic"))
cdsrands<-mel
cdsrands$id<-paste0(cdsrands$rands,"_",cdsrands$localization_cat)

ggplot(mel, aes(localization_cat, rands))+geom_tile(aes(fill=freq_rands))+theme(axis.text.x = element_text(angle = 90))+
  scale_fill_continuous(na.value = 'black')



randsmem<-melt(rands[grepl("membrane", rands$`rands$localization_cat`),], id.vars = "rands$localization_cat", variable.name = "rands", value.name = "freq_rands")
colnames(randsmem)<-paste0(colnames(randsmem),"_membrane")
randscyt<-melt(rands[grepl("cytosol", rands$`rands$localization_cat`),], id.vars = "rands$localization_cat", variable.name = "rands", value.name = "freq_rands")
colnames(randscyt)<-paste0(colnames(randscyt),"_cytosolic")
randsmem<-aggregate(freq_rands_membrane~rands_membrane, data=randsmem,mean)
randscyt<-aggregate(freq_rands_cytosolic~rands_cytosolic, data=randscyt,mean)

locrands<-cbind(randsmem,randscyt)
ggplot(locrands, aes(freq_rands_membrane, freq_rands_cytosolic))+geom_text(aes(label=rands_membrane), size=3)+geom_abline(slope=1, lty=2, colour="grey")
locrands$rands_membrane<-as.character(locrands$rands_membrane)
locrands$labels<-ifelse(locrands$rands_membrane=="CCUG" | locrands$rands_membrane=="CUUC" | locrands$rands_membrane=="CCUC" | locrands$rands_membrane=="CUCU" | locrands$rands_membrane=="AGAA" | locrands$rands_membrane=="UGGA" ,
                       locrands$rands_membrane, NA)
ggplot(locrands, aes(freq_rands_membrane, freq_rands_cytosolic, colour=labels))+geom_text(aes(label=rands_membrane), size=3)+geom_abline(slope=1, lty=2, colour="grey")

###compare between CDS and 3UTR

colnames(utr3rands)<-paste0("utr3_",colnames(utr3rands))
colnames(cdsrands)<-paste0("cds_",colnames(cdsrands))
utr3cds<-merge(utr3rands, cdsrands, by.x="utr3_id", by.y="cds_id")
utr3cds$utr3_rands<-as.character(utr3rands$utr3_rands)
utr3cds$labels<-ifelse(utr3cds$utr3_rands=="CCUG" | utr3cds$utr3_rands=="CUUC" | utr3cds$utr3_rands=="CCUC" | utr3cds$utr3_rands=="CUCU" | utr3cds$utr3_rands=="AGAA" | utr3cds$utr3_rands=="UGGA" ,
                        utr3cds$utr3_rands, NA)

ggplot(utr3cds, aes(cds_freq_rands, utr3_freq_rands))+geom_text(aes(label=utr3_rands), size=3)+geom_abline(slope=1, lty=2, colour="grey")+facet_wrap(~utr3_localization_cat)
ggplot(utr3cds, aes(cds_freq_rands, utr3_freq_rands, colour=labels))+geom_text(aes(label=utr3_rands), size=3)+geom_abline(slope=1, lty=2, colour="grey")+facet_wrap(~utr3_localization_cat)
ggplot(utr3cds, aes(cds_freq_rands, utr3_freq_rands, colour=labels))+geom_text(aes(label=utr3_rands), size=3)+geom_abline(slope=1, lty=2, colour="grey")
utr3_cyto<-subset(utr3rands, utr3_localization_cat=="cytosolic")
cds_mem<-subset(cdsrands, cds_localization_cat=="membrane")
plot<-merge(utr3_cyto, cds_mem, by.x="utr3_rands", by.y="cds_rands")

plot$utr3_rands<-as.character(plot$utr3_rands)
plot$labels<-ifelse(plot$utr3_rands=="CCUG" | plot$utr3_rands=="CUUC" | plot$utr3_rands=="CCUC" | plot$utr3_rands=="CUCU" | plot$utr3_rands=="AGAA" | plot$utr3_rands=="UGGA" ,
                       plot$utr3_rands, NA)
ggplot(plot, aes(cds_freq_rands, utr3_freq_rands))+geom_text(aes(label=utr3_rands), size=3)+geom_abline(slope=1, lty=2, colour="grey")
ggplot(plot, aes(cds_freq_rands, utr3_freq_rands, colour=labels))+geom_text(aes(label=utr3_rands), size=3)+geom_abline(slope=1, lty=2, colour="grey")

utr3_mem<-subset(utr3rands, utr3_localization_cat=="membrane")
cds_cyto<-subset(cdsrands, cds_localization_cat=="cytosolic")
plot<-merge(utr3_mem, cds_cyto, by.x="utr3_rands", by.y="cds_rands")
plot$utr3_rands<-as.character(plot$utr3_rands)
plot$labels<-ifelse(plot$utr3_rands=="CCUG" | plot$utr3_rands=="CUUC" | plot$utr3_rands=="CCUC" | plot$utr3_rands=="CUCU" | plot$utr3_rands=="AGAA" | plot$utr3_rands=="UGGA" ,
                    plot$utr3_rands, NA)

ggplot(plot, aes(cds_freq_rands, utr3_freq_rands, colour=labels))+geom_text(aes(label=utr3_rands), size=3)+geom_abline(slope=1, lty=2, colour="grey")

utr3_cyto<-subset(utr3rands, utr3_localization_cat=="cytosolic")
cds_cyto<-subset(cdsrands, cds_localization_cat=="cytosolic")
plot<-merge(utr3_cyto, cds_cyto, by.x="utr3_rands", by.y="cds_rands")
plot$utr3_rands<-as.character(plot$utr3_rands)
plot$labels<-ifelse(plot$utr3_rands=="CCUG" | plot$utr3_rands=="CUUC" | plot$utr3_rands=="CCUC" | plot$utr3_rands=="CUCU" | plot$utr3_rands=="AGAA" | plot$utr3_rands=="UGGA" ,
                    plot$utr3_rands, NA)

ggplot(plot, aes(cds_freq_rands, utr3_freq_rands, colour=labels))+geom_text(aes(label=utr3_rands), size=3)+geom_abline(slope=1, lty=2, colour="grey")


##cds vs utr3 with no localization
rands <- (Biostrings::oligonucleotideFrequency(seqUtr3, width=4, step = 1, as.prob = F, with.labels = TRUE))
# rands<-data.frame(transcript_id=names(relseq), rands)
rands<-data.frame(transcript_id=names(seqUtr3), rands)

rands<-apply(rands[2:(ncol(rands))], 2, function(x) sum(x))
names_rands<-names(rands)

rands<-data.frame(mer=names(rands), rands)
rands$utr3_freq<-rands$rands/sum(rands$rands)
utr3rands<-rands

rands <- (Biostrings::oligonucleotideFrequency(seqCds, width=4, step = 1, as.prob = F, with.labels = TRUE))
# rands<-data.frame(transcript_id=names(relseq), rands)
rands<-data.frame(transcript_id=names(seqCds), rands)

rands<-apply(rands[2:(ncol(rands))], 2, function(x) sum(x))
names_rands<-names(rands)

rands<-data.frame(mer=names(rands), rands)
rands$cds_freq<-rands$rands/sum(rands$rands)
cdsrands<-rands

utr3cds<-merge(utr3rands, cdsrands, by="mer")
utr3cds$mer<-gsub("T","U", utr3cds$mer)

ggplot(utr3cds, aes(cds_freq, utr3_freq))+geom_text(aes(label=mer), size=3)+geom_abline(slope=1, lty=2, colour="grey")
subset(utr3cds, mer=="CUUC")
utr3cds$labels<-ifelse(utr3cds$mer=="CCUG" | utr3cds$mer=="CUUC" | utr3cds$mer=="CCUC" | utr3cds$mer=="CUCU" | utr3cds$mer=="AGAA" | utr3cds$mer=="UGGA" ,
                       utr3cds$mer, NA)
ggplot(utr3cds, aes(cds_freq, utr3_freq, colour=labels))+geom_text(aes(label=mer), size=3)+geom_abline(slope=1, lty=2, colour="grey")

##which 4mers are most bound - TC signal in 4mers


tcCds<-seqCds[unique(dat$transcript_id)]
dat$tc_from_start_4mer<-ifelse(dat$tc_from_start+3>dat$l_cds , NA, dat$tc_from_start)
dat<-subset(dat, !is.na(tc_from_start_4mer))
#here set position with tc_from_start_4mer
dat<-subset(dat, tc_from_start_4mer-3>1)
dat$tc_4mer<-as.character(Biostrings::subseq(tcCds[as.character(dat$transcript_id)],
                                              start=dat$tc_from_start_4mer-3,end=dat$tc_from_start_4mer))


ag1<-dcast(dat[,c("gene_id","tc_4mer","norm_tc_num1")], gene_id~tc_4mer,sum,value.var="norm_tc_num1")
ag2<-dcast(dat[,c("gene_id","tc_4mer","norm_tc_num2")], gene_id~tc_4mer,sum,value.var="norm_tc_num2")

ag1$tc_4mer_tot<-rowSums(ag1[2:ncol(ag1)],na.rm=T)
ag2$tc_4mer_tot<-rowSums(ag2[2:ncol(ag2)],na.rm=T)
ag1$tc_4mer_tot<-ifelse(ag1$tc_4mer_tot==0,NA,ag1$tc_4mer_tot)
ag2$tc_4mer_tot<-ifelse(ag2$tc_4mer_tot==0,NA,ag2$tc_4mer_tot)

seq_freq <- (Biostrings::oligonucleotideFrequency(tcCds, width=4, step = 1, as.prob = F, with.labels = TRUE))
rownames(seq_freq)<-names(tcCds)
seq_freq<-subset(seq_freq, select=colnames(ag1)[2:(ncol(ag1)-1)])
seq_freq[seq_freq==0]<-NA
# fin1<-ag1[2:(ncol(ag1)-1)]/ag1$tc_codon_tot*1e6
fin1<-ag1[2:(ncol(ag1)-1)]/ag1$tc_4mer_tot/seq_freq*1e6
fin1<-cbind(gene_id=ag1$gene_id, fin1, tc_4mer_tot=ag1$tc_4mer_tot)
# fin1<-subset(fin1, select=colnames(fin1)[!grepl("TGA",colnames(fin1)) & !grepl("TAA",colnames(fin1)) & !grepl("TAG",colnames(fin1))])
heat1<-melt(fin1, measure.vars = colnames(fin1)[2:65] , id.vars="gene_id")

# fin2<-ag2[2:(ncol(ag2)-1)]/ag2$tc_codon_tot*1e6
fin2<-ag2[2:(ncol(ag2)-1)]/ag2$tc_4mer_tot/seq_freq*1e6
fin2<-cbind(gene_id=ag2$gene_id, fin2, tc_4mer_tot=ag2$tc_4mer_tot)
# fin2<-subset(fin2, select=colnames(fin2)[!grepl("TGA",colnames(fin2)) & !grepl("TAA",colnames(fin2)) & !grepl("TAG",colnames(fin2))])
heat2<-melt(fin2, measure.vars = colnames(fin2)[2:65] , id.vars="gene_id")

# ggplot(heat1, aes(gene_id, variable))+geom_tile(aes(fill=log2(value)))+scale_fill_continuous(na.value = 'black')
# ggplot(heat2, aes(gene_id, variable))+geom_tile(aes(fill=log2(value)))+scale_fill_continuous(na.value = 'black')

colnames(fin1)[2:(ncol(fin1))]<-paste0("rep1_",colnames(fin1)[2:(ncol(fin1))])
colnames(fin2)[2:(ncol(fin2))]<-paste0("rep2_",colnames(fin2)[2:(ncol(fin2))])
fin<-merge(fin1, fin2, by="gene_id")

# setwd("~/Google Drive/hdlbp/")
setwd("E:/Google Drive/hdlbp/")
# mas<-read.delim("hdlbp_master_table_with_classes.txt", header=T)
mas<-read.delim("hdlbp_master_table_with_classes_uniq.txt", header=T)
inf<-subset(mas, select=c("gene_id","Symbol",
                          "tpm_cutoff","tc_CDS_norm","localization_cat",
                          "mean_te_293","loc_tar_CDS","tc_CDS_norm_cat",
                          "tc_transcript_norm",
                          "log2FoldChange.mem.cyt.293",
                          "log2FoldChange.mem.cyt.KO.293",
                          "log2FoldChange.ribo.rna.KO.WT"))

fin<-merge(fin, inf, by="gene_id")

ggplot(fin, aes(log2(rep1_tc_4mer_tot), log2(rep2_tc_4mer_tot)))+geom_point()

ave<-(fin[,2:66]+fin[,67:131])/2
colnames(ave)<-gsub("rep1_","",colnames(ave))
ave<-cbind(gene_id=fin$gene_id, ave, fin[,132:ncol(fin)])
colnames(ave)[2:65]<-gsub("T","U",colnames(ave)[2:65])

ave<-subset(ave, tc_4mer_tot>=10) ## thresholdat least 10 T-C per CDS

heat<-melt(ave, measure.vars = colnames(ave)[c(2:65)] , 
           id.vars=c("gene_id","tpm_cutoff","tc_CDS_norm",
                     "localization_cat","mean_te_293","loc_tar_CDS","tc_CDS_norm_cat",
                     "log2FoldChange.mem.cyt.293",
                     "log2FoldChange.mem.cyt.KO.293",
                     "log2FoldChange.ribo.rna.KO.WT",
                     "tc_transcript_norm"))


heat$gene_id<-factor(heat$gene_id, levels=unique(heat$gene_id[order(heat$log2FoldChange.mem.cyt.293, decreasing = T)]), ordered = T)
ggplot(heat, aes(gene_id, variable))+geom_tile(aes(fill=log2(value)))+
  scale_fill_continuous(limits=c(0,20),na.value = 'black')+theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
                                                                 axis.text.x=element_blank(),
                                                                 axis.ticks.x=element_blank())

ggplot(subset(heat, !is.na(localization_cat)), aes(gene_id, variable))+geom_tile(aes(fill=log2(value)))+
  scale_fill_continuous(limits=c(0,20),na.value = 'black')+theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
                                                                 axis.text.x=element_blank(),
                                                                 axis.ticks.x=element_blank())+facet_wrap(~localization_cat, scales="free_x")
ggplot(heat, aes(variable, value))+geom_boxplot()+coord_flip()

###which 4mers are most bound - TC signal in adjacent 4mers


tcCds<-seqCds[unique(dat$transcript_id)]
dat$tc_from_start_4mer<-ifelse(dat$tc_from_start+4>dat$l_cds , NA, dat$tc_from_start)
dat<-subset(dat, !is.na(tc_from_start_4mer))
#here set position with tc_from_start_4mer
# dat<-subset(dat, tc_from_start_4mer-4>1)
dat$tc_4mer<-as.character(Biostrings::subseq(tcCds[as.character(dat$transcript_id)],
                                             start=dat$tc_from_start_4mer+1,end=dat$tc_from_start_4mer+4))


ag1<-dcast(dat[,c("gene_id","tc_4mer","norm_tc_num1")], gene_id~tc_4mer,sum,value.var="norm_tc_num1")
ag2<-dcast(dat[,c("gene_id","tc_4mer","norm_tc_num2")], gene_id~tc_4mer,sum,value.var="norm_tc_num2")

ag1$tc_4mer_tot<-rowSums(ag1[2:ncol(ag1)],na.rm=T)
ag2$tc_4mer_tot<-rowSums(ag2[2:ncol(ag2)],na.rm=T)
ag1$tc_4mer_tot<-ifelse(ag1$tc_4mer_tot==0,NA,ag1$tc_4mer_tot)
ag2$tc_4mer_tot<-ifelse(ag2$tc_4mer_tot==0,NA,ag2$tc_4mer_tot)

seq_freq <- (Biostrings::oligonucleotideFrequency(tcCds, width=4, step = 1, as.prob = F, with.labels = TRUE))
rownames(seq_freq)<-names(tcCds)
seq_freq<-subset(seq_freq, select=colnames(ag1)[2:(ncol(ag1)-1)])
seq_freq[seq_freq==0]<-NA
# fin1<-ag1[2:(ncol(ag1)-1)]/ag1$tc_codon_tot*1e6
fin1<-ag1[2:(ncol(ag1)-1)]/ag1$tc_4mer_tot/seq_freq*1e6
fin1<-cbind(gene_id=ag1$gene_id, fin1, tc_4mer_tot=ag1$tc_4mer_tot)
# fin1<-subset(fin1, select=colnames(fin1)[!grepl("TGA",colnames(fin1)) & !grepl("TAA",colnames(fin1)) & !grepl("TAG",colnames(fin1))])
heat1<-melt(fin1, measure.vars = colnames(fin1)[2:257] , id.vars="gene_id")

# fin2<-ag2[2:(ncol(ag2)-1)]/ag2$tc_codon_tot*1e6
fin2<-ag2[2:(ncol(ag2)-1)]/ag2$tc_4mer_tot/seq_freq*1e6
fin2<-cbind(gene_id=ag2$gene_id, fin2, tc_4mer_tot=ag2$tc_4mer_tot)
# fin2<-subset(fin2, select=colnames(fin2)[!grepl("TGA",colnames(fin2)) & !grepl("TAA",colnames(fin2)) & !grepl("TAG",colnames(fin2))])
heat2<-melt(fin2, measure.vars = colnames(fin2)[2:257] , id.vars="gene_id")

# ggplot(heat1, aes(gene_id, variable))+geom_tile(aes(fill=log2(value)))+scale_fill_continuous(na.value = 'black')
# ggplot(heat2, aes(gene_id, variable))+geom_tile(aes(fill=log2(value)))+scale_fill_continuous(na.value = 'black')

colnames(fin1)[2:(ncol(fin1))]<-paste0("rep1_",colnames(fin1)[2:(ncol(fin1))])
colnames(fin2)[2:(ncol(fin2))]<-paste0("rep2_",colnames(fin2)[2:(ncol(fin2))])
fin<-merge(fin1, fin2, by="gene_id")

# setwd("~/Google Drive/hdlbp/")
setwd("E:/Google Drive/hdlbp/")
# mas<-read.delim("hdlbp_master_table_with_classes.txt", header=T)
mas<-read.delim("hdlbp_master_table_with_classes_uniq.txt", header=T)
inf<-subset(mas, select=c("gene_id","Symbol",
                          "tpm_cutoff","tc_CDS_norm","localization_cat",
                          "mean_te_293","loc_tar_CDS","tc_CDS_norm_cat",
                          "tc_transcript_norm",
                          "log2FoldChange.mem.cyt.293",
                          "log2FoldChange.mem.cyt.KO.293",
                          "log2FoldChange.ribo.rna.KO.WT"))

fin<-merge(fin, inf, by="gene_id")

ggplot(fin, aes(log2(rep1_tc_4mer_tot), log2(rep2_tc_4mer_tot)))+geom_point()

ave<-(fin[,2:258]+fin[,259:515])/2
colnames(ave)<-gsub("rep1_","",colnames(ave))
ave<-cbind(gene_id=fin$gene_id, ave, fin[,516:ncol(fin)])
colnames(ave)[2:257]<-gsub("T","U",colnames(ave)[2:257])

ave<-subset(ave, tc_4mer_tot>=10) ## thresholdat least 10 T-C per CDS

heat<-melt(ave, measure.vars = colnames(ave)[c(2:257)] , 
           id.vars=c("gene_id","tpm_cutoff","tc_CDS_norm",
                     "localization_cat","mean_te_293","loc_tar_CDS","tc_CDS_norm_cat",
                     "log2FoldChange.mem.cyt.293",
                     "log2FoldChange.mem.cyt.KO.293",
                     "log2FoldChange.ribo.rna.KO.WT",
                     "tc_transcript_norm"))


heat$gene_id<-factor(heat$gene_id, levels=unique(heat$gene_id[order(heat$log2FoldChange.mem.cyt.293, decreasing = T)]), ordered = T)
ggplot(heat, aes(gene_id, variable))+geom_tile(aes(fill=log2(value)))+
  scale_fill_continuous(limits=c(0,20),na.value = 'black')+theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
                                                                 axis.text.x=element_blank(),
                                                                 axis.ticks.x=element_blank())

ggplot(subset(heat, !is.na(localization_cat)), aes(gene_id, variable))+geom_tile(aes(fill=log2(value)))+
  scale_fill_continuous(limits=c(0,20),na.value = 'black')+theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
                                                                 axis.text.x=element_blank(),
                                                                 axis.ticks.x=element_blank())+facet_wrap(~localization_cat, scales="free_x")
maim<-aggregate(value~variable+localization_cat, data=heat, mean, na.rm=T)
maimCyt<-subset(maim, localization_cat=="cytosolic")
maimMem<-subset(maim, localization_cat=="membrane")

plot<-merge(maimCyt,maimMem, by="variable")
ggplot(plot, aes(value.x, value.y))+geom_text(aes(label=variable))+geom_abline()
plot$variable<-as.character(plot$variable)
plot$labels<-ifelse(plot$variable=="UGAA" | plot$variable=="CCUG" | plot$variable=="CUUC" | plot$variable=="CCUC" | plot$variable=="CUCU" | plot$variable=="AGAA" | plot$variable=="UGGA" ,
                    plot$variable, NA)
ggplot(plot, aes(value.x, value.y, colour=labels))+geom_text(aes(label=variable))+geom_abline()

cdsAdj<-plot
##TC in 3'UTR- 4mer binding adjacent

dat<-merge(tc, ham, by.x="transcript_id", by.y="transcript")
length(unique(dat$transcript_id))
length(unique(dat$gene_id))

dat$tc_from_start<-dat$tc_stop-dat$cds_stop+1
dat<-subset(dat, tc_from_start<=l_utr3 & tc_from_start>0)

fastapath<-"hg19bt1.transcripts.fa"
seqTrans <- Biostrings::readDNAStringSet(fastapath, format = "fasta", use.names = TRUE)

sub<- seqTrans[ham$transcript]

seqCds <- Biostrings::subseq(sub, start = ham$cds_start, end = ham$cds_stop)
seqUtr3<- Biostrings::subseq(sub, start = ham$cds_stop+1, end = ham$l_tr)
seqUtr5<- Biostrings::subseq(sub, start = 1, end = ham$l_utr5)

tcCds<-seqUtr3[unique(dat$transcript_id)]
dat$tc_from_start_4mer<-ifelse(dat$tc_from_start+4>dat$l_utr3 , NA, dat$tc_from_start)
dat<-subset(dat, !is.na(tc_from_start_4mer))
#here set position with tc_from_start_4mer
# dat<-subset(dat, tc_from_start_4mer-4>1)
dat$tc_4mer<-as.character(Biostrings::subseq(tcCds[as.character(dat$transcript_id)],
                                             start=dat$tc_from_start_4mer+1,end=dat$tc_from_start_4mer+4))


ag1<-dcast(dat[,c("gene_id","tc_4mer","norm_tc_num1")], gene_id~tc_4mer,sum,value.var="norm_tc_num1")
ag2<-dcast(dat[,c("gene_id","tc_4mer","norm_tc_num2")], gene_id~tc_4mer,sum,value.var="norm_tc_num2")

ag1$tc_4mer_tot<-rowSums(ag1[2:ncol(ag1)],na.rm=T)
ag2$tc_4mer_tot<-rowSums(ag2[2:ncol(ag2)],na.rm=T)
ag1$tc_4mer_tot<-ifelse(ag1$tc_4mer_tot==0,NA,ag1$tc_4mer_tot)
ag2$tc_4mer_tot<-ifelse(ag2$tc_4mer_tot==0,NA,ag2$tc_4mer_tot)

seq_freq <- (Biostrings::oligonucleotideFrequency(tcCds, width=4, step = 1, as.prob = F, with.labels = TRUE))
rownames(seq_freq)<-names(tcCds)
seq_freq<-subset(seq_freq, select=colnames(ag1)[2:(ncol(ag1)-1)])
seq_freq[seq_freq==0]<-NA
# fin1<-ag1[2:(ncol(ag1)-1)]/ag1$tc_codon_tot*1e6
fin1<-ag1[2:(ncol(ag1)-1)]/ag1$tc_4mer_tot/seq_freq*1e6
fin1<-cbind(gene_id=ag1$gene_id, fin1, tc_4mer_tot=ag1$tc_4mer_tot)
# fin1<-subset(fin1, select=colnames(fin1)[!grepl("TGA",colnames(fin1)) & !grepl("TAA",colnames(fin1)) & !grepl("TAG",colnames(fin1))])
heat1<-melt(fin1, measure.vars = colnames(fin1)[2:257] , id.vars="gene_id")

# fin2<-ag2[2:(ncol(ag2)-1)]/ag2$tc_codon_tot*1e6
fin2<-ag2[2:(ncol(ag2)-1)]/ag2$tc_4mer_tot/seq_freq*1e6
fin2<-cbind(gene_id=ag2$gene_id, fin2, tc_4mer_tot=ag2$tc_4mer_tot)
# fin2<-subset(fin2, select=colnames(fin2)[!grepl("TGA",colnames(fin2)) & !grepl("TAA",colnames(fin2)) & !grepl("TAG",colnames(fin2))])
heat2<-melt(fin2, measure.vars = colnames(fin2)[2:257] , id.vars="gene_id")

# ggplot(heat1, aes(gene_id, variable))+geom_tile(aes(fill=log2(value)))+scale_fill_continuous(na.value = 'black')
# ggplot(heat2, aes(gene_id, variable))+geom_tile(aes(fill=log2(value)))+scale_fill_continuous(na.value = 'black')

colnames(fin1)[2:(ncol(fin1))]<-paste0("rep1_",colnames(fin1)[2:(ncol(fin1))])
colnames(fin2)[2:(ncol(fin2))]<-paste0("rep2_",colnames(fin2)[2:(ncol(fin2))])
fin<-merge(fin1, fin2, by="gene_id")

# setwd("~/Google Drive/hdlbp/")
setwd("E:/Google Drive/hdlbp/")
# mas<-read.delim("hdlbp_master_table_with_classes.txt", header=T)
mas<-read.delim("hdlbp_master_table_with_classes_uniq.txt", header=T)
inf<-subset(mas, select=c("gene_id","Symbol",
                          "tpm_cutoff","tc_CDS_norm","localization_cat",
                          "mean_te_293","loc_tar_CDS","tc_CDS_norm_cat",
                          "tc_transcript_norm",
                          "log2FoldChange.mem.cyt.293",
                          "log2FoldChange.mem.cyt.KO.293",
                          "log2FoldChange.ribo.rna.KO.WT"))

fin<-merge(fin, inf, by="gene_id")

ggplot(fin, aes(log2(rep1_tc_4mer_tot), log2(rep2_tc_4mer_tot)))+geom_point()

ave<-(fin[,2:257]+fin[,258:513])/2
colnames(ave)<-gsub("rep1_","",colnames(ave))
ave<-cbind(gene_id=fin$gene_id, ave, fin[,514:ncol(fin)])
colnames(ave)[2:256]<-gsub("T","U",colnames(ave)[2:256])

ave<-subset(ave, tc_4mer_tot>=10) ## thresholdat least 10 T-C per CDS

heat<-melt(ave, measure.vars = colnames(ave)[c(2:257)] , 
           id.vars=c("gene_id","tpm_cutoff","tc_CDS_norm",
                     "localization_cat","mean_te_293","loc_tar_CDS","tc_CDS_norm_cat",
                     "log2FoldChange.mem.cyt.293",
                     "log2FoldChange.mem.cyt.KO.293",
                     "log2FoldChange.ribo.rna.KO.WT",
                     "tc_transcript_norm"))


heat$gene_id<-factor(heat$gene_id, levels=unique(heat$gene_id[order(heat$log2FoldChange.mem.cyt.293, decreasing = T)]), ordered = T)
ggplot(heat, aes(gene_id, variable))+geom_tile(aes(fill=log2(value)))+
  scale_fill_continuous(limits=c(0,20),na.value = 'black')+theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
                                                                 axis.text.x=element_blank(),
                                                                 axis.ticks.x=element_blank())

ggplot(subset(heat, !is.na(localization_cat)), aes(gene_id, variable))+geom_tile(aes(fill=log2(value)))+
  scale_fill_continuous(limits=c(0,20),na.value = 'black')+theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
                                                                 axis.text.x=element_blank(),
                                                                 axis.ticks.x=element_blank())+facet_wrap(~localization_cat, scales="free_x")
maim<-aggregate(value~variable+localization_cat, data=heat, mean, na.rm=T)
maimCyt<-subset(maim, localization_cat=="cytosolic")
maimMem<-subset(maim, localization_cat=="membrane")

plot<-merge(maimCyt,maimMem, by="variable")
ggplot(plot, aes(value.x, value.y))+geom_text(aes(label=variable))+geom_abline()
plot$variable<-as.character(plot$variable)
plot$labels<-ifelse(plot$variable=="UGAA" | plot$variable=="CCUG" | plot$variable=="CUUC" | plot$variable=="CCUC" | plot$variable=="CUCU" | plot$variable=="AGAA" | plot$variable=="UGGA" ,
                    plot$variable, NA)
ggplot(plot, aes(value.x, value.y, colour=labels))+geom_text(aes(label=variable))+geom_abline()



##TC in 3'UTR+stop- 4mer binding adjacent

dat<-merge(tc, ham, by.x="transcript_id", by.y="transcript")
length(unique(dat$transcript_id))
length(unique(dat$gene_id))


fastapath<-"hg19bt1.transcripts.fa"
seqTrans <- Biostrings::readDNAStringSet(fastapath, format = "fasta", use.names = TRUE)

sub<- seqTrans[ham$transcript]

tcCds<-sub[unique(dat$transcript_id)]
dat$interval_start<-dat$l_utr5+dat$l_cds-30

dat<-subset(dat, tc_start<interval_start)
#here set position with tc_from_start_4mer
# dat<-subset(dat, tc_from_start_4mer-4>1)
dat$tc_4mer<-as.character(Biostrings::subseq(tcCds[as.character(dat$transcript_id)],
                                             start=dat$tc_start+2,end=dat$tc_start+5))


ag1<-dcast(dat[,c("gene_id","tc_4mer","norm_tc_num1")], gene_id~tc_4mer,sum,value.var="norm_tc_num1")
ag2<-dcast(dat[,c("gene_id","tc_4mer","norm_tc_num2")], gene_id~tc_4mer,sum,value.var="norm_tc_num2")

ag1$tc_4mer_tot<-rowSums(ag1[2:ncol(ag1)],na.rm=T)
ag2$tc_4mer_tot<-rowSums(ag2[2:ncol(ag2)],na.rm=T)
ag1$tc_4mer_tot<-ifelse(ag1$tc_4mer_tot==0,NA,ag1$tc_4mer_tot)
ag2$tc_4mer_tot<-ifelse(ag2$tc_4mer_tot==0,NA,ag2$tc_4mer_tot)

seq_freq <- (Biostrings::oligonucleotideFrequency(tcCds, width=4, step = 1, as.prob = F, with.labels = TRUE))
rownames(seq_freq)<-names(tcCds)
seq_freq<-subset(seq_freq, select=colnames(ag1)[2:(ncol(ag1)-1)])
seq_freq[seq_freq==0]<-NA
# fin1<-ag1[2:(ncol(ag1)-1)]/ag1$tc_codon_tot*1e6
fin1<-ag1[2:(ncol(ag1)-1)]/ag1$tc_4mer_tot/seq_freq*1e6
fin1<-cbind(gene_id=ag1$gene_id, fin1, tc_4mer_tot=ag1$tc_4mer_tot)
# fin1<-subset(fin1, select=colnames(fin1)[!grepl("TGA",colnames(fin1)) & !grepl("TAA",colnames(fin1)) & !grepl("TAG",colnames(fin1))])
heat1<-melt(fin1, measure.vars = colnames(fin1)[2:257] , id.vars="gene_id")

# fin2<-ag2[2:(ncol(ag2)-1)]/ag2$tc_codon_tot*1e6
fin2<-ag2[2:(ncol(ag2)-1)]/ag2$tc_4mer_tot/seq_freq*1e6
fin2<-cbind(gene_id=ag2$gene_id, fin2, tc_4mer_tot=ag2$tc_4mer_tot)
# fin2<-subset(fin2, select=colnames(fin2)[!grepl("TGA",colnames(fin2)) & !grepl("TAA",colnames(fin2)) & !grepl("TAG",colnames(fin2))])
heat2<-melt(fin2, measure.vars = colnames(fin2)[2:257] , id.vars="gene_id")

# ggplot(heat1, aes(gene_id, variable))+geom_tile(aes(fill=log2(value)))+scale_fill_continuous(na.value = 'black')
# ggplot(heat2, aes(gene_id, variable))+geom_tile(aes(fill=log2(value)))+scale_fill_continuous(na.value = 'black')

colnames(fin1)[2:(ncol(fin1))]<-paste0("rep1_",colnames(fin1)[2:(ncol(fin1))])
colnames(fin2)[2:(ncol(fin2))]<-paste0("rep2_",colnames(fin2)[2:(ncol(fin2))])
fin<-merge(fin1, fin2, by="gene_id")

# setwd("~/Google Drive/hdlbp/")
setwd("E:/Google Drive/hdlbp/")
# mas<-read.delim("hdlbp_master_table_with_classes.txt", header=T)
mas<-read.delim("hdlbp_master_table_with_classes_uniq.txt", header=T)
inf<-subset(mas, select=c("gene_id","Symbol",
                          "tpm_cutoff","tc_CDS_norm","localization_cat",
                          "mean_te_293","loc_tar_CDS","tc_CDS_norm_cat",
                          "tc_transcript_norm",
                          "log2FoldChange.mem.cyt.293",
                          "log2FoldChange.mem.cyt.KO.293",
                          "log2FoldChange.ribo.rna.KO.WT"))

fin<-merge(fin, inf, by="gene_id")

ggplot(fin, aes(log2(rep1_tc_4mer_tot), log2(rep2_tc_4mer_tot)))+geom_point()

ave<-(fin[,2:257]+fin[,258:513])/2
colnames(ave)<-gsub("rep1_","",colnames(ave))
ave<-cbind(gene_id=fin$gene_id, ave, fin[,514:ncol(fin)])
colnames(ave)[2:256]<-gsub("T","U",colnames(ave)[2:256])

ave<-subset(ave, tc_4mer_tot>=10) ## thresholdat least 10 T-C per CDS

heat<-melt(ave, measure.vars = colnames(ave)[c(2:257)] , 
           id.vars=c("gene_id","tpm_cutoff","tc_CDS_norm",
                     "localization_cat","mean_te_293","loc_tar_CDS","tc_CDS_norm_cat",
                     "log2FoldChange.mem.cyt.293",
                     "log2FoldChange.mem.cyt.KO.293",
                     "log2FoldChange.ribo.rna.KO.WT",
                     "tc_transcript_norm"))


heat$gene_id<-factor(heat$gene_id, levels=unique(heat$gene_id[order(heat$log2FoldChange.mem.cyt.293, decreasing = T)]), ordered = T)
ggplot(heat, aes(gene_id, variable))+geom_tile(aes(fill=log2(value)))+
  scale_fill_continuous(limits=c(0,20),na.value = 'black')+theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
                                                                 axis.text.x=element_blank(),
                                                                 axis.ticks.x=element_blank())

ggplot(subset(heat, !is.na(localization_cat)), aes(gene_id, variable))+geom_tile(aes(fill=log2(value)))+
  scale_fill_continuous(limits=c(0,20),na.value = 'black')+theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
                                                                 axis.text.x=element_blank(),
                                                                 axis.ticks.x=element_blank())+facet_wrap(~localization_cat, scales="free_x")
maim<-aggregate(value~variable+localization_cat, data=heat, mean, na.rm=T)
maimCyt<-subset(maim, localization_cat=="cytosolic")
maimMem<-subset(maim, localization_cat=="membrane")

plot<-merge(maimCyt,maimMem, by="variable")
ggplot(plot, aes(value.x, value.y))+geom_text(aes(label=variable))+geom_abline()
plot$variable<-as.character(plot$variable)
plot$labels<-ifelse(plot$variable=="UGAA" | plot$variable=="CCUG" | plot$variable=="CUUC" | plot$variable=="CCUC" | plot$variable=="CUCU" | plot$variable=="AGAA" | plot$variable=="UGGA" ,
                    plot$variable, NA)
ggplot(plot, aes(value.x, value.y, colour=labels))+geom_text(aes(label=variable))+geom_abline()

stopUtr3<-plot

adj<-merge(cdsAdj, stopUtr3, by="variable")
colnames(adj)<-c("mer", "loc1", "normTC_cds_cyt","loc2", "normTC_cds_mem","labels1", "loc3", "normTC_utr3stop_cyt", "loc4", "normTC_utr3stop_mem", "labels2" )
ggplot(adj, aes(normTC_cds_mem, normTC_utr3stop_cyt, colour=labels1))+geom_text(aes(label=mer))+geom_abline()+xlim(0,6500)+ylim(0,6500)
ggplot(adj, aes(normTC_cds_cyt, normTC_utr3stop_cyt, colour=labels1))+geom_text(aes(label=mer))+geom_abline()+xlim(0,6500)+ylim(0,6500)
ggplot(adj, aes(normTC_cds_cyt, normTC_utr3stop_mem, colour=labels1))+geom_text(aes(label=mer))+geom_abline()+xlim(0,6500)+ylim(0,6500)
ggplot(adj, aes(normTC_cds_mem, normTC_cds_cyt, colour=labels1))+geom_text(aes(label=mer))+geom_abline()+xlim(0,6500)+ylim(0,6500)
ggplot(adj, aes(normTC_utr3stop_mem, normTC_utr3stop_cyt, colour=labels1))+geom_text(aes(label=mer))+geom_abline()+xlim(0,6500)+ylim(0,6500)

##redo 4mers and crosslinking

library(ggplot2)
library(reshape2)
library(corrplot)
library(Biostrings)
library(dplyr)
# setwd("/Volumes/landthaler/pcp/projects/miha/HDLBP/all_clip_data/mapping_trans/")
# setwd("/Volumes/landthaler/pcp/projects/miha/HDLBP/all_clip_data/reclip/mapping_trans/")
setwd("F:/landthaler/HDLBP/all_clip_data/reclip/mapping_trans/")
len<-read.delim("transcript_cds_utr_lenght.txt", header=T)
rownames(len)<-len$transcript
len$cds_start<-len$l_utr5+1
len$cds_stop<-len$l_utr5+len$l_cds
len$frame_start<-len$l_utr5
len$frame_stop<-len$cds_start


#selecting highest expressed isoform
rsem1<-read.delim("T_293_1.isoforms.results", header=T)
colnames(rsem1)[2:length(colnames(rsem1))]<-paste0(colnames(rsem1)[2:length(colnames(rsem1))],"_","T_293_1")
rsem2<-read.delim("T_293_2.isoforms.results", header=T)
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


tc<-read.delim("reproducible.hdlbp.TCseq.bed", header=F)
colnames(tc)<-c("transcript_id","tc_start","tc_stop","tc_num1","all_reads1","tc_num2","all_reads2","seq")
library(DESeq2)
fac<-estimateSizeFactorsForMatrix(tc[,c("all_reads1","all_reads2")])
tc$norm_tc_num1<-tc$tc_num1/fac[1]
tc$norm_tc_num2<-tc$tc_num2/fac[2]
# ggplot(tc, aes(log2(norm_tc_num1),log2(norm_tc_num2)))+geom_point()+geom_abline(slope=1)
tc<-subset(tc, tc_num1/all_reads1!=1 & tc_num2/all_reads2!=1)

dat<-merge(tc, ham, by.x="transcript_id", by.y="transcript")
length(unique(dat$transcript_id))
length(unique(dat$gene_id))

dat$tc_region<-ifelse(dat$tc_stop>=dat$cds_start & dat$tc_stop<=dat$cds_stop, "cds",
                      ifelse(dat$tc_stop<dat$cds_start, "utr5", 
                             ifelse(dat$tc_stop>dat$cds_stop, "utr3", NA)))
#which position in 4mer
dat$tc_4mer_pos4<-substr(toupper(dat$seq),1, 4)
dat$tc_4mer_pos3<-substr(toupper(dat$seq),2, 5)
dat$tc_4mer_pos2<-substr(toupper(dat$seq),3, 6)
dat$tc_4mer_pos1<-substr(toupper(dat$seq),4, 7)

setwd("E:/Google Drive/hdlbp/")
# mas<-read.delim("hdlbp_master_table_with_classes.txt", header=T)
mas<-read.delim("hdlbp_master_table_with_classes_uniq.txt", header=T)
inf<-subset(mas, select=c("gene_id","Symbol",
                          "tpm_cutoff","tc_CDS_norm","localization_cat",
                          "mean_te_293","loc_tar_CDS","tc_CDS_norm_cat",
                          "tc_transcript_norm",
                          "log2FoldChange.mem.cyt.293",
                          "log2FoldChange.mem.cyt.KO.293",
                          "log2FoldChange.ribo.rna.KO.WT"))

dat<-merge(dat, inf, by="gene_id")


ggplot(subset(dat, !is.na(localization_cat)), 
       aes(tc_4mer_pos4, fill=localization_cat))+geom_bar(position="dodge")+
  facet_wrap(~tc_region, scales = "free_x")+coord_flip()

ggplot(subset(dat, !is.na(localization_cat)), 
       aes(tc_4mer_pos3, fill=localization_cat))+geom_bar(position="dodge")+
  facet_wrap(~tc_region, scales = "free_x")+coord_flip()

ggplot(subset(dat, !is.na(localization_cat)), 
       aes(tc_4mer_pos2, fill=localization_cat))+geom_bar(position="dodge")+
  facet_wrap(~tc_region, scales = "free_x")+coord_flip()

ggplot(subset(dat, !is.na(localization_cat)), 
       aes(tc_4mer_pos1, fill=localization_cat))+geom_bar(position="dodge")+
  facet_wrap(~tc_region, scales = "free_x")+coord_flip()

mer<-melt(dat, measure.vars = colnames(dat)[grepl("tc_4mer_", colnames(dat))],
          id.vars = colnames(dat)[c(1:21,27:ncol(dat))], variable.name = "tc_4mer_pos", value.name = "tc_4mer")


ggplot(subset(mer, !is.na(localization_cat)), 
       aes(tc_4mer, fill=localization_cat))+geom_bar(position="dodge")+
  facet_wrap(~tc_region, scales = "free_x")+coord_flip()+
  theme(text = element_text(size=7))

ggplot(subset(mer, !is.na(localization_cat) & norm_tc_num1>=5 & norm_tc_num2>=5), 
       aes(tc_4mer, fill=localization_cat))+geom_bar(position="dodge")+
  facet_wrap(~tc_region, scales = "free_x")+coord_flip()+
  theme(text = element_text(size=7))

ggplot(subset(mer, !is.na(localization_cat) & TPM_transcript>=10), 
       aes(tc_4mer, fill=localization_cat))+geom_bar(position="dodge")+
  facet_wrap(~tc_region, scales = "free_x")+coord_flip()+
  theme(text = element_text(size=7))

mer$count<-ifelse(is.na(mer$localization_cat), NA,
                  ifelse(mer$norm_tc_num1>1 & mer$norm_tc_num2>1 & mer$TPM_transcript>=10 , 1, NA))
countmers_gene<-dcast(subset(mer, !is.na(count), select=c("transcript_id", "tc_4mer", "count")),transcript_id~tc_4mer,sum,value.var="count")

countmers<-aggregate(count~tc_4mer+tc_region+localization_cat, sum, data=mer)
ggplot(countmers, aes(tc_4mer,count, fill =tc_region))+geom_bar(stat="identity", position="dodge")+facet_wrap(~localization_cat, scales="free")+coord_flip()+
  theme(text = element_text(size=7))

cyt_count<-subset(countmers, localization_cat=="cytosolic")
cyt_count$frac_tot<-cyt_count$count/sum(cyt_count$count)
agr<-aggregate(frac_tot~tc_4mer, sum, data=cyt_count)
agr<-agr[order(agr$frac_tot, decreasing = T),]
agr<-agr[1:20, "tc_4mer"]
rows<-which(cyt_count$tc_4mer %in% agr)
cyt_count<-cyt_count[rows,]
cyt_count$tc_4mer<-factor(cyt_count$tc_4mer, levels=rev(agr))
ggplot(cyt_count, aes(tc_4mer,frac_tot, fill =tc_region))+geom_bar(stat="identity", position="dodge")+coord_flip()+
  theme(text = element_text(size=8))

mem_count<-subset(countmers, localization_cat=="membrane")
mem_count$frac_tot<-mem_count$count/sum(mem_count$count)
agr<-aggregate(frac_tot~tc_4mer, sum, data=mem_count)
agr<-agr[order(agr$frac_tot, decreasing = T),]
agr<-agr[1:20, "tc_4mer"]
rows<-which(mem_count$tc_4mer %in% agr)
mem_count<-mem_count[rows,]
mem_count$tc_4mer<-factor(mem_count$tc_4mer, levels=rev(agr))
ggplot(mem_count, aes(tc_4mer,frac_tot, fill =tc_region))+geom_bar(stat="identity", position="dodge")+coord_flip()+
  theme(text = element_text(size=8))

cyt_count<-subset(countmers, localization_cat=="cytosolic")
cyt_count$frac_tot<-cyt_count$count/sum(cyt_count$count)
mem_count<-subset(countmers, localization_cat=="membrane")
mem_count$frac_tot<-mem_count$count/sum(mem_count$count)

cyt_count$id<-paste0(cyt_count$tc_4mer, "_", cyt_count$tc_region)
mem_count$id<-paste0(mem_count$tc_4mer, "_", mem_count$tc_region)

cyt_mem<-merge(cyt_count, mem_count, by="id")
ggplot(cyt_mem, aes(frac_tot.x, frac_tot.y))+geom_text(aes(label=tc_4mer.x), size=2)+facet_wrap(~tc_region.x, scales="free")+geom_abline(slope=1)

cds_cyt_count<-subset(cyt_count, tc_region=="cds")
utr_cyt_count<-subset(cyt_count, tc_region=="utr3")
cds_mem_count<-subset(mem_count, tc_region=="cds")
utr_mem_count<-subset(mem_count, tc_region=="utr3")

cyt_mem_utr_cds<-merge(utr_cyt_count, cds_mem_count, by="tc_4mer")
ggplot(cyt_mem_utr_cds, aes(frac_tot.x, frac_tot.y))+geom_text(aes(label=tc_4mer), size=2)+geom_abline(slope=1)
ggplot(cyt_mem_utr_cds, aes(frac_tot.x, frac_tot.y))+geom_point()+geom_abline(slope=1)

cyt_mem_utr_utr<-merge(utr_cyt_count, utr_mem_count, by="tc_4mer")
ggplot(cyt_mem_utr_utr, aes(frac_tot.x, frac_tot.y))+geom_text(aes(label=tc_4mer), size=2)+geom_abline(slope=1)
ggplot(cyt_mem_utr_utr, aes(frac_tot.x, frac_tot.y))+geom_text(aes(label=tc_4mer), size=2)+geom_abline(slope=1)

countmers<-aggregate(count~tc_4mer+localization_cat, sum, data=subset(mer, tc_region!="utr5"))
ggplot(countmers, aes(tc_4mer,count, fill =localization_cat))+geom_bar(stat="identity", position="dodge")+coord_flip()+
  theme(text = element_text(size=7))

cyt_count<-subset(countmers, localization_cat=="cytosolic")
cyt_count$frac_tot<-cyt_count$count/sum(cyt_count$count)
cyt_count$loc<-"cytosolic"
mem_count<-subset(countmers, localization_cat=="membrane")
mem_count$frac_tot<-mem_count$count/sum(mem_count$count)
mem_count$loc<-"membrane"
both<-rbind(cyt_count,mem_count)
agr<-aggregate(frac_tot~tc_4mer, sum, data=both)
agr<-agr[order(agr$frac_tot, decreasing = T),]
agr<-agr[1:20, "tc_4mer"]
rows<-which(both$tc_4mer %in% agr)
both<-both[rows,]
both$tc_4mer<-factor(both$tc_4mer, levels=rev(agr))

ggplot(both, aes(tc_4mer,frac_tot, fill =loc))+geom_bar(stat="identity", position="dodge")+coord_flip()+
  theme(text = element_text(size=8))

sca<-merge(cyt_count, mem_count, by="tc_4mer")
ggplot(sca, aes(frac_tot.x, frac_tot.y ))+geom_text(aes(label=tc_4mer))+geom_abline(slope=1)

library(reshape2)

ag1_cds<-dcast(subset(mer, tc_region=="cds", select=c("transcript_id", "tc_4mer", "norm_tc_num1")),transcript_id~tc_4mer,sum,value.var="norm_tc_num1")
ag2_cds<-dcast(subset(mer, tc_region=="cds", select=c("transcript_id", "tc_4mer", "norm_tc_num2")), transcript_id~tc_4mer,sum,value.var="norm_tc_num2")

ag1_utr3<-dcast(subset(mer, tc_region=="utr3", select=c("transcript_id", "tc_4mer", "norm_tc_num1")), transcript_id~tc_4mer,sum,value.var="norm_tc_num1")
ag2_utr3<-dcast(subset(mer, tc_region=="utr3", select=c("transcript_id", "tc_4mer", "norm_tc_num2")), transcript_id~tc_4mer,sum,value.var="norm_tc_num2")

tot1<-aggregate(norm_tc_num1~transcript_id, sum, data=dat)
tot2<-aggregate(norm_tc_num2~transcript_id, sum, data=dat)

len_cds<-data.frame(transcript_id=dat$transcript_id,len=dat$l_cds)
len_cds<-subset(len_cds, !duplicated(transcript_id))
len_utr3<-data.frame(transcript_id=dat$transcript_id,len=dat$l_utr3)
len_utr3<-subset(len_utr3, !duplicated(transcript_id))

ag1_cds<-merge(ag1_cds, tot1, by="transcript_id")
ag1_cds<-merge(ag1_cds, len_cds, by="transcript_id")

ag2_cds<-merge(ag2_cds, tot2, by="transcript_id")
ag2_cds<-merge(ag2_cds, len_cds, by="transcript_id")

ag1_utr3<-merge(ag1_utr3, tot1, by="transcript_id")
ag1_utr3<-merge(ag1_utr3, len_utr3, by="transcript_id")

ag2_utr3<-merge(ag2_utr3, tot2, by="transcript_id")
ag2_utr3<-merge(ag2_utr3, len_utr3, by="transcript_id")


tcCds<-seqCds[unique(ag1_cds$transcript_id)]
seq_freq_cds<- (Biostrings::oligonucleotideFrequency(tcCds, width=4, step = 1, as.prob = F, with.labels = TRUE))
rownames(seq_freq_cds)<-names(tcCds)
seq_freq_cds<-subset(seq_freq_cds, select=colnames(ag1_cds)[2:(ncol(ag1_cds)-2)])
seq_freq_cds[seq_freq_cds==0]<-NA
# fin1_cds<-ag1_cds[2:(ncol(ag1_cds)-1)]/ag1_cds$tc_codon_tot*1e6
# fin1_cds<-ag1_cds[2:(ncol(ag1_cds)-2)]/ag1_cds$norm_tc_num1/ag1_cds$len/seq_freq_cds*1e9
# fin1_cds<-ag1_cds[2:(ncol(ag1_cds)-2)]/ag1_cds$len/seq_freq_cds*1e9
fin1_cds<-ag1_cds[2:(ncol(ag1_cds)-2)]/ag1_cds$len*1e9

seq_freq_cds<- (Biostrings::oligonucleotideFrequency(tcCds, width=4, step = 1, as.prob = F, with.labels = TRUE))
rownames(seq_freq_cds)<-names(tcCds)
seq_freq_cds<-subset(seq_freq_cds, select=colnames(ag2_cds)[2:(ncol(ag2_cds)-2)])
seq_freq_cds[seq_freq_cds==0]<-NA
# fin2_cds<-ag2_cds[2:(ncol(ag2_cds)-1)]/ag2_cds$tc_codon_tot*1e6
# fin2_cds<-ag2_cds[2:(ncol(ag2_cds)-2)]/ag2_cds$norm_tc_num2/ag2_cds$len/seq_freq_cds*1e9
# fin2_cds<-ag2_cds[2:(ncol(ag2_cds)-2)]/ag2_cds$len/seq_freq_cds*1e9
fin2_cds<-ag2_cds[2:(ncol(ag2_cds)-2)]/ag2_cds$len*1e9

# fin2_cds<-ag2_cds[2:(ncol(ag2_cds))]/seq_freq_cds*1e6

fin1_cds<-cbind(transcript_id=ag1_cds$transcript_id, fin1_cds, len=ag1_cds$len, seq_freq=seq_freq_cds)
heat1_cds<-melt(fin1_cds, measure.vars = colnames(fin1_cds)[2:176] , id.vars="transcript_id")
ggplot(heat1_cds, aes(variable, log2(value)))+geom_boxplot()+coord_flip()+theme(text = element_text(size=7))

fin2_cds<-cbind(transcript_id=ag1_cds$transcript_id, fin2_cds, len=ag2_cds$len, seq_freq=seq_freq_cds)
heat2_cds<-melt(fin2_cds, measure.vars = colnames(fin2_cds)[2:176] , id.vars="transcript_id")
ggplot(heat2_cds, aes(variable, log2(value)))+geom_boxplot()+coord_flip()+theme(text = element_text(size=7))


tcUtr3<-seqUtr3[unique(ag1_utr3$transcript_id)]
seq_freq_utr3<- (Biostrings::oligonucleotideFrequency(tcUtr3, width=4, step = 1, as.prob = F, with.labels = TRUE))
rownames(seq_freq_utr3)<-names(tcUtr3)
seq_freq_utr3<-subset(seq_freq_utr3, select=colnames(ag1_utr3)[2:(ncol(ag1_utr3)-2)])
seq_freq_utr3[seq_freq_utr3==0]<-NA
# fin1_utr3<-ag1_utr3[2:(ncol(ag1_utr3)-1)]/ag1_utr3$tc_codon_tot*1e6
# fin1_utr3<-ag1_utr3[2:(ncol(ag1_utr3)-2)]/ag1_utr3$norm_tc_num1/ag1_utr3$len/seq_freq_utr3*1e9
# fin1_utr3<-ag1_utr3[2:(ncol(ag1_utr3)-2)]/ag1_utr3$len/seq_freq_utr3*1e9
fin1_utr3<-ag1_utr3[2:(ncol(ag1_utr3)-2)]/ag1_utr3$len*1e9

seq_freq_utr3<- (Biostrings::oligonucleotideFrequency(tcUtr3, width=4, step = 1, as.prob = F, with.labels = TRUE))
rownames(seq_freq_utr3)<-names(tcUtr3)
seq_freq_utr3<-subset(seq_freq_utr3, select=colnames(ag2_utr3)[2:(ncol(ag2_utr3)-2)])
seq_freq_utr3[seq_freq_utr3==0]<-NA
# fin1_utr3<-ag1_utr3[2:(ncol(ag1_utr3)-1)]/ag1_utr3$tc_codon_tot*1e6
# fin2_utr3<-ag2_utr3[2:(ncol(ag2_utr3)-2)]/ag2_utr3$norm_tc_num2/ag2_utr3$len/seq_freq_utr3*1e9
# fin2_utr3<-ag2_utr3[2:(ncol(ag2_utr3)-2)]/ag2_utr3$len/seq_freq_utr3*1e9
fin2_utr3<-ag2_utr3[2:(ncol(ag2_utr3)-2)]/ag2_utr3$len*1e9

fin1_utr3<-cbind(transcript_id=ag1_utr3$transcript_id, fin1_utr3, len=ag1_utr3$len, seq_freq=seq_freq_utr3)
heat1_utr3<-melt(fin1_utr3, measure.vars = colnames(fin1_utr3)[2:176] , id.vars="transcript_id")
ggplot(heat1_utr3, aes(variable, log2(value)))+geom_boxplot()+coord_flip()

fin2_utr3<-cbind(transcript_id=ag2_utr3$transcript_id, fin2_utr3, len=ag1_utr3$len, seq_freq=seq_freq_utr3)
heat2_utr3<-melt(fin2_utr3, measure.vars = colnames(fin2_utr3)[2:176] , id.vars="transcript_id")
ggplot(heat2_utr3, aes(variable, log2(value)))+geom_boxplot()+coord_flip()


heat1_cds$id<-paste0(heat1_cds$transcript_id,"_",heat1_cds$variable)
heat2_cds$id<-paste0(heat2_cds$transcript_id,"_",heat2_cds$variable)

heat1_utr3$id<-paste0(heat1_utr3$transcript_id,"_",heat1_utr3$variable)
heat2_utr3$id<-paste0(heat2_utr3$transcript_id,"_",heat2_utr3$variable)

# heat<-merge(heat1_cds, heat2_cds, by="id", all=T)
# ggplot(heat, aes(log2(value.x), log2(value.y)))+geom_point(shape=1, size=0.7)+facet_wrap(~variable.x)+geom_abline(slope=1)

heat<-merge(heat1_cds, heat1_utr3, by="id", all=T)
# ggplot(heat, aes(log2(value.x), log2(value.y)))+geom_point(shape=1, size=0.7)+facet_wrap(~variable.x)+geom_abline(slope=1)


# heat<-merge(heat1_cds, heat2_cds, by="id", all=T)


genes<-subset(dat, select=c("transcript_id", "gene_id"))
genes<-subset(genes, !duplicated(transcript_id))
heat<-merge(heat, genes, by.x="transcript_id.x", by.y="transcript_id")
heat<-merge(heat, inf, by="gene_id")
ggplot(heat, aes(log2(value.x), log2(value.y), colour=localization_cat))+geom_point(shape=1, size=0.7)+facet_wrap(~variable.x)+geom_abline(slope=1)

ggplot(subset(heat, !is.na(localization_cat) & variable.x=="ATTT"), aes(log2(value.x), log2(value.y), colour=localization_cat))+geom_point(shape=1, size=1)+facet_wrap(~variable.x)+geom_abline(slope=1)
ggplot(subset(heat, !is.na(localization_cat) & variable.x=="CCTC"), aes(log2(value.x), log2(value.y), colour=localization_cat))+geom_point(shape=1, size=1)+facet_wrap(~variable.x)+geom_abline(slope=1)

ggplot(subset(heat, !is.na(localization_cat) & variable.x=="CTTC"), aes(log2(value.x), log2(value.y), colour=localization_cat))+geom_point(shape=1, size=1)+facet_wrap(~variable.x)+geom_abline(slope=1)
ggplot(subset(heat, !is.na(localization_cat) & variable.x=="CTTC"), 
       aes(log2(value.x), log2(tpm_cutoff), colour=localization_cat))+geom_point(shape=1, size=1)+facet_wrap(~variable.x)+geom_abline(slope=1)
ggplot(subset(heat, !is.na(localization_cat) & variable.x=="CTTT"), 
       aes(log2(value.x/tpm_cutoff), log2(value.y/tpm_cutoff), colour=localization_cat))+geom_point(shape=1, size=1)+facet_wrap(~variable.x)+geom_abline(slope=1)
ggplot(subset(heat, !is.na(localization_cat) ), 
       aes(log2(value.x/tpm_cutoff), log2(value.y/tpm_cutoff), colour=localization_cat))+geom_point(shape=1, size=1)+facet_wrap(~variable.x)+geom_abline(slope=1)


fin1_cds$region<-"cds1"
fin2_cds$region<-"cds2"
fin1_utr3$region<-"utr31"
fin2_utr3$region<-"utr32"

fin<-rbind(fin1_cds, fin2_cds, fin1_utr3, fin2_utr3)
fin<-merge(fin, genes, by.x="transcript_id", by.y="transcript_id")
fin<-merge(fin, inf, by="gene_id")

fin_mel<-melt(fin, measure.vars = colnames(fin)[3:177],
              id.vars = colnames(fin)[c(1,2,354:ncol(fin))])
fin_mel<-subset(fin_mel, !is.na(localization_cat))
fin_mel$region_loc<-paste0(fin_mel$region,"_",fin_mel$localization_cat)

ggplot(subset(fin_mel, !is.na(localization_cat) & variable=="CCTC"), aes(log2(value), colour=region_loc))+geom_density()+facet_wrap(~variable)

rem<-subset(fin_mel, !is.na(localization_cat) & value>2 & tpm_cutoff>0) # how many crosslinks per 4mer
ggplot(subset(rem, tpm_cutoff>=10 & !is.na(localization_cat) & 
                variable=="GTTG"|variable=="CCTC"|variable=="GGCT"|variable=="ATCG"|variable=="CCTC"|variable=="TATA"|variable=="TATA"|variable=="GAAT"|variable=="TGAA"|variable=="CTGA"|variable=="TGGT"|variable=="ATCC"|variable=="TCCT"|variable=="TTCC"|variable=="TCTC"|variable=="CTCT"|variable==""|variable=="CATC"|variable=="CTTT"| variable=="TCTT"| variable=="TGGA"| variable=="CTTC"), 
       aes(variable, log2(value/tpm_cutoff),fill=region))+geom_boxplot(outlier.shape = NA)+facet_wrap(~localization_cat)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+coord_cartesian(ylim=c(7,23))



seqs<-melt(fin, measure.vars = colnames(fin)[179:353],
              id.vars = colnames(fin)[c(1,2,354:ncol(fin))])
seqs<-subset(seqs, !is.na(localization_cat))
seqs$region_loc<-paste0(seqs$region,"_",seqs$localization_cat)

seqs<-aggregate(value~variable+region_loc, sum, data=seqs)
seqs$variable<-gsub(".*\\.", "", seqs$variable)
seqs_reg<-aggregate(value~region_loc, sum, data=seqs)
seqs<-merge(seqs, seqs_reg, by="region_loc")
colnames(seqs)<-c("region_loc", "4mer", "count", "total")
seqs$freq<-seqs$count/seqs$total

cds_utr<- cbind(subset(seqs, region_loc=="cds1_cytosolic" | region_loc=="cds2_cytosolic" | region_loc=="cds1_membrane" | region_loc=="cds2_membrane" ),
             subset(seqs, region_loc=="utr31_cytosolic" | region_loc=="utr32_cytosolic" | region_loc=="utr31_membrane" | region_loc=="utr32_membrane"))
colnames(cds_utr)<-c("region_loc1", "kmer_1", "count1", "total1", "freq_cds", 
                     "region_loc2", "kmer_2", "count2", "total2", "freq_utr3")
cds_utr$loc<-gsub(".*_", "", cds_utr$region_loc1)
cds_utr$rep<-gsub("cds","",gsub("_.*", "", cds_utr$region_loc1))
ggplot(cds_utr, aes(freq_cds,freq_utr3))+geom_text(aes(label=kmer_1), size=2)+facet_wrap(~loc+rep)+geom_abline(slope=1)

cyt_mem<- cbind(subset(seqs, region_loc=="cds1_cytosolic" | region_loc=="cds2_cytosolic" | region_loc=="utr31_cytosolic" | region_loc=="utr32_cytosolic" ),
                subset(seqs, region_loc=="cds1_membrane" | region_loc=="cds2_membrane" | region_loc=="utr31_membrane" | region_loc=="utr32_membrane"))

colnames(cyt_mem)<-c("region_loc1", "kmer_1", "count1", "total1", "freq_cyto", 
                     "region_loc2", "kmer_2", "count2", "total2", "freq_mem")

cyt_mem$reg<-gsub("[0-9].*","", cyt_mem$region_loc1)
cyt_mem$rep<-gsub("cds","",gsub("_.*", "", cyt_mem$region_loc1))
cyt_mem$labels<-ifelse(cyt_mem$kmer_1=="TTCT"| cyt_mem$kmer_1=="TCCT", cyt_mem$kmer_1, NA)
ggplot(subset(cyt_mem, kmer_1!="TTTT"), aes(freq_cyto,freq_mem))+geom_text(aes(label=kmer_1, colour=labels), size=2)+facet_wrap(~reg+rep)+geom_abline(slope=1)

freqs<-aggregate((value/tpm_cutoff)~variable+region_loc, median, data=rem) # or sum???

freqs$kmer_region_loc<-paste0(freqs$variable, "_", freqs$region_loc)
seqs$kmer_region_loc<-paste0(seqs$`4mer`, "_", seqs$region_loc)

djus<-merge(freqs, seqs, by="kmer_region_loc")
colnames(djus)<-c("kmer_region_loc", "kmer", "region_loc1", "binding", "region_loc2", "kmer2", "count", "total", "freq")
ggplot(djus, aes(log2(binding), freq))+geom_text(aes(label=kmer), size=2)+facet_wrap(~region_loc1, ncol=4)
ggplot(djus, aes(log2(binding), freq))+geom_point(shape=1)+facet_wrap(~region_loc1, ncol=4)
djus$kmer<-as.character(djus$kmer)
djus$mark<-ifelse(djus$kmer=="CCTC"|djus$kmer=="CTGA"|djus$kmer=="TGGA"|djus$kmer=="TTCT"|djus$kmer=="CCTC", djus$kmer, NA)
ggplot(djus, aes(log2(binding), freq))+geom_text(aes(label=kmer, colour=mark), size=2)+facet_wrap(~region_loc1, ncol=4)+coord_cartesian(ylim=c(0,0.018))


bck_cds_mem<-subset(cds_utr, region_loc1=="cds1_membrane")
bck_utr_cyt<-subset(cds_utr, region_loc2=="utr31_cytosolic")

bck1<-merge(cyt_mem_utr_cds, bck_cds_mem, by.x="tc_4mer", by.y="kmer_1" )
bck1<-merge(bck1, bck_utr_cyt, by.x="tc_4mer", by.y="kmer_1", all=T )
colnames(bck1)[c(5,10,15,20,26,31)]<-c("bind_utr", "bind_cds", "cds_bck", "none1", "none2", "utr_bck")
ggplot(bck1, aes(bind_utr, bind_cds))+geom_text(aes(label=tc_4mer))

ggplot(bck1, aes(bind_utr, utr_bck))+geom_text(aes(label=tc_4mer))+geom_abline(slope=1)
ggplot(bck1, aes(bind_cds, cds_bck))+geom_text(aes(label=tc_4mer))+geom_abline(slope=1)
bck1$labels<-ifelse(bck1$tc_4mer=="CCCT"|bck1$tc_4mer=="CTCA"|bck1$tc_4mer=="CATC", bck1$tc_4mer, NA)
ggplot(bck1, aes(bind_cds/cds_bck, bind_utr/utr_bck))+geom_text(aes(label=tc_4mer, colour=labels))+geom_abline(slope=.5)

cyt_mem$reg<-gsub("_.*","",cyt_mem$region_loc1)
bck_cyt<-subset(cyt_mem, reg=="cds1" | reg=="utr31")
bck_cyt<-aggregate(freq_cyto~kmer_1, sum, data=bck_cyt)
bck_mem<-subset(cyt_mem, reg=="cds1" | reg=="utr31")
bck_mem<-aggregate(freq_mem~kmer_1, sum, data=bck_mem)

bck1<-merge(sca, bck_cyt, by.x="tc_4mer", by.y="kmer_1" )
bck1<-merge(bck1, bck_mem, by.x="tc_4mer", by.y="kmer_1", all=T )
colnames(bck1)[c(4,8,10,11)]<-c("bind_cyt", "bind_mem", "cyt_bck", "mem_bck")
ggplot(bck1, aes(bind_cyt, bind_mem))+geom_text(aes(label=tc_4mer))

ggplot(bck1, aes(bind_cyt, cyt_bck))+geom_text(aes(label=tc_4mer))+geom_abline(slope=1)
ggplot(bck1, aes(bind_mem, mem_bck))+geom_text(aes(label=tc_4mer))+geom_abline(slope=1)
bck1$labels<-ifelse(bck1$tc_4mer=="CCCT"|bck1$tc_4mer=="CTCA"|bck1$tc_4mer=="CATC", bck1$tc_4mer, NA)
ggplot(bck1, aes(bind_cyt/cyt_bck, bind_mem/mem_bck))+geom_text(aes(label=tc_4mer, colour=labels))+geom_abline(slope=1)



##redo 6mer enrichment in CDS/UTR3 mem/cyt
#run until tcCds

dat$tc_region<-ifelse(dat$tc_stop>=dat$cds_start & dat$tc_stop<=dat$cds_stop, "cds",
                      ifelse(dat$tc_stop<dat$cds_start, "utr5", 
                             ifelse(dat$tc_stop>dat$cds_stop, "utr3", NA)))
ggplot(subset(dat, TPM_transcript>=10), aes(log2(norm_tc_num1),log2(norm_tc_num2)))+geom_point()+facet_wrap(~tc_region)
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

dat<-merge(dat, inf, by="gene_id")

tot<-aggregate(cbind(norm_tc_num1, norm_tc_num2)~transcript_id, data=dat, sum)
totreg<-aggregate(cbind(norm_tc_num1, norm_tc_num2)~transcript_id+tc_region, data=dat, sum)


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
# freqs$labs<-ifelse(grepl("TTCT", freqs$kmer )| grepl("TCTT", freqs$kmer )| grepl("TCCT", freqs$kmer ), freqs$kmer, NA)
labs<-which(freqs$kmer %in% agr)
labs<-freqs[labs,]
labs$labs<-labs$kmer
freqs$labs<-NA
freqs<-rbind(freqs, labs)
ggplot(freqs, aes(freqCds_mem, freqCds_cyt))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()
ggplot(freqs, aes(freqCds_mem, freqUtr3_cyt))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()
ggplot(freqs, aes(freqCdsUtr3_mem, freqCdsUtr3_cyt))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()
ggplot(freqs, aes(freqCds_cyt, freqUtr3_cyt))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()

freqs$labs<-ifelse(grepl("TTCT", freqs$kmer )| grepl("TCTT", freqs$kmer )| grepl("TCCT", freqs$kmer ), freqs$kmer, NA)

ggplot(freqs, aes(freqCds_mem, freqCds_cyt))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()
ggplot(subset(freqs, !is.na(labs)), aes(freqCds_mem, freqCds_cyt))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()
ggplot(subset(freqs, !is.na(labs)), aes(freqCds_mem, freqUtr3_cyt))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()
ggplot(subset(freqs, !is.na(labs)), aes(freqCdsUtr3_mem, freqCdsUtr3_cyt))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()
ggplot(freqs, aes(freqCdsUtr3_mem, freqCdsUtr3_cyt))+geom_text(aes(label=kmer, colour=labs), size=2)+geom_abline()

##redo 7mers and crosslinking

library(ggplot2)
library(reshape2)
library(corrplot)
library(Biostrings)
library(dplyr)
# setwd("/Volumes/landthaler/pcp/projects/miha/HDLBP/all_clip_data/mapping_trans/")
# setwd("/Volumes/landthaler/pcp/projects/miha/HDLBP/all_clip_data/reclip/mapping_trans/")
setwd("D:/landthaler/HDLBP/all_clip_data/reclip/mapping_trans/")
len<-read.delim("transcript_cds_utr_lenght.txt", header=T)
rownames(len)<-len$transcript
len$cds_start<-len$l_utr5+1
len$cds_stop<-len$l_utr5+len$l_cds
len$frame_start<-len$l_utr5
len$frame_stop<-len$cds_start


#selecting highest expressed isoform
rsem1<-read.delim("T_293_1.isoforms.results", header=T)
colnames(rsem1)[2:length(colnames(rsem1))]<-paste0(colnames(rsem1)[2:length(colnames(rsem1))],"_","T_293_1")
rsem2<-read.delim("T_293_2.isoforms.results", header=T)
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


tc<-read.delim("reproducible.hdlbp.TCseq.bed", header=F)
colnames(tc)<-c("transcript_id","tc_start","tc_stop","tc_num1","all_reads1","tc_num2","all_reads2","seq")
library(DESeq2)
fac<-estimateSizeFactorsForMatrix(tc[,c("all_reads1","all_reads2")])
tc$norm_tc_num1<-tc$tc_num1/fac[1]
tc$norm_tc_num2<-tc$tc_num2/fac[2]
# ggplot(tc, aes(log2(norm_tc_num1),log2(norm_tc_num2)))+geom_point()+geom_abline(slope=1)
tc<-subset(tc, tc_num1/all_reads1!=1 & tc_num2/all_reads2!=1)

dat<-merge(tc, ham, by.x="transcript_id", by.y="transcript")
length(unique(dat$transcript_id))
length(unique(dat$gene_id))

dat$tc_region<-ifelse(dat$tc_stop>=dat$cds_start & dat$tc_stop<=dat$cds_stop, "cds",
                      ifelse(dat$tc_stop<dat$cds_start, "utr5", 
                             ifelse(dat$tc_stop>dat$cds_stop, "utr3", NA)))

setwd("E:/Google Drive/hdlbp/")
# mas<-read.delim("hdlbp_master_table_with_classes.txt", header=T)
mas<-read.delim("hdlbp_master_table_with_classes_uniq.txt", header=T)
inf<-subset(mas, select=c("gene_id","Symbol",
                          "tpm_cutoff","tc_CDS_norm","localization_cat",
                          "mean_te_293","loc_tar_CDS","tc_CDS_norm_cat",
                          "tc_transcript_norm",
                          "log2FoldChange.mem.cyt.293",
                          "log2FoldChange.mem.cyt.KO.293",
                          "log2FoldChange.ribo.rna.KO.WT"))

dat<-merge(dat, inf, by="gene_id")

mer<-dat

mer$count<-ifelse(is.na(mer$localization_cat), NA,
                  ifelse(mer$norm_tc_num1>2 & mer$norm_tc_num2>2 & mer$TPM_transcript>=10 , 1, NA))
# countmers_gene<-dcast(subset(mer, !is.na(count), select=c("transcript_id", "seq", "count")),transcript_id~seq,sum,value.var="count")

countmers<-aggregate(count~seq+tc_region+localization_cat, sum, data=mer)
ggplot(countmers, aes(seq,count, fill =tc_region))+geom_bar(stat="identity", position="dodge")+facet_wrap(~localization_cat, scales="free")+coord_flip()+
  theme(text = element_text(size=7))

cyt_count<-subset(countmers, localization_cat=="cytosolic")
cyt_count$frac_tot<-cyt_count$count/sum(cyt_count$count)
agr<-aggregate(frac_tot~seq, sum, data=cyt_count)
agr<-agr[order(agr$frac_tot, decreasing = T),]
agr<-agr[1:20, "seq"]
rows<-which(cyt_count$seq %in% agr)
cyt_count<-cyt_count[rows,]
cyt_count$seq<-factor(cyt_count$seq, levels=rev(agr))
ggplot(cyt_count, aes(seq,frac_tot, fill =tc_region))+geom_bar(stat="identity", position="dodge")+coord_flip()+
  theme(text = element_text(size=8))

mem_count<-subset(countmers, localization_cat=="membrane")
mem_count$frac_tot<-mem_count$count/sum(mem_count$count)
agr<-aggregate(frac_tot~seq, sum, data=mem_count)
agr<-agr[order(agr$frac_tot, decreasing = T),]
agr<-agr[1:20, "seq"]
rows<-which(mem_count$seq %in% agr)
mem_count<-mem_count[rows,]
mem_count$seq<-factor(mem_count$seq, levels=rev(agr))
ggplot(mem_count, aes(seq,frac_tot, fill =tc_region))+geom_bar(stat="identity", position="dodge")+coord_flip()+
  theme(text = element_text(size=8))

cyt_count<-subset(countmers, localization_cat=="cytosolic")
cyt_count$frac_tot<-cyt_count$count/sum(cyt_count$count)
mem_count<-subset(countmers, localization_cat=="membrane")
mem_count$frac_tot<-mem_count$count/sum(mem_count$count)

cyt_count$id<-paste0(cyt_count$seq, "_", cyt_count$tc_region)
mem_count$id<-paste0(mem_count$seq, "_", mem_count$tc_region)

cyt_mem<-merge(cyt_count, mem_count, by="id")
ggplot(cyt_mem, aes(frac_tot.x, frac_tot.y))+geom_text(aes(label=seq.x), size=2)+facet_wrap(~tc_region.x, scales="free")+geom_abline(slope=1)

cds_cyt_count<-subset(cyt_count, tc_region=="cds")
utr_cyt_count<-subset(cyt_count, tc_region=="utr3")
cds_mem_count<-subset(mem_count, tc_region=="cds")
utr_mem_count<-subset(mem_count, tc_region=="utr3")

cyt_mem_utr_cds<-merge(utr_cyt_count, cds_mem_count, by="seq")
ggplot(cyt_mem_utr_cds, aes(frac_tot.x, frac_tot.y))+geom_text(aes(label=seq), size=2)+geom_abline(slope=1)
ggplot(cyt_mem_utr_cds, aes(frac_tot.x, frac_tot.y))+geom_point()+geom_abline(slope=1)

cyt_mem_utr_utr<-merge(utr_cyt_count, utr_mem_count, by="seq")
ggplot(cyt_mem_utr_utr, aes(frac_tot.x, frac_tot.y))+geom_text(aes(label=seq), size=2)+geom_abline(slope=1)
ggplot(cyt_mem_utr_utr, aes(frac_tot.x, frac_tot.y))+geom_text(aes(label=seq), size=2)+geom_abline(slope=1)

countmers<-aggregate(count~seq+localization_cat, sum, data=subset(mer, tc_region!="utr5"))
ggplot(countmers, aes(seq,count, fill =localization_cat))+geom_bar(stat="identity", position="dodge")+coord_flip()+
  theme(text = element_text(size=7))

cyt_count<-subset(countmers, localization_cat=="cytosolic")
cyt_count$frac_tot<-cyt_count$count/sum(cyt_count$count)
cyt_count$loc<-"cytosolic"
mem_count<-subset(countmers, localization_cat=="membrane")
mem_count$frac_tot<-mem_count$count/sum(mem_count$count)
mem_count$loc<-"membrane"
both<-rbind(cyt_count,mem_count)
agr<-aggregate(frac_tot~seq, sum, data=both)
agr<-agr[order(agr$frac_tot, decreasing = T),]
agr<-agr[1:20, "seq"]
rows<-which(both$seq %in% agr)
both<-both[rows,]
both$seq<-factor(both$seq, levels=rev(agr))

ggplot(both, aes(seq,frac_tot, fill =loc))+geom_bar(stat="identity", position="dodge")+coord_flip()+
  theme(text = element_text(size=8))

sca<-merge(cyt_count, mem_count, by="seq")
ggplot(sca, aes(frac_tot.x, frac_tot.y ))+geom_text(aes(label=seq))+geom_abline(slope=1)

mer$tc_4mer<-mer$seq


ag1_cds<-dcast(subset(mer, tc_region=="cds", select=c("transcript_id", "tc_4mer", "norm_tc_num1")),transcript_id~tc_4mer,sum,value.var="norm_tc_num1")
ag2_cds<-dcast(subset(mer, tc_region=="cds", select=c("transcript_id", "tc_4mer", "norm_tc_num2")), transcript_id~tc_4mer,sum,value.var="norm_tc_num2")

ag1_utr3<-dcast(subset(mer, tc_region=="utr3", select=c("transcript_id", "tc_4mer", "norm_tc_num1")), transcript_id~tc_4mer,sum,value.var="norm_tc_num1")
ag2_utr3<-dcast(subset(mer, tc_region=="utr3", select=c("transcript_id", "tc_4mer", "norm_tc_num2")), transcript_id~tc_4mer,sum,value.var="norm_tc_num2")

tot1<-aggregate(norm_tc_num1~transcript_id, sum, data=dat)
tot2<-aggregate(norm_tc_num2~transcript_id, sum, data=dat)

len_cds<-data.frame(transcript_id=dat$transcript_id,len=dat$l_cds)
len_cds<-subset(len_cds, !duplicated(transcript_id))
len_utr3<-data.frame(transcript_id=dat$transcript_id,len=dat$l_utr3)
len_utr3<-subset(len_utr3, !duplicated(transcript_id))

ag1_cds<-merge(ag1_cds, tot1, by="transcript_id")
ag1_cds<-merge(ag1_cds, len_cds, by="transcript_id")

ag2_cds<-merge(ag2_cds, tot2, by="transcript_id")
ag2_cds<-merge(ag2_cds, len_cds, by="transcript_id")

ag1_utr3<-merge(ag1_utr3, tot1, by="transcript_id")
ag1_utr3<-merge(ag1_utr3, len_utr3, by="transcript_id")

ag2_utr3<-merge(ag2_utr3, tot2, by="transcript_id")
ag2_utr3<-merge(ag2_utr3, len_utr3, by="transcript_id")


ag1_cds<-data.frame(transcript_id=ag1_cds$transcript_id, ag1_cds[,2:(ncol(ag1_cds)-2)]/ag1_cds$len)
ag2_cds<-data.frame(transcript_id=ag2_cds$transcript_id, ag2_cds[,2:(ncol(ag2_cds)-2)]/ag2_cds$len)
ag1_utr3<-data.frame(transcript_id=ag1_utr3$transcript_id, ag1_utr3[,2:(ncol(ag1_utr3)-2)]/ag1_utr3$len)
ag2_utr3<-data.frame(transcript_id=ag2_utr3$transcript_id, ag2_utr3[,2:(ncol(ag2_utr3)-2)]/ag2_utr3$len)

fin1_cds<-melt(ag1_cds, measure.vars = colnames(ag1_cds)[2:ncol(ag1_cds)],
               id.vars="transcript_id", value.name = "tc", variable.name = "kmer")

fin2_cds<-melt(ag1_cds, measure.vars = colnames(ag2_cds)[2:ncol(ag2_cds)],
               id.vars="transcript_id", value.name = "tc", variable.name = "kmer")

fin1_utr3<-melt(ag1_utr3, measure.vars = colnames(ag1_utr3)[2:ncol(ag1_utr3)],
               id.vars="transcript_id", value.name = "tc", variable.name = "kmer")

fin2_utr3<-melt(ag1_utr3, measure.vars = colnames(ag2_utr3)[2:ncol(ag2_utr3)],
               id.vars="transcript_id", value.name = "tc", variable.name = "kmer")

fin1_cds$region<-"cds"
fin2_cds$region<-"cds"
fin1_utr3$region<-"utr3"
fin2_utr3$region<-"utr3"

fin1_cds$rep<-"1"
fin2_cds$rep<-"2"
fin1_utr3$rep<-"1"
fin2_utr3$rep<-"2"


fin<-rbind(fin1_cds, fin2_cds, fin1_utr3, fin2_utr3)


eas<-aggregate(tc~kmer+region+rep, data=fin, mean)


eas<-eas[order(eas$tc, decreasing = T),]
