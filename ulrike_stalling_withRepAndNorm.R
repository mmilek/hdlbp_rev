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

setwd("~/Google Drive/koshi_revision/riboprof/")
# setwd("E:/Google Drive/koshi_revision/riboprof/")

setwd("~/hdlbp_git/")
load("ribowaltz_human.RData")
sequences_biost <- Biostrings::readDNAStringSet("~/hdlbp_git/hg19bt1.transcripts.fa", 
                                                format = "fasta", use.names = TRUE)

mas<-read.delim("data/hdlbp_master_table_with_classes_uniq_tsig.txt", header=T)
memTrans<-subset(mas,localization_cat=="membrane" & tpm_cutoff>=10 & gene_biotype=="protein_coding", select="transcript")[,1]
cytTrans<-subset(mas,localization_cat=="cytosolic" & tpm_cutoff>=10 & gene_biotype=="protein_coding", select="transcript")[,1]
tarCDShi<-subset(mas,loc_tar_CDS=="membrane_tc>5.56 & tc<65.26" & tpm_cutoff>=10 & gene_biotype=="protein_coding", select="transcript")[,1]
Trans<-subset(mas,!is.na(localization_cat) & tpm_cutoff>=10 & gene_biotype=="protein_coding", select="transcript")[,1]


fin<-read.delim("data/psite_position_counts_chx.txt", header=T)
# fin<-read.delim("psite_position_counts_nochx.txt", header=T)
colnames(fin)<-gsub("X","", colnames(fin))

normSubset<-subset(fin, pos_from_start >=6 & pos_from_stop<=(-6)) #exclude first and last 2 codons
norm<-colSums(normSubset[,2:7], na.rm=T)

fin<-cbind(fin, fin[,2:7]/norm*1e6)
colnames(fin)[16:ncol(fin)]<-paste0("norm.",colnames(fin)[16:ncol(fin)])

fin<-data.frame(fin,length=width(sequences_biost[as.character(fin$transcript)]))
nt_start<-6 #exclude how many nucleotides from start
nt_stop<-6 #exclude how many nucleotides from stop
fin<-subset(fin,pos_from_start>nt_start&pos_from_stop<(-nt_stop) ) # additionaly filter some codons based on position?
nt_psite<-0 # change to look at a site or e site
fin$psite<-fin$psite+nt_psite 
fin<-subset(fin, length>=(psite+1))
fin$codon_start<-ifelse((fin$psite-1)%%3==0, fin$psite, 
                        ifelse((fin$psite-2)%%3==0,fin$psite-1,
                               fin$psite-2))
fin$psite_codon<-as.character(Biostrings::subseq(sequences_biost[as.character(fin$transcript)],
                                                 start=fin$codon_start,end=fin$codon_start+2))



# df <- subset(fin, pos_from_start%%3==0 & psite>=start_pos & psite<=stop_pos)
df <- subset(fin,  psite>=start_pos & psite<=stop_pos)

thr<-0
df <- subset(df, norm.293_1>=thr & norm.293_2>=thr & norm.guide1_1>=thr & norm.guide1_2>=thr & norm.guide2_1>=thr & norm.guide2_2>=thr)
# df <- subset(df, X293_1>=thr & X293_2>=thr & guide1_1>=thr & guide1_2>=thr & guide2_1>=thr & guide2_2>=thr)

ag_293_1<-dcast(df, transcript~psite_codon,sum,value.var="norm.293_1")
colnames(ag_293_1)<-gsub("T","U",colnames(ag_293_1))
ag_293_2<-dcast(df, transcript~psite_codon,sum,value.var="norm.293_2")
colnames(ag_293_1)<-gsub("T","U",colnames(ag_293_2))
ag_guide1_1<-dcast(df, transcript~psite_codon,sum,value.var="norm.guide1_1")
colnames(ag_guide1_1)<-gsub("T","U",colnames(ag_guide1_1))
ag_guide1_2<-dcast(df, transcript~psite_codon,sum,value.var="norm.guide1_2")
colnames(ag_guide1_2)<-gsub("T","U",colnames(ag_guide1_2))
ag_guide2_1<-dcast(df, transcript~psite_codon,sum,value.var="norm.guide2_1")
colnames(ag_guide2_1)<-gsub("T","U",colnames(ag_guide2_1))
ag_guide2_2<-dcast(df, transcript~psite_codon,sum,value.var="norm.guide2_2")
colnames(ag_guide2_2)<-gsub("T","U",colnames(ag_guide2_2))

transcript_subset<-Trans

cols<-which(ag_293_1$transcript %in% transcript_subset)
ag_293_1<-ag_293_1[cols,]
ag_293_2<-ag_293_2[cols,]
ag_guide1_1<-ag_guide1_1[cols,]
ag_guide1_2<-ag_guide1_2[cols,]
ag_guide2_1<-ag_guide2_1[cols,]
ag_guide2_2<-ag_guide2_2[cols,]

# ag_293_1<-dcast(df, transcript~psite_codon,sum,value.var="X293_1")
# colnames(ag_293_1)<-gsub("T","U",colnames(ag_293_1))
# ag_293_2<-dcast(df, transcript~psite_codon,sum,value.var="X293_2")
# colnames(ag_293_1)<-gsub("T","U",colnames(ag_293_2))
# ag_guide1_1<-dcast(df, transcript~psite_codon,sum,value.var="guide1_1")
# colnames(ag_guide1_1)<-gsub("T","U",colnames(ag_guide1_1))
# ag_guide1_2<-dcast(df, transcript~psite_codon,sum,value.var="guide1_2")
# colnames(ag_guide1_2)<-gsub("T","U",colnames(ag_guide1_2))
# ag_guide2_1<-dcast(df, transcript~psite_codon,sum,value.var="guide2_1")
# colnames(ag_guide2_1)<-gsub("T","U",colnames(ag_guide2_1))
# ag_guide2_2<-dcast(df, transcript~psite_codon,sum,value.var="guide2_2")
# colnames(ag_guide2_2)<-gsub("T","U",colnames(ag_guide2_2))

ag<-list(ag_293_1, ag_293_2, ag_guide1_1,ag_guide1_2,ag_guide2_1,ag_guide2_2)
names(ag)<-colnames(df[,2:7])

cod_aa <- data.frame(codon = c("GCC", "GCG", "GCU", "GCA", 
                               "AGA", "CGG", "AGG", "CGA", "CGC", "CGU", "AAC", "AAU", 
                               "GAC", "GAU", "UGC", "UGU", "CAA", "CAG", "GAG", "GAA", 
                               "GGC", "GGU", "GGA", "GGG", "CAC", "CAU", "AUA", "AUC", 
                               "AUU", "CUG", "CUA", "UUA", "CUU", "UUG", "CUC", "AAA", 
                               "AAG", "AUG", "UUC", "UUU", "CCG", "CCC", "CCU", "CCA", 
                               "AGC", "UCG", "UCU", "UCA", "UCC", "AGU", "UAG", "UAA", 
                               "UGA", "ACA", "ACC", "ACG", "ACU", "UGG", "UAU", "UAC", 
                               "GUA", "GUG", "GUU", "GUC"), aa = c("A", "A", "A", "A", 
                                                                   "R", "R", "R", "R", "R", "R", "N", "N", "D", "D", "C", 
                                                                   "C", "Q", "Q", "E", "E", "G", "G", "G", "G", "H", "H", 
                                                                   "I", "I", "I", "L", "L", "L", "L", "L", "L", "K", "K", 
                                                                   "M", "F", "F", "P", "P", "P", "P", "S", "S", "S", "S", 
                                                                   "S", "S", "*", "*", "*", "T", "T", "T", "T", "W", "Y", 
                                                                   "Y", "V", "V", "V", "V"))


compare<-function (data, threshold=0, condition, reference,site)
{
  # data<-ag
  # threshold<-0
  # condition<-"guide1_"
  # reference<-"X293_"
  # site<-"P"
  
  top<-lapply(data, function(x){row.names(x)<-x$transcript; x})
  top<-lapply(top, function(x) x[ ,-1])
  
  top<-lapply(top, function(x) x[rowSums(x)>threshold,])
  norm<-Map("/", top, lapply(top, rowSums))
  
  Sums<-lapply(norm, colMeans, na.rm=T)
  
  rep1<-(Sums[[paste0(condition,"1")]]-Sums[[paste0(reference,"1")]])/Sums[[paste0(reference,"1")]]
  rep2<-(Sums[[paste0(condition,"2")]]-Sums[[paste0(reference,"2")]])/Sums[[paste0(reference,"2")]]
  freq1.1<-Sums[[paste0(condition,"1")]]
  freq1.2<-Sums[[paste0(reference,"1")]]
  freq2.1<-Sums[[paste0(condition,"2")]]
  freq2.2<-Sums[[paste0(reference,"2")]]
  
  df<-data.frame(codon=names(rep1),rep1, rep2, freq1.1, freq1.2, freq2.1, freq2.2)
  df$av<-rowMeans(df[,2:3])
  df$sd<-apply(df[,2:3],1,sd)
  df$AvgFreq<-rowMeans(df[,4:7])
  df$pos3<-ifelse(substr(df$codon,start = 3,stop=3)=="A"|substr(df$codon,start = 3,stop=3)=="U","A/U","G/C" )
  df<-merge(df, cod_aa, by="codon")
  
  df$codon<-factor(df$codon, levels=df$codon[order(df$av)])
  
  p1<-ggplot(df[df$codon!="UGA"& df$codon!="UAG"&df$codon!="UAA" ,], aes(codon,av))+
    geom_point(aes(size=AvgFreq,colour=pos3))+
    geom_errorbar(aes(ymin=av-sd,ymax=av+sd))+
    theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.2))+
    xlab("")+ylab("mean_codon shift")+
    ggtitle(paste0(condition," vs ",reference," - ",site," site"))+ylim(-2,2)
  
  df$codon<-factor(df$codon, levels=df$codon[order(df$aa)])
  df$codon<-factor(df$codon, levels=df$codon[order(df$codon)])
  
  
  p2<-ggplot(df[df$codon!="UGA"& df$codon!="UAG"&df$codon!="UAA" ,], aes(codon,av))+
    facet_grid(~aa, space="free_x", scales="free_x", switch="x")+
    geom_point(aes(size=AvgFreq,colour=aa))+
    geom_errorbar(aes(ymin=av-sd,ymax=av+sd))+
    scale_colour_discrete(guide = "none")+ 
    theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.2))+
    xlab("")+ylab("mean_codon shift")+
    ggtitle(paste0(condition," vs ",reference," - ",site," site"))+ylim(-.3,.3)
  
  return(p2)
  
}
thres<-0
compare(ag, thres, "guide1_", "X293_", "P")
compare(ag, thres, "guide2_", "X293_", "P")

thres<-0
compare(ag, thres, "guide1_", "X293_", "E")
compare(ag, thres, "guide2_", "X293_", "E")

thres<-0
compare(ag, thres, "guide1_", "X293_", "A")
compare(ag, thres, "guide2_", "X293_", "A")




#done till here


##try to add TC codon crosslink information
### dat variable from hdlbp_codons_backup.R

dat$id<-paste0(dat$transcript_id, "_",dat$codon_start,"_", dat$cds_start,"_", dat$cds_stop)
length(unique(dat$id))

fin$codon_start_cds<-fin$codon_start-fin$start_pos+1
fin$id<-paste0(fin$transcript,"_",fin$codon_start_cds,"_", fin$start_pos, "_", fin$stop_pos)
mer<-merge(dat, fin, by="id", all=F )
nrow(subset(mer, tc_codon==psite_codon))
# df <- subset(mer, pos_from_start%%3==0 & psite>=start_pos & psite<=stop_pos)
# df <- subset(mer,  psite>=start_pos & psite<=stop_pos & tc_num1>2& tc_num2>2)
df <- subset(mer,  psite>=start_pos & psite<=stop_pos & norm_tc_num1>5& norm_tc_num2>5)

thr<-0
df <- subset(df, norm.293_1>=thr & norm.293_2>=thr & norm.guide1_1>=thr & norm.guide1_2>=thr & norm.guide2_1>=thr & norm.guide2_2>=thr)
# df <- subset(df, X293_1>=thr & X293_2>=thr & guide1_1>=thr & guide1_2>=thr & guide2_1>=thr & guide2_2>=thr)

ag_293_1<-dcast(df, transcript~psite_codon,sum,value.var="norm.293_1")
colnames(ag_293_1)<-gsub("T","U",colnames(ag_293_1))
ag_293_2<-dcast(df, transcript~psite_codon,sum,value.var="norm.293_2")
colnames(ag_293_1)<-gsub("T","U",colnames(ag_293_2))
ag_guide1_1<-dcast(df, transcript~psite_codon,sum,value.var="norm.guide1_1")
colnames(ag_guide1_1)<-gsub("T","U",colnames(ag_guide1_1))
ag_guide1_2<-dcast(df, transcript~psite_codon,sum,value.var="norm.guide1_2")
colnames(ag_guide1_2)<-gsub("T","U",colnames(ag_guide1_2))
ag_guide2_1<-dcast(df, transcript~psite_codon,sum,value.var="norm.guide2_1")
colnames(ag_guide2_1)<-gsub("T","U",colnames(ag_guide2_1))
ag_guide2_2<-dcast(df, transcript~psite_codon,sum,value.var="norm.guide2_2")
colnames(ag_guide2_2)<-gsub("T","U",colnames(ag_guide2_2))

# transcript_subset<-Trans
# 
# cols<-which(ag_293_1$transcript %in% transcript_subset)
# ag_293_1<-ag_293_1[cols,]
# ag_293_2<-ag_293_2[cols,]
# ag_guide1_1<-ag_guide1_1[cols,]
# ag_guide1_2<-ag_guide1_2[cols,]
# ag_guide2_1<-ag_guide2_1[cols,]
# ag_guide2_2<-ag_guide2_2[cols,]

# ag_293_1<-dcast(df, transcript~psite_codon,sum,value.var="X293_1")
# colnames(ag_293_1)<-gsub("T","U",colnames(ag_293_1))
# ag_293_2<-dcast(df, transcript~psite_codon,sum,value.var="X293_2")
# colnames(ag_293_1)<-gsub("T","U",colnames(ag_293_2))
# ag_guide1_1<-dcast(df, transcript~psite_codon,sum,value.var="guide1_1")
# colnames(ag_guide1_1)<-gsub("T","U",colnames(ag_guide1_1))
# ag_guide1_2<-dcast(df, transcript~psite_codon,sum,value.var="guide1_2")
# colnames(ag_guide1_2)<-gsub("T","U",colnames(ag_guide1_2))
# ag_guide2_1<-dcast(df, transcript~psite_codon,sum,value.var="guide2_1")
# colnames(ag_guide2_1)<-gsub("T","U",colnames(ag_guide2_1))
# ag_guide2_2<-dcast(df, transcript~psite_codon,sum,value.var="guide2_2")
# colnames(ag_guide2_2)<-gsub("T","U",colnames(ag_guide2_2))

ag<-list(ag_293_1, ag_293_2, ag_guide1_1,ag_guide1_2,ag_guide2_1,ag_guide2_2)
names(ag)<-colnames(df)[grepl("norm\\.", colnames(df))]

cod_aa <- data.frame(codon = c("GCC", "GCG", "GCU", "GCA", 
                               "AGA", "CGG", "AGG", "CGA", "CGC", "CGU", "AAC", "AAU", 
                               "GAC", "GAU", "UGC", "UGU", "CAA", "CAG", "GAG", "GAA", 
                               "GGC", "GGU", "GGA", "GGG", "CAC", "CAU", "AUA", "AUC", 
                               "AUU", "CUG", "CUA", "UUA", "CUU", "UUG", "CUC", "AAA", 
                               "AAG", "AUG", "UUC", "UUU", "CCG", "CCC", "CCU", "CCA", 
                               "AGC", "UCG", "UCU", "UCA", "UCC", "AGU", "UAG", "UAA", 
                               "UGA", "ACA", "ACC", "ACG", "ACU", "UGG", "UAU", "UAC", 
                               "GUA", "GUG", "GUU", "GUC"), aa = c("A", "A", "A", "A", 
                                                                   "R", "R", "R", "R", "R", "R", "N", "N", "D", "D", "C", 
                                                                   "C", "Q", "Q", "E", "E", "G", "G", "G", "G", "H", "H", 
                                                                   "I", "I", "I", "L", "L", "L", "L", "L", "L", "K", "K", 
                                                                   "M", "F", "F", "P", "P", "P", "P", "S", "S", "S", "S", 
                                                                   "S", "S", "*", "*", "*", "T", "T", "T", "T", "W", "Y", 
                                                                   "Y", "V", "V", "V", "V"))

cod_aa<-subset(cod_aa, grepl("U", cod_aa$codon))


compare<-function (data, threshold=0, condition, reference,site)
{
  # data<-ag
  # threshold<-0
  # condition<-"norm.guide1_"
  # reference<-"norm.293_"
  # site<-"P"
  
  top<-lapply(data, function(x){row.names(x)<-x$transcript; x})
  top<-lapply(top, function(x) x[ ,-1])
  
  top<-lapply(top, function(x) x[rowSums(x)>threshold,])
  norm<-Map("/", top, lapply(top, rowSums))
  
  Sums<-lapply(norm, colMeans, na.rm=T)
  
  rep1<-(Sums[[paste0(condition,"1")]]-Sums[[paste0(reference,"1")]])/Sums[[paste0(reference,"1")]]
  rep2<-(Sums[[paste0(condition,"2")]]-Sums[[paste0(reference,"2")]])/Sums[[paste0(reference,"2")]]
  freq1.1<-Sums[[paste0(condition,"1")]]
  freq1.2<-Sums[[paste0(reference,"1")]]
  freq2.1<-Sums[[paste0(condition,"2")]]
  freq2.2<-Sums[[paste0(reference,"2")]]
  
  df<-data.frame(codon=names(rep1),rep1, rep2, freq1.1, freq1.2, freq2.1, freq2.2)
  df$av<-rowMeans(df[,2:3])
  df$sd<-apply(df[,2:3],1,sd)
  df$AvgFreq<-rowMeans(df[,4:7])
  df$pos3<-ifelse(substr(df$codon,start = 3,stop=3)=="A"|substr(df$codon,start = 3,stop=3)=="U","A/U","G/C" )
  df<-merge(df, cod_aa, by="codon")
  
  df$codon<-factor(df$codon, levels=df$codon[order(df$av)])
  
  p1<-ggplot(df[df$codon!="UGA"& df$codon!="UAG"&df$codon!="UAA" ,], aes(codon,av))+
    geom_point(aes(size=AvgFreq,colour=pos3))+
    geom_errorbar(aes(ymin=av-sd,ymax=av+sd))+
    theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.2))+
    xlab("")+ylab("mean_codon shift")+
    ggtitle(paste0(condition," vs ",reference," - ",site," site"))+ylim(-2,2)
  
  df$codon<-factor(df$codon, levels=df$codon[order(df$aa)])
  df$codon<-factor(df$codon, levels=df$codon[order(df$codon)])
  
  
  p2<-ggplot(df[df$codon!="UGA"& df$codon!="UAG"&df$codon!="UAA" ,], aes(codon,av))+
    facet_grid(~aa, space="free_x", scales="free_x", switch="x")+
    geom_point(aes(size=AvgFreq,colour=aa))+
    geom_errorbar(aes(ymin=av-sd,ymax=av+sd))+
    scale_colour_discrete(guide = FALSE)+ 
    theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.2))+
    xlab("")+ylab("mean_codon shift")+
    ggtitle(paste0(condition," vs ",reference," - ",site," site"))+ylim(-1.5,1.5)
  
  return(p2)
  
}
thres<-0
compare(ag, thres, "norm.guide1_", "norm.293_", "P")
compare(ag, thres, "norm.guide2_", "norm.293_", "P")

thres<-0
compare(ag, thres, "guide1_", "X293_", "P")
compare(ag, thres, "norm.guide2_", "norm.293_", "P")

thres<-0
compare(ag, thres, "norm.guide1_", "norm.293_", "E")
compare(ag, thres, "norm.guide2_", "norm.293_", "E")

compare(ag, thres, "norm.guide1_", "norm.293_", "A")
compare(ag, thres, "norm.guide2_", "norm.293_", "A")

##done till here


nrow(subset(mer, tc_num1>5 & tc_num2>5))




length(unique(mer$id))

setwd("/Volumes/landthaler/pcp/projects/miha/ulrike_stalling/")
codon_counts<-function(data, annotation, fastapath, nt_start=0, nt_stop=0,nt_psite=0, transcript_subset) #nt_psite gives how many nucleotides from the psite, asite=3, +1 site=6 etc...
{
  annotation<-anH
  fastapath<-"hg19bt1.transcripts.fa"
  transcript_subset<-Trans
  data<-reads_psite_list[[1]]
  nt_psite=0
  rownames(annotation) <- as.character(annotation$transcript)
  l.transcripts <- rownames(annotation)[which(annotation$l_utr5 > 
                                                0 & annotation$l_cds > 0 & annotation$l_cds%%3 == 0 & 
                                                annotation$l_utr3 > 0)]
  
  filter<-which(l.transcripts %in% transcript_subset)
  l.transcripts<-l.transcripts[filter]
  
  c.transcript <- l.transcripts
  
  sequences_biost <- Biostrings::readDNAStringSet(fastapath, 
                                                  format = "fasta", use.names = TRUE)
  
  fix<-data
  tub<-data.frame(fix,length=width(sequences_biost[as.character(fix$transcript)]))
  tub<-subset(tub,psite_from_start>nt_start&psite_from_stop<(-nt_stop))
  tub$psite<-tub$psite+nt_psite
  tub<-subset(tub, length.1>=(psite+2))
  fix<-tub[,-ncol(tub)]
  
  fix$psite_codon<-as.character(Biostrings::subseq(sequences_biost[as.character(fix$transcript)],
                                                   start=fix$psite,end=fix$psite+2))
  
  sub_sequences_biost <- sequences_biost[c.transcript]
  cds_biost <- Biostrings::subseq(sub_sequences_biost, start = annotation[names(sub_sequences_biost), 
                                                                          "l_utr5"] + 1, end = annotation[names(sub_sequences_biost), 
                                                                                                          "l_utr5"] + annotation[names(sub_sequences_biost), "l_cds"])
  #seq_freq <- (Biostrings::trinucleotideFrequency(cds_biost, step = 3, as.prob = TRUE, with.labels = TRUE, simplify.as = "collapsed")) * 1000
  #names(seq_freq) <- gsub("T", "U", names(seq_freq))
  df <- subset(fix, as.character(transcript) %in% 
                 c.transcript & psite_region == "cds" & psite_from_start%%3 ==  0)
  
  ag<-dcast(df[,c("transcript","psite_codon")], transcript~psite_codon,length,value.var="psite_codon")
  colnames(ag)<-gsub("T","U",colnames(ag))
  
  
  #temp_table <- as.data.frame(table(factor(df$psite_codon, levels = cod_lev)))
  
  return(ag)                                     
}
setwd("/Volumes/landthaler/pcp/projects/miha/ulrike_stalling/")
esite<-lapply(reads_psite_list,codon_counts, anH, "hg19bt1.transcripts.fa", nt_start,nt_stop,nt_psite=-3, Trans)
psite<-lapply(reads_psite_list,codon_counts, anH, "hg19bt1.transcripts.fa", nt_start,nt_stop,nt_psite=0, Trans)
asite<-lapply(reads_psite_list,codon_counts, anH, "hg19bt1.transcripts.fa", nt_start,nt_stop,nt_psite=3, Trans)


# psite<-lapply(reads_psite_list,codon_counts, anH, "hg19bt1.transcripts.fa", nt_start,nt_stop,nt_psite=0, memTrans)
# asite<-lapply(reads_psite_list,codon_counts, anH, "hg19bt1.transcripts.fa", nt_start,nt_stop,nt_psite=3, memTrans)
# 
# psite<-lapply(reads_psite_list,codon_counts, anH, "hg19bt1.transcripts.fa", nt_start,nt_stop,nt_psite=0, cytTrans)
# asite<-lapply(reads_psite_list,codon_counts, anH, "hg19bt1.transcripts.fa", nt_start,nt_stop,nt_psite=3, cytTrans)
# 
# psite<-lapply(reads_psite_list,codon_counts, anH, "hg19bt1.transcripts.fa", nt_start,nt_stop,nt_psite=0, tarCDShi)
# asite<-lapply(reads_psite_list,codon_counts, anH, "hg19bt1.transcripts.fa", nt_start,nt_stop,nt_psite=3, tarCDShi)
# 
# 
# plusOne<-lapply(reads_psite_list,codon_counts, anH, "hg19bt1.transcripts.fa", nt_start,nt_stop,nt_psite=6, memTrans)


#plusTwo<-lapply(reads_psite_list,codon_counts, an, "mm10bt1_vM14.transcripts.fa", nt_start,nt_stop,nt_psite=9)
#plusThree<-lapply(reads_psite_list,codon_counts, an, "mm10bt1_vM14.transcripts.fa", nt_start,nt_stop,nt_psite=12)

cod_aa <- data.frame(codon = c("GCC", "GCG", "GCU", "GCA", 
                               "AGA", "CGG", "AGG", "CGA", "CGC", "CGU", "AAC", "AAU", 
                               "GAC", "GAU", "UGC", "UGU", "CAA", "CAG", "GAG", "GAA", 
                               "GGC", "GGU", "GGA", "GGG", "CAC", "CAU", "AUA", "AUC", 
                               "AUU", "CUG", "CUA", "UUA", "CUU", "UUG", "CUC", "AAA", 
                               "AAG", "AUG", "UUC", "UUU", "CCG", "CCC", "CCU", "CCA", 
                               "AGC", "UCG", "UCU", "UCA", "UCC", "AGU", "UAG", "UAA", 
                               "UGA", "ACA", "ACC", "ACG", "ACU", "UGG", "UAU", "UAC", 
                               "GUA", "GUG", "GUU", "GUC"), aa = c("A", "A", "A", "A", 
                                                                   "R", "R", "R", "R", "R", "R", "N", "N", "D", "D", "C", 
                                                                   "C", "Q", "Q", "E", "E", "G", "G", "G", "G", "H", "H", 
                                                                   "I", "I", "I", "L", "L", "L", "L", "L", "L", "K", "K", 
                                                                   "M", "F", "F", "P", "P", "P", "P", "S", "S", "S", "S", 
                                                                   "S", "S", "*", "*", "*", "T", "T", "T", "T", "W", "Y", 
                                                                   "Y", "V", "V", "V", "V"))

codis<-c("GCC", "GCG", "GCU", "GCA", 
         "AGA", "CGG", "AGG", "CGA", "CGC", "CGU", "AAC", "AAU", 
         "GAC", "GAU", "UGC", "UGU", "CAA", "CAG", "GAG", "GAA", 
         "GGC", "GGU", "GGA", "GGG", "CAC", "CAU", "AUA", "AUC", 
         "AUU", "CUG", "CUA", "UUA", "CUU", "UUG", "CUC", "AAA", 
         "AAG", "AUG", "UUC", "UUU", "CCG", "CCC", "CCU", "CCA", 
         "AGC", "UCG", "UCU", "UCA", "UCC", "AGU", 
         "ACA", "ACC", "ACG", "ACU", "UGG", "UAU", "UAC", 
         "GUA", "GUG", "GUU", "GUC")
threshold<-50
esite_counts<-lapply(esite, function(x){row.names(x)<-x$transcript; x})
esite_counts<-lapply(esite_counts, function(x) x[ , codis])
esite_counts<-lapply(esite_counts, function(x) x[rowSums(x)>threshold,])

# sums<-lapply(esite_counts, function(x) sum(x)/1e6)
esite_counts<-Map("/", esite_counts, lapply(esite_counts, rowSums)) #normalize per transcript, perhaps better to total codon frequency
# esite_counts<-Map("/", esite_counts, sums)
esite_counts<-do.call("rbind", esite_counts)
esite_counts$transcript<-gsub(".*ENST","ENST", row.names(esite_counts))
esite_counts$sample<-gsub("\\.ENST.*","",row.names(esite_counts))

esite_counts<-merge(esite_counts, mas, by="transcript")


mel<-melt(esite_counts, measure.vars = colnames(esite_counts)[2:62], id.vars=c("transcript", "sample", "localization_cat", "loc_tar_CDS"))
ggplot(subset(mel, variable=="CUU") , aes(sample, log2(value), colour=sample))+geom_boxplot()+facet_wrap(~loc_tar_CDS)
ggplot(subset(mel, variable=="CUU") , aes(loc_tar_CDS, log2(value), colour=loc_tar_CDS))+geom_boxplot()+facet_wrap(~sample)

#
ggplot(subset(mel, variable=="CUU") , aes(localization_cat, log2(value), colour=localization_cat))+geom_boxplot()+facet_wrap(~sample)
ggplot(subset(mel, variable=="GAG") , aes(loc_tar_CDS, log2(value), colour=loc_tar_CDS))+geom_boxplot()+facet_wrap(~sample)

ggplot(subset(mel, variable=="CUC") , aes(sample, log2(value), colour=sample))+geom_boxplot()+facet_wrap(~localization_cat)

ggplot(subset(mel, variable=="GAG") , aes(sample, log2(value), colour=sample))+geom_boxplot()

threshold<-50
compare<-function (data, threshold=0, condition, reference,site)
{
  data<-esite
  threshold<-50
  condition<-"IgG_"
  reference<-"nochx.293_"
  site<-"E"
  top<-lapply(data, function(x){row.names(x)<-x$transcript; x})
  top<-lapply(top, function(x) x[ , codis])
  
  
  top<-lapply(top, function(x) x[rowSums(x)>threshold,])
  norm<-Map("/", top, lapply(top, rowSums))
  
  Sums<-lapply(norm, colMeans, na.rm=T)
  
  rep1<-(Sums[[paste0(condition,"1")]]-Sums[[paste0(reference,"1")]])/Sums[[paste0(reference,"1")]]
  rep2<-(Sums[[paste0(condition,"2")]]-Sums[[paste0(reference,"2")]])/Sums[[paste0(reference,"2")]]
  freq1.1<-Sums[[paste0(condition,"1")]]
  freq1.2<-Sums[[paste0(reference,"1")]]
  freq2.1<-Sums[[paste0(condition,"2")]]
  freq2.2<-Sums[[paste0(reference,"2")]]
  
  df<-data.frame(codon=names(rep1),rep1, rep2, freq1.1, freq1.2, freq2.1, freq2.2)
  df$av<-rowMeans(df[,2:3])
  df$sd<-apply(df[,2:3],1,sd)
  df$AvgFreq<-rowMeans(df[,4:7])
  df$pos3<-ifelse(substr(df$codon,start = 3,stop=3)=="A"|substr(df$codon,start = 3,stop=3)=="U","A/U","G/C" )
  df<-merge(df, cod_aa, by="codon")
  
  df$codon<-factor(df$codon, levels=df$codon[order(df$av)])
  
  p1<-ggplot(df[df$codon!="UGA"& df$codon!="UAG"&df$codon!="UAA" ,], aes(codon,av))+
    geom_point(aes(size=AvgFreq,colour=pos3))+
    geom_errorbar(aes(ymin=av-sd,ymax=av+sd))+
    theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.2))+
    xlab("")+ylab("mean_codon shift")+
    ggtitle(paste0(condition," vs ",reference," - ",site," site"))+ylim(-2,2)
  
  df$codon<-factor(df$codon, levels=df$codon[order(df$aa)])
  df$codon<-factor(df$codon, levels=df$codon[order(df$codon)])
  
  
  p2<-ggplot(df[df$codon!="UGA"& df$codon!="UAG"&df$codon!="UAA" ,], aes(codon,av))+
    facet_grid(~aa, space="free_x", scales="free_x", switch="x")+
    geom_point(aes(size=AvgFreq,colour=aa))+
    geom_errorbar(aes(ymin=av-sd,ymax=av+sd))+
    scale_colour_discrete(guide = FALSE)+ 
    theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.2))+
    xlab("")+ylab("mean_codon shift")+
    ggtitle(paste0(condition," vs ",reference," - ",site," site"))+ylim(-2,2)
  
  return(list(p1,p2))
  
}


threshold<-50
compare<-function (data, threshold=0, condition, reference, site, trascript_subset)
{
  data<-psite
  threshold<-100
  condition<-"guide1_"
  reference<-"293_"
  site<-"P"
  transcript_subset<-memTrans # from ip experiment, hdlbp_ip.Rscript
  top<-lapply(data, function(x){row.names(x)<-x$transcript; x})
  top<-lapply(top, function(x) x[ , codis])
  
  top<-lapply(top, function(x) x[rowSums(x)>threshold,])
  
  cond1<-top[[paste0(condition,"1")]]
  cond2<-top[[paste0(condition,"2")]]
  ref1<-top[[paste0(reference,"1")]]
  ref2<-top[[paste0(reference,"2")]]
  
  colnames(cond1)<-paste0(colnames(cond1),"_",condition,"1")
  colnames(cond2)<-paste0(colnames(cond2),"_",condition,"2")
  colnames(ref1)<-paste0(colnames(ref1),"_",reference,"1")
  colnames(ref2)<-paste0(colnames(ref2),"_",reference,"2")
  
  reps<-merge(cond1, ref1, by="row.names")
  reps<-merge(reps, cond2, by.x="Row.names", by.y="row.names")
  reps<-merge(reps, ref2, by.x="Row.names", by.y="row.names")
  
  trans_select<-which(reps$Row.names %in% transcript_subset)
  reps<-reps[trans_select,]
  
  normCond1<-subset(reps, select=colnames(reps)[grepl(paste0(condition,"1"),colnames(reps))])
  normCond1<-normCond1/rowSums(normCond1)*1e6
  normCond2<-subset(reps, select=colnames(reps)[grepl(paste0(condition,"2"),colnames(reps))])
  normCond2<-normCond2/rowSums(normCond2)*1e6
  normRef1<-subset(reps, select=colnames(reps)[grepl(paste0(reference,"1"),colnames(reps))])
  normRef1<-normRef1/rowSums(normRef1)*1e6
  normRef2<-subset(reps, select=colnames(reps)[grepl(paste0(reference,"2"),colnames(reps))])
  normRef2<-normRef2/rowSums(normRef2)*1e6
  
  norm<-cbind(normCond1, normRef1, normCond2, normRef2)
  norm$transcript<-reps$Row.names
  
  rep1<-(colMeans(normCond1)-colMeans(normRef1))/colMeans(normRef1)
  rep2<-(colMeans(normCond2)-colMeans(normRef2))/colMeans(normRef2)
  freq1.cond<-colMeans(normCond1)
  freq1.ref<-colMeans(normRef1)
  freq2.cond<-colMeans(normCond2)
  freq2.ref<-colMeans(normRef2)
  
  df<-data.frame(codon=gsub("_.*","",names(rep1)),rep1, rep2, freq1.cond, freq1.ref, freq2.cond, freq2.ref)
  df$av<-rowMeans(df[,2:3])
  df$sd<-apply(df[,2:3],1,sd)
  df$AvgFreq<-rowMeans(df[,4:7])/1e6
  df$pos3<-ifelse(substr(df$codon,start = 3,stop=3)=="A"|substr(df$codon,start = 3,stop=3)=="U","A/U","G/C" )
  df<-merge(df, cod_aa, by="codon")
  
  sigs<-melt(df, measure.vars = colnames(df)[grepl("freq", colnames(df))], id.vars = c("codon"), variable.name = "sample", value.name = "freq")
  sigs$condition<-gsub(".*\\.","",sigs$sample)
  fm1<-summary(lm(freq~codon+condition, sigs))$coefficients[,4]
  fm1<-p.adjust(fm1, method="fdr")[order(p.adjust(fm1, method="fdr"))]
  fm1<-fm1[-length(fm1)]
  
  sigs<-data.frame(codon=rep(df$codon,2), codon_freq=c(df$rep1, df$rep2), aa=rep(df$aa,2))
  fm1<-summary(lm(codon_freq~codon, sigs))$coefficients[,4]
  fm1<-p.adjust(fm1, method="fdr")[order(p.adjust(fm1, method="fdr"))]
  fm1<-data.frame(codon=as.character(gsub("codon","",names(fm1))), padj=fm1)
  fm1$codon<-ifelse(fm1$codon=="(Intercept)", as.character(df$codon[!(df$codon %in% fm1$codon)]), as.character(fm1$codon))
  df<-merge(df, fm1, by="codon")
  
  df$codon<-factor(df$codon, levels=df$codon[order(df$av)])
  df$codon<-factor(df$codon, levels=df$codon[order(df$aa)])
  df$codon<-factor(df$codon, levels=df$codon[order(df$codon)])
  
  ggplot(df[df$codon!="UGA"& df$codon!="UAG"&df$codon!="UAA" ,], aes(codon,av))+
    facet_grid(~aa, space="free_x", scales="free_x", switch="x")+
    geom_point(aes(size=AvgFreq,colour=aa))+
    geom_errorbar(aes(ymin=av-sd,ymax=av+sd))+
    scale_colour_discrete(guide = FALSE)+ 
    theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.2))+
    xlab("")+ylab("mean_codon shift")+
    ggtitle(paste0(condition," vs ",reference," - ",site," site"))+ylim(-.5,.55)
  
  ggplot(df[df$codon!="UGA"& df$codon!="UAG"&df$codon!="UAA" ,], aes(codon,av))+
    facet_grid(~aa, space="free_x", scales="free_x", switch="x")+
    geom_point(aes(size=AvgFreq,colour=padj))+
    geom_errorbar(aes(ymin=av-sd,ymax=av+sd))+
    scale_color_gradient(high="dodgerblue2", low="orange2")+
    theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.2))+
    xlab("")+ylab("mean_codon shift")+
    ggtitle(paste0(condition," vs ",reference," - ",site," site"))+ylim(-.5,.55)
  
  
  return(list(p2))
  
}


mel<-melt(norm, measure.vars = colnames(norm)[!grepl("transcript", colnames(norm))], id.vars="transcript")
mel$codon<-gsub(";.*","",sub("_", ";", mel$variable))
mel$sample<-gsub(".*;","",sub("_", ";", mel$variable))

ggplot(mel, aes(codon, log2(value)))+geom_boxplot()+facet_wrap(~sample)+coord_flip()
ggplot(subset(mel, codon=="CUU"), aes(sample, log2(value)))+geom_boxplot()

con1<-subset(mel, sample=="HDLBP_2")
re1<-subset(mel, sample=="nochx.293_2")
colnames(con1)<-paste0(colnames(con1),"_hdlbp_2")
colnames(re1)<-paste0(colnames(re1),"_nochx.293_2")
mer<-cbind(con1, re1)

ggplot(mer, aes(log2(value_hdlbp_2), log2(value_nochx.293_2)))+geom_point(shape=1, size=0.6)+facet_wrap(~codon_hdlbp_2)+geom_abline(lty=2, colour="grey")

con1<-subset(mel, sample=="HDLBP_1")
re1<-subset(mel, sample=="IgG_1")
colnames(con1)<-paste0(colnames(con1),"_hdlbp")
colnames(re1)<-paste0(colnames(re1),"_IgG_1")
mer<-cbind(con1, re1)

ggplot(mer, aes(log2(value_hdlbp), log2(value_IgG_1)))+geom_point(shape=1, size=0.6)+facet_wrap(~codon_hdlbp)+geom_abline(lty=2, colour="grey")

con1<-subset(mel, sample=="IgG_2")
re1<-subset(mel, sample=="nochx.293_2")
colnames(con1)<-paste0(colnames(con1),"_IgG_2")
colnames(re1)<-paste0(colnames(re1),"_nochx.293_2")
mer<-cbind(con1, re1)

ggplot(mer, aes(log2(value_IgG_2), log2(value_nochx.293_2)))+geom_point(shape=1, size=0.6)+facet_wrap(~codon_IgG_2)+geom_abline(lty=2, colour="grey")




setwd("/Volumes/landthaler/pcp/projects/miha/ulrike_stalling/")

thres<-50

Comp.g1.ctrl.E<-compare(esite, thres, "guide1_", "293_", "E")
Comp.g2.ctrl.E<-compare(esite, thres, "guide2_", "293_", "E")
nochx.Comp.g1.ctrl.E<-compare(esite, thres, "nochx.guide1_", "nochx.293_", "E")
nochx.Comp.g2.ctrl.E<-compare(esite, thres, "nochx.guide2_", "nochx.293_", "E")
ipInp.Comp.E<-compare(esite, thres, "HDLBP_", "nochx.293_", "E")
ipInp.Comp.hdlbp.g1.E<-compare(esite, thres, "HDLBP_", "guide1_", "E")
ipInp.Comp.hdlbp.g2.E<-compare(esite, thres, "HDLBP_", "guide2_", "E")
ipInp.Comp.igg.inp.E<-compare(esite, thres, "IgG_", "nochx.293_", "E")
ipInp.Comp.igg.g1.E<-compare(esite, thres, "IgG_", "guide1_", "E")
ipInp.Comp.igg.g2.E<-compare(esite, thres, "IgG_", "guide2_", "E")
ipInp.Comp.ig.igg.E<-compare(esite, thres, "HDLBP_", "IgG_", "E")

pdf(paste0("all_esite_HDLBP_stalling_RL29-30_threshold",thres,".pdf"),width=10,height = 5)
invisible(lapply(mget(ls()[grepl("Comp.",ls()) ]),print))
dev.off()

rm(list=ls()[grepl("Comp.",ls()) ])


Comp.g1.ctrl.P<-compare(psite, thres, "guide1_", "293_", "P")
Comp.g2.ctrl.P<-compare(psite, thres, "guide2_", "293_", "P")
nochx.Comp.g1.ctrl.P<-compare(psite, thres, "nochx.guide1_", "nochx.293_", "P")
nochx.Comp.g2.ctrl.P<-compare(psite, thres, "nochx.guide2_", "nochx.293_", "P")
ipInp.Comp.P<-compare(psite, thres, "HDLBP_", "nochx.293_", "P")
ipInp.Comp.hdlbp.g1.P<-compare(psite, thres, "HDLBP_", "guide1_", "P")
ipInp.Comp.hdlbp.g2.P<-compare(psite, thres, "HDLBP_", "guide2_", "P")
ipInp.Comp.igg.inp.P<-compare(psite, thres, "IgG_", "nochx.293_", "P")
ipInp.Comp.igg.g1.P<-compare(psite, thres, "IgG_", "guide1_", "P")
ipInp.Comp.igg.g2.P<-compare(psite, thres, "IgG_", "guide2_", "P")
ipInp.Comp.ig.igg.P<-compare(psite, thres, "HDLBP_", "IgG_", "P")

pdf(paste0("all_psite_HDLBP_stalling_RL29-30_threshold",thres,".pdf"),width=10,height = 5)
invisible(lapply(mget(ls()[grepl("Comp.",ls()) ]),print))
dev.off()

rm(list=ls()[grepl("Comp.",ls()) ])


Comp.g1.ctrl.A<-compare(asite, thres, "guide1_", "293_", "A")
Comp.g2.ctrl.A<-compare(asite, thres, "guide2_", "293_", "A")
nochx.Comp.g1.ctrl.A<-compare(asite, thres, "nochx.guide1_", "nochx.293_", "A")
nochx.Comp.g2.ctrl.A<-compare(asite, thres, "nochx.guide2_", "nochx.293_", "A")
ipInp.Comp.A<-compare(asite, thres, "HDLBP_", "nochx.293_", "A")
ipInp.Comp.hdlbp.g1.A<-compare(asite, thres, "HDLBP_", "guide1_", "A")
ipInp.Comp.hdlbp.g2.A<-compare(asite, thres, "HDLBP_", "guide2_", "A")
ipInp.Comp.igg.inp.A<-compare(asite, thres, "IgG_", "nochx.293_", "A")
ipInp.Comp.igg.g1.A<-compare(asite, thres, "IgG_", "guide1_", "A")
ipInp.Comp.igg.g2.A<-compare(asite, thres, "IgG_", "guide2_", "A")
ipInp.Comp.ig.igg.A<-compare(asite, thres, "HDLBP_", "IgG_", "A")

pdf(paste0("all_asite_HDLBP_stalling_RL29-30_threshold",thres,".pdf"),width=10,height = 5)
invisible(lapply(mget(ls()[grepl("Comp.",ls()) ]),print))
dev.off()




##check the renaming of the HDLBP IP samples
names(esite)[7:10]<-c("HDLBP_1", "IgG_1", "HDLBP_2", "IgG_2")

thres<-50

Comp.g1.ctrl.E<-compare(esite, thres, "guide1_", "293_", "E")
Comp.g2.ctrl.E<-compare(esite, thres, "guide2_", "293_", "E")
nochx.Comp.g1.ctrl.E<-compare(esite, thres, "nochx.guide1_", "nochx.293_", "E")
nochx.Comp.g2.ctrl.E<-compare(esite, thres, "nochx.guide2_", "nochx.293_", "E")
ipInp.Comp.E<-compare(esite, thres, "HDLBP_", "nochx.293_", "E")
ipInp.Comp.hdlbp.g1.E<-compare(esite, thres, "HDLBP_", "guide1_", "E")
ipInp.Comp.hdlbp.g2.E<-compare(esite, thres, "HDLBP_", "guide2_", "E")
ipInp.Comp.igg.inp.E<-compare(esite, thres, "IgG_", "nochx.293_", "E")
ipInp.Comp.igg.g1.E<-compare(esite, thres, "IgG_", "guide1_", "E")
ipInp.Comp.igg.g2.E<-compare(esite, thres, "IgG_", "guide2_", "E")
ipInp.Comp.ig.igg.E<-compare(esite, thres, "HDLBP_", "IgG_", "E")

pdf(paste0("renamed_all_esite_HDLBP_stalling_RL29-30_threshold",thres,".pdf"),width=10,height = 5)
invisible(lapply(mget(ls()[grepl("Comp.",ls()) ]),print))
dev.off()


thres<-20
Comp.g1.ctrl.P<-compare(psite, thres, "guide1_", "293_", "P")
Comp.g1.ctrl.A<-compare(asite, thres, "guide1_", "293_", "A")
Comp.g2.ctrl.P<-compare(psite, thres, "guide2_", "293_", "P")
Comp.g2.ctrl.A<-compare(asite, thres, "guide2_", "293_", "A")
# Comp.g1.ctrl.plusOne<-compare(plusOne, thres, "guide1_", "293_", "plusOne")
# Comp.g2.ctrl.plusOne<-compare(plusOne, thres, "guide2_", "293_", "plusOne")


pdf(paste0("mem_HDLBP_stalling_RL29-30_threshold",thres,".pdf"),width=10,height = 5)
invisible(lapply(mget(ls()[grepl("Comp.",ls()) ]),print))
dev.off()




thres<-20
Comp.g1.ctrl.P<-compare(psite, thres, "guide1_", "293_", "P")
Comp.g1.ctrl.A<-compare(asite, thres, "guide1_", "293_", "A")
Comp.g2.ctrl.P<-compare(psite, thres, "guide2_", "293_", "P")
Comp.g2.ctrl.A<-compare(asite, thres, "guide2_", "293_", "A")
# Comp.g1.ctrl.plusOne<-compare(plusOne, thres, "guide1_", "293_", "plusOne")
# Comp.g2.ctrl.plusOne<-compare(plusOne, thres, "guide2_", "293_", "plusOne")

pdf(paste0("cyto_HDLBP_stalling_RL29-30_threshold",thres,".pdf"),width=10,height = 5)
invisible(lapply(mget(ls()[grepl("Comp.",ls()) ]),print))
dev.off()

thres<-20
Comp.g1.ctrl.P<-compare(psite, thres, "guide1_", "293_", "P")
Comp.g1.ctrl.A<-compare(asite, thres, "guide1_", "293_", "A")
Comp.g2.ctrl.P<-compare(psite, thres, "guide2_", "293_", "P")
Comp.g2.ctrl.A<-compare(asite, thres, "guide2_", "293_", "A")
# Comp.g1.ctrl.plusOne<-compare(plusOne, thres, "guide1_", "293_", "plusOne")
# Comp.g2.ctrl.plusOne<-compare(plusOne, thres, "guide2_", "293_", "plusOne")

pdf(paste0("tarCDShi_HDLBP_stalling_RL29-30_threshold",thres,".pdf"),width=10,height = 5)
invisible(lapply(mget(ls()[grepl("Comp.",ls()) ]),print))
dev.off()

thres<-50
Comp.g1.ctrl.P<-compare(psite, thres, "guide1_", "293_", "P")
Comp.g1.ctrl.A<-compare(asite, thres, "guide1_", "293_", "A")
Comp.g2.ctrl.P<-compare(psite, thres, "guide2_", "293_", "P")
Comp.g2.ctrl.A<-compare(asite, thres, "guide2_", "293_", "A")
Comp.g1.ctrl.plusOne<-compare(plusOne, thres, "guide1_", "293_", "plusOne")
Comp.g2.ctrl.plusOne<-compare(plusOne, thres, "guide2_", "293_", "plusOne")
Comp.g1.ctrl.E<-compare(esite, thres, "guide1_", "293_", "E")
Comp.g2.ctrl.E<-compare(esite, thres, "guide2_", "293_", "E")

pdf(paste0("esite_HDLBP_stalling_RL29-30_threshold",thres,".pdf"),width=10,height = 5)
invisible(lapply(mget(ls()[grepl("Comp.",ls()) ]),print))
dev.off()

example_length_dist_zoom <- rlength_distr(reads, sample="293_1", cl=99)
example_length_dist_zoom[["plot"]]
example_ends_heatmap <- rends_heat(reads, anH, sample="293_1", cl=85,
                                   utr5l = 25, cdsl = 40, utr3l = 25)
head(example_ends_heatmap[["df"]])
example_ends_heatmap[["plot"]]


example_ends_heatmap <- rends_heat(reads, anH, sample="guide1_1", cl=85,
                                   utr5l = 25, cdsl = 40, utr3l = 25)
head(example_ends_heatmap[["df"]])
example_ends_heatmap[["plot"]]


example_ends_heatmap <- rends_heat(reads, anH, sample="guide2_1", cl=85,
                                   utr5l = 25, cdsl = 40, utr3l = 25)
head(example_ends_heatmap[["df"]])
example_ends_heatmap[["plot"]]


#this was replaced because of error in original code
# psite_offset <- psite(reads, flanking = 6, extremity="auto")

##select reads by length
# reads<-lapply(reads, subset, length>30 & length<37)
# reads<-lapply(reads, subset, length>28 & length<32)

psite_offset <- psite_new(reads, flanking = 6, extremity="5end")

head(psite_offset, 10)
reads_psite_list <- psite_info(reads, psite_offset)

head(reads_psite_list[["293_1"]])
psite_cds_list <- psite_per_cds(reads_psite_list,anH, start=0, stop=0) #start is important if you want to exclude some reads from analysis
head(psite_cds_list[["293_1"]])
example_frames_stratified <- frame_psite_length(reads_psite_list, sample="293_1",
                                                region="all", cl=90)
head(example_frames_stratified[["df"]])
example_frames_stratified[["plot"]]
example_frames <- frame_psite(reads_psite_list, sample="293_1", region="all")
head(example_frames[["df"]])
example_frames[["plot"]]
example_metaprofile <- metaprofile_psite(reads_psite_list, anH, sample = "293_1",
                                         utr5l = 20, cdsl = 40, utr3l = 20)
example_metaprofile[["plot"]]

example_metaprofile_31 <- metaprofile_psite(reads_psite_list, anH, sample = "293_1",
                                            length_range = 31, utr5l = 20, cdsl = 40,
                                            utr3l = 20)

example_metaprofile_31[["plot"]]

comparison_df <- list()
comparison_df[["subsample_31nt"]] <- subset(reads_psite_list[["293_1"]], length == 31)
comparison_df[["whole_sample"]] <- reads_psite_list[["293_1"]]
names_list <- list("Only_31" = c("subsample_31nt"),
                   "All" = c("whole_sample"))
example_metaheatmap <- metaheatmap_psite(comparison_df, anH, sample = names_list,
                                         utr5l = 20, cdsl = 40, utr3l = 20, log=F)
example_metaheatmap[["plot"]]


#example codon usage
wt_m1 <- codon_usage_psite_new(reads_psite_list, anH, sample = "guide1_1",
                               fastapath="hg19bt1.transcripts.fa") 
wt_m1[["plot"]]

scatter <- codon_usage_psite_new(reads_psite_list, anH, sample = "293_1",
                                 fastapath="hg19bt1.transcripts.fa", codon_usage = wt_m1[["df"]][,c(1,3)])
scatter[["plot_comparison"]]



wt_m1 <- codon_usage_asite_new(reads_psite_list, anH, sample = "guide1_1",
                               fastapath="hg19bt1.transcripts.fa") 
wt_m1[["plot"]]

scatter <- codon_usage_asite_new(reads_psite_list, anH, sample = "293_1",
                                 fastapath="hg19bt1.transcripts.fa", codon_usage = wt_m1[["df"]][,c(1,3)])
scatter[["plot_comparison"]]
