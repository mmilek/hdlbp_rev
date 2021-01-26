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
load("ribowaltz_human.RData")
# setwd("/Volumes/landthaler/pcp/projects/miha/ulrike_stalling/")
setwd("/Volumes/landthaler/pcp/projects/miha/HDLBP/ribo/meta/considered/bed_chx/")
# setwd("/Volumes/landthaler/pcp/projects/miha/HDLBP/ribo/meta/considered/bed_nochx/")
# setwd("E:/work/hdlbp//ribo/meta/considered/bed/")

reads <- bedtolist(bedfolder=getwd(), annotation=anH)
# names(reads)<-gsub("nochx.","",names(reads))
head(reads[["293_1"]],50)

example_length_dist <- rlength_distr(reads, sample = "293_1", cl=99)
example_length_dist[["plot"]]
example_length_dist <- rlength_distr(reads, sample = "guide1_1", cl=99)
example_length_dist[["plot"]]
example_length_dist <- rlength_distr(reads, sample = "guide2_1", cl=99)
example_length_dist[["plot"]]
example_length_dist <- rlength_distr(reads, sample = "293_2", cl=99)
example_length_dist[["plot"]]
example_length_dist <- rlength_distr(reads, sample = "guide1_2", cl=99)
example_length_dist[["plot"]]
example_length_dist <- rlength_distr(reads, sample = "guide2_2", cl=99)
example_length_dist[["plot"]]


example_ends_heatmap <- rends_heat(reads, anH, sample = "293_1", cl = 85,
                                   utr5l = 25, cdsl = 40, utr3l = 25)
example_ends_heatmap[["plot"]]
example_ends_heatmap <- rends_heat(reads, anH, sample = "guide1_1", cl = 85,
                                   utr5l = 25, cdsl = 40, utr3l = 25)
example_ends_heatmap[["plot"]]
example_ends_heatmap <- rends_heat(reads, anH, sample = "guide2_1", cl = 85,
                                   utr5l = 25, cdsl = 40, utr3l = 25)
example_ends_heatmap[["plot"]]
example_ends_heatmap <- rends_heat(reads, anH, sample = "293_2", cl = 85,
                                   utr5l = 25, cdsl = 40, utr3l = 25)
example_ends_heatmap[["plot"]]
example_ends_heatmap <- rends_heat(reads, anH, sample = "guide1_2", cl = 85,
                                   utr5l = 25, cdsl = 40, utr3l = 25)
example_ends_heatmap[["plot"]]
example_ends_heatmap <- rends_heat(reads, anH, sample = "guide2_2", cl = 85,
                                   utr5l = 25, cdsl = 40, utr3l = 25)
example_ends_heatmap[["plot"]]

reads_length_filtered<-lapply(reads, subset, length>28 & length<32)
# rm(reads)
# gc()

psite_offset <- psite_new(reads_length_filtered, flanking = 6, extremity="5end")


#not sure of length filtered should only be considered here, for now yes
# reads_psite_list <- psite_info(reads, psite_offset)
reads_psite_list <- psite_info(reads_length_filtered, psite_offset)

example_frames_stratified <- frame_psite_length(reads_psite_list, sample = "293_1",
                                                region = "all", cl = 90)
example_frames_stratified[["plot"]]



#make psite tracks
setwd("/Volumes/landthaler/pcp/projects/miha/HDLBP/ribo/meta/considered/")
# for (i in names(reads_psite_list) ) {
#   #i<-"293_1"
#   sub<-reads_psite_list[[i]]
#   sub<-subset(sub, !is.na(psite) & !is.na(psite_from_start) & !is.na(psite_region) )
#   sub$count<-1
#   ag<-aggregate(count~psite+transcript, data=sub, sum)
#   ag$start<-ag$psite-1
#   ag<-ag[,c(2,4,1,3)]
#   ag<-ag[order(ag$transcript, ag$start),]
#   write.table(ag, paste0("l29-31.psite.",i,".bedgraph"),quote=F, sep="\t", row.names=F, col.names = F)
# }


#per position
per_pos<-function(x) {
  #x<-reads_psite_list[[1]]
  x<-subset(x, !is.na(psite) & !is.na(psite_from_start) & !is.na(psite_region) )
  x$count<-1
  ag<-aggregate(count~psite+start_pos+stop_pos+transcript, data=x, sum)
  return(ag)
}

positional<-lapply(reads_psite_list, per_pos)


id<-lapply(positional, function(x) paste0(x$transcript,"_",x$psite,"_",x$start_pos,"_",x$stop_pos))
positional<-Map(cbind, positional, id=id)
positional<-lapply(positional, function(x) x[,-c(1,2,3,4)])


merged<-Reduce(function(...) merge(..., by="id", all=T), positional)

colnames(merged)[2:length(colnames(merged))]<-names(positional)
merged$transcript<-gsub("_.*","",merged$id)
merged$psite<-as.numeric(gsub(".*;","",gsub("_.*","",sub("_",";",merged$id))))
merged$start_pos<-as.numeric(gsub(".*;","",gsub("_.*","",sub("_",";",sub("_",";",merged$id)))))
merged$stop_pos<-as.numeric(gsub(".*_","",merged$id))
merged$pos_from_start<-as.numeric(merged$psite)-as.numeric(merged$start_pos)
merged$pos_from_stop<-as.numeric(merged$psite)-as.numeric(merged$stop_pos)
merged$psite_id_start<-paste0(merged$transcript,"_",merged$pos_from_start)
merged$psite_id_stop<-paste0(merged$transcript,"_",merged$pos_from_stop)

# write.table(merged, "psite_position_counts_chx.txt", quote=F, sep="\t", row.names=F, col.names=T)

thr<-5
nrow(subset(merged, `293_1`>=thr | `293_2`>=thr | guide1_1>=thr | guide1_2>=thr | guide2_1>=thr | guide2_2>=thr))

fin<-subset(merged, `293_1`>=thr | `293_2`>=thr | guide1_1>=thr | guide1_2>=thr | guide2_1>=thr | guide2_2>=thr)

#ggplot(fin, aes(log2(`293_1`), log2(`293_2`)))+geom_point()

# write.table(fin, "psite_position_counts_chx.txt", quote=F, sep="\t", row.names=F, col.names=T)
# write.table(fin, "psite_position_counts_nochx.txt", quote=F, sep="\t", row.names=F, col.names=T)
setwd("/Volumes/landthaler/pcp/projects/miha/HDLBP/ribo/meta/considered/")
fin<-read.delim("psite_position_counts_chx.txt", header=T)
# fin<-read.delim("psite_position_counts_nochx.txt", header=T)
colnames(fin)<-gsub("X","", colnames(fin))

normSubset<-subset(fin, pos_from_start >=6 & pos_from_stop<=(-6)) #exclude first and last 2 codons
norm<-colSums(normSubset[,2:7], na.rm=T)

fin<-cbind(fin, fin[,2:7]/norm*1e6)
colnames(fin)[16:ncol(fin)]<-paste0("norm.",colnames(fin)[16:ncol(fin)])

fin$l_cds<-fin$stop_pos-fin$start_pos+1
fin$start<-0
fin$psite_id<-paste0(fin$transcript,"_",fin$pos_from_start,"_",fin$pos_from_stop)
seqs<-fin[,c("transcript", "start","l_cds")]
seqs$id<-paste0(seqs$transcript,"_",seqs$start,"_",seqs$l_cds)
seqs$end<-seqs$l_cds-1
seqs$minend<-(-seqs$end)
seqs<-subset(seqs, !duplicated(id))
lseqStart<-apply(seqs, 1, function(x) seq(x[2], x[5]))
lseqStop<-apply(seqs, 1, function(x) seq(x[6], x[2]))

dseq<-data.frame(transcript=rep(seqs$transcript, seqs$l_cds), pos_from_start=unlist(lseqStart), pos_from_stop=unlist(lseqStop))

dseq$psite_id<-paste0(dseq$transcript,"_",dseq$pos_from_start,"_",dseq$pos_from_stop)

tust<-merge(dseq, fin, by="psite_id", all.x=T)

#here either take normalized or non-normalized
# allpos<-subset(tust, select=c("psite_id_start","transcript.x","pos_from_start.x","293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2"))
allpos<-subset(tust, select=c("psite_id","transcript.x","norm.293_1","norm.293_2" ,"norm.guide1_1" ,"norm.guide1_2" ,"norm.guide2_1", "norm.guide2_2"))
colnames(allpos)[2:ncol(allpos)]<-gsub("norm.","",colnames(allpos)[2:ncol(allpos)])

values<-allpos[,c("293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2")]
values[is.na(values)]<-0
allpos<-cbind(allpos[,1:2],values)
allpos$pos_from_start<-as.numeric(gsub(".*;","",gsub("_.*","",sub("_",";",allpos$psite_id))))
allpos$pos_from_stop<-as.numeric(gsub(".*_","",allpos$psite_id))

allpos$frame_start<-ifelse(allpos$pos_from_start%%3==0, allpos$pos_from_start/3,
                         ifelse(allpos$pos_from_start%%3==1, ((allpos$pos_from_start)-1)/3, ((allpos$pos_from_start)-2)/3))
allpos$frame_stop<-ifelse(allpos$pos_from_stop%%3==0, allpos$pos_from_stop/3,
                           ifelse(allpos$pos_from_stop%%3==1, ((allpos$pos_from_stop)+2)/3, ((allpos$pos_from_stop)+1)/3))
allpos$id<-paste0(allpos$transcript.x,"_",allpos$frame_start,"_",allpos$frame_stop)
# allpos<-allpos[order(allpos$transcript.x,allpos$pos_from_start),]
# allpos<-allpos[order(allpos$transcript.x,allpos$pos_from_stop, decreasing=T),]

##here one could exclude the out-of frame psites, for now we leave them in and sum every psite per codon
codons<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~id, data=allpos, sum)
codons$frame_start<-as.numeric(gsub(".*;","",gsub("_.*","",sub("_",";",codons$id))))
codons$frame_stop<-as.numeric(gsub(".*_","",codons$id))
codons$transcript<-gsub("_.*","",codons$id)

excluded<-which(codons$frame_start!=0 & codons$frame_start!=1 & codons$frame_stop!=0 & codons$frame_stop!=(-1))
excluded<-codons[excluded,]

normi<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~transcript, data=excluded, mean, na.rm=T)
normi[normi==0]<-NA

minorf<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~transcript, data=excluded, sum, na.rm=T)
thrMinOrf<-5
nrow(subset(minorf, `293_1`>=thrMinOrf & `293_2`>=thrMinOrf ))
nrow(subset(minorf, `293_1`>=thrMinOrf & `293_2`>=thrMinOrf & guide1_1>=thrMinOrf & guide1_2>=thrMinOrf & guide2_1>=thrMinOrf & guide2_2>=thrMinOrf ))
minorf<-subset(minorf, `293_1`>=thrMinOrf & `293_2`>=thrMinOrf )

# tes<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~id, data=tust, sum)
# 
# tes$transcript<-gsub("_.*","",tes$id)
# tes$codon<-as.numeric(gsub(".*_","",tes$id))

mask<-merge(codons, normi, by="transcript")
first<-as.matrix(mask[,3:8])
second<-as.matrix(mask[11:16])
nor<-first/second

nor<-cbind(mask[,c(1,2,9,10)], nor)
colnames(nor)<-gsub("\\.x","",colnames(nor))

minFilter<-which(nor$transcript %in% minorf$transcript)
nor<-nor[minFilter,]

# nor<-subset(nor, `293_1`>0 & `293_2`>0 & guide1_1>0 & guide1_2>0 & guide2_1>0 & guide2_2>0) # exclude all zeros for median - this probably not correct, so left out and did mean instead

avg<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~frame_start, data=nor, mean, na.rm=T)
mel<-melt(avg, measure.vars = c("293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2"), id.vars = "frame_start")
ggplot(mel[mel$frame_start>=2 & mel$frame_start<=500,], aes(frame_start, value, colour=variable))+geom_line()+coord_cartesian(xlim=c(0,300))
avg<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~frame_stop, data=nor, mean, na.rm=T)
mel<-melt(avg, measure.vars = c("293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2"), id.vars = "frame_stop")
ggplot(mel[mel$frame_stop<=(-2) & mel$frame_stop>=(-500),], aes(frame_stop, value, colour=variable))+geom_line()+coord_cartesian(xlim=c(-300,0))




#tsig_annot from ulrike_tsig.R
tsig_annot<-read.delim("tsig_annot.txt", header=T)
memS<-subset(tsig_annot, localization_cat=="membrane" & grepl("SignalP", tsig_annot$tsig) )
memT<-subset(tsig_annot, localization_cat=="membrane" & tsig=="TMhelix-only") 
memN<-subset(tsig_annot, localization_cat=="membrane" & tsig=="mem_notsig") 
cyt<-subset(tsig_annot, localization_cat=="cytosolic" & tsig=="cyt_notsig") 

tesmemS<-which(nor$transcript %in% memS$transcript)
tesmemT<-which(nor$transcript %in% memT$transcript)
tesmemN<-which(nor$transcript %in% memN$transcript)
tescyt<-which(nor$transcript %in% cyt$transcript)

tesnmemS<-nor[tesmemS,]
tesnmemT<-nor[tesmemT,]
tesnmemN<-nor[tesmemN,]
tesncyt<-nor[tescyt,]

length(unique(tesnmemS$transcript))
length(unique(tesnmemT$transcript))
length(unique(tesnmemN$transcript))
length(unique(tesncyt$transcript))

tesnmemS$localization<-"memS"
tesnmemT$localization<-"memT"
tesnmemN$localization<-"memN"
tesncyt$localization<-"cyt"

locn<-rbind(tesnmemS,tesnmemT,tesnmemN,tesncyt)

avg<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~frame_start+localization, data=locn, mean, na.rm=T)
mel<-melt(avg, measure.vars = c("293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2"), id.vars = c("frame_start","localization"))

library(zoo)

ggplot(mel[mel$frame_start>=2 & mel$frame_start<=500 & mel$localization!="memN" ,], aes(frame_start, value, colour=localization))+geom_line(aes(y=rollmean(value, 10, na.pad=T)))+facet_wrap(~variable)+coord_cartesian(xlim=c(0,250))+geom_vline(xintercept=39, lty=2)+geom_vline(xintercept=79, lty=2)+ylab("scaled_psite_coverage")
ggplot(mel[mel$frame_start>=2 & mel$frame_start<=500 & mel$localization!="memN",], aes(frame_start, value, colour=localization))+geom_line(aes(y=rollmean(value, 10, na.pad=T)))+facet_wrap(~variable)+coord_cartesian(xlim=c(0,500))+geom_vline(xintercept=39, lty=2)+geom_vline(xintercept=79, lty=2)+ylab("scaled_psite_coverage")


avg<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~frame_stop+localization, data=locn, mean, na.rm=T)
mel<-melt(avg, measure.vars = c("293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2"), id.vars = c("frame_stop","localization"))

library(zoo)

ggplot(mel[mel$frame_stop<=(-2) & mel$frame_stop>=(-500) & mel$localization!="memN" ,], aes(frame_stop, value, colour=localization))+geom_line(aes(y=rollmean(value, 10, na.pad=T)))+facet_wrap(~variable)+coord_cartesian(xlim=c(-250,0))+ylab("scaled_psite_coverage")
ggplot(mel[mel$frame_stop<=(-2) & mel$frame_stop>=(-500) & mel$localization!="memN" ,], aes(frame_stop, value, colour=localization))+geom_line(aes(y=rollmean(value, 10, na.pad=T)))+facet_wrap(~variable)+coord_cartesian(xlim=c(-500,0))+ylab("scaled_psite_coverage")



#per gene
selected_gene<-"ATP2A2"
trans_id<-which(tsig_annot$Symbol==selected_gene)
trans_id<-tsig_annot[trans_id,]
trans_row<-which(codons$transcript %in% trans_id$transcript)
gene<-codons[trans_row,]
gene<-gene[order(gene$frame_start),]

melito<-melt(gene, id.vars = "frame_start", measure.vars = colnames(gene)[2:7], variable.name = "sample", value.name = "psite_rpm")
ggplot(melito, aes(frame_start, psite_rpm))+geom_line(aes(y=rollmean(psite_rpm, 11, na.pad=T)))+ggtitle(selected_gene)+geom_vline(xintercept = trans_id$start-1,lty=2)+geom_vline(xintercept = trans_id$start-1+40,lty=2)+facet_wrap(~sample)+coord_cartesian(xlim=c(0,trans_id$l_cds/3))


selected_gene<-"DDOST"
trans_id<-which(tsig_annot$Symbol==selected_gene)
trans_id<-tsig_annot[trans_id,]
trans_row<-which(nor$transcript %in% trans_id$transcript)
gene<-nor[trans_row,]
gene<-gene[order(gene$frame_start),]

melito<-melt(gene, id.vars = "frame_start", measure.vars = colnames(gene)[5:10], variable.name = "sample", value.name = "scaled_psite_rpm")
ggplot(melito, aes(frame_start, scaled_psite_rpm))+geom_line(aes(y=rollmean(scaled_psite_rpm, 11, na.pad=T)))+ggtitle(selected_gene)+geom_vline(xintercept = trans_id$start-1,lty=2)+geom_vline(xintercept = trans_id$start-1+40,lty=2)+facet_wrap(~sample)+coord_cartesian(xlim=c(0,trans_id$l_cds/3))


tsig_annot<-tsig_annot[order(tsig_annot$tc_CDS_norm, decreasing = T),]
head(subset(tsig_annot, tpm_cutoff>=200 & loc_tar_CDS =="membrane_tc>5.56 & tc<65.26"),20)




##more cytosolic classes
cytHi<-subset(tsig_annot, localization_cat=="cytosolic" & tsig=="cyt_notsig" & loc_tar_CDS=="cytosolic_tc>0.47 & tc<27.84") 
cytMid<-subset(tsig_annot, localization_cat=="cytosolic" & tsig=="cyt_notsig" & loc_tar_CDS=="cytosolic_tc>0.18 & tc<0.47")
cytLo<-subset(tsig_annot, localization_cat=="cytosolic" & tsig=="cyt_notsig" & loc_tar_CDS=="cytosolic_tc<0.18")


tesmemS<-which(nor$transcript %in% cytHi$transcript)
tesmemT<-which(nor$transcript %in% cytMid$transcript)
tesmemN<-which(nor$transcript %in% cytLo$transcript)

tesnmemS<-nor[tesmemS,]
tesnmemT<-nor[tesmemT,]
tesnmemN<-nor[tesmemN,]

length(unique(tesnmemS$transcript))
length(unique(tesnmemT$transcript))
length(unique(tesnmemN$transcript))

tesnmemS$localization<-"cytHi"
tesnmemT$localization<-"cytMid"
tesnmemN$localization<-"cytLo"

locn<-rbind(tesnmemS,tesnmemT,tesnmemN)

avg<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~frame_start+localization, data=locn, mean, na.rm=T)
mel<-melt(avg, measure.vars = c("293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2"), id.vars = c("frame_start","localization"))

library(zoo)

ggplot(mel[mel$frame_start>=2 & mel$frame_start<=500 & mel$localization!="cytLo",], aes(frame_start, value, colour=localization))+geom_line(aes(y=rollmean(value, 10, na.pad=T)))+facet_wrap(~variable)+coord_cartesian(xlim=c(0,250))+geom_vline(xintercept=39, lty=2)+geom_vline(xintercept=79, lty=2)+ylab("scaled_psite_coverage")
ggplot(mel[mel$frame_start>=2 & mel$frame_start<=500 & mel$localization!="cytLo",], aes(frame_start, value, colour=localization))+geom_line(aes(y=rollmean(value, 10, na.pad=T)))+facet_wrap(~variable)+coord_cartesian(xlim=c(0,500))+geom_vline(xintercept=39, lty=2)+geom_vline(xintercept=79, lty=2)+ylab("scaled_psite_coverage")


avg<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~frame_stop+localization, data=locn, mean, na.rm=T)
mel<-melt(avg, measure.vars = c("293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2"), id.vars = c("frame_stop","localization"))

library(zoo)

ggplot(mel[mel$frame_stop<=(-2) & mel$frame_stop>=(-500) & mel$localization!="cytLo" ,], aes(frame_stop, value, colour=localization))+geom_line(aes(y=rollmean(value, 10, na.pad=T)))+facet_wrap(~variable)+coord_cartesian(xlim=c(-250,0))+ylab("scaled_psite_coverage")
ggplot(mel[mel$frame_stop<=(-2) & mel$frame_stop>=(-500) & mel$localization!="cytLo" ,], aes(frame_stop, value, colour=localization))+geom_line(aes(y=rollmean(value, 10, na.pad=T)))+facet_wrap(~variable)+coord_cartesian(xlim=c(-500,0))+ylab("scaled_psite_coverage")


##more membrane classes with signalP
memHi<-subset(tsig_annot, localization_cat=="membrane" & grepl("SignalP", tsig_annot$tsig) & loc_tar_CDS=="membrane_tc>5.56 & tc<65.26") 
memMid<-subset(tsig_annot, localization_cat=="membrane" & grepl("SignalP", tsig_annot$tsig) & loc_tar_CDS=="membrane_tc>1.66 & tc<5.56") 
memLo<-subset(tsig_annot, localization_cat=="membrane" & grepl("SignalP", tsig_annot$tsig) & loc_tar_CDS=="membrane_tc<1.66") 


tesmemS<-which(nor$transcript %in% memHi$transcript)
tesmemT<-which(nor$transcript %in% memMid$transcript)
tesmemN<-which(nor$transcript %in% memLo$transcript)

tesnmemS<-nor[tesmemS,]
tesnmemT<-nor[tesmemT,]
tesnmemN<-nor[tesmemN,]

length(unique(tesnmemS$transcript))
length(unique(tesnmemT$transcript))
length(unique(tesnmemN$transcript))

tesnmemS$localization<-"memHi"
tesnmemT$localization<-"memMid"
tesnmemN$localization<-"memLo"

locn<-rbind(tesnmemS,tesnmemT,tesnmemN)

avg<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~frame_start+localization, data=locn, mean, na.rm=T)
mel<-melt(avg, measure.vars = c("293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2"), id.vars = c("frame_start","localization"))

library(zoo)

ggplot(mel[mel$frame_start>=2 & mel$frame_start<=500 & mel$localization!="memLo",], aes(frame_start, value, colour=localization))+geom_line(aes(y=rollmean(value, 10, na.pad=T)))+facet_wrap(~variable)+coord_cartesian(xlim=c(0,250))+geom_vline(xintercept=39, lty=2)+geom_vline(xintercept=79, lty=2)+ylab("scaled_psite_coverage")
ggplot(mel[mel$frame_start>=2 & mel$frame_start<=500 & mel$localization!="memLo",], aes(frame_start, value, colour=localization))+geom_line(aes(y=rollmean(value, 10, na.pad=T)))+facet_wrap(~variable)+coord_cartesian(xlim=c(0,500))+geom_vline(xintercept=39, lty=2)+geom_vline(xintercept=79, lty=2)+ylab("scaled_psite_coverage")


avg<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~frame_stop+localization, data=locn, mean, na.rm=T)
mel<-melt(avg, measure.vars = c("293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2"), id.vars = c("frame_stop","localization"))

library(zoo)

ggplot(mel[mel$frame_stop<=(-2) & mel$frame_stop>=(-500) & mel$localization!="memLo" ,], aes(frame_stop, value, colour=localization))+geom_line(aes(y=rollmean(value, 10, na.pad=T)))+facet_wrap(~variable)+coord_cartesian(xlim=c(-250,0))+ylab("scaled_psite_coverage")
ggplot(mel[mel$frame_stop<=(-2) & mel$frame_stop>=(-500) & mel$localization!="memLo" ,], aes(frame_stop, value, colour=localization))+geom_line(aes(y=rollmean(value, 10, na.pad=T)))+facet_wrap(~variable)+coord_cartesian(xlim=c(-500,0))+ylab("scaled_psite_coverage")



##done till here





#older version
dummy<- data.frame(transcript=rep(unique(fin$transcript), each=length(seq(-90,1502))), 
                   pos_from_start=rep(seq(-90,1502), length(unique(fin$transcript))),
                   frame_start=rep(rep(seq(-30,500), each=3),  length(unique(fin$transcript))),
                   pos_from_stop=rep(seq(-1502,90), length(unique(fin$transcript))),
                   frame_stop=rep(rep(seq(-500,30), each=3),  length(unique(fin$transcript))))

dummy$psite_id_start<-paste0(dummy$transcript,"_",dummy$pos_from_start)
dummy$psite_id_stop<-paste0(dummy$transcript,"_",dummy$pos_from_stop)

dummy<-merge(dummy, anH, by="transcript")
dummyMin<-subset(dummy, pos_from_start<0)
dummyPos<-subset(dummy, pos_from_start>=0)
dummyMin<-subset(dummyMin, pos_from_start>=-l_utr5 )
dummyPos<-subset(dummyPos, pos_from_start<=l_cds )

dummy<-rbind(dummyMin,dummyPos)

tust<-merge(dummy, fin, by="psite_id_start", all.x=T)

tust<-subset(tust, select=c("psite_id_start","transcript.x","pos_from_start.x","frame_start","293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2"))
values<-tust[,c("293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2")]
values[is.na(values)]<-0
tust<-cbind(tust[,1:4],values)
tust$id<-paste0(tust$transcript.x,"_",tust$frame_start)

codons<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~id+transcript.x, data=tust, sum)
codons$frame_start<-as.numeric(gsub(".*_","",codons$id))

excluded<-which(codons$frame_start!=0 & codons$frame_start!=1)
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

ggplot(mel, aes(codon, value, colour=variable))+geom_line()+coord_cartesian(xlim=c(-30,500) ,ylim=c(0,10))
ggplot(mel[mel$codon>2,], aes(codon, value, colour=variable))+geom_line()+coord_cartesian(xlim=c(-30,500 ))


setwd("~/Google Drive/hdlbp/")

mas<-read.delim("hdlbp_master_table_with_classes.txt", header=T)
inf<-subset(mas, select=c("gene_id","Symbol",
                          "tpm_cutoff","tc_CDS_norm","localization_cat",
                          "mean_te_293","loc_tar_CDS","tc_CDS_norm_cat",
                          "tc_transcript_norm",
                          "log2FoldChange.mem.cyt.293",
                          "log2FoldChange.mem.cyt.KO.293",
                          "log2FoldChange.ribo.rna.KO.WT",
                          "Transmembrane.helices" ,"Cleavage.site..Signalp." ))


setwd("/Volumes/landthaler/pcp/projects/miha/HDLBP/ribo/meta/considered/")

trans<-read.delim("considered_gene_names.txt", header=T)[,c(1:10,16,18)]
colnames(trans)[10]<-"gene_id"

inf<-merge(inf, trans, by="gene_id")

memS<-subset(inf, localization_cat=="membrane" & Cleavage.site..Signalp.=="SignalP" & Transmembrane.helices=="noTM")
memT<-subset(inf, localization_cat=="membrane" & Cleavage.site..Signalp.=="noSignalP" & Transmembrane.helices=="TMhelix")
memN<-subset(inf, localization_cat=="membrane" & Cleavage.site..Signalp.=="noSignalP" & Transmembrane.helices=="noTM")
cyt<-subset(inf, localization_cat=="cytosolic")

tesmemS<-which(nor$transcript %in% memS$transcript)
tesmemT<-which(nor$transcript %in% memT$transcript)
tesmemN<-which(nor$transcript %in% memN$transcript)
tescyt<-which(nor$transcript %in% cyt$transcript)

tesnmemS<-nor[tesmemS,]
tesnmemT<-nor[tesmemT,]
tesnmemN<-nor[tesmemN,]
tesncyt<-nor[tescyt,]

length(unique(tesnmemS$transcript))
length(unique(tesnmemT$transcript))
length(unique(tesnmemN$transcript))
length(unique(tesncyt$transcript))

tesnmemS$localization<-"memS"
tesnmemT$localization<-"memT"
tesnmemN$localization<-"memN"
tesncyt$localization<-"cyt"

locn<-rbind(tesnmemS,tesnmemT,tesnmemN,tesncyt)

# locn<-subset(locn, rowMeans(locn[,c("293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2")])>30)

avg<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~codon+localization, data=locn, mean, na.rm=T)
mel<-melt(avg, measure.vars = c("293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2"), id.vars = c("codon","localization"))

ggplot(mel[ (mel$localization=="memS" | mel$localization=="cyt"),], aes(codon, value, colour=localization))+geom_line()+coord_cartesian(xlim=c(-30,300 ))+facet_wrap(~variable)
ggplot(mel[mel$codon!=0 & mel$codon!=1 & mel$codon!=2 & mel$codon!=3 & mel$codon!=4,], aes(codon, value, colour=localization))+geom_line()+coord_cartesian(xlim=c(-5,20 ))+facet_wrap(~variable)
ggplot(mel[mel$codon!=0 & mel$codon!=1 ,], aes(codon, value, colour=localization))+geom_line()+coord_cartesian(xlim=c(2,500 ))+facet_wrap(~variable)

ggplot(mel[mel$codon!=0 & mel$codon!=1,], aes(codon, value, colour=localization))+geom_smooth(method="loess")+coord_cartesian(xlim=c(-30,400 ))+facet_wrap(~variable)
ggplot(mel[mel$codon>=5,], aes(codon, value, colour=localization))+geom_smooth(method="loess")+coord_cartesian(xlim=c(-30,300 ))+facet_wrap(~variable)

ggplot(mel[ (mel$localization=="memS" | mel$localization=="cyt"),], aes(codon, value, colour=localization))+geom_smooth()+coord_cartesian(xlim=c(-30,300 ))+facet_wrap(~variable)

ggplot(mel[,], aes(codon, value, colour=localization))+geom_smooth()+coord_cartesian(xlim=c(-30,300 ))+facet_wrap(~variable)




##targets
memS<-subset(inf, loc_tar_CDS=="membrane_tc>10.82 & tc<181.46" & tpm_cutoff>=20 )
memT<-subset(inf, loc_tar_CDS== "membrane_tc>2.81 & tc<10.82"& tpm_cutoff>=20)
memN<-subset(inf, loc_tar_CDS=="membrane_tc<2.81"& tpm_cutoff>=20)
cytS<-subset(inf, loc_tar_CDS=="cytosolic_tc>0.83 & tc<166.4"& tpm_cutoff>=20)
cytT<-subset(inf, loc_tar_CDS=="cytosolic_tc>0.3 & tc<0.83"& tpm_cutoff>=20)
cytN<-subset(inf, loc_tar_CDS=="cytosolic_tc<0.3"& tpm_cutoff>=20)

tesmemS<-which(nor$transcript %in% memS$transcript)
tesmemT<-which(nor$transcript %in% memT$transcript)
tesmemN<-which(nor$transcript %in% memN$transcript)
tescytS<-which(nor$transcript %in% cytS$transcript)
tescytT<-which(nor$transcript %in% cytT$transcript)
tescytN<-which(nor$transcript %in% cytN$transcript)

tesnmemS<-nor[tesmemS,]
tesnmemT<-nor[tesmemT,]
tesnmemN<-nor[tesmemN,]
tesncytS<-nor[tescytS,]
tesncytT<-nor[tescytT,]
tesncytN<-nor[tescytN,]


length(unique(tesnmemS$transcript))
length(unique(tesnmemT$transcript))
length(unique(tesnmemN$transcript))
length(unique(tesncytS$transcript))
length(unique(tesncytT$transcript))
length(unique(tesncytN$transcript))


tesnmemS$localization<-"memS"
tesnmemT$localization<-"memT"
tesnmemN$localization<-"memN"
tesncytS$localization<-"cytS"
tesncytT$localization<-"cytT"
tesncytN$localization<-"cytN"


locn<-rbind(tesnmemS,tesnmemT,tesnmemN,tesncytS, tesncytT, tesncytN)


avg<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~codon+localization, data=locn, mean, na.rm=T)
mel<-melt(avg, measure.vars = c("293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2"), id.vars = c("codon","localization"))

ggplot(mel[ mel$codon>=2,], aes(codon, value, colour=localization))+geom_smooth()+coord_cartesian(xlim=c(-30,500 ))+facet_wrap(~variable)

ggplot(mel[ mel$codon>=2 & (mel$localization=="memS"| mel$localization=="memN" | mel$localization=="memT" | mel$localization=="cytN"),], aes(codon, value, colour=localization))+geom_smooth()+coord_cartesian(xlim=c(-30,500 ))+facet_wrap(~variable)

ggplot(mel[ mel$localization=="memS" & mel$codon>=2,], aes(codon, value, colour=variable))+geom_smooth()+coord_cartesian(xlim=c(-30,500 ))

#rolling mean
library(zoo)
ggplot(mel[ mel$codon>=2 &(mel$localization=="memS"| mel$localization=="memN" | mel$localization=="memT" | mel$localization=="cytN") ,], aes(codon, value, colour=localization))+
  geom_line(aes(y=rollmean(value, 50, na.pad=T)))+coord_cartesian(xlim=c(-30,500))+facet_wrap(~variable)

ggplot(mel[ mel$localization=="memS" & mel$codon>=2,], aes(codon, value, colour=variable))+geom_line(aes(y=rollmean(value, 50, na.pad=T)))+coord_cartesian(xlim=c(-30,500 ))



# from the 3end

dummy<- data.frame(transcript=rep(unique(fin$transcript), each=length(seq(-90,1502))), 
                   pos_from_start=rep(seq(-90,1502), length(unique(fin$transcript))),
                   frame_start=rep(rep(seq(-30,500), each=3),  length(unique(fin$transcript))),
                   pos_from_stop=rep(seq(-1502,90), length(unique(fin$transcript))),
                   frame_stop=rep(rep(seq(-500,30), each=3),  length(unique(fin$transcript))))

dummy$psite_id_start<-paste0(dummy$transcript,"_",dummy$pos_from_start)
dummy$psite_id_stop<-paste0(dummy$transcript,"_",dummy$pos_from_stop)


dummy<-merge(dummy, anH, by="transcript")
dummyMin<-subset(dummy, pos_from_stop<0)
dummyPos<-subset(dummy, pos_from_stop>=0)
dummyMin<-subset(dummyMin, pos_from_stop>=-l_cds )
dummyPos<-subset(dummyPos, pos_from_stop<=l_utr3 )

dummy<-rbind(dummyMin,dummyPos)

tust<-merge(dummy, fin, by="psite_id_stop", all.x=T)

tust<-subset(tust, select=c("psite_id_stop","transcript.x","pos_from_stop.x","frame_stop","293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2"))
values<-tust[,c("293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2")]
values[is.na(values)]<-0
tust<-cbind(tust[,1:4],values)
tust$id<-paste0(tust$transcript.x,"_",tust$frame_stop)

codons<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~id+transcript.x, data=tust, sum)
codons$frame_stop<-as.numeric(gsub(".*_","",codons$id))

excluded<-which(codons$frame_stop!=0 & codons$frame_stop!=(-1))
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


memS<-subset(inf, loc_tar_CDS=="membrane_tc>10.82 & tc<181.46" & tpm_cutoff>=10)
memT<-subset(inf, loc_tar_CDS== "membrane_tc>2.81 & tc<10.82"& tpm_cutoff>=10)
memN<-subset(inf, loc_tar_CDS=="membrane_tc<2.81"& tpm_cutoff>=10)
cytS<-subset(inf, loc_tar_CDS=="cytosolic_tc>0.83 & tc<166.4"& tpm_cutoff>=10)
cytT<-subset(inf, loc_tar_CDS=="cytosolic_tc>0.3 & tc<0.83"& tpm_cutoff>=10)
cytN<-subset(inf, loc_tar_CDS=="cytosolic_tc<0.3"& tpm_cutoff>=10)

tesmemS<-which(nor$transcript %in% memS$transcript)
tesmemT<-which(nor$transcript %in% memT$transcript)
tesmemN<-which(nor$transcript %in% memN$transcript)
tescytS<-which(nor$transcript %in% cytS$transcript)
tescytT<-which(nor$transcript %in% cytT$transcript)
tescytN<-which(nor$transcript %in% cytN$transcript)

tesnmemS<-nor[tesmemS,]
tesnmemT<-nor[tesmemT,]
tesnmemN<-nor[tesmemN,]
tesncytS<-nor[tescytS,]
tesncytT<-nor[tescytT,]
tesncytN<-nor[tescytN,]


length(unique(tesnmemS$transcript))
length(unique(tesnmemT$transcript))
length(unique(tesnmemN$transcript))
length(unique(tesncytS$transcript))
length(unique(tesncytT$transcript))
length(unique(tesncytN$transcript))


tesnmemS$localization<-"memS"
tesnmemT$localization<-"memT"
tesnmemN$localization<-"memN"
tesncytS$localization<-"cytS"
tesncytT$localization<-"cytT"
tesncytN$localization<-"cytN"


locn<-rbind(tesnmemS,tesnmemT,tesnmemN,tesncytS, tesncytT, tesncytN)


avg<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~codon+localization, data=locn, mean, na.rm=T)
mel<-melt(avg, measure.vars = c("293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2"), id.vars = c("codon","localization"))

ggplot(mel, aes(codon, value, colour=localization))+geom_smooth()+facet_wrap(~variable)
ggplot(mel[mel$codon<=0,], aes(codon, value, colour=localization))+geom_smooth()+facet_wrap(~variable)

ggplot(mel[mel$codon>=2,], aes(codon, value, colour=localization))+geom_smooth()+facet_wrap(~variable)


ggplot(mel[  (mel$localization=="memS"| mel$localization=="memN" | mel$localization=="memT" | mel$localization=="cytN"),], aes(codon, value, colour=localization))+geom_smooth()+facet_wrap(~variable)

ggplot(mel[ mel$localization=="memS" & mel$codon>=2,], aes(codon, value, colour=variable))+geom_smooth()+coord_cartesian(xlim=c(-30,500 ))










