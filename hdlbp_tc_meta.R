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
setwd("E:/Google Drive/koshi_revision/riboprof/")
load("ribowaltz_human.RData")

# setwd("/Volumes/landthaler/pcp/projects/miha/HDLBP/all_clip_data/reclip/mapping_trans/")
setwd("D:/landthaler/HDLBP/all_clip_data/reclip/mapping_trans/")

tc<-read.delim("reproducible.hdlbp.TCseq.bed", header=F)
colnames(tc)<-c("transcript_id","tc_start","tc_stop","tc_num1","all_reads1","tc_num2","all_reads2","seq")

# setwd("/Volumes/landthaler/pcp/projects/miha/HDLBP/ribo/meta/considered/")
setwd("D:/landthaler/HDLBP/ribo/meta/considered/")

tsig_annot<-read.delim("tsig_annot.txt", header=T)


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



#with imputation
regionCounts$utr3_vs_cds1<-ifelse(!is.na(regionCounts$utr3.norm.tc_num1) & is.na(regionCounts$cds.norm.tc_num1), 
                                  regionCounts$utr3.norm.tc_num1/min(regionCounts[!is.na(regionCounts$cds.norm.tc_num1), "cds.norm.tc_num1"]),
                                  ifelse(is.na(regionCounts$utr3.norm.tc_num1) & !is.na(regionCounts$cds.norm.tc_num1), 
                                         min(regionCounts[!is.na(regionCounts$utr3.norm.tc_num1), "utr3.norm.tc_num1"])/regionCounts$cds.norm.tc_num1,
                                  ifelse(!is.na(regionCounts$utr3.norm.tc_num1) & !is.na(regionCounts$cds.norm.tc_num1),
                                         regionCounts$utr3.norm.tc_num1/regionCounts$cds.norm.tc_num1, NA)))
regionCounts$utr3_vs_cds2<-ifelse(!is.na(regionCounts$utr3.norm.tc_num2) & is.na(regionCounts$cds.norm.tc_num2), 
                                  regionCounts$utr3.norm.tc_num2/min(regionCounts[!is.na(regionCounts$cds.norm.tc_num2), "cds.norm.tc_num2"]),
                                  ifelse(is.na(regionCounts$utr3.norm.tc_num2) & !is.na(regionCounts$cds.norm.tc_num2), 
                                         min(regionCounts[!is.na(regionCounts$utr3.norm.tc_num2), "utr3.norm.tc_num2"])/regionCounts$cds.norm.tc_num2,
                                         ifelse(!is.na(regionCounts$utr3.norm.tc_num2) & !is.na(regionCounts$cds.norm.tc_num2),
                                                regionCounts$utr3.norm.tc_num2/regionCounts$cds.norm.tc_num2, NA)))

corrplot(cor(regionCounts[2:ncol(regionCounts)], use = "pairwise.complete.obs"), type="upper", method="color")

ggplot(regionCounts, aes(log2(utr3_vs_cds1), log2(utr3_vs_cds2)))+geom_point()


regionCounts<-merge(regionCounts, tsig_annot, by="gene_id")
regionCounts$cds.norm.tc_num1.len<-regionCounts$cds.norm.tc_num1/regionCounts$l_cds*1e3
regionCounts$cds.norm.tc_num2.len<-regionCounts$cds.norm.tc_num2/regionCounts$l_cds*1e3
regionCounts$utr3.norm.tc_num1.len<-regionCounts$utr3.norm.tc_num1/regionCounts$l_utr3*1e3
regionCounts$utr3.norm.tc_num2.len<-regionCounts$utr3.norm.tc_num2/regionCounts$l_utr3*1e3
regionCounts$utr5.norm.tc_num1.len<-regionCounts$utr5.norm.tc_num1/regionCounts$l_utr5*1e3
regionCounts$utr5.norm.tc_num2.len<-regionCounts$utr5.norm.tc_num2/regionCounts$l_utr5*1e3
regionCounts$trans.norm.tc_num1.len<-regionCounts$trans.norm.tc_num1/regionCounts$l_tr*1e3
regionCounts$trans.norm.tc_num2.len<-regionCounts$trans.norm.tc_num2/regionCounts$l_tr*1e3

ggplot(regionCounts,aes(log2(cds.norm.tc_num1.len),log2(cds.norm.tc_num2.len), colour=localization_cat))+geom_point()
ggplot(regionCounts,aes(log2(utr3.norm.tc_num1.len),log2(utr3.norm.tc_num2.len), colour=localization_cat))+geom_point()

ggplot(regionCounts,aes(log2(cds.norm.tc_num1.len),log2(utr3.norm.tc_num1.len), colour=localization_cat))+geom_point()
ggplot(regionCounts,aes(log2(cds.norm.tc_num2.len),log2(utr3.norm.tc_num2.len), colour=localization_cat))+geom_point()
ggplot(regionCounts,aes(log2(cds.norm.tc_num2.len),log2(trans.norm.tc_num2.len), colour=localization_cat))+geom_point()
ggplot(regionCounts,aes(log2(cds.norm.tc_num1.len),log2(trans.norm.tc_num1.len), colour=localization_cat))+geom_point()
ggplot(regionCounts,aes(log2(utr3.norm.tc_num1.len),log2(trans.norm.tc_num1.len), colour=localization_cat))+geom_point()
ggplot(regionCounts,aes(log2(utr3.norm.tc_num2.len),log2(trans.norm.tc_num2.len), colour=localization_cat))+geom_point()

ggplot(regionCounts, aes(log2(utr3.norm.tc_num1.len/cds.norm.tc_num1.len), colour=localization_cat))+stat_ecdf()
ggplot(regionCounts, aes(log2(utr3.norm.tc_num2.len/cds.norm.tc_num2.len), colour=localization_cat))+stat_ecdf()

ggplot(regionCounts, aes(log2(utr3.norm.tc_num2.len/cds.norm.tc_num2.len), log2(trans.norm.tc_num1.len), colour=localization_cat))+geom_point()
ggplot(regionCounts[regionCounts$localization_cat=="membrane",], aes(log2(utr3.norm.tc_num2.len/cds.norm.tc_num2.len), log2(trans.norm.tc_num1.len) ))+geom_point()+geom_point()+coord_cartesian(xlim=c(-10,8))
ggplot(regionCounts[regionCounts$localization_cat=="cytosolic",], aes(log2(utr3.norm.tc_num2.len/cds.norm.tc_num2.len), log2(trans.norm.tc_num1.len), colour=localization_cat))+geom_point()+coord_cartesian(xlim=c(-10,8))


ggplot(regionCounts, aes(log2(utr3_vs_cds1), log2(tc_transcript_norm), colour=localization_cat))+geom_point()

regionCounts$utr3_vs_cds1.norm<-regionCounts$utr3.norm.tc_num1.len/regionCounts$cds.norm.tc_num1.len
regionCounts$utr3_vs_cds2.norm<-regionCounts$utr3.norm.tc_num2.len/regionCounts$cds.norm.tc_num2.len


#get zeros for positional information
seqs<-dat[,c("transcript_id", "l_tr","l_utr5", "l_cds","l_utr3")]
seqs$id<-paste0(seqs$transcript,"_",seqs$l_tr,"_",seqs$l_utr5, "_", seqs$l_cds,"_", seqs$l_utr3)
seqs<-subset(seqs, !duplicated(id))
seqs$start<-1
seqs$min_lutr5<--(seqs$l_utr5)
seqs$stop<--1
seqs$min_lcds<--(seqs$l_cds)
seqs$min_lutr3<--(seqs$l_utr3)
seqs<-subset(seqs, l_utr5!=0)
seqs<-subset(seqs, l_utr3!=0)
seqs<-subset(seqs, l_cds!=0)

lseqStartUtr5<-apply(seqs, 1, function(x) seq(x[7], x[3]))
lseqStopUtr5<-apply(seqs, 1, function(x) seq(x[8], x[9]))
lseqStartCds<-apply(seqs, 1, function(x) seq(x[7], x[4]))
lseqStopCds<-apply(seqs, 1, function(x) seq(x[10], x[9]))
lseqStartUtr3<-apply(seqs, 1, function(x) seq(x[7], x[5]))
lseqStopUtr3<-apply(seqs, 1, function(x) seq(x[11], x[9]))

utr5<-data.frame(transcript=rep(seqs$transcript, seqs$l_utr5), pos_from_start=unlist(lseqStartUtr5), pos_from_stop=unlist(lseqStopUtr5))
cds<-data.frame(transcript=rep(seqs$transcript, seqs$l_cds), pos_from_start=unlist(lseqStartCds), pos_from_stop=unlist(lseqStopCds))
utr3<-data.frame(transcript=rep(seqs$transcript, seqs$l_utr3), pos_from_start=unlist(lseqStartUtr3), pos_from_stop=unlist(lseqStopUtr3))

utr5$tc_id<-paste0(utr5$transcript,"_", utr5$pos_from_start,"_", utr5$pos_from_stop,"_utr5")
cds$tc_id<-paste0(cds$transcript,"_", cds$pos_from_start,"_", cds$pos_from_stop,"_cds")
utr3$tc_id<-paste0(utr3$transcript,"_", utr3$pos_from_start,"_", utr3$pos_from_stop,"_utr3")

dseq<-rbind(utr5,cds,utr3)[,c(1,4)]
rm(utr5, cds, utr3)

newdat<-subset(dat, select=c("tc_id","tc_num1", "tc_num2", "norm.tc_num1", "norm.tc_num2"))
tl<-!(dseq$tc_id %in% newdat$tc_id)
dseq<-dseq[tl,]
dseq<-data.frame(tc_id=dseq$tc_id, tc_num1=rep(0, nrow(dseq)),tc_num2=rep(0, nrow(dseq)),norm.tc_num1=rep(0, nrow(dseq)),norm.tc_num2=rep(0, nrow(dseq)))

newdat<-rbind(newdat, dseq)
rm(dseq, tl)

newdat$tramscript<-gsub("_.*","",newdat$tc_id)
newdat$pos_from_start<-as.numeric(gsub("_.*","",gsub(".*;","",sub("_",";",newdat$tc_id))))
newdat$pos_from_stop<-as.numeric(gsub("_.*","",gsub(".*;","",sub("_",";",sub("_",";",newdat$tc_id)))))
newdat$region<-gsub(".*_","",newdat$tc_id)

total<-aggregate(cbind(norm.tc_num1, norm.tc_num2)~tramscript, data=newdat, sum)

thr<-total[which(total$norm.tc_num1>=5 & total$norm.tc_num2>=5),]

tz<-which(newdat$tramscript %in% thr$tramscript)

newdat<-newdat[tz,]


maxs<-aggregate(cbind(norm.tc_num1, norm.tc_num2)~tramscript, data=newdat, max)

newdat$count<-1

posTrans<-aggregate(count~tramscript, data=newdat, sum)

norm<-maxs[rep(row.names(maxs), posTrans$count),]

newdat<-newdat[order(newdat$tramscript),]

newdat$frac.tc_num1<-newdat$norm.tc_num1/norm$norm.tc_num1
newdat$frac.tc_num2<-newdat$norm.tc_num2/norm$norm.tc_num2


avg<-aggregate(cbind(frac.tc_num1, frac.tc_num2)~pos_from_start+region, data=newdat, mean)
mel<-melt(avg, id.vars = c("pos_from_start", "region"), measure.vars = c("frac.tc_num1", "frac.tc_num2"), variable.name = "sample", value.name = "scaled_tc_rpm")
mel$scaled_tc_rpm<-ifelse(mel$region=="cds" & mel$pos_from_start>1500, NA, mel$scaled_tc_rpm)
mel$scaled_tc_rpm<-ifelse(mel$region=="utr5" & mel$pos_from_start>150, NA, mel$scaled_tc_rpm)
mel$scaled_tc_rpm<-ifelse(mel$region=="utr3" & mel$pos_from_start>800, NA, mel$scaled_tc_rpm)
mel<-subset(mel, !is.na(scaled_tc_rpm))
ggplot(mel, aes(pos_from_start,scaled_tc_rpm))+geom_line()+facet_wrap(~sample+region, scales = "free_x")


library(zoo)

ggplot(mel, aes(pos_from_start,scaled_tc_rpm))+geom_line(aes(y=rollmean(scaled_tc_rpm, 15,na.pad=T)))+facet_wrap(~sample+region, scales = "free_x")


avg<-aggregate(cbind(frac.tc_num1, frac.tc_num2)~pos_from_stop+region, data=newdat, mean)
mel<-melt(avg, id.vars = c("pos_from_stop", "region"), measure.vars = c("frac.tc_num1", "frac.tc_num2"), variable.name = "sample", value.name = "scaled_tc_rpm")
mel$scaled_tc_rpm<-ifelse(mel$region=="cds" & mel$pos_from_stop<(-1500), NA, mel$scaled_tc_rpm)
mel$scaled_tc_rpm<-ifelse(mel$region=="utr5" & mel$pos_from_stop<(-150), NA, mel$scaled_tc_rpm)
mel$scaled_tc_rpm<-ifelse(mel$region=="utr3" & mel$pos_from_stop<(-800), NA, mel$scaled_tc_rpm)
mel<-subset(mel, !is.na(scaled_tc_rpm))
ggplot(mel, aes(pos_from_stop,scaled_tc_rpm))+geom_line()+facet_wrap(~sample+region, scales = "free_x")
library(zoo)
ggplot(mel, aes(pos_from_stop,scaled_tc_rpm))+geom_line(aes(y=rollmean(scaled_tc_rpm, 15,na.pad=T)))+facet_wrap(~sample+region, scales = "free_x")



memHi<-subset(tsig_annot, localization_cat=="membrane" & loc_tar_CDS=="membrane_tc>5.56 & tc<65.26") 
memMid<-subset(tsig_annot, localization_cat=="membrane" & loc_tar_CDS=="membrane_tc>1.66 & tc<5.56") 
memLo<-subset(tsig_annot, localization_cat=="membrane" & loc_tar_CDS=="membrane_tc<1.66") 


tesmemS<-which(newdat$tramscript %in% memHi$transcript)
tesmemT<-which(newdat$tramscript %in% memMid$transcript)
tesmemN<-which(newdat$tramscript %in% memLo$transcript)

tesnmemS<-newdat[tesmemS,]
tesnmemT<-newdat[tesmemT,]
tesnmemN<-newdat[tesmemN,]

length(unique(tesnmemS$tramscript))
length(unique(tesnmemT$tramscript))
length(unique(tesnmemN$tramscript))

tesnmemS$localization<-"memHi"
tesnmemT$localization<-"memMid"
tesnmemN$localization<-"memLo"

locn<-rbind(tesnmemS,tesnmemT,tesnmemN)

avg<-aggregate(cbind(frac.tc_num1, frac.tc_num2)~pos_from_start+localization+region, data=locn, mean, na.rm=T)
mel<-melt(avg, id.vars = c("pos_from_start", "region", "localization"), measure.vars = c("frac.tc_num1", "frac.tc_num2"), variable.name = "sample", value.name = "scaled_tc_rpm")
mel$scaled_tc_rpm<-ifelse(mel$region=="cds" & mel$pos_from_start>1500, NA, mel$scaled_tc_rpm)
mel$scaled_tc_rpm<-ifelse(mel$region=="utr5" & mel$pos_from_start>150, NA, mel$scaled_tc_rpm)
mel$scaled_tc_rpm<-ifelse(mel$region=="utr3" & mel$pos_from_start>800, NA, mel$scaled_tc_rpm)
mel<-subset(mel, !is.na(scaled_tc_rpm))
ggplot(mel, aes(pos_from_start,scaled_tc_rpm, colour=localization))+geom_line()+facet_wrap(~sample+region, scales = "free_x")
ggplot(mel, aes(pos_from_start,scaled_tc_rpm, colour=localization))+geom_line(aes(y=rollmean(scaled_tc_rpm, 33,na.pad=T)))+facet_wrap(~sample+region, scales = "free_x")

avg<-aggregate(cbind(frac.tc_num1, frac.tc_num2)~pos_from_stop+localization+region, data=locn, mean, na.rm=T)
mel<-melt(avg, id.vars = c("pos_from_stop", "region", "localization"), measure.vars = c("frac.tc_num1", "frac.tc_num2"), variable.name = "sample", value.name = "scaled_tc_rpm")
mel$scaled_tc_rpm<-ifelse(mel$region=="cds" & mel$pos_from_stop<(-1500), NA, mel$scaled_tc_rpm)
mel$scaled_tc_rpm<-ifelse(mel$region=="utr5" & mel$pos_from_stop<(-150), NA, mel$scaled_tc_rpm)
mel$scaled_tc_rpm<-ifelse(mel$region=="utr3" & mel$pos_from_stop<(-800), NA, mel$scaled_tc_rpm)
mel<-subset(mel, !is.na(scaled_tc_rpm))
ggplot(mel, aes(pos_from_stop,scaled_tc_rpm, colour=localization))+geom_line()+facet_wrap(~sample+region, scales = "free_x")
ggplot(mel, aes(pos_from_stop,scaled_tc_rpm, colour=localization))+geom_line(aes(y=rollmean(scaled_tc_rpm, 33,na.pad=T)))+facet_wrap(~sample+region, scales = "free_x")




mem<-subset(tsig_annot, localization_cat=="membrane" ) 
cyt<-subset(tsig_annot, localization_cat=="cytosolic" ) 

tesmemS<-which(newdat$tramscript %in% mem$transcript)
tesmemT<-which(newdat$tramscript %in% cyt$transcript)

tesnmemS<-newdat[tesmemS,]
tesnmemT<-newdat[tesmemT,]

length(unique(tesnmemS$tramscript))
length(unique(tesnmemT$tramscript))

tesnmemS$localization<-"mem"
tesnmemT$localization<-"cyt"

locn<-rbind(tesnmemS,tesnmemT)

avg<-aggregate(cbind(frac.tc_num1, frac.tc_num2)~pos_from_start+localization+region, data=locn, mean, na.rm=T)
mel<-melt(avg, id.vars = c("pos_from_start", "region", "localization"), measure.vars = c("frac.tc_num1", "frac.tc_num2"), variable.name = "sample", value.name = "scaled_tc_rpm")
mel$scaled_tc_rpm<-ifelse(mel$region=="cds" & mel$pos_from_start>1500, NA, mel$scaled_tc_rpm)
mel$scaled_tc_rpm<-ifelse(mel$region=="utr5" & mel$pos_from_start>150, NA, mel$scaled_tc_rpm)
mel$scaled_tc_rpm<-ifelse(mel$region=="utr3" & mel$pos_from_start>800, NA, mel$scaled_tc_rpm)
mel<-subset(mel, !is.na(scaled_tc_rpm))
ggplot(mel, aes(pos_from_start,scaled_tc_rpm, colour=localization))+geom_line()+facet_wrap(~sample+region, scales = "free_x")
ggplot(mel, aes(pos_from_start,scaled_tc_rpm, colour=localization))+geom_line(aes(y=rollmean(scaled_tc_rpm, 33,na.pad=T)))+facet_wrap(~sample+region, scales = "free_x")

avg<-aggregate(cbind(frac.tc_num1, frac.tc_num2)~pos_from_stop+localization+region, data=locn, mean, na.rm=T)
mel<-melt(avg, id.vars = c("pos_from_stop", "region", "localization"), measure.vars = c("frac.tc_num1", "frac.tc_num2"), variable.name = "sample", value.name = "scaled_tc_rpm")
mel$scaled_tc_rpm<-ifelse(mel$region=="cds" & mel$pos_from_stop<(-1500), NA, mel$scaled_tc_rpm)
mel$scaled_tc_rpm<-ifelse(mel$region=="utr5" & mel$pos_from_stop<(-150), NA, mel$scaled_tc_rpm)
mel$scaled_tc_rpm<-ifelse(mel$region=="utr3" & mel$pos_from_stop<(-800), NA, mel$scaled_tc_rpm)
mel<-subset(mel, !is.na(scaled_tc_rpm))
ggplot(mel, aes(pos_from_stop,scaled_tc_rpm, colour=localization))+geom_line()+facet_wrap(~sample+region, scales = "free_x")
ggplot(mel, aes(pos_from_stop,scaled_tc_rpm, colour=localization))+geom_line(aes(y=rollmean(scaled_tc_rpm, 33,na.pad=T)))+facet_wrap(~sample+region, scales = "free_x")




##more cytosolic classes
cytHi<-subset(tsig_annot, localization_cat=="cytosolic" & tsig=="cyt_notsig" & loc_tar_CDS=="cytosolic_tc>0.47 & tc<27.84") 
cytMid<-subset(tsig_annot, localization_cat=="cytosolic" & tsig=="cyt_notsig" & loc_tar_CDS=="cytosolic_tc>0.18 & tc<0.47")
cytLo<-subset(tsig_annot, localization_cat=="cytosolic" & tsig=="cyt_notsig" & loc_tar_CDS=="cytosolic_tc<0.18")


tesmemS<-which(newdat$tramscript %in% cytHi$transcript)
tesmemT<-which(newdat$tramscript %in% cytMid$transcript)
tesmemN<-which(newdat$tramscript %in% cytLo$transcript)

tesnmemS<-newdat[tesmemS,]
tesnmemT<-newdat[tesmemT,]
tesnmemN<-newdat[tesmemN,]

length(unique(tesnmemS$tramscript))
length(unique(tesnmemT$tramscript))
length(unique(tesnmemN$tramscript))

tesnmemS$localization<-"cytHi"
tesnmemT$localization<-"cytMid"
tesnmemN$localization<-"cytLo"

locn<-rbind(tesnmemS,tesnmemT,tesnmemN)

avg<-aggregate(cbind(frac.tc_num1, frac.tc_num2)~pos_from_start+localization+region, data=locn, mean, na.rm=T)
mel<-melt(avg, id.vars = c("pos_from_start", "region", "localization"), measure.vars = c("frac.tc_num1", "frac.tc_num2"), variable.name = "sample", value.name = "scaled_tc_rpm")
mel$scaled_tc_rpm<-ifelse(mel$region=="cds" & mel$pos_from_start>1500, NA, mel$scaled_tc_rpm)
mel$scaled_tc_rpm<-ifelse(mel$region=="utr5" & mel$pos_from_start>150, NA, mel$scaled_tc_rpm)
mel$scaled_tc_rpm<-ifelse(mel$region=="utr3" & mel$pos_from_start>800, NA, mel$scaled_tc_rpm)
mel<-subset(mel, !is.na(scaled_tc_rpm))
# ggplot(mel, aes(pos_from_start,scaled_tc_rpm, colour=localization))+geom_line()+facet_wrap(~sample+region, scales = "free_x")
ggplot(mel, aes(pos_from_start,scaled_tc_rpm, colour=localization))+geom_line(aes(y=rollmean(scaled_tc_rpm, 33,na.pad=T)))+facet_wrap(~sample+region, scales = "free_x")

avg<-aggregate(cbind(frac.tc_num1, frac.tc_num2)~pos_from_stop+localization+region, data=locn, mean, na.rm=T)
mel<-melt(avg, id.vars = c("pos_from_stop", "region", "localization"), measure.vars = c("frac.tc_num1", "frac.tc_num2"), variable.name = "sample", value.name = "scaled_tc_rpm")
mel$scaled_tc_rpm<-ifelse(mel$region=="cds" & mel$pos_from_stop<(-1500), NA, mel$scaled_tc_rpm)
mel$scaled_tc_rpm<-ifelse(mel$region=="utr5" & mel$pos_from_stop<(-150), NA, mel$scaled_tc_rpm)
mel$scaled_tc_rpm<-ifelse(mel$region=="utr3" & mel$pos_from_stop<(-800), NA, mel$scaled_tc_rpm)
mel<-subset(mel, !is.na(scaled_tc_rpm))
# ggplot(mel, aes(pos_from_stop,scaled_tc_rpm, colour=localization))+geom_line()+facet_wrap(~sample+region, scales = "free_x")
ggplot(mel, aes(pos_from_stop,scaled_tc_rpm, colour=localization))+geom_line(aes(y=rollmean(scaled_tc_rpm, 33,na.pad=T)))+facet_wrap(~sample+region, scales = "free_x")

##targeting signals
memS<-subset(tsig_annot, localization_cat=="membrane" & grepl("SignalP", tsig_annot$tsig) )
memT<-subset(tsig_annot, localization_cat=="membrane" & tsig=="TMhelix-only") 
memN<-subset(tsig_annot, localization_cat=="membrane" & tsig=="mem_notsig") 
cyt<-subset(tsig_annot, localization_cat=="cytosolic" & tsig=="cyt_notsig") 

tesmemS<-which(newdat$tramscript %in% memS$transcript)
tesmemT<-which(newdat$tramscript %in% memT$transcript)
tesmemN<-which(newdat$tramscript %in% memN$transcript)
tescyt<-which(newdat$tramscript %in% cyt$transcript)

tesnmemS<-newdat[tesmemS,]
tesnmemT<-newdat[tesmemT,]
tesnmemN<-newdat[tesmemN,]
tesncyt<-newdat[tescyt,]

length(unique(tesnmemS$tramscript))
length(unique(tesnmemT$tramscript))
length(unique(tesnmemN$tramscript))
length(unique(tesncyt$tramscript))

tesnmemS$localization<-"memS"
tesnmemT$localization<-"memT"
tesnmemN$localization<-"memN"
tesncyt$localization<-"cyt"

locn<-rbind(tesnmemS,tesnmemT,tesnmemN,tesncyt)

avg<-aggregate(cbind(frac.tc_num1, frac.tc_num2)~pos_from_start+localization+region, data=locn, mean, na.rm=T)
mel<-melt(avg, id.vars = c("pos_from_start", "region", "localization"), measure.vars = c("frac.tc_num1", "frac.tc_num2"), variable.name = "sample", value.name = "scaled_tc_rpm")
mel$scaled_tc_rpm<-ifelse(mel$region=="cds" & mel$pos_from_start>1500, NA, mel$scaled_tc_rpm)
mel$scaled_tc_rpm<-ifelse(mel$region=="utr5" & mel$pos_from_start>150, NA, mel$scaled_tc_rpm)
mel$scaled_tc_rpm<-ifelse(mel$region=="utr3" & mel$pos_from_start>800, NA, mel$scaled_tc_rpm)
mel<-subset(mel, !is.na(scaled_tc_rpm))
ggplot(mel[mel$localization!="memN",], aes(pos_from_start,scaled_tc_rpm, colour=localization))+geom_line()+facet_wrap(~sample+region, scales = "free_x")
ggplot(mel[mel$localization!="memN",], aes(pos_from_start,scaled_tc_rpm, colour=localization))+geom_line(aes(y=rollmean(scaled_tc_rpm, 33,na.pad=T)))+facet_wrap(~sample+region, scales = "free_x")

ggplot(mel[mel$localization!="memN",], aes(pos_from_start,scaled_tc_rpm, colour=localization))+geom_line(aes(y=rollmean(scaled_tc_rpm, 33,na.pad=T)))+facet_wrap(~sample+region, scales = "free_x")+coord_cartesian(xlim=c(0,500))

avg<-aggregate(cbind(frac.tc_num1, frac.tc_num2)~pos_from_stop+localization+region, data=locn, mean, na.rm=T)
mel<-melt(avg, id.vars = c("pos_from_stop", "region", "localization"), measure.vars = c("frac.tc_num1", "frac.tc_num2"), variable.name = "sample", value.name = "scaled_tc_rpm")
mel$scaled_tc_rpm<-ifelse(mel$region=="cds" & mel$pos_from_stop<(-1500), NA, mel$scaled_tc_rpm)
mel$scaled_tc_rpm<-ifelse(mel$region=="utr5" & mel$pos_from_stop<(-150), NA, mel$scaled_tc_rpm)
mel$scaled_tc_rpm<-ifelse(mel$region=="utr3" & mel$pos_from_stop<(-800), NA, mel$scaled_tc_rpm)
mel<-subset(mel, !is.na(scaled_tc_rpm))
ggplot(mel[mel$localization!="memN",], aes(pos_from_stop,scaled_tc_rpm, colour=localization))+geom_line()+facet_wrap(~sample+region, scales = "free_x")
ggplot(mel[mel$localization!="memN",], aes(pos_from_stop,scaled_tc_rpm, colour=localization))+geom_line(aes(y=rollmean(scaled_tc_rpm, 33,na.pad=T)))+facet_wrap(~sample+region, scales = "free_x")


###reanalyze CDS and 3UTR
setwd("~/Google Drive/hdlbp/")
feat<-read.delim("hdlbp_master_table_with_classes_uniq.txt", header=T)

feat$tc_start_norm<-ifelse(is.na(feat$tc_start), NA,
                         ifelse(feat$tpm_cutoff==0 & feat$tc_start>0, feat$tc_start/min(feat[feat$tpm_cutoff!=0,"tpm_cutoff"], na.rm=T),
                                feat$tc_start/feat$tpm_cutoff))

feat$tc_stop_norm<-ifelse(is.na(feat$tc_stop), NA,
                         ifelse(feat$tpm_cutoff==0 & feat$tc_stop>0, feat$tc_stop/min(feat[feat$tpm_cutoff!=0,"tpm_cutoff"], na.rm=T),
                                feat$tc_stop/feat$tpm_cutoff))
thr<-10
quants_mem<-round(quantile(feat[feat$tpm_cutoff>=thr &feat$localization_cat=="membrane","tc_start_norm"], probs=c(0,1/3,2/3,1), na.rm=T),digits=2)
quants_cyt<-round(quantile(feat[feat$tpm_cutoff>=thr &feat$localization_cat=="cytosolic","tc_start_norm"], probs=c(0,1/3,2/3,1), na.rm=T),digits=2)
feat$loc_tar_start<-ifelse(is.na(feat$tc_start_norm) & feat$tpm_cutoff<thr , NA,
                         ifelse(is.na(feat$tc_start_norm) & feat$tpm_cutoff>=thr &feat$localization_cat=="membrane", "nontarget_membrane",
                                ifelse(feat$tpm_cutoff>=thr & feat$tc_start_norm<=quants_mem[2] & feat$tc_start_norm>quants_mem[1] &feat$localization_cat=="membrane", paste0("membrane_tc<",quants_mem[2]),
                                       ifelse(feat$tpm_cutoff>=thr & feat$tc_start_norm<=quants_mem[3] & feat$tc_start_norm>quants_mem[2] &feat$localization_cat=="membrane", paste0("membrane_tc>",quants_mem[2]," & tc<",quants_mem[3]),
                                              ifelse(feat$tpm_cutoff>=thr & feat$tc_start_norm<=quants_mem[4] & feat$tc_start_norm>quants_mem[3] &feat$localization_cat=="membrane", paste0("membrane_tc>",quants_mem[3]," & tc<",quants_mem[4]),
                                                     ifelse(is.na(feat$tc_start_norm) & feat$tpm_cutoff>=thr &feat$localization_cat=="cytosolic", "nontarget_cytosolic",
                                                            ifelse(feat$tpm_cutoff>=thr & feat$tc_start_norm<=quants_cyt[2] & feat$tc_start_norm>quants_cyt[1] &feat$localization_cat=="cytosolic", paste0("cytosolic_tc<",quants_cyt[2]),
                                                                   ifelse(feat$tpm_cutoff>=thr & feat$tc_start_norm<=quants_cyt[3] & feat$tc_start_norm>quants_cyt[2] &feat$localization_cat=="cytosolic", paste0("cytosolic_tc>",quants_cyt[2]," & tc<",quants_cyt[3]),
                                                                          ifelse(feat$tpm_cutoff>=thr & feat$tc_start_norm<=quants_cyt[4] & feat$tc_start_norm>quants_cyt[3] &feat$localization_cat=="cytosolic", paste0("cytosolic_tc>",quants_cyt[3]," & tc<",quants_cyt[4]),NA)))))))))

quants_mem<-round(quantile(feat[feat$tpm_cutoff>=thr &feat$localization_cat=="membrane","tc_stop_norm"], probs=c(0,1/3,2/3,1), na.rm=T),digits=2)
quants_cyt<-round(quantile(feat[feat$tpm_cutoff>=thr &feat$localization_cat=="cytosolic","tc_stop_norm"], probs=c(0,1/3,2/3,1), na.rm=T),digits=2)
feat$loc_tar_stop<-ifelse(is.na(feat$tc_stop_norm) & feat$tpm_cutoff<thr , NA,
                          ifelse(is.na(feat$tc_stop_norm) & feat$tpm_cutoff>=thr &feat$localization_cat=="membrane", "nontarget_membrane",
                                 ifelse(feat$tpm_cutoff>=thr & feat$tc_stop_norm<=quants_mem[2] & feat$tc_stop_norm>quants_mem[1] &feat$localization_cat=="membrane", paste0("membrane_tc<",quants_mem[2]),
                                        ifelse(feat$tpm_cutoff>=thr & feat$tc_stop_norm<=quants_mem[3] & feat$tc_stop_norm>quants_mem[2] &feat$localization_cat=="membrane", paste0("membrane_tc>",quants_mem[2]," & tc<",quants_mem[3]),
                                               ifelse(feat$tpm_cutoff>=thr & feat$tc_stop_norm<=quants_mem[4] & feat$tc_stop_norm>quants_mem[3] &feat$localization_cat=="membrane", paste0("membrane_tc>",quants_mem[3]," & tc<",quants_mem[4]),
                                                      ifelse(is.na(feat$tc_stop_norm) & feat$tpm_cutoff>=thr &feat$localization_cat=="cytosolic", "nontarget_cytosolic",
                                                             ifelse(feat$tpm_cutoff>=thr & feat$tc_stop_norm<=quants_cyt[2] & feat$tc_stop_norm>quants_cyt[1] &feat$localization_cat=="cytosolic", paste0("cytosolic_tc<",quants_cyt[2]),
                                                                    ifelse(feat$tpm_cutoff>=thr & feat$tc_stop_norm<=quants_cyt[3] & feat$tc_stop_norm>quants_cyt[2] &feat$localization_cat=="cytosolic", paste0("cytosolic_tc>",quants_cyt[2]," & tc<",quants_cyt[3]),
                                                                           ifelse(feat$tpm_cutoff>=thr & feat$tc_stop_norm<=quants_cyt[4] & feat$tc_stop_norm>quants_cyt[3] &feat$localization_cat=="cytosolic", paste0("cytosolic_tc>",quants_cyt[3]," & tc<",quants_cyt[4]),NA)))))))))

feat<-merge(feat, regionCounts[,c(1:11,25:29,37,39:ncol(regionCounts))], by="gene_id", all=T)

rel<-subset(feat, tpm_cutoff>=10 & Annotation=="protein_coding" )

ggplot(rel[!is.na(rel$localization_cat),] ,aes(log2(utr3_vs_cds1), log2(utr3_vs_cds2), colour=localization_cat))+geom_point(alpha=0.7, shape=1)
ggplot(rel[!is.na(rel$localization_cat),] ,aes(log2(utr3_vs_cds1), colour=localization_cat))+geom_density()
ggplot(rel[!is.na(rel$localization_cat),] ,aes(log2(utr3_vs_cds1), fill=localization_cat))+geom_histogram(binwidth = 0.3)

ggplot(rel[!is.na(rel$localization_cat),] ,aes(log2(utr3_vs_cds1.norm), log2(utr3_vs_cds2.norm), colour=localization_cat))+geom_point(alpha=0.7, shape=1)
ggplot(rel[!is.na(rel$localization_cat),] ,aes(log2(utr3_vs_cds1.norm), colour=localization_cat))+geom_density()
ggplot(rel[!is.na(rel$localization_cat),] ,aes(log2(utr3_vs_cds1.norm), fill=localization_cat))+geom_histogram(binwidth = 0.3)

ggplot(rel[!is.na(rel$localization_cat),] ,aes(log2(utr3_vs_cds1.norm), log2(tc_transcript_norm), colour=localization_cat))+geom_point()
ggplot(rel[!is.na(rel$localization_cat),] ,aes(log2(utr3_vs_cds1.norm), log2(trans.norm.tc_num1/tpm_cutoff), colour=localization_cat))+geom_point()
plot(log2(rel[!is.na(rel$localization_cat),"utr3_vs_cds1"]),log2(rel[!is.na(rel$localization_cat),"trans.norm.tc_num1"]/rel[!is.na(rel$localization_cat),"tpm_cutoff"]))
identify(log2(rel[!is.na(rel$localization_cat),"utr3_vs_cds1"]),log2(rel[!is.na(rel$localization_cat),"trans.norm.tc_num1"]/rel[!is.na(rel$localization_cat),"tpm_cutoff"]), labels = rel$Symbol)

ggplot(rel[!is.na(rel$localization_cat) & rel$trans.norm.tc_num1>=50,] ,aes(log2(utr3_vs_cds1.norm), log2(trans.norm.tc_num1.len/tpm_cutoff), colour=localization_cat))+geom_point()
subset(rel, trans.norm.tc_num1>=50 & log2(trans.norm.tc_num1/tpm_cutoff)>3 & log2(utr3_vs_cds1)>3, select=c("Symbol", "gene_id"))
subset(rel, trans.norm.tc_num1>=50 & log2(trans.norm.tc_num1/tpm_cutoff)>4 & log2(utr3_vs_cds1)<(-5) & tpm_cutoff>50, select=c("Symbol", "gene_id"))

ggplot(rel[!is.na(rel$localization_cat),] ,aes(log2(utr3_vs_cds2), log2(trans.norm.tc_num2/tpm_cutoff), colour=localization_cat))+geom_point()
ggplot(rel[!is.na(rel$localization_cat),] ,aes(log2(utr3_vs_cds1), log2(tpm_cutoff), colour=localization_cat))+geom_point()

ggplot(rel[!is.na(rel$loc_tar_transcript),], aes(log2(utr3_vs_cds1), log2(trans.norm.tc_num1/tpm_cutoff), colour=loc_tar_transcript))+geom_point(alpha=0.8)
ggplot(rel[!is.na(rel$tc_transcript_cat),], aes(log2(utr3_vs_cds1), log2(trans.norm.tc_num1/tpm_cutoff), colour=tc_transcript_cat))+geom_point(alpha=0.8)

ggplot(rel,aes(log2(utr3_vs_cds1), log2FoldChange.mem.cyt.293, colour=loc_tar_transcript))+geom_point()+ylim(-2.5,6)

ggplot(rel,aes(log2(utr3_vs_cds1), log2FoldChange.mem.cyt.KO.293, colour=loc_tar_transcript))+geom_point()+ylim(-1,1)
ggplot(rel,aes(log2(utr3_vs_cds1), log2FoldChange.ribo.rna.KO.WT, colour=loc_tar_transcript))+geom_point()+ylim(-1,1)

feat$av_utr3_vs_cds<-rowMeans(feat[,c("utr3_vs_cds1", "utr3_vs_cds2")])
feat$av_utr3_vs_cds.norm<-rowMeans(feat[,c("utr3_vs_cds1.norm", "utr3_vs_cds2.norm")])

thr<-10
quants_mem<-round(quantile(feat[feat$tpm_cutoff>=thr &feat$localization_cat=="membrane","av_utr3_vs_cds"], probs=c(0,1/3,2/3,1), na.rm=T),digits=2)
quants_cyt<-round(quantile(feat[feat$tpm_cutoff>=thr &feat$localization_cat=="cytosolic","av_utr3_vs_cds"], probs=c(0,1/3,2/3,1), na.rm=T),digits=2)
feat$loc_utr3_vs_cds<-ifelse(is.na(feat$av_utr3_vs_cds) & feat$tpm_cutoff<thr , NA,
                           ifelse(is.na(feat$av_utr3_vs_cds) & feat$tpm_cutoff>=thr &feat$localization_cat=="membrane", NA,
                                  ifelse(feat$tpm_cutoff>=thr & feat$av_utr3_vs_cds<=quants_mem[2] & feat$av_utr3_vs_cds>quants_mem[1] &feat$localization_cat=="membrane", paste0("membrane_utr3_cds<",quants_mem[2]),
                                         ifelse(feat$tpm_cutoff>=thr & feat$av_utr3_vs_cds<=quants_mem[3] & feat$av_utr3_vs_cds>quants_mem[2] &feat$localization_cat=="membrane", paste0("membrane_utr3_cds>",quants_mem[2]," & tc<",quants_mem[3]),
                                                ifelse(feat$tpm_cutoff>=thr & feat$av_utr3_vs_cds<=quants_mem[4] & feat$av_utr3_vs_cds>quants_mem[3] &feat$localization_cat=="membrane", paste0("membrane_utr3_cds>",quants_mem[3]," & tc<",quants_mem[4]),
                                                       ifelse(is.na(feat$av_utr3_vs_cds) & feat$tpm_cutoff>=thr &feat$localization_cat=="cytosolic", NA,
                                                              ifelse(feat$tpm_cutoff>=thr & feat$av_utr3_vs_cds<=quants_cyt[2] & feat$av_utr3_vs_cds>quants_cyt[1] &feat$localization_cat=="cytosolic", paste0("cytosolic_utr3_cds<",quants_cyt[2]),
                                                                     ifelse(feat$tpm_cutoff>=thr & feat$av_utr3_vs_cds<=quants_cyt[3] & feat$av_utr3_vs_cds>quants_cyt[2] &feat$localization_cat=="cytosolic", paste0("cytosolic_utr3_cds>",quants_cyt[2]," & tc<",quants_cyt[3]),
                                                                            ifelse(feat$tpm_cutoff>=thr & feat$av_utr3_vs_cds<=quants_cyt[4] & feat$av_utr3_vs_cds>quants_cyt[3] &feat$localization_cat=="cytosolic", paste0("cytosolic_utr3_cds>",quants_cyt[3]," & tc<",quants_cyt[4]),NA)))))))))

quants<-round(quantile(feat[feat$tpm_cutoff>=thr,"av_utr3_vs_cds"], probs=c(0,1/3,2/3,1), na.rm=T),digits=2)
feat$av_utr3_vs_cds_cat<-ifelse(is.na(feat$av_utr3_vs_cds) & feat$tpm_cutoff<thr, NA,
                            ifelse(is.na(feat$av_utr3_vs_cds) & feat$tpm_cutoff>=thr, "nontarget",
                                   ifelse(feat$tpm_cutoff>=thr & feat$av_utr3_vs_cds<=quants[2] & feat$av_utr3_vs_cds>quants[1], paste0("tc<",quants[2]),
                                          ifelse(feat$tpm_cutoff>=thr & feat$av_utr3_vs_cds<=quants[3] & feat$av_utr3_vs_cds>quants[2], paste0("tc>",quants[2]," & tc<",quants[3]),
                                                 ifelse(feat$tpm_cutoff>=thr & feat$av_utr3_vs_cds<=quants[4] & feat$av_utr3_vs_cds>quants[3], paste0("tc>",quants[3]," & tc<",quants[4]), NA)))))


quants_mem<-round(quantile(feat[feat$tpm_cutoff>=thr &feat$localization_cat=="membrane","av_utr3_vs_cds.norm"], probs=c(0,1/3,2/3,1), na.rm=T),digits=2)
quants_cyt<-round(quantile(feat[feat$tpm_cutoff>=thr &feat$localization_cat=="cytosolic","av_utr3_vs_cds.norm"], probs=c(0,1/3,2/3,1), na.rm=T),digits=2)
feat$loc_utr3_vs_cds.norm<-ifelse(is.na(feat$av_utr3_vs_cds.norm) & feat$tpm_cutoff<thr , NA,
                             ifelse(is.na(feat$av_utr3_vs_cds.norm) & feat$tpm_cutoff>=thr &feat$localization_cat=="membrane", NA,
                                    ifelse(feat$tpm_cutoff>=thr & feat$av_utr3_vs_cds.norm<=quants_mem[2] & feat$av_utr3_vs_cds.norm>quants_mem[1] &feat$localization_cat=="membrane", paste0("membrane_utr3_cds<",quants_mem[2]),
                                           ifelse(feat$tpm_cutoff>=thr & feat$av_utr3_vs_cds.norm<=quants_mem[3] & feat$av_utr3_vs_cds.norm>quants_mem[2] &feat$localization_cat=="membrane", paste0("membrane_utr3_cds>",quants_mem[2]," & tc<",quants_mem[3]),
                                                  ifelse(feat$tpm_cutoff>=thr & feat$av_utr3_vs_cds.norm<=quants_mem[4] & feat$av_utr3_vs_cds.norm>quants_mem[3] &feat$localization_cat=="membrane", paste0("membrane_utr3_cds>",quants_mem[3]," & tc<",quants_mem[4]),
                                                         ifelse(is.na(feat$av_utr3_vs_cds.norm) & feat$tpm_cutoff>=thr &feat$localization_cat=="cytosolic", NA,
                                                                ifelse(feat$tpm_cutoff>=thr & feat$av_utr3_vs_cds.norm<=quants_cyt[2] & feat$av_utr3_vs_cds.norm>quants_cyt[1] &feat$localization_cat=="cytosolic", paste0("cytosolic_utr3_cds<",quants_cyt[2]),
                                                                       ifelse(feat$tpm_cutoff>=thr & feat$av_utr3_vs_cds.norm<=quants_cyt[3] & feat$av_utr3_vs_cds.norm>quants_cyt[2] &feat$localization_cat=="cytosolic", paste0("cytosolic_utr3_cds>",quants_cyt[2]," & tc<",quants_cyt[3]),
                                                                              ifelse(feat$tpm_cutoff>=thr & feat$av_utr3_vs_cds.norm<=quants_cyt[4] & feat$av_utr3_vs_cds.norm>quants_cyt[3] &feat$localization_cat=="cytosolic", paste0("cytosolic_utr3_cds>",quants_cyt[3]," & tc<",quants_cyt[4]),NA)))))))))

quants<-round(quantile(feat[feat$tpm_cutoff>=thr,"av_utr3_vs_cds.norm"], probs=c(0,1/3,2/3,1), na.rm=T),digits=2)
feat$av_utr3_vs_cds.norm_cat<-ifelse(is.na(feat$av_utr3_vs_cds.norm) & feat$tpm_cutoff<thr, NA,
                                ifelse(is.na(feat$av_utr3_vs_cds.norm) & feat$tpm_cutoff>=thr, "nontarget",
                                       ifelse(feat$tpm_cutoff>=thr & feat$av_utr3_vs_cds.norm<=quants[2] & feat$av_utr3_vs_cds.norm>quants[1], paste0("tc<",quants[2]),
                                              ifelse(feat$tpm_cutoff>=thr & feat$av_utr3_vs_cds.norm<=quants[3] & feat$av_utr3_vs_cds.norm>quants[2], paste0("tc>",quants[2]," & tc<",quants[3]),
                                                     ifelse(feat$tpm_cutoff>=thr & feat$av_utr3_vs_cds.norm<=quants[4] & feat$av_utr3_vs_cds.norm>quants[3], paste0("tc>",quants[3]," & tc<",quants[4]), NA)))))


rel<-subset(feat, tpm_cutoff>=10 & Annotation=="protein_coding" )

ggplot(rel,aes( mean_te_293, colour=loc_utr3_vs_cds))+stat_ecdf()+xlim(0,3)
ggplot(rel,aes( log2FoldChange.mem.cyt.KO.293, colour=loc_utr3_vs_cds))+stat_ecdf()+coord_cartesian(xlim=c(-0.5,0.5))
ggplot(rel,aes( log2FoldChange.mem.cyt.KO.293, colour=loc_utr3_vs_cds.norm))+stat_ecdf()+coord_cartesian(xlim=c(-0.5,0.5))
ggplot(rel,aes( log2FoldChange.mem.cyt.KO.293, colour=av_utr3_vs_cds.norm_cat))+stat_ecdf()+coord_cartesian(xlim=c(-0.5,0.5))


ggplot(rel,aes( log2FoldChange.ribo.rna.KO.WT, colour=loc_utr3_vs_cds))+stat_ecdf()+coord_cartesian(xlim=c(-0.5,0.5))

ggplot(rel,aes( log2FoldChange.ribo.rna.KO.WT, colour=loc_utr3_vs_cds.norm))+stat_ecdf()+coord_cartesian(xlim=c(-0.5,0.5))
ggplot(rel,aes( log2FoldChange.ribo.rna.KO.WT, colour=av_utr3_vs_cds.norm_cat))+stat_ecdf()+coord_cartesian(xlim=c(-0.5,0.5))

ggplot(rel,aes( log2(tpm_cutoff), colour=loc_utr3_vs_cds))+stat_ecdf()
ggplot(rel,aes( log2(tpm_cutoff), colour=loc_utr3_vs_cds.norm))+stat_ecdf()
ggplot(rel,aes( log2FoldChange.tot.KO.WT, colour=loc_utr3_vs_cds))+stat_ecdf()+coord_cartesian(xlim=c(-0.5,0.5))
ggplot(rel,aes(log2(utr3_vs_cds1), log2FoldChange.mem.cyt.293, colour=tc_transcript_cat))+geom_point()+ylim(-2.5,6)
ggplot(rel,aes(log2(tpm_cutoff), log2FoldChange.mem.cyt.293, colour=tc_transcript_cat))+geom_point()+ylim(-2.5,6)
ggplot(rel[rel$localization_cat!="membrane"& rel$loc_tar_CDS!="nontarget_cytosolic",],aes(log2(utr3_vs_cds1), log2(tpm_cutoff), colour=loc_tar_CDS))+geom_point()
ggplot(rel[rel$localization_cat!="membrane",],aes(log2(tpm_cutoff), colour=tc_CDS_norm_cat))+stat_ecdf()
ggplot(rel[rel$localization_cat!="membrane",],aes(log2(tpm_cutoff), colour=loc_utr3_vs_cds))+stat_ecdf()
ggplot(rel[rel$localization_cat!="membrane",],aes(log2FoldChange.tot.KO.WT, colour=loc_utr3_vs_cds))+stat_ecdf()+coord_cartesian(xlim=c(-0.5,0.5))
ggplot(rel[rel$localization_cat!="membrane",],aes(log2FoldChange.cyt.KO.WT, colour=av_utr3_vs_cds_cat))+stat_ecdf()+coord_cartesian(xlim=c(-0.5,0.5))
ggplot(rel[rel$localization_cat!="membrane",],aes(log2FoldChange.tot.KO.WT, colour=av_utr3_vs_cds_cat))+stat_ecdf()+coord_cartesian(xlim=c(-0.5,0.5))
ggplot(rel[rel$localization_cat!="membrane",],aes(log2(tpm_cutoff), colour=av_utr3_vs_cds_cat))+stat_ecdf()

ggplot(rel[rel$localization_cat!="membrane",],aes(log2FoldChange.nuc.cyt.KO.293, colour=loc_utr3_vs_cds))+stat_ecdf()+coord_cartesian(xlim=c(-.5,.5))
ggplot(rel[,],aes(log2FoldChange.nuc.cyt.KO.293, colour=av_utr3_vs_cds_cat))+stat_ecdf()+coord_cartesian(xlim=c(-.5,.5))
ggplot(rel[!is.na(rel$loc_utr3_vs_cds),],aes(log2FoldChange.nuc.cyt.KO.293, log2FoldChange.mem.cyt.KO.293, colour=loc_utr3_vs_cds))+geom_point(alpha=0.6,shape=1)+coord_cartesian(xlim=c(-2,2),ylim=c(-2,2))
ggplot(rel[!is.na(rel$loc_utr3_vs_cds)& rel$localization_cat!="membrane",],aes(log2FoldChange.nuc.cyt.KO.293, log2FoldChange.mem.cyt.KO.293, colour=loc_utr3_vs_cds))+geom_point(alpha=0.6,shape=1)+coord_cartesian(xlim=c(-2,2),ylim=c(-2,2))
ggplot(rel[!is.na(rel$loc_utr3_vs_cds)& rel$localization_cat!="cytosolic",],aes(log2FoldChange.nuc.cyt.KO.293, log2FoldChange.mem.cyt.KO.293, colour=loc_utr3_vs_cds))+geom_point(alpha=0.6,shape=1)+coord_cartesian(xlim=c(-2,2),ylim=c(-2,2))
plot(rel[!is.na(rel$loc_utr3_vs_cds)& rel$localization_cat!="cytosolic","log2FoldChange.nuc.cyt.KO.293"],rel[!is.na(rel$loc_utr3_vs_cds)& rel$localization_cat!="cytosolic","log2FoldChange.mem.cyt.KO.293"])
identify(rel[!is.na(rel$loc_utr3_vs_cds)& rel$localization_cat!="cytosolic","log2FoldChange.nuc.cyt.KO.293"],rel[!is.na(rel$loc_utr3_vs_cds)& rel$localization_cat!="cytosolic","log2FoldChange.mem.cyt.KO.293"],labels=rel$Symbol)

ggplot(rel[,],aes(log2FoldChange.nuc.cyt.293, log2FoldChange.mem.cyt.293, colour=loc_utr3_vs_cds))+geom_point(alpha=0.6,shape=1)+coord_cartesian(xlim=c(-3,6),ylim=c(-3,6))


ggplot(rel[rel$localization_cat!="membrane",],aes(log2FoldChange.mem.nuc.KO.WT, colour=loc_utr3_vs_cds))+stat_ecdf()+coord_cartesian(xlim=c(-.5,.5))

ggplot(rel[rel$localization_cat!="membrane",],aes(log2FoldChange.nuc.tot.2A15.293, colour=loc_utr3_vs_cds))+stat_ecdf()+coord_cartesian(xlim=c(-.5,.5))


ggplot(rel,aes( log2FoldChange.mem.cyt.KO.293, colour=loc_tar_transcript))+stat_ecdf()+coord_cartesian(xlim=c(-0.5,0.5))
ggplot(rel,aes( log2FoldChange.ribo.rna.KO.WT, colour=loc_tar_transcript))+stat_ecdf()+coord_cartesian(xlim=c(-0.5,0.5))

ggplot(rel,aes( log2FoldChange.mem.cyt.KO.293, colour=loc_tar_CDS))+stat_ecdf()+coord_cartesian(xlim=c(-0.5,0.5))
ggplot(rel,aes( log2FoldChange.ribo.rna.KO.WT, colour=loc_tar_CDS))+stat_ecdf()+coord_cartesian(xlim=c(-0.5,0.5))

ggplot(rel,aes( log2FoldChange.mem.cyt.KO.293, colour=loc_tar_stop))+stat_ecdf()+coord_cartesian(xlim=c(-0.5,0.5))

ggplot(rel,aes( log2FoldChange.mem.cyt.293, colour=loc_utr3_vs_cds))+stat_ecdf()+coord_cartesian(xlim=c(-2,0.5))
ggplot(rel,aes( log2FoldChange.mem.cyt.293, colour=av_utr3_vs_cds_cat))+stat_ecdf()

ggplot(rel,aes( log2FoldChange.mem.cyt.KO.293, colour=av_utr3_vs_cds_cat))+stat_ecdf()+coord_cartesian(xlim=c(-0.5,0.5))

ggplot(rel[rel$localization_cat=="membrane",],aes( log2FoldChange.mem.cyt.KO.293, colour=av_utr3_vs_cds_cat))+stat_ecdf()+coord_cartesian(xlim=c(-0.5,1))
ggplot(rel[rel$localization_cat=="cytosolic",],aes( log2FoldChange.mem.cyt.KO.293, colour=av_utr3_vs_cds_cat))+stat_ecdf()+coord_cartesian(xlim=c(-0.5,1))
ggplot(rel[rel$localization_cat=="membrane",],aes( log2FoldChange.ribo.rna.KO.WT, colour=av_utr3_vs_cds_cat))+stat_ecdf()+coord_cartesian(xlim=c(-0.5,1))
ggplot(rel[rel$localization_cat=="cytosolic",],aes( log2FoldChange.ribo.rna.KO.WT, colour=av_utr3_vs_cds_cat))+stat_ecdf()+coord_cartesian(xlim=c(-0.5,1))

ggplot(rel,aes( log2FoldChange.ribo.rna.KO.WT, colour=av_utr3_vs_cds_cat))+stat_ecdf()+coord_cartesian(xlim=c(-0.5,0.5))



#### make list for reporter finding

list<-subset(feat, tpm_cutoff>=10 & Annotation=="protein_coding" & localization_cat=="membrane", 
             select=c("gene_id", "Symbol","tpm_cutoff","transcript","Transmembrane.helices", "Cleavage.site..Signalp.", "log2FoldChange.mem.cyt.KO.293", 
                      "padj.mem.cyt.KO.293" , "log2FoldChange.ribo.rna.KO.WT" ,"padj.ribo.rna.KO.WT","localization_cat",
                      "tc_transcript", "tc_CDS","tc_UTR3","tc_stop",
                      "tc_transcript_norm","tc_CDS_norm","tc_UTR3_norm" ,"tc_stop_norm", 
                      "loc_tar_transcript","loc_tar_CDS" ,"loc_tar_UTR3","loc_tar_stop",
                      "cds.norm.tc_num1", "cds.norm.tc_num2", "utr3.norm.tc_num1", "utr3.norm.tc_num2",
                      "cds.norm.tc_num1.len" ,"cds.norm.tc_num2.len" ,"utr3.norm.tc_num1.len", "utr3.norm.tc_num2.len",
                      "av_utr3_vs_cds","av_utr3_vs_cds.norm","loc_utr3_vs_cds" ))
setwd("~/Google Drive/hdlbp/reporter/")
write.table(list, "mem_hdlbp_reporter_list.txt", quote=F,sep="\t", row.names=F, col.names = T)

list<-subset(feat, tpm_cutoff>=10 & Annotation=="protein_coding" & localization_cat=="cytosolic", 
             select=c("gene_id", "Symbol","tpm_cutoff","transcript","Transmembrane.helices", "Cleavage.site..Signalp.", "log2FoldChange.mem.cyt.KO.293", 
                      "padj.mem.cyt.KO.293" , "log2FoldChange.ribo.rna.KO.WT" ,"padj.ribo.rna.KO.WT","localization_cat",
                      "tc_transcript", "tc_CDS","tc_UTR3","tc_stop",
                      "tc_transcript_norm","tc_CDS_norm","tc_UTR3_norm" ,"tc_stop_norm", 
                      "loc_tar_transcript","loc_tar_CDS" ,"loc_tar_UTR3","loc_tar_stop",
                      "cds.norm.tc_num1", "cds.norm.tc_num2", "utr3.norm.tc_num1", "utr3.norm.tc_num2",
                      "cds.norm.tc_num1.len" ,"cds.norm.tc_num2.len" ,"utr3.norm.tc_num1.len", "utr3.norm.tc_num2.len",
                      "av_utr3_vs_cds","av_utr3_vs_cds.norm","loc_utr3_vs_cds" ))
setwd("~/Google Drive/hdlbp/reporter/")
write.table(list, "cyto_hdlbp_reporter_list.txt", quote=F,sep="\t", row.names=F, col.names = T)

