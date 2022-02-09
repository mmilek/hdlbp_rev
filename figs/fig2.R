library(ggplot2)
library(reshape2)
library(corrplot)

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

ggplot(regionCounts,aes(log2(cds.norm.tc_num1.len),log2(cds.norm.tc_num2.len), colour=localization_cat))+geom_point()
ggplot(regionCounts,aes(log2(utr3.norm.tc_num1.len),log2(utr3.norm.tc_num2.len), colour=localization_cat))+geom_point()

ggplot(regionCounts, aes(log2(utr3.norm.tc_num2.len/cds.norm.tc_num2.len), log2(trans.norm.tc_num1.len), colour=localization_cat))+geom_point()


regionCounts$trans.norm.tc.mean<-rowMeans(regionCounts[,c("trans.norm.tc_num1", "trans.norm.tc_num2")])
regionCounts$trans.norm.tc.mean.exp<- ifelse(regionCounts$tpm_cutoff==0, NA, regionCounts$trans.norm.tc.mean/regionCounts$tpm_cutoff)

regionCounts$cds.norm.tc.mean<-rowMeans(regionCounts[,c("cds.norm.tc_num1", "cds.norm.tc_num2")])
regionCounts$utr3.norm.tc.mean<-rowMeans(regionCounts[,c("utr3.norm.tc_num1", "utr3.norm.tc_num2")])
regionCounts$utr3_vs_cds.mean<-rowMeans(regionCounts[,c("utr3_vs_cds1", "utr3_vs_cds2")])
regionCounts$utr3_vs_cds.norm.mean<-rowMeans(regionCounts[,c("utr3_vs_cds1.norm", "utr3_vs_cds2.norm")])


#fig2b
ggplot(subset(regionCounts, !is.na(localization_cat) & tpm_cutoff>=10 & gene_biotype=="protein_coding"), aes(-log2(utr3_vs_cds.norm.mean), log2(trans.norm.tc.mean.exp), colour=localization_cat))+geom_point(shape=1)+
  scale_colour_manual(values=c("orange2","dodgerblue3"))

#fig2c
ggplot(subset(regionCounts, !is.na(localization_cat) & tpm_cutoff>=10 & gene_biotype=="protein_coding" & tc_CDS_norm_cat!="nontarget"),
         aes(tc_CDS_norm_cat, -log2(utr3_vs_cds.norm.mean),fill=localization_cat))+geom_violin(scale="area",na.rm=T,position=position_dodge())+
    geom_boxplot(width=0.1,na.rm=T, position=position_dodge(width=0.9), outlier.shape = NA)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+scale_fill_manual(values=c("dodgerblue2","orange3"))
  



#for the fig2a, get zeros for positional information (computationally intensive!)
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

newdat$tramscript<-gsub("_.*","",newdat$tc_id)
newdat$pos_from_start<-as.numeric(gsub("_.*","",gsub(".*;","",sub("_",";",newdat$tc_id))))
newdat$pos_from_stop<-as.numeric(gsub("_.*","",gsub(".*;","",sub("_",";",sub("_",";",newdat$tc_id)))))
newdat$region<-gsub(".*_","",newdat$tc_id)

total<-aggregate(cbind(norm.tc_num1, norm.tc_num2)~tramscript, data=newdat, sum)

thr<-total[which(total$norm.tc_num1>=5 & total$norm.tc_num2>=5),]

tz<-which(newdat$tramscript %in% thr$tramscript)

newdat<-newdat[tz,]

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



#fig2a

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





