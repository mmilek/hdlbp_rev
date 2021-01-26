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

# setwd("~/Google Drive/hdlbp/")
setwd("E:/Google Drive/hdlbp/")

mas<-read.delim("hdlbp_master_table_with_classes_uniq.txt", header=T)
inf<-subset(mas, select=c("gene_id","Symbol",
                          "tpm_cutoff","tc_CDS_norm","localization_cat",
                          "mean_te_293","loc_tar_CDS","tc_CDS_norm_cat",
                          "tc_transcript_norm",
                          "log2FoldChange.mem.cyt.293",
                          "log2FoldChange.mem.cyt.KO.293",
                          "log2FoldChange.ribo.rna.KO.WT",
                          "Transmembrane.helices" ,"Cleavage.site..Signalp." ))


# setwd("/Volumes/landthaler/pcp/projects/miha/HDLBP/ribo/meta/considered/")
setwd("D:/landthaler/HDLBP/ribo/meta/considered/")

trans<-read.delim("considered_gene_names.txt", header=T)[,c(1:10,16,18)]
colnames(trans)[10]<-"gene_id"

inf<-merge(inf, trans, by="gene_id")


# setwd("/Volumes/landthaler/pcp/projects/miha/HDLBP/ribo/meta/considered/")
setwd("D:/landthaler/HDLBP/ribo/meta/considered/")

fin<-read.delim("psite_position_counts_chx.txt", header=T)
# fin<-read.delim("psite_position_counts_nochx.txt", header=T)

colnames(fin)[2:3]<-c("293_1", "293_2")



tsig<-read.delim("signalp_tm_positions.txt", header=T) # download this from http://grch37.ensembl.org/biomart/martview careful with id versions
genc<-read.delim("gencode_v19_gene_id_to_gene_name_all.txt", header=F)
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
bedi<-merge(anH, tsig, by.x="transcript", by.y="Transcript.stable.ID.version")
sps<-subset(bedi, !is.na(bedi$Cleavage.site..Signalp..start))
tms<-subset(bedi, !is.na(bedi$Transmembrane.helices.start))

sps$start<-sps$l_utr5+3*(sps$Cleavage.site..Signalp..start-1)
sps$stop<-sps$l_utr5+3*(sps$Cleavage.site..Signalp..end-1)
sps$num<-1

tms$start<-tms$l_utr5+3*(tms$Transmembrane.helices.start-1)
tms$stop<-tms$l_utr5+3*(tms$Transmembrane.helices.end-1)
tms$num<-1

bedi<-rbind(sps,tms)
#write.table(bedi[,c("transcript", "start","stop","tsig")], "tsig.bed", quote=F, sep="\t", row.names=F, col.names = F)

#continue
dms<-aggregate(start~ Transcript.stable.ID.version+Gene.stable.ID.version+tsig, data=tsig, min)

# dms<-subset(dms, !duplicated(Gene.stable.ID.version))

ggplot(dms, aes(start))+geom_histogram(bins=1000)+facet_wrap(~tsig)+coord_cartesian(xlim=c(-10,100))

glock<-which(dms$Transcript.stable.ID.version %in% inf$transcript)

tarsig<-dms[glock,]

ggplot(tarsig, aes(tsig))+geom_bar()

tsig_annot<-merge(inf, tarsig, by.x="gene_id", by.y="Gene.stable.ID.version", all.x=T)
tsig_annot$tsig<-ifelse(is.na(tsig_annot$tsig) & tsig_annot$localization_cat=="membrane", "mem_notsig", 
                        ifelse(is.na(tsig_annot$tsig) & tsig_annot$localization_cat=="cytosolic", "cyt_notsig",tsig_annot$tsig))
ggplot(tsig_annot, aes(tsig))+geom_bar()

nrow(subset(tsig_annot, as.character(transcript)==as.character(Transcript.stable.ID.version)))

man<-subset(tsig_annot, !is.na(tsig))

man$nucstart<-(man$start-1)*3
man$minusstart<-ifelse(man$nucstart==0,0, -man$nucstart)
man$plusstart<-ifelse(man$l_cds>=1502, 1502+man$minusstart,man$l_cds+man$minusstart )

man$dist_stop<-man$l_cds-man$nucstart

man$tsig<-ifelse(man$dist_stop<=150 & man$tsig=="TMhelix-only" , "TailAnchored", man$tsig)
man$tsig<-ifelse(grepl("MT-", man$Symbol), "MitoEncoded", man$tsig)
man$tsig<-ifelse(grepl("SignalP", man$tsig), "SignalP", man$tsig)
# setwd("~/Google Drive/hdlbp/")
setwd("E:/Google Drive/hdlbp/")

mito<-read.delim("mitocarta2.txt", header=T)
man$ens_gene_id<-gsub("\\..*","",man$gene_id)
stud<-merge(man, mito, by.x="ens_gene_id", by.y="gene_id", all.x=T)
stud$tsig<-ifelse(is.na(stud$mito), stud$tsig,
                  ifelse(stud$tsig=="cyt_notsig" & stud$mito=="MitoCarta", "MitoCarta", stud$tsig))

man<-stud

ggplot(man, aes(tsig))+geom_bar()
ggplot(man[man$tpm_cutoff>=10,],aes(tsig,log2(tc_CDS_norm),fill=tsig))+geom_boxplot()+geom_point(aes(fill = tsig), size = 1, shape = 1, position = position_jitterdodge())
ggplot(man[man$tpm_cutoff>=10,],aes(tsig,log2(tc_CDS_norm),fill=tsig))+geom_violin(scale="area",na.rm=T)+geom_boxplot(fill="white", na.rm=T, width=0.3)

ggplot(man[man$tpm_cutoff>=10,],aes(tsig,log2(tc_CDS_norm),fill=tsig))+geom_boxplot()+geom_point(aes(fill = tsig), size = 1, shape = 1, position = position_jitterdodge())

ggplot(man[man$tpm_cutoff>=10,],aes(log2FoldChange.mem.cyt.293,log2(mean_te_293),colour=tsig))+geom_point()



#continue
man<-subset(man, tsig!="TailAnchored")
man<-subset(man, !is.na(minusstart) & !is.na(plusstart))
seqs<-man[,c("transcript", "minusstart","plusstart")]
seqs<-subset(seqs, !is.na(minusstart) & !is.na(plusstart))
lseq<-apply(seqs, 1, function(x) seq(x[2], x[3]))

dseq<-data.frame(transcript=rep(unique(man$transcript), man$plusstart-man$minusstart+1), pos_from_tsig=unlist(lseq))
dseq<-subset(dseq, pos_from_tsig>=(-300))


dseq$psite_id_tsig<-paste0(dseq$transcript,"_",dseq$pos_from_tsig)


df<-merge(fin, tarsig, by.x="transcript", by.y="Transcript.stable.ID.version")

df$pos_from_tsig<-df$pos_from_start-(df$start-1)*3
df$psite_id_tsig<-paste0(df$transcript,"_",df$pos_from_tsig)
df$frame_start<-ifelse(df$pos_from_tsig%%3==0, df$pos_from_tsig/3,
                       ifelse(df$pos_from_tsig%%3==1, ((df$pos_from_tsig)-1)/3, ((df$pos_from_tsig)-2)/3))



tust<-merge(dseq, df, by="psite_id_tsig", all.x=T)
tust$pos_from_tsig<-as.numeric(gsub(".*_","",tust$psite_id_tsig))
tust$frame_start<-ifelse(tust$pos_from_tsig%%3==0, tust$pos_from_tsig/3,
                       ifelse(tust$pos_from_tsig%%3==1, ((tust$pos_from_tsig)-1)/3, ((tust$pos_from_tsig)-2)/3))



tust<-subset(tust, select=c("psite_id_tsig","transcript.x","pos_from_tsig","frame_start","293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2"))
values<-tust[,c("293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2")]
values[is.na(values)]<-0
tust<-cbind(tust[,1:4],values)
tust$id<-paste0(tust$transcript.x,"_",tust$frame_start)

codons<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~id+transcript.x, data=tust, sum)
codons$frame_start<-as.numeric(gsub(".*_","",codons$id))

excluded<-which(codons$frame_start<=40 & codons$frame_start>=20)
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

ggplot(mel, aes(codon, value, colour=variable))+geom_line()+coord_cartesian(xlim=c(-0,500) )
ggplot(mel[mel$codon>2,], aes(codon, value, colour=variable))+geom_line()+coord_cartesian(xlim=c(0,500 ))


##tm or sp
tm<-subset(tsig_annot, tsig=="TMhelix-only" & tpm_cutoff>=10 & localization_cat=="membrane")
sp<-subset(tsig_annot, grepl("SignalP", tsig_annot$tsig) & tpm_cutoff>=10 & localization_cat=="membrane" & !is.na(loc_tar_CDS) & loc_tar_CDS!="nontarget_membrane")
# tm<-subset(tsig_annot, tsig=="TMhelix-only" & tpm_cutoff>=10 )
# sp<-subset(tsig_annot, grepl("SignalP", tsig_annot$tsig) & tpm_cutoff>=10)

tm<-which(nor$transcript %in% tm$transcript)
sp<-which(nor$transcript %in% sp$transcript)

tm<-nor[tm,]
sp<-nor[sp,]

length(unique(tm$transcript))
length(unique(sp$transcript))

tm$localization<-"tm"
sp$localization<-"sp"

locn<-rbind(tm, sp)

avg<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~codon+localization, data=locn, mean, na.rm=T)
mel<-melt(avg, measure.vars = c("293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2"), id.vars = c("codon","localization"))
library(zoo)
ggplot(mel, aes(codon, value, colour=localization))+geom_line(aes(y=rollmean(value, 5, na.pad=T)))+facet_wrap(~variable)


ggplot(mel[ mel$localization=="sp",], aes(codon, value, colour=localization))+geom_line(aes(y=rollmean(value, 5, na.pad=T)))+facet_wrap(~variable)+
  geom_hline(yintercept = 2.5, lty=2)+geom_vline(xintercept = 40, lty=2)+coord_cartesian(ylim=c(0,5),xlim=c(0,250))
ggplot(mel[ mel$localization=="tm",], aes(codon, value, colour=localization))+geom_line(aes(y=rollmean(value, 5, na.pad=T)))+facet_wrap(~variable)+
  geom_hline(yintercept = 2, lty=2)+geom_vline(xintercept = 40, lty=2)+coord_cartesian(ylim=c(0,4),xlim=c(-100,250))


hi_sp<-subset(tsig_annot, grepl("SignalP",tsig_annot$tsig) & tpm_cutoff>=10 & loc_tar_CDS=="membrane_tc>10.82 & tc<181.46")
mid_sp<-subset(tsig_annot, grepl("SignalP",tsig_annot$tsig) & tpm_cutoff>=10 & loc_tar_CDS=="membrane_tc>2.81 & tc<10.82")
lo_sp<-subset(tsig_annot, grepl("SignalP",tsig_annot$tsig) & tpm_cutoff>=10 & loc_tar_CDS=="membrane_tc<2.81")
hi_tm<-subset(tsig_annot, tsig=="TMhelix-only" & tpm_cutoff>=10 & loc_tar_CDS=="membrane_tc>10.82 & tc<181.46")
mid_tm<-subset(tsig_annot, tsig=="TMhelix-only" & tpm_cutoff>=10 & loc_tar_CDS=="membrane_tc>2.81 & tc<10.82")
lo_tm<-subset(tsig_annot, tsig=="TMhelix-only" & tpm_cutoff>=10 & loc_tar_CDS=="membrane_tc<2.81")

# hi_sp<-subset(tsig_annot, grepl("SignalP",tsig_annot$tsig) & tpm_cutoff>=10 & tc_CDS_norm_cat=="tc>1.48 & tc<181.46")
# mid_sp<-subset(tsig_annot, grepl("SignalP",tsig_annot$tsig) & tpm_cutoff>=10 & tc_CDS_norm_cat=="tc>0.41 & tc<1.48")
# lo_sp<-subset(tsig_annot, grepl("SignalP",tsig_annot$tsig) & tpm_cutoff>=10 & tc_CDS_norm_cat=="tc<0.41")
# hi_tm<-subset(tsig_annot, tsig=="TMhelix-only" & tpm_cutoff>=10 & tc_CDS_norm_cat=="tc>1.48 & tc<181.46")
# mid_tm<-subset(tsig_annot, tsig=="TMhelix-only" & tpm_cutoff>=10 & tc_CDS_norm_cat=="tc>0.41 & tc<1.48")
# lo_tm<-subset(tsig_annot, tsig=="TMhelix-only" & tpm_cutoff>=10 & tc_CDS_norm_cat=="tc<0.41")


hi_sp<-which(nor$transcript %in% hi_sp$transcript)
mid_sp<-which(nor$transcript %in% mid_sp$transcript)
lo_sp<-which(nor$transcript %in% lo_sp$transcript)
hi_tm<-which(nor$transcript %in% hi_tm$transcript)
mid_tm<-which(nor$transcript %in% mid_tm$transcript)
lo_tm<-which(nor$transcript %in% lo_tm$transcript)

hi_sp<-nor[hi_sp,]
mid_sp<-nor[mid_sp,]
lo_sp<-nor[lo_sp,]
hi_tm<-nor[hi_tm,]
mid_tm<-nor[mid_tm,]
lo_tm<-nor[lo_tm,]

length(unique(hi_sp$transcript))
length(unique(mid_sp$transcript))
length(unique(lo_sp$transcript))
length(unique(hi_tm$transcript))
length(unique(mid_tm$transcript))
length(unique(lo_tm$transcript))

hi_sp$localization<-"hi_sp"
mid_sp$localization<-"mid_sp"
lo_sp$localization<-"lo_sp"
hi_tm$localization<-"hi_tm"
mid_tm$localization<-"mid_tm"
lo_tm$localization<-"lo_tm"

locn<-rbind(hi_sp, mid_sp, lo_sp, hi_tm, mid_tm, lo_tm)


avg<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~codon+localization, data=locn, mean, na.rm=T)
mel<-melt(avg, measure.vars = c("293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2"), id.vars = c("codon","localization"))
sp<-subset(mel, grepl("hi_sp", localization)| grepl("mid_sp", localization))
tm<-subset(mel, grepl("hi_tm", localization)| grepl("mid_tm", localization))
ggplot(sp, aes(codon, value, colour=localization))+geom_line(aes(y=rollmean(value,5, na.pad=T)))+facet_wrap(~variable)+coord_cartesian(xlim=c(0,300))+geom_vline(xintercept = 40, lty=2)
ggplot(tm, aes(codon, value, colour=localization))+geom_line(aes(y=rollmean(value, 5, na.pad=T)))+facet_wrap(~variable)+coord_cartesian(xlim=c(-100,300))+geom_vline(xintercept = 40, lty=2)

ggplot(mel[mel$localization=="hi_tm",], aes(codon, value, colour=localization))+geom_line(aes(y=rollmean(value, 5, na.pad=T)))+facet_wrap(~variable)+coord_cartesian(xlim=c(-100,300))+geom_vline(xintercept = 40, lty=2)
ggplot(mel[mel$localization=="lo_tm",], aes(codon, value, colour=localization))+geom_line(aes(y=rollmean(value, 5, na.pad=T)))+facet_wrap(~variable)+coord_cartesian(xlim=c(-100,300))+geom_vline(xintercept = 40, lty=2)
ggplot(mel[mel$localization=="hi_sp",], aes(codon, value, colour=localization))+geom_line(aes(y=rollmean(value, 5, na.pad=T)))+facet_wrap(~variable)+coord_cartesian(xlim=c(0,300))+geom_vline(xintercept = 40, lty=2)

ggplot(sp, aes(codon, value, colour=variable))+geom_line(aes(y=rollmean(value, 5, na.pad=T)))+facet_wrap(~variable)+coord_cartesian(xlim=c(0,300))+geom_vline(xintercept = 40, lty=2)

ggplot(tm, aes(codon, value, colour=variable))+geom_line(aes(y=rollmean(value, 5, na.pad=T)))+facet_wrap(~variable)+coord_cartesian(xlim=c(-100,300))+geom_hline(yintercept = 2, lty=2)+geom_vline(xintercept = 40, lty=2)

ggplot(mel[ mel$localization=="tm",], aes(codon, value, colour=localization))+geom_line(aes(y=rollmean(value, 5, na.pad=T)))+facet_wrap(~variable)+geom_hline(yintercept = 2, lty=2)+geom_vline(xintercept = 40, lty=2)


ggplot(mel[grepl("_sp", mel$localization),], aes(codon, value, colour=variable))+geom_line(aes(y=rollmean(value,5, na.pad=T)))+facet_wrap(~variable)+coord_cartesian(xlim=c(0,100))
ggplot(mel[grepl("_tm", mel$localization),], aes(codon, value, colour=variable))+geom_line(aes(y=rollmean(value,5, na.pad=T)))+facet_wrap(~variable)+coord_cartesian(xlim=c(-100,100))

ggplot(mel[grepl("293_1", mel$variable) & grepl("_sp", mel$localization),], aes(codon, value, colour=localization))+geom_line(aes(y=rollmean(value,11, na.pad=T)))
ggplot(mel[grepl("293_1", mel$variable) & grepl("_tm", mel$localization),], aes(codon, value, colour=localization))+geom_line(aes(y=rollmean(value,11, na.pad=T)))+coord_cartesian(xlim=c(-100,300))

ggplot(mel[grepl("_sp", mel$localization),], aes(codon, value, colour=variable))+geom_line(aes(y=rollmean(value,11, na.pad=T)))+facet_wrap(~variable)

ggplot(mel[grepl("_sp", mel$localization),], aes(codon, value, colour=variable))+geom_line(aes(y=rollmean(value,50, na.pad=T)))
ggplot(mel[ grepl("_sp", mel$localization),], aes(codon, value, colour=variable))+geom_smooth()
ggplot(mel[ grepl("_tm", mel$localization),], aes(codon, value, colour=variable))+geom_smooth()


