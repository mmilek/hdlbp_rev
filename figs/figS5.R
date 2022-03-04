library(reshape2)
library(ggplot2)
library(zoo)

#figS5d

dat<-read.delim("data/hdlbp_master_table_with_classes_uniq.txt", header=T)

dat$localization_cat<-as.factor(dat$localization_cat)

dat$labels<-ifelse(dat$Symbol=="HSPA5"|dat$Symbol=="REEP4"|dat$Symbol=="CHST14",
                   dat$Symbol, NA)

ggplot(subset(dat, tpm_cutoff>=10 & Annotation=="protein_coding"),
       aes(log2FoldChange.mem.cyt.293,log2FoldChange.mem.cyt.2A15))+
  geom_point(shape=1)+
  geom_text(aes(label=labels))
  

sd<-subset(dat, tpm_cutoff>=10 & Annotation=="protein_coding" & !is.na(log2FoldChange.mem.cyt.293) & !is.na(log2FoldChange.mem.cyt.2A15),
           select=c("gene_id", "Symbol", "tpm_cutoff", "log2FoldChange.mem.cyt.293", "log2FoldChange.mem.cyt.2A15", "labels"))

# write.table(sd, "source_data/figS5d.txt", quote=F, sep="\t", row.names=F, col.names=T)

# figS5e

PG_summ<-read.delim("data/protein_groups_psilac.txt", header=T)
man<-read.delim("data/hdlbp_master_table_with_classes_uniq_tsig.txt", header=T)
pas<-merge(man, PG_summ, by.x="Symbol", by.y="Gene.names", all.x=T)

pas$mean_l_mem<-rowMeans(pas[,colnames(pas)[grepl("iBAQ.L.Memb", colnames(pas))]])

pas$Ratio.H.M.normalized.Memb_4h_Reverse_1<--pas$Ratio.H.M.normalized.Memb_4h_Reverse_1
pas$Ratio.H.M.normalized.Memb_4h_Reverse_2<--pas$Ratio.H.M.normalized.Memb_4h_Reverse_2

pas$mean_memb_4h<-rowMeans(pas[,colnames(pas)[grepl("Ratio.H.M.normalized.Memb_4h", colnames(pas))]])

pas$target<-ifelse(is.na(pas$tc_transcript_norm_cat), NA,
                   ifelse(pas$tc_transcript_norm_cat=="nontarget", "nontarget", "target"))
pas$target_loc<-ifelse(is.na(pas$target), NA,
                       ifelse(pas$target=="nontarget", "nontarget", 
                              ifelse(pas$target=="target" & pas$localization_cat=="cytosolic", "cyto_target",
                                     ifelse(pas$target=="target" & pas$localization_cat=="membrane", "mem_target", NA))))



int<-subset(pas, mean_memb_4h<0 & target=="target" & log2FoldChange.ribo.rna.KO.WT<0 & localization_cat=="membrane")
int<-subset(int, !duplicated(Symbol))
int$Symbol<-factor(int$Symbol,  levels=int$Symbol[order(int$mean_memb_4h, decreasing = T)])
int<-int[order(int$mean_memb_4h, decreasing = F),]
int<-int[1:60,]

#figS5e
ggplot(int) +geom_point( aes(Symbol, 1, colour=-mean_memb_4h, size= -log2FoldChange.ribo.rna.KO.WT))+coord_flip()+xlab("")+scale_colour_gradient(low="dodgerblue2", high="orange2")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
ggplot(int) +geom_point( aes(Symbol, 1, colour=-mean_memb_4h, size= log2(tc_transcript_norm)))+coord_flip()+xlab("")+scale_colour_gradient(low="dodgerblue2", high="orange2")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
ggplot(int) +geom_point( aes(Symbol, 1, colour=-mean_memb_4h, size= -log2FoldChange.ribo.rna.KO.WT))+coord_flip()+xlab("")+scale_colour_gradient(low="dodgerblue2", high="orange2")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())+facet_wrap(~localization_cat)
ggplot(int) +geom_point( aes(Symbol, 1, colour=-mean_memb_4h, size= log2(tc_transcript_norm)))+coord_flip()+xlab("")+scale_colour_gradient(low="dodgerblue2", high="orange2")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

sd<-subset(int, !is.na(mean_memb_4h),
           select=c("gene_id", "Symbol", "gene_biotype","target_loc", "mean_memb_4h", "log2FoldChange.ribo.rna.KO.WT", "tc_transcript_norm"))
# write.table(sd, "source_data/figS5e.txt", quote=F, sep="\t", row.names=F)

#figS5g begin (same as fig5i)

mas<-read.delim("data/hdlbp_master_table_with_classes_uniq.txt", header=T)
inf<-subset(mas, select=c("gene_id","Symbol",
                          "tpm_cutoff","tc_CDS_norm","localization_cat",
                          "mean_te_293","loc_tar_CDS","tc_CDS_norm_cat",
                          "tc_transcript_norm",
                          "log2FoldChange.mem.cyt.293",
                          "log2FoldChange.mem.cyt.KO.293",
                          "log2FoldChange.ribo.rna.KO.WT",
                          "Transmembrane.helices" ,"Cleavage.site..Signalp." ))

trans<-read.delim("data/considered_gene_names.txt", header=T)[,c(1:10,16,18)]
colnames(trans)[10]<-"gene_id"

inf<-merge(inf, trans, by="gene_id")

fin<-read.delim("data/psite_position_counts_chx.txt", header=T)
colnames(fin)[2:3]<-c("293_1", "293_2")

tsig<-read.delim("data/signalp_tm_positions.txt", header=T) # download this from http://grch37.ensembl.org/biomart/martview careful with id versions
genc<-read.delim("data/gencode_v19_gene_id_to_gene_name_all.txt", header=F)
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

dms<-aggregate(start~ Transcript.stable.ID.version+Gene.stable.ID.version+tsig, data=tsig, min)

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

mito<-read.delim("data/mitocarta2.txt", header=T)
colnames(mito)<-c("gene_id", "mito")
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

tm_list<-subset(tm, select=c("gene_id", "Symbol", "Transcript.stable.ID.version", 
                             "l_utr5", "l_cds", "l_utr3", "start"))
sp_list<-subset(sp, select=c("gene_id", "Symbol", "Transcript.stable.ID.version", 
                             "l_utr5", "l_cds", "l_utr3", "start"))

# make lists of genes for browsing
# write.table(tm_list, "data/tm_list.txt", quote=F, sep="\t", row.names=F)
# write.table(sp_list, "data/sp_list.txt", quote=F, sep="\t", row.names=F)

tm<-which(nor$transcript %in% tm$transcript)
sp<-which(nor$transcript %in% sp$transcript)

tm<-nor[tm,]
sp<-nor[sp,]

length(unique(tm$transcript))
length(unique(sp$transcript))

list1<-data.frame(transcript=unique(sp$transcript), tsig="sp")
list2<-data.frame(transcript=unique(tm$transcript), tsig="tm")
list<-rbind(list1, list2)

# write.table(list, "tm_sp_ribo_tsig_list.txt", quote=F, sep="\t", row.names=F)

tm$localization<-"tm"
sp$localization<-"sp"

locn<-rbind(tm, sp)

avg<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~codon+localization, data=locn, mean, na.rm=T)
mel<-melt(avg, measure.vars = c("293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2"), id.vars = c("codon","localization"))

ggplot(mel, aes(codon, value, colour=localization))+geom_line(aes(y=rollmean(value, 5, na.pad=T)))+facet_wrap(~variable)

ggplot(mel[ mel$localization=="sp",], aes(codon, value, colour=localization))+geom_line(aes(y=rollmean(value, 5, na.pad=T)))+facet_wrap(~variable)+
  geom_hline(yintercept = 2.5, lty=2)+geom_vline(xintercept = 40, lty=2)+coord_cartesian(ylim=c(0,5),xlim=c(0,250))
ggplot(mel[ mel$localization=="tm",], aes(codon, value, colour=localization))+geom_line(aes(y=rollmean(value, 5, na.pad=T)))+facet_wrap(~variable)+
  geom_hline(yintercept = 2, lty=2)+geom_vline(xintercept = 40, lty=2)+coord_cartesian(ylim=c(0,4),xlim=c(-100,250))

# make source data
# write.table(mel, "../../hdlbp_rev/hdlbp_rev/source_data/fig5i.txt", quote=F, sep="\t", row.names = F)

# figS5c this is very computationaly intense! make 
# sure you have enough resources.

fin<-read.delim("data/psite_position_counts_chx.txt", header=T)
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

mask<-merge(codons, normi, by="transcript")
first<-as.matrix(mask[,3:8])
second<-as.matrix(mask[11:16])
nor<-first/second

nor<-cbind(mask[,c(1,2,9,10)], nor)
colnames(nor)<-gsub("\\.x","",colnames(nor))

minFilter<-which(nor$transcript %in% minorf$transcript)
nor<-nor[minFilter,]

avg<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~frame_start, data=nor, mean, na.rm=T)
mel<-melt(avg, measure.vars = c("293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2"), id.vars = "frame_start")
ggplot(mel[mel$frame_start>=2 & mel$frame_start<=500,], aes(frame_start, value, colour=variable))+geom_line()+coord_cartesian(xlim=c(0,300))
avg<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~frame_stop, data=nor, mean, na.rm=T)
mel<-melt(avg, measure.vars = c("293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2"), id.vars = "frame_stop")
ggplot(mel[mel$frame_stop<=(-2) & mel$frame_stop>=(-500),], aes(frame_stop, value, colour=variable))+geom_line()+coord_cartesian(xlim=c(-300,0))

tsig_annot<-read.delim("data/tsig_annot.txt", header=T)
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

# figS5c plotting

avg<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~frame_start+localization, data=locn, mean, na.rm=T)
mel<-melt(avg, measure.vars = c("293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2"), id.vars = c("frame_start","localization"))
mel<-mel[mel$frame_start>=2 & mel$frame_start<=500 & mel$localization!="memN" ,]
# write.table(mel, "source_data/figS5c_start.txt", quote=F, sep="\t", row.names=F)

ggplot(mel[mel$frame_start>=2 & mel$frame_start<=500 & mel$localization!="memN" ,], aes(frame_start, value, colour=localization))+geom_line(aes(y=rollmean(value, 10, na.pad=T)))+facet_wrap(~variable)+coord_cartesian(xlim=c(0,250))+geom_vline(xintercept=39, lty=2)+geom_vline(xintercept=79, lty=2)+ylab("scaled_psite_coverage")
ggplot(mel[mel$frame_start>=2 & mel$frame_start<=500 & mel$localization!="memN",], aes(frame_start, value, colour=localization))+geom_line(aes(y=rollmean(value, 10, na.pad=T)))+facet_wrap(~variable)+coord_cartesian(xlim=c(0,500))+geom_vline(xintercept=39, lty=2)+geom_vline(xintercept=79, lty=2)+ylab("scaled_psite_coverage")

avg<-aggregate(cbind(`293_1`,`293_2` ,guide1_1 ,guide1_2 ,guide2_1, guide2_2)~frame_stop+localization, data=locn, mean, na.rm=T)
mel<-melt(avg, measure.vars = c("293_1","293_2" ,"guide1_1" ,"guide1_2" ,"guide2_1", "guide2_2"), id.vars = c("frame_stop","localization"))
mel<-mel[mel$frame_stop<=(-2) & mel$frame_stop>=(-500) & mel$localization!="memN" ,]
# write.table(mel, "source_data/figS5c_stop.txt", quote=F, sep="\t", row.names=F)

ggplot(mel[mel$frame_stop<=(-2) & mel$frame_stop>=(-500) & mel$localization!="memN" ,], aes(frame_stop, value, colour=localization))+geom_line(aes(y=rollmean(value, 10, na.pad=T)))+facet_wrap(~variable)+coord_cartesian(xlim=c(-250,0))+ylab("scaled_psite_coverage")
ggplot(mel[mel$frame_stop<=(-2) & mel$frame_stop>=(-500) & mel$localization!="memN" ,], aes(frame_stop, value, colour=localization))+geom_line(aes(y=rollmean(value, 10, na.pad=T)))+facet_wrap(~variable)+coord_cartesian(xlim=c(-500,0))+ylab("scaled_psite_coverage")

# figS5a and figS5b - file size limitation, not all
# files provided

library(riboWaltz)
library(data.table)
library(ggplot2)

# prepare annotation for hg19
# anH<-create_annotation(gtfpath = "gencode.v19.primary_assembly.annotation.gtf", dataSource="gencode.v19", organism="Homo sapiens")
# sequences_biost <- Biostrings::readDNAStringSet("hg19.transcripts.fa",
# length(names(sequences_biost) %in% anH$transcript)
# nrow(anH[names(sequences_biost) %in% anH$transcript,])
# length(sequences_biost[as.character(anH$transcript)])
# save.image("ribowaltz_hg38.RData")

# load("ribowaltz_human.RData") # due to size limitations this cannot be provided, but you can build your own annotation as indicated above
anH<-data.table(anH) #higher versions of ribowaltz work with data table

#here transcriptome-mapped read bed files have to be provided (bamToBed)
# setwd("D:/Documents/hdlbp_git/bed_original/")
# 
# files<-list.files(getwd(), pattern="*bed*")
files<-gsub("\\.bed","",files)

#set output dir
# out<-"D:/Documents/hdlbp_git/original_riboQC/"

reads <- bedtolist(getwd(),anH)

psite_offset <- psite(reads, start=TRUE, cl=100, flanking = 6, extremity="auto", 
                      plot=F, plot_dir=paste0(out,"/psite_offset_start"), plot_format="pdf", log_file = F, log_file_dir = paste0(out,"/psite_offset_start/"))

reads_psite_list <- psite_info(reads, psite_offset)

#figS5a
for (i in names(reads)){
  comparison_dt<-list()
  comparison_dt[[paste0(i,"_21nt")]]<-reads_psite_list[[i]][length==21]
  comparison_dt[[paste0(i,"_22nt")]]<-reads_psite_list[[i]][length==22]
  comparison_dt[[paste0(i,"_23nt")]]<-reads_psite_list[[i]][length==23]
  comparison_dt[[paste0(i,"_24nt")]]<-reads_psite_list[[i]][length==24]
  comparison_dt[[paste0(i,"_25nt")]]<-reads_psite_list[[i]][length==25]
  comparison_dt[[paste0(i,"_26nt")]]<-reads_psite_list[[i]][length==26]
  comparison_dt[[paste0(i,"_27nt")]]<-reads_psite_list[[i]][length==27]
  comparison_dt[[paste0(i,"_28nt")]]<-reads_psite_list[[i]][length==28]
  comparison_dt[[paste0(i,"_29nt")]]<-reads_psite_list[[i]][length==29]
  comparison_dt[[paste0(i,"_30nt")]]<-reads_psite_list[[i]][length==30]
  comparison_dt[[paste0(i,"_31nt")]]<-reads_psite_list[[i]][length==31]
  comparison_dt[[paste0(i,"_32nt")]]<-reads_psite_list[[i]][length==32]
  names_list<-list(paste0(i,"_21nt"),
                   paste0(i,"_22nt"),
                   paste0(i,"_23nt"), 
                   paste0(i,"_24nt"),
                   paste0(i,"_25nt"),
                   paste0(i,"_26nt"),
                   paste0(i,"_27nt"),
                   paste0(i,"_28nt"),
                   paste0(i,"_29nt"),
                   paste0(i,"_30nt"),
                   paste0(i,"_31nt"),
                   paste0(i,"_32nt"))
  names(names_list)<-c(paste0(i,"_21nt"),
                       paste0(i,"_22nt"), 
                       paste0(i,"_23nt"),
                       paste0(i,"_24nt"),
                       paste0(i,"_25nt"),
                       paste0(i,"_26nt"),
                       paste0(i,"_27nt"),
                       paste0(i,"_28nt"),
                       paste0(i,"_29nt"),
                       paste0(i,"_30nt"),
                       paste0(i,"_31nt"),
                       paste0(i,"_32nt"))
  example_metaheatmap<-metaheatmap_psite(comparison_dt, anH, sample=names_list,
                                         utr5l=20, cdsl=40, utr3l=20, log=T)
  example_metaheatmap[["plot"]]
  dfwt<-example_metaheatmap[["dt"]]
  # write.table(dfwt, paste0(out,"source_data/",i,".txt"), quote=F, sep="\t", row.names=F)
  ggsave(paste0(out,"/metaheatmap/",i,".pdf" ),width = 10, height = 7, units="in")
}

#figS5b
for (i in names(reads)){
  example_frames<-frame_psite(reads_psite_list, sample=i, region="all")
  example_frames[["plot"]]
  dfwt<-example_frames[["dt"]]
  # write.table(dfwt, paste0(out,"frame_psite/source_data/",i,".txt"), quote=F, sep="\t", row.names=F)
  ggsave(paste0(out,"/frame_psite/",i,".pdf" ),width = 6, height = 3, units="in")
}