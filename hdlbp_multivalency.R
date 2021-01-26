##multivalency


##7mer crosslinked
mer<-dat

mer$count<-ifelse(is.na(mer$localization_cat), NA,
                  ifelse(mer$norm_tc_num1>2 & mer$norm_tc_num2>2 & mer$TPM_transcript>=10 , 1, NA))
mer<-subset(mer, !is.na(count))

#subset for testing
# gr1 <-
#   GRanges(seqnames = mer$transcript_id[1:12], 
#           ranges = IRanges(start=mer$tc_start[1:12], 
#                            end=mer$tc_start[1:12] ),
#           mcols=mer[1:12,c("norm_tc_num1", "norm_tc_num2")])
# as.data.frame(findOverlapPairs(gr1, gr1, maxgap = 29L, select="all"))

gr1 <-
  GRanges(seqnames = mer$transcript_id, 
          ranges = IRanges(start=mer$tc_start, 
                           end=mer$tc_start ),
          mcols=mer[,c("norm_tc_num1", "norm_tc_num2", "localization_cat", "tc_region", "l_tr")])
dfgr<-as.data.frame(findOverlapPairs(gr1, gr1, maxgap = 39L, select="all"))
dfgr$distance<-dfgr$first.X.start-dfgr$second.X.start
max(dfgr$distance)
min(dfgr$distance)

dfgr$tc_seq_pos1_start<-ifelse(dfgr$second.X.start-2<1, NA, dfgr$second.X.start-2)
dfgr$tc_seq_pos1_stop<-ifelse(dfgr$second.X.start+1>dfgr$second.mcols.l_tr, NA, dfgr$second.X.start+1)
dfgr$tc_seq_pos2_start<-ifelse(dfgr$second.X.start-1<1, NA, dfgr$second.X.start-1)
dfgr$tc_seq_pos2_stop<-ifelse(dfgr$second.X.start+2>dfgr$second.mcols.l_tr, NA, dfgr$second.X.start+2)
dfgr$tc_seq_pos3_start<-ifelse(dfgr$second.X.start<1, NA, dfgr$second.X.start)
dfgr$tc_seq_pos3_stop<-ifelse(dfgr$second.X.start+3>dfgr$second.mcols.l_tr, NA, dfgr$second.X.start+3)
dfgr$tc_seq_pos4_start<-ifelse(dfgr$second.X.start+1<1, NA, dfgr$second.X.start+1)
dfgr$tc_seq_pos4_stop<-ifelse(dfgr$second.X.start+4>dfgr$second.mcols.l_tr, NA, dfgr$second.X.start+4)


dfgr$tc_seq_pos1<-toupper(as.character(Biostrings::subseq(sub[as.character(dfgr$second.X.seqnames)],
                                                          start=dfgr$tc_seq_pos1_start,end=dfgr$tc_seq_pos1_stop)))
dfgr$tc_seq_pos2<-toupper(as.character(Biostrings::subseq(sub[as.character(dfgr$second.X.seqnames)],
                                                          start=dfgr$tc_seq_pos2_start,end=dfgr$tc_seq_pos2_stop)))
dfgr$tc_seq_pos3<-toupper(as.character(Biostrings::subseq(sub[as.character(dfgr$second.X.seqnames)],
                                                          start=dfgr$tc_seq_pos3_start,end=dfgr$tc_seq_pos3_stop)))
dfgr$tc_seq_pos4<-toupper(as.character(Biostrings::subseq(sub[as.character(dfgr$second.X.seqnames)],
                                                          start=dfgr$tc_seq_pos4_start,end=dfgr$tc_seq_pos4_stop)))


pfgr<-melt(dfgr, measure.vars = colnames(dfgr)[grepl("pos[0-9]$",colnames(dfgr))],
           id.vars = c("second.X.seqnames","second.mcols.norm_tc_num1", "second.mcols.norm_tc_num2", "second.mcols.localization_cat", "second.mcols.tc_region", "distance"),
           value.name = "kmer", variable.name = "kmer_pos")
#write.table(pfgr, "E:/work/hdlbp/multivalency/pfgr_variable.txt", quote=F, sep="\t", row.names=F)

library(ggplot2)
library(reshape2)
setwd("E:/work/hdlbp/multivalency/")
pfgr<-read.delim("pfgr_variable.txt", header=T)
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

ggplot(subset(freq_hs_plot, second.mcols.tc_region!="utr5"), aes(kmer,freq, fill=second.mcols.tc_region))+
  geom_bar(stat="identity", position="dodge")+facet_wrap(~second.mcols.localization_cat)+coord_flip()+scale_fill_manual(values=c("dodgerblue4", "orange3"))


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

ggplot(sub_pfgr , aes(hval_cat, log2(rowMeans(sub_pfgr[,c("second.mcols.norm_tc_num1","second.mcols.norm_tc_num2")])), fill=hval_cat))+geom_violin(scale = "count")+geom_boxplot(width=0.1, outlier.shape = NA) +theme(axis.text.x = element_text(angle = 45, hjust=1))+
  scale_fill_brewer(palette=4)

head(sub_pfgr)



##most crosslinked 4mers


# pfgr<-melt(dfgr, measure.vars = colnames(dfgr)[grepl("pos[0-9]$",colnames(dfgr))],
           # id.vars = c("second.X.seqnames","second.mcols.norm_tc_num1", "second.mcols.norm_tc_num2", "second.mcols.localization_cat", "second.mcols.tc_region", "distance"),
           # value.name = "kmer", variable.name = "kmer_pos")
pfgr<-read.delim("pfgr_variable.txt", header=T)
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

# xlinked<-c("TTCT","CTTC","TCTT","TTCC","TCCT","CTCT","TTTC","CTTT","TCTC","TTTT") #most crosslinked kmers
xlinked<-unique(freq_hs_plot$kmer)
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


ggplot(subset(plot, ((tc_region=="cds" & localization_cat=="membrane") | (tc_region=="utr3" & localization_cat=="cytosolic") )),
       aes(distance, tc, colour=hval_cat))+geom_line()+facet_wrap(~localization_cat+tc_region+hval_cat, scales="free_y", ncol=5)
ggplot(subset(plot, tc_region!="utr5"),aes(distance, tc, colour=hval_cat))+geom_line()+facet_wrap(~localization_cat+tc_region+hval_cat, scales="free_y", ncol=5)

ggplot(subset(sub_pfgr, ((second.mcols.tc_region=="cds" & second.mcols.localization_cat=="membrane") | (second.mcols.tc_region=="utr3" & second.mcols.localization_cat=="cytosolic") )),
       aes(factor(distance), log2(second.mcols.norm_tc_num1), colour=hval_cat))+geom_boxplot()+facet_wrap(~second.mcols.localization_cat+second.mcols.tc_region+hval_cat, ncol=5)


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


ggplot(subset(plot, ((tc_region=="cds" & localization_cat=="membrane") | (tc_region=="utr3" & localization_cat=="cytosolic") )),
       aes(distance, tc, colour=hval_cat))+geom_line()+facet_wrap(~localization_cat+tc_region+hval_cat, scales="free_y", ncol=5)
ggplot(subset(plot, tc_region!="utr5"),aes(distance, tc, colour=hval_cat))+geom_line()+facet_wrap(~localization_cat+tc_region+hval_cat, scales="free_y", ncol=5)

ggplot(subset(sub_pfgr, ((second.mcols.tc_region=="cds" & second.mcols.localization_cat=="membrane") | (second.mcols.tc_region=="utr3" & second.mcols.localization_cat=="cytosolic") )),
       aes(factor(distance), log2(second.mcols.norm_tc_num1), colour=hval_cat))+geom_boxplot()+facet_wrap(~second.mcols.localization_cat+second.mcols.tc_region+hval_cat, ncol=5)





####multivalency in whole sequences

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


seqs<-seqCds_cyt
seqs<-seqs[width(seqs)>=30]
seqs<-seqs[sample(length(seqs),1000)] #sample 1000 sequences
len<-1:max(width(seqs))

tes<-sapply(width(seqs), function(x) len[len<=(x-29)])
f1<-function(y, i, seqs){
  mtes<-y[[i]]
  mseqs<-seqs[i]
  sapply(mtes, function(x) subseq(mseqs, start=x, width=30))
}

res<-lapply(1:length(tes), f1, y=tes, seqs=seqs)

kmerGroup<-DNAStringSet(c("TCTT", "TTCT", "TTTC",  
                          "CTTT", "CCTT", "CTCT", "CTTC",
                          "TCTC", "TCCT", "TTCC",
                          "CCCT", "CTCC", "CCTC", "TCCC"))

cnt<-lapply(lapply(unlist(res), function(x) vcountPDict(kmerGroup, x)), sum)
df<-data.frame(hval=sapply(cnt, '['), window=unlist(tes), transcript=rep(names(seqs), sapply(tes, max)))

setwd("F:/landthaler/HDLBP/all_clip_data/reclip/mapping_trans/")

# write.table(df, "cdsMem_hval_window.txt", quote=F, sep="\t", row.names = F)
# write.table(df, "utr3Mem_hval_window.txt", quote=F, sep="\t", row.names = F)
# write.table(df, "utr3Cyt_hval_window.txt", quote=F, sep="\t", row.names = F)
# write.table(df, "cdsCyt_hval_window.txt", quote=F, sep="\t", row.names = F)


# setwd("F:/landthaler/HDLBP/all_clip_data/reclip/mapping_trans/")
# setwd("E:/work/hdlbp/multivalency/")
setwd("E:/landthaler/HDLBP/all_clip_data/reclip/mapping_trans/")

cdsMem<-read.delim("cdsMem_hval_window.txt", header=T)
utr3Mem<-read.delim("utr3Mem_hval_window.txt", header=T)
cdsCyt<-read.delim("cdsCyt_hval_window.txt", header=T)
utr3Cyt<-read.delim("utr3Cyt_hval_window.txt", header=T)

cdsMem_bck<-read.delim("bck_cdsMem_hval_window.txt", header=T)
utr3Mem_bck<-read.delim("bck_utr3Mem_hval_window.txt", header=T)
cdsCyt_bck<-read.delim("bck_cdsCyt_hval_window.txt", header=T)
utr3Cyt_bck<-read.delim("bck_utr3Cyt_hval_window.txt", header=T)


cdsMem$region<-"cds"
cdsMem$loc<-"membrane"
utr3Mem$region<-"utr3"
utr3Mem$loc<-"membrane"
cdsCyt$region<-"cds"
cdsCyt$loc<-"cytosolic"
utr3Cyt$region<-"utr3"
utr3Cyt$loc<-"cytosolic"


cdsMem_bck$region<-"cds_bck"
cdsMem_bck$loc<-"membrane_bck"
utr3Mem_bck$region<-"utr3_bck"
utr3Mem_bck$loc<-"membrane_bck"
cdsCyt_bck$region<-"cds_bck"
cdsCyt_bck$loc<-"cytosolic_bck"
utr3Cyt_bck$region<-"utr3_bck"
utr3Cyt_bck$loc<-"cytosolic_bck"

hval<-rbind(cdsMem, utr3Mem, cdsCyt, utr3Cyt, cdsMem_bck, cdsCyt_bck, utr3Cyt_bck, utr3Mem_bck)
# hval<-rbind(cdsMem, utr3Mem, cdsCyt, utr3Cyt)
hval$region_loc<-paste0(hval$region,"_", hval$loc)

thr<-3
sub_hval<-subset(hval, hval>=thr)
tot<-aggregate(hval~region+loc+region_loc, data=sub_hval, sum, na.rm=T)
sums<-aggregate(hval~transcript+region+loc+region_loc, data=sub_hval, sum, na.rm=T)
sums$frac_tot<-sums$hval/sum(hval$hval)
ggplot(sums, aes(frac_tot, colour=region_loc))+geom_density()
ggplot(sums, aes(frac_tot, colour=loc))+stat_ecdf()
ggplot(subset(sums, grepl("cds", region_loc)), aes(frac_tot, colour=region_loc))+stat_ecdf()
ggplot(subset(sums, grepl("utr3", region_loc)), aes(frac_tot, colour=region_loc))+stat_ecdf()
ggplot(subset(sums, !grepl("bck", region_loc)), aes(frac_tot, colour=region_loc))+stat_ecdf()
ggplot(subset(sums, !grepl("bck", region_loc)), aes(frac_tot, colour=region_loc))+geom_density()

ggplot(sums, aes(frac_tot, colour=loc))+geom_density()+coord_cartesian(xlim=c(0,0.0002))
ggplot(sums, aes(frac_tot, colour=loc))+stat_ecdf()

wilcox.test(subset(sums, loc=="membrane", select="frac_tot")[,1],
            subset(sums, loc=="cytosolic", select="frac_tot")[,1])
wilcox.test(subset(sums, loc=="membrane_bck", select="frac_tot")[,1],
            subset(sums, loc=="cytosolic_bck", select="frac_tot")[,1])
wilcox.test(subset(sums, region_loc=="cds_membrane", select="frac_tot")[,1],
            subset(sums, region_loc=="cds_cytosolic", select="frac_tot")[,1])
wilcox.test(subset(sums, region_loc=="cds_bck_membrane_bck", select="frac_tot")[,1],
            subset(sums, region_loc=="cds_bck_cytosolic_bck", select="frac_tot")[,1])

ggplot(sums, aes(region_loc,frac_tot, colour=region_loc))+geom_boxplot()
ggplot(sums, aes(loc,log2(frac_tot), colour=loc))+geom_boxplot()

ggplot(subset(sums, grepl("cds", region_loc)), aes(region_loc,frac_tot, colour=region_loc))+geom_boxplot()



avg<-aggregate(hval~transcript+region+loc+region_loc, data=hval, mean, na.rm=T)
maxs<-aggregate(hval~transcript+region+loc+region_loc, data=hval, max, na.rm=T)

# avg<-aggregate(hval~transcript+region+loc+region_loc, data=sub_hval, mean, na.rm=T)


library(ggplot2)

ggplot(avg, aes(hval, colour=loc))+geom_density()+facet_wrap(~region)
ggplot(subset(avg,grepl("cds",region_loc)), aes(hval, colour=region_loc))+stat_ecdf()
ggplot(subset(avg,grepl("bck",region_loc)), aes(hval, colour=region_loc))+stat_ecdf()
ggplot(avg, aes(hval, colour=loc))+stat_ecdf()

ggplot(avg, aes(hval, colour=loc))+geom_density()
ggplot(avg, aes(loc,hval, colour=loc))+geom_boxplot()
ggplot(avg, aes(region_loc,hval, colour=region_loc))+geom_boxplot()

ggplot(maxs, aes(hval, colour=loc))+geom_density()+facet_wrap(~region)


ggplot(avg, aes(hval, colour=loc))+geom_density()+facet_wrap(~region)

ggplot(avg, aes(hval, colour=region_loc))+geom_density()
ggplot(maxs, aes(hval, colour=region_loc))+geom_density()

quants<-quantile(avg$hval, c(0,1/5, 2/5, 3/5, 4/5, 1))

avg$hval_cat<-ifelse(avg$hval>quants[5], "hval_hi", 
                          ifelse(avg$hval>quants[4], "hval_midhi", 
                                 ifelse(avg$hval>quants[3], "hval_midmid",
                                        ifelse(avg$hval>quants[2], "hval_midlo",
                                               ifelse(avg$hval>=0, "hval_lo", NA)))))
ggplot(avg, aes(hval_cat))+geom_bar()

avg$count<-1
classAvg<-aggregate(count~hval_cat+loc+region+region_loc, data=avg, sum, na.rm=T)
sumAvg<-aggregate(count~loc+region+region_loc, data=avg, sum, na.rm=T)
classAvg<-merge(classAvg, sumAvg, by="region_loc")

classAvg$hval_cat<-factor(classAvg$hval_cat, levels=c("hval_hi", "hval_midhi", "hval_midmid", "hval_midlo", "hval_lo"))

ggplot(classAvg, aes(hval_cat,count.x/count.y, fill=region_loc))+geom_bar(stat="identity", position="dodge")
ggplot(classAvg, aes(hval_cat,count.x/count.y, fill=loc.x))+geom_bar(stat="identity", position="dodge")

head(classAvg)


quants<-quantile(hval$hval, c(0,1/5, 2/5, 3/5, 4/5, 1))

hval$hval_cat<-ifelse(hval$hval>quants[5], "hval_hi", 
                     ifelse(hval$hval>quants[4], "hval_midhi", 
                            ifelse(hval$hval>quants[3], "hval_midlo",
                                   ifelse(hval$hval<=quants[3], "hval_lo",NA))))
                                          
ggplot(hval, aes(hval_cat))+geom_bar()
avg<-aggregate(hval~transcript+region+loc+region_loc+hval_cat, data=hval, mean, na.rm=T)
ggplot(avg, aes(hval_cat,log2(hval), fill=region_loc))+geom_boxplot()
ggplot(avg, aes(log2(hval), colour=region_loc))+geom_density()+facet_wrap(~hval_cat, scales="free_y")+coord_cartesian(xlim=c(0,3))
ggplot(avg, aes(log2(hval), colour=region_loc))+geom_density()+coord_cartesian(xlim=c(0,3))

avg$count<-1
classAvg<-aggregate(count~hval_cat+loc+region+region_loc, data=avg, sum, na.rm=T)
classAvg$hval_cat<-factor(classAvg$hval_cat, levels=c("hval_hi", "hval_midhi", "hval_midmid", "hval_midlo", "hval_lo"))
ggplot(classAvg, aes(hval_cat,count, fill=region_loc))+geom_bar(stat="identity", position="dodge")

head(hval)

ggplot(hval, aes(hval, fill=region_loc))+geom_bar(position="dodge")


threshold<-0
avg<-aggregate(hval~transcript+region+loc+region_loc, 
               data=subset(hval, hval>=threshold), mean, na.rm=T)

ggplot(avg, aes(hval, colour=region_loc))+geom_density()
ggplot(avg, aes(hval, colour=loc))+geom_density()


quants<-quantile(avg$hval, c(0,1/5, 2/5, 3/5, 4/5, 1))

avg$hval_cat<-ifelse(avg$hval>quants[5], "hval_hi", 
                     ifelse(avg$hval>quants[4], "hval_midhi", 
                            ifelse(avg$hval>quants[3], "hval_midmid",
                                   ifelse(avg$hval>quants[2], "hval_midlo",
                                          ifelse(avg$hval>=0, "hval_lo", NA)))))

avg$count<-1
classAvg<-aggregate(count~hval_cat+loc+region+region_loc, data=avg, sum, na.rm=T)
classAvg$hval_cat<-factor(classAvg$hval_cat, levels=c("hval_hi", "hval_midhi", "hval_midmid", "hval_midlo", "hval_lo"))
ggplot(classAvg, aes(hval_cat))+geom_bar()
ggplot(classAvg, aes(hval_cat,count, fill=region_loc))+geom_bar(stat="identity", position="dodge")
cumAvg<-aggregate(count~region_loc, data=classAvg, sum)
ggplot(cumAvg, aes(region_loc,count, fill=region_loc))+geom_bar(stat="identity", position="dodge")


ggplot(classAvg, aes(hval_cat,count, fill=loc))+geom_bar(stat="identity", position="dodge")
ggplot(classAvg, aes(hval_cat,count, fill=region))+geom_bar(stat="identity", position="dodge")


avg$count<-1
classAvg<-aggregate(count~hval_cat+loc+region+region_loc, data=avg, function(x) sum(x, na.rm=T)/sum(avg$count, na.rm=T))
classAvg$hval_cat<-factor(classAvg$hval_cat, levels=c("hval_hi", "hval_midhi", "hval_midmid", "hval_midlo", "hval_lo"))

ggplot(classAvg, aes(hval_cat,count, fill=region_loc))+geom_bar(stat="identity", position="dodge")



hval<-rbind(cdsMem, utr3Mem, cdsCyt, utr3Cyt)
hval$region_loc<-paste0(hval$region,"_", hval$loc)

hval<-subset(hval, hval>=3)

avg<-aggregate(hval~transcript+loc+region_loc, data=hval, sum, na.rm=T)
ggplot(avg, aes(hval, colour=region_loc))+geom_density()

regs<-aggregate(hval~transcript+loc, data=hval, sum, na.rm=T)
ggplot(regs, aes(hval, colour=loc))+geom_density()

### per gene function
setwd("E:/work/hdlbp/multivalency/")
tab<-read.delim("hval_perGene.txt", header=T)

setwd("E:/Google Drive/hdlbp/")
# mas<-read.delim("hdlbp_master_table_with_classes.txt", header=T)
mas<-read.delim("hdlbp_master_table_with_classes_uniq_tsig.txt", header=T)
inf<-subset(mas, select=c("gene_id","Symbol","transcript",
                          "tpm_cutoff","tc_CDS_norm","localization_cat",
                          "mean_te_293","loc_tar_CDS","tc_CDS_norm_cat",
                          "tc_transcript_norm",
                          "log2FoldChange.mem.cyt.293",
                          "log2FoldChange.mem.cyt.KO.293",
                          "log2FoldChange.ribo.rna.KO.WT"))

fus<-merge(mas, tab, by="transcript", all.x=T)
ggplot(fus, aes(hval))+geom_histogram()

fus$hval_cat<-ifelse(is.na(fus$hval), NA, 
                     ifelse(fus$hval >50, ">50",
                            ifelse(fus$hval >10, ">10","<10")))

ggplot(fus, aes(hval_cat))+geom_bar()
ggplot(fus, aes(log2FoldChange.mem.cyt.KO.293, colour=hval_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))

ggplot(fus, aes(log2FoldChange.ribo.rna.KO.WT, colour=hval_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))

# setwd("E:/landthaler/HDLBP/ribo/meta/considered")
# trans<-read.delim("considered_gene_names.txt", header=T)
# subset(trans, gene_name=="TFRC")
