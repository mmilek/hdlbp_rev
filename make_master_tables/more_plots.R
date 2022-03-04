# more data exploration plots

library(reshape2)
library(ggplot2)

fin<-read.delim("data/hdlbp_master_table_with_classes_uniq.txt", header=T)

##change in membrane localization upon KO

subset(fin, tpm_cutoff>=10 & padj.mem.cyt.2A15<0.1 & padj.mem.cyt.3C2<0.1  & log2FoldChange.mem.cyt.293>1.5 & log2FoldChange.mem.cyt.KO.293<(-0.1) & target=="target", 
       select=c("Symbol", "log2FoldChange.mem.cyt.KO.293", "padj.mem.cyt.2A15", "padj.mem.cyt.3C2"))
subset(fin, tpm_cutoff>=10 & padj.mem.cyt.2A15<0.1 & padj.mem.cyt.3C2<0.1  & log2FoldChange.mem.cyt.293>1.5 & log2FoldChange.mem.cyt.KO.293<(-0.1) & target=="nontarget", 
       select=c("Symbol", "log2FoldChange.mem.cyt.KO.293", "padj.mem.cyt.2A15", "padj.mem.cyt.3C2"))


tot_up<-subset(fin, tpm_cutoff>=1 & padj.tot.KO.WT <0.1 & log2FoldChange.tot.KO.WT > 0.5, 
       select=c("Symbol", "padj.tot.KO.WT", "log2FoldChange.tot.KO.WT", "target"))
tot_down<-subset(fin, tpm_cutoff>=1 & padj.tot.KO.WT <0.1 & log2FoldChange.tot.KO.WT < (-0.5), 
               select=c("Symbol", "padj.tot.KO.WT", "log2FoldChange.tot.KO.WT","target"))
nrow(subset(tot_up, target=="target"))
nrow(subset(tot_down, target=="target"))

#write.table(tot_up, "tot_up.txt", quote=F, sep="\t", row.names=F)
#write.table(tot_down, "tot_down.txt", quote=F, sep="\t", row.names=F)

rel<-subset(fin, tpm_cutoff>=10)

ggplot()+ geom_point(aes(rel$log2FoldChange.mem.cyt.293, rel$log2FoldChange.mem.cyt.KO), colour="grey")+
 geom_point(aes(rel[rel$target=="target","log2FoldChange.mem.cyt.293"], rel[rel$target=="target","log2FoldChange.mem.cyt.KO"]), colour="orange")+geom_abline()


plot(rel$log2FoldChange.mem.cyt.293, rel$log2FoldChange.mem.cyt.KO)
identify(rel$log2FoldChange.mem.cyt.293, rel$log2FoldChange.mem.cyt.KO, labels=rel$Symbol)

subset(fin, tpm_cutoff>=1 & log2FoldChange.mem.cyt.KO.293 < (-.5) & log2FoldChange.mem.cyt.KO > 2 & log2FoldChange.mem.cyt.WT > 2 & target=="target", 
       select=c("Symbol", "log2FoldChange.mem.cyt.KO.293","conv_CDS"))

subset(fin, Symbol=="HSPA5", select=grepl("LFQ.*L.*", colnames(fin)))
subset(fin, Symbol=="HSPA5", select=grepl("LFQ.*H.*", colnames(fin)))
subset(fin, Symbol=="HSPA5", select=grepl("Intensity.*L.*", colnames(fin)))
subset(fin, Symbol=="HSPA5", select=grepl("Intensity.*H.*", colnames(fin)))
subset(fin, Symbol=="HSPA5", select=grepl("zero.lfc.*", colnames(fin)))
subset(fin, Symbol=="HSPA5", select=c("Symbol",colnames(fin)[grepl("LFQ.*L.*", colnames(fin))]))
subset(fin, Symbol=="HSPA5", select=c("Symbol",colnames(fin)[grepl("_.*_.*", colnames(fin))| grepl("mem.cyt", colnames(fin))]))
subset(fin, Symbol=="FZD1", select=c("Symbol",colnames(fin)[grepl("_.*_.*", colnames(fin))| grepl("mem.cyt", colnames(fin))]))
subset(fin, Symbol=="MFSD5", select=c("Symbol",colnames(fin)[grepl("_.*_.*", colnames(fin))| grepl("mem.cyt", colnames(fin))]))
subset(fin, Symbol=="PDIA3", select=c("Symbol",colnames(fin)[grepl("_.*_.*", colnames(fin))| grepl("mem.cyt", colnames(fin))]))


fin$localization<-ifelse(fin$tpm_cutoff>=1 & fin$log2FoldChange.mem.cyt.293>1.5 , "membrane", 
                         ifelse(fin$tpm_cutoff>=1 & fin$log2FoldChange.mem.cyt.293<0.5, "cytosol",
                                ifelse(fin$tpm_cutoff>=1 , "intermediate", NA)))

fin$localization_target<-ifelse(fin$tpm_cutoff>=1 & fin$log2FoldChange.mem.cyt.293>1.5 & fin$target=="target", "membrane_target", 
                                ifelse(fin$tpm_cutoff>=1 & fin$log2FoldChange.mem.cyt.293>1.5 & fin$target=="nontarget", "membrane_nontarget",
                                       ifelse(fin$tpm_cutoff>=1 & fin$log2FoldChange.mem.cyt.293<0.5 & fin$target=="target", "cytosol_target" ,
                                              ifelse(fin$tpm_cutoff>=1 & fin$log2FoldChange.mem.cyt.293<0.5 & fin$target=="nontarget", "cytosol_nontarget", 
                                              NA))))


rel<-subset(fin, !is.na(localization_target) & tpm_cutoff>=10)
mel<-melt(rel, measure.vars = colnames(rel)[grepl("T_",colnames(rel))], id.vars = c("Symbol", "target", "Transmembrane.helices", "localization","localization_target"))
ggplot(mel, aes(variable, log2(value), fill=localization))+geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))
ggplot(mel, aes(variable, log2(value), fill=localization_target))+geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))

#does the membrane pool increase? yes, strongly
ggplot(rel, aes(log2FoldChange.mem.tot.2A15.293, colour=localization_target))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(rel, aes(log2FoldChange.mem.tot.3C2.293, colour=localization_target))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(rel, aes(log2FoldChange.mem.tot.2A15.293, colour=localization))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(rel, aes(log2FoldChange.mem.tot.3C2.293, colour=localization))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))

ggplot(rel, aes(log2FoldChange.mem.KO.WT, colour=localization))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))


wilcox.test(rel[!is.na(rel$log2FoldChange.mem.tot.2A15.293) & rel$localization=="cytosol","log2FoldChange.mem.tot.2A15.293"],
            rel[!is.na(rel$log2FoldChange.mem.tot.2A15.293) & rel$localization=="membrane","log2FoldChange.mem.tot.2A15.293"])
wilcox.test(rel[!is.na(rel$log2FoldChange.mem.tot.3C2.293) & rel$localization=="cytosol","log2FoldChange.mem.tot.3C2.293"],
            rel[!is.na(rel$log2FoldChange.mem.tot.3C2.293) & rel$localization=="membrane","log2FoldChange.mem.tot.3C2.293"])
wilcox.test(rel[!is.na(rel$log2FoldChange.mem.KO.WT) & rel$localization=="cytosol","log2FoldChange.mem.KO.WT"],
            rel[!is.na(rel$log2FoldChange.mem.KO.WT) & rel$localization=="membrane","log2FoldChange.mem.KO.WT"])

#does the total pool stay the same? mem mRNAs are slightly lower but distributions overlap
ggplot(rel, aes(log2FoldChange.tot.KO.WT, colour=localization))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(rel, aes(log2FoldChange.tot.KO.WT, colour=localization_target))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))

wilcox.test(rel[!is.na(rel$log2FoldChange.tot.KO.WT) & rel$localization=="cytosol","log2FoldChange.tot.KO.WT"],
            rel[!is.na(rel$log2FoldChange.mem.tot.2A15.293) & rel$localization=="membrane","log2FoldChange.tot.KO.WT"])


#does the cyto pool decrease? yes, mildly
ggplot(rel, aes(log2FoldChange.tot.cyt.KO.293, colour=localization_target))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(rel, aes(log2FoldChange.tot.cyt.KO.293, colour=localization))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))

wilcox.test(rel[!is.na(rel$log2FoldChange.tot.cyt.KO.293) & rel$localization=="cytosol","log2FoldChange.tot.cyt.KO.293"],
            rel[!is.na(rel$log2FoldChange.tot.cyt.KO.293) & rel$localization=="membrane","log2FoldChange.tot.cyt.KO.293"])

#is there an effect on translation? yes
ggplot(rel, aes(log2FoldChange.rpf.KO.WT, colour=localization_target))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(rel, aes(log2FoldChange.rpf.KO.WT, colour=localization))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))

ggplot(rel, aes(log2FoldChange.ribo.rna.KO.WT, colour=localization_target))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))
ggplot(rel, aes(log2FoldChange.ribo.rna.KO.WT, colour=localization))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))

wilcox.test(rel[!is.na(rel$log2FoldChange.ribo.rna.KO.WT) & rel$localization=="cytosol","log2FoldChange.ribo.rna.KO.WT"],
            rel[!is.na(rel$log2FoldChange.ribo.rna.KO.WT) & rel$localization=="membrane","log2FoldChange.ribo.rna.KO.WT"])

wilcox.test(rel[!is.na(rel$log2FoldChange.rpf.KO.WT) & rel$localization=="cytosol","log2FoldChange.ribo.rna.KO.WT"],
            rel[!is.na(rel$log2FoldChange.rpf.KO.WT) & rel$localization=="membrane","log2FoldChange.ribo.rna.KO.WT"])

wilcox.test(rel[!is.na(rel$log2FoldChange.ribo.rna.KO.WT) & rel$localization_target=="cytosol_nontarget","log2FoldChange.ribo.rna.KO.WT"],
            rel[!is.na(rel$log2FoldChange.ribo.rna.KO.WT) & rel$localization_target=="membrane_target","log2FoldChange.ribo.rna.KO.WT"])
wilcox.test(rel[!is.na(rel$log2FoldChange.ribo.rna.KO.WT) & rel$localization_target=="cytosol_target","log2FoldChange.ribo.rna.KO.WT"],
            rel[!is.na(rel$log2FoldChange.ribo.rna.KO.WT) & rel$localization_target=="membrane_target","log2FoldChange.ribo.rna.KO.WT"])
wilcox.test(rel[!is.na(rel$log2FoldChange.ribo.rna.KO.WT) & rel$localization_target=="cytosol_target","log2FoldChange.ribo.rna.KO.WT"],
            rel[!is.na(rel$log2FoldChange.ribo.rna.KO.WT) & rel$localization_target=="membrane_nontarget","log2FoldChange.ribo.rna.KO.WT"])
wilcox.test(rel[!is.na(rel$log2FoldChange.ribo.rna.KO.WT) & rel$localization_target=="cytosol_nontarget","log2FoldChange.ribo.rna.KO.WT"],
            rel[!is.na(rel$log2FoldChange.ribo.rna.KO.WT) & rel$localization_target=="membrane_nontarget","log2FoldChange.ribo.rna.KO.WT"])
wilcox.test(rel[!is.na(rel$log2FoldChange.ribo.rna.KO.WT) & rel$localization_target=="membrane_target","log2FoldChange.ribo.rna.KO.WT"],
            rel[!is.na(rel$log2FoldChange.ribo.rna.KO.WT) & rel$localization_target=="membrane_nontarget","log2FoldChange.ribo.rna.KO.WT"])
wilcox.test(rel[!is.na(rel$log2FoldChange.ribo.rna.KO.WT) & rel$localization_target=="cytosol_target","log2FoldChange.ribo.rna.KO.WT"],
            rel[!is.na(rel$log2FoldChange.ribo.rna.KO.WT) & rel$localization_target=="cytosol_nontarget","log2FoldChange.ribo.rna.KO.WT"])



#is there an effect in pSilac?
ggplot(rel, aes(lfc.lfq.h.ko.293, colour=localization))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))

wilcox.test(rel[!is.na(rel$lfc.lfq.h.ko.293) & rel$localization=="cytosol","lfc.lfq.h.ko.293"],
            rel[!is.na(rel$lfc.lfq.h.ko.293) & rel$localization=="membrane","lfc.lfq.h.ko.293"])

ggplot()+ geom_point(aes(rel$log2FoldChange.mem.cyt.KO.293, rel$log2FoldChange.ribo.rna.KO.WT), colour="grey")+
  geom_point(aes(rel[rel$target=="target","log2FoldChange.mem.cyt.KO.293"], rel[rel$target=="target","log2FoldChange.ribo.rna.KO.WT"]), colour="orange")+
  coord_cartesian(xlim=c(-3,3))

ggplot()+ geom_point(aes(rel$log2FoldChange.mem.cyt.KO.293, rel$log2FoldChange.rpf.KO.WT), colour="grey")+
  geom_point(aes(rel[rel$target=="target","log2FoldChange.mem.cyt.KO.293"], rel[rel$target=="target","log2FoldChange.rpf.KO.WT"]), colour="orange")

plot(rel$log2FoldChange.mem.cyt.293, rel$log2FoldChange.mem.cyt.KO)
identify(rel$log2FoldChange.mem.cyt.293, rel$log2FoldChange.mem.cyt.KO, labels=rel$Symbol)


#is there an effect on cyto?
ggplot(rel, aes(log2FoldChange.nuc.tot.3C2.293, colour=localization_target))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(rel, aes(log2FoldChange.nuc.tot.2A15.293, colour=localization_target))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(rel, aes(log2FoldChange.nuc.cyt.KO.293, colour=localization_target))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))


#how are the mem.cyt nuc.cyt ratios?

ggplot(rel, aes(log2FoldChange.mem.cyt.KO.293, colour=localization_target))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(rel, aes(log2FoldChange.nuc.cyt.KO.293, colour=localization_target))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(rel, aes(log2FoldChange.mem.KO.WT, colour=localization_target))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(rel, aes(log2FoldChange.tot.cyt.KO.293, colour=localization_target))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(rel, aes(log2FoldChange.tot.KO.WT, colour=localization_target))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(rel, aes(log2FoldChange.nuc.KO.WT, colour=localization_target))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))

##define spike in genes
colnames(fin)[grepl("tot", colnames(fin))]
fin$spike_in<-ifelse(fin$tpm_cutoff>=50 & fin$padj.tot.KO.WT >0.9 & fin$log2FoldChange.tot.KO.WT < 0.1 & fin$log2FoldChange.mem.cyt.293 < 0.5 & 
                       fin$target!="target", "spike_in",NA)
nrow(subset(fin, spike_in=="spike_in"))
rel<-subset(fin,tpm_cutoff>=10)
ggplot()+geom_point(aes(log10(fin$baseMean.tot.KO.WT), fin$log2FoldChange.tot.KO.WT), colour="grey")+
  geom_point(aes(log10(fin[fin$spike_in=="spike_in","baseMean.tot.KO.WT"]), 
                 fin[fin$spike_in=="spike_in","log2FoldChange.tot.KO.WT"]), colour="orange")+
  geom_abline(slope=0)

ggplot()+geom_point(aes(log10(fin$baseMean.mem.KO.WT), fin$log2FoldChange.mem.KO.WT), colour="grey")+
  geom_point(aes(log10(fin[fin$spike_in=="spike_in","baseMean.mem.KO.WT"]), 
                 fin[fin$spike_in=="spike_in","log2FoldChange.mem.KO.WT"]), colour="orange")+
  geom_abline(slope=0)
ggplot()+geom_point(aes(log10(fin$baseMean.cyt.KO.WT), fin$log2FoldChange.cyt.KO.WT), colour="grey")+
  geom_point(aes(log10(fin[fin$spike_in=="spike_in","baseMean.cyt.KO.WT"]), 
                 fin[fin$spike_in=="spike_in","log2FoldChange.cyt.KO.WT"]), colour="orange")+
  geom_abline(slope=0)
ggplot()+geom_point(aes(log10(fin$baseMean.nuc.KO.WT), fin$log2FoldChange.nuc.KO.WT), colour="grey")+
  geom_point(aes(log10(fin[fin$spike_in=="spike_in","baseMean.nuc.KO.WT"]), 
                 fin[fin$spike_in=="spike_in","log2FoldChange.nuc.KO.WT"]), colour="orange")+
  geom_abline(slope=0)
ggplot()+geom_point(aes(log10(fin$baseMean.mem.cyt.KO.293), fin$log2FoldChange.mem.cyt.KO.293), colour="grey")+
  geom_point(aes(log10(fin[fin$spike_in=="spike_in","baseMean.mem.cyt.KO.293"]), 
                 fin[fin$spike_in=="spike_in","log2FoldChange.mem.cyt.KO.293"]), colour="orange")+
  geom_abline(slope=0)

nrow(subset(fin, spike_in=="spike_in"))
nrow(subset(fin, spike_in=="spike_in" & target=="target"))
nrow(subset(fin, spike_in=="spike_in" & localization=="cytosol"))
nrow(subset(fin, spike_in=="spike_in" & localization=="membrane"))

spike<-subset(fin, spike_in=="spike_in", select=c("Symbol", "gene_id"))
#write.table(spike, "spike_in_genes.txt", quote=F, sep="\t", row.names=F)

rel<-subset(fin,tpm_cutoff>=10)
ggplot(rel, aes(lfc.lfq.h.ko.293, colour=target))+stat_ecdf()

ggplot()+geom_point(aes(log10(rel$Intensity.H.08.293_2), rel$lfc.lfq.h.ko.293), colour="grey")+
  geom_point(aes(log10(rel[rel$Transmembrane.helices=="TMhelix","Intensity.H.08.293_2"]), 
                 rel[rel$Transmembrane.helices=="TMhelix","lfc.lfq.h.ko.293"]), colour="orange")+
  geom_abline(slope=0)+ylim(-2,2)

ggplot()+geom_point(aes(log10(rel$Intensity.H.08.293_2), rel$lfc.lfq.h.ko.293), colour="grey")+
  geom_point(aes(log10(rel[rel$Cleavage.site..Signalp.=="SignalP","Intensity.H.08.293_2"]), 
                 rel[rel$Cleavage.site..Signalp.=="SignalP","lfc.lfq.h.ko.293"]), colour="orange")+
  geom_abline(slope=0)+ylim(-2,2)

ggplot()+geom_point(aes(log10(rel$Intensity.H.08.293_2), rel$zero.lfc.lfq.h.ko.293), colour="grey")+
  geom_point(aes(log10(rel[rel$Cleavage.site..Signalp.=="SignalP","Intensity.H.08.293_2"]), 
                 rel[rel$Cleavage.site..Signalp.=="SignalP","zero.lfc.lfq.h.ko.293"]), colour="orange")+
  geom_abline(slope=0)

ggplot()+geom_point(aes(log10(rel$Intensity.H.08.293_2), rel$zero.lfc.lfq.h.ko.293), colour="grey")+
  geom_point(aes(log10(rel[rel$Cleavage.site..Signalp.=="SignalP","Intensity.H.08.293_2"]), 
                 rel[rel$Cleavage.site..Signalp.=="SignalP","zero.lfc.lfq.h.ko.293"]), colour="orange")+
  geom_abline(slope=0)


ggplot()+geom_point(aes(log10(rel$Intensity.H.08.293_2), rel$lfc.lfq.l.ko.293), colour="grey")+
  geom_point(aes(log10(rel[rel$Transmembrane.helices=="TMhelix","Intensity.H.08.293_2"]), 
                 rel[rel$Transmembrane.helices=="TMhelix","lfc.lfq.l.ko.293"]), colour="orange")+
  geom_abline(slope=0)+ylim(-2,2)

ggplot()+geom_point(aes(log10(rel$Intensity.H.08.293_2), rel$lfc.lfq.l.ko.293), colour="grey")+
  geom_point(aes(log10(rel[rel$Cleavage.site..Signalp.=="SignalP","Intensity.H.08.293_2"]), 
                 rel[rel$Cleavage.site..Signalp.=="SignalP","lfc.lfq.l.ko.293"]), colour="orange")+
  geom_abline(slope=0)+ylim(-2,2)

nrow(subset(rel, Cleavage.site..Signalp.=="SignalP" & !is.na(lfc.lfq.l.ko.293)))
nrow(subset(rel, Cleavage.site..Signalp.!="SignalP" & !is.na(lfc.lfq.l.ko.293)))
nrow(subset(rel, (Transmembrane.helices=="TMhelix" | Cleavage.site..Signalp.=="SignalP") & !is.na(lfc.lfq.l.ko.293)))
nrow(subset(rel, Transmembrane.helices!="TMhelix" & Cleavage.site..Signalp.!="SignalP" & !is.na(lfc.lfq.l.ko.293)))
nrow(subset(rel, target=="target" & !is.na(lfc.lfq.l.ko.293)))

nrow(subset(rel, target=="target" & Transmembrane.helices=="TMhelix" & Cleavage.site..Signalp.=="SignalP" & !is.na(lfc.lfq.l.ko.293)))
nrow(subset(rel, target=="target" & Transmembrane.helices=="TMhelix" & Cleavage.site..Signalp.=="SignalP" & !is.na(log2FoldChange.mem.cyt.KO.293)))

nrow(subset(rel, target=="target" & !is.na(log2FoldChange.mem.cyt.KO.293)))

nrow(subset(rel, !is.na(log2FoldChange.mem.cyt.KO.293)))
nrow(subset(rel, Transmembrane.helices=="TMhelix" & Cleavage.site..Signalp.=="SignalP" & !is.na(log2FoldChange.mem.cyt.KO.293)))
nrow(subset(rel, Transmembrane.helices=="TMhelix" & Cleavage.site..Signalp.=="SignalP" & !is.na(log2FoldChange.mem.KO.WT)))
nrow(subset(rel, Transmembrane.helices!="TMhelix" & Cleavage.site..Signalp.!="SignalP" & !is.na(log2FoldChange.mem.KO.WT)))

ggplot()+geom_point(aes(rel$log2FoldChange.mem.cyt.2A15, rel$log2FoldChange.mem.cyt.293 ), colour="grey")+
  geom_point(aes(rel[rel$target=="target","log2FoldChange.mem.cyt.2A15"], rel[rel$target=="target","log2FoldChange.mem.cyt.293"] ),colour="orange")+
  geom_text(label=row.names(rel))+geom_abline(slope=1)

ggplot()+geom_point(aes(rel$log2FoldChange.mem.cyt.2A15, rel$log2FoldChange.mem.cyt.293 ), colour="grey")+
  geom_point(aes(rel[rel$Transmembrane.helices=="TMhelix","log2FoldChange.mem.cyt.2A15"], rel[rel$Transmembrane.helices=="TMhelix","log2FoldChange.mem.cyt.293"] ),colour="orange")+
  geom_text(label=row.names(rel))+geom_abline(slope=1)

nrow(subset(rel, log2FoldChange.mem.cyt.WT>2 & !is.na(log2FoldChange.mem.cyt.KO.293)))
nrow(subset(rel, log2FoldChange.mem.cyt.WT>2 & log2FoldChange.mem.cyt.KO>2 &!is.na(log2FoldChange.mem.cyt.KO.293)))
nrow(subset(rel, Cleavage.site..Signalp.=="SignalP" & Transmembrane.helices=="TMhelix" & log2FoldChange.mem.cyt.WT>2 & log2FoldChange.mem.cyt.KO>2 &!is.na(log2FoldChange.mem.cyt.KO.293)))

nrow(subset(rel, Cleavage.site..Signalp.=="SignalP" & Transmembrane.helices=="TMhelix" & log2FoldChange.mem.cyt.WT<2 & log2FoldChange.mem.cyt.KO<2 &!is.na(log2FoldChange.mem.cyt.KO.293)))
nrow(subset(rel, log2FoldChange.mem.cyt.WT<2 & log2FoldChange.mem.cyt.KO<2 &!is.na(log2FoldChange.mem.cyt.KO.293))) # cytosolic
nrow(subset(rel, log2FoldChange.mem.cyt.WT>2 & log2FoldChange.mem.cyt.KO>2 &!is.na(log2FoldChange.mem.cyt.KO.293))) # membrane

nrow(subset(rel, target=="target" & log2FoldChange.mem.cyt.WT<2 & log2FoldChange.mem.cyt.KO<2 &!is.na(log2FoldChange.mem.cyt.KO.293)))
nrow(subset(rel, target=="target" & log2FoldChange.mem.cyt.WT>2 & log2FoldChange.mem.cyt.KO>2 &!is.na(log2FoldChange.mem.cyt.KO.293)))
nrow(subset(rel, target=="target" & !is.na(target) & !is.na(log2FoldChange.mem.cyt.WT) & !is.na(log2FoldChange.mem.cyt.KO) & !is.na(log2FoldChange.mem.cyt.KO.293)))

nrow(subset(rel, target=="target" &  (Transmembrane.helices=="TMhelix" | Cleavage.site..Signalp.=="SignalP" ) & log2FoldChange.mem.cyt.WT>2 & log2FoldChange.mem.cyt.KO>2 &!is.na(log2FoldChange.mem.cyt.KO.293)))
nrow(subset(rel, target=="target" &  Cleavage.site..Signalp.=="SignalP" & log2FoldChange.mem.cyt.WT>2 & log2FoldChange.mem.cyt.KO>2 &!is.na(log2FoldChange.mem.cyt.KO.293)))


library(reshape2)
library(ggplot2)
rel<-subset(fin, tpm_cutoff>=10)
mel<-melt(rel, measure.vars = colnames(rel[,3:34]), id.vars = c("Symbol", "Annotation","target","Transmembrane.helices", 
                                                                "Cleavage.site..Signalp.","Low.complexity..Seg."))
ggplot(mel, aes(variable, log2(value), fill=Transmembrane.helices))+geom_boxplot()



mel<-melt(rel, measure.vars = colnames(rel)[grepl("log2FoldChange.mem.",colnames(rel))], id.vars = c("Symbol", "Annotation","target","Transmembrane.helices", 
                                                                "Cleavage.site..Signalp.","Low.complexity..Seg."))
ggplot(mel, aes(variable, log2(value), fill=Transmembrane.helices))+geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(mel, aes(variable, log2(value), fill=Cleavage.site..Signalp.))+geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


mel<-melt(rel, measure.vars = colnames(rel)[grepl("log2FoldChange.cyt.",colnames(rel))], id.vars = c("Symbol", "Annotation","target","Transmembrane.helices", 
                                                                                                     "Cleavage.site..Signalp.","Low.complexity..Seg."))
ggplot(mel, aes(variable, log2(value), fill=Transmembrane.helices))+geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(mel, aes(variable, log2(value), fill=Cleavage.site..Signalp.))+geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

mel<-melt(rel, measure.vars = colnames(rel)[grepl("log2FoldChange.mem.cyt.",colnames(rel))|grepl("log2FoldChange.nuc.cyt.",colnames(rel))], id.vars = c("Symbol", "Annotation","target","Transmembrane.helices", 
                                                                                                     "Cleavage.site..Signalp.","Low.complexity..Seg."))
ggplot(mel, aes(variable, log2(value), fill=Transmembrane.helices))+geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(rel, aes(log2FoldChange.mem.cyt.WT, colour=target))+stat_ecdf()
ggplot(rel, aes(log2FoldChange.mem.tot.293, colour=target))+stat_ecdf()
ggplot(rel, aes(log2FoldChange.cyt.tot.293, colour=target))+stat_ecdf()


