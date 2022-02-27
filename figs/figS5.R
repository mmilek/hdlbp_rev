#figS5e

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
