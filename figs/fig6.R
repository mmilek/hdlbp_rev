# fig6a

library(ggplot2)
library(matrixTests)

fil<-list.files("data/codons_stalling", pattern="txt", full.names = T)
dat<-lapply(fil, read.delim)
names(dat)<-gsub("\\..*", "", fil)

no_norm<-dat[grepl("no_norm",names(dat))]
norm<-dat[grepl("e_norm",names(dat))]

wt<-norm[grepl("wt",names(norm))]
g1<-norm[grepl("guide1",names(norm))]
g2<-norm[grepl("guide2",names(norm))]

wt_g1<-Map(function(x,y) y$plot_value - x$plot_value, wt, g1)
wt_g2<-Map(function(x,y) y$plot_value - x$plot_value, wt, g2)

wt_g<-Map(cbind, wt, g1, g2, wt_g1=wt_g1, wt_g2=wt_g2)
wt_g<-lapply(wt_g, function(x) x[,c(1,2,6,12,18:20)])

r1<-wt_g[grepl("_1_", names(wt_g))]
r2<-wt_g[grepl("_2_", names(wt_g))]

reps<-Map(cbind, r1, r2)
reps<-lapply(reps, function(x) x[,c(1:7,10:14)])

nam<-c("codon", "aa", "wt_1", "g1_1", "g2_1", 
       "diff_wt_g1_1", "diff_wt_g2_1", "wt_2", "g1_2", "g2_2",
       "diff_wt_g1_2", "diff_wt_g2_2")
reps<-lapply(reps, setNames, nam)

avg_g1<-lapply(reps, function(x) rowMeans(x[,c("diff_wt_g1_1", "diff_wt_g1_2")], na.rm=T))
sd_g1<-lapply(reps, function(x) apply(x[,c("diff_wt_g1_1", "diff_wt_g1_2")], 1,sd, na.rm=T))
avg_g2<-lapply(reps, function(x) rowMeans(x[,c("diff_wt_g2_1", "diff_wt_g2_2")], na.rm=T))
sd_g2<-lapply(reps, function(x) apply(x[,c("diff_wt_g2_1", "diff_wt_g2_2")], 1,sd, na.rm=T))

pval_g1<-lapply(reps, function(x) row_t_equalvar(x[,c("wt_1", "wt_2")],x[,c("g1_1", "g1_2")])$pvalue)
pval_g2<-lapply(reps, function(x) row_t_equalvar(x[,c("wt_1", "wt_2")],x[,c("g2_1", "g2_2")])$pvalue)

padj_g1<-lapply(pval_g1, p.adjust, "fdr")
padj_g2<-lapply(pval_g2, p.adjust, "fdr")

reps<-Map(cbind,reps, avg_g1=avg_g1, sd_g1=sd_g1,
          avg_g2=avg_g2, sd_g2=sd_g2, 
          padj_g1= padj_g1,
          padj_g2= padj_g2)

comp<-do.call("rbind", reps)

comp$site<-gsub("_.*", "",gsub(".*;", "",sub("_",";",sub("_",";", gsub(".*\\/", "", row.names(comp))))))

df<-subset(comp, site=="psite" )

# fig6a
ggplot(df, aes(codon, avg_g1))+
  facet_grid(~aa, space="free_x", scales="free_x", switch="x")+
  scale_colour_discrete(guide = "none")+ 
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.2))+
  geom_point(aes(colour=aa))+geom_errorbar(aes(ymin=avg_g1-sd_g1, ymax=avg_g1+sd_g1))+
  ylab("mean codon shift (KO1-WT)")

ggplot(df, aes(codon, -log10(padj_g1)))+
  facet_grid(~aa, space="free_x", scales="free_x", switch="x")+
  scale_colour_discrete(guide = "none")+ 
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.2))+
  geom_point(aes(colour=aa))+geom_hline(yintercept=1, lty=2)+
  ylab("-log10(Padj)")+coord_cartesian(ylim=c(0,2))

df<-subset(comp, site=="psite" )
ggplot(df, aes(codon, avg_g2))+
  facet_grid(~aa, space="free_x", scales="free_x", switch="x")+
  scale_colour_discrete(guide = "none")+ 
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.2))+
  geom_point(aes(colour=aa))+geom_errorbar(aes(ymin=avg_g2-sd_g2, ymax=avg_g2+sd_g2))+
  ylab("mean codon shift (KO2-WT)")

ggplot(df, aes(codon, padj_g2))+
  facet_grid(~aa, space="free_x", scales="free_x", switch="x")+
  scale_colour_discrete(guide = "none")+ 
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.2))+
  geom_point(aes(colour=aa))+geom_hline(yintercept=1, lty=2)+
  ylab("-log10(Padj)")+coord_cartesian(ylim=c(0,2))

df<-subset(comp, site=="asite" )
ggplot(df, aes(codon, avg_g1))+
  facet_grid(~aa, space="free_x", scales="free_x", switch="x")+
  scale_colour_discrete(guide = "none")+ 
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.2))+
  geom_point(aes(colour=aa))+geom_errorbar(aes(ymin=avg_g1-sd_g1, ymax=avg_g1+sd_g1))+
  ylab("mean codon shift (KO1-WT)")

ggplot(df, aes(codon, padj_g1))+
  facet_grid(~aa, space="free_x", scales="free_x", switch="x")+
  scale_colour_discrete(guide = "none")+ 
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.2))+
  geom_point(aes(colour=aa))+geom_hline(yintercept=1, lty=2)+
  ylab("-log10(Padj)")+coord_cartesian(ylim=c(0,2))

df<-subset(comp, site=="asite" )
ggplot(df, aes(codon, avg_g2))+
  facet_grid(~aa, space="free_x", scales="free_x", switch="x")+
  scale_colour_discrete(guide = "none")+ 
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.2))+
  geom_point(aes(colour=aa))+geom_errorbar(aes(ymin=avg_g2-sd_g2, ymax=avg_g2+sd_g2))+
  ylab("mean codon shift (KO2-WT)")

ggplot(df, aes(codon, padj_g2))+
  facet_grid(~aa, space="free_x", scales="free_x", switch="x")+
  scale_colour_discrete(guide = "none")+ 
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.2))+
  geom_point(aes(colour=aa))+geom_hline(yintercept=1, lty=2)+
  ylab("-log10(Padj)")+coord_cartesian(ylim=c(0,2))

# fig6a
df<-subset(comp, site=="esite" )
ggplot(df, aes(codon, avg_g1))+
  facet_grid(~aa, space="free_x", scales="free_x", switch="x")+
  scale_colour_discrete(guide = "none")+ 
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.2))+
  geom_point(aes(colour=aa))+geom_errorbar(aes(ymin=avg_g1-sd_g1, ymax=avg_g1+sd_g1))+
  ylab("mean codon shift (KO1-WT)")

ggplot(df, aes(codon, padj_g1))+
  facet_grid(~aa, space="free_x", scales="free_x", switch="x")+
  scale_colour_discrete(guide = "none")+ 
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.2))+
  geom_point(aes(colour=aa))+geom_hline(yintercept=1, lty=2)+
  ylab("-log10(Padj)")+coord_cartesian(ylim=c(0,2))

df<-subset(comp, site=="esite" )
ggplot(df, aes(codon, avg_g2))+
  facet_grid(~aa, space="free_x", scales="free_x", switch="x")+
  scale_colour_discrete(guide = "none")+ 
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.2))+
  geom_point(aes(colour=aa))+geom_errorbar(aes(ymin=avg_g2-sd_g2, ymax=avg_g2+sd_g2))+
  ylab("mean codon shift (KO2-WT)")

ggplot(df, aes(codon, padj_g2))+
  facet_grid(~aa, space="free_x", scales="free_x", switch="x")+
  scale_colour_discrete(guide = "none")+ 
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.2))+
  geom_point(aes(colour=aa))+geom_hline(yintercept=1, lty=2)+
  ylab("-log10(Padj)")+coord_cartesian(ylim=c(0,2))

sd<-subset(comp, site=="psite"| site=="esite")
# write.table(sd, "source_data/fig6a.txt", quote=F, sep="\t",row.names = F)





#fig6b
library(ggplot2)

fin<-read.delim("geo_processed_data/processed_data_parclip_trna.txt", header=T)

tcs<-subset(fin, select=c("trna_id","exp_normTCpertrna_hdlbp_mean"))
tcs<-subset(tcs, !duplicated(trna_id))
tcs<-tcs[order(tcs$exp_normTCpertrna_hdlbp_mean, decreasing = T),]
tcs$enrichment_rank<-seq(1,nrow(tcs))
tcs<-tcs[,c(1,3)]
fin<-merge(fin, tcs, by="trna_id")

nc<-read.delim("geo_processed_data/processed_data_ncrna_trna.txt", header=T)
nc$expression_cutoff<-rowMeans(nc[,c("norm.ncrnaseq1_read_count","norm.ncrnaseq2_read_count")])

fin<-merge(fin, nc, by="trna_id")

fin$plot<-as.numeric(fin$enrichment_rank)

fin<-subset(fin, expression_cutoff>=1000)

fin$isotype<-gsub(".*tRNA-","", fin$trna_id)

ggplot(fin, aes(factor(plot),fracTCnumber_hdlbp_mean,color=isotype))+
  geom_errorbar(aes(ymin=fracTCnumber_hdlbp_mean-sd_fracTCpos_hdlbp,ymax=fracTCnumber_hdlbp_mean+sd_fracTCpos_hdlbp))+
  geom_point(size=0.7)+
  theme(axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1))+
  geom_hline(yintercept = 0.3, linetype=2, colour="grey")+
  scale_x_discrete(labels=unique(fin$trna_id)[order(unique(fin$enrichment_rank))])

usg<-subset(fin, select=c("trna_id", "enrichment_rank", "codon_usage","expression_cutoff", "exp_normTCpertrna_hdlbp_mean"))
usg<-subset(usg, !duplicated(trna_id))
usg$trna_id<-factor(usg$trna_id, levels=usg$trna_id[order(usg$enrichment_rank)])

usg$dum<-1
ggplot(usg, aes(trna_id, dum))+geom_tile(aes(fill=codon_usage))+
  theme(axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1))
ggplot(usg, aes(trna_id, dum))+geom_tile(aes(fill=log2(expression_cutoff)))+
  theme(axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1))
ggplot(usg, aes(trna_id, dum))+geom_tile(aes(fill=log2(exp_normTCpertrna_hdlbp_mean)))+
  theme(axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1))

sd<-subset(fin, select=c( "trna_id", "seq_tc_middle_t",             
                          "normTCnumber_hdlbp2", "normTCnumber_hdlbp3",         
                           "fracTCnumber_hdlbp_mean", "sd_fracTCpos_hdlbp"  ,        
                          "exp_normTCpertrna_hdlbp2", "exp_normTCpertrna_hdlbp3",    
                           "exp_normTCpertrna_hdlbp_mean","enrichment_rank" ,            
                           "codon_usage","expression_cutoff", "isotype"))
# write.table(sd, "source_data/fig6b.txt", quote=F, sep="\t", row.names = F)
