#fig6b

fin<-read.delim("geo_processed_data/processed_data_parclip_trna.txt", header=T)

tcs<-subset(fin, select=c("trna_id","exp_normTCpertrna_hdlbp_mean"))
tcs<-subset(tcs, !duplicated(trna_id))
tcs<-tcs[order(tcs$exp_normTCpertrna_hdlbp_mean, decreasing = T),]
tcs$enrichment_rank<-seq(1,nrow(tcs))
tcs<-tcs[,c(1,3)]
fin<-merge(fin, tcs, by="trna_id")

nc<-read.delim("geo_processed_data/processed_data_ncrna_trna.txt", header=T)
nc$expression_cutoff<-rowMeans(nc[,c("norm.ncrnaseq1_read_count","norm.ncrnaseq2_read_count")])
# nc$isotype=gsub(";.*","",fin$id)
# fin$codon<-as.character(reverseComplement(DNAStringSet(gsub(".*;","",fin$id))))
# fin$id<-paste0(fin$id,";",fin$codon)
# fin$id<-paste0(gsub(";.*","-",fin$id),
#                gsub("T","U",gsub(".*-","",sub(";","-",fin$id))))

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
