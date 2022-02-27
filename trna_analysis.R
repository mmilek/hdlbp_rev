# setwd("E:/work/hdlbp/counts_trna_uniq//")
# on the cluster is on fast /fast/AG_Landthaler/miha/hdlbp/reclip/counts_trna
library(ggplot2)
library(reshape2)
library(DESeq2)
library(corrplot)
mas<-read.delim("data/tRNA_hdlbp2_TCseq_filtered.bed", header=F, stringsAsFactors = F)
# mas<-read.delim("tRNA_hdlbp2_TCseq.bed", header=F)
perTrna<-aggregate(V5~V1, data=mas, sum )
mas<-merge(mas, perTrna, by="V1")
mas$fracTCpos<-mas$V5.x/mas$V8
mas$fracTrna<-mas$V5.x/mas$V5.y
mas<-mas[order(mas$fracTCpos, mas$V5.x, decreasing = T),]
mas$trna<-paste0(mas$V1,";",mas$V2,";",mas$V3)
mas<-mas[,c(12,5,7:11)]
colnames(mas)[2:5]<-c("TCnumber","seq","allReadsPerPos","TCpertrna")
colnames(mas)[2:7]<-paste0(colnames(mas)[2:7],"_hdlbp2")
head(mas, 25)

mas2<-read.delim("data/tRNA_hdlbp3_TCseq_filtered.bed", header=F, stringsAsFactors = F)
# mas2<-read.delim("tRNA_hdlbp3_TCseq.bed", header=F)
perTrna<-aggregate(V5~V1, data=mas2, sum )
mas2<-merge(mas2, perTrna, by="V1")
mas2$fracTCpos<-mas2$V5.x/mas2$V8
mas2$fracTrna<-mas2$V5.x/mas2$V5.y
mas2<-mas2[order(mas2$fracTCpos, mas2$V5.x, decreasing = T),]
mas2$trna<-paste0(mas2$V1,";",mas2$V2,";",mas2$V3)
mas2<-mas2[,c(12,5,7:11)]
colnames(mas2)[2:5]<-c("TCnumber","seq","allReadsPerPos","TCpertrna")
colnames(mas2)[2:7]<-paste0(colnames(mas2)[2:7],"_hdlbp3")

mer<-merge(mas, mas2, by="trna")
mer$id<-gsub("-.*","",sub("-",";",gsub(";.*","",gsub(".*tRNA-","",mer$trna))))
normFactors<-estimateSizeFactorsForMatrix(mer[,c("TCnumber_hdlbp2","TCnumber_hdlbp3")])
mer$normTCnumber_hdlbp2<-mer$TCnumber_hdlbp2/normFactors[1]
mer$normTCnumber_hdlbp3<-mer$TCnumber_hdlbp3/normFactors[2]
mer$normTCnumber_hdlbp_mean<-rowMeans(mer[,c("normTCnumber_hdlbp2","normTCnumber_hdlbp3")])
mer$sd_normTCnumber_hdlbp<-apply(mer[,c("normTCnumber_hdlbp2","normTCnumber_hdlbp3")],1, sd)
mer$fracTCnumber_hdlbp_mean<-rowMeans(mer[,c("fracTCpos_hdlbp2","fracTCpos_hdlbp3")])
mer$sd_fracTCpos_hdlbp<-apply(mer[,c("fracTCpos_hdlbp2","fracTCpos_hdlbp3")],1, sd)



files<-list.files(getwd(), pattern="readcounts$")[-c(5,6)]
dat<-lapply(files, read.delim, header=F)
tet<-Reduce(function(...) merge(..., by="V1",all=T), dat)
row.names(tet)<-gsub("-.*","",sub("-",";",gsub(".*tRNA-","",tet[,1])))
tet<-tet[,-1]
colnames(tet)<-gsub("tRNA_","",files)

norm<-sweep(tet,2,colSums(tet,na.rm=T)/1e6,"/") 
colnames(norm)<-paste0("norm.",colnames(norm))
tet<-cbind(tet, norm)


library(Biostrings)
codUs<-read.delim("codonUsage.txt", header=T)
codUs$anti<-as.character(reverseComplement(DNAStringSet(codUs$codon)))
codUs$trna<-paste0(codUs$aa,";",codUs$anti)

tet<-merge(tet, codUs[,c("trna","codon","codonUsage")],by.x="row.names",by.y="trna",all.x=T)


fin<-merge(mer, tet, by.x="id", by.y="Row.names")

tcs<-subset(fin, select=c("id","TCpertrna_hdlbp2","TCpertrna_hdlbp3"))
tcs<-subset(tcs, !duplicated(id))
norm<-sweep(tcs[,2:3],2,colSums(tcs[,2:3])/1e6,"/") 
colnames(norm)<-paste0("norm",colnames(norm))
tcs<-cbind(id=tcs[,1], norm)

normtcs<-tcs[,2:3]/tet[,11:12]
colnames(normtcs)<-paste0("exp_",colnames(normtcs))
tcs<-cbind(tcs, normtcs)

normtcs<-tcs[,2:3]/tet[,c(15,13)]
colnames(normtcs)<-paste0("exptot_",colnames(normtcs))
tcs<-cbind(tcs, normtcs)

tcs$exp_normTCpertrna_hdlbp_mean<-rowMeans(tcs[,c("exp_normTCpertrna_hdlbp2","exp_normTCpertrna_hdlbp3")])
tcs$exp_normTCpertrna_hdlbp_sd<-apply(tcs[,c("exp_normTCpertrna_hdlbp2","exp_normTCpertrna_hdlbp3")],1, sd)
tcs<-tcs[order(tcs$exp_normTCpertrna_hdlbp_mean, decreasing = T),]
tcs$enrichment_rank<-seq(1,nrow(tcs))
fin<-merge(fin, tcs, by="id")

fin$expression_cutoff<-rowMeans(fin[,c("norm.m1_hi_2ome.readcounts","norm.m2_hi_2ome.readcounts")])
fin$isotype=gsub(";.*","",fin$id)
fin$codon<-as.character(reverseComplement(DNAStringSet(gsub(".*;","",fin$id))))
fin$id<-paste0(fin$id,";",fin$codon)
fin$id<-paste0(gsub(";.*","-",fin$id),
               gsub("T","U",gsub(".*-","",sub(";","-",fin$id))))

# fin$id<-factor(fin$id,levels=fin$enrichment_rank, labels=fin$id[fin$enrichment_rank])
fin$plot<-as.numeric(fin$enrichment_rank)
#write.table(fin, "hdlbp_tRNA_perPosition.txt", quote=F, sep="\t", row.names=F)


fin<-subset(fin, expression_cutoff>=1000)

style<- theme(line = element_line(colour="black",size=1,linetype=1)) + theme(text = element_text(family="Arial",size=12, colour="black")) + theme(panel.background = element_rect(fill="white",colour="black",size=1)) + theme(axis.line = element_blank()) + theme(panel.grid.minor = element_blank()) + theme(panel.grid.major = element_blank())+ theme(axis.text = element_text(colour="black")) 
style<- theme(panel.background = element_blank())
ggplot(fin, aes(factor(plot),fracTCnumber_hdlbp_mean,color=isotype))+
  geom_errorbar(aes(ymin=fracTCnumber_hdlbp_mean-sd_fracTCpos_hdlbp,ymax=fracTCnumber_hdlbp_mean+sd_fracTCpos_hdlbp))+
  geom_point(size=0.7)+
  theme(axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1))+
  geom_hline(yintercept = 0.3, linetype=2, colour="grey")+scale_x_discrete(labels=unique(fin$id)[order(unique(fin$enrichment_rank))])



ggplot(fin, aes(plot,fracTCnumber_hdlbp_mean,color=isotype))+
  geom_text(aes(label=seq_hdlbp2), size=3, angle=90)+
  theme(axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1))+
  geom_hline(yintercept = 0.3, linetype=2, colour="grey")

usg<-subset(fin, select=c("id", "enrichment_rank", "codonUsage","expression_cutoff", "exp_normTCpertrna_hdlbp_mean"))
usg<-subset(usg, !duplicated(id))
usg$id<-factor(usg$id, levels=usg$id[order(usg$enrichment_rank)])

usg$dum<-1
p2<-ggplot(usg, aes(id, dum))+geom_tile(aes(fill=codonUsage))+
  theme(axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1))
p3<-ggplot(usg, aes(id, dum))+geom_tile(aes(fill=log2(expression_cutoff)))+
  theme(axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1))
p4<-ggplot(usg, aes(id, dum))+geom_tile(aes(fill=log2(exp_normTCpertrna_hdlbp_mean)))+
  theme(axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1))


print(p2)
print(p3)
print(p4)





