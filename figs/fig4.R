library(corrplot)
library(ggplot2)

#fig4

#fig4d
dat<-read.delim("data/proteinGroups - Exp150+144.txt", header=T)

use<-dat[,c(2,7,11,12,53:57,65,66,42:46)]
colnames(use)[5:9]<-c("lfq.hdlbp1", "lfq.ctrl2", "lfq.hdlbp2", "lfq.ctrl3", "lfq.hdlbp3")

colnames(use)[12:16]<-c("int.hdlbp1", "int.ctrl2", "int.hdlbp2", "int.ctrl3", "int.hdlbp3")

use$Gene.names<-as.character(use$Gene.names)

#compute enrichment per replicate

use$enr1<-ifelse(use$lfq.ctrl2==0, NA, use$lfq.hdlbp1/use$lfq.ctrl2)
use$enr2<-ifelse(use$lfq.ctrl2==0, NA, use$lfq.hdlbp2/use$lfq.ctrl2)
use$enr3<-ifelse(use$lfq.ctrl3==0, NA, use$lfq.hdlbp3/use$lfq.ctrl3)


corrplot(cor(use[,17:19], use="pairwise.complete.obs"), method="color", addCoef.col = "white")
corrplot(cor(use[,5:9], use="pairwise.complete.obs"), method="color", addCoef.col = "white")

pairs(~enr1+enr2+enr3, log="xy", cex=0.5, lower.panel=NULL, data=use)

use$gene<-gsub(";.*","",use$Gene.names)

use$positive<-ifelse(use$gene=="HDLBP", "HDLBP", NA)

ggplot(subset(use, Razor...unique.peptides>=5), aes(log2(lfq.hdlbp1), log2(lfq.ctrl2), colour=positive))+geom_point()+geom_text(aes(label=gene))
ggplot(subset(use, Razor...unique.peptides>=5), aes(log2(lfq.hdlbp2), log2(lfq.ctrl2), colour=positive))+geom_point()+geom_text(aes(label=gene))
ggplot(subset(use, Razor...unique.peptides>=5), aes(log2(lfq.hdlbp3), log2(lfq.ctrl3), colour=positive))+geom_point()+geom_text(aes(label=gene))

thr<-3
use$up1<-ifelse( (use$enr1>=thr& log2(use$lfq.hdlbp1)>=27) | (use$lfq.ctrl2==0 & log2(use$lfq.hdlbp1)>=27),"enr1", NA)
use$up2<-ifelse( (use$enr2>=thr & log2(use$lfq.hdlbp2)>=27) | (use$lfq.ctrl2==0 & log2(use$lfq.hdlbp2)>=27),"enr2", NA)
use$up3<-ifelse( (use$enr3>=thr & log2(use$lfq.hdlbp3)>=27) | (use$lfq.ctrl3==0 & log2(use$lfq.hdlbp3)>=27),"enr3", NA)

#figS4e

ggplot(subset(use, Razor...unique.peptides>=3), aes(log10(lfq.hdlbp1), log10(lfq.ctrl2), colour=up1))+geom_point(shape=1)
ggplot(subset(use, Razor...unique.peptides>=3), aes(log10(lfq.hdlbp2), log10(lfq.ctrl2), colour=up2))+geom_point(shape=1)
ggplot(subset(use, Razor...unique.peptides>=3), aes(log10(lfq.hdlbp3), log10(lfq.ctrl3), colour=up3))+geom_point(shape=1)

use$up<-ifelse(use$up1=="enr1" & use$up2=="enr2" & use$up3=="enr3", "enr", NA)
use$up_labels<-ifelse(use$up=="enr", use$gene,NA)
nrow(subset(use, up=="enr"))

use$mean_enr<-rowMeans(use[,c("enr1", "enr2", "enr3")])
use$mean_lfq<-rowMeans(use[,c("lfq.hdlbp1","lfq.hdlbp2","lfq.hdlbp3")] )
use$imp_enr<-ifelse(is.na(use$mean_enr), max(use$mean_enr, na.rm=T)+2, use$mean_enr)
use<-use[order(use$mean_lfq, decreasing = T),]
use<-subset(use, !duplicated(gene))
use$imp_enr_cat<-ifelse(is.na(use$mean_enr), "imputed", NA)

ggplot(subset(use, Razor...unique.peptides>=5), aes(log10(mean_lfq), log2(mean_enr), colour=up))+geom_point()
ggplot(subset(use, Razor...unique.peptides>=5), aes(log10(mean_lfq), log2(imp_enr), colour=up))+geom_point()

ggplot(subset(use, Razor...unique.peptides>=5), aes(log10(mean_lfq), log2(mean_enr), colour=up))+geom_point()+geom_text(aes(label=up_labels))
ggplot(subset(use, Razor...unique.peptides>=5), aes(log10(mean_lfq), log2(imp_enr), colour=up))+geom_point()+geom_text(aes(label=up_labels))

nrow(subset(use, up=="enr" & Reverse!="+" & Potential.contaminant!="+"& Razor...unique.peptides>=3, select=gene))
subset(use, up=="enr" & Reverse!="+" & Potential.contaminant!="+"& Razor...unique.peptides>=3, select=gene)
int<-subset(use, up=="enr" & Reverse!="+" & Potential.contaminant!="+"& Razor...unique.peptides>=3 )
# write.table(int, "list_top_bioid.txt", quote=F, sep="\t", row.names=F, col.names=T)

int<-int[order(int$mean_lfq, decreasing = T),]
top<-int[1:60,]
ggplot(top,aes(factor(gene, levels=top$gene[order(top$mean_lfq, decreasing = F)]), log10(mean_lfq)))+geom_bar(stat = "identity")+coord_flip()+xlab("")

#fig4d
ggplot(top)+geom_point(aes(factor(gene, levels=top$gene[order(top$mean_lfq, decreasing = F)]), 1, colour=log10(mean_lfq), size=log2(imp_enr)))+coord_flip()+xlab("")+scale_colour_gradient(low="dodgerblue2", high="orange2",limits=c(6,13))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

sd<-subset(top, !is.na(mean_lfq),
           select=c("gene", "mean_lfq", "imp_enr" , "imp_enr_cat" ))

# write.table(sd, "source_data/fig4d.txt", quote=F, sep="\t", row.names = F)




