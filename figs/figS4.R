library(ggplot2)
library(ggpubr)
library(corrplot)
library(VennDiagram)
library(readxl)

#figS4

#figS4d
dat<-read.delim("data/s4d.txt", header=T)
mel<-reshape2::melt(dat[,1:3])
ggbarplot(mel, x = "protein", y = "value",
          ylab= "% total TC transitions", xlab = "protein", color = "protein", fill = "protein",
          add = c("mean_sd"), palette = "jco",
          position = position_dodge(.8)) +
  geom_jitter(aes(protein, value, fill = protein), shape = 21, color = "black",  position = position_jitterdodge(jitter.height = -2, jitter.width = 1.5))


#figS4e
dat<-read.delim("data/proteinGroups - Exp150+144.txt", header=T)

use<-dat[,c(2,7,11,12,53:57,65,66,42:46)]
colnames(use)[5:9]<-c("lfq.hdlbp1", "lfq.ctrl2", "lfq.hdlbp2", "lfq.ctrl3", "lfq.hdlbp3")

colnames(use)[12:16]<-c("int.hdlbp1", "int.ctrl2", "int.hdlbp2", "int.ctrl3", "int.hdlbp3")

use$Gene.names<-as.character(use$Gene.names)
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

sd<-subset(use, Razor...unique.peptides>=3 & Gene.names!="",
           select=c("Gene.names", "Razor...unique.peptides",
                    colnames(use)[grepl("lfq", colnames(use))],
                    "up1", "up2", "up3"))
# write.table(sd, "source_data/figS4e.txt", quote=F, sep="\t", row.names = F)

use$up<-ifelse(use$up1=="enr1" & use$up2=="enr2" & use$up3=="enr3", "enr", NA)
use$up_labels<-ifelse(use$up=="enr", use$gene,NA)
nrow(subset(use, up=="enr"))

use$mean_enr<-rowMeans(use[,c("enr1", "enr2", "enr3")])
use$mean_lfq<-rowMeans(use[,c("lfq.hdlbp1","lfq.hdlbp2","lfq.hdlbp3")] )
use$imp_enr<-ifelse(is.na(use$mean_enr), max(use$mean_enr, na.rm=T)+2, use$mean_enr)
use<-use[order(use$mean_lfq, decreasing = T),]
use<-subset(use, !duplicated(gene))
use$imp_enr_cat<-ifelse(is.na(use$mean_enr), "imputed", NA)


nrow(subset(use, up=="enr" & Reverse!="+" & Potential.contaminant!="+"& Razor...unique.peptides>=3, select=gene))
subset(use, up=="enr" & Reverse!="+" & Potential.contaminant!="+"& Razor...unique.peptides>=3, select=gene)
int<-subset(use, up=="enr" & Reverse!="+" & Potential.contaminant!="+"& Razor...unique.peptides>=3 )
# write.table(int, "list_top_bioid.txt", quote=F, sep="\t", row.names=F, col.names=T)

int<-int[order(int$mean_lfq, decreasing = T),]
top<-int[1:60,]
ggplot(top,aes(factor(gene, levels=gene[order(mean_lfq, decreasing = F)]), log10(mean_lfq)))+geom_bar(stat = "identity")+coord_flip()+xlab("")

int<-int[order(int$mean_enr, decreasing = T),]
top<-int[1:60,]
ggplot(top,aes(factor(gene, levels=gene[order(mean_enr, decreasing = F)]), log2(mean_enr)))+geom_bar(stat = "identity")+coord_flip()

int<-int[order(int$imp_enr, decreasing = T),]
top<-int[1:60,]

top<-subset(int, is.na(mean_enr))
top<-top[order(mean_lfq, decreasing = T),]
top<-top[1:60,]
ggplot(top,aes(factor(gene, levels=gene[order(mean_lfq, decreasing = F)]), log2(mean_lfq)))+geom_bar(stat = "identity")+coord_flip()

int<-int[order(int$Razor...unique.peptides, decreasing = T),]
top<-int[1:60,]
ggplot(top,aes(factor(gene, levels=gene[order(Razor...unique.peptides, decreasing = F)]), log2(Razor...unique.peptides)))+geom_bar(stat = "identity")+coord_flip()


##overlap between replicates

enr1<-subset(use, up1=="enr1")
enr2<-subset(use, up2=="enr2")
enr3<-subset(use, up3=="enr3")

#figS4g

dev.off()
draw.triple.venn(nrow(enr1), nrow(enr2), nrow(enr3),
                 nrow(merge(enr1, enr2, by="gene")),  nrow(merge(enr2, enr3, by="gene")), nrow(merge(enr1, enr3, by="gene")),
                 nrow(merge(enr3, merge(enr1, enr2, by="gene"), by="gene")), euler.d=T,scaled=T, category=c("replicate 1", "replicate 2", "replicate 3"))


enr1<-subset(use, up1=="enr1", select=c("gene","up1"))
enr2<-subset(use, up2=="enr2", select=c("gene","up2"))
enr3<-subset(use, up3=="enr3", select=c("gene","up3"))

sd<-merge(enr3, merge(enr1, enr2, by="gene", all=T), by="gene", all=T)
sd<-sd[,c("gene", "up1", "up2", "up3")]
# write.table(sd, "source_data/figS4g.txt", quote=F, sep="\t", row.names = F)
  
ging<-read_xlsx("data/bioid_published_gingras.xlsx")
ging<-as.data.frame(ging)
colnames(ging)[4]<-"gene"

ging<-ging[order(ging$AvgP, decreasing = T),]

gingtop<-ging[1:300,]

nrow(merge(gingtop, int, by="gene"))


#figS4h
dev.off()
draw.pairwise.venn(nrow(int), nrow(gingtop), 
                   nrow(merge(gingtop, int, by="gene")),   euler.d=T,scaled=T, category=c("this study", "gingras"))

sd<-merge(gingtop, int, by="gene", all=T)
sd<-subset(sd, select=c("gene","AvgP","mean_enr","up" ))
# write.table(sd, "source_data/figS4h.txt", quote=F, sep="\t", row.names = F)
