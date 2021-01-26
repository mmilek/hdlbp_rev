git_dir<-"E:/work/hdlbp/git_rstudio/hdlbp/"
setwd(paste0(git_dir, "/data"))
dat<-read.delim("proteinGroups - Exp150+144.txt", header=T)
library(ggplot2)

use<-dat[,c(2,7,11,12,53:57,65,66,42:46)]
colnames(use)[5:9]<-c("lfq.hdlbp1", "lfq.ctrl2", "lfq.hdlbp2", "lfq.ctrl3", "lfq.hdlbp3")

colnames(use)[12:16]<-c("int.hdlbp1", "int.ctrl2", "int.hdlbp2", "int.ctrl3", "int.hdlbp3")

use$Gene.names<-as.character(use$Gene.names)
use$enr1<-ifelse(use$lfq.ctrl2==0, NA, use$lfq.hdlbp1/use$lfq.ctrl2)
use$enr2<-ifelse(use$lfq.ctrl2==0, NA, use$lfq.hdlbp2/use$lfq.ctrl2)
use$enr3<-ifelse(use$lfq.ctrl3==0, NA, use$lfq.hdlbp3/use$lfq.ctrl3)

library(corrplot)
corrplot(cor(use[,17:19], use="pairwise.complete.obs"), method="color", addCoef.col = "white")
corrplot(cor(use[,5:9], use="pairwise.complete.obs"), method="color", addCoef.col = "white")

pairs(~enr1+enr2+enr3, log="xy", cex=0.5, lower.panel=NULL, data=use)

use$gene<-gsub(";.*","",use$Gene.names)

use$positive<-ifelse(use$gene=="HDLBP", "HDLBP", NA)


ggplot(subset(use, Razor...unique.peptides>=5), aes(log2(lfq.hdlbp1), log2(lfq.ctrl2), colour=positive))+geom_point()+geom_text(aes(label=gene))
ggplot(subset(use, Razor...unique.peptides>=5), aes(log2(lfq.hdlbp2), log2(lfq.ctrl2), colour=positive))+geom_point()+geom_text(aes(label=gene))
ggplot(subset(use, Razor...unique.peptides>=5), aes(log2(lfq.hdlbp3), log2(lfq.ctrl3), colour=positive))+geom_point()+geom_text(aes(label=gene))


ggplot(subset(use, Razor...unique.peptides>=5), aes(log2(enr2), log2(enr3), colour=positive))+geom_point()+geom_text(aes(label=positive))
ggplot(subset(use, Razor...unique.peptides>=5), aes(log2(enr1), log2(enr2), colour=positive))+geom_point()+geom_text(aes(label=positive))
ggplot(subset(use, Razor...unique.peptides>=5), aes(log2(enr1), log2(enr3), colour=positive))+geom_point()+geom_text(aes(label=positive))

ggplot(subset(use, Razor...unique.peptides>=5), aes(log2(enr1), log2(enr2), colour=positive))+geom_point()+geom_text(aes(label=gene))
ggplot(subset(use, Razor...unique.peptides>=5), aes(log2(enr1), log2(enr3), colour=positive))+geom_point()+geom_text(aes(label=gene))
ggplot(subset(use, Razor...unique.peptides>=5), aes(log2(enr2), log2(enr3), colour=positive))+geom_point()+geom_text(aes(label=gene))

thr<-3
use$up1<-ifelse( (use$enr1>=thr& log2(use$lfq.hdlbp1)>=27) | (use$lfq.ctrl2==0 & log2(use$lfq.hdlbp1)>=27),"enr1", NA)
use$up2<-ifelse( (use$enr2>=thr & log2(use$lfq.hdlbp2)>=27) | (use$lfq.ctrl2==0 & log2(use$lfq.hdlbp2)>=27),"enr2", NA)
use$up3<-ifelse( (use$enr3>=thr & log2(use$lfq.hdlbp3)>=27) | (use$lfq.ctrl3==0 & log2(use$lfq.hdlbp3)>=27),"enr3", NA)

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
ggplot(top,aes(factor(gene, levels=top$gene[order(top$mean_lfq, decreasing = F)]), 1))+geom_point(aes(colour=log10(mean_lfq)))
ggplot(top)+geom_point(aes(factor(gene, levels=top$gene[order(top$mean_lfq, decreasing = F)]), 1, colour=log10(mean_lfq), size=log2(imp_enr)))+coord_flip()+xlab("")+scale_colour_gradient(low="dodgerblue2", high="orange2",limits=c(6,13))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

int<-int[order(int$mean_enr, decreasing = T),]
top<-int[1:60,]
ggplot(top,aes(factor(gene, levels=top$gene[order(top$mean_enr, decreasing = F)]), log2(mean_enr)))+geom_bar(stat = "identity")+coord_flip()
ggplot(top,aes(factor(gene, levels=top$gene[order(top$mean_lfq, decreasing = F)]), log2(mean_enr)))+geom_bar(stat = "identity")+coord_flip()
ggplot(top,aes(factor(gene, levels=top$gene[order(top$mean_lfq, decreasing = F)]), log2(imp_enr)))+geom_bar(stat = "identity")+coord_flip()
ggplot(top,aes(factor(gene, levels=top$gene[order(top$mean_lfq, decreasing = F)]), log2(imp_enr), fill=imp_enr_cat))+geom_bar(stat = "identity")+coord_flip()+ylim(0,10)

int<-int[order(int$imp_enr, decreasing = T),]
top<-int[1:60,]
ggplot(top,aes(factor(gene, levels=top$gene[order(top$mean_lfq, decreasing = F)]), log2(imp_enr)))+geom_bar(stat = "identity")+coord_flip()

top<-subset(int, is.na(mean_enr))
top<-top[order(top$mean_lfq, decreasing = T),]
top<-top[1:60,]
ggplot(top,aes(factor(gene, levels=top$gene[order(top$mean_lfq, decreasing = F)]), log2(mean_lfq)))+geom_bar(stat = "identity")+coord_flip()

int<-int[order(int$Razor...unique.peptides, decreasing = T),]
top<-int[1:60,]
ggplot(top,aes(factor(gene, levels=top$gene[order(top$Razor...unique.peptides, decreasing = F)]), log2(Razor...unique.peptides)))+geom_bar(stat = "identity")+coord_flip()

ggplot(use,aes( log2(imp_enr) , log10(mean_lfq) ,colour=up))+geom_point()


## supp table
supp<-subset(use, Reverse!="+" & Potential.contaminant!="+")
supp<-supp[,c(1:3,5:9,17:19,25,27,28)]
colnames(supp)<-c("Majority.protein.IDs","Gene.names",
                  "Razor.unique.peptides","lfq.hdlbp1.dox1",
                  "lfq.ctrl2.nodox2","lfq.hdlbp2.dox2",
                  "lfq.ctrl3.nodox3","lfq.hdlbp3.dox3",
                  "enrichment1","enrichment2",
                  "enrichment3","enriched.yes.no",
                  "mean_enrichment","mean_lfq")
#write.table(supp,"../supp_tables/tableS6_bioid.txt", quote=F,sep="\t",row.names=F)

##overlap between replicates

enr1<-subset(use, up1=="enr1")
enr2<-subset(use, up2=="enr2")
enr3<-subset(use, up3=="enr3")

library(VennDiagram)

draw.triple.venn(nrow(enr1), nrow(enr2), nrow(enr3),
                 nrow(merge(enr1, enr2, by="gene")),  nrow(merge(enr2, enr3, by="gene")), nrow(merge(enr1, enr3, by="gene")),
                           nrow(merge(enr3, merge(enr1, enr2, by="gene"), by="gene")), euler.d=T,scaled=T, category=c("replicate 1", "replicate 2", "replicate 3"))

library(readxl)
ging<-read_xlsx("bioid_published_gingras.xlsx")
ging<-as.data.frame(ging)
colnames(ging)[4]<-"gene"
nrow(merge(ging, int, by="gene"))

ging<-ging[order(ging$SpecSum, decreasing = T),]

gingtop<-ging[1:300,]

nrow(merge(gingtop, int, by="gene"))
ging<-ging[order(ging$AvgP, decreasing = T),]

gingtop<-ging[1:300,]

nrow(merge(gingtop, int, by="gene"))

draw.pairwise.venn(nrow(int), nrow(gingtop), 
                 nrow(merge(gingtop, int, by="gene")),   euler.d=T,scaled=T, category=c("this study", "gingras"))

subset(merge(gingtop, int, by="gene"), select="gene")


##prepare for crapome

colnames(dat)
crap<-dat[,c(7,58:62)]

colnames(crap)<-c("PROTID", "HDLBP1_SPC", "CONTROL2_SPC", "HDLBP2_SPC", "CONTROL3_SPC", "HDLBP3_SPC")
crap<-subset(crap, PROTID!="")
crap<-subset(crap, !duplicated(PROTID))
crap$PROTID<-gsub(";.*","",crap$PROTID)
crap<-subset(crap, !duplicated(PROTID))
write.table(crap,"hdlbp_crapome.txt", quote=F, sep="\t", row.names=F, col.names = T)

pro<-read.delim("hdlbp_prohits_input.txt", header=T)
# pro<-subset(pro, avgSPC>=6 & SP>0.95&FC_A>=2)
# write.table(pro, "hdlbp_prohits_input_filtered.txt", quote=F, sep="\t", row.names=F, col.names = T)

ggplot(pro, aes(log2(FC_A), log2(avgSPC)))+geom_point()
ggplot(subset(pro, SP>=0.95), aes(log2(FC_A), log2(avgSPC)))+geom_point()

pro$reg<-ifelse(pro$SP>=0.95 & pro$avgSPC>=6 & pro$FC_A>=2, "reg", NA)
ggplot(pro, aes( log2(avgSPC),log2(FC_A),colour=reg))+geom_point()

pro<-subset(pro, avgSPC>=6 & SP>0.95&FC_A>=2)

nrow(merge(gingtop, pro, by.x="gene", by.y="prey"))
subset(merge(int, pro, by.x="gene", by.y="prey"), select="gene")


pro<-read.delim("hdlbp_prohits_input_with_neg_controls.txt", header=T)
pro<-subset(pro, avgSPC>=6 & SP>0.95&FC_A>=2)
nrow(merge(int, pro, by.x="gene", by.y="prey"))

# write.table(pro, "hdlbp_prohits_input_filtered.txt", quote=F, sep="\t", row.names=F, col.names = T)
nrow(merge(gingtop, pro, by.x="gene", by.y="prey"))

##check saint
setwd("~/Google Drive/hdlbp/bioid/")

sai<-read.delim("hdlbp_mscount_all_controls.saint.output", header=T)

ggplot(sai, aes(log2(AvgSpec), log2(FoldChange), colour=AvgP ))+geom_point()
ggplot(sai, aes(log2(AvgSpec), log2(FoldChange), colour=BFDR ))+geom_point()

cand<-subset(sai, FoldChange>=3 & AvgP>=0.95 & BFDR<=0.1 & AvgSpec>=20)
cand<-cand[order(cand$AvgSpec, decreasing = T),]
cand<-cand[order(cand$FoldChange, decreasing = T),]

nrow(merge(int, cand, by.x="gene", by.y="Prey"))

