git_dir<-"E:/work/hdlbp/git_rstudio/hdlbp/"

setwd(paste0(git_dir,"/data"))
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(dplyr)
library(ggrepel)
library(VennDiagram)
library(gplots)
library(data.table)







                    #### proteinGroups table load and preparation ####

PG <- fread("proteinGroups-RQ.txt", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T)
PG_noRQ <- fread("proteinGroups-notRQ.txt", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T)

# PG <- fread("second_run_181119/proteinGroups-RQ.txt", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T)
# PG_noRQ <- fread("second_run_181119/proteinGroups-notRQ.txt", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T)



PG$Gene.names <- sapply(strsplit(PG$Gene.names, ";"), "[", 1)
PG$Majority.protein.IDs <- sapply(strsplit(PG$Majority.protein.IDs, ";"), "[", 1)

PG_noRQ$Gene.names <- sapply(strsplit(PG_noRQ$Gene.names, ";"), "[", 1)
PG_noRQ$Majority.protein.IDs <- sapply(strsplit(PG_noRQ$Majority.protein.IDs, ";"), "[", 1)



#### Make sure both tables have same length
PG <- subset(PG, Majority.protein.IDs %in% PG_noRQ$Majority.protein.IDs)
PG_noRQ <- subset(PG_noRQ, Majority.protein.IDs %in% PG$Majority.protein.IDs)





intensities = c("iBAQ.L.Cyto_4h_Forward_1",
                "iBAQ.M.Cyto_4h_Forward_1",
                "iBAQ.H.Cyto_4h_Forward_1",
                
                "iBAQ.L.Cyto_4h_Forward_2",
                "iBAQ.M.Cyto_4h_Forward_2",
                "iBAQ.H.Cyto_4h_Forward_2",
                
                "iBAQ.L.Cyto_4h_Reverse_1",
                "iBAQ.M.Cyto_4h_Reverse_1",
                "iBAQ.H.Cyto_4h_Reverse_1",
                
                "iBAQ.L.Cyto_4h_Reverse_2",
                "iBAQ.M.Cyto_4h_Reverse_2",
                "iBAQ.H.Cyto_4h_Reverse_2",
                
                "iBAQ.L.Memb_4h_Forward_1",
                "iBAQ.M.Memb_4h_Forward_1",
                "iBAQ.H.Memb_4h_Forward_1",
                
                "iBAQ.L.Memb_4h_Forward_2",
                "iBAQ.M.Memb_4h_Forward_2",
                "iBAQ.H.Memb_4h_Forward_2",
                
                "iBAQ.L.Memb_4h_Reverse_1",
                "iBAQ.M.Memb_4h_Reverse_1",
                "iBAQ.H.Memb_4h_Reverse_1",
                
                "iBAQ.L.Memb_4h_Reverse_2",
                "iBAQ.M.Memb_4h_Reverse_2",
                "iBAQ.H.Memb_4h_Reverse_2")


ratios=c("Ratio.H.M.normalized.Cyto_4h_Forward_1",
         "Ratio.H.M.normalized.Cyto_4h_Forward_2",
         "Ratio.H.M.normalized.Cyto_4h_Reverse_1",
         "Ratio.H.M.normalized.Cyto_4h_Reverse_2",
         "Ratio.H.M.normalized.Memb_4h_Forward_1",
         "Ratio.H.M.normalized.Memb_4h_Forward_2",
         "Ratio.H.M.normalized.Memb_4h_Reverse_1",
         "Ratio.H.M.normalized.Memb_4h_Reverse_2")






#### Define Requantified ratios and intensities ####

### Ratios
requantified.ratios <- c()
unscrupulous.ratios <- c()
for (i in ratios) {
  PG[paste0("Requantified.", i)] <- ifelse(is.na(PG_noRQ[i]) &
                                             !is.na(PG[i]), T, F)
  
  requantified.ratios <- c(requantified.ratios, paste0("Requantified.", i))
  unscrupulous.ratios <- c(unscrupulous.ratios, paste0("Unscrupulous.", i))
  rm(i)
}



### intensities
requantified.intensities <- c()
for (i in intensities) {
  PG[paste0("Requantified.", i)] <- ifelse(PG_noRQ[i] == 0 &
                                             PG[i] != 0, T, F)
  
  requantified.intensities <- c(requantified.intensities, paste0("Requantified.", i))
  rm(i)
}









#### Unscrupulous Requantification ####



# Select Unscrupulous Requantified ratios (ratios between two channels when both are requantified)
for (i in 1:(length(ratios))) {
  #H.M ratios 
  PG[, paste0("Unscrupulous.", ratios[i])] <- ifelse((PG[, paste0("Requantified.", intensities[3+(i-1)*3])] &
                                                        PG[, paste0("Requantified.", intensities[2+(i-1)*3])]) &
                                                       !is.na(PG[, ratios[i]]), T, F)
  
  rm(i)
}





# Delet Unscrupulous ratios
for (i in 1:length(ratios)) {
  PG[,ratios[i]] <- ifelse(PG[,unscrupulous.ratios[i]], NA, PG[,ratios[i]])
  
  rm(i)
}





#### Filters ####
#PG <- PG_noRQ
#filter out contaminants, Rev and only identified by site
PG <- subset(PG, Reverse != "+")
PG <- subset(PG, Potential.contaminant != "+")
PG <- subset(PG, Only.identified.by.site != "+")




#transform Intensity and ratios to Log2 (or 10 if you prefer)
PG[c(intensities, ratios)] = log2(PG[c(intensities, ratios)])
# change Inf and NaN values for NA
is.na(PG[c(intensities, ratios)]) <- sapply(PG[c(intensities, ratios)], is.infinite)

is.na(PG[c(intensities, ratios)]) <- sapply(PG[c(intensities, ratios)], is.nan)



# how many proteins were identified (have intensities values not NA) in each L-H group pair
for (i in 1:length(intensities)) {
  cat(intensities[i])
  cat("\t")
  cat("\t")
  cat(nrow(PG[!is.na(PG[intensities[i]]),]))
  cat("\t")
  cat(mean(PG[,intensities[i]], na.rm = T))
  cat("\t")
#  cat(sum(PG[,requantified.intensities[i]], na.rm = T))
  cat("\n")
  rm(i)
}


# how many proteins were identified (have ratios values not NA) in each L-H group pair
for (i in 1:length(ratios)) {
  #group
  cat(ratios[i])
  cat("\t")
  cat("\t")
  #quantified ratios
  cat(nrow(PG[!is.na(PG[ratios[i]]),]))
  cat("\t")
  #ratio mean
  cat(mean(PG[,ratios[i]], na.rm = T))
  cat("\t")
  #Requantified ratios
#  cat(sum(!is.na(PG[ratios[i]]) & PG[, requantified.ratios[i]], na.rm = T))
  cat("\t")
  #Percentage of requantified ratios
#  cat(100*(sum(!is.na(PG[ratios[i]]) & PG[, requantified.ratios[i]], na.rm = T)) / nrow(PG[!is.na(PG[ratios[i]]),]))
  cat("\n")
  rm(i)
}




rm(PG_noRQ, requantified.intensities, requantified.ratios, unscrupulous.ratios)




#### Define PG_summ as working table ####

PG_summ <- subset(PG, select = c("Majority.protein.IDs",
                                 "Gene.names",
                                 intensities, ratios))

PG_summ$Gene.names <- sapply(strsplit(PG_summ$Gene.names, ";"), "[", 1)
PG_summ$Majority.protein.IDs <- sapply(strsplit(PG_summ$Majority.protein.IDs, ";"), "[", 1)

PG_summ<-subset(PG_summ, Gene.names!="") #PG_summ from carlos script
PG_summ<-subset(PG_summ, !duplicated(Gene.names))

setwd(paste0(git_dir,"/data"))
man<-read.delim("hdlbp_master_table_with_classes_uniq_tsig.txt", header=T)
pas<-merge(man, PG_summ, by.x="Symbol", by.y="Gene.names", all.x=T)

mel<-melt(pas, measure.vars = colnames(pas)[grepl("iBAQ.L", colnames(pas))], id.vars = c("Symbol", "localization_cat","tsig"))
plot<-subset(mel, !is.na(localization_cat))
ggplot(plot, aes(variable, value, fill=localization_cat))+geom_boxplot()+coord_flip()
plot<-subset(mel, !is.na(tsig))
ggplot(plot, aes(variable, value, fill=tsig))+geom_boxplot()+coord_flip()


pas$mean_l_mem<-rowMeans(pas[,colnames(pas)[grepl("iBAQ.L.Memb", colnames(pas))]])

mel<-melt(pas, measure.vars = colnames(pas)[grepl("Ratio.H.M.normalized.Cyto", colnames(pas))], id.vars = c("Symbol", "localization_cat","mean_l_mem", "tsig"))
plot<-subset(mel, !is.na(localization_cat))
ggplot(plot, aes(variable, value, fill=localization_cat))+geom_boxplot()+coord_flip()

pas$Ratio.H.M.normalized.Memb_4h_Reverse_1<--pas$Ratio.H.M.normalized.Memb_4h_Reverse_1
pas$Ratio.H.M.normalized.Memb_4h_Reverse_2<--pas$Ratio.H.M.normalized.Memb_4h_Reverse_2
# pas$Ratio.H.M.normalized.Memb_2h_Reverse_1<--pas$Ratio.H.M.normalized.Memb_2h_Reverse_1
# pas$Ratio.H.M.normalized.Memb_2h_Reverse_2<--pas$Ratio.H.M.normalized.Memb_2h_Reverse_2

pas$mean_memb_4h<-rowMeans(pas[,colnames(pas)[grepl("Ratio.H.M.normalized.Memb_4h", colnames(pas))]])

# pas$mean_memb_2h<-rowMeans(pas[,colnames(pas)[grepl("Ratio.H.M.normalized.Memb_2h", colnames(pas))]])

ggplot(subset(pas, !is.na(localization_cat)), aes(mean_memb_4h, colour=localization_cat))+stat_ecdf()+coord_cartesian(xlim=c(-.7,.7))+xlab("H/M KO/WT")
# ggplot(subset(pas, !is.na(localization_cat)), aes(mean_memb_2h, colour=localization_cat))+stat_ecdf()

ks.test(subset(pas, localization_cat=="membrane" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, localization_cat=="cytosolic" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
length(subset(pas,localization_cat=="membrane" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
length(subset(pas,localization_cat=="cytosolic" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])

ggplot(subset(pas, !is.na(tc_transcript_norm_cat)), aes(mean_memb_4h, colour=tc_transcript_norm_cat))+stat_ecdf()+coord_cartesian(xlim=c(-.7,.7))+xlab("H/M KO/WT")
ggplot(subset(pas, !is.na(tc_CDS_norm_cat)), aes(mean_memb_4h, colour=tc_CDS_norm_cat))+stat_ecdf()+coord_cartesian(xlim=c(-.7,.7))+xlab("H/M KO/WT")

ggplot(subset(pas, !is.na(loc_tar_CDS) & localization_cat=="membrane"), aes(mean_memb_4h, colour=loc_tar_CDS))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+xlab("H/M KO/WT")
ggplot(subset(pas, !is.na(loc_tar_CDS) & localization_cat=="cytosolic"), aes(mean_memb_4h, colour=loc_tar_CDS))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))+xlab("H/M KO/WT")

ks.test(subset(pas, tc_transcript_norm_cat=="tc>1.12 & tc<132.32" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
            subset(pas, tc_transcript_norm_cat=="nontarget" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
ks.test(subset(pas, tc_transcript_norm_cat=="tc>1.12 & tc<132.32" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, tc_transcript_norm_cat=="tc>0.3 & tc<1.12" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
ks.test(subset(pas, tc_transcript_norm_cat=="tc>1.12 & tc<132.32" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, tc_transcript_norm_cat=="tc<0.3" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
ks.test(subset(pas, tc_transcript_norm_cat=="nontarget" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, tc_transcript_norm_cat=="tc>0.3 & tc<1.12" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
ks.test(subset(pas, tc_transcript_norm_cat=="nontarget" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, tc_transcript_norm_cat=="tc<0.3" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
ks.test(subset(pas, tc_transcript_norm_cat=="tc>0.3 & tc<1.12" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, tc_transcript_norm_cat=="tc<0.3" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])

ks.test(subset(pas, tc_CDS_norm_cat=="tc>1.04 & tc<65.26" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, tc_CDS_norm_cat=="nontarget" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
ks.test(subset(pas, tc_CDS_norm_cat=="tc>1.04 & tc<65.26" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, tc_CDS_norm_cat=="tc>0.27 & tc<1.04" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
ks.test(subset(pas, tc_CDS_norm_cat=="tc>1.04 & tc<65.26" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, tc_CDS_norm_cat=="tc<0.27" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
ks.test(subset(pas, tc_CDS_norm_cat=="nontarget" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, tc_CDS_norm_cat=="tc>0.27 & tc<1.04" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
ks.test(subset(pas, tc_CDS_norm_cat=="nontarget" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, tc_CDS_norm_cat=="tc<0.27" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
ks.test(subset(pas, tc_CDS_norm_cat=="tc>0.27 & tc<1.04" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, tc_CDS_norm_cat=="tc<0.27" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])



ks.test(subset(pas, tc_transcript_norm_cat=="tc>1.12 & tc<132.32" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1), select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1],
            subset(pas, tc_transcript_norm_cat=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1), select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1])
ks.test(subset(pas, tc_transcript_norm_cat=="tc>1.12 & tc<132.32" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_2), select="Ratio.H.M.normalized.Memb_4h_Forward_2")[,1],
            subset(pas, tc_transcript_norm_cat=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_2), select="Ratio.H.M.normalized.Memb_4h_Forward_2")[,1])
ks.test(subset(pas, tc_transcript_norm_cat=="tc>1.12 & tc<132.32" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_2), select="Ratio.H.M.normalized.Memb_4h_Reverse_2")[,1],
            subset(pas, tc_transcript_norm_cat=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_2), select="Ratio.H.M.normalized.Memb_4h_Reverse_2")[,1])
ks.test(subset(pas, tc_transcript_norm_cat=="tc>1.12 & tc<132.32" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_1), select="Ratio.H.M.normalized.Memb_4h_Reverse_1")[,1],
            subset(pas, tc_transcript_norm_cat=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_1), select="Ratio.H.M.normalized.Memb_4h_Reverse_1")[,1])

length(subset(pas, tc_transcript_norm_cat=="tc>1.12 & tc<132.32" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
length(subset(pas, tc_transcript_norm_cat=="tc>0.3 & tc<1.12" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
length(subset(pas, tc_transcript_norm_cat=="tc<0.3" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
length(subset(pas, tc_transcript_norm_cat=="nontarget" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])

length(subset(pas, tc_CDS_norm_cat=="tc>1.04 & tc<65.26" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
length(subset(pas, tc_CDS_norm_cat=="tc>0.27 & tc<1.04" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
length(subset(pas, tc_CDS_norm_cat=="tc<0.27" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
length(subset(pas, tc_CDS_norm_cat=="nontarget" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])

ggplot(subset(pas, !is.na(tc_transcript_norm_cat)), aes(Ratio.H.M.normalized.Memb_4h_Forward_1, colour=tc_transcript_norm_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1.5,1.5))
ggplot(subset(pas, !is.na(tc_transcript_norm_cat)), aes(Ratio.H.M.normalized.Memb_4h_Forward_2, colour=tc_transcript_norm_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1.5,1.5))
ggplot(subset(pas, !is.na(tc_transcript_norm_cat)), aes(-Ratio.H.M.normalized.Memb_4h_Reverse_1, colour=tc_transcript_norm_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1.5,1.5))
ggplot(subset(pas, !is.na(tc_transcript_norm_cat)), aes(-Ratio.H.M.normalized.Memb_4h_Reverse_2, colour=tc_transcript_norm_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1.5,1.5))

pas$target<-ifelse(is.na(pas$tc_transcript_norm_cat), NA,
                   ifelse(pas$tc_transcript_norm_cat=="nontarget", "nontarget", "target"))
pas$target_loc<-ifelse(is.na(pas$target), NA,
                   ifelse(pas$target=="nontarget", "nontarget", 
                          ifelse(pas$target=="target" & pas$localization_cat=="cytosolic", "cyto_target",
                                 ifelse(pas$target=="target" & pas$localization_cat=="membrane", "mem_target", NA))))

ggplot(subset(pas, gene_biotype=="protein_coding" & !is.na(target)), aes(mean_memb_4h, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))
ggplot(subset(pas, gene_biotype=="protein_coding" & !is.na(target)), aes(Ratio.H.M.normalized.Memb_4h_Forward_1, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))
ggplot(subset(pas, gene_biotype=="protein_coding" & !is.na(target)), aes(Ratio.H.M.normalized.Memb_4h_Forward_2, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))
ggplot(subset(pas, gene_biotype=="protein_coding" & !is.na(target)), aes(-Ratio.H.M.normalized.Memb_4h_Reverse_1, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))
ggplot(subset(pas, gene_biotype=="protein_coding" & !is.na(target)), aes(-Ratio.H.M.normalized.Memb_4h_Reverse_2, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))

ggplot(subset(pas, gene_biotype=="protein_coding" & !is.na(target_loc)), aes(mean_memb_4h, colour=target_loc))+stat_ecdf()+coord_cartesian(xlim=c(-.7,.7))+xlab("mean H/M (KO/WT)")
ggplot(subset(pas, gene_biotype=="protein_coding" & !is.na(target_loc)), aes(Ratio.H.M.normalized.Memb_4h_Forward_1, colour=target_loc))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))
ggplot(subset(pas, gene_biotype=="protein_coding" & !is.na(target_loc)), aes(Ratio.H.M.normalized.Memb_4h_Forward_2, colour=target_loc))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))
ggplot(subset(pas, gene_biotype=="protein_coding" & !is.na(target_loc)), aes(Ratio.H.M.normalized.Memb_4h_Reverse_1, colour=target_loc))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))
ggplot(subset(pas, gene_biotype=="protein_coding" & !is.na(target_loc)), aes(Ratio.H.M.normalized.Memb_4h_Reverse_2, colour=target_loc))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))


length(subset(pas, gene_biotype=="protein_coding" & !is.na(mean_memb_4h) & target=="target", select="mean_memb_4h")[,1])
length(subset(pas, gene_biotype=="protein_coding" & !is.na(mean_memb_4h) & target=="nontarget", select="mean_memb_4h")[,1])

length(subset(pas, gene_biotype=="protein_coding" & !is.na(mean_memb_4h) & target_loc=="mem_target", select="mean_memb_4h")[,1])
length(subset(pas, gene_biotype=="protein_coding" & !is.na(mean_memb_4h) & target_loc=="cyto_target", select="mean_memb_4h")[,1])
length(subset(pas, gene_biotype=="protein_coding" & !is.na(mean_memb_4h) & target_loc=="nontarget", select="mean_memb_4h")[,1])


ks.test(subset(pas, target=="target" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1), select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1],
        subset(pas, target=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1), select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1])
ks.test(subset(pas, target=="target" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_2), select="Ratio.H.M.normalized.Memb_4h_Forward_2")[,1],
        subset(pas, target=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_2), select="Ratio.H.M.normalized.Memb_4h_Forward_2")[,1])
ks.test(subset(pas, target=="target" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_1), select="Ratio.H.M.normalized.Memb_4h_Reverse_1")[,1],
        subset(pas, target=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_1), select="Ratio.H.M.normalized.Memb_4h_Reverse_1")[,1])
ks.test(subset(pas, target=="target" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_2), select="Ratio.H.M.normalized.Memb_4h_Reverse_2")[,1],
        subset(pas, target=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_2), select="Ratio.H.M.normalized.Memb_4h_Reverse_2")[,1])

ks.test(subset(pas, target_loc=="mem_target" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, target_loc=="cyto_target" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
ks.test(subset(pas, target_loc=="mem_target" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, target_loc=="nontarget" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
ks.test(subset(pas, target_loc=="cyto_target" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
        subset(pas, target_loc=="nontarget" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])

ks.test(subset(pas, target_loc=="mem_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1), select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1],
        subset(pas, target_loc=="cyto_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1), select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1])
ks.test(subset(pas, target_loc=="mem_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_2), select="Ratio.H.M.normalized.Memb_4h_Forward_2")[,1],
        subset(pas, target_loc=="cyto_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_2), select="Ratio.H.M.normalized.Memb_4h_Forward_2")[,1])
ks.test(subset(pas, target_loc=="mem_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_1), select="Ratio.H.M.normalized.Memb_4h_Reverse_1")[,1],
        subset(pas, target_loc=="cyto_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_1), select="Ratio.H.M.normalized.Memb_4h_Reverse_1")[,1])
ks.test(subset(pas, target_loc=="mem_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_2), select="Ratio.H.M.normalized.Memb_4h_Reverse_2")[,1],
        subset(pas, target_loc=="cyto_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_2), select="Ratio.H.M.normalized.Memb_4h_Reverse_2")[,1])

ks.test(subset(pas, target_loc=="mem_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1), select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1],
        subset(pas, target_loc=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1), select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1])
ks.test(subset(pas, target_loc=="mem_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_2), select="Ratio.H.M.normalized.Memb_4h_Forward_2")[,1],
        subset(pas, target_loc=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_2), select="Ratio.H.M.normalized.Memb_4h_Forward_2")[,1])
ks.test(subset(pas, target_loc=="mem_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_1), select="Ratio.H.M.normalized.Memb_4h_Reverse_1")[,1],
        subset(pas, target_loc=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_1), select="Ratio.H.M.normalized.Memb_4h_Reverse_1")[,1])
ks.test(subset(pas, target_loc=="mem_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_2), select="Ratio.H.M.normalized.Memb_4h_Reverse_2")[,1],
        subset(pas, target_loc=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_2), select="Ratio.H.M.normalized.Memb_4h_Reverse_2")[,1])

pas$reproducible<-ifelse(!is.na(pas$Ratio.H.M.normalized.Memb_4h_Forward_1) & !is.na(pas$Ratio.H.M.normalized.Memb_4h_Forward_2) &
                           !is.na(pas$Ratio.H.M.normalized.Memb_4h_Reverse_1) & !is.na(pas$Ratio.H.M.normalized.Memb_4h_Reverse_2), "reproducible", "not_reproducible")


ks.test(subset(pas, reproducible=="reproducible" & target_loc=="mem_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1), select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1],
        subset(pas, reproducible=="reproducible" & target_loc=="cyto_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1), select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1])
ks.test(subset(pas, reproducible=="reproducible" & target_loc=="mem_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_2), select="Ratio.H.M.normalized.Memb_4h_Forward_2")[,1],
        subset(pas, reproducible=="reproducible" & target_loc=="cyto_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_2), select="Ratio.H.M.normalized.Memb_4h_Forward_2")[,1])
ks.test(subset(pas, reproducible=="reproducible" & target_loc=="mem_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_1), select="Ratio.H.M.normalized.Memb_4h_Reverse_1")[,1],
        subset(pas, reproducible=="reproducible" & target_loc=="cyto_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_1), select="Ratio.H.M.normalized.Memb_4h_Reverse_1")[,1])
ks.test(subset(pas, reproducible=="reproducible" & target_loc=="mem_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_2), select="Ratio.H.M.normalized.Memb_4h_Reverse_2")[,1],
        subset(pas, reproducible=="reproducible" & target_loc=="cyto_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_2), select="Ratio.H.M.normalized.Memb_4h_Reverse_2")[,1])

ks.test(subset(pas, reproducible=="reproducible" & target_loc=="mem_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1), select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1],
        subset(pas, reproducible=="reproducible" & target_loc=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1), select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1])
ks.test(subset(pas, reproducible=="reproducible" & target_loc=="mem_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_2), select="Ratio.H.M.normalized.Memb_4h_Forward_2")[,1],
        subset(pas, reproducible=="reproducible" & target_loc=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_2), select="Ratio.H.M.normalized.Memb_4h_Forward_2")[,1])
ks.test(subset(pas, reproducible=="reproducible" & target_loc=="mem_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_1), select="Ratio.H.M.normalized.Memb_4h_Reverse_1")[,1],
        subset(pas, reproducible=="reproducible" & target_loc=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_1), select="Ratio.H.M.normalized.Memb_4h_Reverse_1")[,1])
ks.test(subset(pas, reproducible=="reproducible" & target_loc=="mem_target" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_2), select="Ratio.H.M.normalized.Memb_4h_Reverse_2")[,1],
        subset(pas, reproducible=="reproducible" & target_loc=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_2), select="Ratio.H.M.normalized.Memb_4h_Reverse_2")[,1])


ggplot(subset(pas, reproducible=="reproducible" & !is.na(target)), aes(mean_memb_4h, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))
ggplot(subset(pas, reproducible=="reproducible" & !is.na(target)), aes(Ratio.H.M.normalized.Memb_4h_Forward_1, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))
ggplot(subset(pas, reproducible=="reproducible" & !is.na(target)), aes(Ratio.H.M.normalized.Memb_4h_Forward_2, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))
ggplot(subset(pas, reproducible=="reproducible" & !is.na(target)), aes(-Ratio.H.M.normalized.Memb_4h_Reverse_1, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))
ggplot(subset(pas, reproducible=="reproducible" & !is.na(target)), aes(-Ratio.H.M.normalized.Memb_4h_Reverse_2, colour=target))+stat_ecdf()+coord_cartesian(xlim=c(-1,1))

ks.test(subset(pas, reproducible=="reproducible" & target=="target" & !is.na(mean_memb_4h), select=mean_memb_4h)[,1],
        subset(pas, reproducible=="reproducible" & target=="nontarget" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
wilcox.test(subset(pas, reproducible=="reproducible" & target=="target" & !is.na(mean_memb_4h), select=mean_memb_4h)[,1],
        subset(pas, reproducible=="reproducible" & target=="nontarget" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
ks.test(subset(pas, reproducible=="reproducible" & target=="target" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1), select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1],
        subset(pas, reproducible=="reproducible" & target=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1), select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1])
ks.test(subset(pas, reproducible=="reproducible" & target=="target" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_2), select="Ratio.H.M.normalized.Memb_4h_Forward_2")[,1],
        subset(pas, reproducible=="reproducible" & target=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_2), select="Ratio.H.M.normalized.Memb_4h_Forward_2")[,1])
ks.test(subset(pas, reproducible=="reproducible" & target=="target" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_1), select="Ratio.H.M.normalized.Memb_4h_Reverse_1")[,1],
        subset(pas, reproducible=="reproducible" & target=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_1), select="Ratio.H.M.normalized.Memb_4h_Reverse_1")[,1])
ks.test(subset(pas, reproducible=="reproducible" & target=="target" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_2), select="Ratio.H.M.normalized.Memb_4h_Reverse_2")[,1],
        subset(pas, reproducible=="reproducible" & target=="nontarget" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_2), select="Ratio.H.M.normalized.Memb_4h_Reverse_2")[,1])

length(subset(pas, reproducible=="reproducible" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1) & target=="target", select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1])
length(subset(pas, !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1) & target=="target", select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1])

length(subset(pas, reproducible=="reproducible" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1) & target=="target", select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1])


ggplot(subset(pas, !is.na(localization_cat)), aes(Ratio.H.M.normalized.Memb_4h_Forward_1, colour=localization_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1.5,1.5))
ggplot(subset(pas, !is.na(localization_cat)), aes(Ratio.H.M.normalized.Memb_4h_Forward_2, colour=localization_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1.5,1.5))
ggplot(subset(pas, !is.na(localization_cat)), aes(-Ratio.H.M.normalized.Memb_4h_Reverse_1, colour=localization_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1.5,1.5))
ggplot(subset(pas, !is.na(localization_cat)), aes(-Ratio.H.M.normalized.Memb_4h_Reverse_2, colour=localization_cat))+stat_ecdf()+coord_cartesian(xlim=c(-1.5,1.5))



wilcox.test(subset(pas, localization_cat=="membrane" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
            subset(pas, localization_cat=="cytosolic" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])

ks.test(subset(pas, localization_cat=="membrane" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1],
            subset(pas, localization_cat=="cytosolic" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])

length(subset(pas, localization_cat=="membrane" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])
length(subset(pas, localization_cat=="cytosolic" & !is.na(mean_memb_4h), select="mean_memb_4h")[,1])


ks.test(subset(pas, localization_cat=="membrane" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_1), select="Ratio.H.M.normalized.Memb_4h_Reverse_1")[,1],
            subset(pas, localization_cat=="cytosolic" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_1), select="Ratio.H.M.normalized.Memb_4h_Reverse_1")[,1])
ks.test(subset(pas, localization_cat=="membrane" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_2), select="Ratio.H.M.normalized.Memb_4h_Reverse_2")[,1],
            subset(pas, localization_cat=="cytosolic" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_2), select="Ratio.H.M.normalized.Memb_4h_Reverse_2")[,1])
ks.test(subset(pas, localization_cat=="membrane" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1), select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1],
            subset(pas, localization_cat=="cytosolic" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1), select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1])
ks.test(subset(pas, localization_cat=="membrane" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_2), select="Ratio.H.M.normalized.Memb_4h_Forward_2")[,1],
            subset(pas, localization_cat=="cytosolic" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_2), select="Ratio.H.M.normalized.Memb_4h_Forward_2")[,1])

length(subset(pas, localization_cat=="membrane" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_1), select="Ratio.H.M.normalized.Memb_4h_Reverse_1")[,1])
length(subset(pas, localization_cat=="membrane" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_2), select="Ratio.H.M.normalized.Memb_4h_Reverse_2")[,1])
length(subset(pas, localization_cat=="membrane" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1), select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1])
length(subset(pas, localization_cat=="membrane" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_2), select="Ratio.H.M.normalized.Memb_4h_Forward_2")[,1])

length(subset(pas, localization_cat=="cytosolic" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_1), select="Ratio.H.M.normalized.Memb_4h_Reverse_1")[,1])
length(subset(pas, localization_cat=="cytosolic" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_2), select="Ratio.H.M.normalized.Memb_4h_Reverse_2")[,1])
length(subset(pas, localization_cat=="cytosolic" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1), select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1])
length(subset(pas, localization_cat=="cytosolic" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_2), select="Ratio.H.M.normalized.Memb_4h_Forward_2")[,1])

length(subset(pas, reproducible=="reproducible" & localization_cat=="membrane" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_1), select="Ratio.H.M.normalized.Memb_4h_Reverse_1")[,1])
length(subset(pas, reproducible=="reproducible" & localization_cat=="membrane" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_2), select="Ratio.H.M.normalized.Memb_4h_Reverse_2")[,1])
length(subset(pas, reproducible=="reproducible" & localization_cat=="membrane" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1), select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1])
length(subset(pas, reproducible=="reproducible" & localization_cat=="membrane" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_2), select="Ratio.H.M.normalized.Memb_4h_Forward_2")[,1])

length(subset(pas, reproducible=="reproducible" & localization_cat=="cytosolic" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_1), select="Ratio.H.M.normalized.Memb_4h_Reverse_1")[,1])
length(subset(pas, reproducible=="reproducible" & localization_cat=="cytosolic" & !is.na(Ratio.H.M.normalized.Memb_4h_Reverse_2), select="Ratio.H.M.normalized.Memb_4h_Reverse_2")[,1])
length(subset(pas, reproducible=="reproducible" & localization_cat=="cytosolic" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_1), select="Ratio.H.M.normalized.Memb_4h_Forward_1")[,1])
length(subset(pas, reproducible=="reproducible" & localization_cat=="cytosolic" & !is.na(Ratio.H.M.normalized.Memb_4h_Forward_2), select="Ratio.H.M.normalized.Memb_4h_Forward_2")[,1])


ggplot(pas, aes(tsig, log2(utr3_vs_cds.norm.mean)))+geom_boxplot()

       
nrow(subset(pas, tpm_cutoff>=10 & localization_cat=="cytosolic" & target=="target" & gene_biotype=="protein_coding"& !is.na(log2FoldChange.mem.cyt.293)))
nrow(subset(pas, tpm_cutoff>=10 & localization_cat=="membrane" & target=="target" & gene_biotype=="protein_coding" & !is.na(log2FoldChange.mem.cyt.293)))

nrow(subset(pas, tpm_cutoff>=10 & localization_cat=="membrane"  & gene_biotype=="protein_coding" & !is.na(log2FoldChange.mem.cyt.293)))
nrow(subset(pas, tpm_cutoff>=10 & localization_cat=="cytosolic"  & gene_biotype=="protein_coding" & !is.na(log2FoldChange.mem.cyt.293)))


int<-subset(pas, mean_memb_4h<0 & target=="target" & log2FoldChange.ribo.rna.KO.WT<0 & localization_cat=="membrane")
# int<-subset(pas, mean_memb_4h<0 & target=="target" & log2FoldChange.ribo.rna.KO.WT<0 )

# int<-subset(pas, mean_memb_4h<0 & log2FoldChange.ribo.rna.KO.WT<0 & localization_cat=="membrane")
# int<-subset(pas, mean_memb_4h<0 & log2FoldChange.ribo.rna.KO.WT<0 )
int<-subset(int, !duplicated(Symbol))
int$Symbol<-factor(int$Symbol,  levels=int$Symbol[order(int$mean_memb_4h, decreasing = T)])
int<-int[order(int$mean_memb_4h, decreasing = F),]
int<-int[1:60,]

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

ggplot(int, aes(Symbol, mean_memb_4h ))+geom_bar(stat = "identity")+coord_flip()
ggplot(int, aes(Symbol, log2FoldChange.ribo.rna.KO.WT ))+geom_bar(stat = "identity")+coord_flip()

ggplot(subset(pas, tpm_cutoff>=10), aes(mean_memb_4h, log2FoldChange.ribo.rna.KO.WT))+geom_point(shape=1, size=0.5)+coord_cartesian(ylim=c(-1,1))+facet_wrap(~localization_cat)
ggplot(subset(pas, tpm_cutoff>=10), aes(mean_memb_4h, log2FoldChange.ribo.rna.KO.WT))+geom_point(shape=1, size=0.5)+coord_cartesian(ylim=c(-1,1))+facet_wrap(~tc_transcript_norm_cat)
ggplot(subset(pas, tpm_cutoff>=10), aes(mean_memb_4h, log2FoldChange.ribo.rna.KO.WT))+geom_text(aes(label=Symbol), size=2)+facet_wrap(~tc_transcript_norm_cat)


int<-subset(int, select=c("Symbol", "mean_memb_4h", "log2FoldChange.ribo.rna.KO.WT", "tc_transcript_norm"))
colnames(int)<-c("Symbol", "pSILAC_ko_wt", "riboseq_ko_wt", "hdlbp_clip")
# write.table(int, "pSilacList.txt", quote=F, sep="\t", row.names = F)

##supp table
supp<-pas[,c(1:2,8,45,77,90:101,106:109,111,112)]
colnames(supp)<-c("Symbol",
                  "gene_id",
                  "localization",
                  "tarsig",
                  "Majority.protein.IDs",
                  "iBAQ.L.Memb_4h_Forward_1",
                  "iBAQ.M.Memb_4h_Forward_1",
                  "iBAQ.H.Memb_4h_Forward_1",
                  "iBAQ.L.Memb_4h_Forward_2",
                  "iBAQ.M.Memb_4h_Forward_2",
                  "iBAQ.H.Memb_4h_Forward_2",
                  "iBAQ.L.Memb_4h_Reverse_1",
                  "iBAQ.M.Memb_4h_Reverse_1",
                  "iBAQ.H.Memb_4h_Reverse_1",
                  "iBAQ.L.Memb_4h_Reverse_2",
                  "iBAQ.M.Memb_4h_Reverse_2",
                  "iBAQ.H.Memb_4h_Reverse_2",
                  "Ratio.H.M.normalized.Memb_4h_Forward_1",
                  "Ratio.H.M.normalized.Memb_4h_Forward_2",
                  "Ratio.H.M.normalized.Memb_4h_Reverse_1",
                  "Ratio.H.M.normalized.Memb_4h_Reverse_2",
                  "mean.Ratio.H.M.normalized.Memb_4h",
                  "hdlbp.target.yes.no")

setwd(paste0(git_dir, "/supp_tables"))
write.table(supp, "tableS7_psilac.txt", quote=F, sep="\t", row.names = F)
