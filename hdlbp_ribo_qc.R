library(riboWaltz)
library(data.table)
# library(GenomicRanges)
library(ggplot2)
# data(mm81cdna)
# data(reads_list)
# setwd("/Volumes/landthaler/pcp/projects/miha/koshi_riboseq/ribo/mapping_rerun/bed/")
# setwd("E:/Google Drive/koshi_revision/riboprof/")
# setwd("E:/koshi_codon/agami/")
# check annotation gtf and transcript fa
# anH<-create_annotation("gencode.v19.annotation.gtf", dataSource="gencode.v19", organism="Homo sapiens")
# sequences_biost <- Biostrings::readDNAStringSet("hg19bt1.transcripts.fa",
#                                                 format = "fasta", use.names = TRUE)
# length(names(sequences_biost) %in% anH$transcript)
# nrow(an[names(sequences_biost) %in% anH$transcript,])
# length(sequences_biost[as.character(anH$transcript)])

# save.image("ribowaltz_human.RData")

# setwd("~/Google Drive/koshi_revision/riboprof/")
setwd("E:/Google Drive/koshi_revision/riboprof/")
load("ribowaltz_human.RData")
# setwd("/Volumes/landthaler/pcp/projects/miha/ulrike_stalling/")
# setwd("/Volumes/landthaler/pcp/projects/miha/HDLBP/ribo/meta/considered/bed_chx/")
# setwd("/Volumes/landthaler/pcp/projects/miha/HDLBP/ribo/meta/considered/bed_nochx/")
setwd("D:/landthaler/HDLBP/ribo/meta/considered/bed_chx//")
# setwd("D:/landthaler/HDLBP/ribo/hdlbp_ip/bed/")

anH<-data.table(anH) #higher versions of ribowaltz work with data table

files<-list.files(getwd(), pattern="*bed*")
files<-gsub("\\.bed","",files)
files<-files[6]

out<-"E:/work/hdlbp/ribo_qc/"
for (i in files ) {
  # i<-"nochx.293_1"
  reads <- bedtolist(paste0(getwd(),"/",i),anH)
  example_length_dist <- rlength_distr(reads, sample = i, cl=100)
  example_length_dist[["plot"]]
  ggsave(paste0(out,"/all_plots_ip/rlength_dist/",i,".pdf" ),width = 6, height = 5, units="in")
  example_ends_heatmap <- rends_heat(reads, anH, sample = i, cl = 99,
                                     utr5l = 25, cdsl = 40, utr3l = 25, log_colour = T)
  example_ends_heatmap[["plot"]]
  ggsave(paste0(out,"all_plots_ip/rends_heat/",i,".pdf" ),width = 8, height = 8, units="in")
  
  psite_offset <- psite(reads, start=TRUE, cl=100, flanking = 6, extremity="auto", plot=T, plot_dir=paste0(out,"all_plots_ip/psite_offset_start"), plot_format="pdf", log_file = T, log_file_dir = paste0(out,"all_plots_ip/psite_offset_start/",i))
  psite_offset <- psite(reads, start=FALSE, cl=100, flanking = 6, extremity="auto", plot=T, plot_dir=paste0(out,"all_plots_ip/psite_offset_stop"), plot_format="pdf", log_file = T, log_file_dir = paste0(out,"all_plots_ip/psite_offset_stop/",i))
   
  reads_psite_list <- psite_info(reads, psite_offset)
  rm(reads)
  gc()
  example_metaprofile<-metaprofile_psite(reads_psite_list, anH, i,
                                         length_range="all",utr5l=20, cdsl=40, utr3l=20, plot_title = i)
  example_metaprofile[[paste0("plot_",i)]]
  
  ggsave(paste0(out,"all_plots_ip/metaprofile/",i,".pdf" ),width = 8, height = 6, units="in")
  
  comparison_dt<-list()
  comparison_dt[[paste0(i,"_21nt")]]<-reads_psite_list[[i]][length==21]
  comparison_dt[[paste0(i,"_22nt")]]<-reads_psite_list[[i]][length==22]
  comparison_dt[[paste0(i,"_23nt")]]<-reads_psite_list[[i]][length==23]
  comparison_dt[[paste0(i,"_24nt")]]<-reads_psite_list[[i]][length==24]
  comparison_dt[[paste0(i,"_25nt")]]<-reads_psite_list[[i]][length==25]
  comparison_dt[[paste0(i,"_26nt")]]<-reads_psite_list[[i]][length==26]
  comparison_dt[[paste0(i,"_27nt")]]<-reads_psite_list[[i]][length==27]
  comparison_dt[[paste0(i,"_28nt")]]<-reads_psite_list[[i]][length==28]
  comparison_dt[[paste0(i,"_29nt")]]<-reads_psite_list[[i]][length==29]
  comparison_dt[[paste0(i,"_30nt")]]<-reads_psite_list[[i]][length==30]
  comparison_dt[[paste0(i,"_31nt")]]<-reads_psite_list[[i]][length==31]
  comparison_dt[[paste0(i,"_32nt")]]<-reads_psite_list[[i]][length==32]
  names_list<-list(paste0(i,"_21nt"),
                   paste0(i,"_22nt"),
                   paste0(i,"_23nt"), 
                   paste0(i,"_24nt"),
                   paste0(i,"_25nt"),
                   paste0(i,"_26nt"),
                   paste0(i,"_27nt"),
                   paste0(i,"_28nt"),
                   paste0(i,"_29nt"),
                   paste0(i,"_30nt"),
                   paste0(i,"_31nt"),
                   paste0(i,"_32nt"))
  names(names_list)<-c(paste0(i,"_21nt"),
                       paste0(i,"_22nt"), 
                       paste0(i,"_23nt"),
                       paste0(i,"_24nt"),
                       paste0(i,"_25nt"),
                       paste0(i,"_26nt"),
                       paste0(i,"_27nt"),
                       paste0(i,"_28nt"),
                       paste0(i,"_29nt"),
                       paste0(i,"_30nt"),
                       paste0(i,"_31nt"),
                       paste0(i,"_32nt"))
  example_metaheatmap<-metaheatmap_psite(comparison_dt, anH, sample=names_list,
                                         utr5l=20, cdsl=40, utr3l=20, log=T)
  example_metaheatmap[["plot"]]
  ggsave(paste0(out,"all_plots_ip/metaheatmap/",i,".pdf" ),width = 10, height = 7, units="in")
  
  example_frames_stratified<-frame_psite_length(reads_psite_list, sample=i, region="all", cl=100)
  example_frames_stratified[["plot"]]
  
  ggsave(paste0(out,"all_plots_ip/frame_psite_length/",i,".pdf" ),width = 6, height = 7, units="in")
  
  example_frames<-frame_psite(reads_psite_list, sample=i, region="all")
  example_frames[["plot"]]
  
  ggsave(paste0(out,"all_plots_ip/frame_psite/",i,".pdf" ),width = 6, height = 3, units="in")
  
  reads_psite_list<-length_filter(reads_psite_list, "custom", c(29,30,31))
  reads_psite_list<-as.data.frame(reads_psite_list[[i]])
  reads_psite_list<-subset(reads_psite_list, !is.na(psite) & !is.na(psite_from_start) & !is.na(psite_region) )
  reads_psite_list$count<-1
  ag<-aggregate(count~psite+transcript, data=reads_psite_list, sum)
  ag$start<-ag$psite-1
  ag<-ag[,c(2,4,1,3)]
  ag<-ag[order(ag$transcript, ag$start),]
  write.table(ag, paste0(out,"all_plots_ip/tracks/l29-30.psite.",i,".bedgraph"),quote=F, sep="\t", row.names=F, col.names = F)
}



for (i in files ) {
  # i<-"wt293_1"
  reads <- bedtolist(paste0(getwd(),"/",i),anH)
  reads_length<-length_filter(reads, "periodicity", periodicity_threshold = 50)
  rm(reads)
  gc()
  reads_length<-length_filter(reads_length, "custom", c(29,30,31))
  psite_offset <- psite(reads_length, start=FALSE, cl=100, flanking = 6, extremity="auto", plot=F)
  reads_psite_list <- psite_info(reads_length, psite_offset, site=c("psite", "asite", "esite"), fastapath = "D:/landthaler/HDLBP/all_clip_data/reclip/mapping_trans/hg19bt1.transcripts.fa", fasta_genome=F)
  
  cds_cov<-cds_coverage(reads_psite_list, anH, start_nts=6, stop_nts=6)
  cod_cov<-codon_coverage(reads_psite_list, anH, psite=T)
  write.table(cds_cov, "E:/work/hdlbp/ribo_qc/cds_cov_ip/",i,".txt", quote=F, sep="\t", row.names=F)
  write.table(cod_cov, "E:/work/hdlbp/ribo_qc/cds_cov_ip/",i,".txt", quote=F, sep="\t", row.names=F)
}
  
  # setwd("D:/landthaler/HDLBP/ribo/hdlbp_ip/bed/compare/")
  reads <- bedtolist(getwd(),anH)
  reads_length<-length_filter(reads, "periodicity", periodicity_threshold = 50)
  rm(reads)
  gc()
  reads_length<-length_filter(reads_length, "custom", c(30))
  psite_offset <- psite(reads_length, start=FALSE, cl=100, flanking = 6, extremity="auto", plot=F)
  reads_psite_list <- psite_info(reads_length, psite_offset, site=c("psite", "asite", "esite"), fastapath = "D:/landthaler/HDLBP/all_clip_data/reclip/mapping_trans/hg19bt1.transcripts.fa", fasta_genome=F)
  
  rm(reads_length)
  gc()
  
  
  trans<-NULL
  
 
  
  setwd("E:/Google Drive/hdlbp/")
  mas<-read.delim("hdlbp_master_table_with_classes_uniq_tsig.txt", header=T)
  locmem<-subset(mas, gene_biotype=="protein_coding"& tpm_cutoff>=10 & localization_cat=="membrane" & loc_tar_CDS=="membrane_tc>5.56 & tc<65.26", select="transcript")[,1]
  locmem<-subset(mas, gene_biotype=="protein_coding"& tpm_cutoff>=10 & log2FoldChange.ribo.rna.KO.WT<(-0.25) & !is.na(tc_transcript_norm), select="transcript")[,1]
  locmem<-subset(mas, tpm_cutoff>=10 & log2FoldChange.ribo.rna.KO.WT<0 & loc_tar_CDS=="membrane_tc>5.56 & tc<65.26", select="transcript")[,1]
  locmem<-subset(mas, tpm_cutoff>=10 & log2FoldChange.ribo.rna.KO.WT<0 & loc_tar_CDS=="membrane_tc>5.56 & tc<65.26", select="transcript")[,1]
  locmem<-subset(mas, gene_biotype=="protein_coding"& tpm_cutoff>=10 & log2FoldChange.ribo.rna.KO.WT<0 , select="transcript")[,1]
  
  #this kind of works
  locmem<-subset(mas, gene_biotype=="protein_coding"& tpm_cutoff>=10 & tsig=="mem_notsig" & !is.na(tc_transcript_norm), select="transcript")[,1]
  
  
  trans<-locmem
  trans<-NULL
  i<-"wt293_2"
  j<-"guide2_2"
  
  cod_us<-codon_usage_psite(reads_psite_list, anH, sample=c(i,j), site="esite", 
                            fastapath = "D:/landthaler/HDLBP/all_clip_data/reclip/mapping_trans/hg19bt1.transcripts.fa", 
                            fasta_genome=F, transcripts=trans, label_scatter = T, label_number = 64, label_aminoacid = F, frequency_normalization = T)
  cod_us[["plot_comparison"]]
  
  df1<-cod_us[["dt"]][,c(1,2,8,9 )]
  df2<-cod_us[["dt"]][,c(1,2,8,9 )]
  df3<-cod_us[["dt"]][,c(1,2,8,9 )]
  df4<-cod_us[["dt"]][,c(1,2,8,9 )]
  
  df<-cbind(df1, df2[,4],df3[,3:4],df4[,4])
  mel<-melt(df, measure.vars = colnames(df)[3:ncol(df)], id.vars = c("codon", "aa"))
  ggplot(mel, aes(codon, value, colour=variable))+geom_point()+coord_flip()
  
  df$wt1_g1_1<-df$plot_value_guide1_1-df$plot_value_wt293_1
  df$wt1_g2_1<-df$plot_value_guide2_1-df$plot_value_wt293_1
  df$wt2_g1_2<-df$plot_value_guide1_2-df$plot_value_wt293_2
  df$wt2_g2_2<-df$plot_value_guide2_2-df$plot_value_wt293_2
  
  mel<-melt(df, measure.vars = colnames(df)[9:ncol(df)], id.vars = c("codon", "aa"))
  
  ggplot(mel, aes(codon, value, colour=variable))+geom_point()+coord_flip()
  ggplot(mel, aes(codon, value, colour=codon))+geom_boxplot()+coord_flip()
  
  df$mean<-rowMeans(df[,9:ncol(df)])
  df$sd<-apply(df[,9:ncol(df)], 1, sd)
  
  ggplot(df, aes(codon, mean))+
    facet_grid(~aa, space="free_x", scales="free_x", switch="x")+
    scale_colour_discrete(guide = FALSE)+ 
    theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.2))+
    geom_point(aes(colour=aa))+geom_errorbar(aes(ymin=df$mean-df$sd, ymax=df$mean+df$sd))
  
  
  subset(mas, gene_biotype=="protein_coding"& tpm_cutoff>=10 & tsig=="mem_notsig" , select=c("transcript","Symbol"))
  
  ggplot(mas[mas$tpm_cutoff>=10 & mas$gene_biotype=="protein_coding",], aes(log2FoldChange.ribo.rna.KO.WT, colour=tsig))+stat_ecdf()+coord_cartesian(xlim=c(-.5,.5))
  ggplot(mas[mas$tpm_cutoff>=10 & mas$gene_biotype=="protein_coding",], aes(log2FoldChange.ribo.rna.KO.WT, colour=localization_cat))+stat_ecdf()+coord_cartesian(xlim=c(-.5,.5))+facet_wrap(~tc_transcript_norm_cat)
  
  ggplot(mas[mas$tpm_cutoff>=10 & mas$gene_biotype=="protein_coding",], aes(log2(tc_transcript_norm), colour=tsig))+stat_ecdf()
  ggplot(mas[mas$tpm_cutoff>=10 & mas$gene_biotype=="protein_coding",], aes(log2(utr3_vs_cds1), colour=tsig))+stat_ecdf()
  