library(riboWaltz)
library(data.table)
library(ggplot2)

# prepare annotation for hg19

# check annotation gtf and transcript fasta
# anH<-create_annotation(gtfpath = "gencode.v19.primary_assembly.annotation.gtf", dataSource="gencode.v19", organism="Homo sapiens")
# sequences_biost <- Biostrings::readDNAStringSet("hg38star.transcripts.fa",

# save.image("ribowaltz_human.RData")

# load("ribowaltz_human.RData")

# set local directory with stored R session annotation
setwd("D:/Documents/hdlbp_git/")
load("ribowaltz_human.RData")
anH<-data.table(anH) #higher versions of ribowaltz work with data table

# set local directory containing bed files from 
# transcripme mapped bam files (bamToBed) 
setwd("D:/Documents/hdlbp_git/bed_original/")

files<-list.files(getwd(), pattern="*bed*")
files<-gsub("\\.bed","",files)

# Set output directory where QC output files will go
out<-"D:/Documents/hdlbp_git/original_riboQC/"

reads <- bedtolist(getwd(),anH)

psite_offset <- psite(reads, start=TRUE, cl=100, flanking = 6, extremity="auto", 
                      plot=F, plot_dir=paste0(out,"/psite_offset_start"), plot_format="pdf", log_file = F, log_file_dir = paste0(out,"/psite_offset_start/"))

reads_psite_list <- psite_info(reads, psite_offset)

#figS5a
for (i in names(reads)){
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
  
  ggsave(paste0(out,"/metaheatmap/",i,".pdf" ),width = 10, height = 7, units="in")
}

#figS5b
for (i in names(reads)){
  example_frames<-frame_psite(reads_psite_list, sample=i, region="all")
  example_frames[["plot"]]
  
  ggsave(paste0(out,"/frame_psite/",i,".pdf" ),width = 6, height = 3, units="in")
}



for (i in names(reads)){
  example_ends_heatmap <- rends_heat(reads, anH, sample = i, cl = 99,
                                     utr5l = 25, cdsl = 40, utr3l = 25, log_colour = T)
  example_ends_heatmap[["plot"]]
  ggsave(paste0(out,"/rends_heat/",i,".pdf" ),width = 8, height = 8, units="in")
}

for (i in names(reads)){
example_metaprofile<-metaprofile_psite(reads_psite_list, anH, i,
                                       length_range="all",utr5l=20, cdsl=40, utr3l=20, plot_title = i)
example_metaprofile[[paste0("plot_",i)]]

ggsave(paste0(out,"/metaprofile/",i,".pdf" ),width = 8, height = 6, units="in")
}

for (i in names(reads)){

example_frames_stratified<-frame_psite_length(reads_psite_list, sample=i, region="all", cl=100)
example_frames_stratified[["plot"]]

ggsave(paste0(out,"/frame_psite_length/",i,".pdf" ),width = 6, height = 7, units="in")

example_frames<-frame_psite(reads_psite_list, sample=i, region="all")
example_frames[["plot"]]

ggsave(paste0(out,"/frame_psite/",i,".pdf" ),width = 6, height = 3, units="in")
}

for (i in names(reads)){
  example_length_dist <- rlength_distr(reads, sample = i)
  example_length_dist[[paste0("plot_",i)]]
  ggsave(paste0(out,"/length/",i,".pdf" ),width = 6, height = 3, units="in")
}


length_filtered<-length_filter(reads_psite_list, length_filter_mode = "custom", length_range = 28:38)

cds_coverage_example <- cds_coverage(length_filtered, anH)

# write.table(cds_coverage_example, "../novaseq_ribo_count_table_cds.txt", quote=F, sep="\t", row.names=F,col.names = T)

