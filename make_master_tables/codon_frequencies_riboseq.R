# this allows to compute codon frequencies from transcriptome mapped
# bed files (bamToBed converted)

library(riboWaltz)
library(data.table)
library(ggplot2)

# create annotation for transcriptome
# anH<-create_annotation("gencode.v19.annotation.gtf", dataSource="gencode.v19", organism="Homo sapiens")
# sequences_biost <- Biostrings::readDNAStringSet("hg19bt1.transcripts.fa",
                                              # format = "fasta", use.names = TRUE)
# length(names(sequences_biost) %in% anH$transcript)
# nrow(an[names(sequences_biost) %in% anH$transcript,])
# length(sequences_biost[as.character(anH$transcript)])

# save.image("ribowaltz_human.RData")

setwd("D:/Documents/hdlbp_git/")
load("ribowaltz_human.RData")

setwd("D:/Documents/hdlbp_git/bed_original/")

anH<-data.table(anH) #higher versions of ribowaltz work with data table
files<-c("wt293_1", "wt293_2", "guide1_1", "guide1_2", "guide2_1", "guide2_2")

out<-"D:/Documents/hdlbp_git/codons_original/"
for (i in files ) {
  # i<-"wt293_1"
  reads <- bedtolist(paste0(getwd(),"/",i),anH)
  example_length_dist <- rlength_distr(reads, sample = i, cl=100)
  example_length_dist[[paste0("plot_",i)]]
  ggsave(paste0(out,"rlength_dist/",i,".pdf" ),width = 6, height = 5, units="in")
  
  psite_offset <- psite(reads, start=TRUE, cl=100, flanking = 6, extremity="auto", plot=F, plot_dir=paste0(out,"psite_offset_start"), plot_format="pdf", log_file = T, log_file_dir = paste0(out,"psite_offset_start/",i))

  reads_psite_list <- psite_info(reads, psite_offset)
  rm(reads)
  gc()
  example_metaprofile<-metaprofile_psite(reads_psite_list, anH, i,
                                         length_range="all",utr5l=20, cdsl=40, utr3l=20, plot_title = i)
  example_metaprofile[[paste0("plot_",i)]]
  example_frames<-frame_psite(reads_psite_list, sample=i, region="all")
  example_frames[["plot"]]
  
  ggsave(paste0(out,"/frame_psite/",i,".pdf" ),width = 6, height = 3, units="in")
  sites<-c("psite", "asite", "esite")
  for (site in sites){
    cod_us<-codon_usage_psite(reads_psite_list, anH, sample=i, site=site, 
                            fastapath = "D:/Documents/hdlbp_git/hg19bt1.transcripts.fa", 
                            fasta_genome=F, transcripts=NULL, label_scatter = T, label_number = 64, label_aminoacid = F, frequency_normalization = F)
    write.table(cod_us[["dt"]], paste0(out,i,"_",site,"_no_norm.txt"), quote=F, sep="\t", row.names=F) 
    cod_us<-codon_usage_psite(reads_psite_list, anH, sample=i, site=site, 
                              fastapath = "D:/Documents/hdlbp_git/hg19bt1.transcripts.fa", 
                              fasta_genome=F, transcripts=NULL, label_scatter = T, label_number = 64, label_aminoacid = F, frequency_normalization = T)
    write.table(cod_us[["dt"]], paste0(out,i,"_",site,"_norm.txt"), quote=F, sep="\t", row.names=F) 
  }


}


# if you are interested only in a transcript subset

mas<-read.delim("data/hdlbp_master_table_with_classes_uniq.txt", header=T)
inf<-subset(mas, select=c("gene_id","Symbol",
                          "tpm_cutoff","tc_CDS_norm","localization_cat",
                          "mean_te_293","loc_tar_CDS","tc_CDS_norm_cat",
                          "tc_transcript_norm",
                          "log2FoldChange.mem.cyt.293",
                          "log2FoldChange.mem.cyt.KO.293",
                          "log2FoldChange.ribo.rna.KO.WT",
                          "Transmembrane.helices" ,"Cleavage.site..Signalp." ))


trans<-read.delim("data/considered_gene_names.txt", header=T)[,c(1:10,16,18)]
colnames(trans)[10]<-"gene_id"

inf<-merge(inf, trans, by="gene_id")

inf<-subset(inf, tpm_cutoff>=10 & localization_cat=="membrane" )

sp<-read.delim("data/sp_list.txt", header=T)
tm<-read.delim("data/tm_list.txt", header=T)

inf<-sp
colnames(inf)[3]<-"transcript"

inf<-tm
colnames(inf)[3]<-"transcript"


anH<-data.table(anH) #higher versions of ribowaltz work with data table
files<-c("wt293_1", "wt293_2", "guide1_1", "guide1_2", "guide2_1", "guide2_2")

out<-"D:/Documents/hdlbp_git/codons_subset_tm/"
for (i in files ) {
  #i<-"guide2_2"
  reads <- bedtolist(paste0(getwd(),"/",i),anH)
  
  psite_offset <- psite(reads, start=TRUE, cl=100, flanking = 6, extremity="auto", plot=F, plot_dir=paste0(out,"psite_offset_start"), plot_format="pdf", log_file = T, log_file_dir = paste0(out,"psite_offset_start/",i))
  
  reads_psite_list <- psite_info(reads, psite_offset)
  rm(reads)
  gc()
  sites<-c("psite", "asite", "esite")
  for (site in sites){
    cod_us<-codon_usage_psite(reads_psite_list, anH, sample=i, site=site, 
                              fastapath = "D:/Documents/hdlbp_git/hg19bt1.transcripts.fa", 
                              fasta_genome=F, transcripts=inf$transcript, label_scatter = T, label_number = 64, label_aminoacid = F, frequency_normalization = F)
    write.table(cod_us[["dt"]], paste0(out,i,"_",site,"_no_norm.txt"), quote=F, sep="\t", row.names=F) 
    cod_us<-codon_usage_psite(reads_psite_list, anH, sample=i, site=site, 
                              fastapath = "D:/Documents/hdlbp_git/hg19bt1.transcripts.fa", 
                              fasta_genome=F, transcripts=inf$transcript, label_scatter = T, label_number = 64, label_aminoacid = F, frequency_normalization = T)
    write.table(cod_us[["dt"]], paste0(out,i,"_",site,"_norm.txt"), quote=F, sep="\t", row.names=F) 
  }
  
  
}
