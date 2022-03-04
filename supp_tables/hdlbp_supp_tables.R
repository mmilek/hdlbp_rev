# supplemental tables

ribote<-read.delim("geo_processed_data/processed_data_riboseq_te.txt", header=T)

names<-read.delim("data/gencode_v19_gene_id_to_gene_name_all.txt", header=F)
colnames(names)<-c("gene_id", "gene_name")

ribote<-merge(ribote, names, by="gene_id")
ribote<-subset(ribote, select=colnames(ribote)[c(1,20:length(colnames(ribote)))])
ribote<-ribote[,c(1,14,2:7)]

# write.table(ribote, "supp_tables/tableS5_riboseq.txt", quote=F, sep="\t", row.names=F)

# clip

clip<-read.delim("geo_processed_data/processed_data_parclip_pergene.txt", header=T)
# write.table(clip, "supp_tables/tableS2_clip.txt", quote=F, sep="\t", row.names=F)

#rnaseq

rna<-read.delim("geo_processed_data/processed_data_rnaseq_fractionation.txt", header=T)
mas<-read.delim("data/hdlbp_master_table_with_classes_uniq.txt", header=T)
inf<-subset(mas, select=c("gene_id","Annotation","Symbol",
                          "tpm_cutoff","localization_cat",
                          "log2FoldChange.mem.cyt.293"))
colnames(inf)<-c("gene_id", "gene_biotype", "Symbol","tpm_wt", "localization", "log2FoldChange.mem.WT.vs.cyt.WT")
rna<-merge(rna, inf, by="gene_id")
rna<-rna[,c(1,26:31,100:104)]
rna<-rna[,c(1,9,8,11,10,2,3)]
rna<-subset(rna, tpm_wt>=1)
colnames(rna)[c(6,7)]<-c("log2FoldChange.mem.wt.vs.cyt.wt", "padj.mem.wt.vs.cyt.wt")

# write.table(rna, "supp_tables/tableS1_rnaseq.txt", quote=F, sep="\t", row.names=F)

#trna
trna<-read.delim("data/hdlbp_tRNA_perPosition.txt", header=T)
ncrnaseq<-subset(trna, select=c("trna", "m1_hi_2ome.readcounts", "m2_hi_2ome.readcounts", "norm.m1_hi_2ome.readcounts", "norm.m2_hi_2ome.readcounts", "codonUsage"))
ncrnaseq$trna_id<-gsub(";.*", "", ncrnaseq$trna)
ncrnaseq<-subset(ncrnaseq, !duplicated(trna_id))
ncrnaseq<-ncrnaseq[,c(7,2:6)]
colnames(ncrnaseq)<-c("trna_id", "ncrnaseq1_read_count", 
                      "ncrnaseq2_read_count", 
                      "norm.ncrnaseq1_read_count",
                      "norm.ncrnaseq2_read_count",
                      "codon_usage")
# write.table(ncrnaseq, "geo_processed_data/processed_data_ncrna_trna.txt", quote=F, sep="\t", row.names = F)

trna$tc_start<-as.numeric(gsub(";.*", "",
                            gsub(".*sep", "",
                                 sub(";", "sep", trna$trna))))
trna$tc_stop<-as.numeric(gsub(".*;", "", trna$trna))
trna$trna_id<-gsub(";.*", "", trna$trna)
supp<-trna[,c(3:5,9:11,15:16,19,20,39,40,43,45:46,49:51,36)]
supp<-supp[,c(18,16,17,1,4,3,6,2,7,8,9,10,11:13,15,14,19)]
colnames(supp)[16:17]<-c("ncrnaseq_trna_expression", "trna_hdlbp_enrichment_rank")
colnames(supp)[8]<-"seq_tc_middle_t"
# write.table(supp, "supp_tables/tableS7_trna.txt", quote=F, sep="\t", row.names = F)

trnageo<-supp[,1:15]
# write.table(trnageo, "geo_processed_data/processed_data_parclip_trna.txt", quote=F, sep="\t", row.names = F)

#parclip binding sites
pc<-read.delim("geo_processed_data/processed_data_parclip_overlapping_clusters.txt", header=T)
pc<-pc[,c(1:18)]
colnames(pc)[18]<-"region"
# write.table(pc, "geo_processed_data/processed_data_parclip_overlapping_clusters.txt", quote=F, sep="\t", row.names = F)
# write.table(pc, "supp_tables/tableS3_clip_sites.txt", quote=F, sep="\t", row.names = F)


# psilac
# done in psilac_analysis.R

##parclip per gene done in hdlbp_clip.R
##riboseq per gene done in hdlbp_ribo.R
##bioid done in hdlbp_bioid.R


