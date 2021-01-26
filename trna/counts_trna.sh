#!/bin/bash
#$ -V
#$ -l h_vmem=12G,h_rt=12:0:0,data
#$ -o /fast/AG_Landthaler/miha/hdlbp/reclip/counts_trna
#$ -e /fast/AG_Landthaler/miha/hdlbp/reclip/counts_trna

basepath="/fast/AG_Landthaler/miha/hdlbp/reclip/mapping_trna2/"

prefix=""
outpath="/fast/AG_Landthaler/miha/hdlbp/reclip/counts_trna/"
bam=".bam"
genome="/fast/AG_Landthaler/genomes/"

cd $outpath

samtools view -F 16 $basepath/tRNA_$filename$bam | cut -f 3 | uniq -c | awk '{print $2"\t"$1}' > tRNA_"$filename".readcounts

samtools view -bF 0x0010 $basepath/tRNA_$filename$bam | samtools mpileup -d 1000000 -Bf $genome/tRNAgt/extended_uniq_mature_tRNAs.fa /dev/stdin | ./row_mpile_coverage_plus_TC.pl | awk '($3=="T" && $5>0){print $1"\t"($2-1)"\t"$2"\tTC\t"$5"\t+";}' | sort -k1,1 -k2,2n > tRNA_"$filename"_TC.bed

awk '($2>2){print $1,$2-3,$2+4,$4,$5,$6}' OFS="\t" tRNA_"$filename"_TC.bed | bedtools getfasta -fi /fast/AG_Landthaler/genomes/tRNAgt/extended_uniq_mature_tRNAs.fa -bed stdin -tab | awk '(gsub(":","\t",$1))' OFS="\t"| awk '(gsub("-","\t",$2))' OFS="\t" | paste - <(awk '($2>2){print $1,$2-3,$2+4,$4,$5,$6,$2,$3}' OFS="\t" tRNA_"$filename"_TC.bed) | awk '($1==$5 && $2==$6)' | awk '{print $1,$11,$12,$8,$9,$10,$4}' OFS="\t" | multiBamCov -bams $basepath/tRNA_$filename$bam -bed stdin -s -D > tRNA_"$filename"_TCseq.bed

intersectBed -s -a tRNA_"$filename"_TCseq.bed -b tRNA_D.bed -v -nonamecheck > tRNA_"$filename"_TCseq_filtered.bed

echo "track type=bedGraph name=\"${filename}\" db=tRNA_mat_ext" > tRNA_"${filename}"_TC.bedgraph

awk 'BEGIN { OFS="\t" } ($6=="+" && $5>=1){print $1,$2,$3,$5}' tRNA_"$filename"_TCseq_filtered.bed >> tRNA_"${filename}"_TC.bedgraph

echo "track type=bedGraph name=\"${filename}\" db=tRNA_mat_ext" > tRNA_"${filename}"_score.bedgraph

awk 'BEGIN { OFS="\t" } ($6=="+" && $5>=1){print $1,$2,$3,$5/$8}' tRNA_"$filename"_TCseq_filtered.bed >> tRNA_"${filename}"_score.bedgraph

awk '{a[$1]+=$5;}END{for(i in a)print i"\t"a[i];}' tRNA_"$filename"_TCseq_filtered.bed > tRNA_"$filename".TCcounts
