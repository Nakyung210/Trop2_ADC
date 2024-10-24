#!/bin/bash
# Programmer: Nakyung
# This pipeline is for pre-processing ATAC-seq data for differential accessibility analysis.

# Directories
alignment=/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC/data/alignment/
ref=/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC/reference/
src=/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC/src/
peak=/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC/data/peakCall/
count=/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC/data/ATAC_count/

##############################################################################################
# Step1. Count the signal

## Merge bed file
bedtools multiinter -header -names SRR83615{35..44} -i ${peak}SRR83615{35..44}_filtered_idr > ${peak}IDR_merged.bed
awk -F $'\t' 'BEGIN {OFS = FS}{ $2=$2+1; peakid="merged_"++nr;  print peakid,$1,$2,$3,"."}' ${peak}IDR_merged.bed > ${peak}IDR_merged.saf

## Count the signal
for file in SRR83615{35..44}; do echo "featureCounts -p -F SAF -O -a ${peak}IDR_merged.saf -o ${count}${file}.counts ${alignment}${file}.shifted.sorted.bam"; done > ${src}idr_count.txt
cat ${src}idr_count.txt | parallel -j$(nproc)

## Post-processing
echo Peak Chr Start End SRR83615{35..44} | sed -e 's/ /,/g' > ${count}header.csv
for srrid in SRR83615{35..44}; do cat ${count}${srrid}.counts | awk -F"\t" '{OFS=","; print $7}' | sed '1d' | sed '1d' > ${count}${srrid}count.csv; done
cat ${count}SRR8361535.counts | awk -F"\t" '{OFS=","; print $1, $2, $3, $4}' | sed '1d' | sed '1d' > ${count}column_name.csv
paste -d "," ${count}column_name.csv ${count}SRR83615{35..44}count.csv > ${count}count_tmp.csv
cat ${count}header.csv /${count}count_tmp.csv > ${count}GSE124228_PeakCount.csv

rm ${count}header.csv ${count}count_tmp.csv ${count}column_name.csv ${count}SRR83615{35..44}.counts ${count}SRR83615{35..44}.counts.summary

### For the differential accessibility analysis. proceed to differential_accessibility_DESeq2.R

