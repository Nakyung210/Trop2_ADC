#!/bin/bash
# Programmer: Nakyung 
# This pipeline is for pre-processing RNA-seq data obtained from GSE124228. 

# Directories
fq=/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC/data/RNA/fastq/
QC=/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC/data/RNA/QC/
trimmed_fq=/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC/data/RNA/trimmed_fq/
alignment=/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC/data/RNA/alignment/
ref=/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC/reference/gencode.v45/
src=/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC/src/

##############################################################################################
# Step0. Download datasets from the GSE124228. 
cat ${src}SRR_Acc_List.txt | parallel fasterq-dump --outdir ${fq} --temp ${fq} {}

##############################################################################################
# Step1. QC
fastq ${fq}*.fastq -t 16 -o ${QC}
unzip "${QC}*.zip"
for file in `find ${QC} -name summary.txt`; do ls "${file}"; grep -v PASS "${file}"; done

##############################################################################################
# Step2. Trim raw reads
for srrid in SRR8361{592..601}; do trim_galore --paired ${fq}${srrid}_1.fastq ${fq}${srrid}_2.fastq -o ${trimmed_fq}; done

##############################################################################################
# Step3. Mapping and counting
for srrid in SRR8361{592..601}; do
    STAR --runMode alignReads \
         --genomeDir ${ref} \
         --readFilesIn ${trim_fq}"${srrid}"_1_val_1.fq ${trim_fq}"${srrid}"_2_val_2.fq \
         --quantMode GeneCounts \
         --runThreadN 16 \
         --outSAMtype BAM Unsorted \
         --outSAMattributes Standard \
         --outSAMunmapped None \
         --outBAMsortingBinsN 200 \
         --outFileNamePrefix ${alignment}"${srrid}"
done

##############################################################################################
# Step4. Concatanate the read counts
python ${src}concatenate_count.py PE ${alignment} Gene ${alignment}GSE124228_RNA_raw_count.csv