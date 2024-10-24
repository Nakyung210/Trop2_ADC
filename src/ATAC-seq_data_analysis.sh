#!/bin/bash
# Programmer: Nakyung
# This pipeline is for pre-processing ATAC-seq data obtained from GSE124228. 

# Directories
fq=/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC/data/fastq/
QC=/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC/data/QC/
trimmed_fq=/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC/data/trimmed_fq/
alignment=/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC/data/alignment/
ref=/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC/reference/
src=/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC/src/
peak=/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC/data/peakCall/
count=/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC/data/ATAC_count/

##############################################################################################
# Step0. Download datasets from the GSE124228. 
cat ${src}SRR_Acc_List_ATACseq.txt | parallel fasterq-dump --outdir ${fq} --temp ${fq} {}

##############################################################################################
# Step1. QC
fastqc ${fq}*.fastq -t 16 -o ${QC}
unzip "${QC}*.zip"
for file in `find ${QC} -name summary.txt`; do ls "${file}"; grep -v PASS "${file}"; done

##############################################################################################
# Step2. Trim adapter
for srrid in SRR83615{35..44}; do trimmomatic PE -threads 16 -phred33 -trimlog \
 ${trimmed_fq}trimlog.txt isummary ${trimmed_fq}summary.txt \
 ${fq}${srrid}_1.fastq ${fq}${srrid}_2.fastq \
 ${trimmed_fq}${srrid}_paired_1.fq ${trimmed_fq}${srrid}_unpaired_1.fq \
 ${trimmed_fq}${srrid}_paired_2.fq ${trimmed_fq}${srrid}_unpaired_2.fq \
 ILLUMINACLIP:${ref}TruSeq3_PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36; 
done

##############################################################################################
# Step3. Alignment

# Align the reads with Bowtie2
for srrid in SRR83615{35..44}; do bowtie2 -p 16 -x ${hg38} \
 -1 ${trimmed_fq}${srrid}_paired_1.fq -2 ${trimmed_fq}${srrid}_paired_2.fq \
 -X2000 --local --mm --no-mixed --no-discordant -S ${alignment}${srrid}.sam; 
done

# Convert to BAM format from SAM format
for file in ${alignment}*.sam; do \
 samtools view -bS ${file} -o ${file}.bam
 samtools sort -@ 16 ${file}.bam -o ${file}_sorted.bam 
 samtools index -@ 16 ${file}_sorted.bam;
done

##############################################################################################
# Step4. Post-alignment QC

## Remove reads aligned to MT
for srrid in SRR83615{35..44}; do 
 samtools view -h ${aligned}${srrid}_sorted.bam | grep -v chrM | samtools sort -O bam -o ${alignment}${srrid}.rmM.bam; 
done

## Remove PCR duplicates
for file in ${alignment}*.rmM.bam; do MarkDuplicates \
 -I ${file} \
 -O ${file/.rmM.bam/noDup.bam} \
 -M ${alignment}${file}_matrics.txt \
 -REMOVE_DUPLICATES True;
done

## Remove ENCODE blacklist regions with bedtools and samtools
for srrid in SRR83615{35..44}; do bedtools intersect -nonamecheck -v -abam ${alignment}${srrid}noDup.bam -b ${ref}hg38-blacklist.v2.bed > ${alignment}${srrid}.tmp.bam ; done
for srrid in SRR83615{35..44}; do 
 samtools sort -@ 16 -O bam -o ${alignment}${srrid}.rmBL.bam ${alignment}${srrid}.tmp.bam
 samtools index -@ 16 ${alignment}${srrid}.rmBL.bam;
done

## Correct Tn5 shift with deeptools
for srrid in SRR83615{35..44}; do 
 alignmentSieve --numberOfProcessors max --ATACshift --bam ${alignment}${srrid}.rmBL.bam -o ${alignment}${srrid}.shifted.bam
 # If you want to generate BigWig file, please sort and index using commands below
 #samtools sort -@ 16 -O bam ${alignment}${srrid}.shifted.bam -o ${alignment}${srrid}.shifted.sorted.bam
 #samtools index ${alignment}${srrid}.shifted.sorted.bam; 
done

##############################################################################################
# Step5. Generate BigWig files

for srrid in SRR83615{35..44}; do
 samtools sort -@ 16 ${alignment}${srrid}.shifted.bam -o ${alignment}${srrid}.shifted.sorted.bam;
 samtools index -@ 16 ${alignment}${srrid}.shifted.sorted.bam; 
 bamCoverage --numberOfProcessors 16 --binSize 1 --normalizeUsing BPM --effectiveGenomeSize 2805636231 --bam ${alignment}${srrid}.shifted.sorted.bam -o ${alignment}${srrid}.bw ;
done

##############################################################################################
# Step6. Peak calling

## Merge the replicates in each group
samtools merge -@ 16 ${alignment}ARID1A_WT.bam ${alignment}SRR83615{35..38}.shifted.sorted.bam
samtools merge -@ 16 ${alignment}ARID1A_KO.bam ${alignment}SRR83615{39..44}.shifted.sorted.bam

## Call peaks
for file in ARID1A_WT ARID1A_KO ; do echo "macs2 call peak -f BAMPE -g 2805636231 -t ${alignment}${file}.bam --extsize 150 --shift -75 --slocal 5000 --llocal 20000 -B --keep-dup all -p 0.05 -n ${file} --outdir ${peak}" > ${src}PeakCall_MACS2_1.txt 
for file in SRR83615{35..44}; do echo "macs2 call peak -f BAMPE -g 2805636231 -t ${alignment}${file}.bam --extsize 150 --shift -75 --slocal 5000 --llocal 20000 -B --keep-dup all -p 0.01 -n ${file} --outdir ${peak}" > ${src}PeakCall_MACS2_2.txt 
cat ${src}PeakCall_MACS2_1.txt | parallel -j$(nproc)
cat ${src}PeakCall_MACS2_2.txt | parallel -j$(nproc)
for file in $(find ${peak} -name *.narrowPeak); do filename=$(basename ${file} _peak.narrowPeak).sorted.narrowPeak; sort -k8,8nr ${file} > ${peak}${filename}; done

## Check the reproducibility of peaks
for file in SRR83615{35..38}; do echo "idr --sample ${peak}${file}_sorted.narrowPeak ${peak}ARID1A_WT_sorted.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${peak}${file}-idr --plot --log-output-file ${peak}${file}.idr.log"; done > ${src}IDR_WT.txt
for file in SRR83615{39..44}; do echo "idr --sample ${peak}${file}_sorted.narrowPeak ${peak}ARID1A_KO_sorted.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${peak}${file}-idr --plot --log-output-file ${peak}${file}.idr.log"; done > ${src}IDR_KO.txt
cat ${src}IDR_WT.txt | parallel -j$(nproc)
cat ${src}IDR_KO.txt | parallel -j$(nproc)
for file in SRR83615{35..44}; do echo "cat ${peak}${file}-idr | awk '{if($5 >= 540) print $0}' | sort -k1,1 -k2,2n > ${peak}${file}_filtered_idr"; done > ${src}idr_filter.txt
cat ${src}idr_filter.txt | parallel -j$(nproc)

##############################################################################################
# Step7. Obtain signal density

## Run ComputeMatrix
computeMatrix reference-point -S ${alignment}SRR83615{35..42}.bw -R ${ref}Trop2_TSS.bed -a 1000 -b 200 -bs 20 --skipZeros --averageTypeBins sum -out ${count}GSE124228_ATAC_TSS_bin20_mtx
### For the visualization, proceed to Metaplot_TACSTD2.R ###
