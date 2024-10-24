#!/bin/bash
# Programmer: Nakyung
# This pipeline is for pre-processing ATAC-seq data for motif analysis.

# Directories
alignment=/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC/data/alignment/
ref=/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC/reference/
src=/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC/src/
peak=/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC/data/peakCall/
count=/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC/data/ATAC_count/
bed=/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC/data/bed/
fasta=/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC/data/fasta/
FIMO=/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC/data/FIMO/

##############################################################################################
# Step1. Screen the motif

## Run FIMO
cat ${peak}IDR_merged.bed | sed '1d' | awk -F"\t" '{OFS="\t";print $1,$2,$3}' > /${peak}IDR_merged2.bed
cat ${peak}IDR_merged2.bed | xargs -L 1 -P 100 -I {} bash -c 'name=$(echo "{}" | awk "{OFS=\"_\"; print \$1, \$2, \$3}"); echo "{}" > ${bed}line_by_line_IDR/${name}.bed'

export FASTA="${ref}GRCh38.primary_assembly.genome.fa"
export OUT_DIR="${fasta}line_by_line_IDR/"
find ${bed}line_by_line_IDR -name "*.bed" | parallel -j$(nproc) 'name=$(basename {} .bed); bedtools getfasta -fi ${FASTA} -bed {} -name > ${OUT_DIR}/${name}.fasta'
find ${fasta}line_by_line_IDR/ -name *.fasta | while read file; do name=$(basename ${file} .fasta); echo "fimo --o ${FIMO}line_by_line/${name} ${ref}CIS-BP_2.00/Homo_sapiens.meme ${fasta}line_by_line_IDR/${name}.fasta"; done > ${src}fimo_idr.txt

split -l 21000 -d -a 1 --additional-suffix=.txt ${src}fimo_idr.txt ${src}fimo_idr

cat ${src}fimo_idr0.txt | parallel -j$(nproc)
cat ${src}fimo_idr1.txt | parallel -j$(nproc)
cat ${src}fimo_idr2.txt | parallel -j$(nproc)
cat ${src}fimo_idr3.txt | parallel -j$(nproc)
cat ${src}fimo_idr4.txt | parallel -j$(nproc)
cat ${src}fimo_idr5.txt | parallel -j$(nproc)
cat ${src}fimo_idr6.txt | parallel -j$(nproc)
cat ${src}fimo_idr7.txt | parallel -j$(nproc)
cat ${src}fimo_idr8.txt | parallel -j$(nproc)
cat ${src}fimo_idr9.txt | parallel -j$(nproc)

##############################################################################################
# Step2. Post-processing

## Add column name
echo $(cat ${ref}TF_list_GSE124228.txt) | sed -e 's/ /,/g' > ${FIMO}header_idr.csv

## Process the FIMO output
for file in $(find ${FIMO}line_by_line -name fimo.tsv); do fileid=$(echo ${file} | awk -F"/" '{print $9}'); cat ${file} | sed '1d' | awk -F"\t" '{OFS=","; print $2}' | sed 's/(\([^)]*\)).*/\1/' | sort | uniq > ${FIMO}TF_by_line/${fileid}.csv ; done

## Make binary file based on FIMO results
export TF_LIST="${ref}TF_list_GSE124228.txt"
export OUT_DIR="${FIMO}Unique_by_line/"
find ${FIMO}TF_by_line/ -name *.csv | parallel -j$(nproc) '
    fileid=$(basename {} .csv)
    cat ${TF_LIST} | while read line; do
        if grep -qF $line {}; then
            echo "$line,1"
        else
            echo "$line,0"
        fi
    done > ${OUT_DIR}/${fileid}.csv
'
find ${FIMO}Unique_by_line/ -name "*.csv" | parallel -j$(nproc) 'filename=$(basename {} .csv); sed -i "1s/^/TF_list,${filename}\n/" {}'

## Transpose the output file
export OUT_DIR2="${FIMO}Transposed_by_line"
find ${FIMO}Unique_by_line/ -name "*.csv" | parallel -j$(nproc) '
    filename=$(basename {} .csv);
    cat {} | awk -F"," '"'"'{for (i=1; i<=NF; i++)  { a[NR,i] = $i } } NF>p { p = NF };END {OFS=",";for(j=1; j<=p; j++) { str=a[1,j]; for(i=2; i<=NR; i++){ str=str OFS a[i,j]; }; print str}}'"'"' | sed '"'"'1d'"'"' > ${OUT_DIR2}/${filename}.csv
'
## Concatenate the trasposed files
find ${FIMO}Transposed_by_line/ -name "*.csv" -type f | xargs -n 100 cat >> ${FIMO}tmp_idr.csv
cat ${FIMO}header_idr.csv ${FIMO}tmp_idr.csv > ${FIMO}GSE124228_binary_idr.csv

### For the downstream analysis, proceed to motif_differential_accessibilility.R ###

##############################################################################################
# Step3. Narrow down the candidates

cat ${FIMO}chr1_58576426_58577987/fimo.tsv | awk -F"\t" '{if($4 >= 58577099 && $5 <= 58577588) {print $0}}' > ${FIMO}Trop2.tsv
cat ${FIMO}Trop2.tsv | awk -F"\t" '{print $2}' | sed 's/(\([^)]*\)).*/\1/' | sort | uniq > ${FIMO}Trop2_TF.txt
comm -12 ${FIMO}Trop2_TF.txt <(sort ${ref}TF_list_GSE124228.txt) > ${FIMO}Trop2_common.tsv
cat ${FIMO}Trop2_common.tsv | while read line; do cat ${FIMO}Trop2.tsv | grep -w ${line} ; done > ${FIMO}Trop2_common_TF.tsv
cat ${FIMO}Trop2_common_TF.tsv | awk -F"\t" '{print $2}' | sort | uniq > ${FIMO}common_TF_list.txt