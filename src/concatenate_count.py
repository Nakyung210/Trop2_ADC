#!/usr/bin/python3 
# Programmer: Nakyung
# This code is to concatenate all the count tables generated from the STAR counting. 

from pathlib import Path
import csv
import sys

# First argument: paired-end read (PE) or single-end read (SE)
# Second argument: full path to the count files (Please end the path with "/")
# Third argument: is the count TE? Gene?
# Forth argument: full path to the output directory (please include (filename).csv)

if sys.argv[1] == "SE":
    read_end = 1
else:
    read_end = 3

wd = Path(sys.argv[2])

if sys.argv[3] == "Gene":
    countType = "Ensembl_ID"
else:
    countType = sys.argv[3]
    
outputPath = sys.argv[4]

# Concatenate the count files

desired_column = int(read_end)
sample_dict = {}
sample_names = []

for file in wd.glob("*ReadsPerGene.out.tab"):
    file_name = Path(file.stem.replace("ReadsPerGene.out", ""))
    sample_names.append(file_name)
    gene_dict = {}
    with open(file) as tabfile:
        reader = csv.reader(tabfile, delimiter = "\t")
        for row in reader:
            gene_dict[row[0]]=row[desired_column]
    sample_dict[file_name] = gene_dict

with open(outputPath, "wt") as counts_file:
    counts_writer = csv.writer(counts_file)
    counts_writer.writerow([countType]+ sample_names)
    sorted_genes = sorted(list(sample_dict[sample_names[0]].keys()))
    for gene in sorted_genes:
        output = [gene]
        for sample in sample_names:
            output.append(sample_dict[sample][gene])
        counts_writer.writerow(output)

