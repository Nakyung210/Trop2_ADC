setwd("/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC")

library("DESeq2")
library("tidyverse")

# ===== Import peak count data from FeatureCounts ==============================
## merged_8216 is belongs to TACSTD2 region
count_table <- read.csv("./data/count/GSE124228_PeakCount_IDR.csv", 
                        header = TRUE, row.names = NULL)
count <- as.matrix(count_table[,5:ncol(count_table)])
rownames(count) <- count_table$Peak

# Import meta data
meta <- read.csv("./data/meta.csv", header = TRUE, row.names = 1)

# ===== Perform Differential accessibility analysis ============================
dds <- DESeqDataSetFromMatrix(countData = count, colData = meta, design = ~ Condition)
dds$Condition <- relevel(dds$Condition, ref = "ARID1A_WT")
dds <- DESeq(dds)
res <- as.data.frame(results(dds))

res$Peak <- rownames(res)
results <- left_join(res, count_table[,1:4], by = "Peak")

write.csv(results, "./data/DA analysis/GSE124228_DA_IDR.csv")
