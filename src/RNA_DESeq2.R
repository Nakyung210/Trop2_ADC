setwd("/Users/nakyung/Desktop/Bioinformatics/Gilead/RNA/")

library("DESeq2")
library("tidyverse")
library("ggthemes")
library("ggpubr")
library("rtracklayer")
library("ggsignif")

# ===== Import data ============================================================
count_df <- read.csv("./data/GSE124228_RNAseq_raw_counts.csv", header = TRUE,
                  row.names = 1)
meta <- read.csv("./RNA_meta.csv", header = TRUE, row.names = 1)
human_gtf <- readGFF("/Users/nakyung/R/Annotation/gencode.v45.annotation.gtf")

# ===== Perform differential expression analysis ===============================
protein_DE <- function(meta_input, count_df, control, condition, gtf_input, output_path){
  # Preprocess the count file for the DE analysis
  count <- as.matrix(count_df) + 1
  count <- count[,sort(colnames(count))]
  
  # Make DESeq object
  dds <- DESeqDataSetFromMatrix(countData = count, colData = meta_input, 
                                design = as.formula(noquote(paste("~", condition))))
  
  # Add the control and treatment information
  dds$condition <- relevel(dds$condition, ref = control)
  
  # Run differential expression analysis
  dds <- DESeq(dds)
  
  # Get the results
  res <- results(dds)
  res <- res[order(res$padj),]
  
  # Get the normalized count
  normalized_counts <- counts(dds, normalized = TRUE)
  
  # Output the results of differential expression analysis #
  protein_res <- data.frame(res)
  normalized_counts <- as.data.frame(normalized_counts)
  
  # Preprocess GTF file 
  gtf <- subset(gtf_input, gtf_input$type == "gene")
  gtf <- gtf[,c("gene_id", "gene_name")]
  
  # Add gene annotation to the protein_res and normalized_count to the output file
  protein_res$gene_id <- rownames(protein_res)
  normalized_counts$gene_id <- rownames(normalized_counts)
  protein_res_anno <- left_join(protein_res, gtf, by = "gene_id")
  normalized_counts_anno <- left_join(normalized_counts, gtf, by = "gene_id")
  
  # Export results
  write.csv(protein_res_anno, paste0(output_path, "DEanalysis_results.csv"))
  write.csv(normalized_counts_anno, paste0(output_path, "normalized_count.csv"))
  save.image(paste0(output_path, "DESeq2_results.RData"))
}

# input order: meta_input, count_df, control, condition, gtf_input, output_path
protein_DE(meta, count_df, "Control", "condition", human_gtf, "./data/")

# ===== Visualize the differential expression analysis results =================
protein_res <- read.csv("./data/DEanalysis_results.csv", 
                       header = TRUE, row.names = 1)

# Get Trop2's differential expression analysis result
TACSTD2 <- protein_res[protein_res$gene_name == "TACSTD2",]
TACSTD2 <- drop_na(TACSTD2)
stat.test <- data.frame(.y. = "Fold Change", 
                        group1 = "ARID1A WT", group2 = "ARID1A KO",
                        p = TACSTD2$pvalue, p.sig = "***")
TACSTD2.df <- data.frame(genotype = c("ARID1A WT", "ARID1A KO"),
                         FC = c(1, 2^(TACSTD2$log2FoldChange)))

# Visualize with bar plot
TACSTD2.bar <- ggplot(data = TACSTD2.df, 
       aes(x = factor(genotype, levels = c("ARID1A WT", "ARID1A KO")), y = FC,
           fill = genotype)) +
  geom_bar(stat = "identity", width = 0.7, colour = "black", linewidth = 0.3, alpha = 1) +
  geom_signif(comparisons = list(c("ARID1A WT", "ARID1A KO")), 
              annotations = "*** 2.1e-06", size = 0.5, textsize = 6,
              y_position = 1.6, tip_length = 0, vjust = 0.03) +
  scale_y_continuous(limits = c(0, 1.8),
                     breaks = seq(0,1.8, 0.4),
                     expand = expansion(0))+
  scale_fill_manual(values = c("pink", "steelblue")) +
  labs(x = "", y = "Fold Change") +
  theme_par() + 
  theme(legend.position = "none", 
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(face = "bold", size = 15),
        axis.title.y = element_text(face = "bold", size = 18))

ggsave(TACSTD2.bar, filename = "./figure/Trop2 expression bar plot.pdf",
       width = 6, height = 7, dpi = 300)

save.image("./src/DESeq2.RData")




















