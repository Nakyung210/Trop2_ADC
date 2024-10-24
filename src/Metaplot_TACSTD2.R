setwd("/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC")

library("ggplot2")
library("ggthemes")
library("ggpubr")
library("tidyverse")

# ===== Import data ============================================================
read_deeptools_table <- function(file) {
  
  n <- max(count.fields(file, sep = '\t'), na.rm = TRUE)
  x <- readLines(file)
  
  .splitvar <- function(x, sep, n) {
    var <- unlist(strsplit(x, split = sep))
    length(var) <- n
    return(var)
  }
  
  x <- do.call(cbind, lapply(x, .splitvar, sep = '\t', n = n))
  x <- apply(x, 1, paste, collapse = '\t')
  plot_table <- na.omit(read.csv(text = x, sep = '\t')[-1,])  # Remove first row with "gene" label
  
  return(plot_table)
}

bin20 <- read_deeptools_table("./data/metaplot/GSE124228_ATAC_TSS_bin20_mtx_metaplot.txt")

# ===== Visualize with metaplot ================================================
# For the output_dir: Include full path with double quotation like "path/to/save"
meta_plot <- function(input, output_dir){
  # Calculate the average value
  mtx <- input
  mtx$WT <- rowMeans(mtx[,3:6])
  mtx$KO <- rowMeans(mtx[,7:10])
  mtx <- mtx[,c(1,2,11,12)]
  
  # Gather the data
  mtx.long <- mtx %>% pivot_longer(cols = c("WT", "KO"),names_to = "sample",
                                   values_to = "score")
  
  # Visualize the data with metaplot
  metaplot <- ggplot(mtx.long, 
                     aes(x = bins, y = as.numeric(score), color = sample)) +
    geom_line(linewidth = 1.2) +
    scale_x_continuous(breaks = c(1, 10, 60), labels = c("-0.2kb", "TSS", "1kb")) +
    ylab("Normalized signal density") +
    xlab("") +
    theme_par() +
    theme(axis.title.y = element_text(face = "bold", size = 17),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          legend.title = element_text(face = "bold", size = 17),
          legend.text = element_text(size = 16),
          legend.position = c(0.9, 0.8)) +
    scale_color_manual(values = c("WT" = "steelblue", "KO" = "pink")) +
    labs(color = "ARID1A")
  
  # Save the plot
  ggsave(metaplot,filename = output_dir, width = 8, height = 6, dpi = 300)
}

meta_plot(bin20, "./figure/bin20_metaplot.pdf")
save.image("./src/Metaplot_TACSTD2.RData")



