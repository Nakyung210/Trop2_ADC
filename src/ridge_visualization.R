setwd("/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC/")

library("tidyverse")
library("ggthemes")
library("hrbrthemes")
library("RColorBrewer")
library("grid")
library("gridExtra")


# ===== Import data ============================================================
summary_df <- read.csv("./data/FIMO/ridge_summary.csv", row.names = NULL, header = TRUE)
TACSTD2_unique <- read.table("./data/FIMO/common_TF_list.txt")
TACSTD2_unique <- TACSTD2_unique[[1]]
Pearson <- read.csv("/Users/nakyung/Desktop/Bioinformatics/Gilead/TCGA/data/Pearson_correlation.csv", 
                    row.names = NULL, header = TRUE)


# ===== Visualize the ridge regression results =================================
plot_volcano <- function (res, FDR_cutoff, nlabel = 10){
  res <- mutate(res, significance=ifelse(res$FDR<FDR_cutoff, paste0("FDR < ", FDR_cutoff), 
                                         paste0("FDR > ", FDR_cutoff)))
  res = res[!is.na(res$significance),]
  significant_genes <- res %>% filter(significance == paste0("FDR < ", FDR_cutoff))
  top_genes <- significant_genes %>% arrange(FDR) %>% head(nlabel)
  
  
  ggplot(res, aes(cor, -log(FDR))) +
    geom_point(aes(col=significance)) + 
    scale_color_manual(values=c("red", "black")) + 
    ggrepel::geom_text_repel(data=top_genes, aes(label=head(TF,nlabel)), size = 5)+
    scale_y_continuous(limits = c(0,30), 
                       breaks = seq(0,30,10)) +
    labs ( x = "Ridge coefficient", y = expression(-log[10] ~ "FDR", bty = "n"))+
    geom_vline(xintercept = 0.8, linetype = "dotted")+
    geom_vline(xintercept = -0.8, linetype = "dotted")+
    geom_hline(yintercept = -log(0.05), linetype = "dotted")+
    theme_par() +
    theme(legend.position = "none",
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18))
}

volcano <- plot_volcano(summary_df, 0.05, nlabel = 20)
volcano
ggsave(volcano, file = "./figure/volcano.pdf", width = 6.5, height = 7, dpi = 300)


# ===== Visualize the correlation with heatmap =================================
## Preprocess the ridge regression data
data <- summary_df[summary_df$TF %in% TACSTD2_unique,c(1,2,6)]
data <- data %>% filter(cor > 0.8 | cor < -0.8) %>% filter(FDR < 0.05) %>% 
  arrange(desc(cor))
data <- data %>% mutate(label = case_when(data$FDR <= 0.001 ~ "***",
                                          data$FDR <= 0.01 ~ "**",
                                          data$FDR <= 0.05 ~ "*",
                                          data$FDR > 0.05 ~ "ns"))

colnames(data) <- c("TF", "Correlation", "FDR", "label")

saveRDS(data, "./data/FIMO/ridge_visualization.rds")

### Visualize with the heatmap
ridge.plot <- ggplot(data, aes(x = reorder(TF, -Correlation), y = 1, fill = Correlation)) +
  geom_tile(color = "black") +
  scale_fill_distiller(palette = "PiYG") + 
  xlab("    ***       ***        **        **        **         **        **         *         **        **         *         **         *         **         **       ***       ***       ***       ***") +
  ylab("ARID1A") +
  theme_ipsum() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 10, vjust = 33, hjust = 0.007),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.6, 
                                    face = "bold", size = 13),
        plot.margin = unit(c(1,1,0,1), "cm"),
        legend.key.size = unit(0.4,"cm"),
        legend.title = element_text(size = 10)) +
  coord_fixed(ratio = 2) 


## Preprocess the Pearson correlation data
Pearson <- Pearson[Pearson$Column %in% data$TF,]
Pearson <- Pearson %>% mutate(label = 
                                case_when(Pearson$P_value <= 0.001 ~ "***",
                                          Pearson$P_value <= 0.01 ~ "**",
                                          Pearson$P_value <= 0.05 ~ "*",
                                          Pearson$P_value > 0.05 ~ "ns"))

saveRDS(data, "./data/FIMO/ridge_visualization.rds")

### Visualize with the heatmap
pearson.plot <- ggplot(Pearson, aes(x = factor(Column, ARID1A_order),y = 1, fill = Correlation)) +
  geom_tile(color = "black") + 
  scale_fill_distiller(palette = "PiYG") + 
  xlab("     ns       ns       ns        ns       ***       ns        **        ns       **        ns       ***       ns        ***       ns       ns       ***       ns       **        ***") +
  ylab("TACSTD2") +
  theme_ipsum() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                   size = 13, colour = "black"),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 10, vjust = 55, hjust = 0.007),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.9, 
                                    face = "bold", size = 13),
        plot.margin = unit(c(0.2,1,0,1), "cm"),
        legend.key.size = unit(0.4,"cm"),
        legend.title = element_text(size = 10)) +
  coord_fixed(ratio = 2)

## Merge the plots
ridge.g <- ggplotGrob(ridge.plot)
pearson.g <- ggplotGrob(pearson.plot)
g <- rbind(ridge.g, pearson.g, size = "last")
g$widths <- unit.pmax(ridge.g$widths, pearson.g$widths)
grid.newpage()
grid.draw(g)

## Save the figures
ggsave("./figure/ridge_heatmap.tiff", g, width = 8.65, height = 5)
ggsave("./figure/ridge_heatmap.png", g, width = 8.65, height = 5)






























