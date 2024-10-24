setwd("/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC/")

library("tidyverse")
library("glmnet")
library("ridge")

set.seed(12345)

# ===== Import data ============================================================
binary.mtx <- read.csv("./data/FIMO/GSE124228_binary_IDR.csv", header = TRUE, row.names = 1)
DA.res <- read.csv("./data/DA analysis/GSE124228_DA_IDR.csv", header = TRUE, row.names = 1)

# ===== Preprocessing ==========================================================
## Remove NA value
column <- c("log2FoldChange", "Peak", "Chr", "Start", "End")
DA.res <- DA.res %>% select(all_of(column))
na <- DA.res[apply(DA.res, 1, function(x) any(is.na(x))),]
na$peak <- paste0(na$Chr, "_", na$Start, "_", na$End)
DA.res$Start <- DA.res$Start - 1
DA.res$peak <- paste0(DA.res$Chr, "_", DA.res$Start, "_", DA.res$End)
na.list <- na$peak
DA.res <- DA.res[!(DA.res$peak %in% na.list),]

## Obtain log2FoldChange
DA.res <- DA.res[order(DA.res$peak),]
DA.res <- DA.res$log2FoldChange

## Preprocess the binary.matrix ##
binary.mtx <- binary.mtx[!(rownames(binary.mtx) %in% na.list),]
binary.mtx <- binary.mtx[order(rownames(binary.mtx)),]
binary.mtx <- as.matrix(binary.mtx)

# ===== Find optimal lambda ====================================================
## Use five fold cross validation
fit.ridge <- glmnet(binary.mtx, DA.res, alpha = 0)
cv.model <- cv.glmnet(binary.mtx, DA.res, alpha = 0, family = "gaussian", nfolds = 5)
lambda <- cv.model$lambda.min

# ===== Fit to the ridge regression model ======================================
binary.mtx <- as.data.frame(binary.mtx)
model <- linearRidge(DA.res ~ ., data = binary.mtx, lambda = lambda, scaling = "none")
summary <- summary(model)
summary_df <- as.data.frame(summary$summaries$summary1$coefficients)

# ===== Multiple hypothesis correction with BH =================================
colnames(summary_df) <- c("TF", "cor", "SE", "t_value", "p_value")
summary_df$FDR <- p.adjust(summary_df$p_value, method = "BH")








