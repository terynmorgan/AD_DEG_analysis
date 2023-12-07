# Author: Teryn Morgan

# Install neccessary packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")
BiocManager::install("mygene")
BiocManager::install("org.Hs.eg.db")
library(readr)
library(dplyr)
library(GEOquery)
library(mygene)
library(org.Hs.eg.db) 

### IMPORT DATA ### -----------------------------------------------------------
GSE97760 <- read.csv("Data/GSE97760_upreg_DEGs_EntrezIDs.csv") # nrow = 16255
GSE4226 <- read.csv("Data/GSE4226_DEGs_EntrezIDs.csv") # nrow = 373
GSE63060 <- read.csv("Data/GSE63060_DEGs_EntrezIDs.csv") # nrow = 1296
GSE63061 <- read.csv("Data/GSE63061_DEGs_EntrezIDs.csv") # nrow = 720

shared_genes <- read.table("Data/venn_result_upreg_DEGs.txt") # nrow = 1004

### SUBSET EXPRESSION BY SHARED GENES ### -------------------------------------
GSE97760_filtered <- subset(GSE97760, ENTREZ_ID %in% shared_genes$V1)
GSE4226_filtered <- subset(GSE4226, ENTREZ_ID %in% shared_genes$V1)
GSE63060_filtered <- subset(GSE63060, ENTREZ_ID %in% shared_genes$V1)
GSE63061_filtered <- subset(GSE63061, ENTREZ_ID %in% shared_genes$V1)

### GENERATE TOP 100 UP/DOWN REGULATION ### ------------------------------------
# Merge filtered DEG results on ENTREZ_ID
DEG_comb <- rbind(GSE97760_filtered, GSE4226_filtered, GSE63060_filtered, GSE63061_filtered)
DEG_comb <- DEG_comb[order(DEG_comb$AveExpr, rev(DEG_comb$AveExpr), decreasing = TRUE), ]

# Get top 100 upregulated genes
top_100_upreg <- head(DEG_comb, 130) # Top 100 had 6 duplicate ENTREZ_IDs
top_100_upreg <- top_100_upreg[!duplicated(top_100_upreg$ENTREZ_ID),]

### EXPORT TOP 100 ENTREZ_IDS ### ---------------------------------------------
#write.table(top_100_upreg$ENTREZ_ID, "Top100_upreg_DEG_EntrezIDs.txt", sep = ",", row.names = FALSE, quote = FALSE)
#write.csv(top_100_upreg, "Top100_upreg_DEGs.csv", row.names = FALSE)