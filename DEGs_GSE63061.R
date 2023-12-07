# Author: Nick Dibley and Teryn Morgan 

# Install neccessary packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")
BiocManager::install("edgeR")
BiocManager::install("limma")
BiocManager::install("impute")
BiocManager::install("mygene")
BiocManager::install("org.Hs.eg.db")
library(readr)
library(dplyr)
library(GEOquery)
library(mygene)
library(org.Hs.eg.db) 
library(limma)
library(edgeR)
library(impute)

### IMPORT DATA ### -----------------------------------------------------------
counts <- read.csv("Data/GSE63061_series_matrix_noMCI.csv", row.names = 1)
counts <- data.matrix(counts, rownames.force = TRUE)

### LIMMA IMPLMENTATION ### ---------------------------------------------------
# Normalize counts for matrix
counts <- normalizeBetweenArrays(counts)
dge <- DGEList(counts=counts)

# Generate labels for series matrix
condition_ad <- factor(rep("AD", 139))
condition_control <- factor(rep("Control", 134))
condition <- c(condition_control, condition_ad)
design <- model.matrix(~0 + condition)

# Filter DGE list by condition labels
keep <- filterByExpr(dge, design)
dge <- dge[keep,,keep.lib.sizes=FALSE]

# Normalize dge
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log=TRUE, prior.count=3)

# Fit lm model 
fit <- lmFit(logCPM, design)
fit <- treat(fit, lfc=log2(1.2), trend=TRUE)
# Extract DEGs
GSE63061_DEGs <- topTreat(fit, coef=ncol(design), number = 725)

# DEG indices
GSE63061_DEGs_names <- row.names(GSE63061_DEGs)
saveRDS(GSE63061_DEGs_names, file = "GSE63061_DEGs.rds")

### LINK ENTREZ IDS ### --------------------------------
GSE63061_df <- read.csv("Final Project/Stuff for converting DEGs to genes/GSE63061_genes.csv", na.strings=c("","NA"))
GSE63061_filtered <- subset(GSE63061_df, ID %in% GSE63061_DEGs_names) # nrow = 725 

# Add Entrez_Gene_ID to GSE63060_DEGs
GSE63061_DEGs$ENTREZ_ID <- GSE63061_filtered$Entrez_Gene_ID
# Remove rows with NA
GSE63061_DEGs <- subset(GSE63061_DEGs, !is.na(ENTREZ_ID)) # nrow = 720

# Separate into upreg/downreg
GSE63061_downreg <- subset(GSE63061_DEGs, AveExpr < 0) 
GSE63061_upreg <- subset(GSE63061_DEGs, AveExpr > 0) # nrow = 720

### EXPORT FILES ### --------------------------------
#write.csv(GSE63061_DEGs, "Data/GSE63061_DEGs_EntrezIDs.csv", row.names=TRUE)
#write.table(GSE63061_DEGs$ENTREZ_ID, "Data/GSE63061_DEGs_EntrezIDs.txt", sep=",", row.names=FALSE)

