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
counts <- read.csv("R Code/GSE63060_series_matrix_noMCI.csv", row.names = 1)
counts <- data.matrix(counts, rownames.force = TRUE)

### LIMMA IMPLMENTATION ### ---------------------------------------------------
# Normalize counts for matrix
counts <- normalizeBetweenArrays(counts)
dge <- DGEList(counts=counts)

# Generate labels for series matrix
condition_ad <- factor(rep("AD", 145))
condition_control <- factor(rep("Control", 104))
condition <- c(condition_ad, condition_control)
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
GSE63060_DEGs <- topTreat(fit, coef=ncol(design), number = 2057)

# DEG indices
GSE63060_DEGs_names <- row.names(GSE63060_DEGs)
saveRDS(DEGs, file = "GSE63060_DEGs.rds")

### LINK ENTREZ IDS ### --------------------------------
GSE63060_df <- read.csv("Data/GSE63060_genes.csv", na.strings=c("","NA"))

# Filters GSE63060_df by indices in DEGs 
GSE63060_filtered <- subset(GSE63060_df, ID %in% GSE63060_DEGs_names) # nrow = 1709

# Add Entrez_Gene_ID to GSE63060_DEGs
GSE63060_DEGs$ENTREZ_ID <- GSE63060_filtered$Entrez_Gene_ID
# Remove rows with NA
GSE63060_DEGs <- subset(GSE63060_DEGs, !is.na(ENTREZ_ID)) # nrow = 1296

# Separate into upreg/downreg
GSE63060_downreg <- subset(GSE63060_DEGs, AveExpr < 0) # nrow = 0
GSE63060_upreg <- subset(GSE63060_DEGs, AveExpr > 0) # nrow = 1296

### EXPORT FILES ### --------------------------------
#write.csv(GSE63060_DEGs, "Data/GSE63060_DEGs_EntrezIDs.csv", row.names=TRUE)
#write.table(GSE63060_DEGs$ENTREZ_ID, "Data/GSE63060_DEGs_EntrezIDs.txt", sep=",", row.names=FALSE)
