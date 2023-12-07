# Author: Nick Dibley and Teryn Morgan 

# Install neccessary packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("edgeR")
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
counts <- read.csv("Final Project/R Code/GSE97760_series_matrix.csv", row.names = 1)
counts <- data.matrix(counts, rownames.force = TRUE)

### LIMMA IMPLMENTATION ### ---------------------------------------------------
# Normalize counts for matrix
counts <- normalizeBetweenArrays(counts)
dge <- DGEList(counts=counts)

# Generate labels for series matrix
condition <- factor(rep(c("Control", "AD"), each = 10))
condition <- condition[-length(condition)]
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
GSE97760_DEGs <- topTreat(fit, coef=ncol(design), number = 42713)

# DEG indices
GSE97760_DEG_names <- row.names(GSE97760_DEGs)
saveRDS(GSE97760_DEGs_names, file = "GSE97760_DEGs.rds")

### LINK ENTREZ IDS ### --------------------------------
# Get soft file for GSE97760
GSE97760 <- getGEO('GSE97760', destdir=".")
GSE97760 <- getGEO(filename='GPL16699.soft.gz')
GSE97760_df <- Table(GSE97760)
GSE97760_df <- GSE97760_df %>% mutate_all(na_if,"")

# Filter soft file by indices in DEGs 
GSE97760_filtered <- subset(GSE97760_df, ID %in% GSE97760_DEG_names) # nrow = 42713

# Convert GB_ACC (Accession number) to Entrez ID
GSE97760_entrez_ids <- mapIds(org.Hs.eg.db, keys = GSE97760_filtered$GB_ACC,
                             column = "ENTREZID", keytype = "ACCNUM", multiVals = "list")
GSE97760_filtered$ENTREZ_ID <- GSE97760_entrez_ids

GSE97760_DEGs$ENTREZ_ID <- GSE97760_filtered$ENTREZ_ID
# Filtered out NA and NULL values in ENTREZ_ID (42713 -> 31632)
GSE97760_DEGs <- subset(GSE97760_DEGs, !is.na(ENTREZ_ID))
GSE97760_DEGs <- GSE97760_DEGs[!sapply(GSE97760_DEGs$ENTREZ_ID, is.null), ] # nrow =31632

# Collapse ENTREZ_ID lists
GSE97760_DEGs$ENTREZ_ID <- sapply(GSE97760_DEGs$ENTREZ_ID, function(x) paste(x, collapse = ','))

# Separate into upreg/downreg
GSE97760_downreg <- subset(GSE97760_DEGs, AveExpr < 0) # nrow = 15377
GSE97760_upreg <- subset(GSE97760_DEGs, AveExpr > 0) # nrow = 16255

### EXPORT DATA ### ---------------------------------------
#write.csv(GSE97760_DEGs, "Data/GSE97760_DEGs_EntrezIDs.csv", row.names=TRUE)
#write.table(GSE97760_DEGs$ENTREZ_ID, "Data/GSE97760_DEGs_EntrezIDs.txt", sep = ",", row.names = FALSE, quote = FALSE)

#write.csv(GSE97760_upreg, "Data/GSE97760_upreg_DEGs_EntrezIDs.csv", row.names=TRUE)
#write.table(GSE97760_upreg$ENTREZ_ID, "Data/GSE97760_upreg_DEGs_EntrezIDs.txt", sep = ",", row.names = FALSE, quote = FALSE)


