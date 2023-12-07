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
counts <- read.csv("Data/GSE4226_series_matrix.csv", row.names = 1)
counts <- data.matrix(counts, rownames.force = TRUE)

### LIMMA IMPLMENTATION ### ---------------------------------------------------
# Normalize counts for matrix
counts <- normalizeBetweenArrays(counts)
dge <- DGEList(counts=counts)

# Generate labels for series matrix
condition <- factor(rep(c("Control", "AD"), each = ncol(counts)/4, times = 2))
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
GSE4226_DEGs <- topTreat(fit, coef=ncol(design), number = 428)

# DEG indices 
GSE4226_DEG_names <- row.names(GSE4226_DEGs)
saveRDS(GSE4226_DEG_names, file = "GSE4226_DEGs.rds")


### LINK ENTREZ IDS ### -------------------------------------------------------
# Get soft file for GSE4226
GSE4226 <- getGEO('GSE4226', destdir=".")
GSE4226 <- getGEO(filename='GPL1211.soft.gz')
GSE4226_df <- Table(GSE4226)
GSE4226_df <- GSE4226_df %>% mutate_all(na_if,"")

# Filter soft file by indices in DEGs 
GSE4226_filtered <- subset(GSE4226_df, ID %in% GSE4226_DEG_names) # nrow = 428

# Convert GB_ACC (Accession number) to Entrez ID
GSE4226_entrez_ids <- mapIds(org.Hs.eg.db, keys = GSE4226_filtered$GB_ACC,
                             column = "ENTREZID", keytype = "ACCNUM", multiVals = "first")
GSE4226_filtered$ENTREZ_ID <- GSE4226_entrez_ids
GSE4226_DEGs$ENTREZ_ID <- GSE4226_filtered$ENTREZ_ID

# Filtered out NA and NULL values in ENTREZ_ID (428 -> 373)
GSE4226_DEGs <- subset(GSE4226_DEGs, !is.na(ENTREZ_ID))
GSE4226_DEGs <- GSE4226_DEGs[!sapply(GSE4226_DEGs$ENTREZ_ID, is.null), ] # nrow = 373

# Collapse ENTREZ_ID lists
GSE4226_DEGs$ENTREZ_ID <- sapply(GSE4226_DEGs$ENTREZ_ID, function(x) paste(x, collapse = ','))

### EXPORT DATA ### ---------------------------------------
#write.csv(GSE4226_DEGs, "Data/GSE4226_DEGs_EntrezIDs.csv", row.names=TRUE)
#write.table(GSE4226_DEGs$ENTREZ_ID, "Data/GSE4226_DEGs_EntrezIDs.txt", sep=",", quote=FALSE, row.names=FALSE)
