# Author: Teryn Morgan

# Install necessary packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("mygene")
library(org.Hs.eg.db) 
library(clusterProfiler)
library(mygene)

### IMPORT DATA ### -----------------------------------------------------------
top_100_upreg_entrez <- read.csv("Top100_upreg_DEG_EntrezIDs.txt", header=FALSE) 
top_100_upreg_DEGs <- read.csv("Top100_upreg_DEGs.csv") 

hub_genes <- read.table("hub_genes_venn.txt")$V1

### GENERATE RANKED LISTS ### -------------------------------------------------
top_100_upreg_ranked <- top_100_upreg_DEGs[,2]
names(top_100_upreg_ranked) <- as.character(top_100_upreg_DEGs[,7])
top_100_upreg_ranked <- sort(top_100_upreg_ranked, decreasing = TRUE)

# Convert Gene Symbols to Entrez ID
hub_entrez_ids <- mapIds(org.Hs.eg.db, keys = hub_genes,
                             column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
# Extract DEGs for hub genes
hub_DEGs <- subset(top_100_upreg_DEGs, ENTREZ_ID %in% hub_entrez_ids)

# Create ranked list for GSEA
hub_DEGs_ranked <- hub_DEGs[,2]
names(hub_DEGs_ranked) <- as.character(hub_DEGs[,7])
hub_DEGs_ranked <- sort(hub_DEGs_ranked, decreasing = TRUE)

### GO ENRICHMENT ### ---------------------------------------------------------
# Generate 
go_enrichment <- function(ontology_term){
  # Returns enrichment GO categories after FDR control (overrepresented)
  enrich_upreg <- enrichGO(
    gene = top_100_upreg_entrez$V1,
    OrgDb = org.Hs.eg.db,  
    ont = ontology_term,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  return(enrich_upreg)
}

# Generate dot plots for GO enrichments
dotplot(go_enrichment("BP"), color = "p.adjust", showCategory = 10, title ="Upregulated GO Enrichment: Biological Process")
dotplot(go_enrichment("MF"), color = "p.adjust", showCategory = 10, title ="Upregulated GO Enrichment: Molecular Function")
dotplot(go_enrichment("CC"), color = "p.adjust", showCategory = 10, title ="Upregulated GO Enrichment: Cellular Component")

### KEGG ENRICHMENT ### -------------------------------------------------------
# Generates KEGG Enrichment Analysis of a gene set
kegg_enrich <- enrichKEGG(
  gene = top_100_upreg_entrez$x, 
  organism = "hsa", 
  keyType="kegg",
  pvalueCutoff=0.05,
  minGSSize=10,
  maxGSSize=500,
  qvalueCutoff=0.05
)
head(kegg_enrich)
# Opens browser to visualize KEGG pathways
browseKEGG(kegg_enrich, 'hsa04110')

### GSEA HUB GENES ### ------------------------------------------------------
# Generates GSEA of Gene Ontology Enrichment 
GO_enrich_upreg <- gseGO(
  geneList = top_100_upreg_ranked, 
  ont="ALL",
  OrgDb=org.Hs.eg.db,
  keyType="ENTREZID",
  minGSSize=10,
  maxGSSize=500,
  pvalueCutoff=1.0,
  verbose=FALSE,
  scoreType = "pos"
)
head(GO_enrich_upreg) # No enriched genes
dotplot(GO_enrich_upreg, color = "p.adjust", showCategory = 10, title="Upregulated DEG Gene Set Enrichment Analysis of Gene Ontology") 