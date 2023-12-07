# AD_DEG_analysis
AIM: identify and intersect differentially expressed genes (DEGs) from four peripheral blood microarray datasets between Alzheimer’s Disease (AD) patients and controls

Names: Teryn Morgan and Nick Dibley 

**Required Files:**
GSE63060_series_matrix_noMCI.csv —> Series matrix for GEO Dataset GSE63060
GSE4226_series_matrix.csv —> Series matrix for GEO Dataset GSE4226
GSE63061_series_matrix_noMCI.csv —> Series matrix for GEO Dataset GSE63061
GSE97760_series_matrix.csv —> Series matrix for GEO Dataset GSE97760

GSE63060_genes.csv —> Csv file with GSE63060 metadata
GSE63061_genes.csv —> Csv file with GSE63061 metadata
venn_result_upreg_DEGs.txt —> Txt file with consensus DEGs generated from Draw Venn Diagram 

DEGs_GSE97760.R —> R script to generate DEGs for GSE97760
DEGs_GSE63061.R —> R script to generate DEGs for GSE63061
DEGs_GSE63060.R —> R script to generate DEGs for GSE63060
DEGs_GSE4226.R —> R script to generate DEGs for GSE4226
Gene_Set_Interaction.R —> R script to generate consensus DEGs from four GEO Datasets
GO_KEGG_Enrichment.R —> R script to generate GO and KEGG Enrichment on top 100 consensus DEGs

**Required Packages:**
BiocManager, GEOquery, edgeR, limma, impute, my gene, org.Hs.eg.db, clusterProfiler
cytoHubba Plugin in Cytoscape 

**Preprocessing:**
1. Run DEGs_GSE97760.R using GSE97760_series_matrix.csv 
2. Run DEGs_GSE63061.R using GSE63061_series_matrix_noMCI.csv and GSE63061_genes.csv 
3. Run DEGs_GSE63060.R using GSE63060_series_matrix_noMCI.csv and GSE63060_genes.csv 
4. Run DEGs_GSE4226.R using GSE4226_series_matrix.csv 
**Output Files: **
  GSE4226_DEGs_EntrezIDs.txt
  GSE63061_DEGs_EntrezIDs.txt
  GSE97760_DEGs_EntrezIDs.txt
  GSE97760_DEGs_EntrezIDs.txt

**Generate consensus DEGs Entrez IDs:**
1. Upload the following files to the Draw Venn Diagram website:
	GSE4226_DEGs_EntrezIDs.txt
	GSE63061_DEGs_EntrezIDs.txt
	GSE97760_DEGs_EntrezIDs.txt
	GSE97760_DEGs_EntrezIDs.txt
**Output Files:**
  venn_result_upreg_DEGs.txt
  Venn_Result_Upregulated_DEGs.png
  Venn_Result_DEGs.png

**Generate top 100 consensus DEGs:**
1. Run Gene_Set_Interaction.R using the following files: 
	GSE4226_DEGs_EntrezIDs.txt
	GSE63061_DEGs_EntrezIDs.txt
	GSE97760_DEGs_EntrezIDs.txt
	GSE97760_DEGs_EntrezIDs.txt
	venn_result_upreg_DEGs.txt 
**Output Files: **
	Top100_upreg_DEG_EntrezIDs.txt
	Top100_upreg_DEGs.csv

**Generate STRING Interactions:**
1. Using Top100_upreg_DEGs.csv, upload to the STRING database to generate STRING interactions 
**Output File:**
  string_interactions.tsv
  string_normal_image.png


**Cytoscape Network Quantification:**
1. Upload string_interactions.tsv to Cytoscape
2. In the cytoHubba plugin, select the uploaded target network and compute node scores under the following metrics: Degree, MNC, Radially, Stress, and Closeness
3. Download the csv for the top 20 nodes in each metric network
4. Upload the resulting Entrez ID into Draw Venn diagram to generate the consensus hub gene Entrez IDs
**Output Files:**
	hub_genes_venn.txt
  hub_genes_venn_result.ong

**Conduct GO and KEGG Enrichment Analysis:**
1. Run GO_KEGG_Enrichment.R using the following files: 
	Top100_upreg_DEG_EntrezIDs.txt
	Top100_upreg_DEGs.csv
	hub_genes_venn.txt
2. Upload Top100_upreg_DEG_EntrezIDs.txt to Metascape to generate another KEGG enrichment analysis
**Output Files:**
	Enrichment_heatmap_HeatmapSelectedGO_20.pdf
  Enrichment_heatmap_HeatmapSelectedGOTop100.pdf
  GO_Enrichment_CC.png
  GO_Enrichment_BP.png
  GO_Enrichment_MF.png
  KEGG_pathway_enrichment.png
