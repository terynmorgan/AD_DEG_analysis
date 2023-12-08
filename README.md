**Fall 2023 INFO-B 627 Advanced Seminar** <br />
Names: Teryn Morgan and Nick Dibley 

**Required Files:** <br />
GSE63060_series_matrix_noMCI.csv —> Series matrix for GEO Dataset GSE63060<br />
GSE4226_series_matrix.csv —> Series matrix for GEO Dataset GSE4226<br />
GSE63061_series_matrix_noMCI.csv —> Series matrix for GEO Dataset GSE63061<br />
GSE97760_series_matrix.csv —> Series matrix for GEO Dataset GSE97760<br />

GSE63060_genes.csv —> Csv file with GSE63060 metadata<br />
GSE63061_genes.csv —> Csv file with GSE63061 metadata<br />
venn_result_upreg_DEGs.txt —> Txt file with consensus DEGs generated from Draw Venn Diagram <br />

DEGs_GSE97760.R —> R script to generate DEGs for GSE97760<br />
DEGs_GSE63061.R —> R script to generate DEGs for GSE63061<br />
DEGs_GSE63060.R —> R script to generate DEGs for GSE63060<br />
DEGs_GSE4226.R —> R script to generate DEGs for GSE4226<br />
Gene_Set_Interaction.R —> R script to generate consensus DEGs from four GEO Datasets<br />
GO_KEGG_Enrichment.R —> R script to generate GO and KEGG Enrichment on top 100 consensus DEGs<br />

**Required Packages:** <br />
BiocManager, GEOquery, edgeR, limma, impute, my gene, org.Hs.eg.db, clusterProfiler<br />
cytoHubba Plugin in Cytoscape 

**Preprocessing:**
1. Run DEGs_GSE97760.R using GSE97760_series_matrix.csv 
2. Run DEGs_GSE63061.R using GSE63061_series_matrix_noMCI.csv and GSE63061_genes.csv 
3. Run DEGs_GSE63060.R using GSE63060_series_matrix_noMCI.csv and GSE63060_genes.csv 
4. Run DEGs_GSE4226.R using GSE4226_series_matrix.csv<br />
**Output Files:**
  GSE4226_DEGs_EntrezIDs.txt<br />
  GSE63061_DEGs_EntrezIDs.txt<br />
  GSE97760_DEGs_EntrezIDs.txt<br />
  GSE97760_DEGs_EntrezIDs.txt<br />

**Generate consensus DEGs Entrez IDs:**
1. Upload the following files to the Draw Venn Diagram website:
	GSE4226_DEGs_EntrezIDs.txt<br />
	GSE63061_DEGs_EntrezIDs.txt<br />
	GSE97760_DEGs_EntrezIDs.txt<br />
	GSE97760_DEGs_EntrezIDs.txt<br />
**Output Files:** <br />
  venn_result_upreg_DEGs.txt<br />
  Venn_Result_Upregulated_DEGs.png<br />
  Venn_Result_DEGs.png<br />

**Generate top 100 consensus DEGs:**
1. Run Gene_Set_Interaction.R using the following files: 
	GSE4226_DEGs_EntrezIDs.txt<br />
	GSE63061_DEGs_EntrezIDs.txt<br />
	GSE97760_DEGs_EntrezIDs.txt<br />
	GSE97760_DEGs_EntrezIDs.txt<br />
	venn_result_upreg_DEGs.txt<br />
**Output Files:** <br />
	Top100_upreg_DEG_EntrezIDs.txt<br />
	Top100_upreg_DEGs.csv<br />

**Generate STRING Interactions:**
1. Using Top100_upreg_DEGs.csv, upload to the STRING database to generate STRING interactions <br />
**Output File:** <br />
  string_interactions.tsv<br />
  string_normal_image.png

**Cytoscape Network Quantification:**
1. Upload string_interactions.tsv to Cytoscape
2. In the cytoHubba plugin, select the uploaded target network and compute node scores under the following metrics: Degree, MNC, Radially, Stress, and Closeness
3. Download the csv for the top 20 nodes in each metric network
4. Upload the resulting Entrez ID into Draw Venn diagram to generate the consensus hub gene Entrez IDs<br />
**Output Files:** <br />
  hub_genes_venn.txt<br />
  hub_genes_venn_result.ong

**Conduct GO and KEGG Enrichment Analysis:**
1. Run GO_KEGG_Enrichment.R using the following files: 
	Top100_upreg_DEG_EntrezIDs.txt<br />
	Top100_upreg_DEGs.csv<br />
	hub_genes_venn.txt<br />
2. Upload Top100_upreg_DEG_EntrezIDs.txt to Metascape to generate another KEGG enrichment analysis<br />
**Output Files:** <br />
  Enrichment_heatmap_HeatmapSelectedGO_20.pdf<br />
  Enrichment_heatmap_HeatmapSelectedGOTop100.pdf<br />
  GO_Enrichment_CC.png<br />
  GO_Enrichment_BP.png<br />
  GO_Enrichment_MF.png<br />
  KEGG_pathway_enrichment.png
