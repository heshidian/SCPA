# SCPA
A R package for single cell RNA-seq data analysis created by Ruiqing Zhai and Shidian He

### Operation manual

First set up the working directory

````R
setwd("./SCPA/test/pbmc5k")
````

 Load the SCPA 

````R
library(SCPA)
````

 Creating the result file 

````R
SCPA::prepare_1(run)
````

Read the single-cell data file and build the Seurat object

````R
sce <- Seurat::Read10X(data.dir = "./pbmc5k/")
sce <- Seurat::CreateSeuratObject(counts = sce, project = "pbmc5k")
````

 Data quality statistics 

````R
sce <- SCPA::data_quality_2(sce, group = "orig.ident", dir = "./1.Data quality")
````

 Data quality control 

````R
sce <- SCPA::quality_control_3(sce,
                               group = "orig.ident",
                               nFeature_lower = 300,
                               nFeature_upper = 8000,
                               nCount_lower = 3,
                               nCount_upper = 10000,
                               percent_mito_lower = 0,
                               percent_mito_upper = 25,
                               percent_hb_upper = 10,
                               percent_hb_lower = 0,
                               percent_ribo_lower = 3,
                               percent_ribo_upper = 1000,
                               dir = "./2.Quality control")
````

Single cell standard process analysis

````R
sce <- SCPA::standard_flow_4(sce,
                             group = "orig.ident",
                             dim_number = 20,
                             dbs_rate = NULL,
                             remove_doublets = TRUE,
                             nmf_rank_number = 10,
                             Scale_cellcycle = FALSE,
                             resolution_low = 0.1,
                             resolution_up = 1,
                             resolution_by = 0.1,
                             dir = "./3.scRNAflow")
````

Marker gene expression analysis and subgroup difference analysis

````R
sce <- SCPA::find_markers_5(sce,
                            group = "seurat_clusters",
                            genes_list = c("IL7R", "CD3E", "CD8A", "CD3D", "CD2"),
                            dir = "./4.Find markers")
````

Cell  annotation

````R
sce <- SCPA::cell_markers_6(sce,method = "SingleR")
#or
sce <- SCPA::cell_markers_6(sce,method = "personalization",
                            group = "seurat_clusters",
                            cells_markers = c("PTPRC", "CD3D", "CD3E", "CD4", "CD8A", "CD19", "CD79A", "MS4A1","IGHG1", "MZB1", "SDC1", "CD68", "CD163", "CD14", "LAMP3", "IDO1", "IDO2", "CD1E",
"CD1C", "KLRB1", "NCR1", "FGF7", "MME", "ACTA2", "PECAM1", "VWF"),
                            dir = "./5.Cell markers")
sce<-  SCPA::cell_markers_6_personalization(sce,
                                            group = "seurat_clusters",
                                            ClusterID = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16"),
                                            celltype = c("T cells", "T cells", "T cells", "B cells", "T cells", "Monocyte",
                                                         "T cells", "T cells", "B cells", "NK cells","T cells", "T cells",
                                                         "T cells", "T cells", "Monocyte","T cells", "T cells"),
                                            dir = "./5.Cell markers")
````

Cell differentiation analysis and enrichment analysis

````R
sce <- SCPA::cell_genes_7(sce,
                          group_cell_1 = "B cells",
                          group_cell_2 = "T cells",
                          group = "annotation",
                          Two_cell_enrichment = TRUE,
                          dir = "./6.Cell genes"
)
````

 Save the result file 

````R
save(sce,file = "result.RData")
````

