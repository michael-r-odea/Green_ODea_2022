library(dplyr)
library(Seurat)
library(patchwork)

zhong.unfiltered <- read.table("GSE104276_all_pfc_2394_UMI_count_NOERCC.xls", header = TRUE , sep= "\t") # read in count matrix

zhong.unfil <- CreateSeuratObject(counts = zhong.unfiltered, project = "zhong") # create seurat object
zhong.unfil

zhong.unfil[["percent.mt"]] <- PercentageFeatureSet(zhong.unfil, pattern = "^MT.") # calculate percent of features which map to mitochondrial genes (common QC metric)

VlnPlot(zhong.unfil, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, split.by = NULL) # produce QC plots for determining filtering parameters
FeatureScatter(zhong.unfil, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(zhong.unfil, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

zhong <- subset(zhong.unfil, subset = nFeature_RNA > 1000) # remove cells with less than 1000 unique expressed features
zhong

zhong <- NormalizeData(zhong, normalization.method = "LogNormalize", scale.factor = 10000) # log-normalize the data

zhong <- FindVariableFeatures(zhong, selection.method = "vst", nfeatures = 2000) # find variable features

all.genes <- rownames(zhong)
zhong <- ScaleData(zhong, features = all.genes) # scale gene expression

#Linear dimensionality reduction
zhong <- RunPCA(zhong, features = VariableFeatures(object = zhong), npcs = 100) # run PCA
ElbowPlot(zhong, ndims = 100) # create elbow plots to determine number of PCs to use for clustering
ElbowPlot(zhong, ndims = 50) # 13 PCAs looks about right

zhong <- FindNeighbors(zhong, dims = 1:13) # cluster cells with a number of different resolution values

zhong <- FindClusters(object = zhong, resolution = 0.1)
zhong <- FindClusters(object = zhong, resolution = 0.2)
zhong <- FindClusters(object = zhong, resolution = 0.3)
zhong <- FindClusters(object = zhong, resolution = 0.4)
zhong <- FindClusters(object = zhong, resolution = 0.5)
zhong <- FindClusters(object = zhong, resolution = 0.6)
zhong <- FindClusters(object = zhong, resolution = 0.7)
zhong <- FindClusters(object = zhong, resolution = 0.8)
zhong <- FindClusters(object = zhong, resolution = 0.9)
zhong <- FindClusters(object = zhong, resolution = 1.0)
zhong <- FindClusters(object = zhong, resolution = 1.1)
zhong <- FindClusters(object = zhong, resolution = 1.2)
zhong <- FindClusters(object = zhong, resolution = 1.3)
zhong <- FindClusters(object = zhong, resolution = 1.4)
zhong <- FindClusters(object = zhong, resolution = 1.5)
zhong <- FindClusters(object = zhong, resolution = 1.6)
zhong <- FindClusters(object = zhong, resolution = 1.7)
zhong <- FindClusters(object = zhong, resolution = 1.8)
zhong <- FindClusters(object = zhong, resolution = 1.9)
zhong <- FindClusters(object = zhong, resolution = 1.9)
zhong <- FindClusters(object = zhong, resolution = 2.0)

# cluster tree analysis
library(clustree)

clustree(zhong) # create cluster tree plot to help determine proper resolution to use for clustering; using 0.2

# run find clusters again with the chosen resolution
zhong <- FindClusters(object = zhong, resolution = 0.2)

zhong <- RunUMAP(zhong, dims = 1:13) # run UMAP
DimPlot(zhong, reduction = "umap") # create Umap plot (Figure 3F)

zhong.markers <- FindAllMarkers(zhong, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) # find DEGS that define the 7 clusters
zhong.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC) # cluster 5 appears to be microglia

# generate feature plots in Figure 3F (plots were cropped in figure image to emphasize microglia cluster)
FeaturePlot(zhong, features = "TMEM119", order = TRUE)
FeaturePlot(zhong, features = "P2RY12", order = TRUE)
FeaturePlot(zhong, features = "LCP1", order = TRUE)
FeaturePlot(zhong, features = "MRC1", order = TRUE)

microglia.features <- c("TMEM119", "P2RY12", "LCP1", "CX3CR1", "IRF8", "MPEG1", "MRC1", "LYVE1") # microglia markers to test expression (plus  MRC1 and LYVE-1, a lymphatic marker)

DoHeatmap(zhong, features = microglia.features) + NoLegend() # generate heatmap of microglia features expressed across all clusters (not shown in manuscript); confirms cluster 5 is microglia

microglia <- subset(zhong, idents = 5) # subset out microglia
microglia

lymph_microglia <- subset(microglia, subset = MRC1 >= 1.0) # subsets out Mrc1+ microglia
lymph_microglia

other_microglia <- subset(microglia, subset = MRC1 < 1.0) # subsets out Mrc1- microglia
other_microglia

Idents(lymph_microglia) <- 0 # assigns identities to MRc1+ and mrc1- microglia
Idents(other_microglia) <- 1

microglia_combined <- merge(lymph_microglia, y = other_microglia, add.cell.ids = c("lymph", "other"), project = "zhong") # merges Mrc1+ and Mrc1- subsetted datasets together to give one combined microglia dataset with identities based on Mrc1 expression

DoHeatmap(microglia_combined, features = microglia.features) # generate heatmap of microglia marker gene expression between Mrc1+ and Mrc1- microglia (Figure 3H)

microglia_markers <- FindMarkers(microglia_combined, ident.1 = 0, ident.2 = 1, features = microglia.features, logfc.threshold = 0, min.pct = 0) # run DE testing on canonical microglia markers

microglia_markers$p_val <- NULL

library(reactable)

library(htmltools)


## generates DE expression table found in Figure S1D
# Render a bar chart with positive and negative values
bar_chart_pos_neg <- function(label, value, max_value = 2.00, height = "16px",
                              pos_fill = "#02aa5c", neg_fill = "#ff121a") {
  neg_chart <- div(style = list(flex = "1 1 0"))
  pos_chart <- div(style = list(flex = "1 1 0"))
  width <- paste0(abs(value / max_value) * 100, "%")
  
  if (value < 0) {
    bar <- div(style = list(marginLeft = "8px", background = neg_fill, width = width, height = height))
    chart <- div(style = list(display = "flex", alignItems = "center", justifyContent = "flex-end"), label, bar)
    neg_chart <- tagAppendChild(neg_chart, chart)
  } else {
    bar <- div(style = list(marginRight = "8px", background = pos_fill, width = width, height = height))
    chart <- div(style = list(display = "flex", alignItems = "center"), bar, label)
    pos_chart <- tagAppendChild(pos_chart, chart)
  }
  
  div(style = list(display = "flex"), neg_chart, pos_chart)
}


reactable(microglia_markers, columns = list(avg_logFC = colDef(name = "Log Fold Change", format = colFormat(digits = 4), 
                                                               cell = function(value) {
                                                                 label <- round(value, digits = 2)
                                                                 bar_chart_pos_neg(label, value)
                                                               },
                                                               align = "center",
                                                               minWidth = 150), 
                                            pct.1 = colDef(name = "Lymphatic microglia with expression", format = colFormat(percent = TRUE, digits = 1)), 
                                            pct.2 = colDef(name = "Other microglia with expression", format = colFormat(percent = TRUE, digits = 1)), 
                                            p_val_adj = colDef(name = "Mann-Whitney U Test", cell = function(value) {
                                              if (value == 1) {
                                                paste0("p = ", value)
                                              } else {
                                                paste0("p < ", signif(value, digits = 2))
                                              }}, format = colFormat(digits = 5), 
                                              style = function(value) {
                                                if (value <= 0.05) {
                                                  color <- "#008000"
                                                } else {
                                                  color <- "#e00000"
                                                }
                                                list(color = color, fontWeight = "bold")
                                              })), 
          compact = TRUE,
          resizable = TRUE, 
          bordered = TRUE, 
          fullWidth = FALSE, 
          rowStyle = function(index){
            if(index == 1) list(background = "rgba(0, 0, 0, 0.05)")
          })


#export zhong dataset to anndata format for scanpy in python; will be used in "vlnplots.py" to generate stacked violin plot in Figure 3G
microglia_combined$cell_ident <- Idents(microglia_combined)
sceasy::convertFormat(microglia_combined, from="seurat", to="anndata",outFile='zhong_microglia.h5ad')


# subclutering analysis on microglia cluster to determine if Mrc1+ cells are broadly similar to Mrc1- microglia
microglia_combined <- FindVariableFeatures(microglia_combined, selection.method = "vst", nfeatures = 2000) # find variable features
microglia_combined <- RunPCA(microglia_combined, features = VariableFeatures(object = microglia_combined)) # run PCA

ElbowPlot(microglia_combined) # generate elbow plot to help choose number of PCs to use for clustering; 11 looks good
microglia.subclustering <- FindNeighbors(microglia_combined, dims = 1:11)

microglia.subclustering <- FindClusters(object = microglia.subclustering, resolution = 0.1) # cluster cells using different resolution values
microglia.subclustering <- FindClusters(object = microglia.subclustering, resolution = 0.2)
microglia.subclustering <- FindClusters(object = microglia.subclustering, resolution = 0.3)
microglia.subclustering <- FindClusters(object = microglia.subclustering, resolution = 0.4)
microglia.subclustering <- FindClusters(object = microglia.subclustering, resolution = 0.5)
microglia.subclustering <- FindClusters(object = microglia.subclustering, resolution = 0.6)
microglia.subclustering <- FindClusters(object = microglia.subclustering, resolution = 0.7)
microglia.subclustering <- FindClusters(object = microglia.subclustering, resolution = 0.8)
microglia.subclustering <- FindClusters(object = microglia.subclustering, resolution = 0.9)
microglia.subclustering <- FindClusters(object = microglia.subclustering, resolution = 1.0)
microglia.subclustering <- FindClusters(object = microglia.subclustering, resolution = 1.1)
microglia.subclustering <- FindClusters(object = microglia.subclustering, resolution = 1.2)
microglia.subclustering <- FindClusters(object = microglia.subclustering, resolution = 1.3)
microglia.subclustering <- FindClusters(object = microglia.subclustering, resolution = 1.4)
microglia.subclustering <- FindClusters(object = microglia.subclustering, resolution = 1.5)
microglia.subclustering <- FindClusters(object = microglia.subclustering, resolution = 1.6)
microglia.subclustering <- FindClusters(object = microglia.subclustering, resolution = 1.7)
microglia.subclustering <- FindClusters(object = microglia.subclustering, resolution = 1.8)
microglia.subclustering <- FindClusters(object = microglia.subclustering, resolution = 1.9)
microglia.subclustering <- FindClusters(object = microglia.subclustering, resolution = 1.9)
microglia.subclustering <- FindClusters(object = microglia.subclustering, resolution = 2.0)


clustree(microglia.subclustering) # generate clutree plot to help determine resolution to use 

# run find clusters again with the resolution chosen; 0.4
microglia.subclustering <- FindClusters(object = microglia.subclustering, resolution = 0.4)
head(microglia.subclustering[[]])

microglia.subclustering <- RunUMAP(microglia.subclustering, dims = 1:11) # run umap
DimPlot(microglia.subclustering, reduction = "umap") # plot umap (plot not included in manuscript)

DoHeatmap(microglia.subclustering, features = microglia.features) # generate heatmap to determine if canonical microglia gene expression differs between two subclusters; Figure S1C
# heatmap shows cluster 0 exhibts lower expression of microglia markers than cluster 1 (perhaps it is an immature population); however, MRC1+ cells are found in both clusters




