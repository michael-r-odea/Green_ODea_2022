# this script includes the code for the preprocessing and analysis of the Hammond et al (2018) single-cell RNA-sequencing data
# note: file paths are relative to local system and may not apply

# load required packages
library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(egg)
library(tibble)
library(HGC)
library(dendextend)
library(future)
library(ComplexHeatmap)
library(ggpubr)

# import DGE.txt count matrices (obtained from NCBI GEO Series GSE121654); this data includes all samples from E14 to P30; 
data.files <- list.files("./input_data_files/")
setwd("./input_data_files/")
hammond.data <- lapply(data.files, read.table, header = TRUE, row.names = 1, sep = "\t")
setwd("..")

# get sample names from input file names
names(hammond.data) <- stringr::str_replace(data.files, pattern = ".dge.txt", replacement = "")
hammond.names.1 <- names(hammond.data)
names(hammond.data) <- stringr::str_replace(hammond.names.1, pattern = "GSM......._", replacement = "")

sample_ids <- names(hammond.data)
timepoints <- gsub("_.*","",sample_ids) # get sample ages 

# create seurat objects from count matrices
hammond.data.s <- lapply(hammond.data, CreateSeuratObject, project = "Hammond")

# add timepoint and sample_id metadata to each seurat object in list
for (i in seq_along(hammond.data.s)){
  hammond.data.s[[i]]@meta.data[["timepoint"]] <- timepoints[i]
  hammond.data.s[[i]]@meta.data[["sample_id"]] <- sample_ids[i]
}

# Calculate the percentage of UMIs mapping to mitochondrial genes here and store it in percent.mt  
for (i in seq_along(hammond.data.s)){
  hammond.data.s[[i]]@meta.data[["percent.mt"]] <- PercentageFeatureSet(object = hammond.data.s[[i]], pattern = "^mt.")
}

# find highly variable genes
features <- SelectIntegrationFeatures(object.list = hammond.data.s, nfeatures = 3000)

# merge seurat objects
hammond <- merge(x = hammond.data.s[[1]], y = hammond.data.s[2:26])

#hammond$"percent.mt" <- PercentageFeatureSet(hammond, pattern = "^mt.")

# plot qc plots to determine what gene/UMI/and pct mito cutoffs to use in QC filtering: 
p1 <- VlnPlot(hammond, features = "nFeature_RNA", group.by = "sample_id", pt.size = 0) + labs(x = NULL) + 
  ggtitle("Number of genes detected") + NoLegend() + geom_hline(yintercept = 3000, linetype = 'dotted', col = 'red') +
  geom_hline(yintercept = 400, linetype = 'dotted', col = 'red')# cut off 400<x<3,000 genes seems reasonable
p2 <- VlnPlot(hammond, features = "nCount_RNA", group.by = "sample_id", pt.size = 0) + labs(x = NULL) + 
  ggtitle("Number of UMIs") + NoLegend() + geom_hline(yintercept = 10000, linetype = 'dotted', col = 'red') # cut off 10,000 genes seems reasonable
p3 <- VlnPlot(hammond, features = "percent.mt", group.by = "sample_id", pt.size = 0) + labs(x = NULL) + 
  ggtitle("Percent of UMIs mapping to mitochondrial genome") + NoLegend() + geom_hline(yintercept = 3, linetype = 'dotted', col = 'red') # cut off 3% seems reasonable

ggarrange(p1, p2, p3, ncol = 1)

tiff("qc_plots.tiff", res = 300, height = 12, width = 15, units = "in")
ggarrange(p1, p2, p3, ncol = 1)
dev.off()

postscript("qc_plots.ps", height = 12, width = 15)
ggarrange(p1, p2, p3, ncol = 1)
dev.off()


# based on QC plots, filtering out cells with <400 genes or > 3000 unique genes; cells with > 10,000 UMIs, and cells with greater than 3% of reads mapping to mitochondrial genes 
hammond <- subset(hammond, subset = nFeature_RNA > 400 & nFeature_RNA < 3000 & nCount_RNA < 10000 & percent.mt < 3)

# Log-normalize data
hammond <- NormalizeData(hammond) 

# scale all variable features, regressing out technical variables
hammond <- ScaleData(hammond, features = features, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt")) # regressing out UMI counts, genes detected counts, and percent.mt 
hammond <- RunPCA(hammond, features = features, npcs = 50, verbose = TRUE)

# integration with harmony; 
set.seed(59483) # process is stochastic; make sure to set a random seed so integration results are reproducible
hammond <- hammond %>% RunHarmony(c("sample_id"), plot_convergence = TRUE, kmeans_init_nstart=20, kmeans_init_iter_max=100) # integrate across sample_id

# plot Elbow plots before and after Harmony integration
tiff("harmony_elbowplot.tiff", units = "in", width = 5, height = 5, res = 300)
ElbowPlot(hammond, ndims = 50, reduction = "harmony")
dev.off()

tiff("uncorrected_pca_elbowplot.tiff", units = "in", width = 5, height = 5, res = 300) # 17 pc's perhaps?
ElbowPlot(hammond, ndims = 50, reduction = "pca") 
dev.off()

# visualize genes associated with each principal component
tiff("pcs_1_thru_4.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(hammond, dims = 1:4, reduction = "harmony")
dev.off()

tiff("pcs_5_thru_8.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(hammond, dims = 5:8, reduction = "harmony")
dev.off()

tiff("pcs_9_thru_12.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(hammond, dims = 9:12, reduction = "harmony")
dev.off()

tiff("pcs_13_thru_16.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(hammond, dims = 13:16, reduction = "harmony")
dev.off()

tiff("pcs_17_thru_20.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(hammond, dims = 17:20, reduction = "harmony")
dev.off()

tiff("pcs_21_thru_24.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(hammond, dims = 21:24, reduction = "harmony")
dev.off()

tiff("pcs_25_thru_28.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(hammond, dims = 25:28, reduction = "harmony")
dev.off()

tiff("pcs_29_thru_32.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(hammond, dims = 29:32, reduction = "harmony")
dev.off()

tiff("pcs_33_thru_36.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(hammond, dims = 33:36, reduction = "harmony")
dev.off()

tiff("pcs_37_thru_40.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(hammond, dims = 37:40, reduction = "harmony")
dev.off()

# based on the above results; run UMAP with 40 harmony-corrected principal components
hammond <- RunUMAP(hammond, reduction = "harmony", dims = 1:40)

# plot UMAP with sample metadata
DimPlot(hammond, group.by = "sample_id") + theme(aspect.ratio = 1)
DimPlot(hammond, group.by = "timepoint") + theme(aspect.ratio = 1)


# perform hierarchical clustering using HGC
hammond <- FindNeighbors(hammond, reduction = "harmony", dims = 1:40, k.param = 30)
hammond <- FindClusteringTree(hammond, graph.type = "SNN")
tree <- hammond@graphs$ClusteringTree
tree$height = log(tree$height + 1)

# cut tree
initial.k <- cutree(tree, k = 20, order_clusters_as_data = FALSE)
table(initial.k)

# 11 clusters were singlets; subset out subdendrograms and get dendrogram leaves for all other branches, excluding the singlets; 
subtree.list <- get_subdendrograms(as.dendrogram(tree), 12, order_clusters_as_data = FALSE)
dend.without.singlets <- subtree.list[[1]]

test = get_leaves_attr(dend.without.singlets, "label", simplify = FALSE)

# plot dendrogram (without singlets)
dend <- dend.without.singlets %>%
  dendextend::set("branches_k_color", k=5, value = viridis::viridis_pal()(5)) %>%
  dendextend::set("branches_lwd", 3) %>%
  dendextend::set("labels_colors", "white")

tiff("test_subdendrogram.tiff", res = 600, height = 15, width = 20, units = "in")
plot(dend)
dev.off()

# save dendrogram object
saveRDS(dend, "initial_clustering_dendrogram.rds")

# cut dendrogram at k = 11 
k11 <- cutree(dend, k = 11, order_clusters_as_data = FALSE)
table(k11)

# get cluster barcodes 
k11.df <- k11 %>% enframe(name = "barcode", value = "cluster") 
k11.c1 <- k11.df %>% filter(cluster == 1) %>% pull(barcode)
k11.c2 <- k11.df %>% filter(cluster == 2)%>% pull(barcode)
k11.c3 <- k11.df %>% filter(cluster == 3)%>% pull(barcode)
k11.c4 <- k11.df %>% filter(cluster == 4)%>% pull(barcode)
k11.c5 <- k11.df %>% filter(cluster == 5)%>% pull(barcode)
k11.c6 <- k11.df %>% filter(cluster == 6)%>% pull(barcode)
k11.c7 <- k11.df %>% filter(cluster == 7)%>% pull(barcode)
k11.c8 <- k11.df %>% filter(cluster == 8)%>% pull(barcode)
k11.c9 <- k11.df %>% filter(cluster == 9)%>% pull(barcode)
k11.c10 <- k11.df %>% filter(cluster == 10)%>% pull(barcode)
k11.c11 <- k11.df %>% filter(cluster == 11)%>% pull(barcode)

# overlay on UMAP
DimPlot(hammond, cells.highlight = k11.c1) + theme(aspect.ratio = 1) 
DimPlot(hammond, cells.highlight = k11.c2) + theme(aspect.ratio = 1)
DimPlot(hammond, cells.highlight = k11.c3) + theme(aspect.ratio = 1) 
DimPlot(hammond, cells.highlight = k11.c4) + theme(aspect.ratio = 1) 
DimPlot(hammond, cells.highlight = k11.c5) + theme(aspect.ratio = 1) 
DimPlot(hammond, cells.highlight = k11.c6) + theme(aspect.ratio = 1)
DimPlot(hammond, cells.highlight = k11.c7) + theme(aspect.ratio = 1) 
DimPlot(hammond, cells.highlight = k11.c8) + theme(aspect.ratio = 1) 
DimPlot(hammond, cells.highlight = k11.c9) + theme(aspect.ratio = 1) 
DimPlot(hammond, cells.highlight = k11.c9) + theme(aspect.ratio = 1) 
DimPlot(hammond, cells.highlight = k11.c9) + theme(aspect.ratio = 1) 

k11.c1.s <- hammond[,k11.c1]
k11.c1.s$k11_ident <- 1
k11.c2.s <- hammond[,k11.c2]
k11.c2.s$k11_ident <- 2
k11.c3.s <- hammond[,k11.c3]
k11.c3.s$k11_ident <- 3
k11.c4.s <- hammond[,k11.c4]
k11.c4.s$k11_ident <- 4
k11.c5.s <- hammond[,k11.c5]
k11.c5.s$k11_ident <- 5
k11.c6.s <- hammond[,k11.c6]
k11.c6.s$k11_ident <- 6
k11.c7.s <- hammond[,k11.c7]
k11.c7.s$k11_ident <- 7
k11.c8.s <- hammond[,k11.c8]
k11.c8.s$k11_ident <- 8
k11.c9.s <- hammond[,k11.c9]
k11.c9.s$k11_ident <- 9
k11.c10.s <- hammond[,k11.c10]
k11.c10.s$k11_ident <- 10
k11.c11.s <- hammond[,k11.c11]
k11.c11.s$k11_ident <- 11

k11.merged <- merge(k11.c1.s, c(k11.c2.s, k11.c3.s, k11.c4.s, k11.c5.s, k11.c6.s, k11.c7.s, k11.c8.s, k11.c9.s, k11.c10.s, k11.c11.s))
hammond$k11_ident <- k11.merged$k11_ident

rm(list = c("k11.c1.s", "k11.c2.s", "k11.c3.s", "k11.c4.s", "k11.c5.s", "k11.c6.s", "k11.c7.s", "k11.c8.s", "k11.c9.s", "k11.c10.s", "k11.c11.s", 
            "k11.merged"))
gc()

# one cell not assigned: singlet; removing
hammond <- hammond[,!is.na(hammond$k11_ident)] 

# plot clusters on UMAP
DimPlot(hammond, group.by = "k11_ident") + theme(aspect.ratio = 1)


initial.cluster.colors <- c('#e6194b', '#3cb44b', '#4363d8', '#f58231', '#911eb4', 
                            '#f032e6', '#800000', '#000075', '#fabebe', '#008080', '#e6beff')

# plot dendrogram
dend <- dend.without.singlets %>%
  dendextend::set("branches_k_color", k=11, value = initial.cluster.colors) %>%
  dendextend::set("branches_lwd", 3) %>%
  dendextend::set("labels_colors", "white")

dend <- dend %>%
  dendextend::set("branches_k_color", k=11, value = initial.cluster.colors) %>%
  dendextend::set("branches_lwd", 3) %>%
  dendextend::set("labels_colors", "white")

tiff("k11_initial_subdendrogram.tiff", res = 600, height = 15, width = 20, units = "in")
plot(dend)
dev.off()

# rename clusters
Idents(hammond) <- factor(Idents(hammond), levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"))
hammond <- RenameIdents(object = hammond, `11` = "A", `10` = "B", `9` = "C", `8` = "D", `7` = "E", `6` = "F", `3` = "G", 
                        `5` = "H", `4` = "I", 
                        `1` = "J", `2` = "K")

final.cluster.colors <- c('#e6beff', '#008080', '#fabebe', '#000075', '#800000', '#f032e6', '#4363d8', '#911eb4', '#f58231', '#e6194b', '#3cb44b')

# plot UMAP with cluster colors
tiff("initial_clusters_umap.tiff", res = 600, height = 6, width = 6, units = "in")
DimPlot(hammond) + theme(aspect.ratio = 1) + scale_color_manual(values = final.cluster.colors)
dev.off()

postscript("initial_clusters_umap.ps", height = 6, width = 6)
DimPlot(hammond) + theme(aspect.ratio = 1) + scale_color_manual(values = final.cluster.colors)
dev.off()

# run differential expression testing to find cluster markers
plan("multiprocess")
options(future.globals.maxSize= 10000*1024^2)

k11.markers <- FindAllMarkers(hammond)
k11.markers <- k11.markers %>% arrange(cluster, desc(avg_log2FC))

# save cluster markers data frame
saveRDS(k11.markers, "initial_clustering_markers.rds")
write.csv(k11.markers, "initial_clustering_markers.csv")


# plot stacked violin plot with cell type/cluster markers
tiff("initial_clusters_stackedviolinplot.tiff", res = 600, height = 10, width = 10, units = "in")
VlnPlot(hammond, features = c("Tmem119", "P2ry12", "Hexb", "F13a1", "Ccr1", "Ccr2", "Cldn5", "Vtn", "Pecam1", "Neurod6", "Nfib", "Elavl3"), 
        stack = TRUE, flip = TRUE, cols = rev(viridis::viridis_pal(option = "plasma")(12))) + labs(x = NULL, y = "Log-Normalized Expression") + 
  theme(axis.text.x.bottom = element_text(angle = 0, hjust = 0.5)) + NoLegend()
dev.off()

postscript("initial_clusters_stackedviolinplot.ps", height = 10, width = 10)
VlnPlot(hammond, features = c("Tmem119", "P2ry12", "Hexb", "F13a1", "Ccr1", "Ccr2", "Cldn5", "Vtn", "Pecam1", "Neurod6", "Nfib", "Elavl3"), 
        stack = TRUE, flip = TRUE, cols = rev(viridis::viridis_pal(option = "plasma")(12))) + labs(x = NULL, y = "Log-Normalized Expression") + 
  theme(axis.text.x.bottom = element_text(angle = 0, hjust = 0.5)) + NoLegend()
dev.off()

# save clustering dendrogram; 
saveRDS(dend, "initial_dendrogram.rds")
#save Seurat object
saveRDS(hammond, "hammond_initialclustering.rds")

# subcluster analysis: subset out 8 clusters which appear to be either microglia or macrophages on the basis of their marker gene expression
hammond.mgs <- subset(hammond, idents = c("A", "B", "C", "D", "E", "F", "G", "H"))

hammond.mgs 

# to perform subclustering, split object by sample, and find common variable features
hammond.mgs.list <- SplitObject(hammond.mgs, split.by = "sample_id")
features2 <- SelectIntegrationFeatures(object.list = hammond.mgs.list, nfeatures = 2000)

# scale variable features
hammond.mgs <- ScaleData(hammond.mgs, features = features2, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt")) # regressing out UMI counts, genes detected counts, and percent.mt 
# perform PCA
hammond.mgs <- RunPCA(hammond.mgs, features = features2, npcs = 50, verbose = TRUE)

# perform integration with Harmony
set.seed(6839) 
hammond.mgs <- hammond.mgs %>% RunHarmony(c("sample_id"), plot_convergence = TRUE, kmeans_init_nstart=20, kmeans_init_iter_max=1000)

# save Seurat object
saveRDS(hammond.mgs, "hammond_post_subcluster_harmony.rds")

# plot Elbow plots before and after integration
tiff("2harmony_elbowplot.tiff", units = "in", width = 5, height = 5, res = 300)
ElbowPlot(hammond.mgs, ndims = 50, reduction = "harmony")
dev.off()

tiff("2uncorrected_pca_elbowplot.tiff", units = "in", width = 5, height = 5, res = 300) # 17 pc's perhaps?
ElbowPlot(hammond.mgs, ndims = 50, reduction = "pca") 
dev.off()

# plot genes associated with each principal component
tiff("2pcs_1_thru_4.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(hammond.mgs, dims = 1:4, reduction = "harmony")
dev.off()

tiff("2pcs_5_thru_8.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(hammond.mgs, dims = 5:8, reduction = "harmony")
dev.off()

tiff("2pcs_9_thru_12.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(hammond.mgs, dims = 9:12, reduction = "harmony")
dev.off()

tiff("2pcs_13_thru_16.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(hammond.mgs, dims = 13:16, reduction = "harmony")
dev.off()

tiff("2pcs_17_thru_20.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(hammond.mgs, dims = 17:20, reduction = "harmony")
dev.off()

tiff("2pcs_21_thru_24.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(hammond.mgs, dims = 21:24, reduction = "harmony")
dev.off()

tiff("2pcs_25_thru_28.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(hammond.mgs, dims = 25:28, reduction = "harmony")
dev.off()

tiff("2pcs_29_thru_32.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(hammond.mgs, dims = 29:32, reduction = "harmony")
dev.off()

tiff("2pcs_33_thru_36.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(hammond.mgs, dims = 33:36, reduction = "harmony")
dev.off()

tiff("2pcs_37_thru_40.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(hammond.mgs, dims = 37:40, reduction = "harmony")
dev.off()

tiff("2pcs_41_thru_44.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(hammond.mgs, dims = 41:44, reduction = "harmony")
dev.off()

tiff("2pcs_45_thru_48.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(hammond.mgs, dims = 45:48, reduction = "harmony")
dev.off()

tiff("2pcs_49_thru_51.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(hammond.mgs, dims = 49:51, reduction = "harmony")
dev.off()

# Run UMAP
hammond.mgs <- RunUMAP(hammond.mgs, reduction = "harmony", dims = 1:50)

# plot UMAP color by sample variables
DimPlot(hammond.mgs, group.by = "sample_id") + theme(aspect.ratio = 1)
DimPlot(hammond.mgs, group.by = "timepoint") + theme(aspect.ratio = 1)

# plot UMAP overlaid with gene expression
tiff("mrc1_featureplot.tiff", res = 600, height = 6, width = 6, units = "in")
FeaturePlot(hammond.mgs, features = "Mrc1", order = TRUE) + theme(aspect.ratio = 1)
dev.off()

tiff("ms4a7_featureplot.tiff", res = 600, height = 6, width = 6, units = "in")
FeaturePlot(hammond.mgs, features = "Ms4a7", order = TRUE) + theme(aspect.ratio = 1)  
dev.off()

tiff("tmem119_featureplot.tiff", res = 600, height = 6, width = 6, units = "in")
FeaturePlot(hammond.mgs, features = "Tmem119", order = TRUE) + theme(aspect.ratio = 1)
dev.off()

tiff("p2ry12_featureplot.tiff", res = 600, height = 6, width = 6, units = "in")
FeaturePlot(hammond.mgs, features = "P2ry12", order = TRUE) + theme(aspect.ratio = 1)
dev.off()

tiff("lcp1_featureplot.tiff", res = 600, height = 6, width = 6, units = "in")
FeaturePlot(hammond.mgs, features = "Lcp1", order = TRUE) + theme(aspect.ratio = 1)
dev.off()

# perform hierarchical clustering using HGC
hammond.mgs <- FindNeighbors(hammond.mgs, reduction = "harmony", dims = 1:50)
hammond.mgs <- FindClusteringTree(hammond.mgs, graph.type = "SNN")
tree <- hammond.mgs@graphs$ClusteringTree
tree$height = log(tree$height + 1)

## have to throw out singlets that did not cluster
initial.k <- cutree(tree, k = 68, order_clusters_as_data = FALSE) 
table(initial.k)

# 67 singlets; going to get subdendrograms so we can move along clustering everything but the singlets; 
subtree.list <- get_subdendrograms(as.dendrogram(tree), 68, order_clusters_as_data = FALSE)
dend.without.singlets <- subtree.list[[1]]

test = get_leaves_attr(dend.without.singlets, "label", simplify = FALSE) #

# plot dendrogram
dend <- dend.without.singlets %>%
  dendextend::set("branches_k_color", k=21, value = c(ggsci::pal_d3("category20")(20), "#000000")) %>%
  dendextend::set("branches_lwd", 6) %>%
  dendextend::set("labels_colors", "white")

tiff("subdendrogram_microglia_k21.tiff", res = 600, height = 15, width = 20, units = "in")
plot(dend)
dev.off()

# save dendrogram
saveRDS(dend, "initial_clustering_dendrogram_microglia.rds")

# to maximize heterogeneity discoverable by clustering, we cut the dendrogram to produce a large number of clusters
# then, to focus on cell groups that are most robust, we combine small clusters (< approx. 500 cells) with other daughter cluster if they were split into two terminal branches at last division; 
# This yields 17 clusters after merging smaller clusters
k21 <- cutree(dend, k = 21, order_clusters_as_data = FALSE)
table(k21) # get cell counts to identify clusters to merge

k21.df <- k21 %>% enframe(name = "barcode", value = "cluster") 
k21.c1 <- k21.df %>% filter(cluster == 1) %>% pull(barcode)
k21.c2 <- k21.df %>% filter(cluster == 2)%>% pull(barcode)
k21.c3 <- k21.df %>% filter(cluster == 3)%>% pull(barcode)
k21.c4 <- k21.df %>% filter(cluster == 4)%>% pull(barcode)
k21.c5 <- k21.df %>% filter(cluster == 5)%>% pull(barcode)
k21.c6 <- k21.df %>% filter(cluster == 6)%>% pull(barcode)
k21.c7 <- k21.df %>% filter(cluster == 7)%>% pull(barcode)
k21.c8 <- k21.df %>% filter(cluster == 8)%>% pull(barcode)
k21.c9 <- k21.df %>% filter(cluster == 9)%>% pull(barcode)
k21.c10 <- k21.df %>% filter(cluster == 10)%>% pull(barcode)
k21.c11 <- k21.df %>% filter(cluster == 11)%>% pull(barcode)
k21.c12 <- k21.df %>% filter(cluster == 12)%>% pull(barcode)
k21.c13 <- k21.df %>% filter(cluster == 13)%>% pull(barcode)
k21.c14 <- k21.df %>% filter(cluster == 14)%>% pull(barcode)
k21.c15 <- k21.df %>% filter(cluster == 15)%>% pull(barcode)
k21.c16 <- k21.df %>% filter(cluster == 16)%>% pull(barcode)
k21.c17 <- k21.df %>% filter(cluster == 17)%>% pull(barcode)
k21.c18 <- k21.df %>% filter(cluster == 18)%>% pull(barcode)
k21.c19 <- k21.df %>% filter(cluster == 19)%>% pull(barcode)
k21.c20 <- k21.df %>% filter(cluster == 20)%>% pull(barcode)
k21.c21 <- k21.df %>% filter(cluster == 21)%>% pull(barcode)

DimPlot(hammond.mgs, cells.highlight = k21.c1) + theme(aspect.ratio = 1) 
DimPlot(hammond.mgs, cells.highlight = k21.c2) + theme(aspect.ratio = 1) 
DimPlot(hammond.mgs, cells.highlight = k21.c3) + theme(aspect.ratio = 1) 
DimPlot(hammond.mgs, cells.highlight = k21.c4) + theme(aspect.ratio = 1) 
DimPlot(hammond.mgs, cells.highlight = k21.c5) + theme(aspect.ratio = 1) 
DimPlot(hammond.mgs, cells.highlight = k21.c6) + theme(aspect.ratio = 1) 
DimPlot(hammond.mgs, cells.highlight = k21.c7) + theme(aspect.ratio = 1) 
DimPlot(hammond.mgs, cells.highlight = k21.c8) + theme(aspect.ratio = 1) 
DimPlot(hammond.mgs, cells.highlight = k21.c9) + theme(aspect.ratio = 1) 
DimPlot(hammond.mgs, cells.highlight = k21.c10) + theme(aspect.ratio = 1) 
DimPlot(hammond.mgs, cells.highlight = k21.c11) + theme(aspect.ratio = 1) 
DimPlot(hammond.mgs, cells.highlight = k21.c12) + theme(aspect.ratio = 1) 
DimPlot(hammond.mgs, cells.highlight = k21.c13) + theme(aspect.ratio = 1) 
DimPlot(hammond.mgs, cells.highlight = k21.c14) + theme(aspect.ratio = 1)
DimPlot(hammond.mgs, cells.highlight = k21.c15) + theme(aspect.ratio = 1) 
DimPlot(hammond.mgs, cells.highlight = k21.c16) + theme(aspect.ratio = 1)
DimPlot(hammond.mgs, cells.highlight = k21.c17) + theme(aspect.ratio = 1) 
DimPlot(hammond.mgs, cells.highlight = k21.c18) + theme(aspect.ratio = 1) 
DimPlot(hammond.mgs, cells.highlight = k21.c19) + theme(aspect.ratio = 1) 
DimPlot(hammond.mgs, cells.highlight = k21.c20) + theme(aspect.ratio = 1)
DimPlot(hammond.mgs, cells.highlight = k21.c21) + theme(aspect.ratio = 1) 


# reassign cell groups with final cluster labels
k21.c1.s <- hammond.mgs[,k21.c1]
k21.c1.s$k21_ident <- 1
k21.c2.s <- hammond.mgs[,k21.c2]
k21.c2.s$k21_ident <- 2
k21.c3.s <- hammond.mgs[,k21.c3]
k21.c3.s$k21_ident <- 2
k21.c4.s <- hammond.mgs[,k21.c4]
k21.c4.s$k21_ident <- 3
k21.c5.s <- hammond.mgs[,k21.c5]
k21.c5.s$k21_ident <- 4
k21.c6.s <- hammond.mgs[,k21.c6]
k21.c6.s$k21_ident <- 5
k21.c7.s <- hammond.mgs[,k21.c7]
k21.c7.s$k21_ident <- 6
k21.c8.s <- hammond.mgs[,k21.c8]
k21.c8.s$k21_ident <- 7
k21.c9.s <- hammond.mgs[,k21.c9]
k21.c9.s$k21_ident <- 8
k21.c10.s <- hammond.mgs[,k21.c10]
k21.c10.s$k21_ident <- 8
k21.c11.s <- hammond.mgs[,k21.c11]
k21.c11.s$k21_ident <- 8
k21.c12.s <- hammond.mgs[,k21.c12]
k21.c12.s$k21_ident <- 9
k21.c13.s <- hammond.mgs[,k21.c13]
k21.c13.s$k21_ident <- 10
k21.c14.s <- hammond.mgs[,k21.c14]
k21.c14.s$k21_ident <- 11
k21.c15.s <- hammond.mgs[,k21.c15]
k21.c15.s$k21_ident <- 12
k21.c16.s <- hammond.mgs[,k21.c16]
k21.c16.s$k21_ident <- 13
k21.c17.s <- hammond.mgs[,k21.c17]
k21.c17.s$k21_ident <- 13
k21.c18.s <- hammond.mgs[,k21.c18]
k21.c18.s$k21_ident <- 14
k21.c19.s <- hammond.mgs[,k21.c19]
k21.c19.s$k21_ident <- 15
k21.c20.s <- hammond.mgs[,k21.c20]
k21.c20.s$k21_ident <- 16
k21.c21.s <- hammond.mgs[,k21.c21]
k21.c21.s$k21_ident <- 16


k21.merged <- merge(k21.c1.s, c(k21.c2.s, k21.c3.s, k21.c4.s, k21.c5.s, k21.c6.s, k21.c7.s, k21.c8.s, k21.c9.s, k21.c10.s, k21.c11.s, k21.c12.s, k21.c13.s,
                                k21.c14.s, k21.c15.s, k21.c16.s, k21.c17.s, k21.c18.s, k21.c19.s, k21.c20.s, k21.c21.s))
hammond.mgs$final_ident <- k21.merged$k21_ident

rm(list = c("k21.c1.s", "k21.c2.s", "k21.c3.s", "k21.c4.s", "k21.c5.s", "k21.c6.s", "k21.c7.s", "k21.c8.s", "k21.c9.s", "k21.c10.s", "k21.c11.s", 
            "k21.c12.s", "k21.c13.s", "k21.c14.s", "k21.c15.s", "k21.c16.s", "k21.c17.s", "k21.c18.s", "k21.c19.s", "k21.c20.s", "k21.c21.s",
            "k21.merged"))
gc()

# remove singlets unassigned to clusters from dataset
hammond.mgs <- hammond.mgs[,!is.na(hammond.mgs$final_ident)] 
hammond.mgs 

Idents(hammond.mgs) <- hammond.mgs$final_ident
Idents(hammond.mgs) <- factor(Idents(hammond.mgs), levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16))

# plot UMAP overlaid with cluster labels
tiff("hammond_final_clusters.tiff", res = 600, height = 6, width= 6, units = "in")
DimPlot(hammond.mgs, cols = ggsci::pal_d3("category20")(16)) + theme(aspect.ratio = 1)
dev.off()

# save seurat object
saveRDS(hammond.mgs, "hammond_mgs_final_seurat_object.rds")

# differential expression testing for cluster markers
plan("multiprocess")
options(future.globals.maxSize= 10000*1024^2)

final.markers <- FindAllMarkers(hammond.mgs)
final.markers <- final.markers %>% arrange(cluster, desc(avg_log2FC))
saveRDS(final.markers, "final_clustering_markers.rds")
write.csv(final.markers, "final_clustering_markers.csv")

### heatmap of scaled expression of top markers for each cluster
top.marks <- final.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top.marks <- top.marks$gene

plan("sequential")

hammond.mgs <- ScaleData(hammond.mgs, features = rownames(hammond.mgs))

# get average expression of each of these markers within each cluster
exp.mat <- AverageExpression(hammond.mgs, features = top.marks, slot = "scale.data")
exp.mat <- as.matrix(exp.mat$RNA)

# plot heatmap
top.a = HeatmapAnnotation(Cluster = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16"), annotation_name_side = "left", col = list(Cluster = c("1" = ggsci::pal_d3("category20")(16)[1],
                                                                                                                                                                                    "2" = ggsci::pal_d3("category20")(16)[2],
                                                                                                                                                                                    "3" = ggsci::pal_d3("category20")(16)[3],
                                                                                                                                                                                    "4" = ggsci::pal_d3("category20")(16)[4],
                                                                                                                                                                                    "5" = ggsci::pal_d3("category20")(16)[5],
                                                                                                                                                                                    "6" = ggsci::pal_d3("category20")(16)[6],
                                                                                                                                                                                    "7" = ggsci::pal_d3("category20")(16)[7],
                                                                                                                                                                                    "8" = ggsci::pal_d3("category20")(16)[8],
                                                                                                                                                                                    "9" = ggsci::pal_d3("category20")(16)[9],
                                                                                                                                                                                    "10" = ggsci::pal_d3("category20")(16)[10],
                                                                                                                                                                                    "11" = ggsci::pal_d3("category20")(16)[11],
                                                                                                                                                                                    "12" = ggsci::pal_d3("category20")(16)[12],
                                                                                                                                                                                    "13" = ggsci::pal_d3("category20")(16)[13],
                                                                                                                                                                                    "14" = ggsci::pal_d3("category20")(16)[14],
                                                                                                                                                                                    "15" = ggsci::pal_d3("category20")(16)[15],
                                                                                                                                                                                    "16" = ggsci::pal_d3("category20")(16)[16]
)),
show_annotation_name = TRUE,
show_legend = FALSE,
height = unit(10, "mm"))

hmap <- Heatmap(exp.mat, name = "z-scored\navg. expr.",
                col = circlize::colorRamp2(c(-2, -1, 0, 1, 2), viridis::viridis_pal(option = "magma")(5)), # rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu"
                row_names_gp = gpar(fontsize = 4),
                top_annotation = top.a,
                show_column_names = FALSE,
                cluster_columns = FALSE,
                border = TRUE,
                cluster_rows = FALSE,
                #cluster_row_slices = TRUE,
                #split = rep(c("Microglia", "Neuronal", "Macrophage/Monocyte"), each = 30),
                show_row_names = TRUE,
                #right_annotation = gene_ha
)

# save plot
pdf("final_clustering_heatmap.pdf", height = 10, width = 10)
hmap
dev.off()

# after finding cluster markers, cluster 2 appear to be BAMs and cluster 3 are neurons
# next subset out the putative microglia clusters 
hammond.final.mgs <- subset(hammond.mgs, idents = c(1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16))
dim(hammond.final.mgs) 

# subset out microglia into Mrc1+ and Mrc1- 
mrc1.pos <- WhichCells(hammond.final.mgs, expression = `Mrc1` > 0) 
mrc1.neg <- WhichCells(hammond.final.mgs, expression = `Mrc1` == 0) 
length(mrc1.pos)/length(colnames(hammond.final.mgs)) 

mrc1.mgs <- hammond.final.mgs[,mrc1.pos]
Idents(mrc1.mgs) <- "Mrc1+"
other.mgs <- hammond.final.mgs[,mrc1.neg]
Idents(other.mgs) <- "Mrc1-"

combo <- merge(mrc1.mgs, other.mgs)
combo

saveRDS(hammond.final.mgs, "final_seurat_object_mgs_only.rds") # save Seurat object of just microglia
saveRDS(combo, "final_seurat_object_mrc1pos_v_neg_mgs.rds") # save Seurat object of just microglia clusters with Mrc1 identities

# Examine cluster proportions within each sample over time:
table(hammond.final.mgs$final_ident, hammond.final.mgs$timepoint)
prop.df = data.frame("sample" = hammond.final.mgs$sample_id, "timepoint" = hammond.final.mgs$timepoint, "cluster" = hammond.final.mgs$final_ident)

# examine proportion of microglia that are in cluster 10 over time 
clust.10.counts <- prop.df %>% filter(cluster == 10) %>% group_by(timepoint, sample) %>% tally()
total.counts <- prop.df %>% group_by(timepoint, sample) %>% tally()

c10.prop.df <- merge(clust.10.counts, total.counts, by = "sample")
c10.prop.df <- c10.prop.df %>% mutate(prop = n.x/n.y)

c10.prop.df <- c10.prop.df[,-4]
colnames(c10.prop.df) <- c("sample", "timepoint", "n_cluster10", "n_total", "prop_cluster10")

# save cluster 10 microglia proportions data frame
saveRDS(c10.prop.df, "cluster10_proportions.rds")
write.csv(c10.prop.df, "cluster10_proportions.csv")

# examine proportion of mrc1+ microglia over time
mrc1.prop.df = data.frame("sample" = combo$sample_id, "timepoint" = combo$timepoint, "cluster" = combo@active.ident)
mrc1.counts <- mrc1.prop.df %>% filter(cluster == "Mrc1+") %>% group_by(timepoint, sample) %>% tally()
total.mrc1.counts <- mrc1.prop.df %>% group_by(timepoint, sample) %>% tally()

mrc1.prop.df <- merge(mrc1.counts, total.mrc1.counts, by = "sample")
mrc1.prop.df <- mrc1.prop.df %>% mutate(prop = n.x/n.y)

mrc1.prop.df <- mrc1.prop.df[,-4]
colnames(mrc1.prop.df) <- c("sample", "timepoint", "n_cluster10", "n_total", "prop_cluster10")

# save MRC1 proportions data frame
saveRDS(mrc1.prop.df, "mrc1_proportions.rds")
write.csv(mrc1.prop.df, "mrc1_proportions.csv")

# Differential expression testing between mrc1+ and mrc1- microglia 
mgs <- readRDS("final_seurat_object_mrc1pos_v_neg_mgs.rds")

mrc1pos.neg.markers <- FindMarkers(mgs, ident.1 = "Mrc1+", ident.2 = "Mrc1-")
mrc1pos.neg.markers <- mrc1pos.neg.markers %>% arrange(desc(avg_log2FC))

# save marker data frames
saveRDS(mrc1pos.neg.markers, "mrc1_pos_v_neg_all_markers.rds")
write.csv(mrc1pos.neg.markers, "mrc1_pos_v_neg_all_markers.csv")

