### this script provides code for the preprocessing and analysis of the Kracht et al. (2020) single-cell RNA-sequencing data
# note: file paths are relative and may change depdending on local machine

# load required packages
library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(egg)
library(stringr)
library(plyr)
library(EnsDb.Hsapiens.v79)
library(HGC)
library(dendextend)
library(tibble)
library(dplyr)
library(unikn)
library(tibble)
library(ComplexHeatmap)
library(ggpubr)

# generate list of count matrix files (points to a folder containing the _completeCounts.txt files for each batch, obtained from the NCBI GEO Series GSE141862)
setwd("../input_data/GSE141862_RAW/")
data.files <- list.files("../GSE141862_RAW/")

# read in sample metadata (in this case an Excel spreadsheet containing metadata copied from the NCBI GEO repository Series Matrix Files);
# this spreadsheet contained the date, Sample_characteristics_ch1, Sample_title, and geo_accession columns from the original file
sample.metadata <- readxl::read_excel("../../Kracht_seq_R/sample_data.xlsx")
sample.metadata <- sample.metadata[1:239,]

# read in supplementary data table (from Kracht et al. 2020 manuscript; including the gestational age, sex, and sample date for each batch)
other.metadata <- readxl::read_excel("../../aba5906_tables_s1_s10.xlsx", col_names = TRUE)
length(unique(sample.metadata$`!Sample_characteristics_ch1...1`)) # 20 total samples

# read in count matrix files
setwd("../input_data/GSE141862_RAW/")

input.data <- lapply(data.files, read.table, header = TRUE, row.names = 1, sep = "\t")
cell.count <- 0
genes.included <- c()
for(i in 1:length(input.data)){
  cell.count = cell.count + ncol(input.data[[i]]) 
  genes.included <- unique(c(rownames(input.data[[i]]), genes.included))
  colnames(input.data[[i]]) <- paste(paste("run", i, sep = "_"), colnames(input.data[[i]]), sep = "_")
  input.data[[i]]$gene <- rownames(input.data[[i]])
}

setwd("../../Kracht_seq_R/")

# final data: 17,295 cells; 44,200 genes
combined.input.data <- join_all(input.data, by = 'gene', type = 'full')

# 1. Convert from ensembl gene id (from the original count matrices) to human gene.symbol
ensembl.genes <- combined.input.data$gene
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
colnames(geneIDs) <- c("gene_symbol", "gene")

count.matrix <- combined.input.data

count.matrix <- merge(count.matrix, geneIDs, by = "gene", all.x = TRUE)

# examine how successfully gene id's were found for each ensemble id
test = count.matrix$gene
any(is.na(test))

test = count.matrix$gene_symbol
any(is.na(test))

no.genesymbol <- count.matrix[is.na(count.matrix$gene_symbol),] # 1,192 ENSEMBL gene ids have no corresponding gene symbol; excluding from downstream analysis
count.matrix <- count.matrix[!is.na(count.matrix$gene_symbol),] # 43,008 rows now; 

dups <- count.matrix[duplicated(count.matrix$gene_symbol),]$gene_symbol # 1030 duplicates gene symbols (meaning more than one ensemble id refered to the same gene symbol)
length(unique(dups)) # comprised of only 154 genes

# add transcripts counts together for ensemble ids that link to the same gene symbol
count.matrix.nondup <- count.matrix[!duplicated(count.matrix$gene_symbol),]
count.matrix.dups <- count.matrix[duplicated(count.matrix$gene_symbol),]

count.matrix.dups.resolved = ddply(count.matrix.dups, "gene_symbol", numcolwise(sum))

count.matrix.nondup[,1] <- NULL
rownames(count.matrix.nondup) <- count.matrix.nondup[,17296]
count.matrix.nondup[,17296] <- NULL

rownames(count.matrix.dups.resolved) <- count.matrix.dups.resolved[,1]
count.matrix.dups.resolved[,1] <- NULL

count.matrix.final <- rbind(count.matrix.nondup, count.matrix.dups.resolved)

# save final combined count matrices
saveRDS(count.matrix.final, "final_merged_count_matrix.rds")
write.csv(count.matrix.final, "final_merged_count_matrix.csv")

# get sample_id and timepoint info for each; 
colnames(sample.metadata) <- c("ID", "full_time", "sample_title", "GEO_accession")
merged.metadata <- merge(sample.metadata, other.metadata, by = "ID", all = TRUE)

colnames(count.matrix.final)

# make final metadata dataframe with final barcode identifiers for each cell
cell.data.df <- NULL
cell.data.df <- data.frame(ID = character(), full_time = character(), sample_title = character(), 
                           GEO_accession = character(), `Gestational age (weeks + days)` = character(),
                           `Gestational age` = integer(), Sex = character(), ATAC = character(), 
                           `Single cell` = character())
#rownames(cell.data.df) <- ncol(count.matrix.final)
for(i in 1:nrow(merged.metadata)){
  n.repeats <- ncol(input.data[[i]])
  for(x in 1:(n.repeats-1)){
    cell.data.df <- rbind(cell.data.df, merged.metadata[i,])
  }
}

cell.data.df$barcode <- colnames(count.matrix.final)

# save final metadata data frame
saveRDS(cell.data.df, "cell_metadata.rds")
write.csv(cell.data.df, "cell_metadata.csv")

count.matrix.final[is.na(count.matrix.final)] <- 0 # convert all NA's to zero; (not all of the sequencing runs include all genes; add in zero's for undetected genes in a given batch)

# save actual revised final count matrix (with NA's replaced with 0's)
saveRDS(count.matrix.final, "final_final_count_matrix.rds")
write.csv(count.matrix.final, "final_final_count_matrix.csv")


##########################################
# # to save memory; you can clear environment & restart R session; 
# ### read in final count matrix and metadata files
count.matrix.final <- readRDS("../input_data/GSE141862_RAW/final_final_count_matrix.rds")
cell.data.df <- readRDS("../input_data/GSE141862_RAW/cell_metadata.rds")

# Initialize the Seurat object with the raw (non-normalized) data.
mg.data <- CreateSeuratObject(count.matrix.final)

# add cell metadata to seurat object
mg.data$sample_id <- cell.data.df$ID
mg.data$timepoint <- cell.data.df$`Gestational week`
mg.data$subsample_id <- cell.data.df$GEO_accession
mg.data$sex <- cell.data.df$Sex

# add percentage of UMIs per cell mapping to mitochondrial genome (QC metric)
mg.data$percent.mt <- PercentageFeatureSet(object = mg.data, pattern = "^MT.")
mg.data$percent.mt

# split Seurat object into list of seurat objects by sample
mg.list <- SplitObject(mg.data, split.by = "sample_id")

# identify union of highly variable genes across samples
features <- SelectIntegrationFeatures(object.list = mg.list, nfeatures = 3000)

# plot qc plots to determine what gene/UMI/pct. mito cutoffs to use: 
p1 <- VlnPlot(mg.data, features = "nFeature_RNA", group.by = "sample_id", pt.size = 0) + labs(x = NULL) + 
  ggtitle("Number of genes detected") + NoLegend() + geom_hline(yintercept = 3000, linetype = 'dotted', col = 'red') +
  geom_hline(yintercept = 200, linetype = 'dotted', col = 'red')# cut off 200<x<3,000 genes seems reasonable
p2 <- VlnPlot(mg.data, features = "nCount_RNA", group.by = "sample_id", pt.size = 0) + labs(x = NULL) + 
  ggtitle("Number of UMIs") + NoLegend() + geom_hline(yintercept = 90000, linetype = 'dotted', col = 'red') # cut off 90,000 umis
p3 <- VlnPlot(mg.data, features = "percent.mt", group.by = "sample_id", pt.size = 0) + labs(x = NULL) + 
  ggtitle("Percent of UMIs mapping to mitochondrial genome") + NoLegend() + geom_hline(yintercept = 10,
                                                                                       linetype = 'dotted', col = 'red') # cut off 10% seems reasonable

ggarrange(p1, p2, p3, ncol = 1)

tiff("qc_plots.tiff", res = 300, height = 12, width = 15, units = "in")
ggarrange(p1, p2, p3, ncol = 1)
dev.off()

postscript("qc_plots.ps", height = 12, width = 15)
ggarrange(p1, p2, p3, ncol = 1)
dev.off()

# based on QC plots, filtering out cells with <200 genes or > 3000 unique genes; cells with > 90,000 UMIs, and cells with greater than 10% of reads mapping to mitochondrial genes 
mg.data <- subset(mg.data, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & nCount_RNA < 90000 & percent.mt < 10) # 16,203 cells left after filtering
mg.data # seurat object is 42,132 features x 16,203 cells

# clear memory, except for mg.data seurat object and variable features vector
rm(list = setdiff(ls(), c("mg.data", "features")))
gc()

# log-normalize data & scale features;
mg.data <- NormalizeData(mg.data) 
mg.data <- ScaleData(mg.data, features = features, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt")) # regressing out technical variables (UMI counts, genes detected counts, and percent.mt) 

# run PCA 
mg.data <- RunPCA(mg.data, features = features, npcs = 50, verbose = TRUE) 

# again clear memory to save space
rm(list = setdiff(ls(), "mg.data"))
gc()

# Integration of sequencing batches (subsample_id) by Harmony
set.seed(79683) # integration process is stochastic so make sure to set a random seed so results are reproducible
mg.data <- mg.data %>% RunHarmony(c("subsample_id"), plot_convergence = TRUE)

# plot elbow plots before and after integration; 
tiff("harmony_elbowplot.tiff", units = "in", width = 5, height = 5, res = 300)
ElbowPlot(mg.data, ndims = 50, reduction = "harmony")
dev.off()

tiff("uncorrected_pca_elbowplot.tiff", units = "in", width = 5, height = 5, res = 300) 
ElbowPlot(mg.data, ndims = 50, reduction = "pca") 
dev.off()

# also visualize genes contributing to each principal component after integration 
tiff("pcs_1_thru_4.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(mg.data, dims = 1:4, reduction = "harmony")
dev.off()

tiff("pcs_5_thru_8.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(mg.data, dims = 5:8, reduction = "harmony")
dev.off()

tiff("pcs_9_thru_12.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(mg.data, dims = 9:12, reduction = "harmony")
dev.off()

tiff("pcs_13_thru_16.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(mg.data, dims = 13:16, reduction = "harmony")
dev.off()

tiff("pcs_17_thru_20.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(mg.data, dims = 17:20, reduction = "harmony")
dev.off()

tiff("pcs_21_thru_24.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(mg.data, dims = 21:24, reduction = "harmony")
dev.off()

tiff("pcs_25_thru_28.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(mg.data, dims = 25:28, reduction = "harmony")
dev.off()

tiff("pcs_29_thru_32.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(mg.data, dims = 29:32, reduction = "harmony")
dev.off()

# the above plots suggest 15 principal components would yield an informative reduction; 
# next, perform UMAP using 15 harmony-corrected principal components
mg.data <- RunUMAP(mg.data, reduction = "harmony", dims = 1:15)

# visualize UMAP colored by sample variables
DimPlot(mg.data, group.by = "sample_id") + theme(aspect.ratio = 1)
DimPlot(mg.data, group.by = "timepoint") + theme(aspect.ratio = 1)
FeaturePlot(mg.data, features = "timepoint") + theme(aspect.ratio = 1)
FeaturePlot(mg.data, features = "nCount_RNA") + theme(aspect.ratio = 1)
FeaturePlot(mg.data, features = "nFeature_RNA") + theme(aspect.ratio = 1)
FeaturePlot(mg.data, features = "percent.mt") + theme(aspect.ratio = 1)

# save Seurat object
saveRDS(mg.data, "initial_reduction_mg_data.rds") # saved seurat object (original dim reduction; no clutering)

# to identify transcriptional subtypes of microglia (and other cell types), we next perform hierarchical clustering using HGC
mg.data <- FindNeighbors(mg.data, reduction = "harmony", dims = 1:15)
mg.data <- FindClusteringTree(mg.data, graph.type = "SNN")

tree <- mg.data@graphs$ClusteringTree
tree$height = log(tree$height + 1)

cell.labels <- data.frame(Sample = mg.data$sample_id,
                          Time = mg.data$timepoint)


# plot cluster dendrogram
tiff("test_dendrogram.tiff", units = "in", height = 10, width = 25, res = 300)
HGC.PlotDendrogram(tree = tree, k = 7, plot.label = TRUE, labels = cell.labels)
dev.off()

# cutting the tree at 7 clusters
k7 <- cutree(tree, k = 7, order_clusters_as_data = FALSE)
table(k7)

# retrieve cell barcodes for each cluster
k7.df <- k7 %>% enframe(name = "barcode", value = "cluster") 
k7.c1 <- k7.df %>% dplyr::filter(cluster == 1) %>% pull(barcode)
k7.c2 <- k7.df %>% dplyr::filter(cluster == 2)%>% pull(barcode)
k7.c3 <- k7.df %>% dplyr::filter(cluster == 3)%>% pull(barcode)
k7.c4 <- k7.df %>% dplyr::filter(cluster == 4)%>% pull(barcode)
k7.c5 <- k7.df %>% dplyr::filter(cluster == 5)%>% pull(barcode)
k7.c6 <- k7.df %>% dplyr::filter(cluster == 6)%>% pull(barcode)
k7.c7 <- k7.df %>% dplyr::filter(cluster == 7)%>% pull(barcode)

# highlight clusters on UMAP
DimPlot(mg.data, cells.highlight = k7.c1) + theme(aspect.ratio = 1) # erythrocytes
DimPlot(mg.data, cells.highlight = k7.c2) + theme(aspect.ratio = 1) # neurons
DimPlot(mg.data, cells.highlight = k7.c3) + theme(aspect.ratio = 1) # microglia
DimPlot(mg.data, cells.highlight = k7.c4) + theme(aspect.ratio = 1) # microglia
DimPlot(mg.data, cells.highlight = k7.c5) + theme(aspect.ratio = 1) # microglia 
DimPlot(mg.data, cells.highlight = k7.c6) + theme(aspect.ratio = 1) # macrophage/monocytes
DimPlot(mg.data, cells.highlight = k7.c7) + theme(aspect.ratio = 1) # microglia

# reassign cells with cluster labels
k7.c1.s <- mg.data[,k7.c1]
k7.c1.s$k7_ident <- 1
k7.c2.s <- mg.data[,k7.c2]
k7.c2.s$k7_ident <- 2
k7.c3.s <- mg.data[,k7.c3]
k7.c3.s$k7_ident <- 3
k7.c4.s <- mg.data[,k7.c4]
k7.c4.s$k7_ident <- 4
k7.c5.s <- mg.data[,k7.c5]
k7.c5.s$k7_ident <- 5
k7.c6.s <- mg.data[,k7.c6]
k7.c6.s$k7_ident <- 6
k7.c7.s <- mg.data[,k7.c7]
k7.c7.s$k7_ident <- 7

k7.merged <- merge(k7.c1.s, c(k7.c2.s, k7.c3.s, k7.c4.s, k7.c5.s, k7.c6.s, k7.c7.s))
mg.data$k7_ident <- k7.merged$k7_ident

# again clear memory of excess objects
rm(list = c("k7.c1.s", "k7.c2.s", "k7.c3.s", "k7.c4.s", "k7.c5.s", "k7.c6.s", "k7.merged"))
gc()

# plot UMAP overlaid with cluster identities
DimPlot(mg.data, group.by = "k7_ident") + theme(aspect.ratio = 1)

# Differential expression testing to identify marker genes for each cluster
Idents(mg.data) <- mg.data$k7_ident
k7.markers <- FindAllMarkers(mg.data)
k7.markers <- k7.markers %>% dplyr::arrange(cluster, desc(avg_log2FC))

# rename clusters
mg.data <- RenameIdents(mg.data, `3` = "A", `7` = "B", `4` = "C", `5` = "D", `6` = "E", `2` = "F", `1` = "G")

# stacked violin plot highlighting marker genes for clusters & cell identities
tiff("initial_clusters_stackedviolinplot.tiff", res = 600, height = 10, width = 10, units = "in")
VlnPlot(mg.data, features = c("P2RY12", "CX3CR1", "C3", "CSF1R", "AIF1", "F13A1", "LYVE1", "MRC1", 
                              "KIF5A", "MEG3", "NFIB", "HBG2", "HBA2", "HBB"),
        stack = TRUE, flip = TRUE, sort = FALSE, cols = rev(viridis::viridis_pal(option = "plasma")(14))) + labs(x = NULL, y = "Log-Normalized Expression") + 
  theme(axis.text.x.bottom = element_text(angle = 0, hjust = 0.5)) + NoLegend()
dev.off()

postscript("initial_clusters_stackedviolinplot.ps", height = 10, width = 7)
VlnPlot(mg.data, features = c("P2RY12", "CX3CR1", "C3", "CSF1R", "AIF1", "F13A1", "LYVE1", "MRC1", 
                              "KIF5A", "MEG3", "NFIB", "HBG2", "HBA2", "HBB"),
        stack = TRUE, flip = TRUE, sort = FALSE, cols = rev(viridis::viridis_pal(option = "plasma")(14))) + labs(x = NULL, y = "Log-Normalized Expression") + 
  theme(axis.text.x.bottom = element_text(angle = 0, hjust = 0.5)) + NoLegend()
dev.off()

# make color palette
pal <- ggsci::pal_npg()(7)

# make UMAP plots with final color scheme
postscript("mg_dimplot_initial_clusters.ps", height = 6, width = 6)
DimPlot(mg.data, cols = pal) + theme(aspect.ratio = 1)
dev.off()

tiff("mg_dimplot_initial_clusters.tiff", height = 6, width = 6, res = 300, units = "in")
DimPlot(mg.data, cols = pal) + theme(aspect.ratio = 1)
dev.off()

# make color palette for matching UMAP cluster colors to hierarchical clustering dendrogram
adj.pal <- c(ggsci::pal_npg()(7)[7], ggsci::pal_npg()(7)[6], ggsci::pal_npg()(5)[1], ggsci::pal_npg()(7)[3],
             ggsci::pal_npg()(7)[4], ggsci::pal_npg()(5)[5], ggsci::pal_npg()(5)[2])


dend <- as.dendrogram(tree) %>%
  dendextend::set("branches_k_color", k=7, value = adj.pal) %>%
  dendextend::set("branches_lwd", 6) %>%
  dendextend::set("labels_colors", "white")

# plot dendrogram with new colors
tiff("initial_dendrogram.tiff", units = "in", height = 10, width = 10, res = 600)
plot(dend)
dev.off()

# save seurat object; marker data frames; and dendrogram objects
saveRDS(mg.data, "kracht_allcelltypes.rds") # seurat object with clustering results
saveRDS(k7.markers, "initial_clustering_markers.rds") # 7 cluster markers data frame
write.csv(k7.markers, "initial_clustering_markers.csv") # "
saveRDS(tree, "initial_clustering_dendrogram.rds") # clustering dendrogram


mg.subset <- subset(mg.data, idents = c("A", "B", "C", "D", "E")) # subsetting out the five clusters with higher expression of macrophage/microglia markers

# subclustering the subsetted cells

# split object by sample id
mg.subset.list <- SplitObject(mg.subset, split.by = "sample_id")

# find new variable features for subsetted object
features2 <- SelectIntegrationFeatures(object.list = mg.subset.list, nfeatures = 2000)

# clear memory
rm("mg.subset.list")
gc()

# scle new variable features and regress out technical variables
mg.subset <- ScaleData(mg.subset, features = features2, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt")) 

# run PCA
mg.subset <- RunPCA(mg.subset, features = features2, npcs = 50, verbose = TRUE)

# integration with harmony
set.seed(69492)
mg.subset <- mg.subset %>% RunHarmony(c("subsample_id"), plot_convergence = TRUE)

mg.subset 

# plot elbow plots for integration results
tiff("harmony_elbowplot2.tiff", units = "in", width = 5, height = 5, res = 300)
ElbowPlot(mg.subset, ndims = 50, reduction = "harmony")
dev.off()

tiff("uncorrected_pca_elbowplot2.tiff", units = "in", width = 5, height = 5, res = 300) 
ElbowPlot(mg.subset, ndims = 50, reduction = "pca") 
dev.off()

# visualize genes associated with each principal component
tiff("subclustering_pcs_1_thru_4.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(mg.subset, dims = 1:4, reduction = "harmony")
dev.off()

tiff("subclustering_pcs_5_thru_8.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(mg.subset, dims = 5:8, reduction = "harmony")
dev.off()

tiff("subclustering_pcs_9_thru_12.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(mg.subset, dims = 9:12, reduction = "harmony")
dev.off()

tiff("subclustering_pcs_13_thru_16.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(mg.subset, dims = 13:16, reduction = "harmony")
dev.off()

tiff("subclustering_pcs_17_thru_20.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(mg.subset, dims = 17:20, reduction = "harmony")
dev.off()

tiff("subclustering_pcs_21_thru_24.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(mg.subset, dims = 21:24, reduction = "harmony")
dev.off()

tiff("subclustering_pcs_25_thru_28.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(mg.subset, dims = 25:28, reduction = "harmony")
dev.off()

tiff("subclustering_pcs_29_thru_32.tiff", units = "in", height = 10, width = 10, res = 300)
VizDimLoadings(mg.subset, dims = 29:32, reduction = "harmony")
dev.off()

# based on the above plots, we'll use 10 principal components
# run UMAP
mg.subset <- RunUMAP(mg.subset, reduction = "harmony", dims = 1:10)

# plot UMAP colored by metadata features
DimPlot(mg.subset, group.by = "sample_id") + theme(aspect.ratio = 1)
DimPlot(mg.subset, group.by = "timepoint") + theme(aspect.ratio = 1)
FeaturePlot(mg.subset, features = "timepoint") + theme(aspect.ratio = 1)

# overlay UMAP with some microglia and other features for visualization
postscript("p2ry12_featureplot.ps", height = 6, width = 6)
FeaturePlot(mg.subset, features = "P2RY12", order = TRUE) + theme(aspect.ratio = 1)
dev.off()

postscript("csf1r_featureplot.ps", height = 6, width = 6)
FeaturePlot(mg.subset, features = "CSF1R", order = TRUE) + theme(aspect.ratio = 1)
dev.off()

postscript("cx3cr1_featureplot.ps", height = 6, width = 6)
FeaturePlot(mg.subset, features = "CX3CR1", order = TRUE) + theme(aspect.ratio = 1)
dev.off()

#
tiff("p2ry12_featureplot.tiff", height = 6, width = 6, units="in", res = 600)
FeaturePlot(mg.subset, features = "P2RY12", order = TRUE) + theme(aspect.ratio = 1)
dev.off()

tiff("tmem119_featureplot.tiff", height = 6, width = 6, units="in", res = 600)
FeaturePlot(mg.subset, features = "TMEM119", order = TRUE) + theme(aspect.ratio = 1)
dev.off()

tiff("cx3cr1_featureplot.tiff", height = 6, width = 6, units="in", res = 600)
FeaturePlot(mg.subset, features = "CX3CR1", order = TRUE) + theme(aspect.ratio = 1)
dev.off()

tiff("mrc1_featureplot.tiff", height = 6, width = 6, units="in", res = 600)
FeaturePlot(mg.subset, features = "MRC1", order = TRUE) + theme(aspect.ratio = 1)
dev.off()

postscript("mrc1_featureplot.ps", height = 6, width = 6)
FeaturePlot(mg.subset, features = "MRC1") + theme(aspect.ratio = 1)
dev.off()

postscript("ms4a7_featureplot.ps", height = 6, width = 6)
FeaturePlot(mg.subset, features = "MS4A7") + theme(aspect.ratio = 1)
dev.off()

postscript("spi1_featureplot.ps", height = 6, width = 6)
FeaturePlot(mg.subset, features = "SPI1") + theme(aspect.ratio = 1)
dev.off()


postscript("tmem119_featureplot.ps", height = 6, width = 6)
FeaturePlot(mg.subset, features = "TMEM119", order = TRUE) + theme(aspect.ratio = 1)
dev.off()

postscript("lcp1_featureplot.ps", height = 6, width = 6)
FeaturePlot(mg.subset, features = "LCP1") + theme(aspect.ratio = 1)
dev.off()

postscript("f13a1_featureplot.ps", height = 6, width = 6)
FeaturePlot(mg.subset, features = "F13A1") + theme(aspect.ratio = 1)
dev.off()

# perform hierarchical clustering
mg.subset <- FindNeighbors(mg.subset, reduction = "harmony", dims = 1:10)
mg.subset <- FindClusteringTree(mg.subset, graph.type = "SNN")

tree <- mg.subset@graphs$ClusteringTree
tree$height = log(tree$height + 1)

# create new color palette for new clusters
pal <- unikn::usecol(pal_unikn_pair[1:10])

dend <- as.dendrogram(tree) %>%
  dendextend::set("branches_k_color", k=10, pal) %>%
  dendextend::set("branches_lwd", 6) %>%
  dendextend::set("labels_colors", "white")

# plot dendrogram
tiff("subclustering_dendrogram.tiff", units = "in", height = 10, width = 10, res = 600)
plot(dend)
dev.off()

# cut dendrogram at 10 clusters
k10 <- cutree(tree, k = 10, order_clusters_as_data = FALSE)
table(k10)

k10.df <- k10 %>% enframe(name = "barcode", value = "cluster") 
k10.c1 <- k10.df %>% dplyr::filter(cluster == 1) %>% pull(barcode)
k10.c2 <- k10.df %>% dplyr::filter(cluster == 2)%>% pull(barcode)
k10.c3 <- k10.df %>% dplyr::filter(cluster == 3)%>% pull(barcode)
k10.c4 <- k10.df %>% dplyr::filter(cluster == 4)%>% pull(barcode)
k10.c5 <- k10.df %>% dplyr::filter(cluster == 5)%>% pull(barcode)
k10.c6 <- k10.df %>% dplyr::filter(cluster == 6)%>% pull(barcode)
k10.c7 <- k10.df %>% dplyr::filter(cluster == 7)%>% pull(barcode)
k10.c8 <- k10.df %>% dplyr::filter(cluster == 8)%>% pull(barcode)
k10.c9 <- k10.df %>% dplyr::filter(cluster == 9)%>% pull(barcode)
k10.c10 <- k10.df %>% dplyr::filter(cluster == 10)%>% pull(barcode)

# overlay cluster identities on UMAP 
DimPlot(mg.subset, cells.highlight = k10.c1) + theme(aspect.ratio = 1) 
DimPlot(mg.subset, cells.highlight = k10.c2) + theme(aspect.ratio = 1) 
DimPlot(mg.subset, cells.highlight = k10.c3) + theme(aspect.ratio = 1) 
DimPlot(mg.subset, cells.highlight = k10.c4) + theme(aspect.ratio = 1) 
DimPlot(mg.subset, cells.highlight = k10.c5) + theme(aspect.ratio = 1) 
DimPlot(mg.subset, cells.highlight = k10.c6) + theme(aspect.ratio = 1) 
DimPlot(mg.subset, cells.highlight = k10.c7) + theme(aspect.ratio = 1) 
DimPlot(mg.subset, cells.highlight = k10.c8) + theme(aspect.ratio = 1) 
DimPlot(mg.subset, cells.highlight = k10.c9) + theme(aspect.ratio = 1) 
DimPlot(mg.subset, cells.highlight = k10.c10) + theme(aspect.ratio = 1) 

# assign identities to cells 
k10.c1.s <- mg.subset[,k10.c1]
k10.c1.s$k10_ident <- 1
k10.c2.s <- mg.subset[,k10.c2]
k10.c2.s$k10_ident <- 2
k10.c3.s <- mg.subset[,k10.c3]
k10.c3.s$k10_ident <- 3
k10.c4.s <- mg.subset[,k10.c4]
k10.c4.s$k10_ident <- 4
k10.c5.s <- mg.subset[,k10.c5]
k10.c5.s$k10_ident <- 5
k10.c6.s <- mg.subset[,k10.c6]
k10.c6.s$k10_ident <- 6
k10.c7.s <- mg.subset[,k10.c7]
k10.c7.s$k10_ident <- 7
k10.c8.s <- mg.subset[,k10.c8]
k10.c8.s$k10_ident <- 8
k10.c9.s <- mg.subset[,k10.c9]
k10.c9.s$k10_ident <- 9
k10.c10.s <- mg.subset[,k10.c10]
k10.c10.s$k10_ident <- 10

k10.merged <- merge(k10.c1.s, c(k10.c2.s, k10.c3.s, k10.c4.s, k10.c5.s, k10.c6.s, 
                                k10.c7.s, k10.c8.s, k10.c9.s, k10.c10.s))
mg.subset$k10_ident <- k10.merged$k10_ident

# clear memory 
rm(list = c("k10.c1.s", "k10.c2.s", "k10.c3.s", "k10.c4.s", "k10.c5.s", "k10.c6.s", 
            "k10.c7.s", "k10.c8.s", "k10.c9.s", "k10.c10.s", "k10.merged"))
gc()

# plot UMAP with new cluster identities
DimPlot(mg.subset, group.by = "k10_ident", cols = pal) + theme(aspect.ratio = 1) 

Idents(mg.subset) <- mg.subset$k10_ident
Idents(mg.subset) <- factor(Idents(mg.subset), levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))

# differential expression testing to identify cluster markers
k10.markers <- FindAllMarkers(mg.subset)
k10.markers <- k10.markers %>% dplyr::arrange(cluster, desc(avg_log2FC))

# save rsults and subclustered Seurat object
saveRDS(k10.markers, "final_cluster_markers.rds")
write.csv(k10.markers, "final_cluster_markers.csv")
saveRDS(mg.subset, "kracht_final_object.rds") # subclustered seurat object; (with microglia, BAMs, and leukocytes)

# plot stacked violin plots of cluster markers
postscript("final_clusters_stackedviolinplot.ps", height = 10, width = 10)
VlnPlot(mg.data, features = c("P2RY12", "CX3CR1", "C3", "CSF1R", "AIF1", "F13A1", "LYVE1", "MRC1", 
                              "KIF5A", "MEG3", "NFIB", "HBG2", "HBA2", "HBB"),
        stack = TRUE, flip = TRUE, sort = FALSE, cols = rev(viridis::viridis_pal(option = "plasma")(14))) + labs(x = NULL, y = "Log-Normalized Expression") + 
  theme(axis.text.x.bottom = element_text(angle = 0, hjust = 0.5)) + NoLegend()
dev.off()

# plot UMAP 
postscript("mg_dimplot_final_clusters.ps", height = 6, width = 6)
DimPlot(mg.subset, cols = pal) + theme(aspect.ratio = 1)
dev.off()

tiff("mg_dimplot_final_clusters.tiff", height = 6, width = 6, res = 300, units = "in")
DimPlot(mg.subset, cols = pal) + theme(aspect.ratio = 1)
dev.off()

########
# plot heatmap of top 10 marker genes (by log2FC) for each cluster
top.marks <- k10.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top.marks <- top.marks$gene

mg.subset <- ScaleData(mg.subset, features = rownames(mg.subset))
exp.mat <- AverageExpression(mg.subset, features = top.marks, slot = "scale.data")
exp.mat <- as.matrix(exp.mat$RNA)

top.a = HeatmapAnnotation(Cluster = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), annotation_name_side = "left", col = list(Cluster = c("1" = pal[1],
                                                                                                                                                "2" = pal[2],
                                                                                                                                                "3" = pal[3],
                                                                                                                                                "4" = pal[4],
                                                                                                                                                "5" = pal[5],
                                                                                                                                                "6" = pal[6],
                                                                                                                                                "7" = pal[7],
                                                                                                                                                "8" = pal[8],
                                                                                                                                                "9" = pal[9],
                                                                                                                                                "10" = pal[10])),
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

# save heatmap
pdf("final_clustering_heatmap.pdf", height = 10, width = 10)
hmap
dev.off()


# save seurat objects
saveRDS(mg.data, "final_seurat_object.rds") # original all cell types seurat object
saveRDS(mg.subset, "final_seurat_object_clustering2.rds") # subclustering seurat object

