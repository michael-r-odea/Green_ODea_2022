# this script contains the code used for the analysis of bulk RNA-sequencing from Utz et al. (2020), and testing enrichment of the resultant microglia and BAM gene sets in the Hammond et al. and Kracht et al. datasets to identify clusters as microglia or BAMs
# note: file paths are relative to local machine and may change

# load required packages
library(dplyr)
library(ggplot2)
library(Seurat)
library(ComplexHeatmap)
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(DEGreport)
library(tximport)
library(ggplot2)
library(ggrepel)
library(biomaRt)
library(stringr)
library(readxl)
library(ggpubr)
library(viridis)
library(tidyverse)
library(viridis)
library(AUCell)
library(unikn)


## using data from bulk RNA-seq data to identify markers which denote microglia versus BAMs
# dataset: Utz et al. (2020): bulk RNA-seq of BAMS vs Microglia: F4/80hi (CD45+F4/80hiCD11b+Ly6C-Ly6G-) and CD11bhi (CD45+F4/80+CD11bhiLy6C-Ly6G-) 
# data includes samples from several timepoints throughout development (from E10.5, E11.5, E12.5, E14.5, E16.5, and E18.5);
# bulk data was retrieved from SRA BioProject PRJNA612456

## bulk data was first pseudoaligned using Salmon (see salmon_script.sh)

# next, we use DESeq2 to identify differentially expressed genes between microglia and BAMs;

# get salmon transcript count file paths
files <- list.files(path = ".", pattern = "quant.sf", recursive = TRUE)
samples <- sub("/.*", "", files)

# read in SRA metadata table
utz.metadata <- read.table("../SraRunTable.txt", sep = ",", header = TRUE)
utz.metadata <- utz.metadata[match(samples, utz.metadata$GEO_Accession..exp.),]

names(files) <- samples

# convert transcripts isoforms to transcripts
ids <- read.delim(files[1], sep="\t", header=T) 
ids <- as.character(ids[,1])
require(stringr)
ids.strip <- str_replace(ids, "([.][0-9])", "")

# Create a mart object
mart <- useDataset("mmusculus_gene_ensembl", mart = useMart("ENSEMBL_MART_ENSEMBL"))

# Get gene symbol and Ensembl gene IDs
tx2gene <- getBM(
  filters= "ensembl_transcript_id", 
  attributes= c("ensembl_transcript_id", "external_gene_name"),
  values= ids.strip,
  mart= mart)

# import salmon counts files
txi <- tximport(files, type="salmon", txIn = TRUE, txOut = FALSE, tx2gene=tx2gene, importer=read_tsv, ignoreTxVersion=TRUE)
attributes(txi) # check attributes

# process metadata 
rownames(utz.metadata) <- utz.metadata$GEO_Accession..exp.
utz.metadata$Cell_type <- sub("Brain - ", "", utz.metadata$Cell_type)
utz.metadata$Cell_type <- sub("border-associated macrophages", "BAM", utz.metadata$Cell_type)

# Create DEseq2 object
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = utz.metadata,
                                   design = ~Cell_type)

# estimate size factors for normalizing by library size
dds <- estimateSizeFactors(ddsTxi)

# extracting DESeq2-normalized counts 
normalized_counts <- counts(dds, normalized=TRUE)
# save normalized counts as txt file
write.table(normalized_counts, file="bulk_utz_normalized_counts.txt", sep="\t", quote=F, col.names=NA)

# for plotting sample distances, we're performing rlog transformation
rld <- rlog(dds, blind=TRUE)

# extract PCA results from plotPCA deseq2 function
plotPCA(rld, intgroup=c("Cell_type", "Development_stage")) # PC1: 44% variance; PC2: 30% variance 
pca.data <- plotPCA(rld, intgroup=c("Cell_type", "Development_stage"), returnData = TRUE)

# plot PCA 
tiff("pca_plot_celltype_and_timepoint.tiff", units = "in", height = 6, width = 6, res = 300) # principal component 1 is developmental age and principal component 2 is Microglia vs BAM
ggplot(pca.data,aes(x=PC1,y=PC2,color=Development_stage)) + 
  geom_point(aes(shape=Cell_type), size = 3) + 
  scale_color_viridis(option = "magma", discrete = TRUE) +
  theme_pubr() + 
  theme(aspect.ratio = 1, legend.position = "right") + 
  labs(color = "Age", shape = "Cell Type", x = "PC1 (44% variance)", y = "PC2 (30% variance)")
dev.off()


# plot sample distance heatmap
rld_mat <- assay(rld)   
sampleDists <- dist(t(assay(rld)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rld$GEO_Accession..exp.
colnames(sampleDistMatrix) <- rld$GEO_Accession..exp.

# get metadata for heatmap annotations
ha = HeatmapAnnotation(df = data.frame("Cell_type" = rld$Cell_type, "Age" = rld$Development_stage), col = list("Cell_type" = c("Microglia" = "darkred", "BAM" = "navyblue"), 
                                                                                                               "Age" = c("E10.5" = "#000004FF",
                                                                                                                         "E11.5" = "#3B0F70FF",
                                                                                                                         "E12.5" = "#8C2981FF",
                                                                                                                         "E14.5" = "#DE4968FF",
                                                                                                                         "E16.5" = "#FE9F6DFF",
                                                                                                                         "E18.5" = "#FCFDBFFF")))

hmap <- Heatmap(sampleDistMatrix, name = "Distance", col = rev(RColorBrewer::brewer.pal(9, "Blues")), top_annotation = ha)

# plot heatmap
tiff("sample_dists_heatmap.tiff", units = "in", height = 8, width = 9, res = 300)
hmap
dev.off()

# now run DESeq2 to find differentially expressed genes
dds <- DESeq(dds)

# dispersion estimates plot
tiff("dispersion_estimates_plot.tiff", units = "in", height = 6, width = 10, res = 300)
plotDispEsts(dds)
dev.off()

# get results contrasting variable: Cell_Type_Microglia_vs_BAM
# lfc shrunken results
mg.v.bam <- lfcShrink(dds, coef="Cell_type_Microglia_vs_BAM", type="apeglm")
mg.v.bam.results <- mg.v.bam %>% data.frame() %>% filter(abs(log2FoldChange) > 1) %>% filter(padj < 0.05) %>% arrange(desc(baseMean), desc(log2FoldChange))

mg.v.bam.results %>% filter(log2FoldChange > 1) %>% nrow(.) # 574 genes are more highly expressed in microglia in development than BAMs
mg.v.bam.results %>% filter(log2FoldChange < -1) %>% nrow(.) # 325 genes are more highly expressed in BAMs in development than microglia

mg.v.bam.results$gene <- rownames(mg.v.bam.results)

# save results
saveRDS(mg.v.bam.results, "mg_v_bam_deseq2_results.rds")
write.csv(mg.v.bam.results, "mg_v_bam_deseq2_results.csv")

# create heatmap of differentially expressed genes
tmp <- normalized_counts[rownames(mg.v.bam.results),]
tmp <- t(scale(t(tmp)))

ha = HeatmapAnnotation(df = data.frame("Cell_type" = dds$Cell_type, "Age" = dds$Development_stage), col = list("Cell_type" = c("Microglia" = "darkred", "BAM" = "navyblue"), 
                                                                                                               "Age" = c("E10.5" = "#000004FF",
                                                                                                                         "E11.5" = "#3B0F70FF",
                                                                                                                         "E12.5" = "#8C2981FF",
                                                                                                                         "E14.5" = "#DE4968FF",
                                                                                                                         "E16.5" = "#FE9F6DFF",
                                                                                                                         "E18.5" = "#FCFDBFFF")))
# highlight a few notable genes
gene_list = c("Cx3cr1", "Hexb", "Ctsd", "P2ry12", "P2ry13", "Tmem119", "Siglech", "Trem2", "Sall1", "Mpeg1", 
              "F13a1", "Mrc1", "Lyve1", "Cd163", "Siglec1", "Folr2", "Ccl2", "Cd38")

gene_annot = rowAnnotation(genes = anno_mark(at = which(rownames(tmp) %in% gene_list), labels = gene_list))

gene_ha = rowAnnotation(genes = anno_mark(at = match(gene_list, rownames(tmp)), labels = gene_list))

hmap <- Heatmap(tmp, name = "z-scored\nnorm. exp.", col = circlize::colorRamp2(breaks = c(-2, -1, 0, 1, 2), colors = rev(RColorBrewer::brewer.pal(5, "RdBu"))), top_annotation = ha,
                right_annotation = gene_ha,
                show_row_names = FALSE,
                show_column_names = FALSE,
                heatmap_legend_param = list(
                  at = c(-2, -1, 0, 1, 2), 
                  labels = c("-2", "-1", "0", "1", "≥2")
                ))
# plot heatmap
tiff("mg_v_bam_markers_heatmap.tiff", units = "in", height = 8, width = 9, res = 300)
hmap
dev.off()

# extract marker gene lists
mg.markers <- mg.v.bam.results %>% filter(log2FoldChange > 1) %>% rownames(.)
bam.markers <- mg.v.bam.results %>% filter(log2FoldChange < -1) %>% rownames(.)

# read in Hammond data
hammond <- readRDS("hammond_mgs_final_seurat_object.rds")

# find which genes from the MG and BAM datasets were detected in the hammond dataset
mg.markers.hammond <- mg.markers[mg.markers %in% rownames(hammond)]
bam.markers.hammond <- bam.markers[bam.markers %in% rownames(hammond)]

# keep only microglia genes which were expressed in at least 20% of the cells in the dataset (since the vast majority of cells in this dataset are microglia)
genes.to.keep <- Matrix::rowSums(x = hammond@assays$RNA@counts > 0) >= floor(x = 0.20 * ncol(x = hammond@assays$RNA@counts)) 
genes.subset <- rownames(hammond)[genes.to.keep] 

# out of 512 mg markers detected in the hammond dataset; 147 were present in at least 20% of cells in the dataset
mg.markers.hammond.pruned <- mg.markers.hammond[mg.markers.hammond %in% genes.subset]

# to fucus on the most robust expression differences, we use the top 50 microglia genes
mg.markers.hammond.pruned <- mg.markers.hammond[1:50]

# now choosing BAM gene set; we don't set a minimum expression threshold for the BAM gene set because BAMs should represent a much smaller fraction of the dataset)
# we'll use the top 50 DE BAM genes as well
bam.markers.hammond.pruned <- bam.markers.hammond[1:50]

# calculate average expression of the microglia and BAM gene set genes within each cluster
hammond.markers.avg.exp <- as.data.frame(AverageExpression(hammond, features = c(mg.markers.hammond.pruned, bam.markers.hammond.pruned), slot = "data"))
colnames(hammond.markers.avg.exp) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16")

# remove cluster 3, neurons
hammond.markers.avg.exp <- hammond.markers.avg.exp[,c("1", "2", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16")]

# log normalize
hammond.markers.avg.exp <- log1p(hammond.markers.avg.exp)

# scale
z.ham.exp <- t(scale(t(hammond.markers.avg.exp)))

# create heatmap 
ha = HeatmapAnnotation(df = data.frame("Cluster" = colnames(z.ham.exp)), col = list("Cluster" = c("1" = ggsci::pal_d3("category20")(16)[1],
                                                                                                  "2" = ggsci::pal_d3("category20")(16)[2],
                                                                                                  #"3" = ggsci::pal_d3("category20")(16)[3],
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
                                                                                                  "16" = ggsci::pal_d3("category20")(16)[16])))

# highlight certain genes
gene_list = c("Cx3cr1", "Hexb", "Ctsd", "P2ry12", "P2ry13", "Tmem119", "Siglech", "Trem2", "Sall1", "Mpeg1", 
              "F13a1", "Mrc1", "Lyve1", "Cd163", "Siglec1", "Folr2", "Ccl2", "Cd38")

gene_ha = rowAnnotation(genes = anno_mark(at = match(gene_list, rownames(z.ham.exp)), labels = gene_list, labels_gp = gpar(fontsize = 8)))

hmap <- Heatmap(z.ham.exp, name = "z-scored\navg. exp.", col = circlize::colorRamp2(breaks = c(-2, -1, 0, 1, 2), colors = viridis::viridis_pal(option = "magma")(5)), 
                top_annotation = ha,
                show_row_names = FALSE,
                show_column_names = TRUE,
                cluster_rows = TRUE,
                row_split = c(rep("MG", length(mg.markers.hammond.pruned)), rep("BAM", length(bam.markers.hammond.pruned))),
                row_gap = unit(5, "mm"),
                right_annotation = gene_ha,
                heatmap_legend_param = list(
                  at = c(-2, -1, 0, 1, 2), 
                  labels = c("≤-2", "-1", "0", "1", "≥2")
                ))

hmap

# plot heatmap
tiff("mg_v_bam_markers_heatmap.tiff", units = "in", height = 8, width = 6, res = 300)
hmap
dev.off()

# calculate MG and BAM gene set scores for each single-cell
cells_rankings <- AUCell_buildRankings(hammond@assays$RNA@data, nCores=1, plotStats=TRUE)
cells_rankings
cells_AUC <- AUCell_calcAUC(list("MGs" = mg.markers.hammond.pruned, "BAMs" = bam.markers.hammond.pruned), cells_rankings) #

hammond$module1_auc <- as.vector(t(cells_AUC["MGs",]@assays@data@listData[["AUC"]]))
hammond$module2_auc <- as.vector(t(cells_AUC["BAMs",]@assays@data@listData[["AUC"]]))

# extract single-cell scores
tmp.df <- data.frame("Cluster" = hammond@active.ident, "MG_score" = hammond$module1_auc, "BAM_score" = hammond$module2_auc)
tmp.df <- tmp.df %>% filter(Cluster != 3)

# plot cluster gene set score boxplots 
tiff("mg_genes_auc_score_boxplot.tiff", units = "in", height = 4, width = 4, res = 300)
ggboxplot(tmp.df, "Cluster", "MG_score",
          color = "black", fill = "Cluster", notch = TRUE, palette =c(ggsci::pal_d3("category20")(16)[1:2], ggsci::pal_d3("category20")(16)[4:16]),
          outlier.shape = NA) + labs(y = "Gene Set Enrichment Score (AUC)", title = "Microglia Gene Set") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20)) + NoLegend() + 
  coord_flip()
dev.off()

tiff("bam_genes_auc_score_boxplot.tiff", units = "in", height = 4, width = 4, res = 300)
ggboxplot(tmp.df, "Cluster", "BAM_score",
          color = "black", fill = "Cluster", notch = TRUE, palette =c(ggsci::pal_d3("category20")(16)[1:2], ggsci::pal_d3("category20")(16)[4:16]),
          outlier.shape = NA) + labs(y = "Gene Set Enrichment Score (AUC)", title = "BAM Gene Set") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20)) + NoLegend() + 
  coord_flip()
dev.off()

# read in hammond seurat object with microglia labeled as Mrc1+ or Mrc1-
hammond.mrc1.pos.v.neg <- readRDS("final_seurat_object_mrc1pos_v_neg_mgs.rds")

# subset out BAM cells from the main seurat object
hammond.bam <- subset(hammond, idents = 2)

# merge BAMs with purified microglia seurat object 
hammond.2 <- merge(hammond.mrc1.pos.v.neg, hammond.bam)

# calculate average expression of microglia and BAM gene sets in Cluster 2 BAMs, Mrc1+ microglia, and Mrc1- microglia
hammond.markers.avg.exp.2 <- as.data.frame(AverageExpression(hammond.2, features = c(mg.markers.hammond.pruned, bam.markers.hammond.pruned), slot = "data"))
colnames(hammond.markers.avg.exp.2) <- c("2", "Mrc1-", "Mrc1+")

# log normalize
hammond.markers.avg.exp.2 <- log1p(hammond.markers.avg.exp.2)

# scale
z.ham.exp.2 <- t(scale(t(hammond.markers.avg.exp.2)))

# create heatmap
ha = HeatmapAnnotation(df = data.frame("Cluster" = colnames(z.ham.exp.2)), col = list("Cluster" = c("2" = ggsci::pal_d3("category20")(16)[2],
                                                                                                    "Mrc1-" = "purple",
                                                                                                    "Mrc1+" = "darkred")))

gene_list = c("Cx3cr1", "Hexb", "Ctsd", "P2ry12", "P2ry13", "Tmem119", "Siglech", "Trem2", "Sall1", "Mpeg1", 
              "F13a1", "Mrc1", "Lyve1", "Cd163", "Siglec1", "Folr2", "Ccl2", "Cd38")

gene_ha = rowAnnotation(genes = anno_mark(at = match(gene_list, rownames(z.ham.exp.2)), labels = gene_list,  labels_rot = 135, labels_gp = gpar(fontsize = 8)))

hmap <- Heatmap(z.ham.exp.2, name = "z-scored\navg. exp.", col = circlize::colorRamp2(breaks = c(-1, -0.5, 0, 0.5, 1), colors = viridis::viridis_pal(option = "magma")(5)), bottom_annotation = ha,
                show_row_names = FALSE,
                show_column_names = TRUE,
                cluster_rows = TRUE,
                row_split = c(rep("MG", length(mg.markers.hammond.pruned)), rep("BAM", length(bam.markers.hammond.pruned))),
                row_gap = unit(5, "mm"),
                right_annotation = gene_ha,
                heatmap_legend_param = list(
                  at = c(-1, -0.5, 0, 0.5, 1), 
                  labels = c("≤-1", "-0.5", "0", "0.5", "≥1")
                ))

hmap

# plot heatmap
tiff("mg_v_bam_markers_heatmap2.tiff", units = "in", height = 8, width = 4, res = 300)
hmap
dev.off()


# calculate single-cell gene set scores for this new seurat object
cells_rankings.hammond.mrc1.pos.v.neg <- AUCell_buildRankings(hammond.mrc1.pos.v.neg@assays$RNA@data, nCores=1, plotStats=TRUE)
cells_rankings.hammond.mrc1.pos.v.neg
cells_AUC.hammond.mrc1.pos.v.neg <- AUCell_calcAUC(list("MGs" = mg.markers.hammond.pruned, "BAMs" = bam.markers.hammond.pruned), cells_rankings.hammond.mrc1.pos.v.neg) #

hammond.mrc1.pos.v.neg$module1_auc <- as.vector(t(cells_AUC.hammond.mrc1.pos.v.neg["MGs",]@assays@data@listData[["AUC"]]))
hammond.mrc1.pos.v.neg$module2_auc <- as.vector(t(cells_AUC.hammond.mrc1.pos.v.neg["BAMs",]@assays@data@listData[["AUC"]]))

# extract scores
tmp.df.2 <- data.frame("Cluster" = hammond.mrc1.pos.v.neg@active.ident, "MG_score" = hammond.mrc1.pos.v.neg$module1_auc, "BAM_score" = hammond.mrc1.pos.v.neg$module2_auc)
bam.df <- tmp.df %>% filter(Cluster == 2)
tmp.df.3 <- rbind(bam.df, tmp.df.2)
tmp.df.3$Cluster <- factor(tmp.df.3$Cluster, levels = c("Mrc1+", "Mrc1-", "2"))

# plot gene set score boxplots
tiff("mg_genes_auc_score_boxplot_mrc1pos_v_neg.tiff", units = "in", height = 2, width = 4, res = 300)
ggboxplot(tmp.df.3, "Cluster", "MG_score",
          color = "black", fill = "Cluster", notch = TRUE, palette =c("darkred", "purple", ggsci::pal_d3("category20")(16)[2]),
          outlier.shape = NA) + labs(y = "Gene Set Enrichment Score (AUC)", title = "Microglia Gene Set") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20)) + NoLegend() + coord_flip()
dev.off()

tiff("bam_genes_auc_score_boxplot_mrc1pos_v_neg.tiff", units = "in", height = 2, width = 4, res = 300)
ggboxplot(tmp.df.3, "Cluster", "BAM_score",
          color = "black", fill = "Cluster", notch = TRUE, palette =c("darkred", "purple", ggsci::pal_d3("category20")(16)[2]),
          outlier.shape = NA) + labs(y = "Gene Set Enrichment Score (AUC)", title = "BAM Gene Set") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20)) + NoLegend() + coord_flip()
dev.off()


# Kracht analysis: 

# load seurat object
kracht <- readRDS("../Kracht_seq/Kracht_seq_R/kracht_final_object.rds")

# load color palette
pal <- unikn::usecol(pal_unikn_pair[1:10])

# function for converting mouse to human genes
convertMouseGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  return(genesV2)
}

# convert mouse microglia genes list to human genes
mg.markers.orthologs.df <- convertMouseGeneList(mg.markers)
colnames(mg.markers.orthologs.df) <- c("mouse_gene", "human_gene")
mg.markers.orthologs.df <- mg.markers.orthologs.df[match(mg.markers, mg.markers.orthologs.df$mouse_gene),]
mg.markers.orthologs.df.o <- mg.markers.orthologs.df
mg.markers.orthologs.df <- mg.markers.orthologs.df[!(duplicated(mg.markers.orthologs.df$mouse_gene) | duplicated(mg.markers.orthologs.df$mouse_gene, fromLast = TRUE)), ] # remove all many to one orthologs
mg.markers.orthologs.df <- mg.markers.orthologs.df[!(duplicated(mg.markers.orthologs.df$human_gene) | duplicated(mg.markers.orthologs.df$human_gene, fromLast = TRUE)), ] # remove all one to many orthlogs
mg.markers.orthologs <- mg.markers.orthologs.df$human_gene
mg.markers.orthologs <- mg.markers.orthologs[!is.na(mg.markers.orthologs)]
# create list of gene names from this dataframe, only including genes also present in our seurat object
mg.markers.orthologs.kracht <- mg.markers.orthologs[mg.markers.orthologs %in% rownames(kracht)]
# we'll focus on the top 50 genes to examine most robust expression differences
mg.markers.orthologs.kracht.pruned <- mg.markers.orthologs.kracht[1:50]

# convert mouse BAM genes list to human genes
bam.markers.orthologs.df <- convertMouseGeneList(bam.markers)
colnames(bam.markers.orthologs.df) <- c("mouse_gene", "human_gene")
bam.markers.orthologs.df <- bam.markers.orthologs.df[match(bam.markers, bam.markers.orthologs.df$mouse_gene),]
bam.markers.orthologs.df.o <- bam.markers.orthologs.df
bam.markers.orthologs.df <- bam.markers.orthologs.df[!(duplicated(bam.markers.orthologs.df$mouse_gene) | duplicated(bam.markers.orthologs.df$mouse_gene, fromLast = TRUE)), ] # remove all many to one orthologs
bam.markers.orthologs.df <- bam.markers.orthologs.df[!(duplicated(bam.markers.orthologs.df$human_gene) | duplicated(bam.markers.orthologs.df$human_gene, fromLast = TRUE)), ] # remove all one to many orthlogs
bam.markers.orthologs <- bam.markers.orthologs.df$human_gene
bam.markers.orthologs <- bam.markers.orthologs[!is.na(bam.markers.orthologs)]
bam.markers.orthologs.kracht <- bam.markers.orthologs[bam.markers.orthologs %in% rownames(kracht)]
bam.markers.orthologs.kracht.pruned <- bam.markers.orthologs.kracht[1:50]

# calculate average expression of each gene in both gene lists across Kracht dataset clusters
kracht.markers.avg.exp <- as.data.frame(AverageExpression(kracht, features = c(mg.markers.orthologs.kracht.pruned, bam.markers.orthologs.kracht.pruned), slot = "data"))
colnames(kracht.markers.avg.exp) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")

# remove cluster 8 (leukocytes)
kracht.markers.avg.exp <- kracht.markers.avg.exp[,c("1", "2", "3", "4", "5", "6", "7", "9", "10")] 

# log-normalize
kracht.markers.avg.exp <- log1p(kracht.markers.avg.exp)

# Scale
z.kracht.exp <- t(scale(t(kracht.markers.avg.exp)))

# remove NAs
z.kracht.exp <- z.kracht.exp[rowSums(is.na(z.kracht.exp)) != ncol(z.kracht.exp), ]

# create heatmap
ha = HeatmapAnnotation(df = data.frame("Cluster" = colnames(z.kracht.exp)), col = list("Cluster" = c("1" = pal[1],
                                                                                                     "2" = pal[2],
                                                                                                     "3" = pal[3],
                                                                                                     "4" = pal[4],
                                                                                                     "5" = pal[5],
                                                                                                     "6" = pal[6],
                                                                                                     "7" = pal[7],
                                                                                                     #"8" = pal[8],
                                                                                                     "9" = pal[9],
                                                                                                     "10" = pal[10])))
# highlight notable genes
gene_list = c("CX3CR1", "HEXB", "CTSD", "P2RY12", "P2RY13", "TMEM119", "SIGLECH", "TREM2", "SALL1", "MPEG1", 
              "F13A1", "MRC1", "LYVE1", "CD163", "SIGLEC1", "FOLR2", "CCL2", "CD38")

gene_ha = rowAnnotation(genes = anno_mark(at = match(gene_list, rownames(z.kracht.exp)), labels = gene_list, labels_gp = gpar(fontsize = 8)))

hmap <- Heatmap(z.kracht.exp, name = "z-scored\navg. exp.", col = circlize::colorRamp2(breaks = c(-2, -1, 0, 1, 2), colors = viridis::viridis_pal(option = "magma")(5)), 
                top_annotation = ha,
                show_row_names = FALSE,
                show_column_names = TRUE,
                cluster_rows = TRUE,
                row_split = c(rep("MG", length(mg.markers.orthologs.kracht.pruned)), rep("BAM", length(bam.markers.orthologs.kracht.pruned))),
                row_gap = unit(5, "mm"),
                right_annotation = gene_ha,
                heatmap_legend_param = list(
                  at = c(-2, -1, 0, 1, 2), 
                  labels = c("≤-2", "-1", "0", "1", "≥2")
                ))

hmap

# plot heatmap
tiff("mg_v_bam_markers_heatmap_kracht.tiff", units = "in", height = 8, width = 6, res = 300)
hmap
dev.off()

# calculate single-cell gene set scores for both gene sets
cells_rankings.kracht <- AUCell_buildRankings(kracht@assays$RNA@data, nCores=1, plotStats=TRUE)
cells_rankings.kracht
cells_AUC.kracht <- AUCell_calcAUC(list("MGs" = mg.markers.orthologs.kracht.pruned, "BAMs" = bam.markers.orthologs.kracht.pruned), cells_rankings.kracht) #

kracht$module1_auc <- as.vector(t(cells_AUC.kracht["MGs",]@assays@data@listData[["AUC"]]))
kracht$module2_auc <- as.vector(t(cells_AUC.kracht["BAMs",]@assays@data@listData[["AUC"]]))

# extract scores
tmp.df <- data.frame("Cluster" = kracht@active.ident, "MG_score" = kracht$module1_auc, "BAM_score" = kracht$module2_auc)
tmp.df <- tmp.df %>% filter(Cluster != 8) # remove leukocyte cluster

# make boxplots of cluster gene set scores
tiff("kracht_mg_genes_auc_score_boxplot.tiff", units = "in", height = 4, width = 4, res = 300)
ggboxplot(tmp.df, "Cluster", "MG_score",
          color = "black", fill = "Cluster", notch = TRUE, palette =c(pal[1:7], pal[9:10]),
          outlier.shape = NA) + labs(y = "Gene Set Enrichment Score (AUC)", title = "Microglia Gene Set") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20)) + NoLegend() + 
  coord_flip()
dev.off()

tiff("kracht_bam_genes_auc_score_boxplot.tiff", units = "in", height = 4, width = 4, res = 300)
ggboxplot(tmp.df, "Cluster", "BAM_score",
          color = "black", fill = "Cluster", notch = TRUE, palette =c(pal[1:7], pal[9:10]),
          outlier.shape = NA) + labs(y = "Gene Set Enrichment Score (AUC)", title = "BAM Gene Set") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20)) + NoLegend() + 
  coord_flip()
dev.off()

## creating Kracht microglia seurat object labeled with MRC1+ and MRC1- identity 
kracht.mrc1.pos.v.neg2 <- subset(kracht, idents = c(1, 2, 3, 4, 5, 6, 7, 10))

mrc1.pos <- subset(kracht.mrc1.pos.v.neg2, `MRC1` > 0)
dim(mrc1.pos) 
mrc1.neg <- subset(kracht.mrc1.pos.v.neg2, `MRC1` == 0)
dim(mrc1.neg) 
ncol(kracht.mrc1.pos.v.neg2)

ncol(mrc1.pos)/ncol(kracht.mrc1.pos.v.neg2) 

Idents(mrc1.pos) <- "MRC1+"
Idents(mrc1.neg) <- "MRC1-"

kracht.mrc1.pos.v.neg2 <- merge(mrc1.pos, mrc1.neg)

### Differential expression testing between MRC1+ and MRC1- microglia
kracht.new.mrc1.pos.v.neg.markers <- FindMarkers(kracht.mrc1.pos.v.neg2, ident.1 = "MRC1+", ident.2 = "MRC1-")
kracht.new.mrc1.pos.v.neg.markers <- kracht.new.mrc1.pos.v.neg.markers %>% arrange(desc(avg_log2FC))
kracht.new.mrc1.pos.v.neg.markers$gene <- rownames(kracht.new.mrc1.pos.v.neg.markers)

# save DE results
saveRDS(kracht.new.mrc1.pos.v.neg.markers, "kracht_mrc1_pos_v_neg_markers_new.rds")
write_csv(kracht.new.mrc1.pos.v.neg.markers, "kracht_mrc1_pos_v_neg_markers_new.csv")

# create new kracht seurat object, including cluster 9 BAMs and MRC1+ and MRC1- microglia
kracht.bam <- subset(kracht, idents = 9)
kracht.2 <- merge(kracht.mrc1.pos.v.neg2, kracht.bam)

# calculate average expression of MG and BAM gene set genes across each cell group
kracht.markers.avg.exp.2 <- as.data.frame(AverageExpression(kracht.2, features = c(mg.markers.orthologs.kracht.pruned, bam.markers.orthologs.kracht.pruned), slot = "data"))
colnames(kracht.markers.avg.exp.2) <- c("9", "Mrc1-", "Mrc1+")

# log normalize
kracht.markers.avg.exp.2 <- log1p(kracht.markers.avg.exp.2)

# scale
z.kracht.exp.2 <- t(scale(t(kracht.markers.avg.exp.2)))

# remove NAs
z.kracht.exp.2 <- z.kracht.exp.2[rowSums(is.na(z.kracht.exp.2)) != ncol(z.kracht.exp.2), ]

# make heatmap
ha = HeatmapAnnotation(df = data.frame("Cluster" = colnames(z.kracht.exp.2)), col = list("Cluster" = c("9" = pal[9],
                                                                                                       "Mrc1-" = "purple",
                                                                                                       "Mrc1+" = "orange")))
# highlight interesting genes
gene_list = c("CX3CR1", "HEXB", "CTSD", "P2RY12", "P2RY13", "TMEM119", "SIGLECH", "TREM2", "SALL1", "MPEG1", 
              "F13A1", "MRC1", "LYVE1", "CD163", "SIGLEC1", "FOLR2", "CCL2", "CD38")

gene_ha = rowAnnotation(genes = anno_mark(at = match(gene_list, rownames(z.kracht.exp.2)), labels = gene_list,  labels_rot = 135, labels_gp = gpar(fontsize = 8)))

hmap <- Heatmap(z.kracht.exp.2, name = "z-scored\navg. exp.", col = circlize::colorRamp2(breaks = c(-1, -0.5, 0, 0.5, 1), colors = viridis::viridis_pal(option = "magma")(5)), bottom_annotation = ha,
                show_row_names = FALSE,
                show_column_names = TRUE,
                cluster_rows = TRUE,
                row_split = c(rep("MG", length(mg.markers.orthologs.kracht.pruned)), rep("BAM", length(bam.markers.orthologs.kracht.pruned))),
                row_gap = unit(5, "mm"),
                right_annotation = gene_ha,
                heatmap_legend_param = list(
                  at = c(-1, -0.5, 0, 0.5, 1), 
                  labels = c("≤-1", "-0.5", "0", "0.5", "≥1")
                ))

hmap

# save heatmap
tiff("kracht_mg_v_bam_markers_heatmap3.tiff", units = "in", height = 8, width = 4, res = 300)
hmap
dev.off()

# calculate single-cell gene set scores
cells_rankings.kracht.mrc1.pos.v.neg2 <- AUCell_buildRankings(kracht.mrc1.pos.v.neg2@assays$RNA@data, nCores=1, plotStats=TRUE)
cells_rankings.kracht.mrc1.pos.v.neg2
cells_AUC.kracht.mrc1.pos.v.neg2 <- AUCell_calcAUC(list("MGs" = mg.markers.orthologs.kracht.pruned, "BAMs" = bam.markers.orthologs.kracht.pruned), cells_rankings.kracht.mrc1.pos.v.neg2) #

kracht.mrc1.pos.v.neg2$module1_auc <- as.vector(t(cells_AUC.kracht.mrc1.pos.v.neg2["MGs",]@assays@data@listData[["AUC"]]))
kracht.mrc1.pos.v.neg2$module2_auc <- as.vector(t(cells_AUC.kracht.mrc1.pos.v.neg2["BAMs",]@assays@data@listData[["AUC"]]))

# get scores
tmp.df.2 <- data.frame("Cluster" = kracht.mrc1.pos.v.neg2@active.ident, "MG_score" = kracht.mrc1.pos.v.neg2$module1_auc, "BAM_score" = kracht.mrc1.pos.v.neg2$module2_auc)
bam.df <- tmp.df %>% filter(Cluster == 9)

tmp.df.3 <- rbind(bam.df, tmp.df.2)

tmp.df.3$Cluster <- factor(tmp.df.3$Cluster, levels = c("MRC1+", "MRC1-", "9"))

# plot gene set score box plots
tiff("kracht_mg_genes_auc_score_boxplot_mrc1pos_v_neg2.tiff", units = "in", height = 2, width = 4, res = 300)
ggboxplot(tmp.df.3, "Cluster", "MG_score",
          color = "black", fill = "Cluster", notch = TRUE, palette =c("orange", "purple", pal[9]),
          outlier.shape = NA) + labs(y = "Gene Set Enrichment Score (AUC)", title = "Microglia Gene Set") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20)) + NoLegend() + coord_flip()
dev.off()

tiff("kracht_bam_genes_auc_score_boxplot_mrc1pos_v_neg2.tiff", units = "in", height = 2, width = 4, res = 300)
ggboxplot(tmp.df.3, "Cluster", "BAM_score",
          color = "black", fill = "Cluster", notch = TRUE, palette =c("orange", "purple", pal[9]),
          outlier.shape = NA) + labs(y = "Gene Set Enrichment Score (AUC)", title = "BAM Gene Set") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20)) + NoLegend() + coord_flip()
dev.off()

## save lists of genes used for gene sets
saveRDS(mg.markers.orthologs.kracht.pruned, "kracht_mg_markers_used.rds")
saveRDS(bam.markers.orthologs.kracht.pruned, "kracht_bam_markers_used.rds")
write_csv(data.frame("Microglia" = mg.markers.orthologs.kracht.pruned, "BAM" = bam.markers.orthologs.kracht.pruned), "kracht_mg_bam_markers_used.csv")

saveRDS(mg.markers.hammond.pruned, "hammond_mg_markers_used.rds")
saveRDS(bam.markers.hammond.pruned, "hammond_bam_markers_used.rds")
write_csv(data.frame("Microglia" = mg.markers.hammond.pruned, "BAM" = bam.markers.hammond.pruned), "hammond_mg_bam_markers_used.csv")


# test for differential expression of mg and bam gene sets within the Hammond and Kracht datasets between Mrc1+ and Mrc1- microglia
mg.gene.set.de.results.hammond <- FindMarkers(hammond.mrc1.pos.v.neg, ident.1 = "Mrc1+", ident.2 = "Mrc1-", features = mg.markers.hammond.pruned, 
                                              logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, min.cells.feature = 0, min.cells.group = 0)

mg.gene.set.de.results.hammond <- mg.gene.set.de.results.hammond %>% arrange(desc(avg_log2FC))

bam.gene.set.de.results.hammond <- FindMarkers(hammond.mrc1.pos.v.neg, ident.1 = "Mrc1+", ident.2 = "Mrc1-", features = bam.markers.hammond.pruned, 
                                               logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, min.cells.feature = 0, min.cells.group = 0)

bam.gene.set.de.results.hammond <- bam.gene.set.de.results.hammond %>% arrange(desc(avg_log2FC))

mg.gene.set.de.results.kracht <- FindMarkers(kracht.mrc1.pos.v.neg2, ident.1 = "MRC1+", ident.2 = "MRC1-", features = mg.markers.orthologs.kracht.pruned, 
                                             logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, min.cells.feature = 0, min.cells.group = 0)

mg.gene.set.de.results.kracht <- mg.gene.set.de.results.kracht %>% arrange(desc(avg_log2FC))

bam.gene.set.de.results.kracht <- FindMarkers(kracht.mrc1.pos.v.neg2, ident.1 = "MRC1+", ident.2 = "MRC1-", features = bam.markers.orthologs.kracht.pruned, 
                                              logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, min.cells.feature = 0, min.cells.group = 0)

bam.gene.set.de.results.kracht <- bam.gene.set.de.results.kracht %>% arrange(desc(avg_log2FC))

# check for significant differences
mg.gene.set.de.results.hammond %>% filter(avg_log2FC < -0.5) %>% rownames(.) 
bam.gene.set.de.results.hammond %>% filter(avg_log2FC > 0.5) %>% rownames(.) 

mg.gene.set.de.results.kracht %>% filter(avg_log2FC < -0.5) %>% rownames(.) 
bam.gene.set.de.results.kracht %>% filter(avg_log2FC > 0.5) %>% rownames(.)

## Differential expression testing of all microglia versus BAMs, comparing expression of MG and BAM gene set genes
hammond.mgs <- subset(hammond, idents = c(1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16))
hammond.bams <- subset(hammond, idents = 2)

Idents(hammond.mgs) <- "MG"
Idents(hammond.bams) <- "BAM"

hammond.mgs.v.bams <- merge(hammond.mgs, hammond.bams)

mg.gene.set.de.results.hammond.mg.v.bam <- FindMarkers(hammond.mgs.v.bams, ident.1 = "MG", ident.2 = "BAM", features = mg.markers.hammond.pruned, 
                                                       logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, min.cells.feature = 0, min.cells.group = 0)

mg.gene.set.de.results.hammond.mg.v.bam <- mg.gene.set.de.results.hammond.mg.v.bam %>% arrange(desc(avg_log2FC))

bam.gene.set.de.results.hammond.mg.v.bam <- FindMarkers(hammond.mgs.v.bams, ident.1 = "MG", ident.2 = "BAM", features = bam.markers.hammond.pruned, 
                                                        logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, min.cells.feature = 0, min.cells.group = 0)

bam.gene.set.de.results.hammond.mg.v.bam <- bam.gene.set.de.results.hammond.mg.v.bam %>% arrange(desc(avg_log2FC))

mg.gene.set.de.results.hammond.mg.v.bam %>% filter(avg_log2FC > 0.5) %>% rownames(.) 
bam.gene.set.de.results.hammond.mg.v.bam %>% filter(avg_log2FC < -0.5) %>% rownames(.) 

kracht.mgs <- subset(kracht, idents = c(1, 2, 3, 4, 5, 6, 7, 10))
kracht.bams <- subset(kracht, idents = 9)

Idents(kracht.mgs) <- "MG"
Idents(kracht.bams) <- "BAM"

kracht.mgs.v.bams <- merge(kracht.mgs, kracht.bams)

mg.gene.set.de.results.kracht.mg.v.bam  <- FindMarkers(kracht.mgs.v.bams, ident.1 = "MG", ident.2 = "BAM", features = mg.markers.orthologs.kracht.pruned, 
                                                       logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, min.cells.feature = 0, min.cells.group = 0)

mg.gene.set.de.results.kracht.mg.v.bam  <- mg.gene.set.de.results.kracht.mg.v.bam  %>% arrange(desc(avg_log2FC))

bam.gene.set.de.results.kracht.mg.v.bam  <- FindMarkers(kracht.mgs.v.bams, ident.1 = "MG", ident.2 = "BAM", features = bam.markers.orthologs.kracht.pruned, 
                                                        logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, min.cells.feature = 0, min.cells.group = 0)

bam.gene.set.de.results.kracht.mg.v.bam  <- bam.gene.set.de.results.kracht.mg.v.bam  %>% arrange(desc(avg_log2FC))

mg.gene.set.de.results.kracht.mg.v.bam %>% filter(avg_log2FC > 0.5) %>% rownames(.) %>% length(.) 
bam.gene.set.de.results.kracht.mg.v.bam %>% filter(avg_log2FC < -0.5) %>% rownames(.) %>% length(.) 

mg.gene.set.de.results.hammond.mg.v.bam$gene <- rownames(mg.gene.set.de.results.hammond.mg.v.bam)
bam.gene.set.de.results.hammond.mg.v.bam$gene <- rownames(bam.gene.set.de.results.hammond.mg.v.bam)
mg.gene.set.de.results.kracht.mg.v.bam$gene <- rownames(mg.gene.set.de.results.kracht.mg.v.bam)
bam.gene.set.de.results.kracht.mg.v.bam$gene <- rownames(bam.gene.set.de.results.kracht.mg.v.bam)

mg.gene.set.de.results.hammond$gene <- rownames(mg.gene.set.de.results.hammond)
bam.gene.set.de.results.hammond$gene <- rownames(bam.gene.set.de.results.hammond)
mg.gene.set.de.results.kracht$gene <- rownames(mg.gene.set.de.results.kracht)
bam.gene.set.de.results.kracht$gene <- rownames(bam.gene.set.de.results.kracht)

# save differential expression test results
write_csv(mg.gene.set.de.results.hammond.mg.v.bam, "hammond_mg_gene_set_allmgs_vs_bam_de_results.csv")
write_csv(bam.gene.set.de.results.hammond.mg.v.bam, "hammond_bam_gene_set_allmgs_vs_bam_de_results.csv")
write_csv(mg.gene.set.de.results.kracht.mg.v.bam, "kracht_mg_gene_set_allmgs_vs_bam_de_results.csv")
write_csv(bam.gene.set.de.results.kracht.mg.v.bam, "kracht_bam_gene_set_allmgs_vs_bam_de_results.csv")

write_csv(mg.gene.set.de.results.hammond, "hammond_mg_gene_set_mrc1pos_v_neg_de_results.csv")
write_csv(bam.gene.set.de.results.hammond, "hammond_bam_gene_set_mrc1pos_v_neg_de_results.csv")
write_csv(mg.gene.set.de.results.kracht, "kracht_mg_gene_set_mrc1pos_v_neg_de_results.csv")
write_csv(bam.gene.set.de.results.kracht, "kracht_bam_gene_set_mrc1pos_v_neg_de_results.csv")

saveRDS(mg.gene.set.de.results.hammond.mg.v.bam, "hammond_mg_gene_set_allmgs_vs_bam_de_results.rds")
saveRDS(bam.gene.set.de.results.hammond.mg.v.bam, "hammond_bam_gene_set_allmgs_vs_bam_de_results.rds")
saveRDS(mg.gene.set.de.results.kracht.mg.v.bam, "kracht_mg_gene_set_allmgs_vs_bam_de_results.rds")
saveRDS(bam.gene.set.de.results.kracht.mg.v.bam, "kracht_bam_gene_set_allmgs_vs_bam_de_results.rds")

saveRDS(mg.gene.set.de.results.hammond, "hammond_mg_gene_set_mrc1pos_v_neg_de_results.rds")
saveRDS(bam.gene.set.de.results.hammond, "hammond_bam_gene_set_mrc1pos_v_neg_de_results.rds")
saveRDS(mg.gene.set.de.results.kracht, "kracht_mg_gene_set_mrc1pos_v_neg_de_results.rds")
saveRDS(bam.gene.set.de.results.kracht, "kracht_bam_gene_set_mrc1pos_v_neg_de_results.rds")


# Additionally, test for differential expression of Spi1 (pu1) between Mrc1+ and Mrc1- microglia in both datasets
hammond.spi1.results <- FindMarkers(hammond.mrc1.pos.v.neg, features = "Spi1", logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, min.cells.feature = 0, min.cells.group = 0, ident.1 = "Mrc1+", ident.2 = "Mrc1-")
kracht.spi1.results <- FindMarkers(kracht.mrc1.pos.v.neg2, features = "SPI1", logfc.threshold = 0, min.pct = 0, min.diff.pct = 0, min.cells.feature = 0, min.cells.group = 0, ident.1 = "MRC1+", ident.2 = "MRC1-")

# save results
saveRDS(hammond.spi1.results, "hammond_spi1_results.rds")
saveRDS(kracht.spi1.results, "kracht_spi1_results.rds")
write_csv(hammond.spi1.results, "hammond_spi1_results.csv")
write_csv(kracht.spi1.results, "kracht_spi1_results.csv")

# save violin plots of Spi1 expression
tiff("spi1_violinplot_hammond.tiff", units = 'in', height = 4, width = 4, res = 300)
VlnPlot(hammond.mrc1.pos.v.neg, features = "Spi1", cols = c("purple", "darkred"), pt.size = 0) + NoLegend() + labs(x = "", y = "Log-normalized expression") + theme(aspect.ratio = 1)
dev.off()

tiff("spi1_violinplot_kracht.tiff", units = 'in', height = 4, width = 4, res = 300)
VlnPlot(kracht.mrc1.pos.v.neg2, features = "SPI1", cols = c("purple", "darkred")) + NoLegend() + labs(x = "", y = "Log-normalized expression") + theme(aspect.ratio = 1)
dev.off()
