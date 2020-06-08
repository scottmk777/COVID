
library(Seurat)
library(dplyr)
library(Matrix)
library(data.table)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(ggrastr)
library(matrixStats)
library(parallel)
library(reshape2)
library(ComplexHeatmap)
library(slingshot)
library(qlcMatrix)
set.seed(1)

setwd("/path/to/your/files")
srt <- readRDS("unfiltered_srt_object.RDS")

# the unfiltered srt object has all the output of the below code included
# but should anyone decide to reprocess this is the script
srt = NormalizeData(srt, display.progress = FALSE)
srt = FindVariableFeatures(srt, do.plot = F, display.progress = FALSE)
srt = ScaleData(srt, display.progress = FALSE)
srt <- RunPCA(object = srt, verbose = FALSE)
peb <- ElbowPlot(srt, ndims = 50)
srt <- FindNeighbors(object = srt, dims = 1:25)
srt <- FindClusters(object = srt, resolution = 0.4)
srt <- RunUMAP(object = srt, dims = 1:25)

# plot the umap rep
UMAPPlot(srt, label = T) 

# look at some important qc metrics
p_mito <- ggplot(srt@meta.data, aes(x=seurat_clusters, y = mito)) + 
  geom_jitter(height = 0, width = 0.2) +
  geom_boxplot() +
  theme_cowplot()
p_ncount <- ggplot(srt@meta.data, aes(x=seurat_clusters, y = ncount)) + 
  geom_jitter(height = 0, width = 0.2) +
  geom_boxplot() +
  theme_cowplot() + ggtitle("Count RNA Per Cell")
p_nfeat <- ggplot(srt@meta.data, aes(x=seurat_clusters, y = nfeat)) + 
  geom_jitter(height = 0, width = 0.2) +
  geom_boxplot() +
  theme_cowplot() + ggtitle("Unique Feat Per Cell")


# given the expression and awful qc values of cluster 8
# we are going to get rid of it entirely
# the marker genes for the cluster are also suggestive of dead cells 
srt_clean <- subset(srt, idents = "8", invert = T)
# and also get rid of cells that have more than 25% mito RNA
srt_clean <- subset(srt_clean, subset = mito < 25)

# new plot
UMAPPlot(srt_clean, label = T) 

# remove the old object
rm(srt)

#################################
# Correlations

# sub cluster by monocyte related-clusters
srt_sub_clust <- c("3", "4", "5", "10", "11", "16",  "22")
Idents(srt_clean) <- "seurat_clusters"
sub_srt <- subset(srt_clean, idents = srt_sub_clust)

# pull out the gene data
gene_dat <- sub_srt@assays$RNA@data
dim(gene_dat)
# its awfully big
# so we better be careful correlating

# pull out two genes of interest
# S100A12 we found to be important in the Olink analysis
# IFI27 is a very important viral response gene 
sub_genedat <- gene_dat[which(rownames(gene_dat) %in% c("S100A12", "IFI27")),]

# transform it for easy correlations
gene_dat_t <- t(gene_dat)

# it would have been better to swap zeros for NAs here
# but this matrix is very big, and corSparse (from qlcMatrix) is very fast
corres <- corSparse(gene_dat_t, t(sub_genedat))

# configure for later use
corres <- as.data.frame(corres)
corres$gene <- rownames(gene_dat)

# dont want IFI27 correlations anymore
corres$V2 <- NULL

# keep only correlations that are above .3
# just given the number of points these are likely to be sig
# we will confirm that later
corres_sig <- corres[which(abs(corres$V1) > .3),]
colnames(corres_sig)[1] <- "corval"

# pick out the top contenders that correlate with S100A12 (EN-RAGE)
# this is manual but they are the top 5 most and least correlated genes
IFN_genes_cor <- c("HLA-DPA1", "HLA-DPB1", "RPS19", "HLA-DRA", "CD74", 
                    "S100A8", "S100A9", "VCAN", "CD14", "PLBD1")
# subset the cor list
corres_sig_best <- corres_sig[which(corres_sig$gene %in% IFN_genes_cor),]
corres_sig_best <- corres_sig_best[order(corres_sig_best$corval, decreasing = F),]

# factor for plotting
corres_sig_best$gene <- factor(corres_sig_best$gene, levels = corres_sig_best$gene)

# plot the results
p_cor <- ggplot(data=corres_sig_best, aes(x=gene, y=corval)) +
  geom_bar(stat="identity",  fill = "#478A5E") + coord_flip() + theme_cowplot() +
  ggtitle("Correlation EN-RAGE") + ylab("Pearson Correlation") + xlab("")
p_cor




