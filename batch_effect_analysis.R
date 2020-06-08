library(ggplot2)
library(cowplot)
library(ggrastr)
library(matrixStats)
library(parallel)
library(reshape2)
library(ComplexHeatmap)

### To check for batch effect we used 2 methods
# Primarily BEER, to examine the degree of batch effect (there is very little)
# and then harmony to take a look at what happens if we corrected for batch effect


source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

setwd("/path/to/your/files")
srt <- readRDS("final_srt_object.RDS")
Idents(srt) <- "final_clust_review"

meta_uni <- srt@meta.data
meta_uni <- meta_uni[!duplicated(meta_uni$seurat_clusters),]
clust_id <- unique(srt@meta.data$seurat_clusters)

# subset for computational feasibility
srt.small <- subset(srt, downsample = 2000)

# give beer what it needs
expr_dat <- srt.small@assays$RNA@counts
metadat <- srt.small@meta.data
mybeer=BEER(DATA = expr_dat, BATCH = metadat$set, GNUM=30, PCNUM=50, ROUND=1, GN=2000, SEED=1, COMBAT=TRUE, RMG=NULL)   

# plot to look at cors of PCA
# the author of BEER indicate that correlations < 0.7 suggest batch effect
df_plot <- data.frame(cor = mybeer$cor,
                      lcor = mybeer$lcor,
                      pcnum = seq(1:mybeer$PCNUM),
                      clust_id = oneclust_id)

p1 <- ggplot(df_plot, aes(x = cor, y = lcor)) + geom_point() +
  xlab('Rank Correlation of PCs') + ylab('Linear Correlation of PCs') + xlim(0,1) + ylim(0,1) +
  theme_cowplot()

p2 <- ggplot(df_plot, aes(x = pcnum, y = cor)) + geom_point() +
  xlab('PC Number') + ylab('Rank Correlation') +
  theme_cowplot()

pdf("/beer_batch_effect_res.pdf", width = 10, height =5)
grid.arrange(p1, p2, ncol = 2)
dev.off()

