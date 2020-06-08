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
set.seed(1)

setwd("/path/to/your/files")
srt <- readRDS("final_srt_object.RDS")

srt@meta.data$severmod_other <- ifelse(srt@meta.data$Resp %in% c("Severe", "Moderate"), "severe-mod", srt@meta.data$Resp)
srt@meta.data$severmod_other <- ifelse(srt@meta.data$disease %in% c("FluA", "RSV"), "Flu-RSV", srt@meta.data$severmod_other )
table(srt@meta.data$severmod_other)

clust_table <- as.data.frame.matrix(table(srt@meta.data$final_clust_review, srt@meta.data$severmod_other))
clust_ito <- unique(rownames(clust_table))


de_list <- list()
cov_groups <- as.list(c("Conv", "Flu-RSV", "severe-mod"))
cov_one <- cov_pts[1]

Idents(srt) <- "final_clust_review"

delist <- mclapply(cov_groups, function(cov_one){
  print(cov_one)
  i = 1
  de_list <- list()
  cov_one_ind <- which(colnames(clust_table) == cov_one)
  for(x in 1:length(clust_ito)){
    if(clust_table[x, cov_one_ind] > 10 & clust_table$Healthy[x] > 10){
      srt_oneclust <- subset(x = srt, idents = clust_ito[x], invert = FALSE)
      Idents(srt_oneclust) <- "severmod_other"
      srt_oneclust.markers <- FindMarkers(srt_oneclust, ident.1 = cov_one, ident.2 = "Healthy", min.pct = 0.25)
      curclust <- x-1
      srt_oneclust.markers$comp <- paste0(cov_one, "_vs_Healthy_clust_", curclust)
      srt_oneclust.markers$gene <- rownames(srt_oneclust.markers)
      de_list[[i]] <- srt_oneclust.markers
      i <- i+1
    }
  } 
  return(de_list)
}, mc.cores = 48)

saveRDS(file = "deg_wilcox_logfc25_pct25.RDS", delist)

### look at output
delist_df_l <- lapply(delist, function(x) rbindlist(x))
delist_df <- rbindlist(delist_df_l)
colnames(delist_df)[3:4] <- c("Case Group (pct cell expr)", "Healthy (pct cells expr)")

delist_df$clust <- sapply(strsplit(delist_df$comp, "_"), "[", 5)
delist_df$case_group <- sapply(strsplit(delist_df$comp, "_"), "[", 1)

fwrite(file = "deg_list_wilcox_fc25_pct25.csv", delist_df)


