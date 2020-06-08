
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

clust_table <- as.data.frame.matrix(table(srt@meta.data$final_clust_review, srt@meta.data$severmod_other))
clust_ito <- unique(rownames(clust_table))


de_list <- list()
age_gend <- list()
age_gend[[1]] <- c("Age", "Gender_num")
covars <- c(age_gend,
            as.list("Gender_num"),
            list("Age"),
            list("No Covars"))

Idents(srt) <- "final_clust_review"


delist <- mclapply(covars, function(cur_covar){
  print(cur_covar)
  cur_covar_name <- cur_covar
  if(cur_covar == "No Covars"){
    cur_covar <- NULL
  }
  i = 1
  de_list <- list()
  cov_one <- "severe-mod"
  cov_one_ind <- which(colnames(clust_table) == cov_one)
  for(x in 1:25){
    if(clust_table[x, cov_one_ind] > 10 & clust_table$Healthy[x] > 10){
      srt_oneclust <- subset(x = srt, idents = clust_ito[x], invert = FALSE)
      Idents(srt_oneclust) <- "severmod_other"
      
      srt_oneclust.markers <- FindMarkers(srt_oneclust, ident.1 = cov_one, ident.2 = "Healthy", min.pct = 0.25,
                                          test = "MAST", latent.vars = cur_covar)
      curclust <- x-1
      srt_oneclust.markers$comp <- paste0(cov_one, "_vs_Healthy_clust_", curclust)
      srt_oneclust.markers$gene <- rownames(srt_oneclust.markers)
      srt_oneclust.markers$latent_vars <- paste0(cur_covar_name, collapse = "_and_")
      de_list[[i]] <- srt_oneclust.markers
      i <- i+1
    }
  } 
  return(de_list)
}, mc.cores = 48)

saveRDS(file = "deg_mast_logfc25_pct25.RDS", delist)


