library(org.Hs.eg.db)
library(pathview)
library(reshape2)
library(ComplexHeatmap)
library(magrittr)

setwd("/path/to/your/files")
source('/helper_setEnrichment.R')

load('bloodtranscriptionalmodules.RData')
too_few_genes = sapply(modGenes, length) < 5
modGenes[too_few_genes] = NULL

all_genes_in_path = lapply(modGenes, function(x) as.data.frame(paste(x, collapse = ",")))
modgenes_df <- rbindlist(all_genes_in_path)
modgenes_df$path_name = names(modGenes)
colnames(modgenes_df)[1] <- "all_genes_in_path"

                           
BTMIDs = names(modGenes)
modgenes_df$path_name[grep("M85",modgenes_df$path_name)]
m150  = modgenes_df$all_genes_in_path[modgenes_df$path_name == "platelet activation & blood coagulation (M199)"]

### world of genes
all_Genes <- rownames(srt@assays$RNA@data)

# load in the degs from markerfinder
genes_interest <- fread("deg_list_wilcox_fc25_pct25.csv")

# subset to severe-mod vs healthy
genes_interest <- genes_interest[grep("severe",genes_interest$comp),]

#take sig
genes_interest <- genes_interest[which(genes_interest$p_val < 0.05),]

# throwing a lot of immunoglobulin genes
                           # likely because all the plasmablasts are exploding
                           # dumping all their RNA on all the cells
                           # plasmablasts have A LOT of RNA
#genes_interest <- genes_interest[-grep("^IG", genes_interest$gene),]

gene_de_list <- split(as.data.table(genes_interest), by = "clust")
genes_up <- lapply(gene_de_list, function(x){
  genes_interest_oneclust_up <- x[which(x$avg_logFC > 0),]
  num_up_genes = nrow(genes_interest_oneclust_up)
  return(num_up_genes)
})

gene_de_list_enrich_up <- lapply(gene_de_list, function(genes_interest_oneclust){
  genes_interest_oneclust_up <- genes_interest_oneclust[which(genes_interest_oneclust$avg_logFC > 0),]
  de_symb <- genes_interest_oneclust_up$gene
  # run the setenrichment
  pathways_sig = setEnrichment(de_symb, modGenes, all_Genes, BTMIDs)
  enrich <- pathways_sig@enrichmentTable
  enrich$clust <- genes_interest_oneclust_up$clust[1]
  enrich$updown <- "UP Genes"
  enrich <- enrich[which(enrich$p.raw < 0.005),]
  # if more than 10 pathways, cut at 10
  # these are org by increasing pval
  if(nrow(enrich) > 10){
    enrich <- enrich[1:10,]
  }
  return(enrich)
})
all_up_enrich <- rbindlist(gene_de_list_enrich_up)

# see pathways per clust
numpath <- sapply(gene_de_list_enrich_up, function(x) nrow(x))
numpath

# same for down
gene_de_list_enrich_down <- lapply(gene_de_list, function(genes_interest_oneclust){
  genes_interest_oneclust_up <- genes_interest_oneclust[which(genes_interest_oneclust$avg_logFC < 0),]
  sig_genes <- genes_interest_oneclust_up$gene
  de_symb <- sig_genes
  pathways_sig = setEnrichment(de_symb, modGenes, all_Genes, BTMIDs)
  enrich <- pathways_sig@enrichmentTable
  enrich$clust <- genes_interest_oneclust_up$clust[1]
  enrich$updown <- "DOWN Genes"
  enrich <- enrich[which(enrich$p.raw < 0.005),]
  if(nrow(enrich) > 10){
    enrich <- enrich[1:10,]
  }
  return(enrich)
})

numpath <- sapply(gene_de_list_enrich_down, function(x) nrow(x))
numpath
all_down_enrich <- rbindlist(gene_de_list_enrich_down)

# pull overrep sig genes
genes_interest_allclust <- rbindlist(gene_de_list)
genes_interest_allclust_up <- genes_interest_allclust[which(genes_interest_allclust$avg_logFC > 0),]
sig_genes <- genes_interest_allclust_up$gene
de_symb <- sig_genes
pathways_sig = setEnrichment(de_symb, modGenes, all_Genes, BTMIDs)
pthgenes = lapply(pathways_sig@sigGenesInSets, function(x) as.data.frame(paste(x, collapse = ",")))
pthgenes_df <- rbindlist(pthgenes)
pthgenes_df$path_name = names(pthgenes)


colnames(pthgenes_df)[1] <- "all_genes_in_path"
scoregenes <- unlist(strsplit(as.character(pthgenes_df$all_genes_in_path[which_genesind]), "[,]"))


#### Rearranging
# I got this code online and modified it
# I will post a link to the source code when I track it down
# makes pretty ring plots
spot.theme_2 <- list(
  theme_classic(),
  theme(axis.ticks.y=element_blank(), axis.text.y=element_text(size = 10)),
  theme(axis.line=element_blank()),
  theme(text = element_text(size = 10)),
  #theme(legend.position = "none"),
  theme(plot.margin = unit(c(10,10,10,10), "mm")),
  #scale_size_continuous(range = c(0.1, 9)),
  scale_x_discrete(position = "top"))


all_up_enrich <- rbindlist(gene_de_list_enrich_up)
# -log(p) would have worked well too
all_up_enrich$p_inv <- log(1/all_up_enrich$p.raw)
min(all_up_enrich$p_inv)
max(all_up_enrich$p_inv)

all_down_enrich <- rbindlist(gene_de_list_enrich_down)
all_down_enrich$p_inv <- log(1/all_down_enrich$p.raw)
min(all_down_enrich$p_inv)
max(all_down_enrich$p_inv)
all_down_enrich$p_inv <- -all_down_enrich$p_inv


all_up_enrich <- rbind(all_up_enrich, all_down_enrich)
all_up_enrich$pct_in <- all_up_enrich$relevant.genes/all_up_enrich$filtered.size
# drop unnamed pathways
all_up_enrich <- all_up_enrich[-grep("TBA", all_up_enrich$set.name),]
unique(all_up_enrich$clust)


all_up_enrich$clust <- factor(all_up_enrich$clust)
myPalette <- colorRampPalette(brewer.pal(3, "RdBu"))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(-10, 10))

all_up_enrich$set.name_clust <- paste0(all_up_enrich$set.name, all_up_enrich$clust)
all_up_enrich <- all_up_enrich[order(all_up_enrich$p.raw, decreasing = F),]
all_up_enrich <- all_up_enrich[!duplicated(all_up_enrich$set.name_clust),]

all_up_enrich$seurat_clusters <- as.factor(all_up_enrich$clust)
all_up_enrich <- all_up_enrich %>%
  left_join(clust_conv)

all_up_enrich_count <- as.data.frame.matrix(table(all_up_enrich$set.name, all_up_enrich$clust))
all_up_enrich_count$rowsum_count <- rowSums(all_up_enrich_count)
all_up_enrich_count <- all_up_enrich_count[which(all_up_enrich_count$rowsum_count > 1),]

all_up_enrich <- all_up_enrich[which(all_up_enrich$set.name %in% rownames(all_up_enrich_count)),]

## cluster_rows
all_up_enrich_sub <- subset(all_up_enrich, select = c("clust", "set.name", "p_inv"))
all_up_enrich_cast <- dcast(all_up_enrich_sub, formula = clust ~ set.name, value.var = "p_inv", fun.aggregate = mean,
                            fill = 0)
rownames(all_up_enrich_cast) <- paste0("Cluster ", all_up_enrich_cast$clust)
all_up_enrich_cast$clust <- NA
all_up_enrich_cast[is.na(all_up_enrich_cast)] <- 0
heaty <- Heatmap(t(all_up_enrich_cast))
row_order_heaty <- row_order(heaty)

all_up_enrich$set.name <- factor(all_up_enrich$set.name, levels = rev(colnames(all_up_enrich_cast)[row_order_heaty]))


ring.spot   <- ggplot(all_up_enrich, aes(clust, set.name, colour = p_inv)) + 
  spot.theme_2 + ggtitle("All Clust Pathways - Severe & Mod vs Healthy") +
  geom_point(colour = "black",     aes(size = 1.2)) +
  geom_point(colour = "white",     aes(size = 1.0)) +
  geom_point(aes(size = 0.81*pct_in, colour = p_inv)) + ylab("") + xlab("") +
  theme(axis.text.x=element_text(size = 10, 
                                 angle = 90, hjust = 0)) +
  
  scale_color_gradientn(colours=c("navyblue", "white", "darkred"),
                        limits = c(-20,20),
                        oob = scales::squish) +
  labs(size = "% Path is DEG", colour = "Log of 1/P val") #scale_colour_gradient(low = "pink", high = "black")
ring.spot

