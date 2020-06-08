

#### AUTHOR: MICHELE DONATO
## geneSets as a list of geneSets
## gSign is the gene signature
## sigUniverse is the universe of the signature
## an element of geneSets and the signature are character vectors

require(data.table)
require(igraph)
library(network)
library(ndtv)
library(sna)
library(extrafont)

convertIDs <- function( ids, fromKey, toKey, db, ifMultiple=c( "putNA", "useFirst" ) ) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select( 
    db, keys=ids, keytype=fromKey, columns=c(fromKey,toKey) ) )
  if( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]   
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ] }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}

loadfonts()

EnrichmentRes <- setClass("EnrichmentRes",
                          representation = representation(
                            enrichmentTable = 'data.table',
                            sigGenesInSets = 'list'
                          )
)

#' gSign = gene signature in a character vector
#' geneSets = list of character vectors, please use names
#' sigUniverse = gene universe where the signature comes from, character vector
#'
#' test
#' gSign = gene_sig
#' geneSets = keggSets
#' sigUniverse = gene_universe
#' setIDs = keggIDs
setEnrichment <- function(gSign, geneSets, sigUniverse, setIDs = NULL, min_p_genes = 5, min_r_genes = 1){
  
  # min_p_genes = 5
  # min_r_genes = 1
  
  # each set will be filtered; if a gene is not in the reference
  # it will be removed from the set
  geneSetsFiltered = lapply(geneSets, function(x) x[x %in% sigUniverse])
  
  ## S3 type dispatch, hopefully temporary
  if(class(gSign) == 'character'){
    # this is the case where the signature is passed without effect sizes
    
    # make sure that all gSign are in the universe
    gSign = gSign[gSign %in% sigUniverse]
    
    # setSig contains only gene symbols
    setSig = lapply(geneSetsFiltered, function(x) x[x %in% gSign])
    
  }else{
    
    if(is.null(names(gSign)) | !is.numeric(gSign)){
      stop("gSign is not a character vector nor a named vector with numeric
        values.")
    }
    # this is the case where the signature is passed with effect sizes
    
    # filtering gSign by the universe
    gSign = gSign[names(gSign) %in% sigUniverse]
    
    # retrieving effect sizes in each set
    setSig = lapply(geneSetsFiltered, function(x) gSign[names(gSign) %in% x])
    
    # gSign becomes a character vector (legacy)
    gSign = names(gSign)
    
  }
  
  # number of genes in each set that belong to the signature
  sigInSet = sapply(geneSetsFiltered, function(x) sum(x %in% gSign))
  
  # original size of each geneSet
  geneSetSizeOrig = sapply(geneSets, length)
  
  # size of each geneSet (after filtering)
  geneSetSize = sapply(geneSetsFiltered, length)
  
  setEnrichment = data.table('set.name' = names(geneSets),
                             'original.size' = geneSetSizeOrig,
                             'filtered.size' = geneSetSize,
                             'relevant.genes' = sigInSet)
  
  
  if(!is.null(setIDs))
    setEnrichment$set.id = setIDs
  
  set_filter = (setEnrichment$filtered.size >= min_p_genes) & (setEnrichment$relevant.genes >= min_r_genes)
  setEnrichment = setEnrichment[set_filter]
  setSig = setSig[set_filter]
  
  setEnrichment$p.raw = sapply(1:dim(setEnrichment)[1], function(x)
    phyper(
      q = ifelse(setEnrichment[x,relevant.genes] - 1 < 0, 0,
                 setEnrichment[x,relevant.genes] - 1),
      m = setEnrichment[x,filtered.size],
      n = length(sigUniverse) - setEnrichment[x,filtered.size],
      k = length(gSign),
      lower.tail = FALSE
    )
  )
  
  setEnrichment$p.adj = p.adjust(setEnrichment$p.raw, method = 'fdr')
  setOrder <- order(setEnrichment$p.raw)
  setEnrichment <- setEnrichment[setOrder]
  
  return(new(Class = 'EnrichmentRes',
             enrichmentTable = setEnrichment, sigGenesInSets = setSig[setOrder]))
  
}

# get kegg links

enrichmentRes_to_kegg_links_mouse = function(setEnrichmentRes,
                                             n_paths,
                                             gene_sig_es,
                                             organism = 'mmu',
                                             convertGeneIDs = TRUE){
  links = list()
  # i = 1
  for(i in 1:n_paths){
    pid = gsub(organism, '', setEnrichmentRes@enrichmentTable[i, set.id])
    pname = setEnrichmentRes@enrichmentTable[i, set.name]
    genes = gene_sig_es[unlist(gene_sig_es[,1, with = F]) %in% setEnrichmentRes@sigGenesInSets[i][[1]]]
    posGenes = genes[genes[[2]] > 0][[1]]
    negGenes = genes[genes[[2]] < 0][[1]]
    if(length(posGenes) > 0){
      posGenes = convertIDs(posGenes, "SYMBOL", "ENTREZID", org.Mm.eg.db)
    }
    if(length(negGenes) > 0){
      negGenes = convertIDs(negGenes, "SYMBOL", "ENTREZID", org.Mm.eg.db)
    }
    posGenes = posGenes[!is.na(posGenes)]
    negGenes = negGenes[!is.na(negGenes)]
    links[[pname]] = keggLink(posGenes, negGenes, pid, organism)
  }
  return(links)
}

enrichmentRes_to_kegg_links = function(setEnrichmentRes,
                                       n_paths,
                                       gene_sig_es,
                                       organism = 'hsa',
                                       convertGeneIDs = TRUE){
  links = list()
  # i = 1
  for(i in 1:n_paths){
    pid = gsub(organism, '', setEnrichmentRes@enrichmentTable[i, set.id])
    pname = setEnrichmentRes@enrichmentTable[i, set.name]
    genes = gene_sig_es[unlist(gene_sig_es[,1, with = F]) %in% setEnrichmentRes@sigGenesInSets[i][[1]]]
    posGenes = genes[genes[[2]] > 0][[1]]
    negGenes = genes[genes[[2]] < 0][[1]]
    if(length(posGenes) > 0){
      posGenes = convertIDs(posGenes, "SYMBOL", "ENTREZID", org.Hs.eg.db)
    }
    if(length(negGenes) > 0){
      negGenes = convertIDs(negGenes, "SYMBOL", "ENTREZID", org.Hs.eg.db)
    }
    posGenes = posGenes[!is.na(posGenes)]
    negGenes = negGenes[!is.na(negGenes)]
    links[[pname]] = keggLink(posGenes, negGenes, pid, organism)
  }
  return(links)
}



## getters

setGeneric("getEnrichmentTable",
           function(object) {standardGeneric("getEnrichmentTable")}
)

setMethod('getEnrichmentTable', 'EnrichmentRes',
          function(object){
            return(object@enrichmentTable)
          }
)

setGeneric("getSigGenesInSets",
           function(object) {standardGeneric("getSigGenesInSets")}
)

setMethod('getSigGenesInSets', 'EnrichmentRes',
          function(object){
            return(object@sigGenesInSets)
          }
)

## show

setMethod("show", "EnrichmentRes",
          function(object){
            cat("*** Top 10 sets (or max not NA), ranked by p-value ***\n")
            not_na_paths = sum(!is.na(object@enrichmentTable$set.name))
            if(not_na_paths >= 10){
              print(object@enrichmentTable[1:10,])
            }else{
              print(object@enrichmentTable[1:not_na_paths,])
            }
          }
)

## graph helpers

setGeneric("getNetObj",
           function(eobject, sets, nSets) {standardGeneric("getNetObj")}
)

# eobject = btmEnrichmentES[[3]]
# sets = modGenes

setMethod("getNetObj", 'EnrichmentRes',
          function(eobject, sets, nSets){
            # retrieve significant genes in sets
            sigInSets <- getSigGenesInSets(eobject)
            eTable <- getEnrichmentTable(eobject)
            # order the sets in the same order as significant genes in sets list
            sets <- sets[names(sigInSets)]
            # create nodes of the graph
            enodes <- lapply(1:length(sets),
                             function(x) data.table(id = names(sets)[x],
                                                    setSize = length(sets[[x]]),
                                                    de = length(sigInSets[[x]]),
                                                    averageES = mean(sigInSets[[x]]),
                                                    padj = eTable[x,p.adj]
                             )
            )
            enodes = Reduce(rbind, enodes)
            enodes = enodes[1:nSets,]
            # all the pairs, triangular, no diagonal
            allPairs <- combn(nSets, 2)
            eedges = apply(t(allPairs), 1, findEdges,
                           mySets = sets, sigSets = sigInSets)
            eedges = data.table(Reduce(rbind, eedges))
            
            # remove empty intersections
            
            eedges = eedges[sizeInt > 0,]
            edgees = eedges[, from := as.character(from)]
            edgees = eedges[, to := as.character(to)]
            
            eNet = graph.data.frame(eedges, directed = FALSE, vertices = enodes)
            return(eNet)
          }
)

# inds = allPairs[,1]
findEdges <- function(inds, mySets, sigSets){
  ind1 = inds[1]
  ind2 = inds[2]
  from = names(mySets)[ind1]
  to = names(mySets)[ind2]
  sizeInt = length(intersect(mySets[[ind1]], mySets[[ind2]]))
  deInt = length(intersect(sigSets[[ind1]], sigSets[[ind2]]))
  findEdges = data.frame(
    to, from, sizeInt, deInt)
}

setGeneric("ePlot",
           function(eRes, esets, nSets = 10, ...) {standardGeneric("ePlot")}
)


# eRes = btmEnrichmentES[[1]]
# esets = modGenes
# nSets = 10

setMethod('ePlot',
          signature = c('EnrichmentRes', 'list'),
          
          function(eRes, esets, nSets = 10, minDEInt = 1, nodeSizeCex = 50,
                   edgeSizeCex = 5, ...){
            
            eNet = getNetObj(eRes, sets = esets, nSets = nSets)
            
            # removing edges without any DE shared
            eNet <- delete_edges(eNet, which(E(eNet)$deInt < 1)-1)
            # set node colors
            V(eNet)$color = nodeColors(V(eNet)$averageES)
            V(eNet)$size = V(eNet)$de / V(eNet)$setSize * nodeSizeCex
            
            # set edge width
            
            E(eNet)$width = E(eNet)$deInt / E(eNet)$sizeInt * 5
            
            plotlayout = layout.kamada.kawai(eNet)
            
            #     vertexLabels = paste(V(eNet)$name, '\n', V(eNet)$setSize, '/', V(eNet)$de, sep = '')
            #     edgeLabels = paste(E(eNet)$sizeInt, '/', E(eNet)$deInt, sep = "")
            
            vertexLabels = paste(V(eNet)$name, '\n', V(eNet)$de, '/', V(eNet)$setSize, sep = '')
            edgeLabels = paste(E(eNet)$deInt, '/', E(eNet)$sizeInt, sep = "")
            
            # edgeLabels[grep('/0', edgeLabels)] = ""
            
            plot(eNet,
                 vertex.label = vertexLabels,
                 vertex.label.family = "Helvetica",
                 vertex.label.dist = .38, # distance label-node
                 vertex.label.degree = -pi/3, # position of label in radians
                 vertex.label.cex = .7,
                 vertex.label.color = 'black',
                 edge.label = edgeLabels,
                 edge.label.family = "Helvetica",
                 edge.label.cex = 0.6,
                 layout = layout.star,
                 ...
            )
          }
)

#nodesES = V(eNet)$averageES
nodeColors <- function(nodesES, scaleLength = 64){
  # palettes for negative and positives
  negColScale = colorRampPalette(c('blue', 'white'))
  posColScale = colorRampPalette(c('white', 'red'))
  
  #maximum deviation from zero
  maxabs = max(abs(nodesES))
  
  negCols = negColScale(scaleLength) # neg colors
  posCols = posColScale(scaleLength-1) # pos colors
  
  # scalevals tells me where the colors end up being
  # the function in the sapply tells me the closest color to the effect size
  scaleVals = seq(-maxabs, maxabs, length.out = 2 * scaleLength -1)
  myScale = c(negCols,posCols)
  nodeColors = sapply(nodesES, function(x) myScale[which.min(abs(scaleVals-x))])
}



