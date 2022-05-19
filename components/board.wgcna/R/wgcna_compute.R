##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

if(0) {
  load("../../../data/example-data.pgx")
  source("../../base/R/pgx-functions.R")
  source("../../base/R/pgx-cluster.R")
  source("../../base/R/gset-fisher.r")
  load("../../../lib/sysdata.rda",verbose=1)    
  pgx <- ngs
  remove(ngs)

  minmodsize=30;power=6;cutheight=0.25;deepsplit=2;ngenes=1000

}


WgcnaCompute <- function(pgx,
                         ngenes = 1000,
                         minmodsize = 30,
                         power = 6,
                         cutheight = 0.25,
                         deepsplit=2
                         )
{

  require(WGCNA)
  library(WGCNA)
                
  WGCNA::enableWGCNAThreads()
        
  X <- as.matrix(pgx$X)
  dim(X)
  X <- X[order(-apply(X,1,sd,na.rm=TRUE)),]
  X <- X[!duplicated(rownames(X)),]
  datExpr <- t(head(X,ngenes))
    
  ## ngenes     <- input$ngenes
  ## minmodsize <- as.integer(input$minmodsize)
  ## power      <- as.numeric(input$power)
  ## cutheight  <- as.numeric(input$cutheight)
  ## deepsplit  <- as.integer(input$deepsplit)
  
  message("[WgcnaCompute] computing WGCNA modules...")        
  net = WGCNA::blockwiseModules(
    datExpr,
    power = power,
    TOMType = "unsigned",
    minModuleSize = minmodsize,
    reassignThreshold = 0,
    mergeCutHeight = cutheight,
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    deepSplit = deepsplit,
    ## saveTOMs = TRUE, saveTOMFileBase = "WCNA.tom",
    verbose = 3)

  names(net)
  table(net$colors)
  
  ## clean up traits matrix
  datTraits <- pgx$samples
  ## no dates please...
  isdate <- apply(datTraits, 2, is.Date)  ## defined in pgx-functions.R
  datTraits <- datTraits[,!isdate,drop=FALSE]
  dim(datTraits)
  
  ## Expand multi-class discrete phenotypes into binary vectors
  ##datTraits1 <- datTraits
  tr.class <- sapply(type.convert(datTraits, as.is=FALSE),class)
  tr.class
  sel1 <- which(tr.class %in% c("factor"))
  sel2 <- which(tr.class %in% c("integer","numeric"))
  tr1 <- datTraits[,0]
  if(length(sel1)) {
    tr1 <- expandPhenoMatrix(datTraits[,sel1,drop=FALSE], drop.ref=FALSE)
  }
  ## keeping numeric phenotypes
  tr2 <- datTraits[,sel2,drop=FALSE]
  datTraits <- cbind(tr1, tr2)
  dim(datTraits)


  labels2rainbow <- function(net) {
    hc <- net$dendrograms[[1]]
    nc <- length(unique(net$colors))
    n <- length(net$colors)
    ii <- hc$order
    col1 <- labels2colors(net$colors)                
    col.rnk <- rank(tapply(1:n,col1[ii],mean))
    new.col <- rainbow(nc)[col.rnk]
    ## new.col <- heat.colors(nc)[col.rnk]
    names(new.col) <- names(col.rnk)
    new.col["grey"] <- "#AAAAAA"
    new.col
    new.col <- new.col[col1]
        names(new.col) <- net$colors
    new.col
  }
    
  ## get colors of eigengene modules
  me.genes <- tapply( names(net$colors), net$colors, list)
  names(me.genes) <- paste0("ME",names(me.genes))        
  color1 <- labels2rainbow(net)
  me.colors <- color1[!duplicated(color1)]
  names(me.colors) <- paste0("ME",names(me.colors))
  me.colors <- me.colors[names(me.genes)]        
  ##progress$inc(0.4,"")

  ## ---------------------------------------------------------------
  message("[WgcnaCompute] >>> calculating WGCNA clustering...")
  #progress$inc(0.1, "computing dim reductions...")    
  X1 <- t(scale(datExpr))
  dissTOM  <- 1 - WGCNA::TOMsimilarityFromExpr(datExpr, power=power)
  rownames(dissTOM) <- colnames(dissTOM) <- colnames(datExpr)            
  clust <- pgx.clusterBigMatrix(dissTOM, methods=c("umap","tsne","pca"), dims=c(2))
  names(clust)
  if("cluster.genes" %in% names(pgx)) {
    clust[['umap2d']] <- pgx$cluster.genes$pos[['umap2d']][colnames(datExpr),]
  }
  
  ## ---------------------------------------------------------------
  message("[WgcnaCompute] >>> calculating WGCNA module enrichments...")
  ##progress$inc(0,"calculating module enrichment...")
  
  gmt <- getGSETS(grep("HALLMARK|GOBP|^C[1-9]",names(iGSETS),value=TRUE))
  gse <- NULL
  ##bg <- unlist(me.genes)
  bg <- toupper(rownames(pgx$X))
  i=1
  for(i in 1:length(me.genes)) {
    gg <- toupper(me.genes[[i]])
    rr <- gset.fisher( gg, gmt, background=bg, fdr=1 )
    rr <- cbind( module = names(me.genes)[i],
      geneset = rownames(rr), rr)
    rr <- rr[order(rr$p.value),,drop=FALSE]
    if(i==1) gse <- rr
    if(i>1) gse <- rbind(gse, rr)
  }
  rownames(gse) <- NULL
    
  ## construct results object
  out <- list(
    datExpr = datExpr,
    datTraits = datTraits,
    net = net,
    gse = gse,
    clust = clust,
    me.genes = me.genes,
    me.colors = me.colors
  )

  return(out)
}
