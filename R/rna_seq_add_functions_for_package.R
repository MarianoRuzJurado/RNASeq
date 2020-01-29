###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###
###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###
###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###

### alternative run_DESeq2 function
### return list of data.frame with DESeq2 results and corresponding dds object

myRun_DESeq2Function <- function (countDF, targets, cmp, independent = FALSE)
{
  if (class(cmp) != "matrix" & length(cmp) == 2)
    cmp <- t(as.matrix(cmp))
  samples <- as.character(targets$Factor)
  names(samples) <- paste(as.character(targets$SampleName),
                          "", sep = "")
  countDF <- countDF[, names(samples)]
  countDF[is.na(countDF)] <- 0
  deseqDF <- data.frame(row.names = rownames(countDF))
  if (independent == TRUE) {
    loopv <- seq(along = cmp[, 1])
  }
  else {
    loopv <- 1
  }
  for (j in loopv) {
    if (independent == TRUE) {
      subset <- samples[samples %in% cmp[j, ]]
      countDFsub <- countDF[, names(subset)]
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(countDFsub),
                                            colData = data.frame(condition = subset), design = ~condition)
      mycmp <- cmp[j, , drop = FALSE]
    }
    else {
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(countDF),
                                            colData = data.frame(condition = samples), design = ~condition)
      mycmp <- cmp
    }
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    for (i in seq(along = mycmp[, 1])) {
      res <- DESeq2::results(dds, contrast = c("condition",
                                               mycmp[i, ]))
      res[is.na(res[, "padj"]), "padj"] <- 1
      res[is.na(res[, "log2FoldChange"]), "log2FoldChange"] <- 0
      deg <- as.data.frame(res)
      colnames(deg)[colnames(deg) %in% c("log2FoldChange",
                                         "padj")] <- c("logFC", "FDR")
      colnames(deg) <- paste(paste(mycmp[i, ], collapse = "-"),
                             colnames(deg), sep = "_")
      deseqDF <- cbind(deseqDF, deg[rownames(deseqDF),
                                    ])
      mylist <- list(DESeq2_results=deseqDF,DESeq2_Object=dds)
    }
  }
  return(mylist)
}

###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###
###  Parse Deseq2 Data Frame & Object to build a list with all the results                          ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###
###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###

parseDESeq2Res <- function(cmp,collectResults=collectResults,baseDir=baseDir,myRes=myRes,padj =padjDEgenes ,lfc=lFCDEgenes) {
  for (i in 1:nrow(cmp[[1]])) {
   # print(i)
    ctrst <- paste(cmp[[1]][i,1],cmp[[1]][i,2],sep="-")
    #print(ctrst)
   # print(dim(ctrst))
    contrastFolder <- paste(baseDir,"/results/",ctrst,sep="")
    ### add folder to results

    #PathListcontrastFolder[[i]] <- list(contrastFolder)
    if(!dir.exists(contrastFolder)) {
      dir.create(contrastFolder)

    } else print(paste("Directory",ctrst,"exists",sep=" "))

    DESeqFolder <- paste(contrastFolder,"/DESeq2",sep="")

    if(!dir.exists(DESeqFolder)) {
      dir.create(DESeqFolder)
    } else print(paste("Directory",DESeqFolder,"exists",sep=" "))

    #sRes <- myRes$DESeq2_results[,grep(ctrst,colnames(Run_DESeq2frame)),]
    sRes <- myRes$DESeq2_results[,grep(ctrst,colnames(myRes$DESeq2_results)),]
    sRes <- data.frame(sRes,counts(myRes$DESeq2_Object, normalized=TRUE))

    collectResults[["DESeq2"]][["contrasts"]][[ctrst]]<-list()
    collectResults[["DESeq2"]][["contrasts"]][[ctrst]][["folder"]] <- contrastFolder
    collectResults[["DESeq2"]][["contrasts"]][[ctrst]][["DF"]] <- sRes

    sR <- sRes[sRes[,6] < padj & !is.na(sRes[,6]) & abs(sRes[,2])>=lfc ,]

    tL <-list()
    tL[["diffGenes"]]<-rownames(sR)
    tL[["upGenes"]]<-rownames(sR)[sR[,2]>0]
    tL[["downGenes"]]<-rownames(sR)[sR[,2]<0]
    collectResults[["DESeq2"]][["contrasts"]][[ctrst]][["diffGenes"]]<-tL
  }
  return(collectResults)
}
###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###
###  Parse EdgeR Data Frame and count object into results list                                      ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###
###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###

parseEdgeRRes <- function(cmp,myRes,collectResults=collectResults,padj =padjDEgenes ,lfc=lFCDEgenes, baseDir=baseDir) {
  for (i in 1:nrow(cmp[[1]])) {
     print(i)
     ctrst <- paste(cmp[[1]][i,1],cmp[[1]][i,2],sep="-")
     print(ctrst)
  #   #print(dim(ctrst))
    contrastFolder <- paste(baseDir,"/results/",ctrst,sep="")
    print(contrastFolder)
     ### add folder to results
  #
  #   PathListcontrastFolder[[i]] <- list(contrastFolder)
     if(!dir.exists(contrastFolder)) {
       dir.create(contrastFolder)

     } else print(paste("Directory",ctrst,"exists",sep=" "))

     EdgeRFolder <- paste(contrastFolder,"/edgeR",sep="")
  #
     if(!dir.exists(EdgeRFolder)) {
       dir.create(EdgeRFolder)
       print(paste("created folder",EdgeRFolder))
    } else print(paste("Directory",EdgeRFolder,"exists",sep=" "))

     sRes <- myRes[[1]][,grep(ctrst,colnames( myRes[[1]] ))]
     print(head(sRes))
     print(head(myRes[[2]]$fitted.values))
     sRes <- merge(sRes,myRes[[2]]$fitted.values,by="row.names")
     rownames(sRes) <- sRes[,1]
     sRes<-sRes[,-1]
     print(head(sRes))
     #stop()
     collectResults[["edgeR"]][["contrasts"]][[ctrst]]<-list()
     collectResults[["edgeR"]][["contrasts"]][[ctrst]][["folder"]] <- contrastFolder
     collectResults[["edgeR"]][["contrasts"]][[ctrst]][["DF"]] <- sRes

     sR <- sRes[sRes[,5] < padj & !is.na(sRes[,5]) & abs(sRes[,1])>=lfc ,]

     tL <-list()
     tL[["diffGenes"]]<-rownames(sR)
     tL[["upGenes"]]<-rownames(sR)[sR[,1]>0]
     tL[["downGenes"]]<-rownames(sR)[sR[,1]<0]
     collectResults[["edgeR"]][["contrasts"]][[ctrst]][["diffGenes"]]<-tL
  }
  return(collectResults)

}

###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###
###  hacked run_edgeR, gives also a normalized count object   #  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###
###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###

myrun_edgeR <- function (countDF, targets, cmp, independent = TRUE, paired = NULL,
                     mdsplot = "")
{
  if (class(cmp) != "matrix" & length(cmp) == 2)
    cmp <- t(as.matrix(cmp))
  samples <- as.character(targets$Factor)
  names(samples) <- paste(as.character(targets$SampleName),
                          "", sep = "")
  countDF <- countDF[, names(samples)]
  countDF[is.na(countDF)] <- 0
  edgeDF <- data.frame(row.names = rownames(countDF))
  group <- as.character(samples)
  if (independent == TRUE) {
    loopv <- seq(along = cmp[, 1])
  }
  else {
    loopv <- 1
  }
  for (j in loopv) {
    y <- DGEList(counts = countDF, group = group)
    if (independent == TRUE) {
      subset <- samples[samples %in% cmp[j, ]]
      y <- y[, names(subset)]
      y$samples$group <- factor(as.character(y$samples$group))
    }
    keep <- rowSums(cpm(y) > 1) >= 2
    y <- y[keep, ]
    y <- calcNormFactors(y)
    if (length(paired) == 0) {
      design <- model.matrix(~0 + y$samples$group, data = y$samples)
      colnames(design) <- levels(y$samples$group)
    }
    else {
      if (length(paired) > 0 & independent == FALSE)
        stop("When providing values under 'paired' also set independent=TRUE")
      Subject <- factor(paired[samples %in% cmp[j, ]])
      Treat <- y$samples$group
      design <- model.matrix(~Subject + Treat)
      levels(design) <- levels(y$samples$group)
    }
    y <- estimateGLMCommonDisp(y, design, verbose = TRUE)
    y <- estimateGLMTrendedDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    fit <- glmFit(y, design)
    if (independent == TRUE) {
      mycomp <- paste(cmp[j, 1], cmp[j, 2], sep = "-")
    }
    else {
      mycomp <- paste(cmp[, 1], cmp[, 2], sep = "-")
    }
    if (length(paired) == 0)
      contrasts <- makeContrasts(contrasts = mycomp, levels = design)
    for (i in seq(along = mycomp)) {
      if (length(paired) == 0) {
        lrt <- glmLRT(fit, contrast = contrasts[, i])
      }
      else {
        lrt <- glmLRT(fit)
      }
      deg <- as.data.frame(topTags(lrt, n = length(rownames(y))))
      colnames(deg) <- paste(paste(mycomp[i], collapse = "_"),
                             colnames(deg), sep = "_")
      edgeDF <- cbind(edgeDF, deg[rownames(edgeDF), ])
    }
    if (nchar(mdsplot) > 0) {
      pdf(paste("./results/sample_MDS_", paste(unique(subset),
                                               collapse = "-"), ".pdf", sep = ""))
      plotMDS(y)
      dev.off()
    }
  }
  return( list(edgeR_results=edgeDF,edgeR_Object=lrt))
}

###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###
###  function for heatmap Z-score differential & Individual log2f versus control  #  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###
###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###

heatmapbuilder <- function(collectResults){
for (i in names(collectResults )) {
  DETool=names(collectResults[i])
  print(DETool)
  tmp <- collectResults
  #print(tmp)
  CommonResults=collectResults[[DETool]][[1]][[3]]
  #print(CommonResults)

  ### Z-score differential genes
  myMat <- CommonResults+1
  #print(myMat)
  #print(length(tmp))
  if (!length(tmp)==0) {
    for (j in names(tmp[[DETool]]$contrasts)) {
      mynames=names(tmp[[DETool]]$contrasts[j])
      print(mynames)
      sMat <- myMat[tmp[[DETool]]$contrasts[[mynames]]$diffGenes$diffGenes, ]
      ncol <- 20
      mycol <- colorpanel(ncol,"blue","white","red")
      max(abs(t(scale(t(sMat)))))->mx

      mybreaks <- seq(-mx,mx,by=2*mx/ncol)
      outFile <-paste(tmp[[DETool]]$contrasts[[mynames]]$folder,"/",DETool,"/Z-score_diff_genes_",mynames,".pdf",sep="")
      pdf(outFile)
      heatmap.2(t(scale(t(sMat))),
                cexRow=0.1,
                trace=c("none"),
                density.info=c("none"),
                col=mycol,
                breaks=mybreaks,
                Colv="none",
                dendrogram = c("row")
      )
      dev.off()
      outFile<-""
      outFile2 <- paste(tmp[[DETool]]$contrasts[[mynames]]$folder,"/",DETool,"/diff_genes_",mynames,".xls",sep="")
      write.table(tmp[[DETool]]$contrasts[[j]]$DF[tmp[[DETool]]$contrasts[[mynames]]$diffGenes$diffGenes,],file=outFile2,sep="\t",quote=F, col.names = NA) ### use pretty or format( ... , digits=2)
      outFile2 <-""

      ### heatmap with log2Foldchanges individual treatments vs mean control

      controlCondition <- strsplit(j,"-")[[1]][2]
      myTreatmeantCondition <- strsplit(j,"-")[[1]][1]
      mMn <- apply(sMat[,grep(controlCondition,colnames(sMat))],1,mean)
      #print(mMn)
      log2fc <- log2(sMat[,grep(myTreatmeantCondition,colnames(sMat))]/mMn)
      #print(log2fc)
      rL <- list(log2fc)
      #print(rL)
      do.call(cbind, rL) -> l2fcMat
      max(abs(l2fcMat)) ->mx
      mybreaks <- seq(-mx,mx,by=2*mx/ncol)
      pdf(paste(tmp[[DETool]]$contrasts[[j]]$folder,"/",DETool,"/individual_Log2FC_versus_controls_",mynames,".pdf",sep=""))
      heatmap.2(l2fcMat,
                cexRow=0.1,
                trace=c("none"),
                density.info=c("none"),
                col=mycol,
                breaks=mybreaks,
                Colv="none",
                dendrogram = c("row"),
                main=paste(mynames)
      )
      dev.off()

    }
  }
}
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ClusterProfiler GO over-representation test with up and down diff genes     ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


GOoverrepresentation <- function(collectResults,OrgDb,pvalueCutoff,pAdjustMethod,qvalueCutoff,minGSSize,maxGSSize) {

  for (i in names(collectResults )) {

 #  i=names(collectResults[1])
  DETool=names(collectResults[i])
  print(DETool)
  # tmp <- collectResults
  rawcounts = collectResults[[DETool]][[1]][[4]]

  for (j in names(collectResults[[DETool]]$contrast)) {
    # j=names(collectResults[[DETool]]$contrast[1])
   #  print(j)
    ctrst <- names(collectResults[[DETool]]$contrasts[j])
    if (length(collectResults[[DETool]]$contrasts[[ctrst]][["enrichGO_BP"]]) == 0 || myOverwrite == TRUE) {
    print(ctrst)

    collectResults[[DETool]]$contrasts[[ctrst]][["enrichGO_BP"]] <- list()
    collectResults[[DETool]]$contrasts[[ctrst]][["enrichGO_MF"]] <- list()
    collectResults[[DETool]]$contrasts[[ctrst]][["enrichGO_CC"]] <- list()
    as.vector(rownames(rawcounts)) -> egallgenes # universe
    egallgenes = bitr(egallgenes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
   # head(egallgenes)
    as.vector(collectResults[[DETool]]$contrasts[[ctrst]]$diffGenes$upGenes) -> egupgenes
    egupgenes = bitr(egupgenes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
    # head(egupgenes)

    for (k in c("BP", "MF", "CC")) {

      ego <- enrichGO(gene = egupgenes[,2],
                      OrgDb = OrgDb,
                      keyType = "ENTREZID",
                      ont = k,
                      universe = egallgenes[,2],
                      pvalueCutoff = pvalueCutoff,
                      pAdjustMethod = pAdjustMethod,
                      qvalueCutoff = qvalueCutoff,
                      minGSSize = minGSSize,
                      maxGSSize = maxGSSize,
                      readable = TRUE,
                      pool = FALSE)
      # head(ego)

      collectResults[[DETool]]$contrasts[[ctrst]][[paste("enrichGO_",k, sep = "")]] <- ego
      outfile <- ""
      outfile <- paste(collectResults[[DETool]]$contrasts[[ctrst]]$folder, "/", DETool, "/over_representation_",k,".xls",sep = "")
      write.table(collectResults[[DETool]]$contrasts[[ctrst]]$enrichGO_BP, file = outfile,sep="\t",quote = F, col.names = NA)

      }
     } else{print("Already Done!")}
    }
  }
  return(collectResults)
}
