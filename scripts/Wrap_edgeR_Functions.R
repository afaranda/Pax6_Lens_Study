########################         HEADER BLOCK      ###########################
#  File:    Wrap_edgeR_Functions.R                                           #
#  Purpose: Define wrapers to streamline automated edgeR analysis            #
#  Created: Sep 11, 2021                                                     #
#  Author:  Adam Faranda                                                     #
#                                                                            #
##############################################################################

#########      Load Libraries and import expression data          ############
library(edgeR)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(pheatmap)
wd<-getwd()
########  Wrapper function to fit a QLF model to a design matrix  ############
processByDesign <- function(y, design, rob=T, norm="TMM"){
  if(norm == "TMM"){
    y <- y[filterByExpr(y, design),,keep.lib.sizes=F]
    if(norm !="NONE"){y <- calcNormFactors(y,method=norm)}
    y <- estimateDisp(y, design, robust = rob)
    
    fpkm<-rpkm(y, normalized.lib.sizes = T, gene.length = "eu_length")
    for(g in unique(y$samples$group)){
      s<-row.names(y$samples[y$samples$group == g,])
      y$genes[paste0(g,"_Avg_FPKM")]<-apply(fpkm[,s],1,mean)
    }
    
    fit <- glmQLFit(y, design, robust = rob)
    return(list(dge=y, fit=fit, design=design))
    
  } else if(norm == "VOOM-TMM"){
    print(levels(y$samples$group))
    design<-model.matrix(~group, y$samples)
    colnames(design) <- gsub("(\\(|\\))", "", colnames(design))
    colnames(design) <-gsub("group","",colnames(design))
    y <- y[filterByExpr(y, design),,keep.lib.sizes=F]
    y <- calcNormFactors(y)
    
    # NOTE -- FPKM Calculation Fails to account for VOOM precision Weights !!!
    # THIS SHOULD BE FIXED if VOOM is intended for use !!!
    fpkm<-rpkm(y, normalized.lib.sizes = T, gene.length = "eu_length")
    for(g in unique(y$samples$group)){
      s<-row.names(y$samples[y$samples$group == g,])
      y$genes[paste0(g,"_Avg_FPKM")]<-apply(fpkm[,s],1,mean)
    }
    
    v <- voom(y, design=design)
    return(list(v=v, design=design))
  } else if(norm == "VOOM-QNT"){
    print(levels(y$samples$group))
    design<-model.matrix(~group, y$samples)
    colnames(design) <- gsub("(\\(|\\))", "", colnames(design))
    colnames(design) <-gsub("group","",colnames(design))
    y <- y[filterByExpr(y, design),,keep.lib.sizes=F]
    
    # NOTE -- FPKM Calculation Fails to account for VOOM precision Weights !!!
    # THIS SHOULD BE FIXED if VOOM is intended for use !!!
    fpkm<-rpkm(y, normalized.lib.sizes = T, gene.length = "eu_length")
    for(g in unique(y$samples$group)){
      s<-row.names(y$samples[y$samples$group == g,])
      y$genes[paste0(g,"_Avg_FPKM")]<-apply(fpkm[,s],1,mean)
    }
    
    v <- voom(y, design=design, normalize.method ="quantile" )
    return(list(v=v, design=design))
  }
}
############################ Plot Diagnostics ################################
## Generate diagnostic plots. 
diagnostic_plots <- function(
  dge=dge, color_attrib="group", shape_attrib=NULL, 
  respath="LIRTS_DEG_Analysis_results", prefix="DBI_Wildtype"
){
  # Colorblind friendly pallatte from 
  # https://bconnelly.net/
  
  # colors=c(
  #   "#000000", "#E69F00", "#56B4E9", "#009E73",
  #   "#0072B2", "#D55E00", "#CC79A7"
  # )
  colors=c(
    RColorBrewer::brewer.pal(9, name="RdYlBu")[1:3],
    RColorBrewer::brewer.pal(9, name="RdYlBu")[7:9]
  )[c(1,6,2,5,3,4)]
  
  pdf(
    paste(
      respath,"/",prefix,"_BCV_Plot.pdf", sep=""
    ), height=6, width=6
  )  
  plotBCV(dge)                                              # BCV Plot
  dev.off()
  
  # Plot sample projections in first two Principal Components
  pca <-prcomp(t(cpm(dge, log=T)), scale=T)
  fn <- paste(
    respath,"/",prefix,"_PCA_Plot.png", sep=""
  )
  
  if(is.null(color_attrib)){
    ggsave(
      fn,
      autoplot(
        pca, data=dge$samples,size=3
      ), height =4, width=6
    )
  } else if(is.null(shape_attrib)){
    ggsave(
      fn,
      autoplot(
        pca, data=dge$samples,
        colour=color_attrib, size=3
      ) + scale_color_manual(values=colors), 
      height =4, width=6
    )
  } else {
    ggsave(
      fn,
      autoplot(
        pca, data=dge$samples,
        colour=color_attrib, shape=shape_attrib, size=3
      ) + scale_color_manual(values=colors), 
      height =4, width=6
    )
  }
  
  ## Plot sample correlation matrix (using all features)
  
  fpkm<-as.data.frame(                         ## store matrix of FPKM values
    rpkm(
      dge, gene.length="eu_length", 
      log=T, normalized.lib.sizes = T
    )
  )
  colnames(fpkm) <- dge$samples$label
  
  cm <- cor(fpkm, method = "spearman")
  annot <- dge$samples[,c("label", color_attrib, shape_attrib)]
  annot <- annot %>%
    tibble::remove_rownames()%>%
    tibble::column_to_rownames("label")
  
  annot_colors=list(
    color_attrib=colors[1:length(levels(dge$samples[,color_attrib]))],
    shape_attrib=colors[1:length(levels(dge$samples[,shape_attrib]))]
  )
  names(annot_colors$color_attrib)=levels(annot[,color_attrib])
  names(annot_colors$shape_attrib)=levels(annot[,shape_attrib])
  names(annot_colors) <- c(color_attrib, shape_attrib)
  
  pheatmap(
    cm, annotation_col =annot,
    annotation_colors = annot_colors,
    filename = paste0(
      "results/",
      prefix, "_cormat.png"
    ),height = 6, width=8
  )
  
  
  
}

######## Split and Process DGE List based on a set of sample groups ##########
subsetDGEListByGroups<-function(y, groups=c("GR1", "GR2"), norm="TMM"){
  
  # Get samples associated with two groups
  group_subset<-list()
  for(gr in groups){
    group_subset[[gr]] <- row.names(
      y$samples %>% dplyr::filter(group == gr)
    )
  }
  print(group_subset)
  # Extract Subset DGEList 
  y<-y[, unlist(group_subset)]
  
  # Relevel Factors as needed
  for(covariate in colnames(y$samples)){
    if(is.factor(y$samples[covariate])){
      y$samples[covariate] <-droplevels(y$samples[covariate])
    }
  }
  
  # Fix the first group as the reference level / Intercept
  y$samples$group <- relevel(y$samples$group, groups[1])
  
  if(norm == "TMM"){
    design<-model.matrix(~group, y$samples)
    colnames(design) <-gsub("group","",colnames(design))
    y <- y[filterByExpr(y, design),,keep.lib.sizes=F]
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design, robust=T)
    y$genes$com_disp <- y$common.dispersion
    y$genes$trend_disp <- y$trended.dispersion
    y$genes$tag_disp <- y$tagwise.dispersion
    
    fpkm<-rpkm(y, normalized.lib.sizes = T, gene.length = "eu_length")
    for(g in unique(y$samples$group)){
      s<-row.names(y$samples[y$samples$group == g,])
      y$genes[paste0(g,"_Avg_FPKM")]<-apply(fpkm[,s],1,mean)
    }
    
    fit <- glmQLFit(y, design, robust = T)
    return(
      list(
        dge=y, fit=fit, design=design
      )
    )
  } else if (norm == "VOOM-TMM") {
    print(levels(y$samples$group))
    design<-model.matrix(~group, y$samples)
    colnames(design) <- gsub("(\\(|\\))", "", colnames(design))
    colnames(design) <-gsub("group","",colnames(design))
    y <- y[filterByExpr(y, design),,keep.lib.sizes=F]
    y <- calcNormFactors(y)
    
    # NOTE -- FPKM Calculation Fails to account for VOOM precision Weights !!!
    # THIS SHOULD BE FIXED if VOOM is intended for use !!!
    fpkm<-rpkm(y, normalized.lib.sizes = T, gene.length = "eu_length")
    for(g in unique(y$samples$group)){
      s<-row.names(y$samples[y$samples$group == g,])
      y$genes[paste0(g,"_Avg_FPKM")]<-apply(fpkm[,s],1,mean)
    }

    v <- voom(y, design=design)
    return(list(v=v, design=design))
  } else if (norm == "RUVs"){
    
    design<-model.matrix(~group, y$samples)
    colnames(design) <-gsub("group","",colnames(design))
    y <- y[filterByExpr(y, design),,keep.lib.sizes=F]
    controls <- row.names(y)
    
    # Assemble "Samples" matrix required by RUVs
    differences <- matrix(-1,
        nrow = length(levels(y$samples$group)), 
        ncol = max(table(y$samples$group))
    )
    for(gri in 1:length(levels(y$samples$group))){
      level <- levels(y$samples$group)[gri]
      
      s <- which(y$samples$group == level)
      for(j in 1:length(s)){
        differences[gri, j] <- s[j]
      }
    }
    print("Difference Sets for RUVs")
    print(differences)
    
    # Estimate nuisance parameter W_1
    seq <- newSeqExpressionSet(counts = y$counts)
    ruvs <- RUVs(
      seq, cIdx = row.names(y),
      k =1, scIdx=differences
    )
    y$samples$W_1 <- pData(ruvs)[row.names(y$samples), "W_1"]
    design<-model.matrix(~group + W_1, y$samples)
    colnames(design) <-gsub("group","",colnames(design))
    
    y <- y[filterByExpr(y, design),,keep.lib.sizes=F]
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design, robust=T)
    
    ## Add Dispersion estimates to genes table
    y$genes$com_disp <- y$common.dispersion
    y$genes$trend_disp <- y$trended.dispersion
    y$genes$tag_disp <- y$tagwise.dispersion
    
    # NOTE -- FPKM Calculation Fails to account for RUV adjustment !!!
    # THIS may need to be FIXED if RUV is intended for use !!!
    fpkm<-rpkm(y, normalized.lib.sizes = T, gene.length = "eu_length")
    for(g in unique(y$samples$group)){
      s<-row.names(y$samples[y$samples$group == g,])
      y$genes[paste0(g,"_Avg_FPKM")]<-apply(fpkm[,s],1,mean)
    }
    
    fit <- glmQLFit(y, design, robust = T)
    return(
      list(
        dge=y, fit=fit, ruvs=ruvs, design=design
      )
    )
    
  } else if (norm == "RUVr"){
    
    # Fit a typical edgeR model
    z <- y
    design<-model.matrix(~group, z$samples)
    colnames(design) <-gsub("group","",colnames(design))
    z <- z[filterByExpr(z, design),,keep.lib.sizes=F]
    z <- calcNormFactors(z)
    z <- estimateDisp(z, design, robust=T)
    z <- glmQLFit(z, design, robust = T)
    print(head(residuals(z, typ="deviance")))
    
    # Estimate nuisance parameter W_1 using the RUVr method
    seq <- newSeqExpressionSet(counts = z$counts)
    seqUQ <- betweenLaneNormalization(seq, which="upper")
    ruvr <- RUVr(
      seqUQ, cIdx = row.names(z),
      k = 1, residuals = residuals(z, type="deviance")
    )
    y$samples$W_1 <- pData(ruvr)[row.names(y$samples), "W_1"]
    design<-model.matrix(~group + W_1, y$samples)
    colnames(design) <-gsub("group","",colnames(design))
    
    y <- y[filterByExpr(y, design),,keep.lib.sizes=F]
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design, robust=T)
    y$genes$com_disp <- y$common_disp <- y$common.dispersion
    y$genes$trend_disp <- y$trended.dispersion
    y$genes$tag_disp <- y$tagwise.dispersion
    
    # NOTE -- FPKM Calculation Fails to account for RUV adjustment !!!
    # THIS may need to be FIXED if RUV is intended for use !!!
    fpkm<-rpkm(y, normalized.lib.sizes = T, gene.length = "eu_length")
    for(g in unique(y$samples$group)){
      s<-row.names(y$samples[y$samples$group == g,])
      y$genes[paste0(g,"_Avg_FPKM")]<-apply(fpkm[,s],1,mean)
    }
    
    fit <- glmQLFit(y, design, robust = T)
    return(
      list(
        dge=y, fit=fit, ruvr=ruvr, design=design
      )
    )
    
  } else if(norm == "NONE"){
    design<-model.matrix(~group, y$samples)
    colnames(design) <-gsub("group","",colnames(design))
    y <- y[filterByExpr(y, design),,keep.lib.sizes=F]
    y <- estimateDisp(y, design, robust=T)
    ## Leave scaling factors as 1 
    y$genes$com_disp <- y$common.dispersion
    y$genes$trend_disp <- y$trended.dispersion
    y$genes$tag_disp <- y$tagwise.dispersion
    
    fpkm<-rpkm(y, normalized.lib.sizes = T, gene.length = "eu_length")
    for(g in unique(y$samples$group)){
      s<-row.names(y$samples[y$samples$group == g,])
      y$genes[paste0(g,"_Avg_FPKM")]<-apply(fpkm[,s],1,mean)
    }
    
    fit <- glmQLFit(y, design, robust = T)
    return(
      list(
        dge=y, fit=fit, design=design
      )
    )
    
  } else {
    return(y)
  }
}

#### Define function to generate DEG Tables from a fit or DGEList ############
genPairwiseDegTable<-function(y, group1, group2, design){
  print(
    paste(
      group1, "Samples:",
      paste(
        y$samples %>% filter(group == group1) %>% pull("label"),
        collapse = ", "
      )
    )
  )
  print(
    paste(
      group2, "Samples:",
      paste(
        y$samples %>% filter(group == group2) %>% pull("label"),
        collapse = ", "
      )
    )
  )
  
  if(class(y) == "DGEGLM"){
    cn<-sapply(
      colnames(design), 
      function(n){
        ifelse(
          n == group1, -1, 
          ifelse(n == group2, 1, 0
          )
        )
      }
    )
    print(cn)
    return(
      as.data.frame(
        topTags(
          glmQLFTest(y, contrast=cn), n=Inf
        )
      ) %>% 
        dplyr::mutate(
          Test = "QLFTest",
          Group_1 = group1, 
          Group_2 = group2
        ) %>%
        dplyr::select(
          Test, Group_1, Group_2,
          gene_id, logFC, logCPM, PValue, FDR, 
          Avg1 = as.name(paste0(group1, "_Avg_FPKM")),
          Avg2 = as.name(paste0(group2, "_Avg_FPKM"))
        )
    )
  } else if(class(y) == "DGEList") {
    return(
      as.data.frame(
        topTags(
          exactTest(y, pair=c(group1, group2)), n=Inf
        )
      ) %>% 
        dplyr::mutate(
          Test = "ExactTest",
          Group_1 = group1, 
          Group_2 = group2
        ) %>%
        dplyr::select(
          Test, Group_1, Group_2,
          gene_id, logFC, logCPM, PValue, FDR, 
          Avg1 = as.name(paste0(group1, "_Avg_FPKM")),
          Avg2 = as.name(paste0(group2, "_Avg_FPKM"))
        )
    )
  } else if(class(y) == "MArrayLM"){
    cft <- contrasts.fit(
      fit, contrasts=makeContrasts(
        contrasts = paste(group2, group2, sep="-"),
        levels=colnames(design)
      )
    )
    ebs <- eBayes(cft)
    return(
      as.data.frame(
        topTable(cft, n=Inf)
      ) %>% 
        dplyr::mutate(
          Test = "LimmaVoom",
          Group_1 = group1, 
          Group_2 = group2
        ) %>%
        dplyr::select(
          Test, Group_1, Group_2,
          gene_id, logFC, logCPM=AveExper, PValue=P.Value, FDR=adj.P.Val, 
          Avg1 = as.name(paste0(group1, "_Avg_FPKM")),
          Avg2 = as.name(paste0(group2, "_Avg_FPKM"))
        )
    )
  }
}

## Generate a DEG table using coefficients to a design matrix
genDesignCoefDegTable<-function(y, design, coef, group_labels){
  
  # Print Experimental Design and samples associated with each coeffient
  print(design)
  for(i in c(1,coef)){
    print(colnames(design)[i])
    print(y$samples[which(design[,i] == 1),"label"])
  }
  
  if(class(y) == "DGEGLM"){
    print(coef)
    return(
      as.data.frame(
        topTags(
          glmQLFTest(y, coef=coef), n=Inf
        )
      ) %>% 
        dplyr::mutate(
          Test = "QLFTest",
          Group_1 = group_labels[1], 
          Group_2 = group_labels[2],
          Avg1 = 2^((logCPM - (logFC/2)) - log2(eu_length/1000)),
          Avg2 = 2^((logCPM + (logFC/2)) - log2(eu_length/1000))
        ) %>%
        dplyr::select(
          Test, Group_1, Group_2,
          gene_id, logFC, logCPM, PValue, FDR, 
          Avg1, Avg2
        )
    )
  } else if(class(y) == "DGEList") {
    print("Analysis Reqires a fitted model")
    return(NULL)
    
  } else if(class(y) == "MArrayLM"){
    ebs <- eBayes(y)
    return(
      x <- as.data.frame(
        topTable(ebs, n=Inf, coef=coef)
      ) %>% 
        dplyr::mutate(
          Test = "LimmaVoom",
          Group_1 = group_labels[1], 
          Group_2 = group_labels[2],
          Avg1 = 2^((AveExpr - (logFC/2)) - log2(eu_length/1000)),
          Avg2 = 2^((AveExpr + (logFC/2)) - log2(eu_length/1000))
        ) %>%
        dplyr::select(
          Test, Group_1, Group_2,
          gene_id, logFC, logCPM=AveExpr, PValue=P.Value, FDR=adj.P.Val, 
          Avg1, Avg2
        )
    )
  }
}


# Iterate over a set of contrasts, generate DEG Tables and
#  DEG Summary tables for each contrast. 
iterate_edgeR_pairwise_contrasts <- function(
  dge, fit, cntmat=cntmat, df=df, design=design,
  deg=deg,respath="LIRTS_DEG_Analysis_results",
  prefix="DBI"
){
  print(cntmat)
  for (c in colnames(cntmat)) {             
    # Run Exact Tests
    pair<-c(
      names(cntmat[,c])[cntmat[,c] == -1], 
      names(cntmat[,c])[cntmat[,c] == 1]
    )
    print(pair)
    
    deg.et<-genPairwiseDegTable(dge, pair[1], pair[2], design)
    deg.qt<-genPairwiseDegTable(fit, pair[1], pair[2], design)
    
    deg <- bind_rows(
      deg,
      deg.et %>%
        mutate(Samples=prefix) %>%
        tibble::remove_rownames(),
      deg.qt %>%
        mutate(Samples=prefix) %>%
        tibble::remove_rownames()
    )
    
    # Save DEG Tables
    fn<-paste(
      respath,"/",prefix,"_",c,"_Exact_Test_DEG.tsv", sep=""
    )
    write.table(
      deg.et, fn, col.names =T, quote = F, sep="\t", row.names = F
    )
    
    fn<-paste(
      respath,"/",prefix,"_",c,"_QLFTest_DEG.tsv", sep=""
    )
    write.table(
      deg.qt, fn, col.names=T, quote = F, sep="\t", row.names = F
    )
    
    dg<-degSummary(                         # Generate Exact Test Summary tables
      deg.et,
      lfc = 'logFC',
      fdr = 'FDR', 
      Avg1 = pair[1],
      Avg2 = pair[2]
    )%>%
      mutate(Samples=prefix)
    
    dg$contrast<-c
    dg$test<-"Exact Test"
    df<-bind_rows(df, dg)
    
    dg<-degSummary(                         # Generate QLF Test Summary tables 
      deg.qt,
      lfc = "logFC",
      fdr = 'FDR', 
      Avg1 = pair[1],
      Avg2 = pair[2]
    )%>%
      mutate(Samples=prefix)
    dg$contrast<-c
    dg$test<-"QLFTest"
    df<-bind_rows(df, dg)
    
  }
  return(list(df, deg))
}


# Iterate over a set of design matrix coefficients, 
# generate DEG Tables, DEG Summary tables
iterate_edgeR_design_coefficients <- function(
  dge, fit, coefs=c(2,3), df=df, design=design, deg=deg,
  respath="LIRTS_DEG_Analysis_results", prefix="DBI", 
  group_label_list=list(c("ctrl", "trt"), c("mut", "wt")),
  filt="NONE"
){
  print(group_label_list)
  for (c in 1:length(coefs)) {             
    deg.qt<-genDesignCoefDegTable(
      fit, design, coef=coefs[c], 
      group_labels = group_label_list[[c]]
    )
    
    deg <- bind_rows(
      deg, deg.qt %>%
        mutate(
          Filtered= filt,
          Partition=prefix
        ) %>%
        tibble::remove_rownames()
    )
    
    # Save DEG Tables
    fn<-paste(
      respath,"/",prefix,"_",c,"_QLFTest_DEG.tsv", sep=""
    )
    write.table(
      deg.qt, fn, col.names=T, quote = F, sep="\t", row.names = F
    )
    
    dg<-degSummary(                         # Generate QLF Test Summary tables 
      deg.qt,
      lfc = "logFC",
      fdr = 'FDR', 
      Avg1 = "Avg1",
      Avg2 = "Avg2"
    )%>%
      mutate(Samples=prefix)
    dg$contrast<-colnames(design)[coefs[c]]
    dg$test<-"QLFTest"
    df<-bind_rows(df, dg)
    
  }
  return(list(df, deg))
}

################ Extract FPKM Matrix for Features of interest ################
process_selected_features <- function(
  dge, design, genes=NULL, counts=NULL, 
  prefix="Global_Wildtype_Top_Genes",
  respath="results/",
  color_attrib="hours_pcs", 
  shape_attrib = "batch"
){
  
  # colors=c(
  #   "#000000", "#E69F00", "#56B4E9", "#009E73",
  #   "#0072B2", "#D55E00", "#CC79A7"
  # )
  colors=c(
    RColorBrewer::brewer.pal(9, name="RdYlBu")[1:3],
    RColorBrewer::brewer.pal(9, name="RdYlBu")[7:9]
  )[c(1,6,2,5,3,4)]
  
  
  if(!is.null(counts)){
    dge$counts <- counts
  }
  obj <- process_edgeR_ByDesign(y=dge, genes=genes, design=design)
  diagnostic_plots(
    obj$dge, prefix = prefix, 
    color_attrib = color_attrib,
    shape_attrib = shape_attrib
  )
  fpkm<-as.data.frame(
    rpkm(
      obj$dge, gene.length="eu_length", log=T,
      normalized.lib.sizes = T
    )
  )
  colnames(fpkm) <- obj$dge$samples$label
  
  ## Plot samplewise correlation matrix after excluding genes with 
  ## a strong batch response (DBI vs DNA1, DNA2 or DNA3)
  cm <- cor(fpkm, method = "spearman")
  annot <- obj$dge$samples[,c("hours_pcs", "batch", "label")]
  annot <- annot %>%
    tibble::remove_rownames()%>%
    tibble::column_to_rownames("label")
  
  annot_colors=list(
    hours_pcs=colors[1:6],
    batch=colors[2:5]
  )
  names(annot_colors$hours_pcs)=levels(annot$hours_pcs)
  names(annot_colors$batch)=levels(annot$batch)
  
  pheatmap(
    cm, annotation_col =annot,
    annotation_colors = annot_colors,
    filename = paste0(
      respath, prefix, "_cormat.png"
    ),height = 6, width=8
  )
  
  ## Save FPKM matrix for genes with a significant batch dependent logFC
  ## less than 2 (four fold difference) (DBI vs DNA1, DNA2 or DNA3)
  fpkm<-tibble::rownames_to_column(fpkm,"gene_id")
  write.csv(
    fpkm, row.names = F, quote = F,
    paste0(
      "LIRTS_DEG_Analysis_results/",
      prefix,"_TMM-FPKM_Matrix.csv"
    )
  )
  return(obj)
}


