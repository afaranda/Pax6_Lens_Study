########################         HEADER BLOCK      ###########################
#  File:    Generate_Pax6_DEG_Tables.R                                       #
#  Purpose: Construct DEG Tables from the Pax6 Study                         #
#  Created: Sep 11, 2021                                                     #
#  Author:  Adam Faranda                                                     #
#                                                                            #
##############################################################################

#########      Load Libraries and import expression data          ############
library(edgeR)
library(dplyr)
library(ggplot2)
setwd("~/Documents/11Sep2021_Pax6_Study_DEG_Analysis")
wd<-getwd()
data_dir="data"
if(file.exists("data/pax6_master_dgelists.Rdata")){
  load("data/pax6_master_dgelists.Rdata")
} else {
  source("scripts/Prepare_Expression_DGELists.R")
}

source("scripts/Wrap_edgeR_Functions.R")
source("scripts/PrincipalComponents.R")

pax6.deg_master<-data.frame()       # Empty data frame to store all DEG results

########          Generate DEG Tables over All Genomic contrasts         #####
for ( filt in levels(pax6.master$samples$ribo_filter)){
  if (filt != "degnorm"){
    master <- pax6.master
  } else {
    master <- pax6.master.dgn
  }
  # Construct and normalize all-sample DGEList object, fit QLF Model
  dge <- master[,master$samples$ribo_filter == filt ]
  dge$samples$group<-factor(                # Add grouping factor
    paste0(
      dge$samples$genotype, 
      gsub("(i|b|p)","", dge$samples$cell_type)
    ), levels=c("WTE", "WTF", "P6E", "P6F")
  )
  set <- subsetDGEListByGroups(
    dge, groups=as.character(dge$samples$group),
    norm=ifelse(filt=="degnorm", "NONE", "TMM")
  )
  pax6.disp_master <- data.frame(
      Partition = "All",
      Filtered = filt,
      Contrast = "Multiple",
      Disp = set[["dge"]]$common.dispersion
    )
   
  dge$samples$sample <- row.names(dge$samples)
  ggsave(
    filename = "results/P6D_All_Sample_PCA.jpg",
    plotPrinComp(
      cpm(dge, log=T), ft=dge$samples,idCol = 0,
      groupCol = "group",
    ), width=4, height=3
  )
  for(                              # Iterate over relevant contrasts
    cn in list(
      c("WTE", "WTF"),
      c("WTE", "P6E"),
      c("WTF", "P6F"),
      c("P6E", "P6F")
    )
  ){
   
    for(obj in c("dge", "fit")){
      pax6.deg_master <- bind_rows(
        pax6.deg_master, 
        genDegTable(set[[obj]], cn[1],cn[2], set[["design"]]) %>%
          dplyr::mutate(Filtered = filt, Partition = "All")
      ) %>% tibble::remove_rownames()
    }
  }
}

#### Generate DEG Tables over all contrasts, Pairwise Partitions #####
for (filt in c("geno", "ribo", "degnorm") ){
  if (filt != "degnorm"){
    master <- pax6.master      # Use the standard count DEGlist
  } else {
    master <- pax6.master.dgn  # Use DEGNorm adjusted counts
  }
  
  # Construct and normalize all-sample DGEList object, fit QLF Model
  dge <- master[,master$samples$ribo_filter == filt ]
  dge$samples$group<-factor(                # Add grouping factor
    paste0(
      dge$samples$genotype, 
      gsub("(i|b|p)","", dge$samples$cell_type)
    ), levels=c("WTE", "WTF", "P6E", "P6F")
  )
 
  for(
    cn in list(
      c("WTE", "WTF"),
      c("WTE", "P6E"),
      c("WTF", "P6F"),
      c("P6E", "P6F")
    )
  ){
    set <- subsetDGEListByGroups(
      dge, groups=c(cn[1], cn[2]),
      norm=ifelse(filt=="degnorm", "NONE", "TMM")
    )
    pax6.disp_master <- bind_rows(
      pax6.disp_master,
      data.frame(
        Partition = "Pair",
        Filtered = filt,
        Contrast = paste0(cn[2], "v", cn[1]),
        Disp = set[["dge"]]$common.dispersion
      )
    )    
    
    for(i in c("dge", "fit")){
      pax6.deg_master <- bind_rows(
        pax6.deg_master, 
        genDegTable(set[[i]], cn[1],cn[2], set[["design"]]) %>%
          dplyr::mutate(Filtered = filt, Partition = "Pair")
      ) %>% tibble::remove_rownames()
    }
  }
}

############## Generate DEG Tables for factor interaction model ###############
for (filt in c("geno", "ribo", "degnorm") ){
  if (filt != "degnorm"){
    master <- pax6.master      # Use the standard count DEGlist
  } else {
    master <- pax6.master.dgn  # Use DEGNorm adjusted counts
  }
  
  # Construct and normalize all-sample DGEList object, fit QLF Model
  dge <- master[,master$samples$ribo_filter == filt ]
  design <- model.matrix(~genotype * cell_type, dge$samples)
  obj <- processByDesign(y=dge, design=design)
  
  dge$samples$group<-factor(                # Add grouping factor
    paste0(
      dge$samples$genotype, 
      gsub("(i|b|p)","", dge$samples$cell_type)
    ), levels=c("WTE", "WTF", "P6E", "P6F")
  )
 
  for(
    cn in list(
      c("WTE", "WTF"),
      c("WTE", "P6E"),
      c("WTF", "P6F"),
      c("P6E", "P6F")
    )
  ){
    set <- subsetDGEListByGroups(
      dge, groups=c(cn[1], cn[2]),
      norm=ifelse(filt=="degnorm", "NONE", "TMM")
    )
    pax6.disp_master <- bind_rows(
      pax6.disp_master,
      data.frame(
        Partition = "Pair",
        Filtered = filt,
        Contrast = paste0(cn[2], "v", cn[1]),
        Disp = set[["dge"]]$common.dispersion
      )
    )    
    
    for(i in c("dge", "fit")){
      pax6.deg_master <- bind_rows(
        pax6.deg_master, 
        genDegTable(set[[i]], cn[1],cn[2], set[["design"]]) %>%
          dplyr::mutate(Filtered = filt, Partition = "Pair")
      ) %>% tibble::remove_rownames()
    }
  }
}



save(
  pax6.deg_master, 
  file="results/pax6_deg_tables.Rdata"
)



