########################         HEADER BLOCK      ###########################
#  File:    Generate_Pax6_DEG_Tables.R                                       #
#  Purpose: Construct DEG Tables from the Pax6 Study                         #
#  Created: Sep 11, 2021                                                     #
#  Author:  Adam Faranda                                                     #
#                                                                            #
##############################################################################

##########        Load Libraries and setup environment            ############
library(edgeR)
library(dplyr)
library(ggplot2)
library(synapser)
source("scripts/Overlap_Comparison_Functions.R")
source("scripts/Wrap_edgeR_Functions.R")


setwd("~/Documents/11Sep2021_Pax6_Study_DEG_Analysis")
wd<-getwd()
data_dir="data"

# Get Synapse ID's
synLogin()
syn_project <- synFindEntityId("Pax6_Happloinsuficiency_In_The_Lens")
syn_code_dir <- synFindEntityId("code", parent=syn_project)
syn_data_dir <- synFindEntityId("data", parent=syn_project)
syn_deg_dir <- synFindEntityId("DEG_Tables", parent = syn_project)
syn_dgelists <- synFindEntityId(
  "pax6_master_dgelists.Rdata", 
  parent=syn_data_dir
)

# Import Local DGElists or fetch from Synapse
if(file.exists("data/pax6_master_dgelists.Rdata")){
  load("data/pax6_master_dgelists.Rdata")
} else if(!is.null(syn_dgelists)){
  synGet(
    syn_dgelists, downloadLocation=data_dir
  )
  load("data/pax6_master_dgelists.Rdata")
}else{
  source("scripts/Prepare_Expression_DGELists.R")
  syn_dgelists <- synFindEntityId(
    "pax6_master_dgelists.Rdata", 
    parent=syn_data_dir
  )
}

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
 
  diagnostic_plots(
    set$dge, color_attrib = "genotype",
    shape_attrib = "cell_type",
    respath = "results", prefix = paste0("All_",filt=filt)
    
  )
   
  dge$samples$sample <- row.names(dge$samples)
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
        genPairwiseDegTable(set[[obj]], cn[1],cn[2], set[["design"]]) %>%
          dplyr::mutate(Filtered = filt, Partition = "All")
      ) %>% tibble::remove_rownames()
    }
  }
  
  # # Run 2 Way interaction models
  # pax6.deg_master <- iterate_edgeR_design_coefficients(
  #   dge=set$dge, fit=set$fit, deg=pax6.deg_master, filt=filt,
  #   design=design, respath = "results", prefix="2Way",
  #   df=data.frame(), coefs=c(2:4), group_label_list = list(
  #     c("Epi", "Fib"), c("WT", "P6"), c("WTEpi", "P6Fib")
  #   )
  # )[[2]]
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

    for(i in c("dge", "fit")){
      pax6.deg_master <- bind_rows(
        pax6.deg_master, 
        genPairwiseDegTable(set[[i]], cn[1],cn[2], set[["design"]]) %>%
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
  design <- model.matrix(~ cell_type * genotype, dge$samples)
  obj <- processByDesign(
    y=dge, design=design,
    norm=ifelse(filt != "degnorm", "TMM", "NONE")
  )
  
  diagnostic_plots(
    obj$dge, color_attrib = "genotype",
    shape_attrib = "cell_type",
    respath = "results", prefix = paste0("2Way_",filt=filt)
    
  )
  
  pax6.deg_master <- iterate_edgeR_design_coefficients(
    dge=obj$dge, fit=obj$fit, deg=pax6.deg_master, filt=filt,
    design=design, respath = "results", prefix="2Way",
    df=data.frame(), coefs=c(2:4), group_label_list = list(
      c("Epi", "Fib"), c("WT", "P6"), c("REF", "INX")
    )
  )[[2]]
}

save(
  pax6.deg_master, 
  file="results/pax6_deg_tables.Rdata"
)

######################  Push script and data to Synapse ######################
script_path <- "scripts/Generate_Pax6_DEG_Tables.R"
deg_table_path <- "results/pax6_deg_tables.Rdata"
syn_script <- File(
  path=script_path,
  parent=syn_code_dir
)

syn_script <- synStore(
  syn_script,
  used = syn_dgelists
)

syn_deg_tables <- File(
  path=deg_table_path,
  parent=syn_deg_dir
)

syn_deg_tables <- synStore(
  syn_deg_tables,
  executed = syn_script
)


