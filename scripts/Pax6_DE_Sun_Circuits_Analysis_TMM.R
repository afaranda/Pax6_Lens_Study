##############################################################################
#                                                                            #
#  File: Pax6_DE_Sun_Circuits_Analysis_TMM.R                                 #
#  Author: Adam Faranda                                                      #
#  Created: June 18, 2020                                                    #
#  Purpose:                                                                  #
#         Add differential expression estimates to networks curated from     #
#         Sun et. al. 2015.  Prepare a table of node anntations to be used   #
#         with these networks for visualization in Cytoscape                 #
#                                                                            #
##############################################################################

################## Load Libraries and Source Dependencies ####################
library(edgeR)
library(tidyr)
library(tibble)
library(dplyr)
library(tibble)
library(synapser)
options(echo=T)

# Enter Working Directory and Load Raw Data
setwd('/Users/adam/Documents/23Feb2022_Pax6_Study_DEG_Analysis/')
source('scripts/Overlap_Comparison_Functions.R')
source('scripts/Wrap_edgeR_Functions.R')
wd<-getwd()
results<-paste(wd,'results/Pax6_Targets/',sep='/')
if(!dir.exists(results)){
  dir.create(results)
}

data_dir <- paste0(wd,"/data")         ## Local data directory

######################### Setup Synapse Connection ###########################
synLogin()
syn_project <- synFindEntityId("Pax6_Happloinsuficiency_In_The_Lens")
syn_code_dir <- synFindEntityId("code", parent=syn_project)
syn_data_dir <- synFindEntityId("data", parent=syn_project)
syn_fig_dir <- synFindEntityId("Figures", parent=syn_project)
syn_deg_tbl_dir <- synFindEntityId("DEG_Tables", parent = syn_project)
syn_deg_sps_dir <- synFindEntityId("Pax6_Targets", parent = syn_project)

if(is.null(syn_deg_sps_dir)){
  folder <- Folder(
    name="Pax6_Targets",
    parent=syn_project
  )
  syn_deg_sps_dir <- synStore(folder)
} 

## Check for the DEG_Tables folder
syn_deg_tbl_dir <- synapser::synFindEntityId(
  "DEG_Tables",
  syn_project
)

######################## Load in Pax6 Study DGEList ##########################
syn_dgelists <- synFindEntityId(
  "pax6_master_dgelists.Rdata",
  parent=syn_data_dir
)

if(file.exists('data/pax6_master_dgelists.Rdata')){
  load('data/pax6_master_dgelists.Rdata')
} else if(!is.null(syn_dgelists)){
  synGet(
    syn_dgelists, downloadLocation=data_dir
  )
  load('data/pax6_master_dgelists.Rdata')
}else{
  stop("Run Prepare_Expression_DGELists.R first")
}

## Get Ribo Filtered Counts
pax6.master$samples$group <- factor(
  gsub("_ribo_[123]","", pax6.master$samples$label),
  levels=c("WTE", "WTF", "P6E", "P6F")
) 

pax6.master$genes$DESCRIPTION <- gsub(
  " \\[.*\\]", "", pax6.master$genes$DESCRIPTION
)

# Drop unfiltered counts
pax6.master <- pax6.master[,which(pax6.master$samples$ribo_filter=="ribo")]
pax6.master$samples$label <- gsub("_ribo", "", pax6.master$samples$label)
colnames(pax6.master) <- pax6.master$samples$label

## Get DegNorm / Ribo Filtered Counts
pax6.master.dgn$samples$group <- factor(
  gsub("_degnorm_[123]","", pax6.master.dgn$samples$label),
  levels=c("WTE", "WTF", "P6E", "P6F")
) 
pax6.master.dgn$samples$label <- gsub(
  "_degnorm", "", pax6.master.dgn$samples$label
)
colnames(pax6.master.dgn) <- pax6.master.dgn$samples$label

pax6.master.dgn$genes$DESCRIPTION <- gsub(
  " \\[.*\\]", "", pax6.master.dgn$genes$DESCRIPTION
)

######################### Load in Master DEG Table ###########################
syn_deg_master <- synFindEntityId(
  "pax6_deg_tables.Rdata", 
  parent=syn_deg_tbl_dir
)

if(file.exists("results/pax6_deg_tables.Rdata")){
  load("results/pax6_deg_tables.Rdata")
} else if(!is.null(syn_deg_master)){
  synGet(
    syn_deg_master, downloadLocation="results"
  )
  load("results/pax6_deg_tables.Rdata")
}else{
  stop("Run Generate_Pax6_DEG_Tables.R first")
}



########################## Load in Pax6 Targets ##############################
syn_pax6_circuits <- synFindEntityId(
  "Sun_Pax6_Networks.tsv",
  parent = syn_data_dir
)

if(file.exists("data/Sun_Pax6_Networks.tsv")){
  pax6_circuits <- read.table(
    "data/Sun_Pax6_Networks.tsv",
    header=T, sep='\t'
  )
  # names(pax6_circuits) <- c(
  #   "TFactor", "Target", "Direction","PMID"
  # )
  # 
} else if(!is.null(syn_pax6_circuits)){
  synGet(
    syn_pax6_circuits, downloadLocation="data"
  )
  
  pax6_circuits <- read.table(
    "data/Sun_Pax6_Networks.tsv",
  )
  # names(pax6_circuits) <- c(
  #   "TFactor", "Target", "Direction","PMID"
  # )
  # 
}else{
  stop("Prepare TRRUST Targets file First")
}

## Final list of existing data files used by this script
used_files <- list(
 syn_pax6_circuits,
  syn_dgelists,
  syn_deg_master
)
################## Pivot DE Results and extract Node values ##################
result_files <- c()
Sun_Network_DE_Measures <- data.frame(
  Symbol = union(pax6_circuits$Regulator, pax6_circuits$Target)
) %>%
  inner_join(
    pax6.master$genes %>%
      select(
        SYMBOL, gene_id
      ),
    by=c(Symbol ="SYMBOL")
  ) %>%
  left_join(
    pax6.deg_master %>%
      filter(
        Test == "ExactTest" & Group_1 == "WTE" & Group_2 == "WTF" &
          Filtered == "ribo" & Partition == "Pair"
      ) %>%
      select(
        gene_id,
        WTFvsWTE_logFC = logFC,
        WTFvsWTE_FDR = FDR,
        Cell_WTE_FPKM = Avg1,
        Cell_WTF_FPKM = Avg2
      ),
    by="gene_id"
  ) %>%
  left_join(
    pax6.deg_master %>%
      filter(
        Test == "ExactTest" & Group_1 == "P6E" & Group_2 == "P6F" &
          Filtered == "ribo" & Partition == "Pair"
      ) %>%
      select(
        gene_id,
        P6FvsP6E_logFC = logFC,
        P6FvsP6E_FDR = FDR,
        Cell_P6E_FPKM = Avg1,
        Cell_P6F_FPKM = Avg2
      ),
    by="gene_id"
  ) %>%
  left_join(
    pax6.deg_master %>%
      filter(
        Test == "ExactTest" & Group_1 == "WTE" & Group_2 == "P6E" &
          Filtered == "ribo" & Partition == "Pair"
      ) %>%
      select(
        gene_id,
        P6EvsWTE_logFC = logFC,
        P6EvsWTE_FDR = FDR,
        Geno_WTE_FPKM = Avg1,
        Geno_P6E_FPKM = Avg2
      ),
    by="gene_id"
  ) %>%
  left_join(
    pax6.deg_master %>%
      filter(
        Test == "ExactTest" & Group_1 == "WTF" & Group_2 == "P6F" &
          Filtered == "ribo" & Partition == "Pair"
      ) %>%
      select(
        gene_id,
        P6FvsWTF_logFC = logFC,
        P6FvsWTF_FDR = FDR,
        Geno_WTF_FPKM = Avg1,
        Geno_P6F_FPKM = Avg2
      ),
    by="gene_id"
  )
    
Sun_Network_DE_Measures[is.na(Sun_Network_DE_Measures)] <- 0

fn <- "Sun_Pax6_Network_Node_Measurements.csv"
path <- paste(results, fn, sep="/")
write.csv(Sun_Network_DE_Measures, path, row.names = F)
result_files <- append(result_files, path)

######################  Push script and data to Synapse ######################

# Add this script to the code dir
syn_script <- synFindEntityId(
  "Pax6_DE_Sun_Circuits_Analysis_TMM.R",
  parent=syn_code_dir
)

if(is.null(syn_script)){
  syn_script <- File(
    path="scripts/Pax6_DE_Sun_Circuits_Analysis_TMM.R",
    parent=syn_code_dir
  )

  syn_script <- synStore(
    syn_script
  )
}

synSetProvenance(
  syn_script,
  activity = Activity(
    name = "Pax6_Networks_DEG",
    description = "Extract DE Measures for Sun Pax6 Circuits",
    used=used_files
  )
)

syn_act <- Activity(
  name="upload_analysis_results",
  description="upload analysis results"
)
syn_act$executed(syn_script)

for(file in result_files){
  synapse_push <- File(
    path=file,
    parent=syn_deg_sps_dir
  )
  
  synapse_push <- synStore(
    synapse_push
  )

  synSetProvenance(
    synapse_push,
    syn_act
  )
}
