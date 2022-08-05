##############################################################################
#                                                                            #
#  File: Pax6_DE_KEGG_Tables.R                                               #
#  Author: Adam Faranda                                                      #
#  Created: March 3, 2022                                                    #
#  Purpose:                                                                  #
#         Add differential expression estimates to KEGG Pathways enriched in #
#         injury responsive Pax6 LEC DEG for visualization in Cytoscape      #
#                                                                            #
##############################################################################

################## Load Libraries and Source Dependencies ####################
library(edgeR)
library(tidyr)
library(tibble)
library(dplyr)
library(tibble)
library(KEGGgraph)
#library(synapser)
options(echo=T)

# Enter Working Directory and Load Raw Data
setwd('/Users/adam/Documents/23Feb2022_Pax6_Study_DEG_Analysis/')
source('scripts/synapse_reticulate_wrapper.R')
source('scripts/Overlap_Comparison_Functions.R')
source('scripts/Wrap_edgeR_Functions.R')
wd<-getwd()
results<-paste(wd,'results/Pax6_KEGG_Pathways/',sep='/')
if(!dir.exists(results)){
  dir.create(results)
}

data_dir <- paste0(wd,"/data")         ## Local data directory

######################### Setup Synapse Connection ###########################
#synLogin()
syn_project <- synFindEntityId("Pax6_Happloinsuficiency_In_The_Lens")
syn_code_dir <- synFindEntityId("code", parent=syn_project)
syn_data_dir <- synFindEntityId("data", parent=syn_project)
syn_fig_dir <- synFindEntityId("Figures", parent=syn_project)
syn_deg_tbl_dir <- synFindEntityId("DEG_Tables", parent = syn_project)
syn_deg_sps_dir <- synFindEntityId("Pax6_KEGG_Pathways", parent = syn_project)

if(is.null(syn_deg_sps_dir)){
  folder <- Folder(
    name="Pax6_KEGG_Pathways",
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
syn_trrust_targets <- synFindEntityId(
  "TRRUST_Pax6_targets.mouse_03Feb2022.tsv",
  parent = syn_data_dir
)

if(file.exists("data/TRRUST_Pax6_targets.mouse_03Feb2022.tsv")){
  trrust_targets <- read.table(
    "data/TRRUST_Pax6_targets.mouse_03Feb2022.tsv",
  )
  names(trrust_targets) <- c(
    "TFactor", "Target", "Direction","PMID"
  )
  
} else if(!is.null(syn_trrust_targets)){
  synGet(
    syn_trrust_targets, downloadLocation="data"
  )
  
  trrust_targets <- read.table(
    "data/TRRUST_Pax6_targets.mouse_03Feb2022.tsv",
  )
  names(trrust_targets) <- c(
    "TFactor", "Target", "Direction","PMID"
  )
  
}else{
  stop("Prepare TRRUST Targets file First")
}

syn_sun_targets <- synFindEntityId(
  "Sun2015_Pax6_DE_Tagrets.tsv",
  parent = syn_data_dir
)

if(file.exists("data/Sun2015_Pax6_DE_Tagrets.tsv")){
  sun_targets <- read.table(
    "data/Sun2015_Pax6_DE_Tagrets.tsv",
    header=T
  )
  
} else if(!is.null(syn_sun_targets)){
  synGet(
    syn_sun_targets, downloadLocation="data"
  )
  
  sun_targets <- read.table(
    "data/Sun2015_Pax6_DE_Tagrets.tsv",
    header=T
  )
  
}else{
  stop("Prepare TRRUST Targets file First")
}

syn_msig_targets <- synFindEntityId(
  "PAX6_TARGET_GENES_MSigDB_GeneSet.txt",
  parent = syn_data_dir
)

if(file.exists("data/PAX6_TARGET_GENES_MSigDB_GeneSet.txt")){
  msig_targets <- read.table(
    "data/PAX6_TARGET_GENES_MSigDB_GeneSet.txt",
    header=F
  )
} else if(!is.null(syn_msig_targets)){
  synGet(
    syn_msig_targets, downloadLocation="data"
  )
  
  msig_targets <- read.table(
    "data/PAX6_TARGET_GENES_MSigDB_GeneSet.txt",
    header=F
  )
  
}else{
  stop("Prepare MSigDB Targets file First")
}


########################## Load in KEGG Pathways  ############################
syn_mmu_04151 <- synFindEntityId(
  "mmu04151.xml", parent = syn_data_dir
)

if(file.exists("data/mmu04151.xml")){
  mmu04151 <- parseKGML2DataFrame(
    "data/mmu04151.xml"
  )
} else if(!is.null(syn_mmu_04151)){
  synGet(
    syn_mmu_04151, downloadLocation="data"
  )
  
  mmu04151 <- parseKGML2DataFrame(
    "data/mmu04151.xml"
  )
  
}else{
  stop("Download KEGG Data first")
}


syn_mmu_04350 <- synFindEntityId(
  "mmu04350.xml", parent = syn_data_dir
)

if(file.exists("data/mmu04350.xml")){
  mmu04350 <- parseKGML2DataFrame(
    "data/mmu04350.xml"
  )
} else if(!is.null(syn_mmu_04350)){
  synGet(
    syn_mmu_04350, downloadLocation="data"
  )
  
  mmu04350 <- parseKGML2DataFrame(
    "data/mmu04350.xml"
  )
  
}else{
  stop("Download KEGG Data first")
}

######################## Load in Long Injury Table ###########################
syn_injury_deg <- synFindEntityId(
  "Aging_DEG_In_Pax6_LEC_Long_Table.csv",
  parent = synFindEntityId(
    "Pax6_Aging_Analysis",
    parent = syn_project
  )
)

if(
  file.exists(
    "results/Aging_Analysis/Aging_DEG_In_Pax6_LEC_Long_Table.csv"
  )
){
  injury_deg <- read.csv(
    "results/Aging_Analysis/Aging_DEG_In_Pax6_LEC_Long_Table.csv",
    row.names = 1
  )
  
} else if(!is.null(syn_sun_targets)){
  synGet(
    syn_sun_targets, downloadLocation="results/Aging_Analysis"
  )
  
  injury_deg <- read.csv(
    "results/Aging_Analysis/Aging_DEG_In_Pax6_LEC_Long_Table.csv",
    row.names = 1
  )
  
}else{
  stop("Run Pax6_DE_Aging_Analysis_TMM.R first")
}

## Final list of existing data files used by this script
used_files <- list(
  syn_mmu_04151,
  syn_mmu_04350,
  syn_dgelists,
  syn_deg_master,
  syn_injury_deg
)
################## Pivot DE Results and extract Node values ##################
result_files <- c()
for(pathway in c("mmu04151", "mmu04350")){
  kegg <- get(pathway) %>%
    mutate(
      Interaction=paste0(
        toupper(substr(subtype,1,1)),
        substr(subtype,2, nchar(subtype))
      )
    )
  
  kegg_pathway_measures <- data.frame(
    kegg_id = union(
      kegg$from, kegg$to
    )
  ) %>%
    inner_join(
      pax6.master$genes %>%
        dplyr::select(
          SYMBOL, gene_id, ENTREZID
        ) %>% 
        filter(!is.na(ENTREZID))%>%
        mutate(
          kegg_id=paste0("mmu:", ENTREZID)
        ),
      by="kegg_id"
    )  %>% 
    left_join(
      pax6.deg_master %>%
        filter(
          Test == "ExactTest" & Group_1 == "WTE" & Group_2 == "WTF" &
            Filtered == "ribo" & Partition == "Pair"
        ) %>%
        mutate(
          BIO_WTFvsWTE = (
            (Avg1 > 2 | Avg2 > 2) & 
              abs(Avg1 - Avg2) > 2 &
              FDR < 0.05 &
              abs(logFC) > 1
          ),
          ROBUST_WTFvsWTE = (Avg1 > 2 | Avg2 > 2)
        ) %>%
        dplyr::select(
          gene_id,
          WTFvsWTE_logFC = logFC,
          WTFvsWTE_FDR = FDR,
          Cell_WTE_FPKM = Avg1,
          Cell_WTF_FPKM = Avg2,
          ROBUST_WTFvsWTE,
          BIO_WTFvsWTE
        ),
      by="gene_id"
    ) %>%
    left_join(
      pax6.deg_master %>%
        filter(
          Test == "ExactTest" & Group_1 == "P6E" & Group_2 == "P6F" &
            Filtered == "ribo" & Partition == "Pair"
        ) %>%
        mutate(
          BIO_P6FvsP6E = (
            (Avg1 > 2 | Avg2 > 2) & 
              abs(Avg1 - Avg2) > 2 &
              FDR < 0.05 &
              abs(logFC) > 1
          ),
          ROBUST_P6FvsP6E = (Avg1 > 2 | Avg2 > 2)
        ) %>%
        dplyr::select(
          gene_id,
          P6FvsP6E_logFC = logFC,
          P6FvsP6E_FDR = FDR,
          Cell_P6E_FPKM = Avg1,
          Cell_P6F_FPKM = Avg2,
          BIO_P6FvsP6E,
          ROBUST_P6FvsP6E
        ),
      by="gene_id"
    ) %>%
    left_join(
      pax6.deg_master %>%
        filter(
          Test == "ExactTest" & Group_1 == "WTE" & Group_2 == "P6E" &
            Filtered == "ribo" & Partition == "Pair"
        ) %>%
        mutate(
          BIO_P6EvsWTE = (
            (Avg1 > 2 | Avg2 > 2) & 
              abs(Avg1 - Avg2) > 2 &
              FDR < 0.05 &
              abs(logFC) > 1
          ),
          ROBUST_P6EvsWTE = (Avg1 > 2 | Avg2 > 2)
        ) %>%
        dplyr::select(
          gene_id,
          P6EvsWTE_logFC = logFC,
          P6EvsWTE_FDR = FDR,
          Geno_WTE_FPKM = Avg1,
          Geno_P6E_FPKM = Avg2,
          BIO_P6EvsWTE,
          ROBUST_P6EvsWTE
        ),
      by="gene_id"
    ) %>%
    left_join(
      pax6.deg_master %>%
        filter(
          Test == "ExactTest" & Group_1 == "WTF" & Group_2 == "P6F" &
            Filtered == "ribo" & Partition == "Pair"
        )%>%
        mutate(
          BIO_P6FvsWTF = (
            (Avg1 > 2 | Avg2 > 2) & 
              abs(Avg1 - Avg2) > 2 &
              FDR < 0.05 &
              abs(logFC) > 1
          ),
          ROBUST_P6FvsWTF = (Avg1 > 2 | Avg2 > 2)
        ) %>%
        dplyr::select(
          gene_id,
          P6FvsWTF_logFC = logFC,
          P6FvsWTF_FDR = FDR,
          Geno_WTF_FPKM = Avg1,
          Geno_P6F_FPKM = Avg2,
          BIO_P6FvsWTF,
          ROBUST_P6FvsWTF
        ),
      by="gene_id"
    ) %>%
    mutate(
      TRRUST_TARGET = SYMBOL %in% trrust_targets$Target,
      SUN_TARGET = SYMBOL %in% sun_targets$Target,
      MSIG_TARGET = SYMBOL %in% paste0(
        substr(msig_targets$V1, 1,1),       ### NEED BETTER HOMOLOGY MAPPING!
        tolower(
          substr(
            msig_targets$V1, 2, nchar(msig_targets$V1)
          )
        )
      )
    )
  
  kegg_pathway_measures[is.na(kegg_pathway_measures)] <- 0
  print(nrow(kegg_pathway_measures))
  
  fn <- paste0("KEGG_",pathway, "_Network_Node_Measurements.csv")
  file_path <- paste(results, fn, sep="/")
  write.csv(kegg_pathway_measures, file_path, row.names = F, quote=F)
  result_files <- append(result_files, file_path)
  
  fn <- paste0("KEGG_",pathway, "_Network.csv")
  file_path <- paste(results, fn, sep="/")
  write.csv(kegg, file_path, row.names = F, quote = F)
  result_files <- append(result_files, file_path)
}
######################### Advaita Pathway Gene Lists #########################



######################  Push script and data to Synapse ######################

# Add this script to the code dir
syn_script <- synFindEntityId(
  "Pax6_DE_KEGG_Tables.R",
  parent=syn_code_dir
)

if(is.null(syn_script)){
  syn_script <- File(
    path="scripts/Pax6_DE_KEGG_Tables.R",
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
    description = "Extract DE Measures for KEGG Pathways",
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
