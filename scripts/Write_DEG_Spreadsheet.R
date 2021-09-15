##############################################################################
#                                                                            #
#  File: Pairwise_Tests.R                                                    #
#  Author: Adam Faranda                                                      #
#  Created: June 18, 2020                                                    #
#  Purpose: Write Pairwise analysis results to an excel spreadsheet          #
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
setwd('/Users/adam/Documents/11Sep2021_Pax6_Study_DEG_Analysis/')
source('scripts/Excel_Write_Functions.R')
source('scripts/Overlap_Comparison_Functions.R')
wd<-getwd()
results<-paste(wd,'results',sep='/')
data_dir <- paste0(wd,"/data")         ## Local data directory

######################### Setup Synapse Connection ###########################
synLogin()
syn_project <- synFindEntityId("Pax6_Happloinsuficiency_In_The_Lens")
syn_code_dir <- synFindEntityId("code", parent=syn_project)
syn_data_dir <- synFindEntityId("data", parent=syn_project)
syn_fig_dir <- synFindEntityId("Figures", parent=syn_project)
syn_deg_tbl_dir <- synFindEntityId("DEG_Tables", parent = syn_project)
syn_deg_sps_dir <- synFindEntityId("DEG_Spreadsheets", parent = syn_project)


## Check for the DEG_Spreadsheets_Pairwise_Analysis folder
syn_deg_tbl_dir <- synapser::synFindEntityId(
  "DEG_Tables",
  syn_project
)

syn_deg_sps_dir <- synapser::synFindEntityId(
  "DEG_Spreadsheets",
  syn_project
)

syn_deg_master <- synFindEntityId(
  "pax6_deg_tables.Rdata", 
  parent=syn_deg_tbl_dir
)

# Add this script to the code dir
syn_push <- synFindEntityId(
  "Write_DEG_Spreadsheet.R",
  parent=syn_code_dir
)
if(is.null(syn_push)){
  syn_push <- File(
    path="scripts/Write_DEG_Spreadsheet.R",
    parent=syn_code_dir
  )
  
  syn_push <- synStore(
    syn_push
  )
}

# Add this bias plot script to the code dir
syn_bias <- synFindEntityId(
  "Plot_Length_GC_Biases.R",
  parent=syn_code_dir
)
if(is.null(syn_bias)){
  syn_bias <- File(
    path="scripts/Plot_Length_GC_Biases.R",
    parent=syn_code_dir
  )
  
  syn_bias <- synStore(
    syn_bias
  )
}


## Activity definitions for spreadsheet and table provenance
syn_act_sps <- Activity(
  name="upload_deg_spreadsheets",
  description="upload differential expression spreadsheets"
)
syn_act_sps$executed(syn_push)

syn_act_tbl <- Activity(
  name="upload_deg_tables",
  description="upload differential expression tables"
)
syn_act_tbl$executed(syn_push)

syn_act_bpl <- Activity(
  name="generate_bias_plot",
  description="create bias plot for differential expression spreadsheets"
)
syn_act_bpl$executed(syn_bias)

# Fetch files currently stored in DEG Spreadsheets folder
if(!is.null(syn_deg_sps_dir)){
  deg_folder_chl <- as.list(synGetChildren(syn_deg_sps_dir))
  deg_spreadsheets <- sapply(deg_folder_chl, function(s) s[[2]])
  names(deg_spreadsheets) <- sapply(deg_folder_chl, function(s) s[[1]])
}
print(deg_spreadsheets)


########################## Load in Master DGEList ############################
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
  source("scripts/Prepare_Expression_DGELists.R")
  syn_dgelists <- synFindEntityId(
    "pax6_master_dgelists.Rdata", 
    parent=syn_data_dir
  )
}

## Get Ribo Filtered Counts
pax6.master$samples$group <- factor(
  gsub("_ribo_[123]","", pax6.master$samples$label),
  levels=c("WTE", "WTF", "P6E", "P6F")
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
  source("scripts/Generate_Pax6_DEG_Tables.R")
  syn_deg_master <- synFindEntityId(
    "pax6_deg_tables.Rdata", 
    parent=syn_deg_tbl_dir
  )
}

#####################  Setup group and contrast iterators ####################


# Define A list of named contrasts; each element points to a vector with
# a pair of group labels. Positive fold changes will be associated
# with the second group listed. 
contrasts=list(
  WT_Fibers_vs_Epi=c('WTE', 'WTF'),
  P6_Fibers_vs_Epi=c('P6E', 'P6F'),
  Epi_P6_vs_WT=c('WTE', 'P6E'),
  Fibers_P6_vs_WT=c('WTF', 'P6F'),
  `2Way_Fib_vs_Epi`=c("Epi", "Fib"),
  `2Way_P6_vs_WT`=c("WT", "P6"),
  `2Way_Interaction`=c("REF", "INX"),
  DegNorm_WT_Fibers_vs_Epi=c('WTE', 'WTF'),
  DegNorm_P6_Fibers_vs_Epi=c('P6E', 'P6F'),
  DegNorm_Epi_P6_vs_WT=c('WTE', 'P6E'),
  DegNorm_Fibers_P6_vs_WT=c('WTF', 'P6F'),
  DegNorm_2Way_Fib_vs_Epi=c("Epi", "Fib"),
  DegNorm_2Way_P6_vs_WT=c("WT", "P6"),
  DegNorm_2Way_Interaction=c("REF", "INX")
)

# Define A list of named contrasts; each element points to a string which
# briefly describes the contrast
contrast_descriptions<-list(
  WT_Fibers_vs_Epi="Wildtype lens fiber cells vs lens epithelial cells",
  P6_Fibers_vs_Epi="Sey lens fiber cells vs lens epithelial cells",
  Epi_P6_vs_WT="Sey lens epithelial cells vs lens epithelial cells",
  Fibers_P6_vs_WT="Sey lens fiber cells vs lens fiber cells",
  `2Way_Fib_vs_Epi`="2 Way interaction model lens fiber cells vs lens fiber cells",
  `2Way_P6_vs_WT`="2 Way interaction model Sey vs Wildtype",
  `2Way_Interaction`="2 Way interaction model haploinsufficiency and differentiation",
  DegNorm_WT_Fibers_vs_Epi="Wildtype lens fiber cells vs lens epithelial cells",
  DegNorm_P6_Fibers_vs_Epi="Sey lens fiber cells vs lens epithelial cells",
  DegNorm_Epi_P6_vs_WT="Sey lens epithelial cells vs lens epithelial cells",
  DegNorm_Fibers_P6_vs_WT="Sey lens fiber cells vs lens fiber cells",
  DegNorm_2Way_Fib_vs_Epi="2 Way interaction model lens fiber cells vs lens fiber cells",
  DegNorm_2Way_P6_vs_WT="2 Way interaction model Sey vs Wildtype",
  DegNorm_2Way_Interaction="2 Way interaction model haploinsufficiency and differentiation"
)

contrast_match_bias<-list(
  WT_Fibers_vs_Epi="P",
  P6_Fibers_vs_Epi="Sey lens fiber cells vs lens epithelial cells",
  Epi_P6_vs_WT="Sey lens epithelial cells vs lens epithelial cells",
  Fibers_P6_vs_WT="Sey lens fiber cells vs lens fiber cells",
  `2Way_Fib_vs_Epi`="2 Way interaction model lens fiber cells vs lens fiber cells",
  `2Way_P6_vs_WT`="2 Way interaction model Sey vs Wildtype",
  `2Way_Interaction`="2 Way interaction model haploinsufficiency and differentiation",
  DegNorm_WT_Fibers_vs_Epi="Wildtype lens fiber cells vs lens epithelial cells",
  DegNorm_P6_Fibers_vs_Epi="Sey lens fiber cells vs lens epithelial cells",
  DegNorm_Epi_P6_vs_WT="Sey lens epithelial cells vs lens epithelial cells",
  DegNorm_Fibers_P6_vs_WT="Sey lens fiber cells vs lens fiber cells",
  DegNorm_2Way_Fib_vs_Epi="2 Way interaction model lens fiber cells vs lens fiber cells",
  DegNorm_2Way_P6_vs_WT="2 Way interaction model Sey vs Wildtype",
  DegNorm_2Way_Interaction="2 Way interaction model haploinsufficiency and differentiation"
)

# Fetch files currently stored in DEG Spreadsheets folder
if(!is.null(syn_deg_tbl_dir)){
  deg_folder_chl <- as.list(synGetChildren(syn_deg_tbl_dir))
  deg_spreadsheets <- sapply(deg_folder_chl, function(s) s[[2]])
  names(deg_spreadsheets) <- sapply(deg_folder_chl, function(s) s[[1]])
}
print(deg_spreadsheets)

used_files <- list(
  syn_dgelists,
  syn_deg_master
)
        
for(c in names(contrasts)){
  if(grepl("DegNorm",c)){
    master <- pax6.master.dgn
    deg <- pax6.deg_master %>% 
      filter(Filtered == "degnorm")
  } else {
    master <- pax6.master
    deg <- pax6.deg_master %>% 
      filter(Filtered == "ribo")
  }
  
  
  ## Get labels for groups in this contrast
  g1<-contrasts[[c]][1]
  g2<-contrasts[[c]][2]
  
  ## Get samples assigned to groups in this contrast
  if(!grepl("2Way", c)){
    gr<-c(
      master$samples %>% filter(group==g1) %>%pull("label"), 
      master$samples %>% filter(group==g2) %>%pull("label")
    )
    
    ## Subset DEG Table from pax6.deg_master
    degSet <- deg %>%
      filter(
        Partition=="Pair" &
          Test=="ExactTest" &
          Group_1 == g1 &
          Group_2 == g2
      ) %>%
      inner_join(
        master$genes,
        by="gene_id"
      )
    print(nrow(degSet))
    
    ## Get names of length and GC Bias Plots
    gcfn<-paste0(
      "results/LGB_ribo_Pair_ExactTest_",g2,"v",g1,
      "_Stat_GC_Bias_By_Partition.jpg"
    )
    lnfn<-paste0(
      "results/LGB_ribo_Pair_ExactTest_",g2,"v",g1,
      "_Stat_Length_Bias_By_Partition.jpg"
    )
    
    ## Get template for building this spreadsheet
    template<-paste0("scripts/",c,"_Pairwise_Exact_Test_deg_template.xlsx")
    print(template)
    print(file.exists(template))
  } else {
    gr <- colnames(master)
    ## Get names of length and GC Bias Plots
    gcfn<-paste0(
      "results/LGB_ribo_2Way_QLFTest_",g2,"v",g1,
      "_Stat_GC_Bias_By_Partition.jpg"
    )
    lnfn<-paste0(
      "results/LGB_ribo_2Way_QLFTest_",g2,"v",g1,
      "_Stat_Length_Bias_By_Partition.jpg"
    )
    
    ## Subset DEG Table from pax6.deg_master
    degSet <- deg %>%
      filter(
        Partition=="2Way" &
          Test=="QLFTest" &
          Group_1 == g1 &
          Group_2 == g2
      ) %>%
      inner_join(
        master$genes,
        by="gene_id"
      )
    
    print(nrow(degSet))
    
    ## Get template for building this spreadsheet
    template<-paste0("scripts/",c,"_QLF_Interaction_Model_deg_template.xlsx")
    print(template)
    print(file.exists(template))
  }   
 
  if(!file.exists(gcfn)){
    source("scripts/Plot_Length_GC_Biases.R")
  }
 
  ## Get file name for spreadsheet
  fn<-paste0("WDS_",c,"_Differential_Expression.xlsx")
  
  ## Check for existing spreadsheet in Synapse and retrieve, or build a new
  ## spreadsheet if none exists.
  if(fn %in% names(deg_spreadsheets)){
    print(
      paste("Fetching spreadsheet:", fn, "from synapse")
    )
    
    #### synGet can't recognize named vectors This is a stupid band-aid
    spn <- deg_spreadsheets[fn]
    names(spn) <- NULL
    synGet(
      spn,
      downloadLocation=results
    )
  } else {
    print(
      paste("Building spreadsheet:", fn)
    )
    
    # concatenate full path for DEG spreadsheet
    fn<-paste(results, fn, sep="/")
    print(fn)
    
    # Name headers for columns storing group average FPKM estimates.
    names(degSet) <- gsub("Avg1", paste0(g1, "_Avg_FPKM"), names(degSet))
    names(degSet) <- gsub("Avg2", paste0(g2, "_Avg_FPKM"),names(degSet))
    
    Avg1 <- paste(g1, 'Avg_FPKM', sep='_')
    Avg2 <- paste(g2, 'Avg_FPKM', sep='_')
    
    # Specify columns passed to spreadsheet creator
    cols = c(
      'gene_id', 'SYMBOL', 'DESCRIPTION','SEQNAME', 'logFC', 'PValue',
      'FDR', Avg1, Avg2, 'Group_1', 'Group_2'
    )
    print(cols)
    
    createDEGSpreadSheet(
      C1 = c,                          # Name of the contrast
      dg1 = degSet,                    # Data Set for the contrast
      dg1.bioFun = bioSigRNASeq,       # Biological significance filter for dg1
      dg1.fdr = "FDR",                 # Statistic used to filter genes for dg1
      dg1.lfc = "logFC",               # Column in dg1 with log Fold Changes
      dg1.Avg1 = Avg1,                 # Column in dg1 with average value for Group_1
      dg1.Avg2 = Avg2,                 # Column in dg1 with average value for Group_2
      dg1.me = 2,                      # Min. expression for dg1.bioFun
      dg1.x = 23,                      # row number, corner of dg1 Summary table
      dg1.y = 2,                       # col number, corner of dg1 Summary table
      dg1.ds = contrast_descriptions[[c]], # short description for contrast C1 (dg1)
      template = "scripts/deg_template.xlsx",
      descPageName="Data Description", # Name of sheet to write summary tables
      wb = NULL,                       # Optionally pass a workbook object instead.
      pref = "" ,                      # Prefix for output file.
      fname=fn,                        # Manually specify an output file name
      use_lfc = FALSE,                 # Whether to use logFC or Fold_Change
      cols=setdiff(                    # Names of columns to keep in final tables
        cols,
        c('Group_1', 'Group_2')
      ),
      sc_cols=c("PValue", "FDR"),
      extraTables = list(
        Raw_Counts=master[,gr]$counts %>% 
          as.data.frame() %>%
          rownames_to_column("gene_id") %>%
          inner_join(
            master$genes %>% 
              select(
                -pl_length, -pl_gc, -is_principal,
                -GENEBIOTYPE, -tx_id, -gene_id_version
              ), by="gene_id"
          ) %>% 
          relocate(
            gene_id, SYMBOL, DESCRIPTION, SEQNAME, eu_length, eu_gc
          )
      ),
      tableNames=c(                    # Names to use for standard tables
        'All Present Genes', 
        'Statistically Significant',
        'Biologically Significant'
      ),
      img_set=list(
        list(
          fn=gcfn, sr=48, sc=2
        ),
        list(
          fn=lnfn, sr=48, sc=6
        )
      )  
    )
    ## Push bias plots to Synapse
    syn_gc_bias <- File(
      path=gcfn,
      parent =syn_fig_dir
    )
    syn_gc_bias <- synStore(
      syn_gc_bias
    )
    synSetProvenance(
      syn_gc_bias,
      syn_act_bpl
    )
    
    # Add GC Bias plot to list of files used by this script
    used_files <- append(
      used_files,
      syn_gc_bias
    )
    
    syn_ln_bias <- File(
      path=lnfn,
      parent =syn_fig_dir
    )
    syn_ln_bias <- synStore(
      syn_ln_bias
    )
    synSetProvenance(
      syn_ln_bias,
      syn_act_bpl
    )
    
    # Add length Bias plot to list of files used by this script
    used_files <- append(
      used_files,
      syn_ln_bias
    )
    
    ## Push new spreadsheet to Synapse
    syn_spreadsheet <- File(
      path=fn,
      parent=syn_deg_sps_dir
    )
    
    syn_spreadsheet <- synStore(
      syn_spreadsheet
    )
    
    synSetProvenance(
      syn_spreadsheet,
      syn_act_sps
    )
    
    ## Write degSet to a tab separated table for Advaita
    fn <- gsub(
      "_Differential_Expression.xlsx",
      "_DEG_Table.tsv", fn
    )
    write.table(
      degSet, file=fn, row.names = F,
      quote = F, sep="\t"
    )
    
    ## Push DEG table to Synapse
    syn_table <- File(
      path=fn,
      parent=syn_deg_tbl_dir
    )
    syn_table <- synStore(
      syn_table
    )
    synSetProvenance(
      syn_table,
      syn_act_tbl
    )
  }
}

######################  Push script and data to Synapse ######################
synSetProvenance(
  syn_push,
  activity = Activity(
    name = "wrote_deg_spreadsheets",
    description = "collated plots and DEG tables into spreadsheets",
    used=used_files
  )
)

