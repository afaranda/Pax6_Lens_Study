##############################################################################
#                                                                            #
#  File: Pax6_DE_LIRTS_Analysis_TMM.R                                        #
#  Author: Adam Faranda                                                      #
#  Created: June 18, 2020                                                    #
#  Purpose:                                                                  #
#       Identify Pax6 Dependent Genes in LEC that are also subject to injury #
#       dependent differential expression in the LIRTS data. Identify        #
#       the post surgical interval in the LIRTS that is most similar to Pax6 #
#       haploinsufficiency. Identify Pax6 Targets that are differentially    #
#       expressed during mechanical injury, and which of those are also      #
#       differentially expressed in Pax6 haploinsufficient lenses            #
#                                                                            #
##############################################################################

################## Load Libraries and Source Dependencies ####################
library(edgeR)
library(tidyr)
library(tibble)
library(dplyr)
library(tibble)
#library(synapser)
library(ggrepel)
library(EnhancedVolcano)
options(echo=T)

# Enter Working Directory and Load Raw Data
setwd('/Users/adam/Documents/23Feb2022_Pax6_Study_DEG_Analysis/')
source('scripts/Overlap_Comparison_Functions.R')
source('scripts/Wrap_edgeR_Functions.R')
source('scripts/Excel_Write_Functions.R')
source('scripts/synapse_reticulate_wrapper.R')

wd<-getwd()
results<-paste(wd,'results/iSyTE_Analysis',sep='/')
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
syn_deg_sps_dir <- synFindEntityId("Pax6_iSyTE_Analysis", parent = syn_project)

if(is.null(syn_deg_sps_dir)){
  folder <- Folder(
    name="Pax6_iSyTE_Analysis",
    parent=syn_project
  )
  syn_deg_sps_dir <- synStore(folder)
} 

## Check for the DEG_Tables folder
syn_deg_tbl_dir <- synFindEntityId(
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

############################# Load iSyTE Data ################################

syn_isyte <- synFindEntityId(
  "isyte528_long_table.csv",
  parent=syn_data_dir
)

if(file.exists('data/isyte528_long_table.csv')){
  isyte528 <- read.csv('data/isyte528_long_table.csv', row.names = 1)
} else if(!is.null(syn_isyte)){
  synGet(
    syn_isyte, downloadLocation=data_dir
  )
  isyte528 <- read.csv('data/isyte528_long_table.csv', row.names = 1)
}else{
  stop("Run Prepare_Expression_DGELists.R first")
}

## Deduplicate iSyTE Table -- this does not impact 
## other spreadsheets. 
isyte528 %>%
  group_by(Platform, Interval, MGI.symbol) %>%
  filter(
    (
      Lens_Expression + WEB_Expression == max(
        Lens_Expression + WEB_Expression
      )
    )
  ) %>% group_by() %>%
  left_join(
    pax6.master$genes %>%
      dplyr::select(SYMBOL, gene_id),
    by=c(MGI.symbol = "SYMBOL")
  ) -> isyte528

isyte_528_P56 <- isyte528 %>% 
  filter(
    !is.na(gene_id) & Interval == "P56" & p_value < 0.05 & fold_change > 2
  ) %>%
  pull("gene_id")


## Final list of existing data files used by this script
used_files <- list(
  syn_dgelists,
  syn_isyte,
  syn_deg_master
)


#####################  Setup group and contrast iterators ####################

### Compare the Pax6 vs WT LEC to the pairwise exact test for
### each post surgical interval and report:
###
## Analyses by study
###   1.) Contingency tables for Overall Differential Expression
###   2.) Contingency tables for directional dependence
###   3.) Fold Change Plot.

isyte_enrichment <- function(
  pax6_deg = data.frame(),
  isyte = data.frame(),
  result_label = "DNA1_6_Hours_Post_Injury"
){

  ### Run the Fisher's Exact Test for Overrepresentation
  universe <- union(
    pax6_deg$gene_id,
    isyte$gene_id
  )
  
  ## Join Tables for uniform analysis and add indicator columns
  test_data <- data.frame(
    gene_id = universe
  ) %>%
    left_join(
      pax6_deg %>%
        select(
          gene_id, pax6_logFC = logFC, pax6_FDR = FDR,
          pax6_Avg1=Avg1, pax6_Avg2=Avg2
        ),
      by="gene_id"
    ) %>%
    left_join(
      isyte %>%
        select(
          gene_id, isyte_fold_change=fold_change, isyte_pvalue= p_value,
          isyte_Lens_Expression = Lens_Expression, 
          isyte_WEB_Expression=WEB_Expression,
          Platform, Interval
        ),
      by="gene_id"
    ) %>%
    replace_na(replace=list(rep(0,11))) %>%
    mutate(
      Contrast = result_label,
      IS_PAX6 = (
        abs(pax6_logFC) > 1 & pax6_FDR < 0.05 &
          (pax6_Avg1 > 2 | pax6_Avg2 > 2) &
          (abs(pax6_Avg1 - pax6_Avg2) > 2) &
          !is.na(pax6_logFC)
      ),
      UP_PAX6 = (
        pax6_logFC > 1 & pax6_FDR < 0.05 &
          (pax6_Avg1 > 2 | pax6_Avg2 > 2) &
          (abs(pax6_Avg1 - pax6_Avg2) > 2) &
          !is.na(pax6_logFC)
      ),
      DN_PAX6 = (
        pax6_logFC < -1 & pax6_FDR < 0.05 &
          (pax6_Avg1 > 2 | pax6_Avg2 > 2) &
          (abs(pax6_Avg1 - pax6_Avg2) > 2) &
          !is.na(pax6_logFC)
      ),
      IS_ISYTE = (
        !is.na(isyte_fold_change)
      ),
      IS_ISYTE_DE = (
        abs(isyte_fold_change) > 2 & isyte_pvalue < 0.05 & 
          !is.na(isyte_fold_change)
      ),
      UP_ISYTE = (
        isyte_fold_change > 2 & isyte_pvalue < 0.05  & 
          !is.na(isyte_fold_change)
      )
    ) 
  
  print(nrow(test_data))
  ## Prepare Overall Enrichment Table
  print("Overall overrepresentation -- all 528 iSyTE")
  contingency <- table(
    IS_ISYTE=test_data$IS_ISYTE,
    IS_PAX6=test_data$IS_PAX6
  )
  print(contingency)
  print(fisher.test(contingency))
  
  print("Overall overrepresentation -- Lens Enriched at P56")
  contingency <- table(
    UP_ISYTE =test_data$UP_ISYTE,
    IS_PAX6=test_data$IS_PAX6
  )
  print(contingency)
  print(fisher.test(contingency))
  
  print("Directional Dependence")
  contingency <- table(
    test_data %>%
      filter(IS_PAX6 & IS_ISYTE) %>%
      select(UP_ISYTE, UP_PAX6)
  )
  print(sum(contingency))
  print(contingency)
  print(fisher.test(contingency))
  
  return(
    test_data %>%
      filter(IS_PAX6 & IS_ISYTE)
  )
}


# Define A list of named contrasts; each element points to a vector with
# a pair of group labels. Positive fold changes will be associated
# with the second group listed. 
contrasts=list(
  Pax6vsWT_LEC_X_P56=c('WTE', 'P6E', 'P56'),
  Pax6vsWT_LFC_X_P56=c('WTF', 'P6F', 'P56'),
  WT_LECvsLFC_X_P56=c('WTE', 'WTF', 'P56'),
  P6_LECvsLFC_X_P56=c('P6E', 'P6F', 'P56')
)

########### Calculate Fisher's Exact Test results for each contrast  #########
result_files <- c()

pax6_isyte_deg_table <- data.frame()
fn <- "iSyTE_Enrichment_In_Pax6_DEG.txt"
path <- paste(results, fn, sep="/")
result_files <- append(result_files, path)
sink(path)
for(c in names(contrasts)){
  pax6_deg <- pax6.deg_master %>% filter(
    Partition == "Pair",
    Filtered == "ribo",
    Group_1 == contrasts[[c]][1],
    Group_2 == contrasts[[c]][2],
    Test == "ExactTest"
  )%>%
    left_join(
      pax6.master$genes %>%
        select(
          gene_id, SYMBOL, DESCRIPTION,
          IS_ISYTE_P56, IS_ZONULE, IS_TRRUST_PAX6_TARGET,
          IS_SUN_PAX6_TARGET, IS_SUN_PAX6_LENS_PEAK,
          IS_SUN_PAX6_FOREBRAIN_PEAK
        ),
      by="gene_id"
    ) 
  print(paste("All ", nrow(pax6_deg)))
  pax6_deg <- union(
    pax6_deg %>% 
      filter(SYMBOL !="") %>%
      group_by(SYMBOL) %>%
      filter((Avg1 + Avg2) == max(Avg1 + Avg2)),
    pax6_deg %>% filter(SYMBOL =="")
  ) %>% group_by() %>% as.data.frame()
  print(paste("Deduped ", nrow(pax6_deg)))
  
  isyte <- isyte528 %>%
    filter(
      Platform == "affy430",
      Interval == contrasts[[c]][3]
    )
  print(c)
  pax6_isyte_deg_table <-bind_rows(
    pax6_isyte_deg_table,
    isyte_enrichment(pax6_deg, isyte, result_label=c)
  )
}
sink()

fn <- "ISYTE_GENES_In_Pax6_Long_Table.csv"
path <- paste(results, fn, sep="/")
write.csv(pax6_isyte_deg_table, path)
result_files <- append(result_files, path)


pax6.deg_master<- pax6.deg_master %>%
  mutate(
    is_p56 = gene_id %in% isyte_528_P56,
    is_bio = abs(Avg1 - Avg2) > 2 & (Avg1 > 2 | Avg2 > 2)
  )

fn <- "Pax6_iSyTE_DEG_Counts.csv"
path <- paste(results, fn, sep="/")
result_files <- append(result_files, path)
pax6.deg_master %>% 
  filter(Filtered == "ribo" & Partition == "Pair") %>%
  group_by(Test, Group_1, Group_2) %>%
  summarise(
    Total_OBS = n() ,#length(unique(hs_symbol)),
    Total_DEG = sum(abs(logFC)>1 & FDR < 0.05 & is_bio),
    DEG_In_p56 = sum(is_p56 & abs(logFC)>1 & FDR < 0.05 & is_bio),
    UP_In_p56 = sum(is_p56 & (logFC > 1 & FDR < 0.05 & is_bio)),
    DOWN_In_p56= sum(is_p56 & (logFC < -1 & FDR < 0.05 & is_bio)),
    NON_DEG_In_p56=sum(is_p56 & !(abs(logFC)>1 & FDR < 0.05 & is_bio)),
    NO_OBS_In_p56=length(setdiff(isyte_528_P56, gene_id))
  ) %>% write.csv(path)

############ Generate Enhanced Volcano Plot Highlighting iSyTE DEG #############

pax6_deg <- pax6.deg_master %>%
  filter(Filtered == "ribo" & Partition == "Pair" & Test == "ExactTest") %>%
  filter(Group_1 == "WTE" & Group_2 == "P6E") %>%
  filter(Avg1 >2 | Avg2 > 2) %>%
  inner_join(
    pax6.master$genes %>%
      select(gene_id, SYMBOL),
    by="gene_id"
  ) %>% mutate(
    IS_ISYTE = (gene_id %in% isyte_528_P56 ) & ( abs(logFC) > 1) & (FDR < 0.05)
  ) %>% filter(IS_ISYTE)

fn <- "ISYTE_GENES_DE_IN_PAX6_LEC.jpg"
path <- paste(results, fn, sep="/")
result_files <- append(result_files, path)

p <- EnhancedVolcano(
  pax6_deg %>% filter(IS_ISYTE),
  x = "logFC", y="FDR", lab=pax6_deg$SYMBOL,
  selectLab = pax6_deg %>% filter(IS_ISYTE) %>% pull(SYMBOL),
  drawConnectors = T, max.overlaps = 1000,
  pCutoff = 5*10^-2, 
  title = "Lens Preferred Genes influenced by Pax6 haploinsufficiency in LEC",
  subtitle = ""
)
ggsave(path, plot=p, width=12, height=7)

pax6_deg <- pax6.deg_master %>%
  filter(Filtered == "ribo" & Partition == "Pair" & Test == "ExactTest") %>%
  filter(Group_1 == "WTF" & Group_2 == "P6F") %>%
  filter(Avg1 >2 | Avg2 > 2) %>%
  inner_join(
    pax6.master$genes %>%
      select(gene_id, SYMBOL),
    by="gene_id"
  ) %>% mutate(
    IS_ISYTE = (gene_id %in% isyte_528_P56 ) & ( abs(logFC) > 1) & (FDR < 0.05)
  ) %>% filter(IS_ISYTE)

fn <- "ISYTE_GENES_DE_IN_PAX6_LFC.jpg"
path <- paste(results, fn, sep="/")
result_files <- append(result_files, path)

p <- EnhancedVolcano(
  pax6_deg %>% filter(IS_ISYTE),
  x = "logFC", y="FDR", lab=pax6_deg$SYMBOL,
  selectLab = pax6_deg %>% filter(IS_ISYTE) %>% pull(SYMBOL),
  drawConnectors = T, max.overlaps = 1000,
  pCutoff = 5*10^-2, 
  title = "Lens Preferred Genes influenced by Pax6 haploinsufficiency in LFC",
  subtitle = ""
)
ggsave(path, plot=p, width=12, height=7)


######################  Push script and data to Synapse ######################

# Add this script to the code dir
syn_script <- synFindEntityId(
  "Pax6_DE_iSyTE_Analysis.R",
  parent=syn_code_dir
)

if(is.null(syn_script)){
  syn_script <- File(
    path="scripts/Pax6_DE_iSyTE_Analysis.R",
    parent=syn_code_dir
  )

  syn_script <- synStore(
    syn_script
  )
}

synSetProvenance(
  syn_script,
  activity = Activity(
    name = "Pax6_iSyTE_DEG",
    description = "Lens preferred Genes that are Pax6 Dependent",
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
