##############################################################################
#                                                                            #
#  File: Pax6_DE_Zonule_Analysis_TMM.R                                       #
#  Author: Adam Faranda                                                      #
#  Created: June 18, 2020                                                    #
#  Purpose:  Identify Zonular proteins differentially regulated in Pax6      #
#            haploinsufficient lenses compared to wildtype                   #
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
results<-paste(wd,'results/Zonule_Analysis',sep='/')
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
syn_deg_sps_dir <- synFindEntityId("Pax6_Zonule_Analysis", parent = syn_project)
syn_hpco_gene_meta <- synFindEntityId(
  "Gene_Metadata",
  synFindEntityId("Human_PCO_Study_Data_Set_1")
)

if(is.null(syn_deg_sps_dir)){
  folder <- Folder(
    name="Pax6_Zonule_Analysis",
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

syn_zonules <- synFindEntityId(
  "Zonule_Proteins_Mapped_Ensembl_ID.txt",
  parent=syn_hpco_gene_meta
)

syn_homology <- synFindEntityId(
  "human_mouse_homology_groups.csv",
  parent=syn_hpco_gene_meta
)

if(file.exists('data/Zonule_Proteins_Mapped_Ensembl_ID.txt')){
  zonules <- read.table(
    'data/Zonule_Proteins_Mapped_Ensembl_ID.txt',
    header=T, quote="", sep="\t"
  ) %>%
    distinct(
      Gene.stable.ID, Gene.name
    )
} else if(!is.null(syn_zonules)){
  synGet(
    syn_zonules, downloadLocation=data_dir
  )
  zonules <- read.table(
    'data/Zonule_Proteins_Mapped_Ensembl_ID.txt',
    header=T, quote="", sep="\t"
  ) %>%
    distinct(
      Gene.stable.ID, Gene.name
    )
}else{
  stop("Run Prepare_Expression_DGELists.R first")
}

syn_homology <- synFindEntityId(
  "human_mouse_homology_groups.csv",
  parent=syn_hpco_gene_meta
)

if(file.exists('data/human_mouse_homology_groups.csv')){
  homology <- read.csv(
    'data/human_mouse_homology_groups.csv', row.names = 1
  ) 
} else if(!is.null(syn_homology)){
  synGet(
    syn_homology, downloadLocation=data_dir
  )
  homology <- read.csv(
    'data/human_mouse_homology_groups.csv', row.names = 1
  ) 
}else{
  stop("Run Prepare_Expression_DGELists.R first")
}

zonules <- inner_join(
  zonules,
  homology,
  by="Gene.stable.ID"
)

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

### For Combined Targets, Use the observation from Sun 2015 
combined_targets <- bind_rows(
  sun_targets %>%
    mutate(
      Reference=ifelse(
        Target %in% intersect(
          sun_targets$Target, 
          trrust_targets$Target
        ),
        "Both","Sun_2015"
      ) 
    ) %>%
    select(
      Reference,
      Target,
      Direction
    ),
  trrust_targets %>%
    mutate(
      Reference="TRRUST"
    ) %>%
    select(
      Reference,
      Target,
      Direction
    )%>% filter(
      !Target %in% intersect(sun_targets$Target, trrust_targets$Target)
    )
)

## Final list of existing data files used by this script
used_files <- list(
  syn_sun_targets,
  syn_trrust_targets,
  syn_dgelists,
  syn_zonules,
  syn_homology,
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

zonule_enrichment <- function(
  pax6_deg = data.frame(),
  zonules = data.frame(),
  homology = data.frame(),
  result_label = "DNA1_6_Hours_Post_Injury"
){
  
  ### Run the Fisher's Exact Test for Overrepresentation
  
  ### Define the universe as the set of all genes included in
  ### both the Pax6 DEG Table, and the Mouse genes encoding Zonular
  ### Proteins. 
  
  universe <- union(
    pax6_deg$gene_id,
    zonules$Mouse.gene.stable.ID
  )
  
  ## Join Tables for uniform analysis and add indicator columns
  ## Total number of records is the total number of included genes
  test_data <- data.frame(
    gene_id = universe
  ) %>%
    left_join(
      pax6_deg, 
      by = "gene_id"
    ) %>%
    mutate(
      IS_PAX6 = (
        abs(logFC) > 1 & FDR < 0.05 &
          (Avg1 > 2 | Avg2 > 2) &
          (abs(Avg1 - Avg2) > 2) &
          !is.na(logFC)
      ),
      UP_PAX6 = (
        logFC > 1 & FDR < 0.05 &
          (Avg1 > 2 | Avg2 > 2) &
          (abs(Avg1 - Avg2) > 2) &
          !is.na(logFC)
      ),
      DN_PAX6 = (
        logFC < -1 & FDR < 0.05 &
          (Avg1 > 2 | Avg2 > 2) &
          (abs(Avg1 - Avg2) > 2) &
          !is.na(logFC)
      )
    ) %>% left_join(
      homology %>%
        distinct(Mouse.gene.stable.ID, Hom_Group),
      by=c(gene_id="Mouse.gene.stable.ID")
    ) %>%
    mutate(
      IS_ZONULE = Hom_Group %in% zonules$Hom_Group,
      Hom_Group = ifelse(is.na(Hom_Group), 0, Hom_Group)
    ) %>%
    left_join(
      pax6.master$genes %>%
        select(
          gene_id, SYMBOL, DESCRIPTION
        ),
      by="gene_id"
    )
  
  print(nrow(test_data))
  ## Prepare Overall Enrichment Table
  print("Overall overrepresentation")
  contingency <- test_data %>%
    group_by(Hom_Group) %>%
    summarize(
      IS_ZONULE = any(IS_ZONULE),
      IS_PAX6 = any(IS_PAX6)
    ) %>% select(IS_ZONULE, IS_PAX6) %>% table()
  print(contingency)
  print(fisher.test(contingency))
  
  print("Overrepresentation -- Upregulated")
  contingency <- test_data %>%
    group_by(Hom_Group) %>%
    summarize(
      IS_ZONULE = any(IS_ZONULE),
      UP_PAX6 = any(UP_PAX6)
    ) %>% select(IS_ZONULE,UP_PAX6) %>% table()
  print(contingency)
  print(fisher.test(contingency))
  
  print("Overrepresentation -- Downregulated")
  contingency <- test_data %>%
    group_by(Hom_Group) %>%
    summarize(
      IS_ZONULE = any(IS_ZONULE),
      DN_PAX6 = any(DN_PAX6)
    ) %>% select(IS_ZONULE,DN_PAX6,) %>% table()
  print(contingency)
  print(fisher.test(contingency))
  
  return(
    test_data %>%
      left_join(
        zonules %>%
          select(Hom_Group, Human_SYMBOL = Gene.name),
        by = "Hom_Group"
      ) %>% distinct() %>%
      group_by(Hom_Group) %>%
      filter(IS_ZONULE & any(!is.na(logFC))) %>%
      mutate(
        PAX6_TARGET = SYMBOL %in% combined_targets$Target
      )
  )
}

zonules %>% 
  group_by(Hom_Group) %>%
  filter(n() > 1)

# Define A list of named contrasts; each element points to a vector with
# a pair of group labels. Positive fold changes will be associated
# with the second group listed. 
contrasts=list(
  Pax6vsWT_LEC_X_ZON=c('WTE', 'P6E'),
  Pax6vsWT_LFC_X_ZON=c('WTF', 'P6F'),
  WT_LECvsLFC_X_ZON=c('WTE', 'WTF'),
  P6_LECvsLFC_X_ZON=c('P6E', 'P6F')
)

########### Calculate Fisher's Exact Test results for each contrast  #########
result_files <- c()

pax6_zonule_deg_table <- data.frame()
fn <- "Zonule_Enrichment_In_Pax6_DEG.txt"
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
  )
  print(c)
  pax6_zonule_deg_table <-bind_rows(
    pax6_zonule_deg_table,
   x <- zonule_enrichment(
     pax6_deg, zonules, homology=homology, result_label=c
   )
  )
}
sink()

fn <- "ZONULE_GENES_In_Pax6_Long_Table.csv"
path <- paste(results, fn, sep="/")
write.csv(pax6_zonule_deg_table, path)
result_files <- append(result_files, path)

pax6.deg_master <- pax6.deg_master %>% left_join(
  homology %>%
    distinct(Mouse.gene.stable.ID, Hom_Group),
  by=c(gene_id="Mouse.gene.stable.ID")
) %>%
  mutate(
    IS_ZONULE = Hom_Group %in% zonules$Hom_Group,
    is_bio = abs(Avg1 - Avg2) > 2 & (Avg1 > 2 | Avg2 > 2)
  )

fn <- "Pax6_Zonule_DEG_Counts.csv"
path <- paste(results, fn, sep="/")
# pax6.deg_master %>%
#   filter(Filtered == "ribo" & Partition == "Pair") %>%
#   group_by(Test, Group_1, Group_2) %>%
#   summarise(
#     Total_OBS = n() ,#length(unique(hs_symbol)),
#     Total_DEG = sum(abs(logFC)>1 & FDR < 0.05 & is_bio),
#     DEG_In_ZONULE = sum(IS_ZONULE & abs(logFC)>1 & FDR < 0.05 & is_bio),
#     UP_In_ZONULE = sum(IS_ZONULE & (logFC > 1 & FDR < 0.05 & is_bio)),
#     DOWN_In_ZONULE= sum(IS_ZONULE & (logFC < -1 & FDR < 0.05 & is_bio)),
#     NON_DEG_In_ZONULE=sum(IS_ZONULE & !(abs(logFC)>1 & FDR < 0.05 & is_bio)),
#     NO_OBS_In_ZONULE=length(setdiff(zonules$Mouse.gene.stable.ID, gene_id))
#   ) %>% write.csv(path)

pax6.deg_master %>%
  filter(Filtered == "ribo" & Partition == "Pair") %>%
  group_by(Test, Group_1, Group_2, Hom_Group) %>%
  filter(
    rank(FDR) == 1 & rank(desc(abs(logFC))) == 1
  ) %>%
  filter(row_number() == 1) %>%
  group_by(Test, Group_1, Group_2) %>%
  summarise(
    Total_OBS = n() ,#length(unique(hs_symbol)),
    Total_DEG = sum(abs(logFC)>1 & FDR < 0.05 & is_bio),
    DEG_In_ZONULE = sum(IS_ZONULE & abs(logFC)>1 & FDR < 0.05 & is_bio),
    UP_In_ZONULE = sum(IS_ZONULE & (logFC > 1 & FDR < 0.05 & is_bio)),
    DOWN_In_ZONULE= sum(IS_ZONULE & (logFC < -1 & FDR < 0.05 & is_bio)),
    NON_DEG_In_ZONULE=sum(IS_ZONULE & !(abs(logFC)>1 & FDR < 0.05 & is_bio)),
    NO_OBS_In_ZONULE=length(setdiff(zonules$Mouse.gene.stable.ID, gene_id))
  ) %>% write.csv(path)
result_files <- append(result_files, path)

############ Generate Enhanced Volcano Plot Highlighting Zonules #############

pax6_deg <- pax6.deg_master %>%
  filter(Filtered == "ribo" & Partition == "Pair" & Test == "ExactTest") %>%
  filter(Group_1 == "WTE" & Group_2 == "P6E") %>%
  filter(Avg1 >2 | Avg2 > 2) %>%
  inner_join(
    pax6.master$genes %>%
      select(gene_id, SYMBOL),
    by="gene_id"
  ) %>% mutate(
    IS_ZONULE = (gene_id %in% zonules$Mouse.gene.stable.ID ) & ( abs(logFC) > 1) & (FDR < 0.05)
  )%>% filter(IS_ZONULE)

fn <- "ZONULE_GENES_DE_IN_PAX6_LEC.jpg"
path <- paste(results, fn, sep="/")
result_files <- append(result_files, path)

p <- EnhancedVolcano(
  pax6_deg %>% filter(IS_ZONULE),
  x = "logFC", y="FDR", lab=pax6_deg$SYMBOL,
  selectLab = pax6_deg %>% filter(IS_ZONULE) %>% pull(SYMBOL),
  drawConnectors = T, max.overlaps = 1000,
  pCutoff = 5*10^-2, 
  title = "Zonular Proteins influenced by Pax6 haploinsufficiency",
  subtitle = ""
)
ggsave(path, plot=p, width=12, height=7)

######################  Push script and data to Synapse ######################

# Add this script to the code dir
syn_script <- synFindEntityId(
  "Pax6_DE_Zonule_Analysis.R",
  parent=syn_code_dir
)

if(is.null(syn_script)){
  syn_script <- File(
    path="scripts/Pax6_DE_Zonule_Analysis.R",
    parent=syn_code_dir
  )

  syn_script <- synStore(
    syn_script
  )
}

synSetProvenance(
  syn_script,
  activity = Activity(
    name = "Pax6_Zonule_DEG",
    description = "Genes encoding Zonular Proteins that are Pax6 Dependent",
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
