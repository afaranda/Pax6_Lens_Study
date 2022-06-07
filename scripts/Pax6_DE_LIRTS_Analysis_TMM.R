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
library(synapser)
library(ggrepel)
options(echo=T)

# Enter Working Directory and Load Raw Data
setwd('/Users/adam/Documents/23Feb2022_Pax6_Study_DEG_Analysis/')
source('scripts/Overlap_Comparison_Functions.R')
source('scripts/Wrap_edgeR_Functions.R')
source('scripts/Excel_Write_Functions.R')

wd<-getwd()
results<-paste(wd,'results/LIRTS_Analysis',sep='/')
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
syn_deg_sps_dir <- synFindEntityId("Pax6_LIRTS_Analysis", parent = syn_project)

if(is.null(syn_deg_sps_dir)){
  folder <- Folder(
    name="Pax6_LIRTS_Analysis",
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

########################## Load LIRTS Deg Tables #############################

### For now, use a set of previously generated DEG Tables
### Next week: Integrate the LIRTS Deg Table script with this project

syn_lirts_dgelists <- synFindEntityId(
  "LIRTS_master_dgelist.Rdata",
  parent=syn_data_dir
)

if(file.exists('data/LIRTS_master_dgelist.Rdata')){
  load('data/LIRTS_master_dgelist.Rdata')
} else if(!is.null(syn_lirts_dgelists)){
  synGet(
    syn_lirts_dgelists, downloadLocation=data_dir
  )
  load('data/LIRTS_master_dgelist.Rdata')
}else{
  stop("Run Prepare_Expression_DGELists.R first")
}

syn_lirts_deg_master <- synFindEntityId(
  "LIRTS_Master_DEG_Table.Rdata", 
  parent=syn_deg_tbl_dir
)

if(file.exists("results/")){
  load("results/LIRTS_Master_DEG_Table.Rdata")
} else if(!is.null(syn_lirts_deg_master)){
  synGet(
    syn_lirts_deg_master, downloadLocation="results"
  )
  load("results/LIRTS_Master_DEG_Table.Rdata")
}else{
  stop("Run Generate_LIRTS_DEG_Tables.R first")
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

########################### Load KEGG Pathway data ###########################
kegg_mmu <- ROntoTools::keggPathwayGraphs(
  organism = "mmu"
)

mmu04060_entrez <- as.numeric(
  gsub(
    "mmu:", "",
    kegg_mmu[['path:mmu04060']]@nodes
  )
)

mmu04060_genes <- AnnotationDbi::select(
  org.Mm.eg.db::org.Mm.eg.db,
  columns = c("SYMBOL", "GENENAME"),
  keys=as.character(mmu04060_entrez),
  keytype = "ENTREZID"
)

mmu04151_entrez <- as.numeric(
  gsub(
    "mmu:", "",
    kegg_mmu[['path:mmu04151']]@nodes
  )
)

mmu04151_genes <- AnnotationDbi::select(
  org.Mm.eg.db::org.Mm.eg.db,
  columns = c("SYMBOL", "GENENAME"),
  keys=as.character(mmu04151_entrez),
  keytype = "ENTREZID"
)
## Final list of existing data files used by this script
used_files <- list(
  syn_sun_targets,
  syn_trrust_targets,
  syn_dgelists,
  syn_deg_master,
  syn_lirts_dgelists,
  syn_lirts_deg_master
)


#####################  Setup group and contrast iterators ####################

### Compare the Pax6 vs WT LEC to the pairwise exact test for
### each post surgical interval and report:
###
## Analyses by study
###   1.) Contingency tables for Overall Differential Expression
###   2.) Contingency tables for directional dependence
###   3.) Fold Change Plot.

compare_deg <- function(
  pax6_deg = data.frame(),
  injury_deg = data.frame(),
  result_label = "DNA1_6_Hours_Post_Injury"
){

  ### Run the Fisher's Exact Test for Overrepresentation
  universe <- union(
    pax6_deg$gene_id,
    injury_deg$gene_id
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
      injury_deg %>%
        select(
          gene_id, injury_logFC = logFC, injury_FDR = FDR,
          injury_Avg1=Avg1, injury_Avg2=Avg2
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
      IS_INJURY = (
        abs(injury_logFC) > 1 & injury_FDR < 0.05 &
          (injury_Avg1 > 2 | injury_Avg2 > 2) &
          (abs(injury_Avg1 - injury_Avg2) > 2) &
          !is.na(injury_logFC)
      ),
      UP_INJURY = (
        injury_logFC > 1 & injury_FDR < 0.05 &
          (injury_Avg1 > 2 | injury_Avg2 > 2) &
          (abs(injury_Avg1 - injury_Avg2) > 2) &
          !is.na(injury_logFC)
      )
    ) %>%
    left_join(
      pax6.master$genes %>%
        select(
          gene_id, SYMBOL, DESCRIPTION
        ),
      by="gene_id"
    )
  
  ## Prepare Overall Enrichment Table
  print("Overall overrepresentation")
  contingency <- table(
    test_data$IS_PAX6,
    test_data$IS_INJURY
  )
  print(contingency)
  print(fisher.test(contingency))
  
  print("Directional Dependence")
  contingency <- table(
    test_data %>%
      filter(IS_PAX6 & IS_INJURY) %>%
      select(UP_PAX6, UP_INJURY)
  )
  print(sum(contingency))
  print(contingency)
  print(fisher.test(contingency))
  
  return(
    test_data %>%
      filter(IS_PAX6 & IS_INJURY) %>%
      mutate(
        PAX6_TARGET = SYMBOL %in% combined_targets$Target
      )
  )
}


# Define A list of named contrasts; each element points to a vector with
# a pair of group labels. Positive fold changes will be associated
# with the second group listed. 
contrasts=list(
  DBI_WT24vs0H=c('WT_0H_DBI', 'WT_24H_DBI', 'DBI_Wildtype'),
  DBI_WT48vs0H=c('WT_0H_DBI', 'WT_48H_DBI', 'DBI_Wildtype'),
  DNA1_WT6vs0H=c('WT_0H_DNA1', 'WT_6H_DNA1', 'DNA1_Wildtype'),
  DNA1_WT24vs0H=c('WT_0H_DNA1', 'WT_24H_DNA1', 'DNA1_Wildtype'),
  DNA2_WT120vs0H=c('WT_0H_DNA2', 'WT_120H_DNA2', 'DNA2_Wildtype'),
  DNA3_WT72vs0H=c('WT_0H_DNA3', 'WT_72H_DNA3', 'DNA2_Wildtype'),
  LFC_P6vsWT = c('WTF', 'P6F')
)


####### Iterate over contrasts and extract DE Results for Pax6 Targets #######
########### Calculate Fisher's Exact Test results for each contrast  #########
result_files <- c()
pax6_deg <- pax6.deg_master %>% filter(
  Partition == "Pair",
  Filtered == "ribo",
  Group_1 == "WTE",
  Group_2 == "P6E",
  Test == "ExactTest"
)

pax6_injury_deg_table <- data.frame()
fn <- "LIRTS_Enrichment_In_Pax6_DEG.txt"
path <- paste(results, fn, sep="/")
result_files <- append(result_files, path)
sink(path)
for(c in names(contrasts)){
  inj <- res[[2]] %>%
    filter(
      Test == "ExactTest",
      Group_1 == contrasts[[c]][1],
      Group_2 == contrasts[[c]][2],
      Samples == contrasts[[c]][3]
    )
  if(c == "LFC_P6vsWT"){
    inj <- pax6.deg_master %>% filter(
      Partition == "Pair",
      Filtered == "ribo",
      Group_1 == "WTF",
      Group_2 == "P6F",
      Test == "ExactTest"
    )
  }
  print(c)
  pax6_injury_deg_table <-bind_rows(
    pax6_injury_deg_table,
    compare_deg(pax6_deg, inj, result_label=c)
  )
}

sink()

fn <- "LIRTS_DEG_In_Pax6_LEC_Long_Table.csv"
path <- paste(results, fn, sep="/")
write.csv(pax6_injury_deg_table, path)
result_files <- append(result_files, path)
###################### Setup Intersection Spreadsheets #######################
for(c in names(contrasts)){
  inj <- res[[2]] %>%
    filter(
      Test == "ExactTest",
      Group_1 == contrasts[[c]][1],
      Group_2 == contrasts[[c]][2],
      Samples == contrasts[[c]][3]
    )
  
  if(c == "LFC_P6vsWT"){
    inj <- pax6.deg_master %>% filter(
      Partition == "Pair",
      Filtered == "ribo",
      Group_1 == "WTF",
      Group_2 == "P6F",
      Test == "ExactTest"
    )
  }
  
  print(c)
  
  C1 <- "P6vsWT"
  C2 <- paste0(
    gsub(
      "_DBI","",
      gsub(
        "_DNA1","",
        gsub(
          "_DNA2","",
          gsub(
            "_DNA3","",
            contrasts[[c]][2]
          )
        )
      )  
    ) ,"vs",
    gsub(
      "_DBI","",
      gsub(
        "_DNA1","",
        gsub(
          "_DNA2","",
          gsub(
            "_DNA3","",
            contrasts[[c]][1]
          )
        )
      )  
    )
  )
  template <- paste0("templates/Pax6_LEC_",c,"_overlap_template.xlsx")
  fn <- paste0("Pax6_HS_vs_WT_X_",c,"DEG_Comparison.xlsx")
  path <- paste(results, fn, sep="/")
  write.csv(pax6_injury_deg_table, path)
  result_files <- append(result_files, path)
  createMethodComparisonSpreadsheet(
    C1 = C1, C2 = C2, template =template, #"templates/overlap.xlsx",
    dg1 = pax6_deg, dg2 = inj, pref = "PCO", fname = path,
    dg2.me = 2, dg1.ds = "Pax6 HS vs WT LEC",
    dg2.ds = c, unlog = T,
    dg2.bioFun=bioSigRNASeq, idc = "gene_id",
    annot=pax6.master$genes %>% select(gene_id, SYMBOL, DESCRIPTION), #rnc=comp.meta[[c]][["rnc"]]
  )
}

######################       Build Pivoted Table        ######################

#### THIS TABLE IS WRONG
fn <- "LIRTS_DEG_In_Pax6_LEC_Injury_Resp_Summary.csv"
path <- paste(results, fn, sep="/")
write.csv(pax6_injury_deg_table, path)
result_files <- append(result_files, path)

data.frame(
  gene_id = pax6_injury_deg_table %>%
    filter(IS_PAX6) %>%
    pull(gene_id) %>% unique()
) %>% 
  left_join(
    pax6_injury_deg_table %>% 
      select(
        gene_id, SYMBOL, DESCRIPTION, PAX6_TARGET,
        matches("^pax6", ignore.case=F)) %>%
      distinct(), 
    by="gene_id"
  ) %>% 
  left_join(
    pax6_injury_deg_table %>% 
      mutate(
        Injury_Response = ifelse(
          UP_INJURY, "Upregulated",
          ifelse(
            IS_INJURY, "Downregulated",
            "Not DE"
          )
        ),
        Contrast = factor(
          Contrast, levels=c(
            "DNA1_WT6vs0H", "DNA1_WT24vs0H", "DBI_WT24vs0H",
            "DBI_WT48vs0H","DNA2_WT120vs0H", "DNA3_WT72vs0H",
            "LFC_P6vsWT"
          )
        )
      ) %>% 
      select(gene_id, Contrast, Injury_Response) %>% 
      pivot_wider(
        id_cols = "gene_id", 
        names_from="Contrast", 
        values_from = "Injury_Response", 
        values_fill = "Not Observed"
      ),
    by="gene_id"
  ) -> pivoted 
pivoted %>% write.csv(path)
################### Cytokine-Cytokine Receptor Comparison ####################

ccr_genes <-c(
  "Il1rn", "Tnfsf8", "Ccl5", "Lepr", "Il17re", "Tnfrsf11b", "Ccl7", "Ccl17", 
  "Cxcl14", "Eda2r", "Tgfbr2", "Il6ra", "Cd40", "Gdf15", "Bmpr1b", "Lif",
  "Tgfb1", "Clcf1", "Il17rb", "Ccl2", "Cx3cl1", "Tnfrsf19", "Bmp6",
  "Tnfrsf12a", "Ackr4", "Tnfrsf21", "Osmr", "Il3ra", "Csf1", "Relt", "Cxcl12",
  "Ghr", "Tnfrsf11a", "Lifr", "Tnfrsf25", "Ctf1", "Tgfb3", "Tnfrsf1b", 
  "Gdf10", "Il17rc", "Il10rb"
)

ccr_genes <- pax6.master$genes %>%
  filter(SYMBOL %in% ccr_genes) %>%
  select(gene_id, SYMBOL)


ccr_genes <- inner_join(
  ccr_genes, 
  pax6.deg_master %>%
    filter(
      Test == "ExactTest" & 
        Group_1 == "WTE" & 
        Group_2 == "P6E" & 
        Filtered == "ribo" &
        Partition == "Pair" &
        abs(logFC) > 1 &
        FDR  < 0.05
    ) %>%
    select(
      gene_id,
      Pax6_logFC = logFC
    ) 
)

ccr_genes <- inner_join(
  ccr_genes, 
  res[[2]] %>%
    filter(
      Test == "ExactTest" & 
        Group_1 == "WT_0H_DBI" & 
        Group_2 == "WT_24H_DBI" & 
        Samples == "DBI_Wildtype" &
        abs(logFC) > 1 &
        FDR  < 0.05
    ) %>%
    select(
      gene_id,
      Injury_logFC = logFC
    ) 
)

fn <- "Cytokine_Receptor_Pax6_X_LIRTS_Focused.jpg"
path <- paste(results, fn, sep="/")
result_files <- append(result_files, path)
ggsave(
  path,
  ggplot(
    ccr_genes, aes(x=Injury_logFC, y=Pax6_logFC)
  ) + 
    geom_point() + 
    geom_label_repel(aes(label=SYMBOL)) +
    xlab("Log Fold Change 24 vs 0 Hours After Injury") +
    ylab("Log Fold Change Sey vs Wildtype Lens Epithelium") +
    theme(
      axis.text = element_text(size=12),
      axis.title = element_text(size=12)
    )
)

################### Fibrotic Mediators / Cytokines table #####################
inflammation_fibrosis_genes <-c(
  "Ccl5",  "Ccl7", "Ccl17","Tgfbr2", "Gdf15", "Tgfb1", "Ccl2",
  "Csf1", "Tnfrsf25", "Tgfb3", "Itgb8", "Itga1", "Col1a1",
  "Fn1", "Itga5", "Thbs1", "Thbs2"
)
fn <- "Inflammation_and_Fibrosis_Gene_Table_Pax6_X_LIRTS.csv"
path <- paste(results, fn, sep="/")
write.csv(pax6_injury_deg_table, path)
result_files <- append(result_files, path)
#pivoted %>% filter(SYMBOL %in% inflammation_fibrosis_genes) %>%
#write.csv(path)

########################## Plot All Genes comparison #########################

## for each deg table:
## get log fold changes and FDR values of all genes expressed
## at biologically significant levels in at least one condition
## and inner join DEG Tables on gene_id. Add gene symbols, and indictor column
## to flag genes that should get a call out.
pathways <- list(PI3K = mmu04151_genes, Cytokine = mmu04060_genes)
for(pw in names(pathways)){
  for(c in names(contrasts)){
    inj <- res[[2]] %>%
      filter(
        Test == "ExactTest",
        Group_1 == contrasts[[c]][1],
        Group_2 == contrasts[[c]][2],
        Samples == contrasts[[c]][3]
      )
    if(c == "LFC_P6vsWT"){
      inj <- pax6.deg_master %>% filter(
        Partition == "Pair",
        Filtered == "ribo",
        Group_1 == "WTF",
        Group_2 == "P6F",
        Test == "ExactTest"
      )
    }
    
    print(c)
    df <- inner_join(
      pax6_deg %>%
        filter(Avg1 > 2 | Avg2 > 2) %>%
        filter(FDR < 0.05)  %>%
        select(gene_id, Pax6_logFC = logFC, pax6_fdr=FDR) ,
      inj %>%
        filter(Avg1 > 2 | Avg2 > 2) %>%
        filter(FDR < 0.05) %>%
        select(gene_id, Injury_logFC = logFC, inj_fdr=FDR),
      by="gene_id"
    ) %>% inner_join(
      pax6.master$genes %>%
        select(gene_id, SYMBOL),
      by="gene_id"
    ) %>% mutate(
      SYMBOL = ifelse(SYMBOL %in% pathways[[pw]]$SYMBOL, SYMBOL, "")
    ) %>% rowwise() %>%
      mutate(
        DIST = sqrt(sum(c(Injury_logFC^2), Pax6_logFC^2))
      ) %>% group_by() %>%
      mutate(
        DIST = (DIST / max(DIST))
      ) %>% rowwise() %>%
      mutate(
        ALPHA = ifelse(SYMBOL == "", max(c(DIST, 0.05)), 1)
      ) %>%
      group_by() %>%
      mutate(
        SIZE = ifelse(SYMBOL == "",1, 2),
        COLOR = ifelse(
          (Pax6_logFC > 0 & Injury_logFC > 0), "darkred",ifelse(
            (Pax6_logFC < 0 & Injury_logFC < 0), "blue",
            "black"
          )
        )
      ) %>%
      mutate(
        COLOR = ifelse(SYMBOL == "", "black", COLOR)
      )
    
    fn <- paste0(pw,"_Pax6_X_",c,"_All_Genes.jpg")
    path <- paste(results, fn, sep="/")
    result_files <- append(result_files, path)
    
    p <- ggplot(df,aes(x=Injury_logFC, y=Pax6_logFC)) +
      geom_point(alpha = df$ALPHA, size = df$SIZE, colour = df$COLOR) +
      geom_label_repel(aes(label=SYMBOL), max.overlaps = 1000) +
      xlab("Log Fold Change 24 vs 0 Hours After Injury") +
      ylab("Log Fold Change Sey vs Wildtype Lens Epithelium") +
      theme(
        axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.position = "none"
      )
    print(nrow(df))
    ggsave( path, p, height=8, width = 8)
  }
}
  ## Use ggplot2 to to construct the following



######################  Push script and data to Synapse ######################

# Add this script to the code dir
syn_script <- synFindEntityId(
  "Pax6_DE_LIRTS_Analysis_TMM.R",
  parent=syn_code_dir
)

if(is.null(syn_script)){
  syn_script <- File(
    path="scripts/Pax6_DE_LIRTS_Analysis_TMM.R",
    parent=syn_code_dir
  )

  syn_script <- synStore(
    syn_script
  )
}

synSetProvenance(
  syn_script,
  activity = Activity(
    name = "Pax6_LIRTS_DEG",
    description = "Genes that are both Pax6 and Injury Dependent",
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
