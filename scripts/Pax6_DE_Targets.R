##############################################################################
#                                                                            #
#  File: Pax6_DE_Targets.R                                                   #
#  Author: Adam Faranda                                                      #
#  Created: June 18, 2020                                                    #
#  Purpose: Identify validated Pax6 targets that are DE in this study        #
#                                                                            #
##############################################################################

################## Load Libraries and Source Dependencies ####################
library(edgeR)
library(tidyr)
library(tibble)
library(dplyr)
library(tibble)
#library(synapser)
options(echo=T)

# Enter Working Directory and Load Raw Data
setwd('/Users/adam/Documents/23Feb2022_Pax6_Study_DEG_Analysis/')
source('scripts/synapse_reticulate_wrapper.R')
source('scripts/Overlap_Comparison_Functions.R')
source('scripts/Wrap_edgeR_Functions.R')
wd<-getwd()
results<-paste(wd,'results/Pax6_Targets',sep='/')
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
syn_deg_sps_dir <- synFindEntityId("Pax6_Targets", parent = syn_project)

if(is.null(syn_deg_sps_dir)){
  folder <- Folder(
    name="Pax6_Targets",
    parent=syn_project
  )
  syn_deg_sps_dir <- synStore(folder)
} 

## Check for the DEG_Tables folder
syn_deg_tbl_dir <- synFindEntityId(
  "DEG_Tables",
  syn_project
)

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

## Define design matrix for interaction model. 
design <- model.matrix(~ cell_type * genotype, pax6.master$samples)

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

### Genes considered PAX6 targets in Human by MSigDB
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

## Peak Annotations for lens from Sun et al. 2015
syn_sun_lens_peaks <- synFindEntityId(
  "Lens_Peak_Annotations.tsv", 
  parent = syn_data_dir
)

if(file.exists("data/Lens_Peak_Annotations.tsv")){
  sun_lens_peaks <- read.table(
    "data/Lens_Peak_Annotations.tsv",
    header=T, sep="\t"
  )
} else if(!is.null(syn_msig_targets)){
  synGet(
    syn_sun_lens_peaks, downloadLocation="data"
  )
  
  sun_lens_peaks <- read.table(
    "data/Lens_Peak_Annotations.tsv",
    header=T, sep="\t"
  )
  
}else{
  stop("Get Sun Lens Peaks First")
}

## Forebrain Peak Annotations for lens from Sun et al. 2015
syn_sun_forebrain_peaks <- synFindEntityId(
  "Forebrain_Peak_Annotations.tsv", 
  parent = syn_data_dir
)

if(file.exists("data/Forebrain_Peak_Annotations.tsv")){
  sun_forebrain_peaks <- read.table(
    "data/Forebrain_Peak_Annotations.tsv",
    header=T, sep="\t"
  )
} else if(!is.null(syn_msig_targets)){
  synGet(
    syn_sun_forebrain_peaks, downloadLocation="data"
  )
  
  sun_forebrain_peaks <- read.table(
    "data/Forebrain_Peak_Annotations.tsv",
    header=T, sep="\t"
  )
  
}else{
  stop("Get Sun Forebrain Peaks First")
}

### Update Col4a3bp to Cert1 in sSun Peaks
sun_lens_peaks$mm9_Symbol <- gsub(
  "Col4a3bp", "Cert1", sun_lens_peaks$mm9_Symbol
)

sun_lens_peaks$mm10_Symbol <- gsub(
  "Col4a3bp", "Cert1", sun_lens_peaks$mm10_Symbol
)

sun_forebrain_peaks$mm9_Symbol <- gsub(
  "Col4a3bp", "Cert1", sun_forebrain_peaks$mm9_Symbol
)

sun_forebrain_peaks$mm10_Symbol <- gsub(
  "Col4a3bp", "Cert1", sun_forebrain_peaks$mm10_Symbol
)






### Many of the symbols missing in the mm10 column were lost during
### Mapping.  In most cases this is because the RefSeq ID has been
### deprecated, either due to new data about the locus or because
### the specific transcript could not be validated (eg Auts2).  
### In this analysis, the transcript is irrelevant thus we will
### replace missing mm10 symbols with their mm9 counterparts.
sun_lens_peaks <- sun_lens_peaks %>%
  rowwise() %>%
  mutate(
    mm10_Symbol = ifelse(mm10_Symbol == "", mm9_Symbol, mm10_Symbol)
  ) %>% group_by() %>%
  filter(!grepl("\\|", mm10_Symbol))

sun_forebrain_peaks <- sun_forebrain_peaks %>%
  rowwise() %>%
  mutate(
    mm10_Symbol = ifelse(mm10_Symbol == "", mm9_Symbol, mm10_Symbol)
  ) %>% group_by() %>%
  filter(!grepl("\\|", mm10_Symbol))





### Collapse to genes and count peaks of each category
sun_lens_peaks %>%
  group_by(mm10_Symbol, Peak) %>%
  filter(row_number() == 1) %>%
  group_by(mm10_Symbol) %>%
  summarise(
    Adjacent_Peaks = n(),
    Promoter_Peaks = sum(grepl("prot", Type)),
    Exon_Peaks = sum(grepl("exon", Type)),
    Intron_Peaks = sum(grepl("intron", Type)),
    Distal_Peaks = sum(grepl("distal", Type))
  ) %>% group_by() %>%
  rename(Target = mm10_Symbol) %>% 
  as.data.frame() -> sun_lens_peaks

sun_forebrain_peaks %>%
  group_by(mm10_Symbol, Peak) %>%
  filter(row_number() == 1) %>%
  group_by(mm10_Symbol) %>%
  summarise(
    Adjacent_Peaks = n(),
    Promoter_Peaks = sum(grepl("prot", Type)),
    Exon_Peaks = sum(grepl("exon", Type)),
    Intron_Peaks = sum(grepl("intron", Type)),
    Distal_Peaks = sum(grepl("distal", Type))
  ) %>% group_by() %>%
  rename(Target = mm10_Symbol) %>% 
  as.data.frame() -> sun_forebrain_peaks

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
  syn_deg_master,
  syn_sun_lens_peaks,
  syn_sun_forebrain_peaks
)


#####################  Setup group and contrast iterators ####################


# Define A list of named contrasts; each element points to a vector with
# a pair of group labels. Positive fold changes will be associated
# with the second group listed. 
contrasts=list(
  WT_Fibers_vs_Epi=c('WTE', 'WTF'),
  P6_Fibers_vs_Epi=c('P6E', 'P6F'),
  Epi_P6_vs_WT=c('WTE', 'P6E'),
  Fibers_P6_vs_WT=c('WTF', 'P6F')#,
  # `2Way_Fib_vs_Epi`=c("Epi", "Fib"),
  # `2Way_P6_vs_WT`=c("WT", "P6"),
  # `2Way_Interaction`=c("REF", "INX")#,
  # DegNorm_WT_Fibers_vs_Epi=c('WTE', 'WTF'),
  # DegNorm_P6_Fibers_vs_Epi=c('P6E', 'P6F'),
  # DegNorm_Epi_P6_vs_WT=c('WTE', 'P6E'),
  # DegNorm_Fibers_P6_vs_WT=c('WTF', 'P6F'),
  # DegNorm_2Way_Fib_vs_Epi=c("Epi", "Fib"),
  # DegNorm_2Way_P6_vs_WT=c("WT", "P6"),
  # DegNorm_2Way_Interaction=c("REF", "INX")
)

# Define A list of named contrasts; each element points to a string which
# briefly describes the contrast
contrast_descriptions<-list(
  WT_Fibers_vs_Epi="Wildtype lens fiber cells vs lens epithelial cells",
  P6_Fibers_vs_Epi="Sey lens fiber cells vs lens epithelial cells",
  Epi_P6_vs_WT="Sey lens epithelial cells vs lens epithelial cells",
  Fibers_P6_vs_WT="Sey lens fiber cells vs lens fiber cells"#,
  # `2Way_Fib_vs_Epi`="2 Way interaction model lens fiber cells vs lens fiber cells",
  # `2Way_P6_vs_WT`="2 Way interaction model Sey vs Wildtype",
  # `2Way_Interaction`="2 Way interaction model haploinsufficiency and differentiation",
  # DegNorm_WT_Fibers_vs_Epi="Wildtype lens fiber cells vs lens epithelial cells",
  # DegNorm_P6_Fibers_vs_Epi="Sey lens fiber cells vs lens epithelial cells",
  # DegNorm_Epi_P6_vs_WT="Sey lens epithelial cells vs lens epithelial cells",
  # DegNorm_Fibers_P6_vs_WT="Sey lens fiber cells vs lens fiber cells",
  # DegNorm_2Way_Fib_vs_Epi="2 Way interaction model lens fiber cells vs lens fiber cells",
  # DegNorm_2Way_P6_vs_WT="2 Way interaction model Sey vs Wildtype",
  # DegNorm_2Way_Interaction="2 Way interaction model haploinsufficiency and differentiation"
)

####### Iterate over contrasts and extract DE Results for Pax6 Targets #######
########### Calculate Fisher's Exact Test results for each contrast  #########
result_files <- c()
for(c in names(contrasts)){
  if(grepl("DegNorm",c)){
    master <- pax6.master.dgn
    deg <- pax6.deg_master %>% 
      filter(Filtered == "degnorm")
    filt <- "degnorm"
  } else {
    master <- pax6.master
    deg <- pax6.deg_master %>% 
      filter(Filtered == "ribo")
    filt <- "ribo"
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
    
  } else {
    gr <- colnames(master)
    
    obj <- processByDesign(
      y=master, design=design,
      norm=ifelse(filt != "degnorm", "TMM", "NONE")
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
        obj$dge$genes %>%
          select(
            gene_id, SYMBOL,
            matches("FPKM")
          ),
        by="gene_id"
      )
    
    print(nrow(degSet))
  }   
  
  ## Deduplicate symbols based on feature with maximum expression
  degSet <- union(
    degSet %>% 
      filter(SYMBOL !="") %>%
      group_by(SYMBOL) %>%
      filter((Avg1 + Avg2) == max(Avg1 + Avg2)),
    degSet %>% filter(SYMBOL =="")
  ) %>% group_by() %>% as.data.frame()
  print(paste("Deduped ", nrow(degSet)))
  
  ## Add Bio Sig Flag
  degSet <- degSet %>%
    mutate(
      IS_DEG = (
        abs(logFC) > 1 & FDR < 0.05
      ),
      IS_BIO = (
        abs(logFC) > 1 & FDR < 0.05 & 
          (Avg1 > 2 | Avg2 > 2) & 
          (abs(Avg1 - Avg2)>2)
      )
    )
  
  # Name headers for columns storing group average FPKM estimates.
  names(degSet) <- gsub("Avg1", paste0(g1, "_Avg_FPKM"), names(degSet))
  names(degSet) <- gsub("Avg2", paste0(g2, "_Avg_FPKM"),names(degSet))
  
  print(c)
  
  #### Tabulate Pax6 Peaks reported For Lens in Sun et. al 2015
  deg_sun_lens_peaks <- inner_join(
    sun_lens_peaks,
    degSet,
    by=c("Target"="SYMBOL")
  ) %>%
    filter(abs(logFC) > 1 & FDR < 0.05)
  
  ## Get file name for the table
  fn<-paste0("P6T_",c,"_Sun_2015_Lens_Peaks.csv")
  path<-paste(results, fn, sep="/")
  write.csv(deg_sun_lens_peaks, file=path, row.names = F)
  result_files <- append(result_files, path)
  
  ## Run Enrichment tests on Sun_2015 Targets
  fn<-paste0("P6T_",c,"_Sun_2015_Lens_Peak_Enrichment.txt")
  path<-paste(results, fn, sep="/")
  result_files <- append(result_files, path)
  
  sink(path)
  print("Statistically Significant")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% sun_lens_peaks$Target)
    ) %>%
    select( IS_DEG, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
   contingency, alternative = "two.sided"
  ) %>% print()
  
  print("Statistically Positive")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% sun_lens_peaks$Target),
      IS_DEG = IS_DEG & logFC > 1
    ) %>%
    
    select( IS_DEG, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
    contingency,alternative = "two.sided"
  ) %>% print()
  
  print("Statistically Negative")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% sun_lens_peaks$Target),
      IS_DEG = IS_DEG & logFC < -1
    ) %>%
    select(IS_DEG, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
    contingency,alternative = "two.sided"
  ) %>% print()
  
  
  print("Biologically Significant")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% sun_lens_peaks$Target)
    ) %>%
    select( IS_BIO, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
    contingency, alternative = "two.sided"
  ) %>% print()
  
  print("Biologically Positive")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% sun_lens_peaks$Target),
      IS_BIO = IS_BIO & logFC > 1
    ) %>%
    select( IS_BIO, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
    contingency,alternative = "two.sided"
  ) %>% print()
  
  print("Biologically Negative")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% sun_lens_peaks$Target),
      IS_BIO = IS_BIO & logFC < -1
    ) %>%
    select( IS_BIO, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
    contingency,alternative = "two.sided"
  ) %>% print()
  
  sink()
  
  #### Tabulate Pax6 Peaks reported For Forebrain in Sun et. al 2015
  deg_sun_forebrain_peaks <- inner_join(
    sun_forebrain_peaks,
    degSet,
    by=c("Target"="SYMBOL")
  ) %>%
    filter(abs(logFC) > 1 & FDR < 0.05)
  
  ## Get file name for the table
  fn<-paste0("P6T_",c,"_Sun_2015_Forebrain_Peaks.csv")
  path<-paste(results, fn, sep="/")
  write.csv(deg_sun_forebrain_peaks, file=path, row.names = F)
  result_files <- append(result_files, path)
  
  ## Run Enrichment tests on Sun_2015 Targets
  fn<-paste0("P6T_",c,"_Sun_2015_Forebrain_Peak_Enrichment.txt")
  path<-paste(results, fn, sep="/")
  result_files <- append(result_files, path)
  
  sink(path)
  print("Statistically Significant")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% sun_forebrain_peaks$Target)
    ) %>%
    select( IS_DEG, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
    contingency, alternative = "two.sided"
  ) %>% print()
  
  print("Statistically Positive")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% sun_forebrain_peaks$Target),
      IS_DEG = IS_DEG & logFC > 1
    ) %>%
    
    select( IS_DEG, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
    contingency,alternative = "two.sided"
  ) %>% print()
  
  print("Statistically Negative")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% sun_forebrain_peaks$Target),
      IS_DEG = IS_DEG & logFC < -1
    ) %>%
    select(IS_DEG, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
    contingency,alternative = "two.sided"
  ) %>% print()
  
  
  print("Biologically Significant")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% sun_forebrain_peaks$Target)
    ) %>%
    select( IS_BIO, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
    contingency, alternative = "two.sided"
  ) %>% print()
  
  print("Biologically Positive")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% sun_forebrain_peaks$Target),
      IS_BIO = IS_BIO & logFC > 1
    ) %>%
    select( IS_BIO, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
    contingency,alternative = "two.sided"
  ) %>% print()
  
  print("Biologically Negative")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% sun_forebrain_peaks$Target),
      IS_BIO = IS_BIO & logFC < -1
    ) %>%
    select( IS_BIO, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
    contingency,alternative = "two.sided"
  ) %>% print()
  
  sink()
  
  #### Tabulate Targets reported in Sun et. al 2015
  deg_sun_targets <- inner_join(
    sun_targets,
    degSet,
    by=c("Target"="SYMBOL")
  ) %>%
    filter(abs(logFC) > 1 & FDR < 0.05)
  
  ## Get file name for the table
  fn<-paste0("P6T_",c,"_Sun_2015_Targets.csv")
  path<-paste(results, fn, sep="/")
  write.csv(deg_sun_targets, file=path, row.names = F)
  result_files <- append(result_files, path)
  
  ## Run Enrichment tests on Sun_2015 Targets
  fn<-paste0("P6T_",c,"_Sun_2015_Target_Enrichment.txt")
  path<-paste(results, fn, sep="/")
  result_files <- append(result_files, path)
  
  sink(path)
  print("Statistically Significant")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% sun_targets$Target)
    ) %>%
    select( IS_DEG, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
    contingency,alternative = "two.sided"
  ) %>% print()
  
  print("Statistically Positive")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% sun_targets$Target),
      IS_DEG = IS_DEG & logFC > 1
    ) %>%
    
    select( IS_DEG, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
    contingency,alternative = "two.sided"
  ) %>% print()
  
  print("Statistically Negative")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% sun_targets$Target),
      IS_DEG = IS_DEG & logFC < -1
    ) %>%
    select(IS_DEG, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
    contingency,alternative = "two.sided"
  ) %>% print()
  
  
  print("Biologically Significant")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% sun_targets$Target)
    ) %>%
    select( IS_BIO, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
    contingency,alternative = "two.sided"
  ) %>% print()
  
  print("Biologically Positive")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% sun_targets$Target),
      IS_BIO = IS_BIO & logFC > 1
    ) %>%
    select( IS_BIO, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
    contingency,alternative = "two.sided"
  ) %>% print()
  
  print("Biologically Negative")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% sun_targets$Target),
      IS_BIO = IS_BIO & logFC < -1
    ) %>%
    select( IS_BIO, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
    contingency,alternative = "two.sided"
  ) %>% print()
  
  sink()
  
  #### Tabulate Targets Recognized by TRRUST
  deg_trrust_targets <- inner_join(
    trrust_targets,
    degSet,
    by=c("Target"="SYMBOL")
  ) %>%
    filter(abs(logFC) > 1 & FDR < 0.05)
  
  ## Get file name for the table
  fn<-paste0("P6T_",c,"_TRRUST_Targets.csv")
  path<-paste(results, fn, sep="/")
  write.csv(deg_trrust_targets, file=path, row.names = F)
  result_files <- append(result_files, path)
  
  ## Run Enrichment tests on TRRUST Targets
  fn<-paste0("P6T_",c,"_TRRUST_Target_Enrichment.txt")
  path<-paste(results, fn, sep="/")
  result_files <- append(result_files, path)
  
  sink(path)
  print("Statistically Significant")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% trrust_targets$Target)
    ) %>%
    select( IS_DEG, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
    contingency,alternative = "two.sided"
  ) %>% print()
  
  print("Statistically Positive")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% trrust_targets$Target),
      IS_DEG = IS_DEG & logFC > 1
    ) %>%
    select( IS_DEG, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
    contingency,alternative = "two.sided"
  ) %>% print()
  
  print("Statistically Negative")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% trrust_targets$Target),
      IS_DEG = IS_DEG & logFC < -1
    ) %>%
    select( IS_DEG, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
    contingency,alternative = "two.sided"
  ) %>% print()
  
  print("Biologically Significant")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% trrust_targets$Target)
    ) %>%
    select( IS_BIO, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
    contingency,alternative = "two.sided"
  ) %>% print()
  
  print("Biologically Positive")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% trrust_targets$Target),
      IS_BIO = IS_BIO & logFC > 1
    ) %>%
    select( IS_BIO, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
    contingency,alternative = "two.sided"
  ) %>% print()
  
  print("Biologically Negative")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% trrust_targets$Target),
      IS_BIO = IS_BIO & logFC < -1
    ) %>%
    select( IS_BIO, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
    contingency,alternative = "two.sided"
  ) %>% print()
  sink()
  
  
  #### Tabulate Targets Recognized by TRRUST and Sun 2015
  deg_combined_targets <- inner_join(
    combined_targets,
    degSet,
    by=c("Target"="SYMBOL")
  ) %>%
    filter(abs(logFC) > 1 & FDR < 0.05)
  
  ## Get file name for the table
  fn<-paste0("P6T_",c,"_Combined_Targets.csv")
  path<-paste(results, fn, sep="/")
  write.csv(deg_combined_targets, file=path, row.names = F)
  result_files <- append(result_files, path)
  
  ## Run Enrichment tests on Sun_2015 Targets
  fn<-paste0("P6T_",c,"_Combined_Target_Enrichment.txt")
  path<-paste(results, fn, sep="/")
  result_files <- append(result_files, path)
  
  sink(path)
  print("Statistically Significant")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% combined_targets$Target)
    ) %>%
    select( IS_DEG, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
    contingency,alternative = "two.sided"
  ) %>% print()
  
  print("Statistically Positive")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% combined_targets$Target),
      IS_DEG = IS_DEG & logFC > 1
    ) %>%
    select( IS_DEG, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
    contingency,alternative = "two.sided"
  ) %>% print()
  
  print("Statistically Negative")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% combined_targets$Target),
      IS_DEG = IS_DEG & logFC < -1
    ) %>%
    select( IS_DEG, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
    contingency,alternative = "two.sided"
  ) %>% print()
  
  print("Biologically Significant")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% combined_targets$Target)
    ) %>%
    select( IS_BIO, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
    contingency,alternative = "two.sided"
  ) %>% print()
  
  print("Biologically Positive")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% combined_targets$Target),
      IS_BIO = IS_BIO & logFC > 1
    ) %>%
    select( IS_BIO, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
    contingency,alternative = "two.sided"
  ) %>% print()
  
  print("Biologically Negative")
  print(c)
  contingency <- degSet %>%
    mutate(
      IS_PAX6 = (SYMBOL %in% combined_targets$Target),
      IS_BIO = IS_BIO & logFC < -1
    ) %>%
    select( IS_BIO, IS_PAX6) %>% table()
  
  print(contingency)
  fisher.test(
    contingency,alternative = "two.sided"
  ) %>% print()
  
  sink()
}

######################       Build Pivoted Table        ######################
combined_targets %>%
  left_join(
    pax6.deg_master %>%
      filter(
        Group_1 == "WTE", Group_2 == "WTF",
        Test == "ExactTest", Filtered == "ribo",
        Partition == "Pair"
      ) %>%
      inner_join(
        pax6.master$genes %>%
          select(gene_id, SYMBOL, DESCRIPTION),
        by="gene_id"
      ) %>%
      select(
        WTE_FPKM = Avg1, WTF_FPKM = Avg2,
        WTF_vs_WTE_LogFC = logFC,
        WTF_vs_WTE_FDR = FDR, SYMBOL
      ),
    by=c(Target="SYMBOL")
  ) %>%
  left_join(
    pax6.deg_master %>%
      filter(
        Group_1 == "P6E", Group_2 == "P6F",
        Test == "ExactTest", Filtered == "ribo",
        Partition == "Pair"
      ) %>%
      inner_join(
        pax6.master$genes %>%
          select(gene_id, SYMBOL),
        by="gene_id"
      ) %>%
      select(
        P6E_FPKM = Avg1, P6F_FPKM = Avg2,
        P6F_vs_P6E_LogFC = logFC,
        P6F_vs_P6E_FDR = FDR, SYMBOL
      ),
    by=c(Target="SYMBOL")
  )%>%
  rowwise() %>%
  filter(
    sum(across(matches("FPKM"), ~ is.na(.x))) < 4
  ) -> deg_combined_by_genotype

## Write Combination By Genotype Table
fn<-paste0(
  "P6T_Combined_Targets_Differentiation_By_Genotype.csv"
)
path<-paste(results, fn, sep="/")
write.csv(
  deg_combined_by_genotype, file=path, row.names = F
)
result_files <- append(result_files, path)

fn<-paste0(
  "P6T_Combined_Targets_Upregulation_Defect.csv"
)
path<-paste(results, fn, sep="/")
write.csv(
  deg_combined_by_genotype %>%
    filter(
      WTF_vs_WTE_LogFC  > 1 &
        WTF_vs_WTE_FDR < 0.05 &
        P6F_vs_P6E_LogFC <= 1
    ), file=path, row.names = F)
result_files <- append(result_files, path)

fn<-paste0(
  "P6T_Combined_Targets_Downregulation_Defect.csv"
)
path<-paste(results, fn, sep="/")
write.csv(
  deg_combined_by_genotype %>%
    filter(
      WTF_vs_WTE_LogFC  < -1 &
        WTF_vs_WTE_FDR < 0.05 &
        P6F_vs_P6E_LogFC >= -1
    ), file=path, row.names = F)
result_files <- append(result_files, path)

################### Generate Summary tables for Target DE ####################
fn <- "Pax6_Targets_DEG_Counts.csv"
path <- paste(results, fn, sep="/")

bind_rows(
  pax6.deg_master %>% 
    filter(
      Test == "ExactTest" & Partition == "Pair" & Filtered == "ribo"
    ) %>% inner_join(pax6.master$genes, by="gene_id") %>%
    rowwise() %>%
    mutate(IS_BIO =(Avg1>2 | Avg2>2) & abs(Avg1 - Avg2) > 2) %>%
    group_by(Group_1, Group_2) %>%
    summarize(
      Total_OBS = n() ,#length(unique(hs_symbol)),
      Total_DEG = sum(abs(logFC)>1 & FDR < 0.05 & IS_BIO),
      TARGET_OBS = sum(IS_TRRUST_PAX6_TARGET),
      DEG_TARGET = sum(abs(logFC)>1 & FDR < 0.05 & IS_BIO & IS_TRRUST_PAX6_TARGET),
      UP_TARGET = sum(logFC>1 & FDR < 0.05 & IS_BIO & IS_TRRUST_PAX6_TARGET),
      DOWN_TARGET = sum(logFC < -1 & FDR < 0.05 & IS_BIO & IS_TRRUST_PAX6_TARGET),
      NON_DEG_TARGET = sum(!(abs(logFC)>1 & FDR < 0.05 & IS_BIO) & IS_TRRUST_PAX6_TARGET),
      NO_OBS_TARGET = length(setdiff(trrust_targets$Target,SYMBOL))
    ) %>% group_by() %>% mutate(
      TARGET_SET = "TRRUST"
    ) %>% as.data.frame(),
  
  pax6.deg_master %>% 
    filter(
      Test == "ExactTest" & Partition == "Pair" & Filtered == "ribo"
    ) %>% inner_join(pax6.master$genes, by="gene_id") %>%
    rowwise() %>%
    mutate(IS_BIO =(Avg1>2 | Avg2>2) & abs(Avg1 - Avg2) > 2) %>%
    group_by(Group_1, Group_2) %>%
    summarize(
      Total_OBS = n() ,#length(unique(hs_symbol)),
      Total_DEG = sum(abs(logFC)>1 & FDR < 0.05 & IS_BIO),
      TARGET_OBS = sum(IS_SUN_PAX6_TARGET),
      DEG_TARGET = sum(abs(logFC)>1 & FDR < 0.05 & IS_BIO & IS_SUN_PAX6_TARGET),
      UP_TARGET = sum(logFC>1 & FDR < 0.05 & IS_BIO & IS_SUN_PAX6_TARGET),
      DOWN_TARGET = sum(logFC < -1 & FDR < 0.05 & IS_BIO & IS_SUN_PAX6_TARGET),
      NON_DEG_TARGET = sum(!(abs(logFC)>1 & FDR < 0.05 & IS_BIO) & IS_SUN_PAX6_TARGET),
      NO_OBS_TARGET = length(setdiff(sun_targets$Target,SYMBOL))
    ) %>% group_by()  %>% mutate(
      TARGET_SET = "SUN_TARGETS"
    ) %>% as.data.frame(),
  
  pax6.deg_master %>% 
    filter(
      Test == "ExactTest" & Partition == "Pair" & Filtered == "ribo"
    ) %>% inner_join(pax6.master$genes, by="gene_id") %>%
    rowwise() %>%
    mutate(IS_BIO =(Avg1>2 | Avg2>2) & abs(Avg1 - Avg2) > 2) %>%
    group_by(Group_1, Group_2) %>%
    summarize(
      Total_OBS = n() ,#length(unique(hs_symbol)),
      Total_DEG = sum(abs(logFC)>1 & FDR < 0.05 & IS_BIO),
      TARGET_OBS = sum(IS_SUN_PAX6_LENS_PEAK),
      DEG_TARGET = sum(abs(logFC)>1 & FDR < 0.05 & IS_BIO & IS_SUN_PAX6_LENS_PEAK),
      UP_TARGET = sum(logFC>1 & FDR < 0.05 & IS_BIO & IS_SUN_PAX6_LENS_PEAK),
      DOWN_TARGET = sum(logFC < -1 & FDR < 0.05 & IS_BIO & IS_SUN_PAX6_LENS_PEAK),
      NON_DEG_TARGET = sum(!(abs(logFC)>1 & FDR < 0.05 & IS_BIO) & IS_SUN_PAX6_LENS_PEAK),
      NO_OBS_TARGET = length(setdiff(sun_lens_peaks$Target,SYMBOL))
    ) %>% group_by() %>% mutate(
      TARGET_SET = "SUN_LENS_PEAKS"
    ) %>% as.data.frame(),
  pax6.deg_master %>% 
    filter(
      Test == "ExactTest" & Partition == "Pair" & Filtered == "ribo"
    ) %>% inner_join(pax6.master$genes, by="gene_id") %>%
    #group_by(Group_1, Group_2, SYMBOL) %>%
    #filter(Avg1+Avg2==max(Avg1+Avg2))%>%
    rowwise() %>%
    mutate(IS_BIO =(Avg1>2 | Avg2>2) & abs(Avg1 - Avg2) > 2) %>%
    group_by(Group_1, Group_2) %>%
    summarize(
      Total_OBS = n() ,#length(unique(hs_symbol)),
      Total_DEG = sum(abs(logFC)>1 & FDR < 0.05 & IS_BIO),
      TARGET_OBS = sum(IS_SUN_PAX6_FOREBRAIN_PEAK),
      DEG_TARGET = sum(abs(logFC)>1 & FDR < 0.05 & IS_BIO & IS_SUN_PAX6_FOREBRAIN_PEAK),
      UP_TARGET = sum(logFC>1 & FDR < 0.05 & IS_BIO & IS_SUN_PAX6_FOREBRAIN_PEAK),
      DOWN_TARGET = sum(logFC < -1 & FDR < 0.05 & IS_BIO & IS_SUN_PAX6_FOREBRAIN_PEAK),
      NON_DEG_TARGET = sum(!(abs(logFC)>1 & FDR < 0.05 & IS_BIO) & IS_SUN_PAX6_FOREBRAIN_PEAK),
      NO_OBS_TARGET = length(setdiff(sun_forebrain_peaks$Target,SYMBOL))
    ) %>% group_by() %>% mutate(
      TARGET_SET = "SUN_FOREBRAIN_PEAKS"
    ) %>% as.data.frame()
)  %>% write.csv(path)
result_files <- append(result_files, path)
######################  Push script and data to Synapse ######################

# Add this script to the code dir
syn_script <- synFindEntityId(
  "Pax6_DE_Targets.R",
  parent=syn_code_dir
)

if(is.null(syn_script)){
  syn_script <- File(
    path="scripts/Pax6_DE_Targets.R",
    parent=syn_code_dir
  )
  
  syn_script <- synStore(
    syn_script
  )
}

synSetProvenance(
  syn_script,
  activity = Activity(
    name = "extract_pax6_targets",
    description = "extracted DE Pax6 Targets",
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
