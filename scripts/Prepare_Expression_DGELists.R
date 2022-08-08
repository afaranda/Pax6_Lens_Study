##############################################################################
#                                                                            #
#  File:    Prepare_Expression_DGELists.R                                    #
#  Purpose: Import expression measurements and generate DGEList opjects      #
#  Created: Feb 2, 2021                                                      #
#  Author:  Adam Faranda                                                     #
#                                                                            #
##############################################################################

# Load libraries, set working directory
library(edgeR)
library(dplyr)
#library(synapser)
setwd("~/Documents/23Feb2022_Pax6_Study_DEG_Analysis/")
source("scripts/synapse_reticulate_wrapper.R")
wd<-getwd()
data_dir="data"

# Import table of gene annotations 
# Import Annotations, append to "Gene Length" Table
fn<-paste(data_dir,'Gene_Annotations.csv', sep='/')
lt<-read.csv(fn, row.names = 1)
if(!"ENTREZID" %in% names(lt)){
  library(org.Mm.eg.db)
  lt<-merge(
    lt, AnnotationDbi::select(
      org.Mm.eg.db, keys=unique(lt$SYMBOL), 
      columns = c("ENTREZID"), 
      keytype = "SYMBOL"
    ), by='SYMBOL',all.x=T
  )
  row.names(lt)<-lt$gene_id
  detach(package:org.Mm.eg.db, unload = T)
  detach(package:AnnotationDbi, unload=T)
  write.csv(lt, fn)
} 


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
syn_hpco_gene_meta <- synFindEntityId(
  "Gene_Metadata",
  synFindEntityId("Human_PCO_Study_Data_Set_1")
)


##############################################################################
#                   Import metadata to annotate genes table                  #
##############################################################################


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

isyte528 %>%
  filter(Platform == "affy430") %>%
  rowwise() %>%
  mutate(
    IS_ISYTE_P56 = (
      Interval=="P56" & fold_change > 1 & p_value < 0.05
    )
  ) %>% as.data.frame() -> isyte528
############################# Load Zonule Data ###############################

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

######################### Prepare Pax6 Peak data ##############################

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


##############################################################################
#                        Create DGEList for Pax6 Study                       #
##############################################################################

# Build phenotype / attribute table for samples
samples<-data.frame(file=list.files(data_dir, pattern="GeneCount"))%>%
  mutate(
    genotype = factor(
      ifelse(grepl("(plus_minus|Het)", file),"P6", "WT"),
        levels = c("WT", "P6")
    ),
    cell_type = factor(
      ifelse(grepl("LE", file), "E", "F"),
      levels = c("E", "F")
    ),
    ribo_filter = factor(
      ifelse(grepl("_rf_", file), "ribo", "geno"),
      levels = c("geno", "ribo", "degnorm")
    )
  ) %>% 
  group_by(genotype, cell_type, ribo_filter) %>%
  mutate(
    sample=paste(
      genotype, cell_type,"_", 
      ribo_filter,"_", row_number(), 
      sep="")
  ) %>% 
  arrange(ribo_filter, cell_type, genotype) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("sample")
  samples$label <- row.names(samples)

# Aggregate Pax6 Study Samples
pax6.master<-readDGE(
  path=data_dir, 
  files=samples$file,
  labels=row.names(samples),
  header=FALSE
)

# Assemble DGEList Object
pax6.master<-DGEList(
  counts=pax6.master$counts[row.names(lt), ],
  samples=samples[row.names(pax6.master$samples),],
  genes=lt
)
load(paste0(data_dir,"/","DegNorm_Analysis_Results_ribo.Rdata"))
 
samples <- pax6.master$samples %>%
  filter(ribo_filter=="ribo") %>%
  select(-group, -lib.size, -norm.factors) %>%
  mutate(
    label=gsub("_ribo_", "_degnorm_", label),
    ribo_filter="degnorm"
  )

row.names(samples) <- gsub(
  "_rf_GeneCount.txt", "_sorted_rf_alignment.bam",
  samples$file
)

################ Add Annotation flags to gene metadata table #################

pax6.master$genes$IS_ISYTE_P56 <- pax6.master$genes$SYMBOL %in% (
  isyte528[isyte528$IS_ISYTE_P56,"MGI.symbol"]
)

pax6.master$genes$IS_ZONULE <- pax6.master$genes$gene_id %in% ( 
  zonules$Gene.stable.ID
)

pax6.master$genes$IS_TRRUST_PAX6_TARGET <- pax6.master$genes$SYMBOL %in% (
  trrust_targets$Target
)

pax6.master$genes$IS_SUN_PAX6_TARGET <- pax6.master$genes$SYMBOL %in% (
  sun_targets$Target
)

pax6.master$genes$IS_SUN_PAX6_LENS_PEAK <- pax6.master$genes$SYMBOL %in% (
  sun_lens_peaks$Target
)

pax6.master$genes$IS_SUN_PAX6_FOREBRAIN_PEAK <- (
  pax6.master$genes$SYMBOL %in% sun_lens_peaks$Target
)

#pax6.master$genes$SUN_LENS_PROMOTER <- 0
#pax6.master$genes$SUN_LENS_EXON <- 0
#pax6.master$genes$SUN_LENS_INTRON <- 0
#pax6.master$genes$SUN_LENS_DISTAL <- 0



####################  Create master for DegNorm counts #######################
counts <- degnorm_results$counts_normed[,row.names(samples)]
genes <- pax6.master$genes[row.names(counts),]
pax6.master.dgn <- DGEList(
  counts=counts,
  genes = genes,
  samples=samples
)
colnames(pax6.master.dgn) <- pax6.master.dgn$samples$label


### Save counts and DI's extracted from the DegNorm results object
write.csv(
  counts,
  file=paste0(data_dir,"/","Pax6_Lens_DegNorm_Counts.csv")
)

write.csv(
  degnorm_results$DI,
  file=paste0(data_dir,"/","Pax6_Lens_DegNorm_DI.csv")
)

rm(degnorm_results)

## Create DGEList Objects
save(
  list=c("pax6.master", "pax6.master.dgn"), 
  file="data/pax6_master_dgelists.Rdata"
)

######################  Push script and data to Synapse ######################
script_path <- "scripts/Prepare_Expression_DGELists.R"
dgelists_path <- "data/pax6_master_dgelists.Rdata"

### Upload Count Files and DegNorm Counts / DI Matrices
### used bythis script
syn_count_files <- list()
for(f in c(
  list.files(data_dir, pattern="_GeneCount.txt"),
  "Pax6_Lens_DegNorm_Counts.csv",
  "Pax6_Lens_DegNorm_DI.csv")
){
  count_file <- File(
    path=paste0(data_dir,"/",f),
    parent=syn_data_dir
  )
  syn_count_files <- append(
    syn_count_files,
    synStore(
      count_file
    )
  )
}

## Upload this script
syn_script <- File(
  path=script_path,
  parent=syn_code_dir
)

syn_script <- synStore(
  syn_script,
  used = syn_count_files
)

## Upload this script's results
syn_dgelists <- File(
  path=dgelists_path,
  parent=syn_data_dir
)

syn_dgelists <- synStore(
  syn_dgelists,
  executed = syn_script
)


