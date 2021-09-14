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
setwd("~/Documents/11Sep2021_Pax6_Study_DEG_Analysis")
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
# syn_deg_dir <- synFindEntityId("DEG_Tables", parent = syn_project)
# syn_dgelists <- synFindEntityId(
#   "pax6_master_dgelists.Rdata", 
#   parent=syn_data_dir
# )


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
      ifelse(grepl("LE", file), "Epi", "Fib"),
      levels = c("Epi", "Fib")
    ),
    ribo_filter = factor(
      ifelse(grepl("_rf_", file), "ribo", "geno"),
      levels = c("geno", "ribo", "degnorm")
    )
  ) %>% 
  group_by(genotype, cell_type, ribo_filter) %>%
  mutate(
    sample=paste(
      genotype, cell_type, 
      ribo_filter, row_number(), 
      sep="_")
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
  #tibble::rownames_to_column("label") %>%
  select(-group, -lib.size, -norm.factors) %>%
  mutate(
    label=gsub("_ribo_", "_degnorm_", label),
    ribo_filter="degnorm"
  )

row.names(samples) <- gsub(
  "_rf_GeneCount.txt", "_sorted_rf_alignment.bam",
  samples$file
)

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


