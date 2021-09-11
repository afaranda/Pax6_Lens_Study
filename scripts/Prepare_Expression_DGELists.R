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
rm(degnorm_results)

## Create DGEList Objects
save(
  list=c("pax6.master", "pax6.master.dgn"), 
  file="data/pax6_master_dgelists.Rdata"
)
