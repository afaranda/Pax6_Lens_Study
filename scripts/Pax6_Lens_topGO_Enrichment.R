##############################################################################
#                                                                            #
#  File: Pax6_Lens_topGO_Enrichment.R                                        #
#  Author: Adam Faranda                                                      #
#  Created: June 18, 2020                                                    #
#  Purpose: Analyze DEG tables to identify enriched GO Terms and prepare     #
#           report elements including:                                       #
#                Tables of top enriched terms for each desired contrast      #
#                                                                            #
#                Bar-and-dot plots summarizing transcript abundance for      #
#                all genes and / or selected genes in a target category      #
#                across age dependent contrasts                              #
#                                                                            #
#                Tables of Statistics for all genes and / or selected genes  #
#                accross age dependent contrasts                             #
##############################################################################

################## Load Libraries and Source Dependencies ####################
library(DBI)
library(edgeR)
library(tidyr)
library(tibble)
library(AnnotationDbi)
library(clusterProfiler)
library(topGO)
library(ROntoTools)
library(dplyr)
library(tibble)
library(synapser)
options(echo=T)

# Enter Working Directory and Load Raw Data
setwd('~/Documents/11Sep2021_Pax6_Study_DEG_Analysis/')
# source('scripts/Excel_Write_Functions.R')
# source('scripts/Length_GC_Bias_Plot_Functions.R')
# source('scripts/Overlap_Comparison_Functions.R')
wd<-getwd()
results<-paste(wd,'results/topGO_Enrichment',sep='/')
if(!dir.exists(results)){
  dir.create(results)
}
data_dir <- paste0(wd,"/data")         ## Local data directory


################## Load in Master DGELists and DEG tables ####################

load('data/pax6_master_dgelists.Rdata')
load('results/pax6_deg_tables.Rdata')

#####################  Setup group and contrast iterators ####################


# Define A list of named contrasts; each element points to a vector with
# a pair of group labels. Positive fold changes will be associated
# with the second group listed. 
contrasts=list(
  WT_Fibers_vs_Epi=c('WTE', 'WTF', 'ribo', 'Pair'),
  P6_Fibers_vs_Epi=c('P6E', 'P6F', 'ribo', 'Pair'),
  Epi_P6_vs_WT=c('WTE', 'P6E','ribo', 'Pair'),
  Fibers_P6_vs_WT=c('WTF', 'P6F','ribo', 'Pair'),
  `2Way_Fib_vs_Epi`=c("Epi", "Fib", 'ribo', 'All'),
  `2Way_P6_vs_WT`=c("WT", "P6", 'ribo', 'All'),
  `2Way_Interaction`=c("REF", "INX", 'ribo', 'All') #,
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
contrast_descriptions=list(
  WT_Fibers_vs_Epi="Wildtype lens fiber cells vs lens epithelial cells",
  P6_Fibers_vs_Epi="Sey lens fiber cells vs lens epithelial cells",
  Epi_P6_vs_WT="Sey lens epithelial cells vs lens epithelial cells",
  Fibers_P6_vs_WT="Sey lens fiber cells vs lens fiber cells",
  `2Way_Fib_vs_Epi`="2 Way interaction model lens fiber cells vs lens fiber cells",
  `2Way_P6_vs_WT`="2 Way interaction model Sey vs Wildtype",
  `2Way_Interaction`="2 Way interaction model haploinsufficiency and differentiation"#,
  # DegNorm_WT_Fibers_vs_Epi="Wildtype lens fiber cells vs lens epithelial cells",
  # DegNorm_P6_Fibers_vs_Epi="Sey lens fiber cells vs lens epithelial cells",
  # DegNorm_Epi_P6_vs_WT="Sey lens epithelial cells vs lens epithelial cells",
  # DegNorm_Fibers_P6_vs_WT="Sey lens fiber cells vs lens fiber cells",
  # DegNorm_2Way_Fib_vs_Epi="2 Way interaction model lens fiber cells vs lens fiber cells",
  # DegNorm_2Way_P6_vs_WT="2 Way interaction model Sey vs Wildtype",
  # DegNorm_2Way_Interaction="2 Way interaction model haploinsufficiency and differentiation"
)

##### Iterate over contrasts and generate top DEG Tables
for( cn in names(contrasts)){
  
  ## Get Selection Criteria for contrast of interest
  g1 <- contrasts[[cn]][1]
  g2 <- contrasts[[cn]][2]
  norm <- contrasts[[cn]][3]
  part <- contrasts[[cn]][4]
  
  cat(contrast_descriptions[[cn]])
  cat("\n")
  cat(
    paste(
      paste("Group 1 Label:", g1), " --- ",
      paste("Group 2 Label:", g2)
    )
  )
  cat("\n")
  
  ## Get DE results for contrast of interest
  dg1 <- pax6.deg_master %>%
    dplyr::filter(
      Group_1 == g1 & Group_2 == g2 & 
        Test == "ExactTest" & Filtered == norm &
        Partition == part
    )
  cat(paste("Measured Genes:", norm, nrow(dg1), "\n"))
  
  
  ## Get named vector of statistically significant genes
  ag.stat <- ifelse(dg1[,"FDR"] < 0.05 & abs(dg1[,'logFC']) > 1, 1, 0)
  names(ag.stat)<-dg1[,"gene_id"]
  
  ## Get named vector of statistically significant genes
  ag.bio <- ifelse(
    dg1[,"FDR"] < 0.05 & abs(dg1[,'logFC']) > 1 & 
      abs(dg1$Avg1 - dg1$Avg2) > 2 & (dg1$Avg1 > 2 | dg1$Avg2 > 2), 1, 0
  )
  names(ag.bio)<-dg1[,"gene_id"]
  
  cat(
    "Statistically significant: ", sum(ag.stat),"\n",
    "Biologically significant: ", sum(ag.bio)
  )
  cat("\n\n")
  
  ## Iterate over filtering thresholds
  for(ag in list(list("stat", ag.stat), list("bio", ag.bio))){
    
    ## Iterate over target ontologies
    for(ont in c("BP", "MF", "CC")){
      
      ## Check whether enrichment tables have already been generated
      check <- c(
        paste0(
          cn,"_",ag[[1]],"_",norm,"_",part,"_topGO_",
          ont, "_classic_Enrichment.csv"
        ),
        paste0(
          cn,"_",ag[[1]],"_",norm,"_",part,"_topGO_",
          ont, "_elim_Enrichment.csv"
        ),
        paste0(
          cn,"_",ag[[1]],"_",norm,"_",part,"_topGO_",
          ont, "_weight_Enrichment.csv"
        )
      )
      if(! all(check %in% list.files(results))){
        GOdata <-new(
          "topGOdata", ontology=ont, allGenes = ag[[2]],
          geneSelectionFun = function(allScore){
            return(allScore > 0)
          }, annot=annFUN.org,
          mapping="org.Mm.eg.db",
          ID="ensembl"
        )
        #print(GOdata.stat)
        res.F.C<-runTest(GOdata, algorithm = "classic", statistic = "fisher")
        res.F.E<-runTest(GOdata, algorithm = "elim", statistic = "fisher")
        res.F.W<-runTest(GOdata, algorithm = "weight", statistic = "fisher")
        
        for(scoring in list(
          c("classic", "elim"),
          c("elim", "weight"),
          c("weight", "classic")
        )){
          fn <- paste(
            cn,"_",ag[[1]],"_",norm,"_",part,"_topGO_",ont,
            "_",scoring[1],"_Enrichment.csv", sep=""
          )
          if(!file.exists(fn)){
            
            allRes<-GenTable(
              GOdata, classic=res.F.C,
              elim=res.F.E, weight=res.F.W,
              orderBy=scoring[1], ranksOf=scoring[2],
              topNodes=100, numChar=100
            )
            print(paste("Scoring Method:", scoring[1]))
            print(allRes)
            
            write.csv(
              allRes, file=paste(results,fn, sep="/"),
              row.names = FALSE
            )
          } else {
            print("File Exists")
          }
        }
      } else {
        print("Files Exist")
      }
    }
  }
} 


################### Generate Gene Tables for selected Terms ####################
# mydb <- homology_sqlite <- "data/DEG_Homology_DB.sqlite" ## Set DB Path
# mydb <- dbConnect(
#   RSQLite::SQLite(),
#   homology_sqlite
# )
# 
# de_cross_tissues_dgn <- dbGetQuery(
#   mydb,
#   "SELECT * FROM hs_human_aging_accross_celltypes_degnorm;"
# )
# 
# de_cross_tissues_tmm <- dbGetQuery(
#   mydb,
#   "SELECT * FROM hs_human_aging_accross_celltypes_ribo;"
# )
# 
# ag <- rep(1, nrow(de_cross_tissues_dgn))
# names(ag) <- de_cross_tissues_dgn$ensembl_id
# GOdata <-new(
#   "topGOdata", ontology="BP", allGenes = ag,
#   geneSelectionFun = function(allScore){
#     return(allScore > 0)
#   }, annot=annFUN.org,
#   mapping="org.Hs.eg.db",
#   ID="ensembl"
# )


# term <- "GO:2000812"
# term <- "GO:0002088"
# term <- "GO:1990349"
# term <- "GO:0006000"
# term <- "GO:1905383"
# term <- "GO:0015793"
# term <- "GO:0006833"
# term <- "GO:0150094"
# term <- "GO:0007413"  ## Axonal Fasciulation
# term <- "GO:0097485"  ## Neuron Projection Guidance
# term <- "GO:0008333"  ## Endosome to lysosome transport
# term <- "GO:0048812"  ## Neuron Projection
# term <- "GO:0042552"  ## Myelination
# term <- "GO:0042775"  ## Mitochondrial ATP Synthesis Coupled electron transport
# term <- ""
# de_cross_tissues_dgn %>%
#   filter(ensembl_id %in% genesInTerm(GOdata, term)[[1]])  %>% View()
# 
# de_cross_tissues_dgn %>%
#   filter(grepl("ANK", symbol))  %>% View()
# 
# term <- "GO:0015250"
# 
# GOdata <-new(
#   "topGOdata", ontology="MF", allGenes = ag,
#   geneSelectionFun = function(allScore){
#     return(allScore > 0)
#   }, annot=annFUN.org,
#   mapping="org.Hs.eg.db",
#   ID="ensembl"
# )
# 
# 
# term <- "GO:0015250"
# term <- "GO:0004364"  ## Glutathione S-Transferase Activity
# de_cross_tissues %>%
#   filter(ensembl_id %in% genesInTerm(GOdata, term)[[1]]) %>% View()
# 
# ## Generate Tables for selected Gene Sets
# ### Genes associated with Alzheimers
# genes <- c(
#   "APOE", "CLU", "ABCA7", "SORL1", "CR1", "CD33", "MS4A",
#   "TREM2",  "CD2AP", "PICALM" ,"EPHA1", "HLA-DRB5",
#   "INPP5D","MEF2C","HLA-DRB1", "CASS4","PTK2B", "NME8", "ZCWPW1",
#   "CELF1","FERMT2","SLC24A4","DSG2", "PLD3","UNC5C","AKAP9","ADAM10",
#   "BIN1", "MME", "BACE1", "BACE2", "PSEN1", "PSEN2", "APP", "MMEL1", 
#   "IDE", "ECE1", "ECE2", "NEDD9"
# )
# 
# de_cross_tissues <- dbGetQuery(
#   mydb,
#   "SELECT * FROM hs_human_aging_accross_celltypes_degnorm;"
# )
# 
# de_cross_tissues %>%
#   filter(symbol %in% genes) %>% write.csv("results/topGO_Enrichment/Alzheimers_Genes_DegNorm.csv", row.names = F)
# 
# de_cross_tissues <- dbGetQuery(
#   mydb,
#   "SELECT * FROM hs_human_aging_accross_celltypes_ribo;"
# )
# 
# de_cross_tissues %>%
#   filter(symbol %in% genes)  %>% write.csv("results/topGO_Enrichment/Alzheimers_Genes_Ribo.csv", row.names = F)
# 
# 
# ### Genes associated with Neuron Architecture
# genes <- c(
#   "CNTN4", "TMOD1", "CAPRIN2", "NRCAM",
# )
# 
# de_cross_tissues <- dbGetQuery(
#   mydb,
#   "SELECT * FROM hs_human_aging_accross_celltypes_degnorm_bio;"
# )
# 
# de_cross_tissues %>%
#   filter(symbol %in% genes) %>% View()#%>% write.csv("results/topGO_Enrichment/Alzheimers_Genes_DegNorm.csv", row.names = F)
# 
# de_cross_tissues <- dbGetQuery(
#   mydb,
#   "SELECT * FROM hs_human_aging_accross_celltypes_ribo;"
# )
# 
# de_cross_tissues %>%
#   filter(symbol %in% genes)  %>% write.csv("results/topGO_Enrichment/Alzheimers_Genes_Ribo.csv", row.names = F)
# # Genes associated with the pathway "Neuroactive Ligand Receptor"
# genes <- c(
#   "GABRA3", "CHRNA2", "S1PR1", "HRH1", "NPFFR1",
#   "VIPR1", "APLN", "MC1R", "C3", "ADM", "S1PR3", 
#   "GABRE", "GRIN3A", "NPY1R", "ADRB2", "P2RY6", 
#   "PRLR", "GRIN2A", "GRIA4", "GRIN2C", "EDN1", "MTNR1A",
#   "GRIK2", "GPR35","GRM3", "TBXA2R", "GRIN2D", "F2R", "NDRG2"
#   
# )
# de_cross_tissues %>%
#   filter(symbol %in% genes) %>% View()
# # 
# # 
# # 
# 
# dg1 <- pax6.deg_master %>%
#   dplyr::filter(
#     #Group_1 == "YC0" & Group_2 == "AC0" & 
#       Test == "ExactTest" & Filtered == "degnorm"
#   )
# ag<-ifelse(dg1[,"FDR"] < 0.05 & abs(dg1[,'logFC']) > 1, 1, 0)
# names(ag)<-dg1[,"gene_id"]
# GOdata<-new(
#   "topGOdata", ontology="BP", allGenes = ag,
#   geneSelectionFun = function(allScore){
#     return(allScore > 0)
#   }, annot=annFUN.org,
#   mapping="org.Hs.eg.db",
#   ID="ensembl"
# )
# 
# # Run Enrichment Tests #######################################################
# res.F.C<-runTest(GOdata, algorithm = "classic", statistic = "fisher")
# res.F.E<-runTest(GOdata, algorithm = "elim", statistic = "fisher")
# res.F.W<-runTest(GOdata, algorithm = "weight", statistic = "fisher")
# 
# dg1.allRes<-GenTable(
#   GOdata, classFisher=res.F.C, 
#   elimFisher=res.F.E, weightFisher=res.F.W,
#   orderBy="classFisher", ranksOf="classFisher",
#   topNodes =100, numChar=100
# )
# 
# ### Get Genes associated with a specific term
# pax6.master.dgn$genes[genesInTerm(GOdata, "GO:0034330")[[1]],c("gene_id","SYMBOL", "DESCRIPTION")]
# ####### Try Cluster Profiler Methods
# library(org.Hs.eg.db)
# library(AnnotationHub)
# ah <- AnnotationHub()
# ens <- ah[['AH95744']]
# library(clusterProfiler)
# 
# data(geneList, package="DOSE")
# gene <- names(geneList)[abs(geneList) > 2]
# 
# # Entrez gene ID
# head(gene)
# 
# 
# 
# ggo <-enrichGO(
#   gene = dg1 %>%
#     filter(abs(logFC) > 1 & FDR < 0.05) %>%
#     #filter(abs(Avg1 - Avg2) > 2) %>%
#     pull(gene_id),
#   OrgDb    = org.Hs.eg.db,
#   keyType = "ENSEMBL",
#   ont      = "BP",
#   readable = TRUE
# )
# 
# nrow(ggo@result %>% dplyr::select( -geneID) %>% filter(qvalue< 0.05))
# ggo@result  %>% filter(qvalue< 0.05)%>% dplyr::select( qvalue)
# 
# ggo@result %>% dplyr::select( -geneID) %>% filter(qvalue< 0.05)
# 
# 
# comp <- inner_join(
#   ggo@result %>% select(ID, clusterProfiler = pvalue),
#   dg1.allRes %>% select(ID = GO.ID, topGO=classFisher) %>% mutate(
#     topGO = p.adjust(topGO, method = "none")
#   ),
#   by="ID"
# )
# 
# 
# comp <- inner_join(
#   ggo@result %>% select(ID, clusterProfiler = GeneRatio) %>%
#     mutate(clusterProfiler = as.integer(sub("\\/.*", "", clusterProfiler))),
#   dg1.allRes %>% select(ID = GO.ID, topGO=Significant),
#   by="ID"
# )
# 
# 
# comp <- inner_join(
#   ggo@result %>% select(ID, clusterProfiler = BgRatio) %>%
#     mutate(clusterProfiler = as.integer(sub("\\/.*", "", clusterProfiler))),
#   dg1.allRes %>% select(ID = GO.ID, topGO=Annotated),
#   by="ID"
# )
# 
# ###### ROntoTools Analysis  ######
# library(ROntoTools)
# kpg <- keggPathwayGraphs("hsa", verbose = FALSE)
# kpn <- keggPathwayNames("hsa")
# 
# kpg <- setEdgeWeights(
#   kpg, edgeTypeAttr = "subtype",
#   edgeWeightByType = list(
#     activation = 1, inhibition = -1,
#     expression = 1, repression = -1),defaultWeight = 0
# )
# 
# 
# load(system.file("extdata/E-GEOD-21942.topTable.RData", package = "ROntoTools"))
# fc <- top$logFC[top$adj.P.Val <= .01]
# names(fc) <- top$entrez[top$adj.P.Val <= .01]
# pv <- top$P.Value[top$adj.P.Val <= .01]
# names(pv) <- top$entrez[top$adj.P.Val <= .01]
# head(fc)
# ref <- top$entrez
# 
# kpg <- setNodeWeights(kpg, weights = alphaMLG(pv), defaultWeight = 1)
# peRes <- pe(x=fc, ref=ref, graphs=kpg)
# summary(peRes, kpn=kpn)
# head(genes)
# 
# x <- data.frame(
#   FPKM = c(10,20,15, 12, 8, 10, 6, 7,5, 11,18,21, 13, 9, 11, 2, 1,5),
#   CellType=rep(rep(c("Rh", "Eq", "Cf"), each =3), times=2),
#   Age=rep(c("Young", "Aged"), each =9)
# )
# x$CellType <- as.factor(x$CellType)
# levels(x$CellType) <- c("Central LEC", "Equatorial LEC", "Cortical Fibers")
# ################## Plotting Function for Gene accross age and cell type
# 
# plot_single_gene_abundance <- function(fpkm_data){
#   ## Accepts a long format data frame with the columns Age, CellType 
#   ## and FPKM.  Age and CellType should be ordered factors.  Returns
#   ## A ggplot object that can be printed or saved
#   
#   ## Calculate group means for bars
#   group_means <- fpkm_data %>%
#     dplyr::group_by(Age, CellType) %>%
#     summarise(FPKM = mean(FPKM))
#   
#   ## Generate Bar Plot for Group means
#   p <- ggplot(group_means) + geom_bar(
#     mapping=aes(fill=CellType, y=FPKM, x=Age),
#     stat = "identity", position="dodge", alpha=0.8
#   ) + 
#     scale_fill_grey() + 
#     ## Add FPKM values for individual samples
#     geom_point(
#       data = fpkm_data,
#       aes(
#         x=Age,
#         y=FPKM,
#         group = CellType
#       ),
#       size=4,
#       position=position_dodge(width=0.9)
#     ) + 
#     theme(panel.background = element_rect(fill="white"))
#   return(p)
# }
# p <- plot_single_gene_abundance(x)
# 
# 
# 
# y <- bind_rows(
#   x %>% mutate(
#     gene_id ="FN1", 
#     FPKM = FPKM *1
#   ),
#   x %>% mutate(
#     gene_id ="ACTA2", 
#     FPKM = FPKM *10
#   ),
#   x %>% mutate(
#     gene_id ="TNC", 
#     FPKM = FPKM *50
#   ),
#   x %>% mutate(
#     gene_id ="CXCL1", 
#     FPKM = FPKM *100
#   )
# )
# 
# fpkm_data <- y
# group_means <- fpkm_data %>%
#   dplyr::group_by(Age, CellType, gene_id) %>%
#   summarise(FPKM = mean(FPKM))
# 
# 
# 
# p <- ggplot(
#   group_means
# ) + geom_bar(
#   mapping=aes(fill=CellType, y=FPKM, x=Age),
#   stat = "identity", position="dodge", alpha=0.85
# ) + scale_fill_grey() + geom_point(
#   data = fpkm_data,
#   aes(
#     x=Age,
#     y=FPKM,
#     group = CellType
#   ),
#   size=2,
#   position=position_dodge(width=0.9)
# )+ facet_grid(
#   rows = vars(gene_id), scale="free"
# ) + scale_color_grey() + theme(panel.background = element_rect(fill="white"))
# 
# 
# 
# 
# 
# 
# 
# ggarrange(p,q)
# 
# 
# library(ggridges)
# library(ggplot2)
# 
# # Diamonds dataset is provided by R natively
# #head(diamonds)
# 
# # basic example
# ggplot(diamonds, aes(x = price, y = cut, fill = cut)) +
#   geom_density_ridges() +
#   theme_ridges() + 
#   theme(legend.position = "none")
# 
# 


