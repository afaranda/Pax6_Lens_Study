########################         HEADER BLOCK      ###########################
#  File:    Plot_Length_GC_Biases.R                                          #
#  Purpose: Evaluate Intersections between DEG Tables                        #
#  Created: Feb 2, 2021                                                      #
#  Author:  Adam Faranda                                                     #
#                                                                            #
##############################################################################

#########      Load Libraries and import expression data          ############
library(edgeR)
library(dplyr)
library(ggplot2)
library(venn)
library(clusterProfiler)
setwd("~/Documents/11Sep2021_Pax6_Study_DEG_Analysis/")
wd<-getwd()
data_dir="data"
source("scripts/Overlap_Comparison_Functions.R")


if(file.exists("results/pax6_deg_tables.Rdata")){
  load("results/pax6_deg_tables.Rdata")
  load("data/pax6_master_dgelists.Rdata")
  pax6.deg_master <- pax6.deg_master %>% mutate(
    Partition = factor(Partition),
    Test = factor(Test),
    Filtered = factor(Filtered)
  )
} else {
  source("scripts/Generate_Pax6_DEG_Tables.R")
  pax6.deg_master <- pax6.deg_master %>% mutate(
    Partition = factor(Partition),
    Test = factor(Test),
    Filtered = factor(Filtered)
  )
}

genes<-pax6.master$genes

########### Define dplyr Query to summarize bias in DEG Analyses #############
bias_summary<-function(
  deg_master, genes, minExp=0, minDiff=0, minLFC =1, maxSig = 0.05
){
  deg_master%>%
    filter(Avg1 > minExp | Avg2 > minExp) %>%
    filter(abs(Avg1 - Avg2) > minDiff) %>%
    inner_join(
      genes %>%
        select(gene_id, eu_length, eu_gc), 
      by="gene_id"
    ) %>%
    group_by(Test, Group_1, Group_2, Filtered, Partition) %>%
    summarize(
      Length_Pearson = cor(log2(eu_length), logFC),
      GC_Pearson = cor(log2(eu_gc), logFC),
      Length_Spearman = cor(log2(eu_length), logFC, method = "spearman"),
      GC_Spearman = cor(log2(eu_gc), logFC, method = "spearman"),
      Total_DEG = sum(abs(logFC) > minLFC & FDR < maxSig),
      UP_DEG = sum(logFC >  minLFC & FDR < maxSig),
      DOWN_DEG = sum(logFC < - minLFC & FDR < maxSig)
    ) %>% group_by() %>% as.data.frame()
  
}

####################### Define Functions to Plot biases ######################
length_bias_plot <- function(deg, fn){
  Contrast <- paste0(
    unique(deg$Group_2,), "v",
    unique(deg$Group_1)
  )
  
  mn<-paste("Length Bias in",Contrast)
  yl <-paste(                                    # Construct y axis label
    "Log2 Fold Change in",  unique(deg$Group_2,),          
    "vs", unique(deg$Group_1,)
  )
  
  png(fn, width=6, height = 5, units="in", res=1200)
  plot(
    x=log(deg$eu_length,2), y=deg$logFC, main=mn,
    xlab = "Log2 Gene Length (exon-union)",
    ylab = yl,
    col = ifelse(
      deg[,"Avg1"] == 0,
      "red",
      ifelse(deg[,"Avg2"] == 0, "red", "black")
    )
  )

  # Add Correlation Coefficients to tests
  ct=cor.test(deg$logFC, log(deg$eu_length,2), method="spearman")
  rho=paste("rho:", round(ct$estimate,3))
  sig=paste("p value:", signif(ct$p.value,3), sep="")
  abline(lm(deg$logFC~log(deg$eu_length,2)), col="red", lwd=2)
  text(7,max(deg$logFC)-1,rho)
  text(7,max(deg$logFC)-3,sig)
  dev.off()
}

gc_bias_plot <- function(deg, fn){
  Contrast <- paste0(
    unique(deg$Group_2,), "v",
    unique(deg$Group_1)
  )
  
  mn<-paste("GC Bias in",Contrast)
  yl <-paste(                                    # Construct y axis label
    "Log2 Fold Change in",  unique(deg$Group_2,),          
    "vs", unique(deg$Group_1,)
  )
  
  png(fn, width=6, height = 5, units="in", res=1200)
  plot(
    x=deg$eu_gc, y=deg$logFC, main=mn,
    xlab = "Fractional GC content (exon-union)",
    ylab = yl,
    col = ifelse(
      deg[,"Avg1"] == 0,
      "red",
      ifelse(deg[,"Avg2"] == 0, "red", "black")
    )
  )
  
  # Add Correlation Coefficients to tests
  ct=cor.test(deg$logFC, deg$eu_gc, method="spearman")
  rho=paste("rho:", round(ct$estimate,3))
  sig=paste("p value:", signif(ct$p.value,3), sep="")
  abline(lm(deg$logFC~deg$pl_gc), col="red", lwd=2)
  text(0.35,max(deg$logFC)-1,rho)
  text(0.35,max(deg$logFC)-3,sig)
  dev.off()
}

################# Iterate over studies and generate bias plots ################

# for(
#   study in list(
#     Pax6=list(
#       deg=pax6_deg_master %>% 
#       filter(
#         Filtered == 'ribo', 
#         Test == "ExactTest",
#         Partition == "Pair",
#       ) %>%
#       mutate(
#         Contrast = paste0(
#           Group_2, "v", Group_1
#         )
#       ) %>% inner_join(
# #         genes, by="gene_id"
#       ),
#       title="Pax6_Study"
#     )
# )
for(
  filt in unique(pax6.deg_master$Filtered)
){
  tests <- unique(
    pax6.deg_master %>% 
      filter(Filtered == filt) %>% 
      pull("Test")
  )
  for(test in tests){
    parts <- unique(
      pax6.deg_master %>% 
        filter(Filtered == filt & Test==test) %>% 
        pull("Partition")
    )
    for(part in parts){
      cns <- unique(
        pax6.deg_master %>%
          filter(
            Filtered == filt,
            Test == test,
            Partition == part
          ) %>%
          mutate(Contrast=paste0(Group_2,"v",Group_1)) %>%
          pull("Contrast")
      )
      print(paste(filt, test, part))
      for(cn in cns){
        deg <- pax6.deg_master%>%
          mutate(Contrast=paste0(Group_2,"v",Group_1)) %>%
          filter(
            Filtered == filt,
            Test == test,
            Partition == part,
            Contrast == cn
          ) %>%
          inner_join(
            genes,
            by="gene_id"
          )
        
        print(
          paste(
            cn,
            nrow(deg)
          )
        )
        filename = paste0(
          "results/LGB_",
          filt,"_",part,"_",test,"_",cn,
          "_Stat_Length_Bias_By_Partition.jpg"
        )
        print(filename)
        length_bias_plot(
          deg, 
          fn=filename
        )
        
        filename = paste0(
          "results/LGB_",
          filt,"_",part,"_",test,"_",cn,
          "_Stat_GC_Bias_By_Partition.jpg"
        )
        
        gc_bias_plot(
          deg, 
          fn=filename
        )
      }
    }
  }
}    

