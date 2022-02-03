################################################################################
# File: Overlap_Comparison_Functions.R                                         #
# Author: Adam Faranda                                                         #
# Created: June 27, 2019                                                       #
# Purpose: Compare differential expression results from various ocular tissues #
#          in pairwise contrasts between wildtype and pax6 heterozygous mice   #
#                                                                              #
################################################################################

# Setup Workspace
library('openxlsx')
library('dplyr')
wd<-getwd()

# Given a data frame and a vector of column names; return a vector
# of each column's corresponding number in the data frame. 
colNum<-function(df, n){
  if(is.character(n)){
    positions<-c()
    for(i in n){
      positions<-c(
        positions,
        grep(i, colnames(df))
      )
    }
  } else {
    positions<-n
  }
  positions
}

# Function recieves a data frame and a list of name pairs, returns a data frame
# with names updated based on the list
dfRename<-function(df, name_pairs){
  if (!is.null(names(name_pairs))){
    for(np in names(name_pairs)){
      names(df)[grep(name_pairs[np], names(df))]<-np
    }
  }
  df
}

dfSubname<-function(df, name_pairs, pf="^", po="$"){
  if (!is.null(names(name_pairs))){
    for(np in names(name_pairs)){
      r<-paste(pf,name_pairs[np], po,sep="")
      names(df)<-gsub(r, np, names(df))
    }
  }
  df
}

# Query for Pairwise overlap between two Deglists -- assumes both tables
# have the same column headers, and that "MGI.symbol" is a unique id in both
# lists
query<-function(
  dg1=ss_master, dg2=ss_master, id_col="MGI.symbol",
  cols=c("logFC", "PValue", "FDR", "Group_1", "Group_2", "Avg1", "Avg2")
){
 dg1<-dg1[,c(id_col, cols)]
 dg2<-dg2[,c(id_col, cols)]
 
 for(c in cols){
   names(dg1)[grep(c, names(dg1))]<-paste("dg1.",c, sep="")
   names(dg2)[grep(c, names(dg1))]<-paste("dg2.",c, sep="")
 }
 out<-full_join(dg1, dg2, by=id_col)
 out
}


# Function takes a set of DEGs tablutes Total, Up and Down Genes
# For statistical and biological significance
degSummary<-function(
  df, lfc_min=1, fdr_max=0.05, 
  Avg1="Avg1", Avg2="Avg2",
  lfc="logFC", fdr="FDR", 
  bioFun = bioSigRNASeq, minExp = 2
){
  names(df)[grep(lfc, names(df))]<-"lfc"
  names(df)[grep(fdr, names(df))]<-"fdr"
  names(df)[grep(Avg1, names(df))]<-"Avg1"
  names(df)[grep(Avg2, names(df))]<-"Avg2"
  
  dg<-data.frame(
    criteria = c("Statistically Significant", "Biologically Significant"),
    Total = c(
      nrow(df %>% filter(abs(lfc) > lfc_min, fdr < fdr_max)), 
      nrow(
      	bioFun(df, lfc ="lfc", stat="fdr", a1="Avg1", a2="Avg2", 
  		minExp= minExp, maxStat=fdr_max, minLfc=lfc_min)
      )
    ),
    Up = c(
      nrow(df %>% filter(lfc > lfc_min, fdr < fdr_max)), 
      nrow(
      	bioFun(df, lfc ="lfc", stat="fdr", a1="Avg1", a2="Avg2", 
  		minExp= minExp, maxStat=fdr_max, minLfc=lfc_min) %>%
        filter(lfc > lfc_min)
      )
    ),
    Down = c(
      nrow(df %>% filter(lfc < -lfc_min, fdr < fdr_max)),
      nrow(
      	bioFun(df, lfc ="lfc", stat="fdr", a1="Avg1", a2="Avg2", 
  		minExp= minExp, maxStat=fdr_max, minLfc=lfc_min) %>%
        filter(lfc < lfc_min)
      ) 
    )
  )
  dg
}

# Extract directional subsets of statistically significant genes
subsetTables<-function(
  df,                       # Data frame with a joined pair of results
  id_col="MGI.symbol",      # Unique Identifier for this gene
  Contrast_1 = "LE",        # Name of the first contrast in df
  Contrast_2 = "PCO",       # Name of the second contrast in df
  dg1="dg1",                # Prefix for contrast 1
  dg2="dg2",                # Prefix for contrast 2
  lfc="logFC",              # Column with log 2 fold change values
  pvl="p_value",            # Column with p value for pairwise test
  fdr="FDR",                # Column with FDR values
  g1 = "Group_1",           # Column with Group_1 label
  g2 = "Group_2",           # Column with Group_2 label
  a1 = "Avg1",              # Column with average values for Group_1
  a2 = "Avg2",              # Column with average values for Group_2
  stat = T,                 # Whether to use 'Stat' or 'Bio' naming scheme
  unlog = T,                # Whether to report absolute or log2 fold changes
  descname = F,             # Use original, or descriptive attribute names
  annot = NULL,             # Optionally provide table (keyed on ID)
  dropGroup = T,            # Drop or keep group label columns
  unit = "Avg"              # Aggregated unit for groupwise abundance
){
  # Standardize column headers
  cols<-c(lfc=lfc, pvl=pvl, fdr=fdr, g1=g1, g2=g2, a1=a1, a2=a2)
  df<-dfSubname(df, cols, pf="")
  
  prf<-c(dg1=dg1, dg2=dg2)
  df<-dfSubname(df, prf, po="")
  
  if(unlog){
    df$dg1.lfc<-ifelse(df$dg1.lfc >= 0, 2^df$dg1.lfc, -1/2^df$dg1.lfc)
    df$dg2.lfc<-ifelse(df$dg2.lfc >= 0, 2^df$dg2.lfc, -1/2^df$dg2.lfc)
  }
  if(!is.null(annot)){
    if(any(grepl(id_col, names(annot)))){
       if(length(unique(annot[,id_col])) == nrow(annot)){
        acols<-setdiff(names(annot), id_col)
        dcols<-setdiff(names(df), id_col)
        df<-left_join(df, annot, by=id_col)
        df<-data.frame(df, stringsAsFactors = F)
        df<-df[, c(id_col, acols, dcols)]
      }
    }
  }
  # Subset results
  if (stat){
    tables<-list(
      `Stat Sig Intersection`=df,
      `SS Dn C1 Up C2` = df %>% filter(dg1.lfc < 0, dg2.lfc > 0 ),
      `SS Dn C1 Dn C2` = df %>% filter(dg1.lfc < 0, dg2.lfc < 0 ),
      `SS Up C1 Up C2` = df %>% filter(dg1.lfc > 0, dg2.lfc > 0 ),
      `SS Up C1 Dn C2` = df %>% filter(dg1.lfc > 0, dg2.lfc < 0 )
    )
    
  # Use if the data tables submitted via df are biologically significant
  } else {
    tables<-list(
      `Bio Sig Intersection`=df,
      `BS Dn C1 Up C2`= df %>% filter(dg1.lfc < 0, dg2.lfc > 0 ),
      `BS Dn C1 Dn C2`= df %>% filter(dg1.lfc < 0, dg2.lfc < 0 ),
      `BS Up C1 Up C2`= df %>% filter(dg1.lfc > 0, dg2.lfc > 0 ),
      `BS Up C1 Dn C2`= df %>% filter(dg1.lfc > 0, dg2.lfc < 0 )
    )
  }
  
  # Rename Contrasts
  names(tables)<-gsub("C1", Contrast_1, names(tables))
  names(tables)<-gsub("C2", Contrast_2, names(tables))
  if(descname){
    cols<-c(
      dg1.g1 = paste(Contrast_1, df$dg1.g1[1], sep="_"),
      dg1.g2 = paste(Contrast_1, df$dg1.g2[1], sep="_"),
      dg2.g1 = paste(Contrast_2, df$dg2.g1[1], sep="_"),
      dg2.g2 = paste(Contrast_2, df$dg2.g2[1], sep="_"),
      dg1.a1 = paste(Contrast_1, df$dg1.g1[1], unit, sep="_"),
      dg1.a2 = paste(Contrast_1, df$dg1.g2[1], unit, sep="_"),
      dg2.a1 = paste(Contrast_2, df$dg2.g1[1], unit, sep="_"),
      dg2.a2 = paste(Contrast_2, df$dg2.g2[1], unit, sep="_"),
      dg1.lfc = ifelse(
        unlog, paste(Contrast_1, "Fold_Change",sep="_"),
        paste(Contrast_1, lfc,sep="_")
      ),
      dg2.lfc=ifelse(
        unlog, paste(Contrast_2, "Fold_Change",sep="_"),
        paste(Contrast_2, lfc,sep="_")
      ),
      dg1.pvl = paste(Contrast_1, pvl,sep="_"),
      dg2.pvl = paste(Contrast_2, pvl,sep="_"),
      dg1.fdr = paste(Contrast_1, fdr,sep="_"),
      dg2.fdr = paste(Contrast_2, fdr,sep="_")
    )
    
    for (t in names(tables)){
      if(dropGroup){
        sloc<-names(cols)
        names(sloc)<-cols
        keep<-setdiff(
          names(tables[[t]]), 
          sloc[c(grep("\\.g1$", sloc), grep("\\.g2$", sloc))]
        )
        sloc<-sloc[-c(grep("\\.g1$", sloc), grep("\\.g2$", sloc)) ]
        # print(keep)
        tables[[t]]<-tables[[t]][, keep]
        tables[[t]]<-dfSubname(tables[[t]], sloc)
        
      } else {
        sloc<-names(cols)
        names(sloc)<-cols
        tables[[t]]<-dfSubname(tables[[t]], sloc)
      }
    }
    
  } else {
    for(t in names(tables)){
      if(dropGroup){
        sloc<-names(cols)
        frp<-names(prf)
        names(sloc)<-cols
        names(frp)<-prf
        
        drop<-c(
          grep("\\.g1$", names(tables[[t]])), 
          grep("\\.g2$", names(tables[[t]]))
        )
        keep<-names(tables[[t]])[-drop]
        # print(keep)
        sloc<-sloc[-c(grep("g1$", sloc), grep("g2$", sloc))]
        tables[[t]]<-dfSubname(tables[[t]], sloc)
        tables[[t]]<-dfSubname(tables[[t]], frp)
        
      } else {
        sloc<-names(cols)
        frp<-names(prf)
        names(sloc)<-cols
        names(frp)<-prf
        tables[[t]]<-dfSubname(tables[[t]], sloc)
        tables[[t]]<-dfSubname(tables[[t]], frp)
      }
    }
  }
  tables<-append(
    tables, list(Contrasts=c(Contrast_1=Contrast_1, Contrast_2=Contrast_2))
  )
  tables
}

################################################################

# Extract directional subsets of statistically significant genes
compareHits<-function(
  df,                       # Data frame with a joined pair of results
  id_col="MGI.symbol",      # Unique Identifier for this gene
  Contrast_1 = "LE",        # Name of the first contrast in df
  Contrast_2 = "PCO",       # Name of the second contrast in df
  dg1="dg1",                # Prefix for contrast 1
  dg2="dg2",                # Prefix for contrast 2
  dg1.unit="Avg_FPKM",      # Units of expression measure for dg1
  dg2.unit="Avg_FPKM",      # Units of expression measure for dg2
  lfc="logFC",              # Column with log 2 fold change values
  pvl="p_value",            # Column with p value for pairwise test
  fdr="FDR",                # Column with FDR values
  g1 = "Group_1",           # Column with Group_1 label
  g2 = "Group_2",           # Column with Group_2 label
  a1 = "Avg1",              # Column with average values for Group_1
  a2 = "Avg2",              # Column with average values for Group_2
  stat = T,                 # Whether to use 'Stat' or 'Bio' naming scheme
  unlog = T,                # Whether to report absolute or log2 fold changes
  descname = F,             # Use original, or descriptive attribute names
  annot = NULL,             # Optionally provide table (keyed on ID)
  dropGroup = T,            # Drop or keep group label columns
  dg1.lfcmin = 1,           # Fold Change threshold in dg1
  dg2.lfcmin = 1,           # Fold Change threshold in dg2
  minDiff = 2,              # Minimum difference for biological significance
  minAvg  = 2               # Minimum Abundance for biologuical significance
){
  # print(length(stat))

  # Standardize column headers
  cols<-c(lfc=lfc, pvl=pvl, fdr=fdr, g1=g1, g2=g2, a1=a1, a2=a2)
  df<-dfSubname(df, cols, pf="")
  
  prf<-c(dg1=dg1, dg2=dg2)
  df<-dfSubname(df, prf, po="")
  
  if(!is.null(annot)){
    # print(head(annot))
    if(any(grepl(id_col, names(annot)))){
      print("found ID")
      if(length(unique(annot[,id_col])) == nrow(annot)){
        acols<-setdiff(names(annot), id_col)
        dcols<-setdiff(names(df), id_col)
        df<-left_join(df, annot, by=id_col)
        df<-data.frame(df, stringsAsFactors = F)
        df<-df[, c(id_col, acols, dcols)]
      } else {
        print("Length Mismatch")
      }
    }
  }
  # Subset results
  # print(length(stat))
  if (stat){
    tables<-list(
      `Stat Sig Intersection`= df %>%
        filter(abs(dg1.lfc) > dg1.lfcmin & dg1.fdr < 0.05)%>%
        filter(abs(dg2.lfc) > dg2.lfcmin & dg2.fdr < 0.05),
      `SS Dn C1 Up C2` = df %>% 
        filter(
          dg1.lfc < 0-dg1.lfcmin & dg1.fdr < 0.05, 
          dg2.lfc > dg2.lfcmin & dg2.fdr < 0.05),
      
      `SS Dn C1 Dn C2` = df %>% 
        filter(
          dg1.lfc < 0-dg1.lfcmin & dg1.fdr < 0.05, 
          dg2.lfc < 0-dg2.lfcmin & dg2.fdr < 0.05),
      
      `SS Up C1 Up C2` =df %>% 
        filter(dg1.lfc > dg1.lfcmin & dg1.fdr < 0.05, 
               dg2.lfc > dg2.lfcmin & dg2.fdr < 0.05),
      
      `SS Up C1 Dn C2` = df %>% 
        filter(dg1.lfc > dg1.lfcmin & dg1.fdr < 0.05, 
               dg2.lfc < 0-dg2.lfcmin & dg2.fdr < 0.05),
      
      `SS Up C1 Only`  = df %>% 
        filter(dg1.lfc > dg1.lfcmin & dg1.fdr < 0.05) %>%
        filter(!((abs(dg2.lfc) > dg2.lfcmin & dg2.fdr < 0.05)) | is.na(dg2.lfc)),
      
      `SS Dn C1 Only`  = df %>% 
        filter(dg1.lfc < 0-dg1.lfcmin & dg1.fdr < 0.05) %>% 
        filter(!((abs(dg2.lfc) > dg2.lfcmin & dg2.fdr < 0.05)) | is.na(dg2.lfc)),    
      
      `SS Up C2 Only`  = df %>% 
        filter(dg2.lfc > dg2.lfcmin & dg2.fdr < 0.05) %>%
        filter(!((abs(dg1.lfc) > dg1.lfcmin & dg1.fdr < 0.05)) | is.na(dg1.lfc)),
      
      `SS Dn C2 Only`  = df %>% 
        filter(dg2.lfc < 0-dg2.lfcmin & dg2.fdr < 0.05) %>%
        filter(!((abs(dg1.lfc) > dg1.lfcmin & dg1.fdr < 0.05)) | is.na(dg1.lfc))
      
    )
    
    # Use if the data tables submitted via df are biologically significant
  } else {
    dg <- df
    df <- df %>% 
      filter(abs(dg1.a1 - dg1.a2) > minDiff) %>%
      filter(dg1.a1 > minAvg | dg1.a2 > minAvg) %>%
      filter(abs(dg2.a1 - dg2.a2) > minDiff) %>%
      filter(dg2.a1 > minAvg | dg2.a2 > minAvg)
    
    tables<-list(
      # `Bio Sig Intersection`=df,
      # `BS Dn C1 Up C2`= df %>% filter(dg1.lfc < 0, dg2.lfc > 0 ),
      # `BS Dn C1 Dn C2`= df %>% filter(dg1.lfc < 0, dg2.lfc < 0 ),
      # `BS Up C1 Up C2`= df %>% filter(dg1.lfc > 0, dg2.lfc > 0 ),
      # `BS Up C1 Dn C2`= df %>% filter(dg1.lfc > 0, dg2.lfc < 0 )
      
      
      `Bio Sig Intersection`=df %>%
        filter(abs(dg1.lfc) > dg1.lfcmin & dg1.fdr < 0.05)%>%
        filter(abs(dg2.lfc) > dg2.lfcmin & dg2.fdr < 0.05),
      `BS Dn C1 Up C2` = df %>% 
        filter(
          dg1.lfc < 0-dg1.lfcmin & dg1.fdr < 0.05, 
          dg2.lfc > dg2.lfcmin & dg2.fdr < 0.05),
      
      `BS Dn C1 Dn C2` = df %>% 
        filter(
          dg1.lfc < 0-dg1.lfcmin & dg1.fdr < 0.05, 
          dg2.lfc < 0-dg2.lfcmin & dg2.fdr < 0.05),
      
      `BS Up C1 Up C2` =df %>% 
        filter(dg1.lfc > dg1.lfcmin & dg1.fdr < 0.05, 
               dg2.lfc > dg2.lfcmin & dg2.fdr < 0.05),
      
      `BS Up C1 Dn C2` = df %>% 
        filter(dg1.lfc > dg1.lfcmin & dg1.fdr < 0.05, 
               dg2.lfc < 0-dg2.lfcmin & dg2.fdr < 0.05),
      
      `BS Up C1 Only`  = df %>% 
        filter(dg1.lfc > dg1.lfcmin & dg1.fdr < 0.05, 
               !(abs(dg2.lfc) > dg2.lfcmin & dg2.fdr < 0.05)),
      
      `BS Dn C1 Only`  = df %>% 
        filter(dg1.lfc < 0-dg1.lfcmin & dg1.fdr < 0.05, 
               !(abs(dg2.lfc) > dg2.lfcmin & dg2.fdr < 0.05)),
      
      `BS Up C2 Only`  = df %>% 
        filter(!(abs(dg1.lfc) > dg1.lfcmin & dg1.fdr < 0.05), 
               dg2.lfc > dg2.lfcmin & dg2.fdr < 0.05),
      
      `BS Dn C2 Only`  = df %>% 
        filter(!(abs(dg1.lfc) > dg1.lfcmin & dg1.fdr < 0.05), 
               dg2.lfc < 0-dg2.lfcmin & dg2.fdr < 0.05),
      `All Observations` = dg
    )
  }
  
  # Rename Contrasts
  names(tables)<-gsub("C1", Contrast_1, names(tables))
  names(tables)<-gsub("C2", Contrast_2, names(tables))
  if(descname){
    
    dg1.g1<-df %>% filter(!is.na(dg1.g1)) %>% pull(dg1.g1) %>% unique()
    dg1.g2<-df %>% filter(!is.na(dg1.g2)) %>% pull(dg1.g2) %>% unique()
    dg2.g1<-df %>% filter(!is.na(dg2.g1)) %>% pull(dg2.g1) %>% unique()
    dg2.g2<-df %>% filter(!is.na(dg2.g2)) %>% pull(dg2.g2) %>% unique()
    
    cols<-c(
      dg1.g1 = paste(Contrast_1, dg1.g1, sep="_"),
      dg1.g2 = paste(Contrast_1, dg1.g2, sep="_"),
      dg2.g1 = paste(Contrast_2, dg2.g1, sep="_"),
      dg2.g2 = paste(Contrast_2, dg2.g2, sep="_"),
      dg1.a1 = paste(Contrast_1, dg1.g1, dg1.unit, sep="_"),
      dg1.a2 = paste(Contrast_1, dg1.g2, dg1.unit, sep="_"),
      dg2.a1 = paste(Contrast_2, dg2.g1, dg2.unit, sep="_"),
      dg2.a2 = paste(Contrast_2, dg2.g2, dg2.unit, sep="_"),
      dg1.lfc = ifelse(
        unlog, paste(Contrast_1, "Fold_Change",sep="_"),
        paste(Contrast_1, lfc,sep="_")
      ),
      dg2.lfc=ifelse(
        unlog, paste(Contrast_2, "Fold_Change",sep="_"),
        paste(Contrast_2, lfc,sep="_")
      ),
      dg1.pvl = paste(Contrast_1, pvl,sep="_"),
      dg2.pvl = paste(Contrast_2, pvl,sep="_"),
      dg1.fdr = paste(Contrast_1, fdr,sep="_"),
      dg2.fdr = paste(Contrast_2, fdr,sep="_")
    )
    
    for (t in names(tables)){
      if(dropGroup){
        sloc<-names(cols)
        names(sloc)<-cols
        keep<-setdiff(
          names(tables[[t]]), 
          sloc[c(grep("\\.g1$", sloc), grep("\\.g2$", sloc))]
        )
        sloc<-sloc[-c(grep("\\.g1$", sloc), grep("\\.g2$", sloc)) ]
        # print(keep)
        tables[[t]]<-tables[[t]][, keep]
        tables[[t]]<-dfSubname(tables[[t]], sloc)
        
      } else {
        sloc<-names(cols)
        names(sloc)<-cols
        if(unlog){
          tables[[t]]$dg1.lfc<-ifelse(
            tables[[t]]$dg1.lfc >= 0, 2^tables[[t]]$dg1.lfc, 
            -1/2^tables[[t]]$dg1.lfc
          )
          tables[[t]]$dg2.lfc<-ifelse(
            tables[[t]]$dg2.lfc >= 0, 
            2^tables[[t]]$dg2.lfc, 
            -1/2^tables[[t]]$dg2.lfc
          )
        }
        tables[[t]]<-dfSubname(tables[[t]], sloc)
      }
    }
    
  } else {
    for(t in names(tables)){
      if(dropGroup){
        sloc<-names(cols)
        frp<-names(prf)
        names(sloc)<-cols
        names(frp)<-prf
        
        drop<-c(
          grep("\\.g1$", names(tables[[t]])), 
          grep("\\.g2$", names(tables[[t]]))
        )
        keep<-names(tables[[t]])[-drop]
        # print(keep)
        sloc<-sloc[-c(grep("g1$", sloc), grep("g2$", sloc))]
        if(unlog){
          tables[[t]]$dg1.lfc<-ifelse(
            tables[[t]]$dg1.lfc >= 0, 2^tables[[t]]$dg1.lfc, 
            -1/2^tables[[t]]$dg1.lfc
          )
          tables[[t]]$dg2.lfc<-ifelse(
            tables[[t]]$dg2.lfc >= 0, 
            2^tables[[t]]$dg2.lfc, 
            -1/2^tables[[t]]$dg2.lfc
          )
        }
        tables[[t]]<-dfSubname(tables[[t]], sloc)
        tables[[t]]<-dfSubname(tables[[t]], frp)
        
      } else {
        sloc<-names(cols)
        frp<-names(prf)
        names(sloc)<-cols
        names(frp)<-prf
        if(unlog){
          tables[[t]]$dg1.lfc<-ifelse(
            tables[[t]]$dg1.lfc >= 0, 2^tables[[t]]$dg1.lfc, 
            -1/2^tables[[t]]$dg1.lfc
          )
          tables[[t]]$dg2.lfc<-ifelse(
            tables[[t]]$dg2.lfc >= 0, 
            2^tables[[t]]$dg2.lfc, 
            -1/2^tables[[t]]$dg2.lfc
          )
        }
        tables[[t]]<-dfSubname(tables[[t]], sloc)
        tables[[t]]<-dfSubname(tables[[t]], frp)
      }
    }
  }
  tables<-append(
    tables, list(Contrasts=c(Contrast_1=Contrast_1, Contrast_2=Contrast_2))
  )
  tables
}
################################################################

# Tabulate Directional Intersections between the two data sets
# Recieves set of tables generated by "subsetTables()" returns a
# data frame with row counts for each directional intersect
tabulateOverlap<-function(tables, rename=F){
  stat.Intersect<-data.frame(
    C1_UP=c(nrow(tables[[4]]), nrow(tables[[5]])),
    C1_Down=c(nrow(tables[[2]]), nrow(tables[[3]]))
  )
  row.names(stat.Intersect)<-c("C2_Up", "C2_Down")
  if(rename){
    names(stat.Intersect)<-gsub(
      "C1", tables[['Contrasts']][1], names(stat.Intersect)
    )
    row.names(stat.Intersect)<-gsub(
      "C2", tables[['Contrasts']][2], row.names(stat.Intersect)
    )
  }
  stat.Intersect
}

# Biosig Filter -- filter a pair of deg-lists for biologically significant
# genes in both contrasts (only  applicable for a pair of RNASeq contrasts)
bioSig<-function(ar){
  ar %>%
    filter(dg1.Avg1 > 2 | dg1.Avg2 > 2) %>%
    filter(abs(dg1.Avg1 - dg1.Avg2) > 2) %>%
    filter(dg2.Avg1 > 2 | dg2.Avg2 > 2) %>%
    filter(abs(dg2.Avg1 - dg2.Avg2) > 2)
}


bioSigRNASeq<-function(
  df, lfc = "logFC", stat="FDR", a1="Avg1", a2="Avg2", 
  minExp=2, maxStat=0.05, minLfc=1
  ){
    cols<-c(lfc=lfc, stat=stat, a1=a1, a2=a2)
    df<-dfSubname(df, cols)
    df<-df %>% 
      filter(abs(lfc) >= minLfc & stat <= maxStat) %>% 
      filter(a1 >= minExp | a2 >= minExp) %>%
      filter(abs( a1 - a2) > minExp)
    
    sloc<-names(cols)
    names(sloc)<-cols
    df<-dfSubname(df, sloc)
    df
}

bioSigArray<-function(
  df, lfc = "logFC", stat="FDR", a1="Avg1", a2="Avg2", 
  minExp=0, maxStat=0.05, minLfc=1
){
  cols<-c(lfc=lfc, stat=stat, a1=a1, a2=a2)
  df<-dfSubname(df, cols)
  df<-df %>% 
    filter(abs(lfc) >= minLfc & stat <= maxStat) %>% 
    filter(a1 >= minExp | a2 >= minExp)
  sloc<-names(cols)
  names(sloc)<-cols
  df<-dfSubname(df, sloc)
  df
}


bioSigNone<-function(
  df, lfc = "logFC", stat="FDR", a1="Avg1", a2="Avg2", 
  minExp=0, maxStat=0.05, minLfc=1
){
  # Pass-through function that can handle the same parameters
  # as the other sig filters; returns the same data that it is
  # passed
  df
}



# Returns true if all values in a vector have the same sign, false otherwise
dxn<-function(x){
  x<-x[!is.na(x)]
  for(i in sign(x)){
    if(any(i != sign(x))){
      return(FALSE)
    }
  }
  return(TRUE)
}


# Given a grouped data set, returns the value of f1 that corresponds
# to the row in f2 matching the "key" value -- this function can be used
# to pivot data in a "summarize" query

pivot<-function(f1, f2, key){ 
  if(length(grep(key, f2)) == 1){
    return(
      nth(f1, grep(key, f2))
    )
  } else {
    return (NA)
  }
}

# General Purpose "Unlog" function. Given a data frame and a range or list of
# columns, convert log2 fold changes in the given columns to fold changes. 
unlog<-function(x, cols){
  for(c in cols){
    x[,c]<-ifelse(x[,c] >= 0, 2^x[,c], -1/2^x[,c])
  }
  x
}


# Function to select one row from a set of duplicates.  For any
# duplicate symbol, the record with the greatest absolute logFC is
# retained and all others are discarded
uniqueMaxLfc<-function(
  df, idc="MGI.symbol", lfc="logFC", fdr="FDR", fdr_min=0.05){
  names(df)[grep(idc, names(df))]<-"idc"
  names(df)[grep(lfc, names(df))]<-"lfc"
  names(df)[grep(fdr, names(df))]<-"fdr"
  
  df <- bind_rows(  
    df %>% 
      group_by(idc) %>% filter(min(fdr) > fdr_min) %>%
      group_by(idc) %>% filter(abs(lfc) == max(abs(lfc))) %>% 
      filter(row_number() == 1),
    df %>% 
      group_by(idc) %>% filter(min(fdr) <= fdr_min ) %>%
      group_by(idc) %>% filter(fdr <= fdr_min) %>%
      group_by(idc) %>% filter(abs(lfc) ==  max(abs(lfc)))
  )
  
  names(df)[grep("idc", names(df))]<-idc
  names(df)[grep("lfc", names(df))]<-lfc
  names(df)[grep("fdr", names(df))]<-fdr
  return (data.frame(df, stringsAsFactors=F))
}

# Helper function convert fold change to Log2 fold change
logify<-function(x, base=2){
  if(x >= 0){
    return (log(x, base))
  }
  else if( x < 0){
    return( log( 1/abs(x), base ))
  }
}

# Function to detect duplicates in a data frame and report them based
# on one or more columns given as unique ID's
printDups<-function(df, idc='MGI.symbol'){
	df$pk<-apply(data.frame(df[,idc]), 1, paste, collapse = "")
	dups<-unique(df[duplicated(df$pk), 'pk'])
	dg<-df[df$pk %in% dups,]
	dg<-dg[order(dg$pk),setdiff(names(dg), 'pk')]
	dg
}


# Function to select one row from a set of duplicates.  For any
# duplicate symbol, the record with the greatest total abundance based
# on the sum of group averages is selected for downstream analysis
uniqueTotalExp<-function(
  df, idc="MGI.symbol", av1="Avg1", av2="Avg2"){
  names(df)[grep(paste("^",idc,"$", sep=""), names(df))]<-"idc"
  names(df)[grep(paste("^",av1,"$", sep=""), names(df))]<-"av1"
  names(df)[grep(paste("^",av2,"$", sep=""), names(df))]<-"av2"
    

  df <- bind_rows(  
    df %>% 
      group_by(idc) %>% filter(n() == 1),
    df %>% 
      group_by(idc) %>% filter(n() > 1 ) %>%
      group_by(idc) %>% filter(av1 + av2 == max(av1 + av2)) %>%
      group_by(idc) %>% filter(row_number() == 1)
  )
  
  names(df)[grep("idc", names(df))]<-idc
  names(df)[grep("av1", names(df))]<-av1
  names(df)[grep("av2", names(df))]<-av2
  return (data.frame(df, stringsAsFactors=F))
}

# vennIntersections: Given a joined pair of DEG tables, this function
# tabulates the DEGs detected in both tables, or only in one table
# with total, and directional partitions. It returns a four column
# data frame. 

# vennIntersections<-function(
#   df,                       # Data frame with a joined pair of results
#   id_col="MGI.symbol",      # Unique Identifier for this gene
#   Contrast_1 = "LE",        # Name of the first contrast in df
#   Contrast_2 = "PCO",       # Name of the second contrast in df
#   dg1="dg1",                # Prefix for contrast 1
#   dg2="dg2",                # Prefix for contrast 2
#   lfc="logFC",              # Column with log 2 fold change values
#   pvl="p_value",            # Column with p value for pairwise test
#   fdr="FDR",                # Column with FDR values
#   g1 = "Group_1",           # Column with Group_1 label
#   g2 = "Group_2",           # Column with Group_2 label
#   a1 = "Avg1",              # Column with average values for Group_1
#   a2 = "Avg2",              # Column with average values for Group_2
#   stat = T,                 # Whether to use 'Stat' or 'Bio' naming scheme
#   lfcmin = 1,               # Fold Change threshold
#   minDiff = 2,              # Minimum difference for biological significance
#   minAvg  = 2               # Minimum Abundance for biologuical significance
# ){
#   # Standardize column headers
#   cols<-c(lfc=lfc, pvl=pvl, fdr=fdr, g1=g1, g2=g2, a1=a1, a2=a2)
#   df<-dfSubname(df, cols, pf="")
#   
#   prf<-c(dg1=dg1, dg2=dg2)
#   df<-dfSubname(df, prf, po="")
#   
#   df <- df %>% 
#     filter(abs(dg1.a1 - dg1.a2) > minDiff) %>%
#     filter(dg1.a1 > minAvg | dg1.a2 > minAvg) %>%
#     filter(abs(dg2.a1 - dg2.a2) > minDiff) %>%
#     filter(dg2.a1 > minAvg | dg2.a2 > minAvg)
#   
#   Venns<-data.frame(
#     Partition = c('Total', 'Increased', 'Decreased'),
#     C1 = c(
#       nrow(df %>% filter((abs(dg1.lfc) > lfcmin & dg1.fdr < 0.05) & !(abs(dg2.lfc) > lfcmin & dg2.fdr < 0.05))),
#       nrow(df %>% filter((dg1.lfc > lfcmin & dg1.fdr < 0.05) & !(abs(dg2.lfc) > lfcmin & dg2.fdr < 0.05))),
#       nrow(df %>% filter((dg1.lfc < 0-lfcmin & dg1.fdr < 0.05) & !(abs(dg2.lfc) > lfcmin & dg2.fdr < 0.05)))),
#     Both = c(
#       nrow(df %>% filter((abs(dg1.lfc) > lfcmin & dg1.fdr < 0.05) & (abs(dg2.lfc) > lfcmin & dg2.fdr < 0.05))), 
#       nrow(df %>% filter((dg1.lfc > lfcmin & dg1.fdr < 0.05) & (dg2.lfc > lfcmin & dg2.fdr < 0.05))),
#       nrow(df %>% filter((dg1.lfc < 0-lfcmin & dg1.fdr < 0.05) & (dg2.lfc < 0-lfcmin & dg2.fdr < 0.05)))),
#     C2 = c(
#       nrow(df %>% filter(!(abs(dg1.lfc) > lfcmin & dg1.fdr < 0.05) & (abs(dg2.lfc) > lfcmin & dg2.fdr < 0.05))),
#       nrow(df %>% filter(!(abs(dg1.lfc) > lfcmin & dg1.fdr < 0.05) & (dg2.lfc > lfcmin & dg2.fdr < 0.05))),
#       nrow(df %>% filter(!(abs(dg1.lfc) > lfcmin & dg1.fdr < 0.05) & (dg2.lfc < 0-lfcmin & dg2.fdr < 0.05))))
#   )
#   names(Venns)[grep('C1', names(Venns))]<-Contrast_1
#   names(Venns)[grep('C2', names(Venns))]<-Contrast_2
#   return(Venns)
# }


vennIntersections<-function(
  dg1, dg2,                     # Data frame with DEG results to compare
  dg1.idc="gene_id",            # Unique Identifier for this gene
  dg1.lbl = "LE",               # Name of the first contrast (dg1)
  dg1.lfc="logFC",              # Column with log 2 fold change values
  dg1.pvl="PValue",             # Column with p value for pairwise test
  dg1.fdr="FDR",                # Column with FDR values
  dg1.g1 = "Group_1",           # Column with Group_1 label
  dg1.g2 = "Group_2",           # Column with Group_2 label
  dg1.a1 = "Avg1",              # Column with average values for Group_1
  dg1.a2 = "Avg2",              # Column with average values for Group_2
  dg2.idc="gene_id",            # Unique Identifier for this gene
  dg2.lbl = "PCO",              # Name of the second contrast (dg2)
  dg2.lfc="logFC",              # Column with log 2 fold change values
  dg2.pvl="PValue",             # Column with p value for pairwise test
  dg2.fdr="FDR",                # Column with FDR values
  dg2.g1 = "Group_1",           # Column with Group_1 label
  dg2.g2 = "Group_2",           # Column with Group_2 label
  dg2.a1 = "Avg1",              # Column with average values for Group_1
  dg2.a2 = "Avg2",              # Column with average values for Group_2
  stat = T,                     # Whether to use 'Stat' or 'Bio' naming scheme
  dg1.lfcmin = 1,               # Fold Change threshold in dg1
  dg2.lfcmin = 1,               # Fold Change threshold in dg2
  fdrmax = 0.05,                # Maximum FDR adjusted p value
  minDiff = 2,                  # Minimum difference for biological significance
  minAvg = 2                    # Minimum Abundance for biologuical significance
){
  # Standardize table columns and headers
  dg1.cols<-c(
    idc=dg1.idc, lfc=dg1.lfc, fdr=dg1.fdr, 
    g1=dg1.g1, g2=dg1.g2, a1=dg1.a1, a2=dg1.a2
  )
  dg1<-dfSubname(dg1, dg1.cols) %>%
    dplyr::select(
      idc, lfc, fdr, 
      g1, g2, a1, a2
    ) %>%
    filter(abs(a1 - a2) > minDiff) %>%
    filter(a1 > minAvg | a2 > minAvg) %>%
    filter(abs(lfc) > dg1.lfcmin & fdr < fdrmax) 
  print(names(dg1))
  
  dg2.cols<-c(
    idc=dg2.idc, lfc=dg2.lfc, fdr=dg2.fdr, 
    g1=dg2.g1, g2=dg2.g2, a1=dg2.a1, a2=dg2.a2
  )
  dg2<-dfSubname(dg2, dg2.cols) %>%
    dplyr::select(
      idc, lfc, fdr, 
      g1, g2, a1, a2
    ) %>% 
    filter(abs(a1 - a2) > minDiff) %>%
    filter(a1 > minAvg | a2 > minAvg) %>%
    filter(abs(lfc) > dg2.lfcmin & fdr < fdrmax)
  print(names(dg2))
  
  # Tally Genes that intersect with the other data set  
  dg1<-dg1 %>%
    mutate(
      Inx.all = idc %in% (dg2 %>% pull(idc)),
      Inx.up = idc %in% (dg2 %>% filter(lfc >= 0) %>% pull(idc)),
      Inx.down = idc %in% (dg2 %>% filter(lfc < 0) %>% pull(idc))
    )
  
  dg2<-dg2 %>%
    mutate(
      Inx.all = idc %in% (dg1 %>% pull(idc)),
      Inx.up = idc %in% (dg1 %>% filter(lfc >= 0) %>% pull(idc)),
      Inx.down = idc %in% (dg1 %>% filter(lfc < 0) %>% pull(idc))
    )
  
  # Summarize Intersections
  venn<-bind_rows(
    dg1 %>%
      summarise(Partition="All", C1=sum(!Inx.all), Both=sum(Inx.all)),
    dg1 %>%
      filter(lfc > 0) %>%
      summarise(Partition="Increasing", C1=sum(!Inx.all), Both=sum(Inx.up)),
    dg1 %>%
      filter(lfc < 0) %>%
      summarise(Partition="Decreasing", C1=sum(!Inx.all), Both=sum(Inx.down))
  ) %>%
    inner_join(
      bind_rows(
        dg2 %>%
          summarise(Partition="All", C2=sum(!Inx.all)),
        dg2 %>%
          filter(lfc > 0) %>%
          summarise(Partition="Increasing", C2=sum(!Inx.all)),
        dg2 %>%
          filter(lfc < 0) %>%
          summarise(Partition="Decreasing", C2=sum(!Inx.all))
      ), by="Partition"
    )%>%
    dplyr::rename(
      !!as.name(dg1.lbl):="C1",
      !!as.name(dg2.lbl):="C2"
    )
  
    print(venn)
    # rename(
    #   C1=dg1.lbl,
    #   C2=dg2.lbl
    # )
    # 
  
  
  print(paste("Rows in dg1:", nrow(dg1), "Up:", sum(dg1$Inx.all)))
  print(paste("Rows in dg2:", nrow(dg2), "Up:", sum(dg2$Inx.all)))
  return(
    venn
  )
}

# val1<-data.frame(
#   gene_id=paste0("ENSMUS", 1:10000),           # 10000 Features
#   logFC=c(
#     rep(5, 250),                              # 250 Features lfc == 5
#     rep(1, 250),                              # 250 Features lfc == 1
#     rep(0, 9000),                             # 9000 Features lfc == 0
#     rep(-1, 250),                             # 250 Features lfc == -1
#     rep(-5, 250)                              # 250 Features lfc == -5
#   ),
#   PValue=rep(0,10000),
#   FDR = c(
#     rep(0.049, 500),                # First and last 500 rows have FDR < 0.05
#     rep(0.05, 9000),
#     rep(0.049, 500)
#   ),
#   Avg1= c(
#     rep(5, 100),                    # First and last 500 rows have FDR < 0.05
#     rep(0.5, 100),
#     rep(0.3, 50),
#     rep(10, 9500),
#     rep(1.5, 50),
#     rep(2.5, 100),
#     rep(25, 100)
#   ),
#   Avg2= c(
#     rep(25, 100),
#     rep(2.5, 100),
#     rep(1.5, 50),
#     rep(10, 9500),
#     rep(0.3, 50),
#     rep(0.5, 100),
#     rep(5, 100)
#   )
# )
# 
# # ## Duplicate "val1", but rearange feature names such that half of the significant
# ## genes in this have a change in the same direction as they do in "val1"
# val2<-val1
# val2$gene_id<-c(
#   paste("ENSMUS", seq(1,10000,2), sep=""),
#   paste("ENSMUS", seq(2,10000,2), sep="")
# )
# 
# ## Duplicate "val1", but rearange feature names such that half of the 
# ## the genes signficantly upregulated in "val1" are downregulated in "val3"
# val3<-val1
# val3$gene_id<-c(
#   paste("ENSMUS", seq(1,10000,2), sep=""),
#   paste("ENSMUS", seq(10000,2,-2), sep="")
# )
# 
# ## Duplicate "val1", but rearange feature names such that half of the 
# ## the DEG in "val4" are concordant with "val4" and half are discordant
# 
# mix<-function(a,b){
#   x<-vector()
#   for(i in 1:length(a)){
#     x<-c(x, a[i], b[i])
#   }
#   x
# }
# 
# val4<-val1
# val4$gene_id<-mix(
#   paste("ENSMUS", seq(1,10000,2), sep=""),
#   paste("ENSMUS", seq(10000,2,-2), sep="")
# )
# 
# ## Duplicate "val1", but rearange feature names such that half of the 
# ## the intersecting DEG in "val5" are concordant with "val1" and half are discordant
# ## there are 50 DEG that are only measured in val1, and 50 that are only measured in
# ## val5
# 
# mix<-function(a,b){
#   x<-vector()
#   for(i in 1:length(a)){
#     x<-c(x, a[i], b[i])
#   }
#   x
# }
# 
# val5<-val1
# val5$gene_id<-mix(
#   paste("ENSMUS", seq(51,10050,2), sep=""),
#   paste("ENSMUS", seq(10050,52,-2), sep="")
# )
# 
# degSummary(val1)
# degSummary(val2)
# degSummary(val3)
# degSummary(val4)
# degSummary(val5)
# 
# q<-query(
#   val1 %>% filter(abs(logFC) > 1 & FDR < 0.05), 
#   val2 %>% filter(abs(logFC) > 1 & FDR < 0.05),
#   id_col = "gene_id",
#   cols=setdiff(names(val1), "gene_id")
# )
# h<-compareHits(q, id_col="gene_id")
# tabulateOverlap(h)
# 
# q<-query(
#   val1 %>% filter(abs(logFC) > 1 & FDR < 0.05), 
#   val3 %>% filter(abs(logFC) > 1 & FDR < 0.05),
#   id_col = "gene_id",
#   cols=setdiff(names(val1), "gene_id")
# )
# h<-compareHits(q, id_col="gene_id")
# tabulateOverlap(h)
# 
# q<-query(
#   val1 %>% filter(abs(logFC) > 1 & FDR < 0.05), 
#   val4 %>% filter(abs(logFC) > 1 & FDR < 0.05),
#   id_col = "gene_id",
#   cols=setdiff(names(val1), "gene_id")
# )
# h<-compareHits(q, id_col="gene_id")
# tabulateOverlap(h)
# 
# q<-query(
#   val1 %>% filter(abs(logFC) > 1 & FDR < 0.05), 
#   val5 %>% filter(abs(logFC) > 1 & FDR < 0.05),
#   id_col = "gene_id",
#   cols=setdiff(names(val1), "gene_id")
# )
# h<-compareHits(q, id_col="gene_id", unlog = T)
# tabulateOverlap(h)
# 
# fj<-full_join(
#   val1 %>% filter(abs(logFC) > 1 & FDR < 0.05), 
#   val5 %>% filter(abs(logFC) > 1 & FDR < 0.05),
#   by="gene_id"
# )
# 








