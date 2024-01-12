# Take raw data from rhAmp assay run on a CFX96 instrument
# Input data should be csv format
# B. Sutherland, Liam Surry (2024-01-11)

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path) # take main directory
setwd(current.path)
rm(current.path)

# Identify the files
filenames.vec <- list.files(path = "02_input_data/", pattern = ".csv", full.names = T)

# Isolate only the column-purified plates
filenames.vec <- head(sort(filenames.vec), n = 5) # The first five plates are column purified

#### 01. Reading in and viewing the data ####
filename.FN <- NULL; short_filename.FN <- NULL; data.df <- NULL
data_FAM.df <- NULL; data_VIC.df <- NULL; data_wide.df <- NULL
data.list <- list()
for(i in 1:length(filenames.vec)){
  
  # Set the filename for this loop
  filename.FN <- filenames.vec[i]
  
  # Make a short filename for use in saving out
  filename.FN
  short_filename.FN <- gsub(pattern = "_202.*", replacement = "", x = filename.FN)
  short_filename.FN <- gsub(pattern = "02_input_data//", replacement = "", x =  short_filename.FN)
  
  # Reporting
  print(paste0("Working on plate ", short_filename.FN))
  
  # Read in the datafile
  data.df <- read.csv(file = filename.FN, header = T)
  dim(data.df)
  head(data.df)
  #str(data.df)
  
  # Subset to only the required columns
  data.df <- data.df[, c("Well", "Fluor", "Content", "Cq")]
  head(data.df)
  
  # What dyes are present? 
  unique(data.df$Fluor)
  
  # Go from long form to horizontal form by separating VIC and FAM into individual df
  data_FAM.df <- data.df[data.df$Fluor=="FAM",]
  data_VIC.df <- data.df[data.df$Fluor=="VIC",]
  dim(data.df)
  dim(data_FAM.df)
  dim(data_VIC.df)
  head(data_VIC.df)
  head(data_FAM.df)
  
  # Combine VIC and FAM back into a single, wide df
  data_wide.df <- merge(x = data_FAM.df, y = data_VIC.df, by = "Well"
                        , suffixes = c(".fam", ".vic")
  )
  head(data_wide.df)
  data_wide.df$geno <- NA
  head(data_wide.df)
  
  #### 02. Inspecting difference between the dyes ####
  # Calculate the difference between VIC and FAM detections
  data_wide.df$diff <- data_wide.df$Cq.fam - data_wide.df$Cq.vic
  head(data_wide.df)
  # note: when one of the dyes is missing (due to genotype, or to a lack of detection)
  #    , then the diff will be NA
  
  # Plot the difference between the dyes, as well as scatterplot both dyes together
  #   note: only considers when both dyes are detected (no homozygotes, no missing data)
  pdf(file = paste0("03_results/Cq_diffs_", short_filename.FN, ".pdf" ), width = 10.5, height = 5)
  par(mfrow=c(1,2))
  plot(data_wide.df$diff, ylab = "Cq difference (FAM - VIC)"
       , main = short_filename.FN
       , las = 1
  )
  plot(data_wide.df$Cq.fam, data_wide.df$Cq.vic
       , ylab = "Cq. value (VIC)"
       , xlab = "Cq. value (FAM)"
       , main = short_filename.FN
       , las = 1
       
  )
  dev.off()
  
  ## Save out the dataframe to a list
  data.list[[short_filename.FN]] <- data_wide.df
  
  
}

names(data.list)
data_wide.df <- data.list[["OCV23_rhAmp_plate_02"]]
head(data_wide.df)

# Plate 2 has both false positive VIC and false positive FAM

#### 03. Correcting false-positives (VIC) ####
# Correct for false positive detections of VIC
head(data_wide.df)
data_wide.df$Cq.vic.corr <- NA
str(data_wide.df$Cq.vic.corr)

for(j in 1:nrow(data_wide.df)){
  
  if(!is.na(data_wide.df$diff[j])){
    
    if(data_wide.df$diff[j] < -5){
      
      print("False positive detected, leaving value as NA")
      
    }else{
      
      data_wide.df$Cq.vic.corr[j] <- data_wide.df$Cq.vic[j]
      
    }
    
  }else{
    
    data_wide.df$Cq.vic.corr[j] <- data_wide.df$Cq.vic[j]
    
  }
  
}


head(data_wide.df)

# Loop over the dataframe, get the difference between Cqs and the derived genotype
for(j in 1:nrow(data_wide.df)){
  
  if(is.na(data_wide.df$Cq.fam[j]) & is.na(data_wide.df$Cq.vic.corr[j])){
    
    data_wide.df$geno[j] <- "no.geno"
    
  }else if(is.na(data_wide.df$Cq.fam[j]) & !is.na(data_wide.df$Cq.vic.corr[j])){
    
    data_wide.df$geno[j] <- "homozyg.alt"
    
  }else if(!is.na(data_wide.df$Cq.fam[j]) & is.na(data_wide.df$Cq.vic.corr[j])){
    
    data_wide.df$geno[j] <- "homozyg.ref"
    
  }else if(!is.na(data_wide.df$Cq.fam[j]) & !is.na(data_wide.df$Cq.vic.corr[j])){
    
    data_wide.df$geno[j] <- "heterozyg"
    
  }
  
}

table(data_wide.df$geno)

