# Data input and formatting of raw rhAmp assay data (csv format)
# B. Sutherland, Liam Surry, VIU (initialized 2024-01-11)

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path) # take main directory
setwd(current.path)

# Set user variables
sample_info.FN <- "00_archive/CHRM8_AAB24_Liam_Surry_2024-06-07.csv"

## Load libraries
#install.packages("dplyr")
library("dplyr")

# Identify the files
filenames.vec <- list.files(path = "02_input_data/", pattern = ".csv", full.names = T)


#### Read in the data, reduce columns, convert to wide format ####
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
  #dim(data.df)
  #head(data.df)
  #str(data.df)
  
  print(paste0("This plate has ", length(unique(data.df$Well)), " wells."))
  
  # Subset to only the required columns
  data.df <- data.df[, c("Well", "Fluor", "Content", "Cq")]
  #head(data.df)
  
  # Update NaN to NA
  data.df$Cq[data.df$Cq=="NaN"] <- NA
  
  # What dyes are present? 
  print(paste0("Fluorophores present: "))
  print(unique(data.df$Fluor))
  
  # Convert from long format to wide format
  data_FAM.df <- data.df[data.df$Fluor=="FAM",]
  data_VIC.df <- data.df[data.df$Fluor=="VIC",]

  # Combine VIC and FAM back into a single, wide df
  data_wide.df <- merge(x = data_FAM.df, y = data_VIC.df, by = "Well"
                        , suffixes = c(".fam", ".vic")
                        )
  
  # Add empty genotype column
  data_wide.df$geno <- NA
  
  # Add full identifier
  data_wide.df$full.id <- NA
  data_wide.df$full.id <- paste0(short_filename.FN, "__", data_wide.df$Well)
  
  
  #### 02. Inspecting difference between the dyes ####
  # Calculate FAM - VIC
  data_wide.df$diff <- data_wide.df$Cq.fam - data_wide.df$Cq.vic
  head(data_wide.df)
  # note: when one of the dyes is missing (due to genotype, or to a lack of detection)
  #    , then the diff column will be an NA
  
  # When both dyes are present, plot FAM-VIC
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
       #, main = short_filename.FN
       , las = 1
       
  )
  dev.off()
  
  ## Save out the dataframe to a list
  data.list[[short_filename.FN]] <- data_wide.df
  
  
}

names(data.list)

head(data.list[[1]])


## The data is now prepared for analysis, go to 01_scripts/rhamp_02_data_checking.R

