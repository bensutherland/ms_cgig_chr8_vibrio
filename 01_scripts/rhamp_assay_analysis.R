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


#### 01. Read in the data, reduce columns, convert to wide format ####
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


#### 03. Join data from all plates ####
names(data.list)
head(data.list[[1]])


library("dplyr")
all_plates.df <- dplyr::bind_rows(data.list)
dim(all_plates.df)
length(unique(all_plates.df$full.id)) # 394

## All data is present in all_plates.df
## note: this assumes a constant cutoff for FP designation, not a plate-specific cutoff

head(all_plates.df)

##### 03.2 Data checking #####
# Before any correction, how much missing data is there? 
table(is.na(all_plates.df$Cq.fam) & is.na(all_plates.df$Cq.vic))
noCq_wells.df  <- all_plates.df[is.na(all_plates.df$Cq.fam) & is.na(all_plates.df$Cq.vic), ] 
write.csv(x = noCq_wells.df, file = "03_results/no_Cq_wells.csv", row.names = F)

#### 04. Correcting false-positives (both dyes) ####
# Correct for false positive detections of either dye
#   based on difference (i.e., FAM-VIC Cqs)
#   note: we are not correcting FP homozygotes yet

head(all_plates.df)
all_plates.df$Cq.vic.corr <- NA
all_plates.df$Cq.fam.corr <- NA
head(all_plates.df)

# Set FP cutoffs
FP_cutoff_fam.val <-  4
FP_cutoff_vic.val <- -4


# Loop across df to correct Cq vals that are outside FP cutoffs
for(j in 1:nrow(all_plates.df)){
  
  # If the call has both vic and fam (i.e., a 'diff' value), correct FPs as needed
  if(!is.na(all_plates.df$diff[j])){
    
    # If the call is outside of the cutoff for VIC, do not retain the value
    if(all_plates.df$diff[j] < FP_cutoff_vic.val){
      
      print(paste0(j, "- VIC false positive detected, leaving VIC value as NA"))
      
      # Retain FAM
      all_plates.df$Cq.fam.corr[j] <- all_plates.df$Cq.fam[j]
      
    # If the call is outside of the cutoff for FAM, do not retain the value
    }else if(all_plates.df$diff[j] > FP_cutoff_fam.val){
      
      print(paste0(j, "- FAM false positive detected, leaving FAM value as NA"))
      
      # Retain VIC
      all_plates.df$Cq.vic.corr[j] <- all_plates.df$Cq.vic[j]
    
    # If the call is inside both cutoffs, retain the value
    }else{
      
      # Retain VIC
      all_plates.df$Cq.vic.corr[j] <- all_plates.df$Cq.vic[j]
      
      # Retain FAM
      all_plates.df$Cq.fam.corr[j] <- all_plates.df$Cq.fam[j]
    }
      
  # If the diff is NA, we can't remove potential FPs, so keep both values as is
    # TODO: set hard cutoff as well (?)
    
  }else{
    
    all_plates.df$Cq.vic.corr[j] <- all_plates.df$Cq.vic[j]
    all_plates.df$Cq.fam.corr[j] <- all_plates.df$Cq.fam[j]
    
  }
  
}


head(all_plates.df, n = 20)


#### 05. Calling genotypes ####

# Loop over the dataframe, get the difference between Cqs and the derived genotype
for(j in 1:nrow(all_plates.df)){
  
  # If both Cq are NA, it is an uncalled sample
  if(is.na(all_plates.df$Cq.fam[j]) & is.na(all_plates.df$Cq.vic.corr[j])){
    
    all_plates.df$geno[j] <- "no.geno"
    
  # If VIC is present & FAM is NA, it is a homozygous alternate
  }else if(is.na(all_plates.df$Cq.fam[j]) & !is.na(all_plates.df$Cq.vic.corr[j])){
    
    all_plates.df$geno[j] <- "homozyg.alt"
    
  # If FAM is present & VIC is NA, it is a homozygous reference
  }else if(!is.na(all_plates.df$Cq.fam[j]) & is.na(all_plates.df$Cq.vic.corr[j])){
    
    all_plates.df$geno[j] <- "homozyg.ref"
    
  # If both dyes are present, after correction, it is a true heterozygote
  }else if(!is.na(all_plates.df$Cq.fam[j]) & !is.na(all_plates.df$Cq.vic.corr[j])){
    
    all_plates.df$geno[j] <- "heterozyg"
    
  }
  
}

table(all_plates.df$geno)

head(all_plates.df)


