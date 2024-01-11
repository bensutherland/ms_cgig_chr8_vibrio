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

filename.FN <- filenames.vec[1]


data.df <- read.csv(file = filename.FN, header = T)
dim(data.df)
head(data.df)
str(data.df)

data.df <- data.df[, c("Well", "Fluor", "Content", "Cq")]
head(data.df)

unique(data.df$Fluor)

# Go from long form to horizontal form
data_FAM.df <- data.df[data.df$Fluor=="FAM",]
data_VIC.df <- data.df[data.df$Fluor=="VIC",]
dim(data.df)
dim(data_FAM.df)
dim(data_VIC.df)
head(data_VIC.df)
head(data_FAM.df)

data_wide.df <- merge(x = data_FAM.df, y = data_VIC.df, by = "Well"
                      , suffixes = c(".fam", ".vic")
                      )
head(data_wide.df)
data_wide.df$geno <- NA
head(data_wide.df)

data_wide.df$diff <- data_wide.df$Cq.fam - data_wide.df$Cq.vic
head(data_wide.df)

filename.FN
short_filename.FN <- gsub(pattern = "_2023.*", replacement = "", x = filename.FN)
short_filename.FN <- gsub(pattern = "02_input_data//", replacement = "", x =  short_filename.FN)

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

