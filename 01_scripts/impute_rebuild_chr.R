# Rebuild the full datafile after ai2 that has been sep by chr
# B. Sutherland (2024-07-22)

### Front Matter ####
# Clean space
# rm(list=ls())

## Install and load packages
#install.packages("rstudioapi")
#install.packages("data.table")
library("rstudioapi")
library("data.table")
library("tidyr")

# Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

## Info
# sessionInfo()

# Set user variables
input_files.vec <- list.files(path = "12_impute_impute/", pattern = ".genotypes")

#### Read in imputed chr data ####
imputed_chr.list <- list(); input.FN <- NULL
for(i in 1:length(input_files.vec)){
  
  # Set filename
  input.FN <- paste0("12_impute_impute/", input_files.vec[i])
  
  # Reporting
  print(paste0("Reading in file ", input.FN))
  
  # Read in
  imputed_chr.list[[input_files.vec[i]]] <- fread(file = input.FN)
  
}


#### Read in pre-imputed chr data for the marker names ####

# Also bring in mnames
# load(file = "12_impute_impute/prepared_matrices_for_imputation.RData")
# input.dat$mname


#### THE FOLLOWING TO BE DELETED WHEN THE ABOVE WAS RUN CORRECTLY ####
# Set variables
input.FN             <- "12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic_ai2.txt"
chr_indicator.string <- "NC_047"

#### 01. Import data ####
input.dat <- fread(file = input.FN, sep = "\t", header = T)
input.dat[1:5,1:5]
dim(input.dat)

#### /END/ THE FOLLOWING TO BE DELETED WHEN THE ABOVE WAS RUN CORRECTLY ####

# Keep the marker names only for those chrs that were collected
input.dat <- input.dat[grep(pattern = chr_indicator.string, x = input.dat$mname), ]
input.dat[1:5,1:5]
input.dat$mname


#### Apply the mname to each slot of the imputed list ####

names(imputed_chr.list)

chr.oi <- NULL; mname.vec <- NULL; imputed_named.df <- NULL
for(i in 1:length(imputed_chr.list)){
  
  # Determine the chr of interest
  chr.oi <- names(imputed_chr.list)[i]
  chr.oi <- gsub(pattern = "ai2_input_", replacement = "", x = chr.oi)
  chr.oi <- gsub(pattern = ".genotypes", replacement = "", x = chr.oi)
  print(paste0("Working on ", chr.oi))
  
  # obtain mnames for this chr
  mname.vec <- input.dat$mname[grep(pattern = chr.oi, x = input.dat$mname)]
  head(mname.vec)
  mname.vec <- c("ind", mname.vec)
  head(mname.vec)
  length(mname.vec)
  
  imputed_chr.list[[i]][1:5,1:5]
  print(ncol(imputed_chr.list[[i]])==length(mname.vec))
  
  colnames(imputed_chr.list[[i]]) <- mname.vec
  
  # Extract and build
  imputed_named.df <- cbind(imputed_named.df, imputed_chr.list[[i]])
  imputed_named.df[1:5,1:5]
  print(dim(imputed_named.df))
  
}

#imputed_named.df[1:5,1:5]

# Then write out
### NOTE: this will not work, because the space is separating NC* and position ###
fwrite(x = imputed_named.df, file = "13_impute_compare/imputed_named_not_transposed.txt", sep = " "
       , quote = F, col.names = T, row.names = T
)

all_chr.df <- fread(input = "13_impute_compare/imputed_named_not_transposed.txt", sep = " ")


# Then transpose
imputed_named_t.df <- as.data.frame(imputed_named.df)
imputed_named_t.df <- t(imputed_named_t.df)

fwrite(x = imputed_named_t.df, file = "13_impute_compare/imputed_named_t.txt", sep = " "
       , quote = F, col.names = T, row.names = T
)
