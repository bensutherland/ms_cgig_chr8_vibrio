# Convert geno matrix to inputs needed for AlphaImpute2
# requires that 01_scripts/bcf_to_ai2.sh has already been run
# B. Sutherland (2024-07-19)

### Front Matter ####
# Clean space
# rm(list=ls())

## Install and load packages
#install.packages("rstudioapi")
#install.packages("data.table")
library("rstudioapi")
library("data.table")

# Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

## Info
# sessionInfo()

# Set variables
input.FN <- "12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic_ai2_full.txt"

#### 01. Import data ####
input.dat <- fread(file = input.FN, sep = "\t", header = T)
input.dat[1:5,1:5]

# Find the unique chr present
chr.df <- input.dat$mname
chr.df <- as.data.frame(chr.df)
head(chr.df)
library(tidyr)
chr.df <- separate(data = chr.df, col = "chr.df", into = c("chr", "pos"), sep = " ", remove = F)
head(chr.df)

chr_present <- unique(chr.df$chr)
chr_present <- chr_present[grep(pattern = "NC_047", x = chr_present)]
chr_present 

matrix.list <- list()
for(i in 1:length(chr_present)){
  
  print(paste0("Working on chr ", chr_present[i]))
  
  matrix.list[[chr_present[i]]] <- input.dat[grep(pattern = chr_present[i], x = input.dat$mname), ]
  
  print("Number records for this chr: ")
  print(nrow(matrix.list[[chr_present[i]]]))
  
}

length(matrix.list)

# Next:
# - ensure it is in the correct format, need to be transposed? 
# - write out each element of the matrix
# - build bash script to run ai2 per chr


## R solution
# header.df <- read.delim(file = "12_impute_impute/temp_ind_names.txt", header = F, sep = "\t")
# header.df <- as.data.frame(header.df)
# header.df <- c("mname", header.df$V1)
# header.df <- as.data.frame(header.df)
# header.df <- t(header.df)
# input.dat <- rbind(header.df, input.dat)


#### 02. Prepare genotypes file ####
genos.df <- extract.gt(input.vcf)
dim(genos.df)
genos.df[1:5,1:5]

# Convert all missing data (NA) to "9
genos.df[is.na(genos.df)] <- "9"
genos.df[1:10,1:10]

# Convert to alt allele dosage
genos.df <- gsub(pattern = "0/0", replacement = "0", x = genos.df)
genos.df <- gsub(pattern = "0/1", replacement = "1", x = genos.df)
genos.df <- gsub(pattern = "1/1", replacement = "2", x = genos.df)

# Convert to df
genos.df <- as.data.frame(genos.df)
genos.df[1:10, 1:10]

# Transpose
genos_t.df <- t(genos.df)
genos_t.df <- as.data.frame(genos_t.df)
genos_t.df[1:10,1:4]

# Make sure all loci for an individual are what are expected
#unique(genos_t.df[which(rownames(genos_t.df)=="55-41F"),])

# Create a vector of individual names
indiv_names.df <- rownames(genos_t.df)
indiv_names.df <- as.data.frame(indiv_names.df)
colnames(indiv_names.df) <- "ind"
head(indiv_names.df)

# Combine ind names and genotypes
genos_t.df <- cbind(indiv_names.df, genos_t.df)
genos_t.df[1:10,1:4]

# Write output
data.table::fwrite(x = genos_t.df, file = "12_impute_impute/genos.txt", quote = F, sep = " ", row.names = F, col.names = F)


#### 03. Prepare pedigree file ####
pedigree.df <- rownames(genos_t.df)
pedigree.df <- as.data.frame(pedigree.df)
pedigree.df$sire <- NA
pedigree.df$dame <- NA

head(pedigree.df)

write.table(x = pedigree.df, file = "12_impute_impute/pedigree.txt", sep = " ", quote = F
            , row.names = F, col.names = F
)

# Next, annotate the pedigree then go to AlphaImpute2 for imputation
