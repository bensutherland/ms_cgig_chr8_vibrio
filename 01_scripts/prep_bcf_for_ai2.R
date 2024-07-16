# Convert BCF file to input needed for AlphaImpute2
# B. Sutherland (2024-07-12)

### Front Matter ####
# Clean space
# rm(list=ls())

## Install and load packages
#install.packages("rstudioapi")
#install.packages("vcfR")
#install.packages("data.table")
library("rstudioapi")
library("vcfR")
library("adegenet")
library("data.table")

# Set working directory to the ms_scallop_popgen repo
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

## Info
# sessionInfo()

# Set variables
#input.FN <- "12_impute_impute/all_inds_wgrs_and_panel_NC_047567_1.vcf" # subset dataset
input.FN <- "12_impute_impute/all_inds_wgrs_and_panel.vcf"


#### 01. Import data ####
input.vcf <- vcfR::read.vcfR(file = input.FN)
input.vcf


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

# Temporary fix, drop triallelic site
genos.df <- gsub(pattern = "1/2", replacement = "9", x = genos.df)
genos.df <- gsub(pattern = "0/2", replacement = "9", x = genos.df)
genos.df <- gsub(pattern = "2/2", replacement = "9", x = genos.df)

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
