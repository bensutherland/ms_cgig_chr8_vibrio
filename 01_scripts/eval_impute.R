# Read in AlphaImpute output file, evaluate imputation
# B. Sutherland (2024-07-15)

### Front Matter ####
# Clean space
# rm(list=ls())

## Install and load packages
#install.packages("rstudioapi")
#install.packages("vcfR")
#install.packages("data.table")
#install.packages("dartR")
library("rstudioapi")
library("vcfR")
library("adegenet")
library("data.table")
library("dartR")

# Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

## Info
# sessionInfo()

# Set variables
input.FN     <- "12_impute_impute/genos_imputed_converted.txt"
input_vcf.FN <- "../../CHR8_wgrs/wgrs_workflow/05_genotyping/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1._NC047567.1.vcf" 

# Read in the imputed data
imputed.df <- fread(file = "12_impute_impute/genos_imputed_converted.txt", sep = "\t")
dim(imputed.df) # 775,020
imputed.df[1:5,1:5]

# Read in VCF file
my_vcf <- read.vcfR(file = input_vcf.FN) # n = 230,329 variants
my_vcf

genos_10X.df <- extract.gt(x = my_vcf, element = "GT")
dim(genos_10X.df)
genos_10X.df[1:5,1:5]

# Now need to make rows and cols match for the following
imputed.df[1:5,1:5]
library(tidyr)
imputed.df <- separate(data = imputed.df, col = "mnames.df", into = c("chr", "pos", "marker.id"), sep = "__")
imputed.df[1:5,1:5]
mnames.df <- paste0(imputed.df$chr, "_", imputed.df$pos)
mnames.df <- as.data.frame(mnames.df)
colnames(mnames.df) <- "mname"
imputed_w_names.df <- cbind(mnames.df, imputed.df)
imputed_w_names.df[1:5,1:20]
imputed_w_names.df <- imputed_w_names.df[grep(pattern = "mname|ASY2", x = colnames(imputed_w_names.df))]
imputed_w_names.df[1:5,1:5]


genos_10X.df[1:5,1:5]
genos_10X.df <- as.data.frame(genos_10X.df)
mnames.df <- rownames(genos_10X.df)
mnames.df <- as.data.frame(mnames.df)
head(mnames.df)
colnames(mnames.df) <- "mname"

genos_10X_w_names.df <- cbind(mnames.df, genos_10X.df)
genos_10X_w_names.df[1:5,1:5]

# NEXT: make colnames match, then convert to 0, 1, 2 format, then compare
