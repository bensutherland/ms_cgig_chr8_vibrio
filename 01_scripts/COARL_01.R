# Import VCF from wgrs_workflow and perform general population statistics
# B. Sutherland (2024-02-09)


### Front Matter ####
# Clean space
# rm(list=ls())

# Prior to running the following, source simple_pop_stats and choose Pacific oyster

## Install and load packages
#install.packages("rstudioapi")
#install.packages("adegenet")
#install.packages("vcfR")

library("rstudioapi")
library("adegenet")
library("vcfR")

## Info
# sessionInfo()

# Set variables
input.FN <- "../ms_cgig_chr8/02_input_data/mpileup_calls_filt_AF_0.05_LD.0.5.50kb_subset_0.01.vcf"


#### 01. Load data ####
my_data.vcf <- read.vcfR(file = input.FN)

# convert vcf to genlight
my_data.gl <- vcfR2genlight(x = my_data.vcf)
my_data.gl

# convert vcf to genind
my_data.gid <- vcfR2genind(x = my_data.vcf)
my_data.gid


#### 02. Basic characteristics ####
indNames(my_data.gid) # use this to link to the breeding map
pop(my_data.gid) <- rep(x = "OFR", times = nInd(my_data.gid))

#### 03. Basic analyses ####
pca_from_genind(data = my_data.gid, PCs_ret = 4, plot_eigen = T, retain_pca_obj = T
                , plot_allele_loadings = F
                )





