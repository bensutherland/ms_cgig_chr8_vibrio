# Import VCF from wgrs_workflow and perform general population statistics
# B. Sutherland (2024-02-09)

### Front Matter ####
# Clean space
# rm(list=ls())

# Prior to running the following, source simple_pop_stats and choose Pacific oyster

## Install and load packages (not included in sps)
#install.packages("vcfR")

library("vcfR")

## Info
# sessionInfo()

# Set variables
input.FN <- "../ms_cgig_chr8/02_input_data/mpileup_calls_filt_AF_0.05_LD.0.5.50kb_subset_0.01.vcf"
indiv.FN <- "../ms_cgig_chr8/00_archive/COARL2_parental_genotyping_label_map.txt"


#### 01. Load data ####
# read in vcf
my_data.vcf <- read.vcfR(file = input.FN)

# convert vcf to genlight
my_data.gl <- vcfR2genlight(x = my_data.vcf)
my_data.gl

# convert vcf to genind
my_data.gid <- vcfR2genind(x = my_data.vcf)
my_data.gid



#### 02. Preparing data ####
# Rename samples to sample identifier
inds.df <- indNames(my_data.gid) # what are the indiv names and order in the data? 
inds.df <- as.data.frame(inds.df)
head(inds.df)

# Separate into component parts of the name
inds.df <- separate(data = inds.df, col = "inds.df", into = c("project", "sample.num", "file.info")
                     , sep = "-", remove = T
                   )
head(inds.df)

# Create name that will match the sample identifier df
inds.df <- paste0(inds.df$project, "_", inds.df$sample.num)
inds.df <- as.data.frame(inds.df)
head(inds.df) 

# Data checking
cbind(indNames(my_data.gid), inds.df) # this remains in the order of the data

# Read in conversion file (tube name to sample name)
conversion.df <- read.table(file = indiv.FN, header = T, sep = "\t")
head(conversion.df)

# Convert from tube_label to sample_ID
head(inds.df)       # this is the ordered identifiers from the data
dim(inds.df)
head(conversion.df) # this is the correspondence to the true sample IDs
dim(conversion.df)

# Remove other project info from the conversion file (note: not clear why merge doesn't drop them)
conversion.df <- conversion.df[grep(pattern = "COARL", x = conversion.df$tube_label), ]

ordered_conversion.df <- merge(x = inds.df, y = conversion.df
                          , by.x = "inds.df", by.y = "tube_label", sort = F
                          #, all.x = T, all.y = F
                          ) # important to not sort

head(conversion.df) # ordered identifiers 
dim(conversion.df)

indNames(my_data.gid) <- conversion.df$sample_ID
#pop(my_data.gid) <- indNames(my_data.gid)

pop(my_data.gid) <- gsub(pattern = "_.*", replacement = "", x = indNames(my_data.gid))

# pop(my_data.gid) <- rep(x = "OFR", times = nInd(my_data.gid))


#### 03. Basic analyses ####
pca_from_genind(data = my_data.gid, PCs_ret = 4, plot_eigen = T, retain_pca_obj = T
                , plot_allele_loadings = F
                )

# Run relatedness analysis

# Plot MAF distribution (note: LD-filtered)


#### 04. Bring in cross information, per cross, per locus, infer offspring AF ####
# 


