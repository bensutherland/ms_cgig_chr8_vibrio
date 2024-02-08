# Calling genotypes
# B. Sutherland, Liam Surry, VIU (initialized 2024-01-11)

# Note: requires that 01_scripts/rhamp_03_call_genos.R has been run

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

## Load libraries

## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path) # take main directory
setwd(current.path)
rm(current.path)

# Load data
load(file = "03_results/prepared_rhamp_results.RData")

# Set user variables
sample_info.FN <- "00_archive/OSU_CHR8_VC_Mapping_Family114-117_2024-01-21_MFrhAmpDNAid.txt"

### Load sample info data
annot.df <- read.delim(file = sample_info.FN, header = T, sep = "\t")
head(annot.df)


### Format full geno data
head(all_plates.df)

# Reduce the cols
data.df  <-    all_plates.df[, 
                             c("Well", "Cq.fam", "Cq.fam.corr", "Cq.vic", "Cq.vic.corr", "diff", "geno", "full.id")]

data.df <- as.data.frame(data.df)
head(data.df)


#### Connect the full.id in the data to the full.id in the annotation file
data_annot.df <- merge(x = data.df, y = annot.df, by = "full.id", all.x = T)
head(data_annot.df)


write.table(x = data_annot.df, file = "03_results/rhAmp_full_results.txt", quote = F, sep = "\t", row.names = F)



# Next: go to 01_scripts/rhamp_04_well_to_sample.R

# # Separate full ID into components
# data.df    <-  separate(data = data.df, col = "full.id"
#                      , into = c("plate.id", "well.id")
#                      , sep = "__", remove = T)
# 
# head(data.df)



