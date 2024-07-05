# Update the provided VCF with the contig IDs and coordinates provided with the panel
#  requires the input VCF and the contig/position information from Sutherland et al. 2024 (G3), Additional File S1 (see README)
#  initialized 2024-06-13
#  Ben J. G. Sutherland (VIU)

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

# Load libraries
library(vcfR)

## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path) # take main directory
setwd(current.path)
rm(current.path)

# User set variables
chr_info.FN <- "00_archive/additional_file_S1_amp_panel_design_info.txt"
vcf.FN <- "02_input_data/G0923-21-VIUN_converted.vcf.gz"

# Load chr info and prepare marker ID, contig name, and position of SNP
chr_info.df <- read.delim(file = chr_info.FN)
head(chr_info.df, n = 2)
chr_info.df$marker_ID <- as.character(chr_info.df$marker_ID) # ensure marker ID is a character
str(chr_info.df)
chr_info.df <- chr_info.df[,c("marker_ID", "chr", "start")] # retain only the needed cols
head(chr_info.df)

# Load input VCF
my_vcf <- read.vcfR(file = vcf.FN)

# # Inspect
# test.df <- my_vcf@fix
# test.df <- as.data.frame(test.df)
# sort(table(test.df$ID, useNA = "ifany"), decreasing = T)
# head(test.df)
# rm(test.df)
# 
# # Remove any locus that does not have a marker name
# my_vcf <- my_vcf[!is.na(my_vcf@fix[,"ID"]), ]
# my_vcf@gt[1:10,1:10]
# geno.df <- my_vcf@gt
# View(geno.df)


#### 02. Replace missing info with chr info ####
# Pull out the fix section of the VCF (to replace the 0 and 0 for CHROM and POS)
fix_section <- my_vcf@fix
fix_section <- as.data.frame(fix_section)
head(fix_section)

# Combine the chr info with the VCF fix section
replacement_fix <- merge(x = fix_section, y = chr_info.df, by.x = "ID", by.y = "marker_ID", sort = F, all.x = T)
dim(replacement_fix)
head(replacement_fix)

# Reorder cols
replacement_fix <- replacement_fix[,c("chr", "start", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")]
head(replacement_fix)

# Rename as per the fix section
colnames(replacement_fix) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
head(replacement_fix)

# Make matrix
replacement_fix.mat <- as.matrix(replacement_fix)

# Replace existing fix section with the complete information
my_vcf@fix <- replacement_fix.mat
my_vcf@fix[1:5,]

# Write out compressed VCF file
output.FN <- paste0("03_results/", gsub(pattern = "\\.vcf\\.gz", x =  basename(vcf.FN), replacement = "_annot.vcf.gz"))
print(paste0("Your output will be a compressed VCF with the name ", output.FN))
vcfR::write.vcf(x = my_vcf, file = output.FN)

# Now run snplift to transfer coordinates to chr assembly genome
