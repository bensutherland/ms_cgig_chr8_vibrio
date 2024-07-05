# Create working VCF in correct orientation
# B. Sutherland
# 2024-07-05

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

# Load libraries
#install.packages("tidyr")
require(tidyr)

## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path) # take main directory
setwd(current.path)
rm(current.path)

# Set variable names
source_map.FN <- "02_input_data/G0923-21-VIUN.map"
correct_bedfile.FN <- "../amplitargets/cgig/current/WGAG22008_BJS_OYRv01_Hotspot_A.bed"

# Read in plink map file
map.df <- read.table(file = source_map.FN, header = F, sep = "\t")
colnames(map.df) <- c("chr", "mname", "allele1", "allele2")
head(map.df)

# Read in hotspot bed file
bed.df <- read.table(file = correct_bedfile.FN, header = F, sep = "\t", skip = 1)
head(bed.df)
dim(bed.df)
colnames(bed.df) <- c("contig", "start", "stop", "mname", "spacer", "strand", "details", "type")
bed.df <- as.data.frame(bed.df)
head(bed.df)

# Separate details filed to get ref and obs
bed.df <- separate(data = bed.df, col = "details", into = c("ref", "alt", "anchor"), sep = ";", remove = F)
head(bed.df)
bed.df <- separate(data = bed.df, col = "ref", into = c("ref", "ref.allele"), sep = "=", remove = T)
bed.df <- separate(data = bed.df, col = "alt", into = c("alt", "alt.allele"), sep = "=", remove = T)
head(bed.df)

# Combine the map file and hotspot file, keeping in the order of the map
all.df <- merge(x = map.df, y = bed.df, by = "mname", all.x = T, sort = F)
head(all.df)
dim(all.df)
unique(all.df$ref.allele)
unique(all.df$alt.allele)

# # Fix map
# #fixed_map.df <- cbind(paste0(all.df$contig, "__", all.df$start), all.df$mname, all.df$ref.allele, all.df$alt.allele)
# #fixed_map.df <- cbind(all.df$contig, all.df$mname, all.df$ref.allele, all.df$alt.allele)
# #fixed_map.df <- cbind(rep(0, times = nrow(all.df)), all.df$mname, all.df$ref.allele, all.df$alt.allele)
# fixed_map.df <- cbind(all.df$mname, all.df$start, all.df$ref.allele, all.df$alt.allele)
# fixed_map.df <- as.data.frame(fixed_map.df)
# head(fixed_map.df)
# 
# # Write out
# write.table(x = fixed_map.df, file = "G0923-21-VIUN_fixed_alleles.map", quote = F, sep = "\t"
#             , col.names = F, row.names = F
#             )

# Create plink file to set the reference allele
my_ref_allele.df <- cbind(all.df$mname, all.df$ref.allele)
my_ref_allele.df <- as.data.frame(my_ref_allele.df)
head(my_ref_allele.df)

# Write out
write.table(x = my_ref_allele.df, file = "02_input_data/G0923-21-VIUN_set_allele.txt", quote = F, sep = "\t"
            , col.names = F, row.names = F
)

# The written file can be used to reorient the plink map and ped to VCF file with the correct reference allele

# End
