# Determine phase shifts
# 2024-07-15
# B. Sutherland

library(data.table)
library(vcfR)

setwd("~/Documents/cgig/CHR8_impute/ms_cgig_chr8")

# Get VCF file to have locus names
input.FN <- "04_impute_panel/wgrs_filtered_parent_loci_amp_panel_parent_loci_amp_panel_offspring_loci_NC_047567.1.vcf"
my_vcf <- read.vcfR(file = input.FN)
geno_info.df <- my_vcf@fix
geno_info.df <- as.data.frame(geno_info.df)
mname <- paste0(geno_info.df$CHROM, "__", geno_info.df$POS, "__", geno_info.df$ID)
length(mname)
head(mname)

which(!is.na(geno_info.df$ID)==TRUE)


geno  <- fread(file = "NC047567.genotypes", sep = " ")
geno[1:40,1:10]

mom.geno    <- geno[geno$V1=="65-4F", ]
dad.geno    <- geno[geno$V1=="58-33M", ]
offspr.geno <- geno[grep(pattern = "_114_", x = geno$V1), ] 

phase <- fread(file = "NC047567.haplotypes", sep = " ")
phase[1:40,1:10]

mom.phase    <- phase[phase$V1=="65-4F", ]
dad.phase    <- phase[phase$V1=="58-33M", ]
offspr.phase <- phase[grep(pattern = "_114_", x = phase$V1), ] 

mom.geno[ ,c(1,47290:47310)]
mom.phase[,c(1,47290:47310)]

dad.geno[ ,c(1,47290:47310)]
dad.phase[,c(1,47290:47310)]

offspr.geno[1:4,c(1,47290:47310)]
offspr.phase[1:8,c(1,47290:47310)]


subset.df <- rbind(mom.geno, mom.phase, dad.geno, dad.phase, offspr.geno, offspr.phase)
colnames(subset.df) <- c("ind", mname)

fwrite(x = subset.df, file = "subset_data.txt", sep = "\t", quote = F)





