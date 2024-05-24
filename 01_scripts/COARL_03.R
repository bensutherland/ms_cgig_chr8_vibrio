# Prepare data for GWAS
# B. Sutherland (initialized 2024-02-14)

# Requires that 01_scripts/COARL_02.R was already run
# Clear workspace, launch simple_pop_stats

# Install/ load additional libraries
#install.packages("vcfR")
#devtools::install_github('kaustubhad/fastman', build_vignettes = TRUE)
#install.packages("missMethods")
library(vcfR)
library(missMethods)
library(fastman)

# Global variables
input_AF.FN <- "03_results/per_family_inferred_allele_frequency_data.RData"
date <- format(Sys.time(), "%Y-%m-%d_%Hh%M")

# User-set variables
#pheno_of_interest <- "dw_size_mean"
pheno_of_interest <- "dw_per_d_tot"
#pheno_of_interest <- "dw_per_d_perf"
#pheno_of_interest <- "sw_per_d_perf"
#pheno_of_interest <- "dw_minus_sw_size_mean"
#pheno_of_interest <- "sw_size_mean"

# Prepare an output folder
output.dir <- paste0("03_results/gemma_run_", pheno_of_interest, "_", date)
dir.create(path = output.dir)


#### 01. Prepare marker info and genotypes ####
# Load allele frequency data
load(input_AF.FN)
loci.df <- geno_sim
dim(loci.df)
loci.df[1:5,1:5]

# Load marker information that matches genotypes
my_vcf <- read.vcfR(file = genos.FN)
my_vcf_info.df <- my_vcf@fix # obtain details on markers
my_vcf_info.df <- as.data.frame(my_vcf_info.df)
my_vcf_info.df[1:5,1:5]
my_vcf_info.df$mname <- paste0(my_vcf_info.df$CHROM, "_", my_vcf_info.df$POS) # create identifier
my_vcf_info.df <- my_vcf_info.df[,c("mname", "REF", "ALT")] # only keep the required columns
head(my_vcf_info.df)
my_vcf_info.df$mname <- gsub(pattern = "\\.", replacement = "_",x = my_vcf_info.df$mname) # update the identifier to match genotypic data
head(my_vcf_info.df)

# Check concordance between the marker info and the genotypes
nrow(my_vcf_info.df)
nrow(loci.df)
length(intersect(x = my_vcf_info.df$mname, y = colnames(loci.df))) # OK (some filtered out)
colnames(loci.df) <- gsub(pattern = "\\.1", replacement = "", x = colnames(loci.df)) # Remove the allele detail from colnames
length(intersect(x = my_vcf_info.df$mname, y = colnames(loci.df))) # OK (some filtered out)
loci.df[1:5,1:5]
head(my_vcf_info.df)

# Prepare geno data
gemma_geno_input.df <- t(loci.df)
gemma_geno_input.df[1:5,1:5]
gemma_geno_input.df <- as.data.frame(gemma_geno_input.df)
gemma_geno_input.df$loci <- rownames(gemma_geno_input.df)

# Combine the marker info and the genotypes
geno <- merge(x = my_vcf_info.df, y = gemma_geno_input.df, by.x = "mname", by.y = "loci")
geno[1:5,1:5] # BIMBAM Format

# Write out geno info
write.table(x = geno, file = paste0(output.dir, "/gemma_geno.txt")
            , sep = "\t"
            , row.names = F, col.names = F
            , quote = F
            )


#### 02. Prepare phenotype information ####
head(cross_and_pheno.df)
cross_and_pheno.df$family.id <- paste0("family_", cross_and_pheno.df$family)
head(cross_and_pheno.df)

geno_order.df <- rownames(loci.df)
geno_order.df <- as.data.frame(geno_order.df)
head(geno_order.df)

cross_and_pheno_ordered.df <-        merge(x = geno_order.df, y = cross_and_pheno.df
                                            , by.x = "geno_order.df", by.y = "family.id"
                                            , sort = F
                                           )
head(cross_and_pheno_ordered.df)
dim(cross_and_pheno_ordered.df)
dim(cross_and_pheno.df)

# New calculation
cross_and_pheno_ordered.df$dw_minus_sw_size_mean <- cross_and_pheno_ordered.df$dw_size_mean - cross_and_pheno_ordered.df$sw_size_mean


### TODO: this is where we could keep more than one pheno of interest if we wanted ###
pheno <- cross_and_pheno_ordered.df[, pheno_of_interest]

write.table(x = pheno
            , file = paste0(output.dir, "/gemma_pheno_", pheno_of_interest, ".txt")
            , quote = F, sep = "\t", row.names = F, col.names = F
            )


### Marker annotation info
# first column is SNP id, second column bp position, third column is chr number
annot.df <- my_vcf@fix
annot.df <- as.data.frame(annot.df)
head(annot.df)
annot.df <- annot.df[,c("CHROM", "POS")]
head(annot.df)
annot.df$mname <- paste0(annot.df$CHROM, "_", annot.df$POS)
annot.df$mname <- gsub(pattern = "\\.", replacement = "_", x = annot.df$mname)
head(annot.df)

# Keep only the records that are in the output file
annot.df <- annot.df[annot.df$mname %in% geno$mname, ]
dim(annot.df)

head(annot.df)
annot.df <- annot.df[,c("mname", "POS", "CHROM")] # put in required order

# ADJUST TO CHR IDS IF WANT (based on Konstantin Divilov's GWAS.R script)
# https://github.com/kdivilov/Aquaculture_2023/blob/main/GWAS.R
annot.df$CHROM <- gsub(pattern = "NC_047559.1", replacement = "Chr07", x = annot.df$CHROM)
annot.df$CHROM <- gsub(pattern = "NC_047560.1", replacement = "Chr01", x = annot.df$CHROM)
annot.df$CHROM <- gsub(pattern = "NC_047561.1", replacement = "Chr09", x = annot.df$CHROM)
annot.df$CHROM <- gsub(pattern = "NC_047562.1", replacement = "Chr06", x = annot.df$CHROM)
annot.df$CHROM <- gsub(pattern = "NC_047563.1", replacement = "Chr03", x = annot.df$CHROM)
annot.df$CHROM <- gsub(pattern = "NC_047564.1", replacement = "Chr02", x = annot.df$CHROM)
annot.df$CHROM <- gsub(pattern = "NC_047565.1", replacement = "Chr04", x = annot.df$CHROM)
annot.df$CHROM <- gsub(pattern = "NC_047566.1", replacement = "Chr05", x = annot.df$CHROM)
annot.df$CHROM <- gsub(pattern = "NC_047567.1", replacement = "Chr10", x = annot.df$CHROM)
annot.df$CHROM <- gsub(pattern = "NC_047568.1", replacement = "Chr08", x = annot.df$CHROM)
#table(annot.df$CHROM)

annot.df <- annot.df[with(annot.df, order(annot.df$CHROM)), ]
head(annot.df)
write.table(x = annot.df, file = paste0(output.dir, "/gemma_geno_annot.txt")
            , sep = "\t", row.names = F, col.names = F, quote = F
            )

# Reporting
print(paste0("Use the output data in ", output.dir, " with gemma to analyze GWAS."))
#TODO: the below is an issue, because it overwrites previous, should either add pheno to name or drop
#save.image(file = "03_results/finalized_gemma_inputs_completed.RData")

# Next: use gemma to analyze data

