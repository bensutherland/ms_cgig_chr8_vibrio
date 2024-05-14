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

# User set variables
input_AF.FN <- "03_results/per_family_inferred_allele_frequency_data.RData"
plink_map.FN <- "../ms_cgig_chr8/02_input_data/myplink.map"
input_VCF.FN <- "../ms_cgig_chr8/02_input_data/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.001_AF_0.05_LD0.5w50kb_subset0.001.vcf"

pheno_of_interest <- "dw_size_mean"



# Load data from previous step
load(input_AF.FN)
dim(loci.df)
loci.df[1:5,1:5]
#loci.df <- loci.df[, grep(pattern = "rep", x = colnames(loci.df), invert = T)]

# Read in VCF
my_vcf <- read.vcfR(file = input_VCF.FN)
str(my_vcf@fix)
my_vcf_info.df <- my_vcf@fix
my_vcf_info.df <- as.data.frame(my_vcf_info.df)
my_vcf_info.df[1:5,1:5]
my_vcf_info.df$mname <- paste0(my_vcf_info.df$CHROM, "_", my_vcf_info.df$POS)
my_vcf_info.df <- my_vcf_info.df[,c("mname", "REF", "ALT")]
head(my_vcf_info.df)
my_vcf_info.df$mname <- gsub(pattern = "\\.", replacement = "_",x = my_vcf_info.df$mname)
head(my_vcf_info.df)

# Check concordance
nrow(my_vcf_info.df)
nrow(loci.df)
length(intersect(x = my_vcf_info.df$mname, y = loci.df$loci)) # OK (some filtered out)

colnames(loci.df) <- gsub(pattern = "_a2_ppn", replacement = "", x = colnames(loci.df))
loci.df[1:5,1:5]
head(my_vcf_info.df)

# Combine the annotation info to the VCF info
geno <- merge(x = my_vcf_info.df, y = loci.df, by.x = "mname", by.y = "loci")
geno[1:5,1:5]

write.table(x = geno, file = "03_results/gemma_geno.txt", sep = "\t"
            , row.names = F, col.names = F
            , quote = F
            )



#### FAMILY INFO AND PHENOS
head(cross_and_pheno.df)
cross_and_pheno.df$family.id <- paste0("family_", cross_and_pheno.df$family)
head(cross_and_pheno.df)

geno_order.df <- colnames(loci.df)[2:ncol(loci.df)]
geno_order.df <- as.data.frame(geno_order.df)
head(geno_order.df)

cross_and_pheno_ordered.df <-        merge(x = geno_order.df, y = cross_and_pheno.df
                                            , by.x = "geno_order.df", by.y = "family.id"
                                            , sort = F)
head(cross_and_pheno_ordered.df)
pheno <- cross_and_pheno_ordered.df[,pheno_of_interest]
write.table(x = pheno, file = paste0("03_results/gemma_pheno_", pheno_of_interest, ".txt")
            , quote = F, sep = "\t", row.names = F, col.names = F
            )


### Annotation info
# first column is SNP id, second column bp position, third column is chr number
# annot.df <- geno$mname
# annot.df <- as.data.frame(annot.df)
# colnames(annot.df)[1] <- "mname"
# head(annot.df)
# annot.df <- separate(data = annot.df, col = "mname", into = c("chr", "pos"), sep = "_1_")

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

# ADJUST TO CHR IDS IF WANT (based on Konstantin's GWAS.R script)
annot.df$CHROM <- gsub(pattern = "NC_047559.1", replacement = "Chr7", x = annot.df$CHROM)
annot.df$CHROM <- gsub(pattern = "NC_047560.1", replacement = "Chr1", x = annot.df$CHROM)
annot.df$CHROM <- gsub(pattern = "NC_047561.1", replacement = "Chr9", x = annot.df$CHROM)
annot.df$CHROM <- gsub(pattern = "NC_047562.1", replacement = "Chr6", x = annot.df$CHROM)
annot.df$CHROM <- gsub(pattern = "NC_047563.1", replacement = "Chr3", x = annot.df$CHROM)
annot.df$CHROM <- gsub(pattern = "NC_047564.1", replacement = "Chr2", x = annot.df$CHROM)
annot.df$CHROM <- gsub(pattern = "NC_047565.1", replacement = "Chr4", x = annot.df$CHROM)
annot.df$CHROM <- gsub(pattern = "NC_047566.1", replacement = "Chr5", x = annot.df$CHROM)
annot.df$CHROM <- gsub(pattern = "NC_047567.1", replacement = "Chr10", x = annot.df$CHROM)
annot.df$CHROM <- gsub(pattern = "NC_047568.1", replacement = "Chr8", x = annot.df$CHROM)
#table(annot.df$CHROM)

annot.df <- annot.df[with(annot.df, order(annot.df$CHROM)), ]
head(annot.df)
write.table(x = annot.df, file = "03_results/gemma_geno_annot.txt"
            , sep = "\t", row.names = F, col.names = F, quote = F
            )


