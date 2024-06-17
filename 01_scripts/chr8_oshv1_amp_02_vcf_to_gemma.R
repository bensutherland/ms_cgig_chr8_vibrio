# Prepare gemma files based on input VCF and phenotypes
#  requires (...)
#  initialized 2024-06-14
#  Ben J. G. Sutherland (VIU)

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

# Load libraries
library(vcfR)
#devtools::install_github('kaustubhad/fastman',build_vignettes = TRUE)
library(fastman)
#install.packages("missMethods")
library(missMethods)

## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path) # take main directory
setwd(current.path)
rm(current.path)

# User set variables
phenos.FN       <- "00_archive/qcat992_sample_mort_pheno_2024-06-17.txt"
vcf.FN          <- "02_input_data/G0923-21-VIUN_annot_snplift_to_roslin.vcf"
inds_to_keep.FN <- "../simple_pop_stats/03_results/retained_inds.txt"
loci_to_keep.FN <- "../simple_pop_stats/03_results/retained_loci.txt"

# Load pheno info
pheno.df <- read.delim2(file = phenos.FN, header = T, sep = "\t")
pheno.df <- as.data.frame(pheno.df)
head(pheno.df)
pheno.df <- separate(data = pheno.df, col = "tube_label", into = c("assay", "family", "replicate", "indiv"), sep = "_", remove = F)
table(pheno.df$family, useNA = "ifany") # see if any NAs
table(pheno.df$mort_surv, useNA = "ifany") # see if any NAs
table(paste0(pheno.df$family, "__", pheno.df$mort_surv))

# Load input VCF
my_vcf <- read.vcfR(file = vcf.FN)

# # Are there any individuals that do not have a phenotype? 
# genotyped_indivs <- colnames(my_vcf@gt)[2:ncol(my_vcf@gt)]
# pheno.df[!pheno.df$tube_label %in% genotyped_indivs, "tube_label"]
# pheno.df[!genotyped_indivs %in% pheno.df$tube_label, "tube_label"]


# Load keep info
inds_to_keep <- read.table(file = inds_to_keep.FN)
inds_to_keep.df <- as.data.frame(inds_to_keep)
head(inds_to_keep.df)
colnames(inds_to_keep.df) <- "inds_to_keep"
nrow(inds_to_keep.df) # 173

lost_sample <- as.data.frame("ASY2_114_R4_10")
colnames(lost_sample) <- "lost_sample"
inds_to_keep <- setdiff(x = inds_to_keep.df$inds_to_keep, y = lost_sample$lost_sample)
length(inds_to_keep) # 172
rm(inds_to_keep.df)

loci_to_keep <- read.table(file = loci_to_keep.FN)
loci_to_keep.df <- as.data.frame(loci_to_keep)
head(loci_to_keep.df)
colnames(loci_to_keep.df) <- "loci_to_keep"


# Filter to only keep specific loci and inds
my_vcf.filt <- my_vcf[, c("FORMAT", inds_to_keep)]

# which(as.character(loci_to_keep.df$loci_to_keep) %in% my_vcf@fix[,"ID"])
# which(my_vcf@fix[,"ID"] %in% as.character(loci_to_keep.df$loci_to_keep))
# length(which(my_vcf@fix[,"ID"] %in% as.character(loci_to_keep.df$loci_to_keep))) # 261 (some diff?)

my_vcf.filt <- my_vcf.filt[which(my_vcf@fix[,"ID"] %in% as.character(loci_to_keep.df$loci_to_keep)), ]

my_vcf.filt
my_vcf <- my_vcf.filt


#### Prepare gemma inputs ####
## Extract genotypes
geno <- extract.gt(x = my_vcf, element = "GT")
geno[1:5, 1:5]

## Convert marker names to chr and pos info
marker_names <- rownames(geno) # pull rownames
marker_names <- as.data.frame(marker_names)
head(marker_names)

positional_info.df <- my_vcf@fix # pull chr pos info
positional_info.df <- as.data.frame(positional_info.df)
head(positional_info.df)

# Combine ordered marker names in VCF with pos info
marker_names_w_positional.df <- merge(x = marker_names, y = positional_info.df
                                      , by.x = "marker_names", by.y = "ID", all.x = T, sort = F
                                      )
head(marker_names_w_positional.df)
tail(marker_names_w_positional.df)

# Create new ID, using chr and pos info
new_id <- paste0(marker_names_w_positional.df$CHROM, "__", marker_names_w_positional.df$start)
rownames(geno) <- new_id # use the new ID as the row names
rm(new_id)
geno[1:5, 1:5]

## Convert genotypes to numeric values
geno[geno=="0/0"] = 0
geno[geno=="0/1"] = 1
geno[geno=="1/1"] = 2
geno = t(geno) # rows = samples, cols = loci
mode(geno) = "numeric"

geno[1:5,1:5]

# this is the genotypes file needed for gemma


## Obtain vector of family identities
head(rownames(geno))
var_family <- rep(NA, nrow(geno))
var_family[grep(pattern = "_114_", x = rownames(geno))] <- "F114"
var_family[grep(pattern = "_115_", x = rownames(geno))] <- "F115"
var_family[grep(pattern = "_116_", x = rownames(geno))] <- "F116"
var_family[grep(pattern = "_117_", x = rownames(geno))] <- "F117"
table(var_family)

# Retain as covariate
gwascovar = model.matrix(~as.factor(var_family))


# Create phenotype holding mort/survival info
indiv.df <- rownames(geno)
indiv.df <- as.data.frame(indiv.df)
dim(indiv.df)

# Do all genotyped individuals have a phenotype? 
setdiff(x = indiv.df$indiv.df, y = pheno.df$tube_label) # if any are present, need to remove at the start of this script and start over
# this was done for ASY2_114_R4_10

# Merge
indiv.df <- merge(x = indiv.df, y = pheno.df, by.x = "indiv.df", by.y = "tube_label", all.x = T, sort = F) 
# NOTE: if any NAs, this will not work
head(indiv.df) # in the same order as the geno df
tail(indiv.df)
dim(indiv.df)
table(indiv.df$mort_surv, useNA = "ifany") # should all have phenos
table(paste0(indiv.df$family, "__", indiv.df$mort_surv))

indiv.df$mort_surv <- gsub(pattern = "S", replacement = "1", x = indiv.df$mort_surv)
indiv.df$mort_surv <- gsub(pattern = "M", replacement = "0", x = indiv.df$mort_surv)
var_status <- as.numeric(indiv.df$mort_surv)
var_status

# GWAS variables
gwaspheno <-  var_status
str(gwascovar) # already present
gwasgeno = t(geno) 
gwasgeno[1:5,1:5] # preview
gwasgeno <- cbind(rownames(gwasgeno),"X","Y",gwasgeno)

gwasanno <- rownames(gwasgeno)
gwasanno <- as.data.frame(gwasanno)
gwasanno <- separate(data = gwasanno, col = "gwasanno", into = c("chr", "pos"), sep = "__", remove = F)
head(gwasanno)


write.table(x = gwaspheno, file = "03_results/gwaspheno.txt", row.names = F, col.names = F)
write.table(x = gwascovar, file = "03_results/gwascovar.txt", row.names = F, col.names = F)
write.table(x = gwasgeno, file = "03_results/gwasgeno.txt", row.names = F, col.names = F, quote = F)
write.table(x = gwasanno, file = "03_results/gwasanno.txt", row.names = F, col.names = F, quote = F)




