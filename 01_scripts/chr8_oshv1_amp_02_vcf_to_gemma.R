# Prepare gemma files based on input VCF and phenotypes
#  requires (...)
#  initialized 2024-06-14
#  Ben J. G. Sutherland (VIU), incl. code dev by Konstantin Divilov

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

# Source simple_pop_stats for functions (but use the ms_cgig_chr8 setwd command below)

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
# inds_to_keep.FN <- "../simple_pop_stats/03_results/retained_inds.txt"
# loci_to_keep.FN <- "../simple_pop_stats/03_results/retained_loci.txt"

max_missing <- 0.3

#### 01. Load pheno info ####
pheno.df <- read.delim2(file = phenos.FN, header = T, sep = "\t")
pheno.df <- as.data.frame(pheno.df)
head(pheno.df)
pheno.df <- separate(data = pheno.df, col = "tube_label", into = c("assay", "family", "replicate", "indiv"), sep = "_", remove = F)
table(pheno.df$family, useNA = "ifany") # see if any NAs
table(pheno.df$mort_surv, useNA = "ifany") # see if any NAs
table(paste0(pheno.df$family, "__", pheno.df$mort_surv))


#### 02. Load VCF file and do filtering, summary analysis ####
# Load input VCF
my_vcf <- read.vcfR(file = vcf.FN)
my_vcf

# Convert to genind for sps functions
obj <- vcfR2genind(x = my_vcf)
obj
head(indNames(x = obj))

# Annotate population based on family
indiv.df <- as.data.frame(indNames(x = obj))
colnames(indiv.df) <- "indiv"
head(indiv.df)
indiv.df <- separate(data = indiv.df, col = "indiv", into = c("assay", "family", "replicate", "sample_number"), sep = "_", remove = F)
head(indiv.df)

# Characterize missing data
percent_missing_by_ind(df = obj)
head(missing_data.df)
missing_data.df$pop <- NA
missing_data.df[grep(pattern = "_114_", x = missing_data.df$ind), "pop"] <- "F114"
missing_data.df[grep(pattern = "_115_", x = missing_data.df$ind), "pop"] <- "F115"
missing_data.df[grep(pattern = "_116_", x = missing_data.df$ind), "pop"] <- "F116"
missing_data.df[grep(pattern = "_117_", x = missing_data.df$ind), "pop"] <- "F117"
head(missing_data.df)
tail(missing_data.df)

head(pheno.df) # add pheno info
missing_data.df <- merge(x = missing_data.df, y = pheno.df, by.x = "ind", by.y = "tube_label", all.x = T, sort = F)
head(missing_data.df)

# Plot missing data
# Plot missing data by individual, colour by pop
pdf(file = "03_results/geno_rate_by_ind_pop.pdf", width = 9, height = 6)
par(mfrow=c(1,1))
missing_data.df <- missing_data.df[with(missing_data.df, order(missing_data.df$pop)), ]
plot(100 * (1 - missing_data.df$ind.per.missing), ylab = "Genotyping rate (%)"
     , col = as.factor(missing_data.df$pop)
     , las = 1
     , xlab = "Individual"
     , ylim = c(0,100)
     , pch=19
     , cex = 0.8
)
abline(h = 50, lty = 3)
abline(h = 70, lty = 2)

legend("bottomleft", legend = unique(missing_data.df$pop)
       , fill = as.factor(unique(missing_data.df$pop))
       , cex = 0.8
       , bg = "white"
)
dev.off()

# Plot by mortality
head(missing_data.df)
missing_data.df$mort_col <- NA
missing_data.df$mort_col[grep(pattern = "M", x = missing_data.df$mort_surv)] <- "red"
missing_data.df$mort_col[grep(pattern = "S", x = missing_data.df$mort_surv)] <- "lightblue"
missing_data.df$mort_col[is.na(missing_data.df$mort_surv)] <- "black"

pdf(file = "03_results/geno_rate_by_ind_pheno.pdf", width = 9, height = 6)
plot(100 * (1 - missing_data.df$ind.per.missing), ylab = "Genotyping rate (%)"
     , col = missing_data.df$mort_col
     , las = 1
     , xlab = "Individual"
     , ylim = c(0,100)
     , pch=19
     , cex = 0.8
)
abline(h = 50, lty = 3)
abline(h = 70, lty = 2)

legend(x = 15, y = 30, legend = c("m", "s", "NA")
       , fill = unique(missing_data.df$mort_col)
       , cex = 0.8
       #, bg = ""
      )
dev.off()

# Summarize samples before any filtering is done
table(paste0(pheno.df$family, "__", pheno.df$mort_surv))


# Filter based on missing data
head(missing_data.df)
filtered_data.df <- missing_data.df[missing_data.df$ind.per.missing < max_missing, ]
table(paste0(filtered_data.df$pop, "__", filtered_data.df$mort_surv))

keep <- missing_data.df[missing_data.df$ind.per.missing < max_missing, "ind"]

length(keep)
nInd(obj)

obj.filt <- obj[(keep)]
obj.filt


##### Loci - missing data #####
# Filter loci based on missing data
obj.df <- genind2df(obj.filt)
obj.df[1:5,1:5]
obj.df <- t(obj.df)
obj.df[1:5,1:5]
obj.df <- obj.df[2:nrow(obj.df),] # remove pop row
obj.df[1:5,1:5]
dim(obj.df)
str(obj.df)

obj.df <- as.data.frame(obj.df)
dim(obj.df)
str(obj.df)
obj.df[1:5,1:5] # See top left of file
obj.df[(dim(obj.df)[1]-5):dim(obj.df)[1], (dim(obj.df)[2]-5):dim(obj.df)[2]] # See bottom right of file

# Add collector col
obj.df$marker.per.missing <- NA

for(i in 1:(nrow(obj.df))){
  
  # Per marker                      sum all NAs for the marker, divide by total number markers
  obj.df$marker.per.missing[i] <-  (sum(is.na(obj.df[i,]))-1) / (ncol(obj.df)-1) 
  
}


# Plot missing data by marker
pdf(file = "03_results/geno_rate_by_marker.pdf", width = 9, height = 6)
plot(100 * (1- obj.df$marker.per.missing), xlab = "Marker", ylab = "Genotyping rate (%)", las = 1
     , ylim = c(0,100)
     #, pch = 1
     #, cex = plot_cex
)
abline(h = 50
       #, col = "grey60"
       , lty = 3)
dev.off()

# Filter markers by genotyping rate
keep <- rownames(obj.df[obj.df$marker.per.missing < max_missing, ])

# How many loci will be removed? 
nLoc(obj.filt)
nLoc(obj.filt) - length(keep)

# Drop loci from genind
obj.all.filt <- obj.filt[, loc=keep]

# Rename back to obj
obj <- obj.all.filt
obj



#obj.bck <- obj
drop_loci(drop_monomorphic = T) ### TODO: need to require an input object identifier
# Drops 174 monomorphic loci
obj <- obj_filt
obj

population <- rep(NA, times = nInd(obj))
population[grep(pattern = "114", x = indNames(obj))] <- "F114"
population[grep(pattern = "115", x = indNames(obj))] <- "F115"
population[grep(pattern = "116", x = indNames(obj))] <- "F116"
population[grep(pattern = "117", x = indNames(obj))] <- "F117"
population

pop(obj) <- population

## Run PCA
pca_from_genind(data = obj, PCs_ret = 4, plot_eigen = T, plot_allele_loadings = F, retain_pca_obj = T)

# Save info
loci_to_keep <- locNames(obj)
length(loci_to_keep)
inds_to_keep <- indNames(obj)
length(inds_to_keep)

# Additional step: do any of the indiv set to be retained not have phenotypic info? 
inds_to_keep %in% pheno.df$tube_label # there is a missing one
inds_to_keep <- inds_to_keep[inds_to_keep %in% pheno.df$tube_label]
length(inds_to_keep)

##### Filter VCF ####
# Filter to only keep specific inds
my_vcf.filt <- my_vcf[, c("FORMAT", inds_to_keep)]

# Filter to only keep specific loci
my_vcf.filt <- my_vcf.filt[which(my_vcf@fix[,"ID"] %in% as.character(loci_to_keep)), ]

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

#### Imputation ####
# Create family-specific genotype matrices
geno_F114 = geno[grep("_114_", rownames(geno)), ]
dim(geno_F114)
geno_F115 = geno[grep("_115_", rownames(geno)), ]
dim(geno_F115)
geno_F116 = geno[grep("_116_", rownames(geno)), ]
dim(geno_F116)
geno_F117 = geno[grep("_117_", rownames(geno)), ]
dim(geno_F117)

# Run family-specific mean imputation on genotypes
geno_F114_impute <- impute_mean(geno_F114)
geno_F115_impute <- impute_mean(geno_F115)
geno_F116_impute <- impute_mean(geno_F116)
geno_F117_impute <- impute_mean(geno_F117)


# Reconstruct full genotype matrix
geno_imputed <- rbind(  geno_F114_impute
                      , geno_F115_impute
                      , geno_F116_impute
                      , geno_F117_impute
                      )

# this is the genotypes file needed for gemma
geno <- geno_imputed

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




