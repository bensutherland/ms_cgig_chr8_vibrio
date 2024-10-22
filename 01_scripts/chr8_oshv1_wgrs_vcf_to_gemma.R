# Read in whole-genome resequence VCF file (empirical, 10X data), characterize, then prepare GWAS gemma files
#  requires wgrs_workflow-based VCF file (in 02_input_data), and phenotype file (in 00_archive)
#  note: all code and directories listed are within the simple_pop_stats repository
#  note: all output will go into simple_pop_stats/03_results
#  initialized 2024-10-07
#  Ben J. G. Sutherland (VIU), incl. code dev by Konstantin Divilov


#### 00. Front Matter ####
# Clear space
# rm(list=ls())

## Install packages
# install.packages("vcfR")
# install.packages("data.table")
# install.packages("tidyr")

## Load libraries
library(vcfR)
library(data.table)
library(tidyr)

# Set working directory to the ms_scallop_popgen repo
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

## User set variables
phenos.FN  <- "00_archive/G0923-21-VIUN_SampleInventory_V2_recd_2024-08-16.txt"
vcf.FN     <- "02_input_data/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_offspring_only_rename.vcf"

impute             <- FALSE
family_sep_outputs <- TRUE

# Obtain date
date <- format(Sys.time(), "%Y-%m-%d")


#### 01. Load phenotype info ####
# Load and prepare phenotype file
pheno.df       <- read.delim2(file = phenos.FN, header = T, sep = "\t")
pheno.df       <- as.data.frame(pheno.df)
pheno.df$DPE   <- as.numeric(pheno.df$DPE) # convert DPE to numeric
pheno.df$indiv <- gsub(pattern = "_", replacement = "-", x = pheno.df$indiv) # for matching sample IDs in VCF file
pheno.df       <- separate(data = pheno.df, col = "indiv"
                     , into = c("assay", "family", "replicate", "sample.id")
                     , sep = "-", remove = F) # Obtain details about each sample
head(pheno.df)

# Summarize available phenotypes
table(pheno.df$family, useNA = "ifany")                                        # indivs per family
table(pheno.df$survival_state, useNA = "ifany")                                # summarize number inds w/ each survival state
table(paste0(pheno.df$family, "__", pheno.df$survival_state), useNA = "ifany") # summarize mortality state by family
#table(paste0(pheno.df$family, "__", pheno.df$survival_state, "__", pheno.df$DPE), useNA = "ifany") # view by day as well

# Remove samples with missing phenotypes
pheno.df <- pheno.df[!is.na(pheno.df$DPE), ]
nrow(pheno.df)

# Convert survivor DPE to day 17 (one day after the trial ended)
max(pheno.df$DPE, na.rm = T)
pheno.df[pheno.df$survival_state=="S", "DPE"] <- 17
head(pheno.df, n = 20)

# Show mortality by family in retained samples
pdf(file = "03_results/DPE_by_family.pdf", width = 7.5, height = 4.5)
boxplot(pheno.df$DPE ~ as.factor(pheno.df$family), las = 1, ylab = "DPE"
        , xlab = "Family"
)
dev.off()


#### 02. Load genotypes ####
# Load input VCF
my_vcf <- read.vcfR(file = vcf.FN)


# #### 03. Prepare GEMMA inputs ####
# Extract genotypes
geno <- extract.gt(x = my_vcf, element = "GT")
geno[1:5, 1:5]

# Extract positional info
positional_info.df <- my_vcf@fix
positional_info.df <- as.data.frame(positional_info.df)
head(positional_info.df)

# Use positional info as the locus.id
locus.id <- paste0(positional_info.df$CHROM, "__", positional_info.df$POS)
rownames(geno) <- locus.id # use the new ID as the row names
rm(locus.id)
geno[1:5, 1:5]

# Convert genotypes to numeric allele dosage
geno[geno=="0/0"] = 0
geno[geno=="0/1"] = 1
geno[geno=="1/1"] = 2
geno[1:5, 1:5]

geno = t(geno) # rows = samples, cols = loci
geno[1:5, 1:5]

mode(geno) = "numeric"
geno[1:5,1:5]

# Optional imputation 
if(impute==TRUE){

  #### Imputation ####
  # Create family-specific genotype matrices
  geno_F114 = geno[grep("-114-", rownames(geno)), ]
  dim(geno_F114)
  geno_F115 = geno[grep("-115-", rownames(geno)), ]
  dim(geno_F115)
  geno_F116 = geno[grep("-116-", rownames(geno)), ]
  dim(geno_F116)
  geno_F117 = geno[grep("-117-", rownames(geno)), ]
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

}else if(impute!=TRUE){
 
  print("Not imputing, expect missing data.")
   
}


# Obtain vector of family identities
head(rownames(geno))
var_family <- rep(NA, nrow(geno))
var_family[grep(pattern = "-114-", x = rownames(geno))] <- "F114"
var_family[grep(pattern = "-115-", x = rownames(geno))] <- "F115"
var_family[grep(pattern = "-116-", x = rownames(geno))] <- "F116"
var_family[grep(pattern = "-117-", x = rownames(geno))] <- "F117"
table(var_family)

# Retain as covariate
gwascovar = model.matrix(~as.factor(var_family))

# Obtain the order of sample names for the geno object and put it into a df
indiv.df <- rownames(geno)
indiv.df <- as.data.frame(indiv.df)
colnames(indiv.df) <- "indiv"
head(indiv.df)
nrow(indiv.df)

## Merge with earlier-developed pheno df
# Checking
setdiff(x = indiv.df$indiv, y = pheno.df$indiv) # if any inds are present without phenotypes, something has gone wrong in filters above, start over

# Merge (sort = F is needed)
indiv.df <- merge(x = indiv.df, y = pheno.df, by = "indiv", all.x = T, sort = F)
# NOTE: if any NAs, this will not work

# Checking
head(indiv.df) # in the same order as the geno df
dim(indiv.df)
geno[1:5,1:5]
tail(indiv.df)
table(paste0(indiv.df$family, "__", indiv.df$survival_state))


# Create pheno objects
indiv.df$survival_state_gemma <- indiv.df$survival_state
indiv.df$survival_state.gemma <- gsub(pattern = "S", replacement = "1", x = indiv.df$survival_state_gemma)
indiv.df$survival_state.gemma <- gsub(pattern = "M", replacement = "0", x = indiv.df$survival_state.gemma)

pheno_survival_state <- as.numeric(indiv.df$survival_state.gemma)
pheno_survival_state
  

## Retain phenotype for DPE
pheno_DPE <- as.numeric(indiv.df$DPE)
pheno_DPE  

## Prepare genotypes
gwasgeno = t(geno) 
gwasgeno[1:5,1:5] # preview
gwasgeno <- cbind(rownames(gwasgeno),"X","Y",gwasgeno)

# Write outputs
write.table(x = pheno_survival_state, file = "03_results/gwas_pheno_survival_state.txt", row.names = F, col.names = F)
write.table(x = pheno_DPE, file = "03_results/gwas_pheno_DPE.txt", row.names = F, col.names = F)
write.table(x = gwascovar, file = "03_results/gwas_covar.txt", row.names = F, col.names = F)
#fwrite(x = gwasgeno, file = "03_results/gwas_geno.txt", sep = " ", col.names = F, quote = F) # problem: NAs are empty
write.table(x = gwasgeno, file = "03_results/gwas_geno.txt", row.names = F, col.names = F, quote = F)


# Also creating output as family-separated data? 
if(family_sep_outputs==TRUE){
  
  # Make an output directory
  output.dir <- paste0("03_results/family_sep_gemma_inputs_", date)
  dir.create(path = output.dir)
  
  # Genotypes
  F114_gwasgeno <- gwasgeno[, c(1:3, grep(pattern = "\\-114\\-", x = colnames(gwasgeno)))] # 1:3 gives the first three added cols
  F114_gwasgeno[1:5,1:5]
  dim(F114_gwasgeno)
  write.table(x = F114_gwasgeno, file = paste0(output.dir, "/gwas_geno_F114.txt"), row.names = F, col.names = F, quote = F)
  
  F115_gwasgeno <- gwasgeno[, c(1:3, grep(pattern = "\\-115\\-", x = colnames(gwasgeno)))] # 1:3 gives the first three added cols
  F115_gwasgeno[1:5,1:5]
  dim(F115_gwasgeno)
  write.table(x = F115_gwasgeno, file = paste0(output.dir, "/gwas_geno_F115.txt"), row.names = F, col.names = F, quote = F)
  
  F116_gwasgeno <- gwasgeno[, c(1:3, grep(pattern = "\\-116\\-", x = colnames(gwasgeno)))] # 1:3 gives the first three added cols
  F116_gwasgeno[1:5,1:5]
  dim(F116_gwasgeno)
  write.table(x = F116_gwasgeno, file = paste0(output.dir, "/gwas_geno_F116.txt"), row.names = F, col.names = F, quote = F)
  
  F117_gwasgeno <- gwasgeno[, c(1:3, grep(pattern = "\\-117\\-", x = colnames(gwasgeno)))] # 1:3 gives the first three added cols
  F117_gwasgeno[1:5,1:5]
  dim(F117_gwasgeno)
  write.table(x = F117_gwasgeno, file = paste0(output.dir, "/gwas_geno_F117.txt"), row.names = F, col.names = F, quote = F)
  
  # Save out per-family phenotypes too
  pheno_survival_state_F114 <- as.numeric(indiv.df[indiv.df$family=="114", "survival_state.gemma"])
  length(pheno_survival_state_F114)
  write.table(x = pheno_survival_state_F114, file = paste0(output.dir, "/gwas_pheno_survival_state_F114.txt"), row.names = F, col.names = F)
  
  pheno_survival_state_F115 <- as.numeric(indiv.df[indiv.df$family=="115", "survival_state.gemma"])
  length(pheno_survival_state_F115)
  write.table(x = pheno_survival_state_F115, file = paste0(output.dir, "/gwas_pheno_survival_state_F115.txt"), row.names = F, col.names = F)
  
  pheno_survival_state_F116 <- as.numeric(indiv.df[indiv.df$family=="116", "survival_state.gemma"])
  length(pheno_survival_state_F116)
  write.table(x = pheno_survival_state_F116, file = paste0(output.dir, "/gwas_pheno_survival_state_F116.txt"), row.names = F, col.names = F)
  
  pheno_survival_state_F117 <- as.numeric(indiv.df[indiv.df$family=="117", "survival_state.gemma"])
  length(pheno_survival_state_F117)
  write.table(x = pheno_survival_state_F117, file = paste0(output.dir, "/gwas_pheno_survival_state_F117.txt"), row.names = F, col.names = F)
  
}

# Next: go to GEMMA for analysis (see README)
