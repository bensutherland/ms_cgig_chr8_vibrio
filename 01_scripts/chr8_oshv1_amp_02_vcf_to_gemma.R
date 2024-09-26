# Read in de novo genotyped amp panel VCF, characterize, then prepare GWAS gemma files
#  requires amplitools-based BCF file, renamed by impute_workflow, and phenotype file (in 00_archive)
#   convert to VCF file, and put in 02_input_data (see filenames in 'User set variables' below)
#  note: all code and directories listed are within the simple_pop_stats repository
#  initialized 2024-06-14
#  Ben J. G. Sutherland (VIU), incl. code dev by Konstantin Divilov


#### 00. Front Matter ####
# Clear space
# rm(list=ls())

# Source simple_pop_stats ( https://github.com/bensutherland/simple_pop_stats )

# Load additional libraries to those loaded in simple_pop_stats
#devtools::install_github('kaustubhad/fastman',build_vignettes = TRUE)
library(fastman)
#install.packages("missMethods")
library(missMethods)

## User set variables
# Filenames
phenos.FN       <- "00_archive/G0923-21-VIUN_SampleInventory_V2_recd_2024-08-16.txt"
vcf.FN          <- "02_input_data/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_offspring_only_rename.vcf"

# Variables
max_missing <- 0.3
impute <- TRUE
filter_by_GR <- FALSE # want to filter by GR? 

# Phenotype variable
pheno_of_interest <- "survival_state"
#pheno_of_interest <- "DPE"


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
table(paste0(pheno.df$family, "__", pheno.df$survival_state, "__", pheno.df$DPE), useNA = "ifany") # view by day as well

# Remove samples with missing phenotypes
pheno.df <- pheno.df[!is.na(pheno.df$DPE), ]
nrow(pheno.df)

# Convert survivor DPE to day 17 (one day after the trial ended)
max(pheno.df$DPE, na.rm = T)
pheno.df[pheno.df$survival_state=="S", "DPE"] <- 17
head(pheno.df, n = 20)


#### 02. Load genotypes ####
# Load input VCF
my_vcf <- read.vcfR(file = vcf.FN)
my_vcf

# Convert to genind for simple_pop_stats functions
obj <- vcfR2genind(x = my_vcf)
obj
head(indNames(x = obj)) # indiv names

# Only retain those individuals with phenotypic information
obj <- obj[indNames(obj) %in% pheno.df$indiv,]
obj


#### 03. Update population ID ####
# Annotate population based on family
population <- indNames(x = obj)
population[grep(pattern = "-114-", x = population)] <- "F114"
population[grep(pattern = "-115-", x = population)] <- "F115"
population[grep(pattern = "-116-", x = population)] <- "F116"
population[grep(pattern = "-117-", x = population)] <- "F117"
population
table(population)

# Update obj pop ID
pop(obj) <- population
rm(population)

# Create df w/ info for later matching
population.df <- cbind(indNames(obj), as.character(obj@pop))
colnames(population.df) <- c("indiv", "family")
head(population.df)


#### 04. Per individual missing data ####
##### 04.a. Characterize missing data #####
percent_missing_by_ind(df = obj)
head(missing_data.df)

# Add population and survivorship info
missing_data.df <- merge(x = missing_data.df, y = pheno.df, by.x = "ind", by.y = "indiv", sort = F)
head(missing_data.df)
nrow(missing_data.df) 

# Add a colour for mortality
missing_data.df$mort_col <- NA
missing_data.df$mort_col[grep(pattern = "M", x = missing_data.df$survival_state)] <- "red"
missing_data.df$mort_col[grep(pattern = "S", x = missing_data.df$survival_state)] <- "black"

# Add a shape for mortality
missing_data.df$mort_shape <- NA
missing_data.df$mort_shape[grep(pattern = "M", x = missing_data.df$survival_state)] <- 1
missing_data.df$mort_shape[grep(pattern = "S", x = missing_data.df$survival_state)] <- 16


##### 04.b. Plot missing data #####
# Plot missing data by individual, colour by family
pdf(file = "03_results/geno_rate_by_ind.pdf", width = 9, height = 5)
plot(100 * (1 - missing_data.df$ind.per.missing), ylab = "Genotyping rate (%)"
     , col = as.factor(missing_data.df$family)
     , las = 1
     , xlab = "Individual"
     , ylim = c(0,100)
     , pch=missing_data.df$mort_shape
     , cex = 0.8
     #, xlim = c(0,250)
)
abline(h = (100 * (1 - max_missing)), lty = 2)

legend("bottomleft", legend = unique(missing_data.df$family)
       , fill = as.factor(unique(missing_data.df$family))
       , cex = 0.8
       , bg = "white"
)
dev.off()


##### 04.c. Filter based on missing data #####
filtered_data.df <- missing_data.df[missing_data.df$ind.per.missing < max_missing, ]
nrow(filtered_data.df)

# View summary
table(filtered_data.df$family)
table(paste0(filtered_data.df$family, "__", filtered_data.df$survival_state))

# Identify which inds should be kept
keep <- missing_data.df[missing_data.df$ind.per.missing < max_missing, "ind"]
length(keep)

# Filter obj to only keep inds w/ GR > cutoff
obj <- obj[(keep)]
obj


#### 05. Per locus missing data ####
##### 05.a. Characterize missing data #####
percent_missing_by_locus(df = obj)
head(missing_data_loci.df)

##### 05.b. Plot missing data #####
pdf(file = "03_results/geno_rate_by_marker.pdf", width = 9, height = 5)
hist(missing_data_loci.df$perc.missing, breaks = 20
     , las = 1
     , xlab = "Per locus missing (%)"
     , main = ""
     )
dev.off()

##### 05.c. Filter based on missing data #####
# Filter markers by genotyping rate
keep <- missing_data_loci.df[missing_data_loci.df$perc.missing < max_missing, "locus.id"]
length(keep)

# How many loci will be removed? 
nLoc(obj)
nLoc(obj) - length(keep)

# Retain only those loci above the cutoff
obj <- obj[, loc=keep]


#### 06. Drop monomorphic loci ####
drop_loci(df = obj, drop_monomorphic = T)
obj <- obj_filt # rename output back to original
obj


#### 07. Run PCA ####
pca_from_genind(data = obj
                , PCs_ret = 4
                , plot_eigen = T
                , plot_allele_loadings = F
                , retain_pca_obj = T
                )
# Results are saved out to PDF file


#### 08. Filter VCF based on above filters ####
if(filter_by_GR==TRUE){
  
  print("Filtering VCF file by per-sample genotyping rate")
  my_vcf <- my_vcf[, c("FORMAT", indNames(obj))]
  print(my_vcf)
  
}else{
  
  print("Not filtering VCF file by per-sample genotyping rate")
  print("Removing inds from VCF file that did not have phenotypes")
  my_vcf <- my_vcf[, c("FORMAT", pheno.df$indiv)]
  print(my_vcf)
  
}

# Filter to only keep specific loci
# It is possible to filter based on the fix col, would need to populate the ID column via paste, then use similar to the following
#my_vcf[which(paste0(my_vcf@fix[,"ID"] %in% as.character(locNames(obj))), ]
# note: not implementing 


#### 09. Prepare GEMMA inputs ####
# Note: the following assumes that the order of the genotypes is retained to the fix section, 
#  which it should be, assuming you are drawing from the same VCF file

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


##### 09.b. Mean imputing (optional) #####
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


# Create vector with phenotypic information
indiv.df <- rownames(geno)
indiv.df <- as.data.frame(indiv.df)
colnames(indiv.df) <- "indiv"
dim(indiv.df)
head(indiv.df)

## Merge with earlier-developed pheno df
# Checking
setdiff(x = indiv.df$indiv, y = pheno.df$indiv) # if any inds are present without phenotypes, something has gone wrong in filters above, start over

# Merge (sort = F is needed)
indiv.df <- merge(x = indiv.df, y = pheno.df, by = "indiv", all.x = T, sort = F)
# NOTE: if any NAs, this will not work

# Checking
head(indiv.df) # in the same order as the geno df
geno[1:5,1:5]
tail(indiv.df)
table(paste0(indiv.df$family, "__", indiv.df$survival_state))


#### TODO: MOVE UP and SAVE OUT####
# Show mortality by family in retained samples
boxplot(pheno.df$DPE ~ as.factor(pheno.df$family), las = 1, ylab = "DPE"
        , xlab = "Family"
        )
#### /END/ TODO: MOVE ####

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
write.table(x = gwasgeno, file = "03_results/gwas_geno.txt", row.names = F, col.names = F, quote = F)

# Next: go to GEMMA for analysis (see README)

