# Read in AlphaImpute outputs (imputed or 10X, all chromosomes), and compare to evaluate imputation
# B. Sutherland (2024-07-25)

### Front Matter ####
# Clean space
# rm(list=ls())

## Install and load packages
#install.packages("rstudioapi")
#install.packages("data.table")
library("rstudioapi")
library("data.table")
library(tidyr)
library(dplyr)

# Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

## Info
# sessionInfo()

# Set variables
offspring_imputed_ai2.FN     <- "13_impute_compare/all_chr_combined.txt" # imputed
offspring_10X_ai2.FN         <- "13_impute_compare/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_ai2.txt" # 10X genotypes

#### 01. Load data ####
# Read in imputed data
imputed.df <- fread(file = offspring_imputed_ai2.FN, sep = "\t")
dim(imputed.df)
imputed.df <- as.data.frame(imputed.df) # convert to df
imputed.df[1:5,1:5]

# Read in empirical data (10X)
empirical.df <- fread(file = offspring_10X_ai2.FN, sep = "\t")
dim(empirical.df)
empirical.df <- as.data.frame(empirical.df) # convert to df
empirical.df[1:5,1:5]


#### 02. Prepare data for matching ####
## Make sample names match for imputed data
head(colnames(imputed.df), n = 20)

# Remove parents and grandparents (if present)
imputed.df <- imputed.df[, grep(pattern = "mname|ASY2", x = colnames(imputed.df))] # remove parents, keep mname and ASY2 inds
dim(imputed.df)

# Remove 10X data (if present; this was used to support imputation in one evaluation)
imputed.df <- imputed.df[, grep(pattern = "fastq.gz", x = colnames(imputed.df), invert = T)] # remove parents, keep mname and ASY2 inds
dim(imputed.df)

# Remove string '_ReAMP' from the end of sample names
colnames(imputed.df) <- gsub(pattern = "_ReAMP", replacement = "", x = colnames(imputed.df))
head(colnames(imputed.df))
table(duplicated(colnames(imputed.df))) # any duplicates?


## Make sample names match for empirical data
#colnames(empirical.df)
head(colnames(empirical.df))

# Remove everything after the first underscore
colnames.df <- gsub(pattern = "_.*", replacement = "", x = colnames(empirical.df))
colnames.df <- as.data.frame(colnames.df)
colnames.df <- colnames.df[2:nrow(colnames.df), ] # drop mname for now
colnames.df <- as.data.frame(colnames.df)
colnames.df <- separate(data = colnames.df, col = "colnames.df", into = c("assay", "fam", "rep", "ind", "noise"), sep = "-", remove = T)
colnames.df$name <- paste0(colnames.df$assay, "_", colnames.df$fam, "_", colnames.df$rep, "_", colnames.df$ind)
colnames(empirical.df) <- c("mname", colnames.df$name)
rm(colnames.df)
table(duplicated(colnames(empirical.df))) # there are two cols that are duplicated! 
colnames(empirical.df)[duplicated(colnames(empirical.df))]
# Remove the duplicated columns
empirical.df <- empirical.df[, !duplicated(colnames(empirical.df))]

#### 03. Matching ####
## What columns match between the two
keep.cols <- intersect(x = colnames(imputed.df), y = colnames(empirical.df))
keep.cols
length(keep.cols)

# Keep only the keep cols from imputed
imputed.df <- imputed.df[, colnames(imputed.df) %in% keep.cols]
dim(imputed.df)

# Keep only the keep cols from empirical
empirical.df <- empirical.df[, colnames(empirical.df) %in% keep.cols]
dim(empirical.df)


#### 04. Clean up mname
imputed.df$mname <- gsub(pattern = " ", replacement = "__", x = imputed.df$mname)
imputed.df[1:5,1:5]
empirical.df$mname <- gsub(pattern = " ", replacement = "__", x = empirical.df$mname)
empirical.df[1:5,1:5]


#### 05. Add unique identifier to separate the two datasets
colnames(imputed.df) <- paste0(colnames(imputed.df), "_imputed")
colnames(empirical.df) <- paste0(colnames(empirical.df), "_empirical")


#### 06. Merge the two datasets
all_data.df <- merge(x = imputed.df, y = empirical.df, by.x = "mname_imputed", by.y = "mname_empirical")
dim(all_data.df)
all_data.df[1:5,1:5]
colnames(all_data.df) <- gsub("mname_imputed", replacement = "mname", x = colnames(all_data.df)) # clean up mname colname

# Sort by colname, then put mname first
all_data.df <- all_data.df[ , order(colnames(all_data.df))]
all_data.df <- all_data.df %>% 
  select("mname", everything())
all_data.df[1:5,1:5]


#### 07. Per chromosome, evaluate concordance ####
# Identify chr
chr <- unique(gsub(pattern = "__.*", replacement = "", x = all_data.df$mname))
chr

# Identify samples
samples.vec <- unique(gsub(pattern = "_empirical|_imputed", replacement = "", x = colnames(all_data.df)))
samples.vec <- samples.vec[grep(pattern = "mname", x = samples.vec, invert = T)]


# Per chromosome
subset_data.df <- NULL; result.df <- NULL; soi <- NULL ; score <- NULL
for(c in 1:length(chr)){
  
  print(paste0("Working on chr ", chr[c]))
  
  # Subset to the target chr
  subset_data.df <- all_data.df[grep(pattern = chr[c], x = all_data.df$mname), ]
  print(paste0("Number of records in this chr: ", nrow(subset_data.df)))
  
  # Loop to evaluate concordance per sample
  # Set up an empty df to fill
  result.df <- matrix(data = NA, nrow = length(samples.vec), ncol = 5)
  result.df <- as.data.frame(result.df)
  colnames(result.df) <- c("sample", "num.loci", "num.match", "num.missing", "prop.match")
  
  soi <- NULL ; score <- NULL
  for(i in 1:length(samples.vec)){
    
    soi <- samples.vec[i]
    score <- sum(subset_data.df[, paste0(soi, "_empirical")] == subset_data.df[, paste0(soi, "_imputed")])
    num_missing <- sum(subset_data.df[, paste0(soi, "_empirical")]==9) + sum(subset_data.df[, paste0(soi, "_imputed")]==9)
    prop_corr <- score / (nrow(subset_data.df) - num_missing)
    
    # Store results per sample
    result.df[i,"sample"] <- soi
    result.df[i,"num.loci"] <- nrow(subset_data.df)
    result.df[i,"num.match"] <- score
    result.df[i,"num.missing"] <- num_missing
    result.df[i,"prop.match"] <- prop_corr
    
    
  }
  
  # Summarize (printout)
  print(paste0("Mean ppn of typed loci concordant per sample: ", round(mean(result.df$prop.match, na.rm = T), digits = 3)))
  print(paste0("Mean number of typed loci concordant per sample: ", round(mean(result.df$num.match, na.rm = T), digits = 3)))
  
  # Write out
  write.table(x = result.df, file = paste0("13_impute_compare/concord_eval_", chr[c], "_comparison.txt")
              , sep = "\t", row.names = F
              )
  
  pdf(file = paste0("13_impute_compare/concord_eval_", chr[c], "_hist.pdf"), width = 6.5, height = 3.5)
  hist(x = result.df$prop.match, main = "", breaks = 10
       , xlab = paste0("Per sample proportion of typed loci concordant (", chr[c], ")")
       , las = 1
       )
  #text(x = 0.67, y = 30, paste0(nrow(subset_data.df), " loci"))
  dev.off()
              
}




# Write out all info
write.table(x = all_data.df, file = "13_impute_compare/all_loci_data_comparison.txt", sep = "\t", row.names = F)



# end #
