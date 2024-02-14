# Infer offspring genotypes based on cross information
# B. Sutherland (2024-02-11)

# Requires that 01_scripts/COARL_02.R was already run
# Clear workspace, launch simple_pop_stats

load("03_results/prepared_data.RData")

#### 04. Bring in cross information, per cross, per locus, infer offspring AF ####
# read in crosses info
cross_and_pheno.df <- read.table(file = cross.FN, header = T, sep = "\t")
head(cross_and_pheno.df)
dim(cross_and_pheno.df)

# match to the genos
cross_and_pheno.df$Dam <- str_pad(string = cross_and_pheno.df$Dam, width = 2, side = "left", pad = "0")
cross_and_pheno.df$Sire <- str_pad(string = cross_and_pheno.df$Sire, width = 2, side = "left", pad = "0")

cross_and_pheno.df$Dam  <- paste0("female_", cross_and_pheno.df$Dam)
cross_and_pheno.df$Sire <- paste0("male_", cross_and_pheno.df$Sire)

head(cross_and_pheno.df)

# Create genotype dataframe
geno.df <- obj$tab
geno.df[1:5,1:5]

### Infer offspring allele frequencies for each family in the cross
# Create locus dataframe
loci <- locNames(obj)
loci.df <- as.data.frame(loci)
rm(loci)
head(loci.df)

# Set nulls
family <- NULL; dam <- NULL; column_of_interest <- NULL

for(i in 1:nrow(cross_and_pheno.df)){
  
  print(i)
  
  # family
  family <- cross_and_pheno.df[i, "Family"]
  
  # dam
  dam  <- cross_and_pheno.df[i, "Dam"]
  
  # sire
  sire <- cross_and_pheno.df[i, "Sire"]
  
  # Reporting
  print(paste0("Family ", family, " is comprised of ", dam, " and ", sire))
  
  # Create an empty column for inferred numbers of alt alleles per locus
  #loci.df <- cbind(loci.df, rep(NA, times = nrow(loci.df)))
  
  # Set the column name for the family
  column_of_interest <- paste0("family_", family, "_a2_ppn")
  colnames(loci.df)[i+1] <- column_of_interest
  
  # Reporting
  print("Starting to fill estimated proportions")
  print(head(loci.df))
  
  # Set nulls
  num.alt <- NULL; mname <- NULL
  
  for(m in 1:nrow(loci.df)){
    
    # Marker name for this iteration
    mname <- paste0(loci.df[m, "loci"], ".1")
    
    # Extract the number of alt alleles for each parent (order does not matter)
    loci.df[m, column_of_interest] <- paste0(geno.df[dam, mname], "__", geno.df[sire, mname])
    
  }
  
  # What are the unique genotypes?
  unique(loci.df[,column_of_interest])
  
  ## Fix the proportion data
  # if any NA, cannot estimate
  loci.df[grep(pattern = "NA", x = loci.df[, column_of_interest]), column_of_interest] <- "NA"
  
  # What are the unique genotypes, now that the NA-containing ones are removed?
  unique(loci.df[,column_of_interest])
  
  ### Inference logic ###
  # No alt alleles in either parent
  loci.df[, column_of_interest] <- gsub(pattern = "0__0", replacement = "0", x = loci.df[, column_of_interest])
  
  # One parent has one alt allele
  loci.df[, column_of_interest] <- gsub(pattern = "0__1|1__0", replacement = "0.25", x = loci.df[, column_of_interest])
  
  # Both parents have one alt allele
  loci.df[, column_of_interest] <- gsub(pattern = "1__1", replacement = "0.5", x = loci.df[, column_of_interest])
  
  # One parent has two alt allele, one parent has no alt alleles
  loci.df[, column_of_interest] <- gsub(pattern = "0__2|2__0", replacement = "0.5", x = loci.df[, column_of_interest])
  
  # One parent has two alt alleles and one parent has one alt allele
  loci.df[, column_of_interest] <- gsub(pattern = "2__1|1__2", replacement = "0.75", x = loci.df[, column_of_interest])
  
  # Both parents have two alt alleles
  loci.df[, column_of_interest] <- gsub(pattern = "2__2", replacement = "1", x = loci.df[, column_of_interest])
  
}


loci.df[1:10, ]

save.image(file = "03_results/per_family_inferred_allele_frequency_data.RData")

# Next step will be to conduct the GWAS
