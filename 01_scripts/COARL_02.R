# Infer offspring genotypes based on cross information
# B. Sutherland (2024-02-11)

# Requires that 01_scripts/COARL_02.R was already run, and still in environment

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

# Create loci dataframe
loci <- locNames(obj)
loci.df <- as.data.frame(loci)
rm(loci)
head(loci.df)


family <- NULL; dam <- NULL; loci.df <- NULL
# full: nrow(cross_and_pheno.df)
for(i in 1:2){
  
  # family
  family <- cross_and_pheno.df[i, "Family"]
  
  # dam
  dam  <- cross_and_pheno.df[i, "Dam"]
  
  # sire
  sire <- cross_and_pheno.df[i, "Sire"]
  
  loci.df <- cbind(loci.df, rep(NA, times = nrow(loci.df)))
  
  colnames(loci.df)[i+1] <- paste0("family_", family, "_a2_count")
  head(loci.df)
  
  num.alt <- NULL; mname <- NULL
  for(m in 1:nrow(loci.df)){
    
    mname <- paste0(loci.df[m, "loci"], ".1")
    
    # pull the marker number of alt alleles for dam and sire combined
    loci.df[m,i+1] <- geno.df[dam, mname] + 
                           geno.df[sire, mname]
    
  }
}


loci.df[1:10, ]

# Next: convert from number of alt alleles to AF
unique(loci.df[,"family_1_a2_count"])

# Finally, GWAS
