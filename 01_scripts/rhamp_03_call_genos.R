# Calling genotypes
# B. Sutherland, Liam Surry, VIU (initialized 2024-01-11)

# Note: requires that 01_scripts/rhamp_02_data_checking.R has been run
#        and remains in environment

#### Call genotypes ####
# Based on dyes, call geno
for(j in 1:nrow(all_plates.df)){
  
  # Both dyes are NA, uncalled
  if(is.na(all_plates.df$Cq.fam[j]) & is.na(all_plates.df$Cq.vic.corr[j])){
    
    all_plates.df$geno[j] <- "no.geno"
    
  # VIC present & FAM = NA, genotype: homozygous alternate
  }else if(is.na(all_plates.df$Cq.fam[j]) & !is.na(all_plates.df$Cq.vic.corr[j])){
    
    all_plates.df$geno[j] <- "homo.alt"
    
  # FAM present & VIC = NA, genotype: homozygous reference 
  }else if(!is.na(all_plates.df$Cq.fam[j]) & is.na(all_plates.df$Cq.vic.corr[j])){
    
    all_plates.df$geno[j] <- "homo.ref"
    
  # VIC present & FAM present, genotype: heterozygous
  }else if(!is.na(all_plates.df$Cq.fam[j]) & !is.na(all_plates.df$Cq.vic.corr[j])){
    
    all_plates.df$geno[j] <- "het"
    
  }
  
}

# Genos
table(all_plates.df$geno)

write.table(x = table(all_plates.df$geno), file = "03_results/rhAmp_genos.txt"
            , quote = F, row.names = F, col.names = F
            )


head(all_plates.df)

# Next: go to 01_scripts/rhamp_04_well_to_sample.R
