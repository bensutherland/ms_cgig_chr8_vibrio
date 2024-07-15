# Calling genotypes
# B. Sutherland, Liam Surry, VIU (initialized 2024-01-11)

# Note: requires that 01_scripts/rhamp_02_data_checking.R has been run
#        and remains in environment

#### Call genotypes ####
# Based on dyes, call geno
for(j in 1:nrow(all_plates.df)){
  
  # If both dyes are NA, set value as uncalled
  if(is.na(all_plates.df$Cq.fam.corr[j]) & is.na(all_plates.df$Cq.vic.corr[j])){
    
    all_plates.df$geno[j] <- "no.geno"
    
  # If VIC is present & FAM is NA, set genotype: homozygous alternate
  }else if(is.na(all_plates.df$Cq.fam.corr[j]) & !is.na(all_plates.df$Cq.vic.corr[j])){
    
    all_plates.df$geno[j] <- "homo.alt"
    
  # If FAM is present & VIC is NA, set genotype: homozygous reference 
  }else if(!is.na(all_plates.df$Cq.fam.corr[j]) & is.na(all_plates.df$Cq.vic.corr[j])){
    
    all_plates.df$geno[j] <- "homo.ref"
    
  # If VIC is present & FAM is present, set genotype: heterozygous
  }else if(!is.na(all_plates.df$Cq.fam.corr[j]) & !is.na(all_plates.df$Cq.vic.corr[j])){
    
    all_plates.df$geno[j] <- "het"
    
  }
  
}

# Genos
table(all_plates.df$geno)

write.table(x = all_plates.df, file = "03_results/rhAmp_genos.txt"
            , quote = F, row.names = F, col.names = F
            )


head(all_plates.df)

save.image(file = "03_results/prepared_rhamp_results.RData")

# Next: go to 01_scripts/rhamp_04_well_to_sample.R
