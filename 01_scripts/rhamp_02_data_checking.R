# Data analysis of formatted rhAmp data
# B. Sutherland, Liam Surry, Ally Gignac VIU (initialized 2024-01-11)

# Note: requires that 01_scripts/rhamp_01_read_in_and_format.R has been run
#        and remains in environment

# Note: currently assumes a set cutoff for FP designation (not plate-specific)
# Note: no correction applied to homozygotes (# TODO: set hard cutoff as well (?))

# Set variables
# FP cutoffs
FP_cutoff_fam.val <-  (-1)
FP_cutoff_vic.val <-  (-15)


#### Join data from all plates ####
# View data
names(data.list)
head(data.list[[1]])

# Combine all plates by column names
all_plates.df <- dplyr::bind_rows(data.list)
dim(all_plates.df)
print(paste0("Total number of unique wells in dataset: ", length(unique(all_plates.df$full.id))))


#### Missing data check ####
# How much missing data? (i.e., no call for both)
print("Missing values (including neg. controls): ")
table(is.na(all_plates.df$Cq.fam) & is.na(all_plates.df$Cq.vic))

# Save out complete missing data samples
noCq_wells.df  <- all_plates.df[is.na(all_plates.df$Cq.fam) & is.na(all_plates.df$Cq.vic), ]
write.csv(x = noCq_wells.df, file = "03_results/rhAmp_no_Cq_wells.csv", row.names = F)


#### Correct false-positive heterozygotes (both dyes) ####
# Set variables
all_plates.df$Cq.vic.corr <- NA
all_plates.df$Cq.fam.corr <- NA

# Correct Cq vals outside of FP het cutoffs
corrected <- NULL
for(j in 1:nrow(all_plates.df)){
  
  # If the call has both VIC and FAM (i.e., a 'diff' value), consider for correction
  if(!is.na(all_plates.df$diff[j])){
    
    # VIC: remove values outside the cutoff
    if(all_plates.df$diff[j] < FP_cutoff_vic.val){
      
      print(paste0("[",j,"]", "- VIC false positive detected, leaving VIC value as NA"))
      print(paste0("diff: ", all_plates.df$diff[j]))
      
      # FAM: keep
      all_plates.df$Cq.fam.corr[j] <- all_plates.df$Cq.fam[j]
      
      # Retain name of corrected well
      corrected <- c(corrected, paste0("VIC_FP: ", all_plates.df[j, "Well"]))
      
    # FAM: remove values outside the cutoff
    }else if(all_plates.df$diff[j] > FP_cutoff_fam.val){
      
      print(paste0("[",j,"]", "- FAM false positive detected, leaving FAM value as NA"))
      print(paste0("diff: ", all_plates.df$diff[j]))
      
      # VIC: retain
      all_plates.df$Cq.vic.corr[j] <- all_plates.df$Cq.vic[j]
      
      # Retain name of corrected well
      corrected <- c(corrected, paste0("FAM_FP: ", all_plates.df[j, "Well"]))
      
    # Genotype call difference within range: keep heterozygote
    }else{
      
      # Retain VIC
      all_plates.df$Cq.vic.corr[j] <- all_plates.df$Cq.vic[j]
      
      # Retain FAM
      all_plates.df$Cq.fam.corr[j] <- all_plates.df$Cq.fam[j]
    }
    
    # If the diff is NA, we can't remove potential FPs, so keep both values as is
    
    
  }else{
    
    all_plates.df$Cq.vic.corr[j] <- all_plates.df$Cq.vic[j]
    all_plates.df$Cq.fam.corr[j] <- all_plates.df$Cq.fam[j]
    
  }
  
}

write.table(x = corrected, file = "03_results/rhAmp_corrected_wells.txt"
            , row.names = F, col.names = F, quote = F
            )

#head(all_plates.df, n = 20)

# Next: go to 01_scripts/rhamp_03_call_genos.R
