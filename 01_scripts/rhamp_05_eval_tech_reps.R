# Evaluate tech reps
# B. Sutherland, Liam Surry, VIU (initialized 2024-02-08)

# Note: requires that 01_scripts/rhamp_04_well_to_sample.R has been run

# Here is the current dataframe
head(data_annot.df)

# Evaluate technical reproducibility at the level of the genotype
data_annot.df <- data_annot.df[order(data_annot.df$sample), ]
head(data_annot.df)
head(data_annot.df, n = 100)

samples <- unique(data_annot.df$sample)
length(samples)

# Remove controls
samples <- samples[which(samples!="Neg")]

# Remove no genos
dim(data_annot.df)
data_annot.df.bck <- data_annot.df

# Remove rows without a geno call
data_annot.df <- data_annot.df[data_annot.df$geno!="no.geno", ]
dim(data_annot.df)

# Make a new dataframe
results.df <- as.data.frame(samples)
head(results.df)
rm(samples)

# Make an empty column to fill with the most frequently viewed genotype
results.df$majority.geno <- rep(NA, times = nrow(results.df))

# Make an empty column to fill with the percentage of calls as the most frequent geno
results.df$ppn.corr <- rep(NA, times = nrow(results.df))

# Make an empty column to fill with the number of total calls (with a genotype)
results.df$num.calls <- rep(NA, times = nrow(results.df))

head(results.df)

# Loop over samples
sample <- NULL; slice <- NULL; majority.geno <- NULL
for(i in 1:nrow(results.df)){
  
  # Sample of interest
  sample <- results.df[i,"samples"]
  
  # Reporting
  print(paste0("Working on sample ", sample))
  
  # Take the data for the target sample
  slice <- data_annot.df[data_annot.df$sample==sample, ]
  
  # Define the majority geno
  majority.geno <- names(head(sort(table(slice$geno), decreasing = TRUE), n = 1))
  
  # Fill in the value for the majority geno
  results.df[results.df$samples==sample, "majority.geno"] <- majority.geno

  # Fill in the value for the proportion correct
  results.df[results.df$samples==sample, "ppn.corr"]      <- round(x = nrow(slice[slice$geno==majority.geno, ]) / nrow(slice), digits = 3)
  
  # Fill in the value for the number calls
  results.df[results.df$samples==sample, "num.calls"]   <- nrow(slice)

  
}


results.df

write.table(x = results.df, file = "03_results/rhAmp_per_sample_genotype_summary.txt", quote = F, sep = "\t", row.names = F)

