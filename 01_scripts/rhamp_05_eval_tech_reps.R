# Evaluate tech reps
# B. Sutherland, Liam Surry, VIU (initialized 2024-02-08)

# Note: requires that 01_scripts/rhamp_04_well_to_sample.R has been run

# Here is the current dataframe
head(data_annot.df)

# Evaluate technical reproducibility at the level of the genotype
data_annot.df <- data_annot.df[order(data_annot.df$DNA.id), ]
head(data_annot.df)
head(data_annot.df, n = 100)
View(data_annot.df)
