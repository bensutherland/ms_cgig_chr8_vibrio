# Calling genotypes
# B. Sutherland, Liam Surry, VIU (initialized 2024-01-11)

# Note: requires that 01_scripts/rhamp_03_call_genos.R has been run
#        and remains in environment

# Set user variables
sample_info.FN <- ""


head(all_plates.df)

# Reduce the cols
data.df  <-    all_plates.df[, 
                             c("Well", "Cq.fam", "Cq.fam.corr", "Cq.vic", "Cq.vic.corr", "diff", "geno", "full.id")]

data.df <- as.data.frame(data.df)
head(data.df)

# # Separate full ID into components
# data.df    <-  separate(data = data.df, col = "full.id"
#                      , into = c("plate.id", "well.id")
#                      , sep = "__", remove = T)
# 
# head(data.df)



