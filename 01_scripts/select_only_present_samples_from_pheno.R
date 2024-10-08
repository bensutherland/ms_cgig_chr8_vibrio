setwd("~/Documents/cgig/CHR8_impute/simple_pop_stats_10X_GWAS")
present_samples <- read.delim("03_results/samples.txt", header = F)
head(present_samples)

pheno.df <- read.delim(file = "00_archive/G0923-21-VIUN_SampleInventory_V2_recd_2024-08-16.txt")
head(pheno.df)
pheno.df$indiv <- gsub(pattern = "_", replacement = "-", x = pheno.df$indiv)

pheno.df <- pheno.df[(pheno.df$indiv %in% present_samples$V1), ]
write.table(x = pheno.df, file = "03_results/present_sample_phenos.txt", quote = F, sep = "\t", col.names = T, row.names = F)
