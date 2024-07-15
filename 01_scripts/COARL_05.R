# Post gemma run
# B. Sutherland (initialized 2024-05-15)

# Requires that 01_scripts/COARL_04.R was already run, and gemma was run
# Clear workspace, launch simple_pop_stats

library(ggplot2)

# User set variables
gemma_result_dir <- "03_results/gemma_run_dw_size_mean_2024-05-16_14h48/output"
gemma_result.FN <- "gwas_dw_size.assoc.txt"
pheno <- "dw_size_mean"

# Load gemma result
gemma_gwas <- read.table(paste0(gemma_result_dir, "/", gemma_result.FN), header = T)

# Sort by chr
gemma_gwas <- gemma_gwas[with(gemma_gwas, order(gemma_gwas$chr)), ]

head(gemma_gwas)

# Global variables
input_AF.FN <- "03_results/per_family_inferred_allele_frequency_data.RData"
date <- format(Sys.time(), "%Y-%m-%d_%Hh%M")

load(input_AF.FN)

## Plotting
# This contains the geno sim info
geno_sim[1:5, 1:5]
#geno_sim.bck <- geno_sim
colnames(geno_sim) <- gsub(pattern = "\\.1", replacement = "", x = colnames(geno_sim))
geno_sim[1:5, 1:5]

length(colnames(geno_sim) %in% gemma_gwas$rs)

# Obtain pheno data and give a good identifier
cross_and_pheno.df$family.id <- paste0("family_", cross_and_pheno.df$family)


#### Which are the top loci? ####
# Limit by chr
target_chr <- "Chr05"
#target_chr <- "Chr02"

gemma_gwas_subset.df <- gemma_gwas[gemma_gwas$chr==target_chr, ]
gemma_gwas_subset.df <- gemma_gwas_subset.df[with(gemma_gwas_subset.df, order(gemma_gwas_subset.df$p_wald, decreasing = F)), ]
gemma_gwas_subset.df <- gemma_gwas_subset.df[, c("chr", "rs", "ps", "af", "p_wald")]
head(gemma_gwas_subset.df, n = 25)

# What is the p-value threshold? 
print("To be considered significant, the p_wald is less than... ") 
print(0.05 / ncol(geno_sim))

write.table(x = gemma_gwas_subset.df, file = paste0("03_results/", "gwas_results_subset_", pheno, "_", target_chr, ".txt")
            , quote = F, sep = "\t", row.names = F, col.names = T
)

## Pick a locus
## Chr05
#plot_target <- "NC_047566_1_8411918"
#plot_target <- "NC_047566_1_8414668"
#plot_target <- "NC_047566_1_8796621"
#plot_target <- "NC_047566_1_8810784"
#plot_target <- "NC_047566_1_8181404"
#plot_target <- "NC_047566_1_7937836"

#plot_target <- "NC_047566_1_9593741"

## Pick a locus
## Chr02
#plot_target <- "NC_047564_1_46866720"
#plot_target <- "NC_047564_1_39587146"
#plot_target <- "NC_047564_1_41473163"
#plot_target <- "NC_047564_1_44629652"
#plot_target <- "NC_047564_1_47054066"
#plot_target <- "NC_047564_1_50899201"
#plot_target <- "NC_047564_1_58133796"

#plot_target <- "NC_047564_1_46680698"




# Obtain the genotypes of the target in the inferred dataset
plot_target.df <- geno_sim[, which(colnames(geno_sim) %in% plot_target)]
plot_target.df <- as.data.frame(plot_target.df)
colnames(plot_target.df) <- "locus_alt"
plot_target.df$family.id <- rownames(plot_target.df)
head(plot_target.df)

# Combine the genotype to the phenotype
head(cross_and_pheno.df)
plot.df <- merge(x = plot_target.df, y = cross_and_pheno.df, by = "family.id")
nrow(plot.df)==nrow(plot_target.df) # to make sure no rows were lost in the combining
head(plot.df)

# Add the AF
head(gemma_gwas_subset.df)
af <- gemma_gwas_subset.df[gemma_gwas_subset.df$rs==plot_target, "af"]

## base R version
# boxplot(plot.df$dw_size_mean ~ plot.df$locus_alt, main = paste0(target_chr, ", ", plot_target))
# points(x = plot.df$locus_alt, y = plot.df$dw_size_mean, col = "red")
#table(plot.df$locus_alt)

# Plot
p <- ggplot(data = plot.df, aes(x = locus_alt, y = dw_size_mean, group = locus_alt, color=locus_alt)) + geom_boxplot()
p <- p + theme(legend.position="none")
p <- p + ggtitle(paste0(target_chr, ", ", plot_target, ", AF=", af))
p
p_w_points <- p +  geom_jitter()
p_w_points

# Save out
pdf(file = paste0("03_results/", pheno, "_plot_alt_freq_", target_chr, "_", plot_target, ".pdf")
    , width = 7, height = 4.5)
print(p)
dev.off()

# Save out
pdf(file = paste0("03_results/", pheno, "_plot_alt_freq_", target_chr, "_", plot_target, "_w_points.pdf")
    , width = 7, height = 4.5)
print(p_w_points)
dev.off()






