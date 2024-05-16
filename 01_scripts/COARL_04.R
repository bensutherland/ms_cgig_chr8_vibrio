# Post gemma run
# B. Sutherland (initialized 2024-05-15)

# Requires that 01_scripts/COARL_03.R was already run, and gemma was run
# Clear workspace, launch simple_pop_stats

# Load libraries
library(missMethods)
library(fastman)


# User set variables
gemma_result_dir <- "03_results/gemma_run_dw_size_mean_2024-05-16_14h48/output"
gemma_result.FN <- "gwas_dw_size.assoc.txt"
pheno <- "dw_size_mean"

# gemma_result_dir <- "03_results/gemma_run_dw_per_d_tot_2024-05-15_11h00/output"
# gemma_result.FN <- "gwas_dw_per_d_tot.assoc.txt"

# gemma_result_dir <- "03_results/gemma_run_dw_per_d_perf_2024-05-15_11h21/output"
# gemma_result.FN <- "gemma_gwas_prepared.assoc.txt"

# gemma_result_dir <- "03_results/gemma_run_sw_per_d_perf_2024-05-15_11h30/output"
# gemma_result.FN <- "gemma_gwas_prepared.assoc.txt"
# pheno <- "sw_size_mean"

# gemma_result_dir <- "03_results/gemma_run_dw_minus_sw_size_mean_2024-05-15_11h39/output"
# gemma_result.FN <- "gemma_gwas_prepared.assoc.txt"
# pheno <- "dw_minus_sw_size_mean"

# gemma_result_dir <- "03_results/gemma_run_sw_size_mean_2024-05-15_11h54/output"
# gemma_result.FN <- "gemma_gwas_prepared.assoc.txt"
# pheno <- "sw_size_mean"


p_val <- "p_wald"
# p_val <- "p_score"
#p_val <- "p_lrt"


# Load gemma result
gemma_gwas = read.table(paste0(gemma_result_dir, "/", gemma_result.FN), header = T)
head(gemma_gwas)

# Sort by chr
gemma_gwas <- gemma_gwas[with(gemma_gwas, order(gemma_gwas$chr)), ]

# Only keep chr? 


# Plot pval distributions
pdf(file = paste0(gemma_result_dir, "/", "pval_histogram.pdf"), width = 10.5, height = 6.5)
par(mfrow=c(3,1))
hist(x = gemma_gwas$p_wald, breaks = 200, main = "p_wald")
hist(x = gemma_gwas$p_score, breaks = 200, main = "p_score")
hist(x = gemma_gwas$p_lrt, breaks = 200, main = "p_lrt")
dev.off()

ylim_val <- ((-log10(0.05/nrow(gemma_gwas))) + 1)

pdf(file = paste0(gemma_result_dir, "/", "Manhattan_plot.pdf"), width = 10.5, height = 6.5)
par(mfrow = c(1,1), mar = c(5,4,4,2)+0.1, mgp = c(3,1,0))
fastman(gemma_gwas
        , chr = "chr"
        , bp = "ps"
        , p=p_val
        , genomewideline = -log10(0.05/nrow(gemma_gwas))
        , suggestiveline = NULL
        , cex=1
        , cex.lab=1
        , cex.axis=1
        , main=pheno
        , xlab = ""
        , ylim = c(0,ylim_val)
        )
dev.off()

# Compare p-vals
par(mfrow=c(2,1))
plot(x = gemma_gwas$p_wald, y = gemma_gwas$p_lrt)
plot(x = gemma_gwas$p_wald, y = gemma_gwas$p_score)

