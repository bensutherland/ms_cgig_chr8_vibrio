# Plot multiple fastman results together
# B. Sutherland (VIU)
# 2024-10-22
library(fastman)

## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path) # take main directory
setwd(current.path)
rm(current.path)


# Set filenames
hotspot_only_data <- "~/Documents/cgig/CHR8_amp/ms_cgig_chr8/03_results/output_pheno_survival_state/post_Manhattan_plotting.RData"
denovo_panel_data <- "~/Documents/cgig/CHR8_amp/ms_cgig_chr8/03_results_denovo_genos/output_pheno_survival_state_GR70_noimpute/post_Manhattan_plotting.RData"
ai2_data          <- "~/Documents/cgig/CHR8_impute/impute_workflow_v.0.7/07_GWAS/ai2_imputed_survival_state_pheno/post-plot.RData"
fi3_data          <- "~/Documents/cgig/CHR8_impute/impute_workflow_v.0.7/07_GWAS/fi3_imputed_survival_state_pheno/post-plot.RData"

# Output plot
output.FN <- "03_results/multipanel_Manhattan_plot.pdf"

# Load hotspot data
load(hotspot_only_data)

# Save gemma output data
gemma_gwas_hotspot <- gemma_gwas
colnames(gemma_gwas_hotspot)[grep(pattern = "pos.true", x = colnames(gemma_gwas_hotspot))] <- "pos"
gemma_gwas_hotspot <- gemma_gwas_hotspot[grep(pattern = "^Chr", x = gemma_gwas_hotspot$chr), ] # only keep loci on chr

# Load denovo data
load(denovo_panel_data)

# Save gemma output data
gemma_gwas_denovo_panel <- gemma_gwas
colnames(gemma_gwas_denovo_panel)[grep(pattern = "pos.true", x = colnames(gemma_gwas_denovo_panel))] <- "pos"
gemma_gwas_denovo_panel <- gemma_gwas_denovo_panel[grep(pattern = "^Chr", x = gemma_gwas_denovo_panel$chr), ] # only keep loci on chr


# Load ai2 data
load(ai2_data)

# Save gemma output data
gemma_gwas_ai2 <- gemma_gwas


# Load fi3 data
load(fi3_data)

# Save gemma output data
gemma_gwas_fi3 <- gemma_gwas


#### Plot ####
pdf(file = output.FN, height = 7, width =  6)

par(mfrow = c(4,1), mar = c(5,6,2,2) +0.1, mgp = c(3,1,0))

# hotspot
fastman(m = gemma_gwas_hotspot, chr = "chr", bp = "pos", p = "p_wald"
        , genomewideline = -log10(0.05/nrow(gemma_gwas_hotspot))
        , suggestiveline = NULL
        , cex = 0.7, cex.lab = 1, cex.axis = 1
        #, ylim = c(0,10)
        , maxP = plot_maxP
)

mtext(text = "A", side = 2, line = 3
      , at = (max(-log10(gemma_gwas_hotspot$p_wald)) + ( 0.1  * max(-log10(gemma_gwas_hotspot$p_wald))))
      , las = 1)

# denovo panel
par(mar = c(5,6,2,2) +0.1, mgp = c(3,1,0))
fastman(m = gemma_gwas_denovo_panel, chr = "chr", bp = "pos", p = "p_wald"
        , genomewideline = -log10(0.05/nrow(gemma_gwas_denovo_panel))
        , suggestiveline = NULL
        , cex = 0.7, cex.lab = 1, cex.axis = 1
        #, ylim = c(0,10)
        , maxP = plot_maxP
        
)

mtext(text = "B", side = 2, line = 3
      , at = (max(-log10(gemma_gwas_denovo_panel$p_wald)) + ( 0.1  * max(-log10(gemma_gwas_denovo_panel$p_wald))))
      , las = 1)


# ai2
par(mar = c(5,6,2,2) +0.1, mgp = c(3,1,0))
fastman(m = gemma_gwas_ai2, chr = "chr", bp = "pos", p = "p_wald"
        , genomewideline = -log10(0.05/nrow(gemma_gwas_ai2))
        , suggestiveline = NULL
        , cex = 0.7, cex.lab = 1, cex.axis = 1
        #, ylim = c(0,10)
        , maxP = plot_maxP
        
)

mtext(text = "C", side = 2, line = 3
      , at = (max(-log10(gemma_gwas_ai2$p_wald)) + ( 0.1  * max(-log10(gemma_gwas_ai2$p_wald))))
      , las = 1)


# fi3
par(mar = c(5,6,2,2) +0.1, mgp = c(3,1,0))
fastman(m = gemma_gwas_fi3, chr = "chr", bp = "pos", p = "p_wald"
        , genomewideline = -log10(0.05/nrow(gemma_gwas_fi3))
        , suggestiveline = NULL
        , cex = 0.7, cex.lab = 1, cex.axis = 1
        #, ylim = c(0,10)
        , maxP = plot_maxP
        
)

mtext(text = "D", side = 2, line = 3
      , at = (max(-log10(gemma_gwas_fi3$p_wald)) + ( 0.1  * max(-log10(gemma_gwas_fi3$p_wald))))
      , las = 1)


dev.off()

