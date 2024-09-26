# Analyze gemma output 
#  requires GEMMA run in ms_cgig_chr8/03_results as per README.md
#  initialized 2024-06-17
#  Ben J. G. Sutherland (VIU), incl. code dev by Konstantin Divilov

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

# Load libraries
library(vcfR)
#devtools::install_github('kaustubhad/fastman',build_vignettes = TRUE)
library(fastman)
#install.packages("missMethods")
library(missMethods)
library(tidyr)
library(ggplot2)

## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path) # take main directory
setwd(current.path)
rm(current.path)

# User set variables
#gemma_output.FN       <- "03_results/output_pheno_survival_state/gwas_all_fam_covar.assoc.txt"
gemma_output.FN       <- "03_results/output_pheno_DPE/gwas_all_fam_covar.assoc.txt"

#### 01. Load GEMMA results ####
# Read in GEMMA output
gemma_gwas <- read.table(file = gemma_output.FN, header = T)
gemma_gwas <- as.data.frame(gemma_gwas)
head(gemma_gwas)

# Separate marker name into chr and pos
gemma_gwas <- separate(data = gemma_gwas, col = "rs", into = c("chr.true", "pos.true"), sep = "__", remove = F)
gemma_gwas$pos.true <- as.numeric(gemma_gwas$pos.true)

# convert to linkage group (LG), based on https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_902806645.1/
# LG to Chr info can be found in Sup File in https://doi.org/10.1093/gigascience/giab020
gemma_gwas$chr <- gemma_gwas$chr.true
gemma_gwas$chr <- gsub(pattern = "NC_047559.1", replacement = "Chr07", x = gemma_gwas$chr)
gemma_gwas$chr <- gsub(pattern = "NC_047560.1", replacement = "Chr01", x = gemma_gwas$chr)
gemma_gwas$chr <- gsub(pattern = "NC_047561.1", replacement = "Chr09", x = gemma_gwas$chr)
gemma_gwas$chr <- gsub(pattern = "NC_047562.1", replacement = "Chr06", x = gemma_gwas$chr)
gemma_gwas$chr <- gsub(pattern = "NC_047563.1", replacement = "Chr03", x = gemma_gwas$chr)
gemma_gwas$chr <- gsub(pattern = "NC_047564.1", replacement = "Chr02", x = gemma_gwas$chr)
gemma_gwas$chr <- gsub(pattern = "NC_047565.1", replacement = "Chr04", x = gemma_gwas$chr)
gemma_gwas$chr <- gsub(pattern = "NC_047566.1", replacement = "Chr05", x = gemma_gwas$chr)
gemma_gwas$chr <- gsub(pattern = "NC_047567.1", replacement = "Chr10", x = gemma_gwas$chr)
gemma_gwas$chr <- gsub(pattern = "NC_047568.1", replacement = "Chr08", x = gemma_gwas$chr)

gemma_gwas <- gemma_gwas[with(gemma_gwas, order(gemma_gwas$chr)), ]

hist(gemma_gwas$p_wald, breaks = 20)


# Output plot filename
plot.FN <- gsub(pattern = "03_results\\/", replacement = "", x = gemma_output.FN)
plot.FN <- gsub(pattern = "\\/gwas_all_fam_covar.assoc.txt", replacement = "", x = plot.FN)

pdf(file = paste0("03_results/GWAS_plots/Manhattan_", plot.FN, ".pdf"), width = 9, height = 5)
par(mfrow = c(1,1), mar = c(5,4,4,2) +0.1, mgp = c(3,1,0))
fastman(m = gemma_gwas, chr = "chr", bp = "pos.true", p = "p_wald"
        , genomewideline = -log10(0.05/nrow(gemma_gwas))
        , suggestiveline = -log10(0.05)
        , cex = 1, cex.lab = 1, cex.axis = 1
        #, ylim = c(0,10)
        )
dev.off()

# end

