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
<<<<<<< Updated upstream
input_folder <- "03_results/output_pheno_survival_state"
#input_folder <- "03_results/output_pheno_DPE"
input_folder <- "03_results_denovo_genos/output_pheno_survival_state_GR70_noimpute" # denovo variants

# Build output filename
gemma_output.FN  <- paste0(input_folder, "/gwas_all_fam_covar.assoc.txt")

=======
#gemma_output.FN       <- "03_results/output_pheno_survival_state/gwas_all_fam_covar.assoc.txt"
# gemma_output.FN       <- "03_results/output_pheno_DPE/gwas_all_fam_covar.assoc.txt"
#gemma_output.FN <- "03_results/family_sep_gemma_inputs_2024-10-11/F114/output/result.assoc.txt"
#gemma_output.FN <- "03_results/family_sep_gemma_inputs_2024-10-11/F115/output/result.assoc.txt"
#gemma_output.FN <- "03_results/family_sep_gemma_inputs_2024-10-11/F116/output/result.assoc.txt"
gemma_output.FN <- "03_results/family_sep_gemma_inputs_2024-10-11/F117/output/result.assoc.txt"
>>>>>>> Stashed changes

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


# Output plot filename
plot.FN <- gsub(pattern = "03_results\\/", replacement = "", x = gemma_output.FN)
plot.FN <- gsub(pattern = "\\/gwas_all_fam_covar.assoc.txt", replacement = "", x = plot.FN)
plot.FN <- gsub(pattern = "\\/result.assoc.txt", replacement = "", x = plot.FN)
plot.FN <- gsub(pattern = "\\/", replacement = "_", x = plot.FN)

# Plot histogram of pvals
pdf(file = paste0(plot.FN, "/pval_hist.pdf"), width = 7, height = 4.5)
hist(gemma_gwas$p_wald, breaks = 20)
dev.off()

# Create directory for output (if does not already exist)
dir.create(path = "03_results/GWAS_plots")

<<<<<<< Updated upstream
# Manhattan plot
output_folder <- gsub(pattern = "\\/.*", replacement = "", x = input_folder)
name_details <- gsub(pattern = ".*\\/", replacement = "", x = input_folder)
pdf(file = paste0(output_folder, "/GWAS_plots/Manhattan_", name_details, ".pdf"), width = 9, height = 3.5)
=======

pdf(file = paste0("03_results/GWAS_plots/Manhattan_", plot.FN, ".pdf"), width = 9, height = 5)
>>>>>>> Stashed changes
par(mfrow = c(1,1), mar = c(5,4,4,2) +0.1, mgp = c(3,1,0))
fastman(m = gemma_gwas, chr = "chr", bp = "pos.true", p = "p_wald"
        , genomewideline = -log10(0.05/nrow(gemma_gwas))
        , suggestiveline = -log10(0.05)
        , cex = 1, cex.lab = 1, cex.axis = 1
        #, ylim = c(0,10)
        )
dev.off()


##### End matter ####
save.image(file = paste0(input_folder, "/", "post_Manhattan_plotting.RData"))


# end
