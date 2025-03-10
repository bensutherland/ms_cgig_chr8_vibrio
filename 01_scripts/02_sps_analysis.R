# simple_pop_stats analysis component of the OCV23 RADseq analysis
# B. Sutherland
# Initialized 2023-10-23
# Requires first running "ms_cgig_chr8_vibrio/01_scripts/01_sps_char_and_filt.R"

# Clear space and source simple_pop_stats_start.R

# Load packages
#install.packages("ggpubr")
require("ggpubr")

#### 01. Load Data ####
load(file = "03_results/post_all_filters.RData") # loaded from prerequisite script above

# Data is present in
obj


#### 02. Multivariate analysis ####
# PCA from genind
pca_from_genind(data = obj
                , PCs_ret = 4
                , plot_eigen = TRUE
                , plot_allele_loadings = FALSE
                , retain_pca_obj = TRUE
                , colour_file = "00_archive/pop_cols.csv"
)

#### 0.3 Export ####
# Write out object
save.image(file = "03_results/post_all_filters_post_multivariate.RData")

# Go to "01_scripts/03_GWAS.R
