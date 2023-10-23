# simple_pop_stats analysis component of the Yesso scallop RADseq analysis
# B. Sutherland
# Initialized 2023-10-23
# Requires first running "ms_cgig_chr8/01_scripts/02_sps_char_and_filt.R"



#install.packages("ggpubr")
require("ggpubr")

#### 01. Load Data ####
load(file = "03_results/post_all_filters.RData") # loaded from prerequisite script above

# Data is present in
obj



## Multivariate
# For an unknown reason, sps currently requires a
#   manual sourcing of the function to properly use the retain_pca_obj function
#source("~/Documents/pyes/simple_pop_stats/01_scripts/utilities/pca_from_genind.r", echo=TRUE)

# PCA from genind
pca_from_genind(data = obj
                , PCs_ret = 4
                , plot_eigen = TRUE
                , plot_allele_loadings = TRUE
                , retain_pca_obj = TRUE
                , colour_file = "00_archive/formatted_cols.csv"
)

file.copy(from = "03_results/pca_scores_per_sample.txt", to = "03_results/pca_scores_per_sample_sibs_incl.txt")

## Prepare an eigenvalue plot for inset
num_eigenvals <- 10
vals.df <- as.data.frame(pca.obj$eig[1:num_eigenvals])
colnames(vals.df)[1] <- "vals"
vals.df$pc <- seq(1:num_eigenvals)
vals.df
colnames(vals.df) <- c("PVE", "PC")

# Express eigenvalues as a percentage of total variation explained
tot.var <- sum(pca.obj$eig)
vals.df$PVE <- vals.df$PVE/tot.var *100

# Barplot
eig.plot <- ggplot(data = vals.df, aes(x=PC, y=PVE)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x=element_blank() #remove x axis labels
        , axis.ticks.x=element_blank() #remove x axis ticks
        , panel.background = element_blank()
  )


## Plot
# Remove legend pc1 v pc2
pc1_v_pc2.plot  <- pc1_v_pc2.plot + theme(legend.position = "none")
pc1_v_pc2.plot  <- pc1_v_pc2.plot + annotation_custom(ggplotGrob(eig.plot)
                                                      , xmin = 1, xmax = 5
                                                      , ymin = -10, ymax = -3.5
)

# Legend inside panel second plot
pc3_v_pc4.plot <- pc3_v_pc4.plot + theme(legend.justification = c(1,0), legend.position = c(1,0)
                                         , legend.background = element_rect(colour = "black", fill = "white", linetype = "solid")
)

final.figure <- ggarrange(pc1_v_pc2.plot, pc3_v_pc4.plot
                          , labels = c("A", "B")
                          , ncol = 2, nrow = 1
)

pdf(file = "03_results/pca_composite_figure_w_close_kin.pdf", width = 12, height = 6.5)
print(final.figure)
dev.off()


# STILL TO DO: 
# Make separate populations as separate datasets
# FST comparison between dead/alive
# Plot in Manhattan plot
# Possibly statistic test for p-value



# single SNP per locus analysis is complete
