# Read in AlphaImpute output file
# B. Sutherland (2024-07-15)

### Front Matter ####
# Clean space
# rm(list=ls())

## Install and load packages
#install.packages("rstudioapi")
#install.packages("vcfR")
#install.packages("data.table")
library("rstudioapi")
library("vcfR")
library("adegenet")
library("data.table")

# Set working directory to the ms_scallop_popgen repo
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

## Info
# sessionInfo()

# Set variables
input.FN <- "ai2test.genotypes"

imputed.mat <- fread(file = input.FN, sep = " ", header = F, stringsAsFactors = F)

dim(imputed.mat)
str(imputed.mat)
imputed.mat[1:5,1:5]
imputed.df <- as.data.frame(imputed.mat)
imputed.df[1:5,1:5]

imputed_t.df <- t(imputed.df)
dim(imputed_t.df)
imputed_t.df[1:5,1:5]

colnames(imputed_t.df) <- imputed_t.df[1,] # use the first row as colnames
imputed_t.df <- imputed_t.df[2:nrow(imputed_t.df),]        # then drop the first row, leaving only genotypes
imputed_t.df[1:5,1:5]
imputed_t.mat <- as.matrix(imputed_t.df)
imputed_t.mat[1:5,1:5]

# Convert back to standard format
imputed_t.mat <- gsub(pattern = "\\b0\\b", replacement = "0/0", x = imputed_t.mat)
imputed_t.mat <- gsub(pattern = "\\b1\\b", replacement = "0/1", x = imputed_t.mat)
imputed_t.mat <- gsub(pattern = "\\b2\\b", replacement = "1/1", x = imputed_t.mat)

# info.df <- rep(NA, times = nrow(imputed_t.mat))
# info.df <- as.data.frame(info.df)
# colnames(info.df) <- "FORMAT"
# head(info.df)
# 
# imputed_t_w_fmt.mat <- cbind(info.df, imputed_t.mat)
# dim(imputed_t_w_fmt.mat)
# imputed_t_w_fmt.mat[1:5,1:5]

#unique(imputed_t.mat[,which(colnames(imputed_t.mat)=="55-41F")])

# Read in the contributing VCF file
input_vcf.FN <- "04_impute_panel/wgrs_filtered_parent_loci_amp_panel_parent_loci_amp_panel_offspring_loci_NC_047567.1.vcf"

#### 01. Import VCF ####
input.vcf <- vcfR::read.vcfR(file = input_vcf.FN)
input.vcf

dim(input.vcf@gt)
input.vcf@gt[1:5,1:5]
dim(input.vcf@gt)
input.vcf@fix[1:5,]

geno_info.df <- input.vcf@fix
geno_info.df <- as.data.frame(geno_info.df)

mname <- paste0(geno_info.df$CHROM, "__", geno_info.df$POS, "__", geno_info.df$ID)
length(mname)

imputed_t.mat[1:5,1:5]
rownames(imputed_t.mat) <- mname

imputed_t.mat[1:5,1:5]

data.mat <- t(imputed_t.mat)
data.mat[1:5,1:5]
colnames(data.mat) <- gsub(pattern = "\\.", replacement = "_", x = colnames(data.mat))

data.gid <- df2genind(X = data.mat, sep = "/")

# source sps fn
# Set working directory to the ms_scallop_popgen repo
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

pca_from_genind(data = data.gid, PCs_ret = 4, plot_eigen = F, plot_allele_loadings = F
                , plot_ellipse = F, retain_pca_obj = T
                )


# # Replace the gt 
# dim(imputed_t_w_fmt.mat)
# imputed_t_w_fmt.mat <- as.matrix(imputed_t_w_fmt.mat)
# 
# input.vcf@gt <- imputed_t_w_fmt.mat
# 
# data.gl <- vcfR2genlight(x = input.vcf)
# data.gl$gen
# 
# pca1 <- glPca(x = data.gl, nf = 4, parallel = T)
# 
# 
# 
