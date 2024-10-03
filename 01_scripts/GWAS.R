# Prepare inputs for a GWAS from stacks_workflow to GEMMA
# K. Divilov, edits by B. Sutherland, L. Surry
# initialized 2023-12-29

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

## Install packages
# install.packages("vcfR")
# install.packages("fastman")
# install.packages("missMethods")
# install.packages("rstudioapi")
# install.packages("tidyr")
# install.packages("norm")

# Load libraries
library(vcfR)
library(fastman)
library(missMethods)
library(rstudioapi)
library(tidyr)
library(norm)

# Set working directory to the ms_scallop_popgen repo
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

# User-set variables
VCF.FN <- "02_input_data/populations.snps_single-SNP_per_tag_2023-10-23.vcf"
map.FN <- "02_input_data/populations.plink_2023-10-23.map"
interp.FN <- "02_input_data/sample_interp_2024-07-18.csv"

impute_type <- "EM" # either "mean" or "EM" (i.e., expectation-maximization)


#### 01. Read in data ####
# Read in genotype data
vcf = read.vcfR(file = VCF.FN)
vcf

# Read in map data
map = read.table(file = map.FN, header = F)
#colnames(map) <- c("scaffold", "mname", "unkn", "pos")
head(map)


#### 02. Prepare marker names and CHR info ####
# Add linkage group (LG) info to map file, based on https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_902806645.1/
map$LG = NA
map$LG[map$V1=="NC_047559.1"] = 1
map$LG[map$V1=="NC_047560.1"] = 2
map$LG[map$V1=="NC_047561.1"] = 3
map$LG[map$V1=="NC_047562.1"] = 4
map$LG[map$V1=="NC_047563.1"] = 5
map$LG[map$V1=="NC_047564.1"] = 6
map$LG[map$V1=="NC_047565.1"] = 7
map$LG[map$V1=="NC_047566.1"] = 8
map$LG[map$V1=="NC_047567.1"] = 9
map$LG[map$V1=="NC_047568.1"] = 10

# Convert LG to chr (info can be found in Sup File in https://doi.org/10.1093/gigascience/giab020 )
map$Chr = NA
map$Chr[map$LG==1] = "Chr7"
map$Chr[map$LG==2] = "Chr1"
map$Chr[map$LG==3] = "Chr9"
map$Chr[map$LG==4] = "Chr6"
map$Chr[map$LG==5] = "Chr3"
map$Chr[map$LG==6] = "Chr2"
map$Chr[map$LG==7] = "Chr4"
map$Chr[map$LG==8] = "Chr5"
map$Chr[map$LG==9] = "Chr10"
map$Chr[map$LG==10] = "Chr8"

head(map)

# Create marker name (chr_pos)
map$SNPname = paste(map$Chr, map$V4, sep="_")
head(map)


#### 03. Prepare genotypes ####
# Extract genotypes from vcf
geno = extract.gt(vcf, element = "GT")
geno[1:5, 1:20]

# Change SNP names to "Chr[1-10]_[location in bp]"
rownames(geno) = map$SNPname

# Remove SNPs on contigs, i.e., not on Chr1-10
dim(geno)
geno = geno[-grep("NA_", rownames(geno)), ]
dim(geno)

# missingness sanity check
plot(rowSums(is.na(geno))/ncol(geno)) # missing data by locus
plot(colSums(is.na(geno))/nrow(geno)) # missing data by sample

# Change genotypes to numeric values
geno[1:5, 1:5]
geno[geno=="0/0"] = 0
geno[geno=="0/1"] = 1
geno[geno=="1/1"] = 2

# Transpose
geno = t(geno)

# Convert to numeric
mode(geno) = "numeric"

# Remove parental samples
geno = geno[-grep("F0",rownames(geno)),]
dim(geno)

# Create family-specific genotype matrices
geno_F114 = geno[grep("F114",rownames(geno)),]
geno_F115 = geno[grep("F115",rownames(geno)),]
geno_F116 = geno[grep("F116",rownames(geno)),]
geno_F117 = geno[grep("F117",rownames(geno)),]
geno_OFR6.10 = geno[grep("OFR6.10",rownames(geno)),]

# Family-specific imputation
if(impute_type == "mean"){
  
  print("Using imputation type 'mean'")
  
  # Run family-specific mean imputation on genotypes
  geno_F114_impute = impute_mean(geno_F114)
  geno_F115_impute = impute_mean(geno_F115)
  geno_F116_impute = impute_mean(geno_F116)
  geno_F117_impute = impute_mean(geno_F117)
  geno_OFR6.10_impute = impute_mean(geno_OFR6.10)
  
}else if(impute_type == "EM"){
  
  print("Using imputation type 'EM'")
  
  # Run family-specific EM imputation on genotypes
  geno_F114_impute = impute_EM(geno_F114)
  geno_F115_impute = impute_EM(geno_F115)
  geno_F116_impute = impute_EM(geno_F116)
  geno_F117_impute = impute_EM(geno_F117)
  geno_OFR6.10_impute = impute_EM(geno_OFR6.10)
  
}


# Reconstruct full genotype matrix
geno = rbind(geno_F114_impute,
             geno_F115_impute,
             geno_F116_impute,
             geno_F117_impute,
             geno_OFR6.10_impute)


#### 04. Prepare GEMMA inputs ####
#create variable holding family information
var_family = rep(NA,nrow(geno))
var_family[grep("F114",rownames(geno))] = "F114"
var_family[grep("F115",rownames(geno))] = "F115"
var_family[grep("F116",rownames(geno))] = "F116"
var_family[grep("F117",rownames(geno))] = "F117"
var_family[grep("OFR6.10",rownames(geno))] = "OFR6.10"

#create variable holding live/dead information
var_status = rep(NA,nrow(geno))
var_status[grep("Alive",rownames(geno))] = "1"
var_status[grep("Dead",rownames(geno))] = "0"
var_status = as.numeric(var_status)

# Also create continuous day-of-death variable to be used as phenotype
# Read in interp data
interp.df <- read.table(file = interp.FN, header = T, sep = ",")
interp.df <- as.data.frame(interp.df)
head(interp.df)
colnames(interp.df) <- c("dna.id", "sample", "family", "day.sampled", "state")
head(interp.df)

geno[1:5,1:5] # genotype matrix
indiv.df <- rownames(geno) # obtain the sample IDs from the genotype matrix
indiv.df <- as.data.frame(indiv.df)
head(indiv.df) # order of the samples in the genotype matrix
indiv.df <- separate(data = indiv.df, col = "indiv.df", into = c("family.state", "dna.id"), sep = "_", remove = F)
head(indiv.df)

# Combine interp with the ordered individual names
length(intersect(x = indiv.df$dna.id, y = interp.df$dna.id))
nrow(indiv.df)
all.df <- merge(x = indiv.df, y = interp.df, by = "dna.id", all.x = T, sort = F)
dim(all.df)
head(all.df)
head(indiv.df)
tail(all.df)
tail(indiv.df)

# Correct day.sampled for survivors
all.df$day.sampled[all.df$state=="Alive"] # which 'day sampled' are found in the Alive records
all.df$day.sampled[all.df$state=="Dead"]  # which 'day sampled' are found in the Dead records
# set alive state as day 7
all.df$day.sampled[all.df$state=="Alive"] <- "D7"
# Correct naming issue
all.df$day.sampled <- gsub(pattern = "Day 4", replacement = "D4", x = all.df$day.sampled)
table(all.df$day.sampled)


#GWAS with all families
gwaspheno = var_status
gwascovar = model.matrix(~as.factor(var_family))
gwasgeno = t(geno)
gwasgeno = cbind(rownames(gwasgeno),"X","Y",gwasgeno)
gwasanno = cbind(rownames(gwasgeno),
                 sapply(strsplit(rownames(gwasgeno),"_"), `[`, 2),
                 sapply(strsplit(sapply(strsplit(rownames(gwasgeno),"_"), `[`, 1),"Chr"), `[`, 2),
                 0)

write.table(gwaspheno, "03_results/gwaspheno.txt",row.names = F,col.names = F)
write.table(gwascovar, "03_results/gwascovar.txt",row.names = F,col.names = F)
write.table(gwasanno, "03_results/gwasanno.txt",row.names = F,col.names = F,quote = F)
write.table(gwasgeno, "03_results/gwasgeno.txt",row.names = F,col.names = F,quote = F)

# GWAS additional phenotype
gwaspheno2 <- as.numeric(gsub(pattern = "D", replacement = "", x = all.df$day.sampled))
write.table(gwaspheno2, "03_results/gwaspheno2.txt",row.names = F, col.names = F)

#### 05. Instructions for GEMMA (run in terminal) ####
##Run in command-line
# cd 03_results
##gemma-0.98.5 available at https://github.com/genetics-statistics/GEMMA
# ./gemma-0.98.5 -g gwasgeno.txt -p gwaspheno.txt -gk -maf 0.05 -o gwas_allfam
# ./gemma-0.98.5 -g gwasgeno.txt -p gwaspheno.txt -k output/gwas_allfam.cXX.txt -n 1 -c gwascovar.txt -a gwasanno.txt -maf 0.05 -lmm 4 -o gwas_allfam_covar
# ./gemma-0.98.5 -g gwasgeno.txt -p gwaspheno.txt -k output/gwas_allfam.cXX.txt -n 1 -a gwasanno.txt -maf 0.05 -lmm 4 -o gwas_allfam_nocovar

## Jump to plotting section #


#### Unused variations (subsets) of GWAS ####
# #GWAS with only MBP families
# gwaspheno = var_status[-which(var_family=="OFR6.10")]
# gwascovar = model.matrix(~as.factor(var_family[-which(var_family=="OFR6.10")]))
# gwasgeno = t(geno[-which(var_family=="OFR6.10"),])
# gwasgeno = cbind(rownames(gwasgeno),"X","Y",gwasgeno)
# gwasanno = cbind(rownames(gwasgeno),
#                  sapply(strsplit(rownames(gwasgeno),"_"), `[`, 2),
#                  sapply(strsplit(sapply(strsplit(rownames(gwasgeno),"_"), `[`, 1),"Chr"), `[`, 2),
#                  0)
# 
# write.table(gwaspheno,"gwaspheno.txt",row.names = F,col.names = F)
# write.table(gwascovar,"gwascovar.txt",row.names = F,col.names = F)
# write.table(gwasanno,"gwasanno.txt",row.names = F,col.names = F,quote = F)
# write.table(gwasgeno,"gwasgeno.txt",row.names = F,col.names = F,quote = F)


#./gemma-0.98.5 -g gwasgeno.txt -p gwaspheno.txt -gk -maf 0.05 -o gwas_mbpfam
#./gemma-0.98.5 -g gwasgeno.txt -p gwaspheno.txt -k output/gwas_mbpfam.cXX.txt -n 1 -c gwascovar.txt -a gwasanno.txt -maf 0.05 -lmm 4 -o gwas_mbpfam_covar
#./gemma-0.98.5 -g gwasgeno.txt -p gwaspheno.txt -k output/gwas_mbpfam.cXX.txt -n 1 -a gwasanno.txt -maf 0.05 -lmm 4 -o gwas_mbpfam_nocovar


# #GWAS with only F114
# gwaspheno = var_status[which(var_family=="F114")]
# gwasgeno = t(geno[which(var_family=="F114"),])
# gwasgeno = cbind(rownames(gwasgeno),"X","Y",gwasgeno)
# gwasanno = cbind(rownames(gwasgeno),
#                  sapply(strsplit(rownames(gwasgeno),"_"), `[`, 2),
#                  sapply(strsplit(sapply(strsplit(rownames(gwasgeno),"_"), `[`, 1),"Chr"), `[`, 2),
#                  0)
# 
# write.table(gwaspheno,"gwaspheno.txt",row.names = F,col.names = F)
# write.table(gwasanno,"gwasanno.txt",row.names = F,col.names = F,quote = F)
# write.table(gwasgeno,"gwasgeno.txt",row.names = F,col.names = F,quote = F)
# 
# #./gemma-0.98.5 -g gwasgeno.txt -p gwaspheno.txt -n 1 -a gwasanno.txt -maf 0.05 -lm 2 -o gwas_F114
# 
# #GWAS with only F115
# gwaspheno = var_status[which(var_family=="F115")]
# gwasgeno = t(geno[which(var_family=="F115"),])
# gwasgeno = cbind(rownames(gwasgeno),"X","Y",gwasgeno)
# gwasanno = cbind(rownames(gwasgeno),
#                  sapply(strsplit(rownames(gwasgeno),"_"), `[`, 2),
#                  sapply(strsplit(sapply(strsplit(rownames(gwasgeno),"_"), `[`, 1),"Chr"), `[`, 2),
#                  0)
# 
# write.table(gwaspheno,"gwaspheno.txt",row.names = F,col.names = F)
# write.table(gwasanno,"gwasanno.txt",row.names = F,col.names = F,quote = F)
# write.table(gwasgeno,"gwasgeno.txt",row.names = F,col.names = F,quote = F)
# 
# #./gemma-0.98.5 -g gwasgeno.txt -p gwaspheno.txt -n 1 -a gwasanno.txt -maf 0.05 -lm 2 -o gwas_F115
# 
# #GWAS with only F116
# gwaspheno = var_status[which(var_family=="F116")]
# gwasgeno = t(geno[which(var_family=="F116"),])
# gwasgeno = cbind(rownames(gwasgeno),"X","Y",gwasgeno)
# gwasanno = cbind(rownames(gwasgeno),
#                  sapply(strsplit(rownames(gwasgeno),"_"), `[`, 2),
#                  sapply(strsplit(sapply(strsplit(rownames(gwasgeno),"_"), `[`, 1),"Chr"), `[`, 2),
#                  0)
# 
# write.table(gwaspheno,"gwaspheno.txt",row.names = F,col.names = F)
# write.table(gwasanno,"gwasanno.txt",row.names = F,col.names = F,quote = F)
# write.table(gwasgeno,"gwasgeno.txt",row.names = F,col.names = F,quote = F)
# 
# #./gemma-0.98.5 -g gwasgeno.txt -p gwaspheno.txt -n 1 -a gwasanno.txt -maf 0.05 -lm 2 -o gwas_F116
# 
# #GWAS with only F117
# gwaspheno = var_status[which(var_family=="F117")]
# gwasgeno = t(geno[which(var_family=="F117"),])
# gwasgeno = cbind(rownames(gwasgeno),"X","Y",gwasgeno)
# gwasanno = cbind(rownames(gwasgeno),
#                  sapply(strsplit(rownames(gwasgeno),"_"), `[`, 2),
#                  sapply(strsplit(sapply(strsplit(rownames(gwasgeno),"_"), `[`, 1),"Chr"), `[`, 2),
#                  0)
# 
# write.table(gwaspheno,"gwaspheno.txt",row.names = F,col.names = F)
# write.table(gwasanno,"gwasanno.txt",row.names = F,col.names = F,quote = F)
# write.table(gwasgeno,"gwasgeno.txt",row.names = F,col.names = F,quote = F)
# 
# #./gemma-0.98.5 -g gwasgeno.txt -p gwaspheno.txt -n 1 -a gwasanno.txt -maf 0.05 -lm 2 -o gwas_F117
# 
# #GWAS with only OFR6.10
# gwaspheno = var_status[which(var_family=="OFR6.10")]
# gwasgeno = t(geno[which(var_family=="OFR6.10"),])
# gwasgeno = cbind(rownames(gwasgeno),"X","Y",gwasgeno)
# gwasanno = cbind(rownames(gwasgeno),
#                  sapply(strsplit(rownames(gwasgeno),"_"), `[`, 2),
#                  sapply(strsplit(sapply(strsplit(rownames(gwasgeno),"_"), `[`, 1),"Chr"), `[`, 2),
#                  0)
# 
# write.table(gwaspheno,"gwaspheno.txt",row.names = F,col.names = F)
# write.table(gwasanno,"gwasanno.txt",row.names = F,col.names = F,quote = F)
# write.table(gwasgeno,"gwasgeno.txt",row.names = F,col.names = F,quote = F)
# 
# #./gemma-0.98.5 -g gwasgeno.txt -p gwaspheno.txt -n 1 -a gwasanno.txt -maf 0.05 -lm 2 -o gwas_OFR6.10


#### 06. Plot GEMMA output ####
# All family, with covariate
# gemma_gwas2 = read.table("03_results/output/gwas_allfam_covar.assoc.txt")
# gemma_gwas = read.table("03_results/output/gwas_allfam_covar.assoc.txt", skip = 1)
# names(gemma_gwas) = gemma_gwas2[1,]

# Binary phenotype
gemma_output <- read.table(file = "03_results/output/gwas_allfam_covar.assoc.txt", header = T)
gemma_gwas   <- gemma_output
head(gemma_gwas)

nind <- length(colnames(gwasgeno)[colnames(gwasgeno)!=""])

pdf(file = "03_results/Manhattan_plot_all_fam_w_covar_binary_pheno.pdf", width = 6.5, height = 4.5)
fastman(gemma_gwas,
        chr = "chr",
        bp = "ps",
        p="p_lrt",
        genomewideline = -log10(0.05/nrow(gemma_gwas)),
        suggestiveline = NULL,
        cex=1.5,cex.lab=1.5,cex.axis=1,
        ylim=c(0,6),
        #main= paste0("All families with fixed covariate, dead/alive pheno (n = ", nind, ")")
        )
dev.off()

# Days-to-death phenotype
gemma_output <- read.table(file = "03_results/output/gwas_allfam_covar_pheno_day_to_death.assoc.txt", header = T)
gemma_gwas   <- gemma_output
head(gemma_gwas)

nind <- length(colnames(gwasgeno)[colnames(gwasgeno)!=""])

pdf(file = "03_results/Manhattan_plot_all_fam_w_covar_day-to-death_pheno.pdf", width = 6.5, height = 4.5)
fastman(gemma_gwas,
        chr = "chr",
        bp = "ps",
        p="p_lrt",
        genomewideline = -log10(0.05/nrow(gemma_gwas)),
        suggestiveline = NULL,
        cex=1.5,cex.lab=1.5,cex.axis=1,
        ylim=c(0,6),
        #main= paste0("All families with fixed covariate, days-to-death pheno (n = ", nind, ")")
)
dev.off()






#### Unused plotting ####
# gemma_gwas2 = read.table("output/gwas_allfam_nocovar.assoc.txt")
# gemma_gwas = read.table("output/gwas_allfam_nocovar.assoc.txt",skip = 1)
# names(gemma_gwas) = gemma_gwas2[1,]
# 
# fastman(gemma_gwas,
#         chr = "chr",
#         bp = "ps",
#         p="p_lrt",
#         genomewideline = -log10(0.05/nrow(gemma_gwas)),
#         suggestiveline = NULL,
#         cex=1.5,cex.lab=1.5,cex.axis=1,
#         ylim=c(0,6),
#         main="All familes + no covariable (n=165)")
# 
# 
# gemma_gwas2 = read.table("output/gwas_mbpfam_covar.assoc.txt")
# gemma_gwas = read.table("output/gwas_mbpfam_covar.assoc.txt",skip = 1)
# names(gemma_gwas) = gemma_gwas2[1,]
# 
# fastman(gemma_gwas,
#         chr = "chr",
#         bp = "ps",
#         p="p_lrt",
#         genomewideline = -log10(0.05/nrow(gemma_gwas)),
#         suggestiveline = NULL,
#         cex=1.5,cex.lab=1.5,cex.axis=1,
#         ylim=c(0,6),
#         main="MBP familes + fixed covariable (n=139)")
# 
# 
# 
# gemma_gwas2 = read.table("output/gwas_mbpfam_nocovar.assoc.txt")
# gemma_gwas = read.table("output/gwas_mbpfam_nocovar.assoc.txt",skip = 1)
# names(gemma_gwas) = gemma_gwas2[1,]
# 
# fastman(gemma_gwas,
#         chr = "chr",
#         bp = "ps",
#         p="p_lrt",
#         genomewideline = -log10(0.05/nrow(gemma_gwas)),
#         suggestiveline = NULL,
#         cex=1.5,cex.lab=1.5,cex.axis=1,
#         ylim=c(0,6),
#         main="MBP familes + no covariable (n=139)")
# 
# 
# 
# gemma_gwas2 = read.table("output/gwas_F114.assoc.txt")
# gemma_gwas = read.table("output/gwas_F114.assoc.txt",skip = 1)
# names(gemma_gwas) = gemma_gwas2[1,]
# 
# fastman(gemma_gwas,
#         chr = "chr",
#         bp = "ps",
#         p="p_lrt",
#         genomewideline = -log10(0.05/nrow(gemma_gwas)),
#         suggestiveline = NULL,
#         cex=1.5,cex.lab=1.5,cex.axis=1,
#         ylim=c(0,6),
#         main="F114 (n=32)")
# 
# 
# 
# gemma_gwas2 = read.table("output/gwas_F115.assoc.txt")
# gemma_gwas = read.table("output/gwas_F115.assoc.txt",skip = 1)
# names(gemma_gwas) = gemma_gwas2[1,]
# 
# fastman(gemma_gwas,
#         chr = "chr",
#         bp = "ps",
#         p="p_lrt",
#         genomewideline = -log10(0.05/nrow(gemma_gwas)),
#         suggestiveline = NULL,
#         cex=1.5,cex.lab=1.5,cex.axis=1,
#         ylim=c(0,6),
#         main="F115 (n=37)")
# 
# 
# gemma_gwas2 = read.table("output/gwas_F116.assoc.txt")
# gemma_gwas = read.table("output/gwas_F116.assoc.txt",skip = 1)
# names(gemma_gwas) = gemma_gwas2[1,]
# 
# fastman(gemma_gwas,
#         chr = "chr",
#         bp = "ps",
#         p="p_lrt",
#         genomewideline = -log10(0.05/nrow(gemma_gwas)),
#         suggestiveline = NULL,
#         cex=1.5,cex.lab=1.5,cex.axis=1,
#         ylim=c(0,6),
#         main="F116 (n=36)")
# 
# 
# 
# gemma_gwas2 = read.table("output/gwas_F117.assoc.txt")
# gemma_gwas = read.table("output/gwas_F117.assoc.txt",skip = 1)
# names(gemma_gwas) = gemma_gwas2[1,]
# 
# fastman(gemma_gwas,
#         chr = "chr",
#         bp = "ps",
#         p="p_lrt",
#         genomewideline = -log10(0.05/nrow(gemma_gwas)),
#         suggestiveline = NULL,
#         cex=1.5,cex.lab=1.5,cex.axis=1,
#         ylim=c(0,6),
#         main="F117 (n=34)")
# 
# 
# 
# gemma_gwas2 = read.table("output/gwas_OFR6.10.assoc.txt")
# gemma_gwas = read.table("output/gwas_OFR6.10.assoc.txt",skip = 1)
# names(gemma_gwas) = gemma_gwas2[1,]
# 
# fastman(gemma_gwas,
#         chr = "chr",
#         bp = "ps",
#         p="p_lrt",
#         genomewideline = -log10(0.05/nrow(gemma_gwas)),
#         suggestiveline = NULL,
#         cex=1.5,cex.lab=1.5,cex.axis=1,
#         ylim=c(0,6),
#         main="OFR6.10 (n=26)")



