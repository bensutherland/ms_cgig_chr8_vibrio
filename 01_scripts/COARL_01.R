# Import VCF from wgrs_workflow and perform general population statistics
# B. Sutherland (2024-02-09)

### Front Matter ####
# Clean space
# rm(list=ls())

# Prior to running the following, source simple_pop_stats and choose Pacific oyster

## Install and load packages (not included in sps)
#install.packages("vcfR")

library("vcfR")

## Info
# sessionInfo()

# Set variables
genos.FN <- "../ms_cgig_chr8/02_input_data/mpileup_calls_filt_AF_0.05_LD.0.5.50kb_subset_0.01.vcf"
indiv.FN <- "../ms_cgig_chr8/00_archive/COARL2_parental_genotyping_label_map.txt"
cross.FN <- "../ms_cgig_chr8/00_archive/COARL2_crosses_and_phenos.txt"


#### 01. Load genotype data and format ####
# read in vcf
my_data.vcf <- read.vcfR(file = genos.FN)
my_data.vcf

# convert vcf to genind
my_data.gid <- vcfR2genind(x = my_data.vcf)
my_data.gid

## Simplify sample names and rename from tube label 
inds.df <- indNames(my_data.gid) # what are the indiv names and order in the data? 
inds.df <- as.data.frame(inds.df)
colnames(inds.df) <- "fastq.name"
head(inds.df)

# separate into component parts of fastq name
inds.df <- separate(data = inds.df, col = "fastq.name", into = c("project", "sample.num", "file.info")
                     , sep = "-", remove = T
                   )
head(inds.df)

# create matchable name (i.e., COARL2_01)
inds.df <- paste0(inds.df$project, "_", inds.df$sample.num)
inds.df <- as.data.frame(inds.df)
head(inds.df) 

# Data checking
#cbind(indNames(my_data.gid), inds.df) # this remains in the order of the data

#### 02. Load correspondence and rename ####
# read in correspondence file
conversion.df <- read.table(file = indiv.FN, header = T, sep = "\t")
head(conversion.df)

# Convert from tube_label to sample_ID
head(inds.df)       # this is the ordered identifiers from the data
dim(inds.df)
head(conversion.df) # this is the correspondence to the true sample IDs
dim(conversion.df)  # note there are extra lines in addition to COARL project

# Remove other project info from the conversion file
conversion.df <- conversion.df[grep(pattern = "COARL", x = conversion.df$tube_label), ]

ordered_conversion.df <- merge(x = inds.df, y = conversion.df
                          , by.x = "inds.df", by.y = "tube_label", sort = F
                          ) # important to not sort

head(conversion.df) # ordered identifiers 
tail(conversion.df)
dim(conversion.df)

# Use the ordered identifier to rename samples
indNames(my_data.gid) <- conversion.df$sample_ID

# Individual sex
sex <- gsub(pattern = "_.*", replacement = "", x = indNames(my_data.gid))
sex <- as.character(sex)

# use sex as population
pop(my_data.gid) <- sex


#### 03. Basic analyses ####
pca_from_genind(data = my_data.gid, PCs_ret = 4, plot_eigen = T, retain_pca_obj = T
                , plot_allele_loadings = F, width = 15, height = 15, plot_label = T
                )


# Run relatedness analysis
relatedness_calc(data = my_data.gid, datatype = "SNP")
datatype <- "SNP"
relatedness_plot(file = "03_results/kinship_analysis_2024-02-11.Rdata", same_pops = FALSE, plot_by = "names"
                 , plot_by_group = F)


# Plot MAF distribution (note: keep in mind that the data is LD-filtered)
my_data.gid
maf_filt(data = my_data.gid, maf = 0.01)
obj <- obj_maf_filt


# Plot
pdf(file = paste0("03_results/MAF_hist.pdf"), width = 7, height = 4)
hist(myFreq
     #, proba=T # note: does not sum to 1, not worth using
     , col="grey", xlab = "MAF"
     , main = ""
     #, ylim = c(0, 1250)
     , ylab = "Number of loci"
     , las = 1
     , breaks = 50
)
text(x = 0.4, y = 6000, labels = paste("n = ", length(myFreq), " loci", sep = "" ))
abline(v = 0.01, lty = 2)

dev.off()

# Save out the MAF calculation as a table
myFreq <- round(myFreq, digits = 3)
write.table(x = myFreq, file = "03_results/allele_freq.txt"
            , sep = "\t", quote = F
            , row.names = T, col.names = F
)


save.image(file = "03_results/prepared_data.RData")


#### 04. Bring in cross information, per cross, per locus, infer offspring AF ####
# 


