# Import VCF from wgrs_workflow and perform general population statistics
# B. Sutherland (initialized 2024-02-09)

### Front Matter ####
# prior to running the following, clear the workspace, source simple_pop_stats and choose Pacific oyster

## Install and load packages
#install.packages("vcfR")
library("vcfR")

# Set variables
genos.FN <- "../wgrs_workflow_COARL_v.0.2/05_genotyping/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.001_AF_0.05_LD0.5w50kb_subset0.05.vcf"
indiv.FN <- "../accessory_files/COARL2_parental_genotyping_label_map_2024-02-12.txt"
cross.FN <- "../accessory_files/COARL2_60Fert_AvgPhenotypeData_2024_04_03.txt"

# User set variables
#plot_by <- "sex"
plot_by <- "family"


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


#### 02. Load indiv label and sample name correspondence and rename ####
# read in correspondence file
conversion.df <- read.table(file = indiv.FN, header = T, sep = "\t")
head(conversion.df)
unique(conversion.df$Oyster.ID)

# Remove samples that are not part of the project from the conversion table
conversion.df <- conversion.df[grep(pattern = "COARL2", x =  conversion.df$Tube.label), ]

# Remove the parenthetical information from Oyster.ID, as we will use this for population info
conversion.df$Oyster.ID <- gsub(pattern = "\\(.*", replacement = "", x = conversion.df$Oyster.ID)
table(conversion.df$Oyster.ID)

# Convert from tube_label to sample_ID
head(inds.df)       # this is the ordered identifiers from the data
dim(inds.df)
head(conversion.df) # this is the correspondence to the true sample IDs
dim(conversion.df)

# Combine the VCF samples with the sample information
ordered_conversion.df <- merge(x = inds.df, y = conversion.df
                          , by.x = "inds.df", by.y = "Tube.label"
                          , sort = F   # important to not sort
                          ) 

head(cbind(inds.df, ordered_conversion.df)) # order of first two cols should be retained
tail(cbind(inds.df, ordered_conversion.df)) # order of first two cols should be retained
dim(ordered_conversion.df)

# Use the ordered identifier to rename samples
indNames(my_data.gid) <- ordered_conversion.df$Sample.ID

# Add individual sex info as column
sex <- gsub(pattern = "_.*", replacement = "", x = indNames(my_data.gid))
sex <- as.character(sex)
ordered_conversion.df$sex <- sex
rm(sex)
head(ordered_conversion.df)

# Set population with sex or family
if(plot_by=="sex"){
 
  pop(my_data.gid) <- ordered_conversion.df$sex
  
}else if(plot_by=="family"){
  
  pop(my_data.gid) <- ordered_conversion.df$Oyster.ID
  
}

unique(pop(my_data.gid))


#### 03. Basic exploratory analyses ####
pca_from_genind(data = my_data.gid
                , PCs_ret = 4, plot_eigen = T, retain_pca_obj = T
                , plot_allele_loadings = F
                , width = 15, height = 15, plot_label = T
                , colour_file = NULL
                , plot_ellipse = F
                )


# Run relatedness analysis
relatedness_calc(data = my_data.gid, datatype = "SNP")
datatype <- "SNP"

# Get filename for the latest kinship analysis
kinship.FN <- head(sort(list.files(path = "03_results/", pattern = "kinship_analysis", full.names = T)), n = 1)

relatedness_plot(file = kinship.FN, same_pops = FALSE, plot_by = "names"
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

# Move to 01_scripts/COARL_02.R
