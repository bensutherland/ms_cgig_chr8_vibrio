# simple_pop_stats characterizing and filtering for the OCV23 RADseq analysis
# B. Sutherland
# Initialized 2023-10-23
# Requires: 
#  - VCF file from Stacks in the input folder (set in user-set variables)
#  - simple_pop_stats repository

# Prior to running the following, source simple_pop_stats and choose Pacific oyster

#### 01. Load Data ####
# User-set variables
VCF.FN <- "02_input_data/populations.snps.vcf" # standard analysis
datatype <- "SNP" # simple_pop_stats requirement, normally assigned by load_genepop() when input is a genepop

# Read in genotype data
vcf <- read.vcfR(file = VCF.FN)
vcf

# Convert to genind file
obj <- vcfR2genind(x = vcf)

# Add population attribute based on sample names
assign_pops.df <- indNames(obj)
assign_pops.df <- as.data.frame(assign_pops.df)
colnames(assign_pops.df) <- "indiv"
assign_pops.df <- separate(data = assign_pops.df, col = "indiv", into = c("pop", "ind_identifier"), sep = "-", remove = F)
head(assign_pops.df)
unique(assign_pops.df$pop) # OK

pop(obj) <- assign_pops.df$pop

table(pop(obj))


#### 02. Prepare colour file ####
# Define population colours
pops_in_genepop.df <- unique(pop(obj))
pops_in_genepop.df <- as.data.frame(pops_in_genepop.df)
colnames(pops_in_genepop.df) <- "collection"
pops_in_genepop.df

# Designate colours
pops_in_genepop.df$colour <- c("black", "darkgreen", "orange", "blue", "purple", "red")

# Save out
write.csv(x = pops_in_genepop.df, file = "00_archive/pop_cols.csv", quote = F, row.names = F)


#### 03. Characterize data and filter ####
characterize_genepop(obj)

##### 03.1. Individuals - missing data #####
percent_missing_by_ind(df = obj)
head(missing_data.df)

## Plot
# Add colours
missing_data.df <- merge(x = missing_data.df, y = pops_in_genepop.df, by.x = "pop", by.y = "collection", all.x = T
                      , sort = F
)

# Plot missing data by individual
pdf(file = "03_results/geno_rate_by_ind.pdf", width = 7, height = 4)
plot(1 - missing_data.df$ind.per.missing, ylab = "Genotyping rate (%)"
     , col = missing_data.df$colour
     , las = 1
     , pch = 16
     , cex = 0.8
     , xlab = "Individual"
     , ylim = c(0,1)
)

abline(h = 0.7, lty = 3)

legend("bottomleft"
       , legend = pops_in_genepop.df$collection
       , fill = pops_in_genepop.df$colour
       , cex = 0.85
       , bg = "white"
)
dev.off()


## Filter individuals based on GR
#   keep indiv with >=70% genotyping rate (% missing < 0.3)
keep <- missing_data.df[missing_data.df$ind.per.missing < 0.3, "ind"]
print(paste0("Retaining ", length(keep), " of the total ", nInd(obj), " individuals"))

# Retain only the keep indiv
obj.filt <- obj[(keep)]
obj.filt
table(pop(obj.filt))

# Rename
obj <- obj.filt

## Retain names of retained indiv and loci
inds <- indNames(obj)
loci <- locNames(obj)

write.table(x = inds, file = "03_results/retained_individuals.txt", sep = "\t", quote = F
            , row.names = F, col.names = F
)

write.table(x = loci, file = "03_results/retained_loci.txt", sep = "\t", quote = F
            , row.names = F, col.names = F
)

# Percent missing post filters? 
mean(missing_data.df[missing_data.df$ind.per.missing < 0.3, "ind.per.missing"]) # 0.042
sd(missing_data.df[missing_data.df$ind.per.missing < 0.3, "ind.per.missing"])   # 0.038

# Summarize number survived and number dead post filters
table(gsub(pattern = "_.*", replacement = "", x = indNames(obj)))



##### 03.2 Loci by HWE and excess Hobs #####
# Not running HOBS or HWE filter, due to family effects expected


##### 03.3 Post-indiv missing data filter allele freq calculations #####
maf_filt(data = obj, maf = 0.01)
#myFreq
obj <- obj_maf_filt # Rename

# Plot
pdf(file = paste0("03_results/maf_hist_post_filter.pdf"), width = 6, height = 4)
hist(myFreq
     #, proba=T # note: does not sum to 1, not worth using
     , col="grey", xlab = "Minor allele frequency (MAF)"
     , main = ""
     #, ylim = c(0, 2500)
     , ylab = "Number of loci"
     , las = 1
     , breaks = 20
)
text(x = 0.4, y = 1000, labels = paste("n = ", length(myFreq), " loci", sep = "" ))
dev.off()

# Save out the MAF calculation as a table
myFreq <- round(myFreq, digits = 3)
write.table(x = myFreq, file = "03_results/allele_freq_retained_loci.txt"
            , sep = "\t", quote = F
            , row.names = T, col.names = F
)

# Exploration
table(myFreq < 0.01)
table(myFreq < 0.1) 


##### Clean up space ####
rm(obj.filt)
rm(obj_maf_filt)

gc()

#### 0.4 Export ####
# Write out object
save.image(file = "03_results/post_all_filters.RData")

# Go to "01_scripts/03_sps_analysis.R
