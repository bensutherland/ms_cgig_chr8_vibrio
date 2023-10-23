# simple_pop_stats characterizing and filtering for the OCV23 RADseq analysis
# B. Sutherland
# Initialized 2023-10-23
# Requires running "ms_cgig_chr8/01_scripts/01_import_plink_to_genind.R" first

# Prior to running the following, source simple_pop_stats and choose Pacific oyster

#### 01. Load Data ####
load(file = "../ms_cgig_chr8/03_results/prepared_genind.RData") # loaded from prerequisite script above
datatype <- "SNP" # normally assigned by load_genepop() when input is a genepop

obj <- my.data.gid
obj

#### 02. Prepare Data ####
characterize_genepop(obj)

### Prepare Colours ###
## Define population colours
pops_in_genepop.df <- unique(pop(obj))
pops_in_genepop.df <- as.data.frame(pops_in_genepop.df)
colnames(pops_in_genepop.df) <- "collection"
pops_in_genepop.df

# Manually designate colours
pop_colours <- matrix(c("OSU.FO", "F114", "F115", "F116", "F117", "OFR6.10", "black", "darkgreen", "yellow", "blue", "purple", "magenta"), nrow = 6, ncol = 2)
colnames(pop_colours) <- c("collection", "colour")
pop_colours

# Connect colours to present populations
colours <- merge(x = pops_in_genepop.df, y =  pop_colours, by = "collection"
                 #, sort = F
                 , all.x = T
)
colours

# Clean space
rm(pop_colours)
rm(pops_in_genepop.df)

# Save out colours to use later
colours
write.csv(x = colours, file = "00_archive/formatted_cols.csv", quote = F, row.names = F)


#### 03. Characterize missing data (indiv) and filter ####
##### 03.1 Individuals - missing data #####
percent_missing_by_ind(df = obj)
head(missing_data.df)

# What is the average missing data, prior to removals
mean(missing_data.df$ind.per.missing) # 0.052
sd(missing_data.df$ind.per.missing)   # 0.072

### Plot per-individual missing data ###
# Provide population IDs to missing data, based on names
missing_data.df$pop <- rep(x = NA, times = nrow(missing_data.df))

# Provide population based on the individual name
missing_data.df$pop[grep(pattern = "F0-F0", x = missing_data.df$ind)] <- "OSU.FO"
missing_data.df$pop[grep(pattern = "F114", x = missing_data.df$ind)] <- "F114"
missing_data.df$pop[grep(pattern = "F115", x = missing_data.df$ind)] <- "F115"
missing_data.df$pop[grep(pattern = "F116", x = missing_data.df$ind)] <- "F116"
missing_data.df$pop[grep(pattern = "F117", x = missing_data.df$ind)] <- "F117"
missing_data.df$pop[grep(pattern = "OFR6.10", x = missing_data.df$ind)] <- "OFR6.10"

table(missing_data.df$pop)

# Integrate colours into dataframe for plotting. Note: don't sort, as it is still in the same order as the obj
colours
plot_cols.df <- merge(x = missing_data.df, y = colours, by.x = "pop", by.y = "collection", all.x = T
                      , sort = F
)

# Plot missing data by individual
pdf(file = "03_results/geno_rate_by_ind.pdf", width = 8, height = 5)
plot(1 - plot_cols.df$ind.per.missing, ylab = "Genotyping rate (%)"
     , col = plot_cols.df$colour
     , las = 1
     , xlab = "Individual"
     , ylim = c(0,1)
)

abline(h = 0.7, lty = 3)

legend("bottomleft", legend = unique(plot_cols.df$pop)
       , fill = unique(plot_cols.df$colour)
       , cex = 1.0
       , bg = "white"
)
dev.off()


## Filter individuals
# Identify which samples to retain based on genotyping rate
#  to keep inds with >=70% genotyping rate (% missing < 0.3)
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


##### 03.2 Loci by HWE and excess Hobs #####
# Not running HOBS or HWE filter, due to family effects expected


##### 03.3 Post-indiv missing data filter allele freq calculations #####
maf_filt(data = obj, maf = 0.01)
myFreq
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


# ##### 03.4 Variants per pop #####
# obj.sep <- seppop(x = obj, drop = TRUE) # Note: here "drop" is necessary to discard alleles that are no longer present in the subset of the data
# drop_loci(df = obj.sep$BC, drop_monomorphic = T)
# drop_loci(df = obj.sep$JPN, drop_monomorphic = T)
# drop_loci(df = obj.sep$VIU, drop_monomorphic = T)

##### Clean up space ####
rm(my.data.gid)
rm(obj_filt)
rm(obj_maf_filt)
rm(obj.filt)
gc()

#### 0.4 Export ####
# Write out object
save.image(file = "03_results/post_all_filters.RData")

# Go to "01_scripts/03_sps_analysis.R
