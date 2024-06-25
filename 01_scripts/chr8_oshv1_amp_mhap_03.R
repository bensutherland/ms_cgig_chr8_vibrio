# Analyze mhap data as genepop, viewing general trends
# Requires that a genepop has been generated using amplitools
# B. Sutherland
# 2024-06-25

# Source simple_pop_stats initiator

# Read in data
obj <- read.genepop(file = "02_input_data/genepop.gen", ncode = 2)
obj

# See indiv names
head(indNames(x = obj))

# Generate popmap, then manually edit
generate_popmap(df = obj)

# Annotate samples
annotate_from_popmap(df = obj, popmap.FN = "00_archive/my_data_ind-to-pop_annot.txt", convert_to_alt_ID = TRUE)
obj <- obj_annot
head(indNames(obj))

# Characterize missing data
percent_missing_by_ind(df = obj)
head(missing_data.df, n = 10)
missing_data.df$pop <- NA
missing_data.df[grep(pattern = "_114_", x = missing_data.df$ind), "pop"] <- "F114"
missing_data.df[grep(pattern = "_115_", x = missing_data.df$ind), "pop"] <- "F115"
missing_data.df[grep(pattern = "_116_", x = missing_data.df$ind), "pop"] <- "F116"
missing_data.df[grep(pattern = "_117_", x = missing_data.df$ind), "pop"] <- "F117"
missing_data.df[grep(pattern = "_wgrs", x = missing_data.df$ind), "pop"] <- "F0_wgrs"
missing_data.df[grep(pattern = "_amp", x = missing_data.df$ind), "pop"] <- "F0_amp"

head(missing_data.df)
tail(missing_data.df)
unique(missing_data.df$pop)
table(missing_data.df$pop)

# Plot missing data
# Plot missing data by individual
pdf(file = "03_results/geno_rate_by_ind.pdf", width = 9, height = 6)
plot(100 * (1 - missing_data.df$ind.per.missing), ylab = "Genotyping rate (%)"
     , col = as.factor(missing_data.df$pop)
     , las = 1
     , xlab = "Individual"
     , ylim = c(-50,100)
     #, pch=plot_pch
     , cex = 0.8
)

abline(h = 50, lty = 3)
abline(h = 70, lty = 2)

legend("bottomleft", legend = unique(missing_data.df$pop)
       , fill = as.factor(unique(missing_data.df$pop))
       , cex = 0.6
       , bg = "white"
)
dev.off()

# Filter based on missing data
max_missing <- 0.3
head(missing_data.df)

keep <- missing_data.df[missing_data.df$ind.per.missing < max_missing, "ind"]

length(keep)
nInd(obj)

obj.filt <- obj[(keep)]
obj.filt

# Samples remaining after filters
table(pop(obj.filt))


##### Loci - missing data #####
# Filter loci based on missing data
obj.df <- genind2df(obj.filt)
obj.df[1:5,1:5]
obj.df <- t(obj.df)
obj.df[1:5,1:5]
obj.df <- obj.df[2:nrow(obj.df),] # remove pop row
obj.df[1:5,1:5]
dim(obj.df)
str(obj.df)

obj.df <- as.data.frame(obj.df)
dim(obj.df)
str(obj.df)
obj.df[1:5,1:5] # See top left of file
obj.df[(dim(obj.df)[1]-5):dim(obj.df)[1], (dim(obj.df)[2]-5):dim(obj.df)[2]] # See bottom right of file

# Add collector col
obj.df$marker.per.missing <- NA

for(i in 1:(nrow(obj.df))){
  
  # Per marker                      sum all NAs for the marker, divide by total number markers
  obj.df$marker.per.missing[i] <-  (sum(is.na(obj.df[i,]))-1) / (ncol(obj.df)-1) 
  
}


# Plot missing data by marker
pdf(file = "03_results/geno_rate_by_marker.pdf", width = 9, height = 6)
plot(100 * (1- obj.df$marker.per.missing), xlab = "Marker", ylab = "Genotyping rate (%)", las = 1
     , ylim = c(0,100)
     #, pch = 1
     #, cex = plot_cex
)
abline(h = 50
       #, col = "grey60"
       , lty = 3)
dev.off()

# Filter markers by genotyping rate
keep <- rownames(obj.df[obj.df$marker.per.missing < max_missing, ])

# How many loci will be removed? 
nLoc(obj.filt)
nLoc(obj.filt) - length(keep)

# Drop loci from genind
obj.all.filt <- obj.filt[, loc=keep]

# Rename back to obj
obj <- obj.all.filt
obj

# Drop monomorphs
drop_loci(drop_monomorphic = T) ### TODO: need to require an input object identifier
obj <- obj_filt
obj

##### Analysis ####
pca_from_genind(data = obj, PCs_ret = 4, plot_eigen = F, plot_allele_loadings = F, retain_pca_obj = T)

#### This ends the preliminary analysis ####
