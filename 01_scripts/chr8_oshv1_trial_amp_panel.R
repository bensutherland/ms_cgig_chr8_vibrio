# Analyze amplicon panel data for the CHR8 OsHV-1 trial at CATC
# B. Sutherland, initialized 2024-04-30

# Source amplitools initiator
# Copy the two input files into amplitools/02_input_data

# Prepare genotype block
proton_to_genepop(neg_control = "BLANK")

# Prepare the genepop
#./01_scripts/format_genepop.sh 02_input_data/prepped_matrices/<filename>.txt

# Copy the results to simple_pop_stats
# cp ./02_input_data/prepped_genepops/*.gen ../simple_pop_stats/02_input_data/


### B. Filter in simple_pop_stats ####
# Open and source simple_pop_stats/01_scripts/simple_pop_stats_start.R, then choose '9-Pacific oyster'
## note: change variable 'on_network' to FALSE

# NOTE: may need a multimapper file

# load_genepop(datatype = "SNP") # simple_pop_stats/02_input_data/G0923-21-VIUN_File1_gen_data.gen
# obj.file1 <- obj

load_genepop(datatype = "SNP") # simple_pop_stats/02_input_data/G0923-21-VIUN_File2_gen_data.gen
obj.file2 <- obj

rm(obj)

# Rename
# simplify_names(df = obj.file1, format = "amplitools")
# obj.file1 <- obj_simplified

simplify_names(df = obj.file2, format = "amplitools")
obj.file2 <- obj_simplified

# indNames(obj.file1)
indNames(obj.file2) <- gsub(pattern = "_ReAMP", replacement = "", x = indNames(obj.file2))

# length(intersect(x = indNames(obj.file1), y = indNames(obj.file2))) # 218
# setdiff(x = indNames(obj.file1), y = indNames(obj.file2)) # no extra in file 1
# setdiff(x = indNames(obj.file2), y = indNames(obj.file1)) # 22 extra in file 2

# Note: it appears that File 2 is the superior file, so we should just work with that file for now. 
obj <- obj.file2

# Aside, calculate MAF, without any filters applied
maf_filt(data = obj, maf = 0.0001)

write.table(x = myFreq, file = paste0("03_results/myFreq_", nLoc(obj_filt), "_loci.txt")
            , quote = F, sep = "\t"
            )

generate_popmap(df = obj)
# manually enter, then save out as *_annot.txt

# Annotate samples
annotate_from_popmap(df = obj, popmap.FN = "00_archive/my_data_ind-to-pop_annot.txt")
obj <- obj_annot

# Characterize missing data
percent_missing_by_ind(df = obj)
head(missing_data.df)
missing_data.df$pop <- NA
missing_data.df[grep(pattern = "_114_", x = missing_data.df$ind), "pop"] <- "F114"
missing_data.df[grep(pattern = "_115_", x = missing_data.df$ind), "pop"] <- "F115"
missing_data.df[grep(pattern = "_116_", x = missing_data.df$ind), "pop"] <- "F116"
missing_data.df[grep(pattern = "_117_", x = missing_data.df$ind), "pop"] <- "F117"
head(missing_data.df)
tail(missing_data.df)

# Plot missing data
# Plot missing data by individual
pdf(file = "03_results/geno_rate_by_ind.pdf", width = 9, height = 6)
plot(100 * (1 - missing_data.df$ind.per.missing), ylab = "Genotyping rate (%)"
     , col = as.factor(missing_data.df$pop)
     , las = 1
     , xlab = "Individual"
     , ylim = c(0,100)
     #, pch=plot_pch
     , cex = 0.8
)

abline(h = 50, lty = 3)
abline(h = 70, lty = 2)

legend("bottomleft", legend = unique(missing_data.df$pop)
       , fill = as.factor(unique(missing_data.df$pop))
       , cex = 0.8
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

## TODO: a per locus filter, but write a function for it ##

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



#obj.bck <- obj
drop_loci(drop_monomorphic = T) ### TODO: need to require an input object identifier
# Drops 174 monomorphic loci
obj <- obj_filt
obj

##### Analysis ####
pca_from_genind(data = obj, PCs_ret = 4, plot_eigen = T, plot_allele_loadings = F, retain_pca_obj = T)

# Write out loci to keep
keep_loci.vec <- locNames(obj)
write.table(x = keep_loci.vec, file = "03_results/retained_loci.txt", sep = "\t", quote = F, row.names = F, col.names = F)

# Write out indiv to keep
keep_inds.vec <- indNames(obj)
write.table(x = keep_inds.vec, file = "03_results/retained_inds.txt", sep = "\t", quote = F, row.names = F, col.names = F)
