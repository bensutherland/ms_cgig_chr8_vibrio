## Analysis script for OCP23 
# uses amplitools
# Ben Sutherland and Liam Surry, VIU
# initialized 2023-08-24

# If you are running from raw data, do step 1. If you've already prepared a genepop, go to step 2. 

#### 01. Prepare genepop in amplitools ####
# Open and source amplitools/01_scripts/00_initiator.R

# Convert proton to genepop
proton_to_genepop(hotspot_only = TRUE, neg_control = "Blank")

# Prepare the genepop
#./01_scripts/format_genepop.sh 02_input_data/prepped_matrices/<filename>.txt

# Copy the results to simple_pop_stats
# cp ./02_input_data/prepped_genepops/R_2023_07_26_12_44_23_user_GSS5PR-0268-78-Ampseq_Oyster_20230725_gen_data.gen ../simple_pop_stats/02_input_data/


#### 02. Filter in simple_pop_stats ####
# Open and source simple_pop_stats/01_scripts/simple_pop_stats_start.R, then choose '9-Pacific oyster'
## note: change variable 'on_network' to FALSE

#### 01. Load data ####
load_genepop(datatype = "SNP")
## note: input file is "02_input_data/R_2023_07_26_12_44_23_user_GSS5PR-0268-78-Ampseq_Oyster_20230725_gen_data.gen"

#### 02. Prepare data ####
##### 02.1 Manually assign population names based on samples present #####
#### TODO: AMPLITOOLS PREP FUNCTION: CREATE INDLIST ####
## Create a list of individuals for manual addition of population
indiv.df <- as.data.frame(indNames(obj))
colnames(indiv.df) <- "indiv"
head(indiv.df)

# Separate components of indiv ID (i.e., run, barcode, sample)
indiv.df <- separate(data = indiv.df, col = "indiv", into = c("run", "barcode", "indiv"), sep = "__", remove = T)
head(indiv.df)

# Use reduced indiv name as indname in genind
indNames(obj) <- indiv.df$indiv

# How many samples from each run? 
table(indiv.df$run) # 190, only a single run here

# Prepare to write out a clean text file to add pop attribute per sample
indiv.df <- as.data.frame(indiv.df[, "indiv"])
colnames(indiv.df) <- "indiv"
indiv.df$pop <- NA # Add dummy column to fill manually
head(indiv.df)

# Write out empty file to provide pop names
write.table(x = indiv.df, file = "02_input_data/my_data_ind-to-pop.txt"
            , sep = "\t", col.names = T, row.names = F
            , quote = F
)

### /END/ AMPLITOOLS PREP FUNCTION: CREATE INDLIST ###


# In folder above, *manually annotate* the output file above
# , save with "_annot.txt" appended, populate with pop names (no spaces)

### TODO: AMPLITOOLS CONNECT POPS ####
# Load annotated df
indiv_annot.df <- read.table(file = "02_input_data/my_data_ind-to-pop_annot.txt"
                             , header = T, sep = "\t"
                             #, quote = F
)

## Update population names
# Merge with the population annotation, do not sort
indiv_annot_in_order.df <- merge(x = indiv.df, indiv_annot.df, by = "indiv"
                                 , all.x = T, sort = FALSE # very necessary line
)

head(indiv_annot_in_order.df)

# Observe order remained same as (x) above
head(cbind(indiv_annot_in_order.df, indiv.df), n = 10)
tail(cbind(indiv_annot_in_order.df, indiv.df), n = 10)

### TODO: add data-check in this step###
# # Write a little test to be sure
# test <- cbind(indNames(obj), indiv_annot_in_order.df$indiv)
# table(test[,1] == test[ ,2])
# 
# # test <- cbind(indNames(obj), sort(indiv_annot_in_order.df$indiv))
# # table(test[,1] == test[ ,2])
# ## /END/ ##

# Assign the pop IDs to the genind
pop(obj) <- indiv_annot_in_order.df$pop.y
table((pop(obj)))


##### 02.2 Add in population colours #####
## Population colours
colours <- matrix(c("F1", "F0", "OAR", "darkgreen", "purple", "black"), nrow = 3, ncol = 2)
colnames(colours) <- c("my.pops", "my.cols")
colours

# Save out colours to be used downstream
colnames(x = colours) <- c("collection", "colour")
write.csv(x = colours, file = "00_archive/formatted_cols.csv", quote = F, row.names = F)


#### 03. Characterize missing data (indiv and loci) and filter ####
### This should be a function (SPS)
# Set variables to use for both plots
plot_width <- 8
plot_height <- 5
plot_cex <- 0.85
plot_pch <- 16

##### 03.1 Individuals - missing data #####
percent_missing_by_ind(df = obj)
head(missing_data.df)

# Add pop IDs to the missing data df
missing_data.df <- merge(x = missing_data.df, y = indiv_annot.df, by.x = "ind", by.y = "indiv", all.x = T)
head(missing_data.df)

# Add colours to the missing data df
colours
plot_cols.df <- merge(x = missing_data.df, y = colours, by.x = "pop", by.y = "collection", all.x = T
                      , sort = F
)
head(plot_cols.df)

# Plot missing data by individual
pdf(file = "03_results/geno_rate_by_ind.pdf", width = plot_width, height = plot_height)
plot(100 * (1 - plot_cols.df$ind.per.missing), ylab = "Genotyping rate (%)"
     , col = plot_cols.df$my.cols
     , las = 1
     , xlab = "Individual"
     , ylim = c(0,100)
     , pch=plot_pch
     , cex = plot_cex
)

abline(h = 50, lty = 3)

legend("bottomleft", legend = unique(plot_cols.df$pop)
       , fill = unique(plot_cols.df$my.cols)
       , cex = 0.7
       , bg = "white"
)
dev.off()

# Filter individuals by missing data
obj.df <- missing_data.df
head(obj.df)

keep <- obj.df[obj.df$ind.per.missing < 0.5, "ind"]

length(keep)
nInd(obj)

obj.filt <- obj[(keep)]
obj.filt

# Samples remaining after filters
table(pop(obj.filt))


##### 03.2 Loci - missing data #####
### TODO: SPS function
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
pdf(file = "03_results/geno_rate_by_marker.pdf", width = plot_width, height = plot_height)
plot(100 * (1- obj.df$marker.per.missing), xlab = "Marker", ylab = "Genotyping rate (%)", las = 1
     , ylim = c(0,100)
     , pch = plot_pch
     , cex = plot_cex
)
abline(h = 50
       #, col = "grey60"
       , lty = 3)
dev.off()

# Filter markers by genotyping rate
keep <- rownames(obj.df[obj.df$marker.per.missing < 0.5, ])

# How many loci will be removed? 
nLoc(obj.filt)
nLoc(obj.filt) - length(keep)

# Drop loci from genind
obj.all.filt <- obj.filt[, loc=keep]

# Rename back to obj
obj <- obj.all.filt


##### 03.3 Drop monomorphic loci #####
drop_loci(df = obj, drop_monomorphic = TRUE) # drops monomorphic markers

obj <- obj_filt


##### 03.4 Post-QC info collection #####
obj

## View the ind or loc names
inds <- indNames(obj)
loci <- locNames(obj)

# Save out which individuals have passed the filters
write.table(x = inds, file = "03_results/retained_individuals.txt", sep = "\t", quote = F
            , row.names = F, col.names = F
)

write.table(x = loci, file = "03_results/retained_loci.txt", sep = "\t", quote = F
            , row.names = F, col.names = F
)


##### 03.5 per marker stats and filters #####
# # MAF information
# maf_filt(data = obj, maf = 0.01)
# head(myFreq)
# 
# pdf(file = "03_results/maf_freq.pdf", width = 7, height = 5)
# hist(myFreq, breaks = 20, main = "", xlab = "MAF", las = 1)
# dev.off()
# 
# obj <- obj_maf_filt

## Per locus statistics
per_locus_stats(data = obj)
head(per_loc_stats.df)


# Plot Fst by Hobs
pdf(file = "03_results/per_locus_Fst_v_Hobs.pdf", width = 8, height = 5) 
plot(per_loc_stats.df$Fst, per_loc_stats.df$Hobs
     , las = 1
     , xlab = "Per locus FST"
     , ylab = "Per locus HOBS"
     , pch = 16
     , cex = 0.85
)
dev.off()

table(per_loc_stats.df$Hobs > 0.5) 

# Not dropping any loci based on HWP or HOBS at this time, given the non-random selection of individuals

# Remove specific individuals
keep <- indNames(obj)[grep(pattern = "_2$", x = indNames(obj), invert = T)]
obj <- obj[(keep)]
table(pop(obj))


#### 04. Analysis ####
####### Convert genepop to Rubias format #####
# Need to create a tab-delim stock code file in format of e.g., 
## row 1: collection	repunit
## row 2: boundary_bay	lower_mainland

# Here we will just create a df based on existing populations where collection = repunit
stock_code.df <- as.data.frame(unique(pop(obj)))
colnames(stock_code.df) <- "collection"
stock_code.df$repunit <- stock_code.df$collection
stock_code.df

# Write it out
write_delim(x = stock_code.df, file = "00_archive/stock_code.txt", delim = "\t", col_names = T)
micro_stock_code.FN <- "00_archive/stock_code.txt"
# this is for annotate_rubias(), for an unknown reason it requires the name micro_stock_code.FN

## Convert genepop to rubias
datatype <- "SNP" # required for genepop_to_rubias_SNP
as.data.frame(cbind(indNames(obj), as.character(pop(obj)))) # Note: BR27 should be VIU_F0 [confirmed]
obj # the current analysis object

# All filters applied
genepop_to_rubias_SNP(data = obj, sample_type = "reference", custom_format = TRUE, micro_stock_code.FN = micro_stock_code.FN)
print("Your output is available as '03_results/rubias_output_SNP.txt")
file.copy(from = "03_results/rubias_output_SNP.txt", to = "../amplitools/03_results/cgig_all_rubias.txt", overwrite = T)


#### Parentage Analysis ####
# Clear space, and launch amplitools initiator (i.e., 01_scripts/00_initiator.R)

# Using this output, move to "amplitools/01_scripts/ckmr_from_rubias.R"
ckmr_from_rubias(input.FN = "03_results/cgig_all_rubias.txt", parent_pop = "F0"
                 , offspring_pop = "F1", cutoff = 5
)

# Generate report from the output of CKMR-sim
prep_report(relationship = "PO")



# Plot the output results
graph_relatives(input.FN = "03_results/po_F0_vs_F1_pw_logl_5.txt", logl_cutoff = 5
                , drop_string = "G00", directed = F, plot_width = 5, plot_height = 5
)

graph_relatives(input.FN = "03_results/fs_offsp_F1_pw_logl_5.txt", logl_cutoff = 5
                , drop_string = "G00", directed = F, plot_width = 5, plot_height = 5
)

graph_relatives(input.FN = "03_results/fs_parent_F0_pw_logl_5.txt", logl_cutoff = 5
                , drop_string = "G00", directed = F, plot_width = 5, plot_height = 5
)


