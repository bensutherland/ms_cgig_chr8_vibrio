# Compare genotypes between technologies
# B. Sutherland (VIU)
# 2024-07-08

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

# Load libraries
#install.packages("vcfR")
require(vcfR)
require(tidyr)

## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path) # take main directory
setwd(current.path)
rm(current.path)

# Set variable names
panel_records.FN <- "../compare_amp_panel_parents_and_wgrs_parents/0002.vcf" # records from panel shared in both
wgrs_records.FN  <- "../compare_amp_panel_parents_and_wgrs_parents/0003.vcf" # records from wgrs shared in both

# Read in data
panel.vcf <- read.vcfR(file = panel_records.FN)
wgrs.vcf  <- read.vcfR(file = wgrs_records.FN)

panel.vcf
wgrs.vcf

# Extract genotype portion
panel.df <- extract.gt(panel.vcf)
wgrs.df  <- extract.gt(wgrs.vcf)

panel.df[1:5,1:5]
wgrs.df[1:5,1:5]

## Update indiv names to match the two datatypes
# panel
panel_indiv_names.df <- colnames(panel.df)
panel_indiv_names.df <- as.data.frame(panel_indiv_names.df)
colnames(panel_indiv_names.df) <- "indiv"
panel_indiv_names.df$indiv <- gsub(pattern = "OCP_063_1", replacement = "55-41F_p", x = panel_indiv_names.df$indiv)
panel_indiv_names.df$indiv <- gsub(pattern = "OCP_068_2", replacement = "65-19M_p", x = panel_indiv_names.df$indiv)
panel_indiv_names.df$indiv <- gsub(pattern = "OCP_080_1", replacement = "65-8F_p", x = panel_indiv_names.df$indiv)
panel_indiv_names.df$indiv <- gsub(pattern = "OCP_135_1", replacement = "65-4F_p", x = panel_indiv_names.df$indiv)
panel_indiv_names.df$indiv <- gsub(pattern = "OCP_138_2", replacement = "79-1M_p", x = panel_indiv_names.df$indiv)
panel_indiv_names.df$indiv <- gsub(pattern = "OCP_148_2", replacement = "79-13M_p", x = panel_indiv_names.df$indiv)
panel_indiv_names.df$indiv <- gsub(pattern = "OCP_171_1", replacement = "58-33M_p", x = panel_indiv_names.df$indiv)
panel_indiv_names.df$indiv <- gsub(pattern = "OCP_177_1", replacement = "58-9F_p", x = panel_indiv_names.df$indiv)
colnames(panel.df) <- panel_indiv_names.df$indiv

# wgrs
wgrs_indiv_names.df <- colnames(wgrs.df)
wgrs_indiv_names.df <- as.data.frame(wgrs_indiv_names.df)
colnames(wgrs_indiv_names.df) <- "indiv"
wgrs_indiv_names.df <- separate(data = wgrs_indiv_names.df, col = "indiv", into = c("fam", "ind", "suffix"), sep = "-", remove = T)
wgrs_indiv_names.df$full_name <- paste0(wgrs_indiv_names.df$fam, "-", wgrs_indiv_names.df$ind, "_g")
colnames(wgrs.df) <- wgrs_indiv_names.df$full_name

# Inspect
panel.df[1:5,]
wgrs.df[1:5,1:5]

# Combine the dfs based on shared rownames (loci)
all.df <- cbind(panel.df, wgrs.df)
dim(all.df)
all.df[1:5,1:5]

# Order so that same indivs w diff tech are beside each other
all.df <- all.df[, order(colnames(all.df))]
all.df[1:5,]

# Remove non-hotspot loci (all hotspot loci have mnames that consist of a single numeric string)
all_hotspot.df <- all.df[grep(pattern = "_", x = rownames(all.df), invert = T), ]
all_hotspot.df[1:10,1:8]

# Rename df
dat.df <- all_hotspot.df

## Loop to tally and calculate proportions correct
# Create empty df
eval.df <- NULL # output df name
eval.df$mname <- rownames(dat.df)
eval.df$na.tally <- NA
eval.df$tally <- NA
eval.df$pair.tally <- NA
eval.df <- as.data.frame(eval.df)
head(eval.df)

# Define unique indivs
unique_inds.vec <- unique(gsub(pattern = "_g|_p", replacement = "", x = colnames(dat.df)))

# Set nulls
moi <- NULL; tally <- NULL; na.tally <- NULL; pair.tally <- NULL; 
panel.call <- NULL; genome.call <- NULL

# Loop
for(i in 1:length(eval.df$mname)){
  
  ## Debugging
  #print(i)
  
  # Select the mname for this iteration
  moi <- eval.df$mname[i]
  
  # Reset tallies to zero for this mname
  tally <- 0
  na.tally <- 0
  pair.tally <- 0
  
  # Evaluate across indivs
  for(j in 1:length(unique_inds.vec)){
    
    ## Debugging
    #print(j)
    
    # Determine call by panel and wgrs for this indiv and mname
    panel.call  <- dat.df[moi, paste0(unique_inds.vec[j], "_p")]
    genome.call <- dat.df[moi, paste0(unique_inds.vec[j], "_g")]
    
    # If neither call is NA
    if(!is.na(panel.call) & !is.na(genome.call)){
      
      # Record that a pair is evaluable
      pair.tally <- pair.tally + 1
      
      # If the two calls are concordant
      if(panel.call==genome.call){
        
        # Record that a concordant call was made
        tally <- tally + 1
        
      }else{
        
        # Keep concordant tally constant (no addition)
        tally <- tally
        
      }
    
    # If an NA is observed  
    }else{
      
      # Record the number of NAs observed
      na.tally <- na.tally + sum(is.na(c(panel.call, genome.call)))
      
    }
     
  }
  
  # Record the tallies for the mname
  eval.df[i,"tally"]         <- tally
  eval.df[i,"na.tally"]      <- na.tally
  eval.df[i,"pair.tally"]    <- pair.tally
    
}
  

head(eval.df)

# Add a percentage correct (i.e., how many of the present pairs were found to be concordant?)
eval.df$percent.corr <- eval.df$tally / eval.df$pair.tally
head(eval.df)

pdf(file = "03_results/histogram_concordant_calls.pdf", width = 7.5, height = 4.75)
par(mfrow=c(1,2))
hist(eval.df$percent.corr
     , xlab = "Proportion concordant (no NAs) per marker", las = 1
     , main = ""
)
hist(eval.df$tally
     , xlab = "Number concordant calls per marker", las = 1
     , main = ""
     )
dev.off()


write.table(x = eval.df, file = "03_results/concordance_eval.txt", sep = "\t", row.names = T, col.names = NA)
write.table(x = dat.df, file = "03_results/concordance_calls.txt", sep = "\t", row.names = T, col.names = NA)

# end
