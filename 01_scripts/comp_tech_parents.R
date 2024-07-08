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
panel_records.FN <- "../compare_amp_panel_parents_and_wgrs_parents/0002.vcf"
wgrs_records.FN  <- "../compare_amp_panel_parents_and_wgrs_parents/0003.vcf"

# Read in data
panel.vcf <- read.vcfR(file = panel_records.FN)
wgrs.vcf  <- read.vcfR(file = wgrs_records.FN)

panel.vcf
wgrs.vcf

#TODO: this is where we could limit to only hotspot loci

panel.df <- extract.gt(panel.vcf)
wgrs.df  <- extract.gt(wgrs.vcf)

panel.df[1:5,1:5]
wgrs.df[1:5,1:5]

# Update colnames to make match
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

wgrs_indiv_names.df <- colnames(wgrs.df)
wgrs_indiv_names.df <- as.data.frame(wgrs_indiv_names.df)
colnames(wgrs_indiv_names.df) <- "indiv"
wgrs_indiv_names.df <- separate(data = wgrs_indiv_names.df, col = "indiv", into = c("fam", "ind", "suffix"), sep = "-", remove = T)
wgrs_indiv_names.df$full_name <- paste0(wgrs_indiv_names.df$fam, "-", wgrs_indiv_names.df$ind, "_g")
colnames(wgrs.df) <- wgrs_indiv_names.df$full_name


panel.df[1:5,]
wgrs.df[1:5,1:5]

# Combine the two dfs
all.df <- cbind(panel.df, wgrs.df)
dim(all.df)
all.df[1:5,1:5]
all.df <- all.df[, order(colnames(all.df))]
all.df[1:5,]


# Drop non-hotspot loci
all_hotspot.df <- all.df[grep(pattern = "_", x = rownames(all.df), invert = T), ]
all_hotspot.df[1:10,1:8]


dat.df <- all_hotspot.df

# Loop to build tally
eval.df <- NULL
eval.df$mname <- rownames(dat.df)
eval.df <- as.data.frame(eval.df)
head(eval.df)
eval.df$na.tally <- NA
eval.df$tally <- NA

# unique parents
unique_inds.vec <- unique(gsub(pattern = "_g|_p", replacement = "", x = colnames(dat.df)))

moi <- NULL
for(i in 1:length(eval.df$mname)){
  
  print(i)
  
  moi <- eval.df$mname[i]
  tally <- 0
  na.tally <- 0
  panel.call <- NULL
  genome.call <- NULL
  
  for(ind in 1:length(unique_inds.vec)){
    
    print(ind)
    
    panel.call  <- dat.df[moi, which(colnames(dat.df)==paste0(unique_inds.vec[ind], "_p"))]
    genome.call <- dat.df[moi, paste0(unique_inds.vec[ind], "_g")]
    
    if(!is.na(panel.call) & !is.na(genome.call)){
      
      if(panel.call==genome.call){
        
        tally <- tally + 1
        
      }else{
        
        tally <- tally
        
      }
      
    }else{
      
      na.tally <- na.tally + sum(is.na(c(panel.call, genome.call)))
      
    }
     
  }
  
  eval.df[i,"na.tally"] <- na.tally
  eval.df[i,"tally"]    <- tally
    
}
  
head(eval.df)

pdf(file = "03_results/histogram_of_concordance_tally.pdf", width = 6.5, height = 4.5)
hist(eval.df$tally)
dev.off()


write.table(x = eval.df, file = "03_results/concordance_eval.txt", sep = "\t", row.names = T, col.names = NA)
write.table(x = dat.df, file = "03_results/concordance_input_data.txt", sep = "\t", row.names = T, col.names = NA)

