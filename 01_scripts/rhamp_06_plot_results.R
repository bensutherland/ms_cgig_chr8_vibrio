#Setting Up
setwd("C:/Users/ereth/Desktop/BIOL_491/rhAmpprocessing/03_results/")
getwd()
install.packages("ggplot2")
library(ggplot2)
install.packages("tidyr")
library(tidyr)
install.packages("dplyr")
library(dplyr)

#import data
rhamp.df <- read.csv(file = "sample_day_of_death_and_DNA_ID.csv")
rhamp_to_geno_summary <- read.csv(file = "rhAmp_per_sample_genotype_summary.csv")
#merge data frames
rhamp_merged <- merge(rhamp.df, rhamp_to_geno_summary, by = "samples")

#Split sample column into family and day of death
rhamp_merged <- separate(rhamp_merged, col = family.ID, into = c("Sample_family", "day_of_death"), sep = "_D")
head(rhamp_merged)


#table(data.separate$day_of_death)
#table(paste0(data.separate$Sample_family, "_", data.separate$day_of_death))


#Remove no genos
#rhamp_separate_day[2,8]
rhamp_no_geno <- rhamp_merged[!is.na(rhamp_merged$day_of_death) & rhamp_merged$majority.geno!="no.geno",]
dim(rhamp_merged)

#Reassign day_of_death call "6_ALIVE as arbitrary day of death 7 
rhamp_no_geno$day_of_death[rhamp_no_geno$day_of_death == "6_ALIVE"] <- 7
head(rhamp_no_geno)

#Assign day_of_death values to numerical 
rhamp_no_geno$day_of_death <- as.numeric(as.character(rhamp_no_geno$day_of_death))

#Plotting: basic boxplot
boxplot(rhamp_no_geno$day_of_death ~ rhamp_no_geno$majority.geno, xlab = "Genotype", ylab = "Day_of_Death")

#Plotting without survivors 
rhamp_no_survivors<- rhamp_no_geno[!is.na(rhamp_no_geno$day_of_death) & rhamp_no_geno$day_of_death !="7",]
boxplot(rhamp_no_survivors$day_of_death ~ rhamp_no_survivors$majority.geno, xlab = "Genotype", ylab = "Day_of_Death(no survivors)")
#Violin plot representation 
rhamp.violin.no.jitter<- ggplot(data= rhamp_no_geno)+ geom_violin(aes(x=majority.geno, y= day_of_death, group = majority.geno, fill = majority.geno, colour = majority.geno)) + theme_classic() + xlab("Genotype") + ylab("Day of Death") + labs(fill = "Genotype", colour = "Genotype") + scale_y_continuous(labels = c("3", "4", "5", "6", "Survivors"))
rhamp.violin.no.jitter
#Now add jitter + make transparent
rhamp.violin.jitter<- ggplot(data= rhamp_no_geno)+ geom_violin(aes(x=majority.geno, y= day_of_death, group = majority.geno, fill = majority.geno, colour = majority.geno), alpha = 0.2) + geom_jitter(aes(x=majority.geno, y= day_of_death, group = majority.geno, fill = majority.geno, colour = majority.geno)) + theme_classic() + xlab("Genotype") + ylab("Day of Death") + labs(fill = "Genotype", colour = "Genotype") + scale_y_continuous(labels = c("3", "4", "5", "6", "Survivors"), breaks = c(3, 4, 5, 6, 7))
rhamp.violin.jitter

#add in line at day 6 for survivors 
rhamp.violin.no.jitter.line<- ggplot(data= rhamp_no_geno)+ geom_violin(aes(x=majority.geno, y= day_of_death, group = majority.geno, fill = majority.geno, colour = majority.geno)) + theme_classic() + xlab("Genotype") + ylab("Day of Death") + labs(fill = "Genotype", colour = "Genotype") + geom_hline(aes(yintercept = 6.01), colour = "black")+ annotate(geom="text", x=2, y=6.3, label="Survivors", color = "black")  + scale_y_continuous(labels = c("3", "4", "5", "6", "Survivors"))
rhamp.violin.no.jitter.line
#Per family plotting 
#Separate data frame into families 
rhamp_families <- separate(rhamp_no_geno, col = Sample_family, into = c("Project", "family"), sep = "_F")
head(rhamp_families)
#Separate into indiivudal families:
rhamp_family_114 <-subset(rhamp_families, family == "114")
rhamp_family_115 <-subset(rhamp_families, family == "115")
rhamp_family_116 <-subset(rhamp_families, family == "116")
rhamp_family_117 <-subset(rhamp_families, family == "117")

#Looping to make box plot and violin plot for every family.
families <- c(114, 115, 116, 117)
for (f in families) {
  rhamp_family <- subset(rhamp_families, family == f)
  boxplot(rhamp_no_geno$day_of_death ~ rhamp_no_geno$majority.geno, xlab = "Genotype", ylab = "Day_of_Death")
  tiff(paste0("boxplot_", f, ".tiff"), compression = "lzw")
  def.off()
  
  rhamp.violin <- ggplot(data= rhamp_family)+ geom_violin(aes(x=majority.geno, y= day_of_death, group = majority.geno, fill = majority.geno, colour = majority.geno)) + theme_classic() + xlab("Genotype") + ylab("Day of Death") + labs(fill = "Genotype", colour = "Genotype") + scale_y_continuous(labels = c("3", "4", "5", "6", "Survivors"))
  
  tiff(paste0("violin_", f, ".tiff"), compression = "lzw")
  def.off()
}


rhamp.violin.no.jitter<- ggplot(data= rhamp_no_geno)+ geom_violin(aes(x=majority.geno, y= day_of_death, group = majority.geno, fill = majority.geno, colour = majority.geno)) + theme_classic() + xlab("Genotype") + ylab("Day of Death") + labs(fill = "Genotype", colour = "Genotype") + scale_y_continuous(labels = c("3", "4", "5", "6", "Survivors"))

#Skeleton Code
####Remove NA day of death (controls)
####data_subset <- subset(data.separate, day_of_death!= "NA")
####remove no-geno calls
####data_subset2 <- subset(data_subset, geno!="no.geno")
#gogole rename function for renaming 6_alive
#geom jitter (points over box or violin plot)
#remove no.geno
#put line between 6 and 6_alive to indicate survival, 
#maybe include survival plot -> would have to summarize based on genotype for each day 



#violin plot
violin_allfamilies <- ggplot(data = rhamp.separate, aes(x = geno, y = day_of_death, group = geno, fill = geno, colour = geno)) + geom_violin(aes(x = geno, y = day_of_death, group = geno, fill = geno, colour = geno)) + geom_jitter(aes(x= geno, y= day_of_death))
violin_allfamilies

head(data_subset2)

data_families <- separate(data_subset2, col = Sample_family, into = c("Project", "family"), sep = "_F")
head(data_families)
#Separate Individual families

#Individual violin plots
violin_F114 <- ggplot(data = data_subset_114, aes(x = geno, y = day_of_death, group = geno, fill = geno, colour = geno)) + geom_violin(aes(x = geno, y = day_of_death, group = geno, fill = geno, colour = geno)) + geom_jitter(aes(x= geno, y= day_of_death))
violin_F114
violin_F115 <- ggplot(data = data_subset_115, aes(x = geno, y = day_of_death, group = geno, fill = geno, colour = geno)) + geom_violin(aes(x = geno, y = day_of_death, group = geno, fill = geno, colour = geno)) + geom_jitter(aes(x= geno, y= day_of_death))
violin_F115
violin_F116 <- ggplot(data = data_subset_116, aes(x = geno, y = day_of_death, group = geno, fill = geno, colour = geno)) + geom_violin(aes(x = geno, y = day_of_death, group = geno, fill = geno, colour = geno)) + geom_jitter(aes(x= geno, y= day_of_death))
violin_F116
violin_F117 <- ggplot(data = data_subset_117, aes(x = geno, y = day_of_death, group = geno, fill = geno, colour = geno)) + geom_violin(aes(x = geno, y = day_of_death, group = geno, fill = geno, colour = geno))+ geom_jitter(aes(x= geno, y= day_of_death))
violin_F117

#boxplot2.0
boxplottest <- boxplot(data_subset2)
