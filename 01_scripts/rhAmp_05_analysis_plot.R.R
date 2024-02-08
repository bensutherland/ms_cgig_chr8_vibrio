#Get yourself setup
setwd("C:/Users/ereth/Desktop/BIOL_491/rhAmpprocessing")
getwd()
install.packages("ggplot2")
library(ggplot2)
install.packages("tidyr")
library(tidyr)
install.packages("dplyr")
library(dplyr)
#import data
data.df <- read.csv(file = "rhAmp_full_results.csv")
head(data.df)

data.separate <- separate(data.df, col = Sample, into = c("Sample_family", "day_of_death"), sep = "_D")
View(data.separate)
#Remove NA day of death (controls)
data_subset <- subset(data.separate, day_of_death!= "NA")
#remove no-geno calls
data_subset2 <- subset(data_subset, geno!="no.geno")

#gogole rename function for renaming 6_alive
#geom jitter (points over box or violin plot)
#remove no.geno
#put line between 6 and 6_alive to indicate survival, 
#maybe include survival plot -> would have to summarize based on genotype for each day 
Jitter<- ggplot(data= data_subset2)+ geom_violin(aes(x=geno, y= day_of_death)) + geom_jitter(aes(x= geno, y= day_of_death))
Jitter

#violin plot
violin_allfamilies <- ggplot(data = data_subset2, aes(x = geno, y = day_of_death, group = geno, fill = geno, colour = geno)) + geom_violin(aes(x = geno, y = day_of_death, group = geno, fill = geno, colour = geno)) + geom_jitter(aes(x= geno, y= day_of_death))
violin_allfamilies

head(data_subset2)

data_families <- separate(data_subset2, col = Sample_family, into = c("Project", "family"), sep = "_F")
head(data_families)
#Separate Individual families
data_subset_114 <-subset(data_families, family == "114")
data_subset_115 <-subset(data_families, family == "115")
data_subset_116 <-subset(data_families, family == "116")
data_subset_117 <-subset(data_families, family == "117")
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
