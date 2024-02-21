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
install.packages("gridExtra")
library(gridExtra)
install.packages("gridGraphics")
library(gridGraphics)
install.packages("cowplot")
library(cowplot)
# create an empty list to store the plots
families <- c(114, 115, 116, 117)
boxplots <- list()
violins <- list()

# set the layout of the graphics device for the box plots
par(mfrow = c(2, 2))

for (f in families) {
  # subset the data for each family
  rhamp_family <- subset(rhamp_families, family == f)
  
  # save and plot the box plot for each family 
  png(paste0("boxplot_", f, ".png"), width = 6, height = 4, units = "in", res = 300) 

 
  boxplot <- ggplot(data = rhamp_family) + geom_boxplot(aes(x = majority.geno, y = day_of_death, fill = majority.geno)) + labs(x = "Genotype", y = "Day of Death") + theme_classic() + xlab("Genotype") + ylab("Day of Death") + labs(fill = "Genotype", colour = "Genotype") + scale_y_continuous(labels = c("3", "4", "5", "6", "Survivors")) + annotate("text", x = Inf, y = Inf, label = paste0("F", f), hjust = 1.4, vjust = 1, size = 5, fontface = "bold") 
  
  print(boxplot)
  graphics.off() 
  boxplots[[f]] <- boxplot
  
  # save and plot the violin plot for each family
  png(paste0("violin_", f, ".png"), width = 6, height = 4, units = "in", res = 300) 
  rhamp.violin <- ggplot(data= rhamp_family)+ geom_violin(aes(x=majority.geno, y= day_of_death, group = majority.geno, fill = majority.geno, colour = majority.geno)) + theme_classic() + xlab("Genotype") + ylab("Day of Death") + labs(fill = "Genotype", colour = "Genotype") + scale_y_continuous(labels = c("3", "4", "5", "6", "Survivors")) + annotate("text", x = Inf, y = Inf, label = paste0("F", f), hjust = 1.4, vjust = 1, size = 5, fontface = "bold")
  print(rhamp.violin) 
  graphics.off() 
  
  # add the violin plot to the list
  violins[[f]] <- rhamp.violin
}

#Install purr
if (!require(purrr)) {
  install.packages("purrr")
  library(purrr)
}
# remove NULL values from the lists
boxplots <- purrr::compact(boxplots)
violins <- purrr::compact(violins)

#Combine box plots into one
tiff("combined_boxplots.tiff", width = 11, height = 9, units = "in", res = 300)
grid.arrange(grobs = boxplots, nrow = 2, ncol = 2)
graphics.off()

# Combine all violin plots into one plot and save it as a tiff file
tiff("combined_violins.tiff", width = 11, height = 9, units = "in", res = 300)
grid.arrange(grobs = violins, nrow = 2, ncol = 2)
graphics.off()


#######Determine proportion of each type of genotype per family 
#Make new data frame
rhamp.prop <- as.data.frame(rhamp_families)
#Calculate proportions 
rhamp.prop <- rhamp.prop %>%
  group_by(family, majority.geno) %>%
  summarise(count = n()) %>%
  mutate(prop = count / sum(count))
#Make bar plot
geno.prop <- ggplot(rhamp.prop, aes(x = family, y = prop, fill = majority.geno)) + geom_bar(stat = "identity", position = "dodge") + geom_text(aes(label = paste0(round(prop*100, 0), "%")), position = position_dodge(width = 0.9), vjust = -0.25) + scale_y_continuous(labels = scales::percent, limits = c(0, 1)) + ylab("Proportion in Family") + xlab("Family") + labs(fill = "Genotype") + theme_classic()
geno.prop

######Make bar plot of %Survival and % Mortality for all families including non-mapping + control 
mortalitybarplot <- read.csv("mortalitybarplot.csv")
#Convert % to proportions
mortalitybarplot$survival <- mortalitybarplot$survival / 100
mortalitybarplot$death <- mortalitybarplot$death / 100

#Switch F118 to Control
mortalitybarplot$Family[mortalitybarplot$Family == "118"] <- "Control"

#convert to long format
mort_long <- reshape2::melt(mortalitybarplot, id.vars = "Family", variable.name = "status", value.name = "percent")

#Relabel family as chr
mort_long$Family <- as.character(mort_long$Family)

#Make plot
mort_bar <- ggplot(mort_long, aes(x = Family, y = percent, fill = status)) + geom_bar(stat = "identity") + labs(x = "Family", y = "Proportion", fill = "Status") + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_discrete(name = "Status", labels = c("Survival(%)", "Mortality(%)"))

mort_bar





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
