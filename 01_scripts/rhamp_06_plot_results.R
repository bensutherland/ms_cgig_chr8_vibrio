# Plot genotype and mortality results
# Liam Surry, Ben Sutherland
# Initialized 2024-02-07

#### 00. Front matter ####
# Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

# Install and load packages
#install.packages("ggplot2")
#install.packages("tidyr")
#install.packages("dplyr")
#install.packages("gridExtra")
#install.packages("gridGraphics")
#install.packages("cowplot")
#install.packages("ggpubrr")

library(ggplot2)
library(tidyr)
library(dplyr)
library(gridGraphics)
library(gridExtra)
library(cowplot)
library(ggpubr)

# Set user variables
sample_id_interp.FN <- "00_archive/sample_day_of_death_and_DNA_ID.csv"
rhamp_to_geno_summary.FN <- "03_results/rhAmp_per_sample_genotype_summary.csv"


#### 01. Load and prep data ####
sample_id_interp.df <- read.csv(file = sample_id_interp.FN)
sample_id_interp.df <- as.data.frame(sample_id_interp.df)
head(sample_id_interp.df)
dim(sample_id_interp.df)

rhamp_to_geno_summary.df <- read.csv(file = rhamp_to_geno_summary.FN)
rhamp_to_geno_summary.df <- as.data.frame(rhamp_to_geno_summary.df)
head(rhamp_to_geno_summary.df)
dim(rhamp_to_geno_summary.df)

# Merge input files
rhamp_merged.df <- merge(x = sample_id_interp.df, y = rhamp_to_geno_summary.df, by = "samples")
dim(rhamp_merged.df)
head(rhamp_merged.df)

# Split sample info column (family.ID) into component parts (note: assumes _D separates)
rhamp_merged.df <- separate(data = rhamp_merged.df, col = "family.ID"
                             , into = c("sample_family", "day_of_death")
                             , sep = "_D"
                           )
head(rhamp_merged.df)
unique(rhamp_merged.df$day_of_death)
# Update alive to a arbitrary later-date of 7
rhamp_merged.df$day_of_death <- gsub(pattern = "6_ALIVE", replacement = "7", x = rhamp_merged.df$day_of_death)
table(rhamp_merged.df$day_of_death)

# Assign day_of_death values to numerical 
rhamp_merged.df$day_of_death <- as.numeric(rhamp_merged.df$day_of_death)

# table(paste0(rhamp_merged.df$sample_family, "_", rhamp_merged.df$day_of_death))

# Ensure all samples have genotypes
unique(rhamp_merged.df$majority.geno)


#### 02. Plot day of death by genotype per-family ####
# Change genotypes to number of alternate alleles
rhamp_merged.df$num_alt <- rhamp_merged.df$majority.geno
rhamp_merged.df$num_alt <- gsub(pattern = "homo.ref", replacement = "0", x = rhamp_merged.df$num_alt)
rhamp_merged.df$num_alt <- gsub(pattern = "het", replacement = "1", x = rhamp_merged.df$num_alt)
rhamp_merged.df$num_alt <- gsub(pattern = "homo.alt", replacement = "2", x = rhamp_merged.df$num_alt)
head(rhamp_merged.df)
str(rhamp_merged.df)
rhamp_merged.df$num_alt <- as.numeric(rhamp_merged.df$num_alt)

# Subset families
head(rhamp_merged.df)

# Separate data frame into families 
rhamp_merged.df  <- separate(rhamp_merged.df, col = sample_family, into = c("project", "family"), sep = "_F", remove = F)

#add a new column "mortality" to each data frame in the families list, encoding “died” as 1 and “survived” as 0 based on the day_of_death column.
rhamp_merged.df <- rhamp_merged.df %>%
  mutate(mortality = ifelse(`day_of_death` %in% 3:6, 1, 0))

#Subset families
rhamp_family_114 <- subset(rhamp_merged.df, family == "114")
rhamp_family_115 <- subset(rhamp_merged.df, family == "115")
rhamp_family_116 <- subset(rhamp_merged.df, family == "116")
rhamp_family_117 <- subset(rhamp_merged.df, family == "117")


# Plot with ggplot2
F114.plot <-  ggplot(rhamp_family_114, aes(x = as.factor(num_alt), y= day_of_death)) + 
                  geom_boxplot(fill="gray") +
                  geom_jitter() + 
                  labs(x= "Number alternate alleles", y = "Days to death") +
                  theme(axis.text = element_text(size = 12), 
                        axis.title = element_text(size = 14)) + 
                  theme_classic()
                
F115.plot  <-   ggplot(rhamp_family_115, aes(x = as.factor(num_alt), y= day_of_death)) + 
                  geom_boxplot(fill="gray") +
                  geom_jitter() + 
                  labs(x= "Number alternate alleles", y = "Days to death") +
                  theme(axis.text = element_text(size = 12), 
                        axis.title = element_text(size = 14))+ 
  theme_classic()
                
F116.plot  <- ggplot(rhamp_family_116, aes(x = as.factor(num_alt), y= day_of_death)) + 
                  geom_boxplot(fill="gray") +
                  geom_jitter() + 
                  labs(x= "Number alternate alleles", y = "Days to death") +
                  theme(axis.text = element_text(size = 12), 
                        axis.title = element_text(size = 14))+ 
  theme_classic()
                
F117.plot   <-  ggplot(rhamp_family_117, aes(x = as.factor(num_alt), y= day_of_death)) + 
                  geom_boxplot(fill="gray") +
                  geom_jitter() + 
                  labs(x= "Number alternate alleles", y = "Days to death") +
                  theme(axis.text = element_text(size = 12), 
                        axis.title = element_text(size = 14))+ 
  theme_classic()


                
final_fig <- ggarrange(F114.plot, F115.plot, F116.plot, F117.plot
                , labels = c("A", "B", "C", "D")
                , ncol = 2, nrow = 2)

pdf(file = "03_results/boxplot_days_to_death_by_num_alt_alleles.pdf", width = 8, height = 5.5)
print(final_fig)
dev.off()

# Base R boxplot and linear model stats
par(mfrow=c(1,1))
boxplot(rhamp_merged.df$day_of_death[grep(pattern = "F114", x = rhamp_merged.df$sample_family)] 
        ~ rhamp_merged.df$num_alt[grep(pattern = "F114", x = rhamp_merged.df$sample_family)]
        , ylab = "day to death", xlab = "Number alt. allele"
        , las = 1
        )

mod <- lm(formula = rhamp_merged.df$day_of_death[grep(pattern = "F114", x = rhamp_merged.df$sample_family)] 
                  ~ rhamp_merged.df$num_alt[grep(pattern = "F114", x = rhamp_merged.df$sample_family)])
summary(mod)


boxplot(rhamp_merged.df$day_of_death[grep(pattern = "F115", x = rhamp_merged.df$sample_family)] 
        ~ rhamp_merged.df$num_alt[grep(pattern = "F115", x = rhamp_merged.df$sample_family)]
        , ylab = "day to death", xlab = "Number alt. allele"
        , las = 1
        
)
mod <- lm(formula = rhamp_merged.df$day_of_death[grep(pattern = "F115", x = rhamp_merged.df$sample_family)] 
          ~ rhamp_merged.df$num_alt[grep(pattern = "F115", x = rhamp_merged.df$sample_family)])
summary(mod)


boxplot(rhamp_merged.df$day_of_death[grep(pattern = "F116", x = rhamp_merged.df$sample_family)] 
        ~ rhamp_merged.df$num_alt[grep(pattern = "F116", x = rhamp_merged.df$sample_family)]
        , ylab = "day to death", xlab = "Number alt. allele"
        , las = 1
        
)
mod <- lm(formula = rhamp_merged.df$day_of_death[grep(pattern = "F116", x = rhamp_merged.df$sample_family)] 
          ~ rhamp_merged.df$num_alt[grep(pattern = "F116", x = rhamp_merged.df$sample_family)])
summary(mod)

mod <- aov(formula = rhamp_merged.df$day_of_death[grep(pattern = "F116", x = rhamp_merged.df$sample_family)] 
          ~ as.factor(rhamp_merged.df$num_alt[grep(pattern = "F116", x = rhamp_merged.df$sample_family)]))
summary(mod)

boxplot(rhamp_merged.df$day_of_death[grep(pattern = "F117", x = rhamp_merged.df$sample_family)] 
        ~ rhamp_merged.df$num_alt[grep(pattern = "F117", x = rhamp_merged.df$sample_family)]
        , ylab = "day to death", xlab = "Number alt. allele"
        , las = 1
        
)

mod <- lm(formula = rhamp_merged.df$day_of_death[grep(pattern = "F117", x = rhamp_merged.df$sample_family)] 
          ~ rhamp_merged.df$num_alt[grep(pattern = "F117", x = rhamp_merged.df$sample_family)])
summary(mod)


#### 03. Mortality (dead or alive) as function of # individuals versus number of alternate alleles#### 

# Make list of all families
families <- list(rhamp_family_114, rhamp_family_115, rhamp_family_116, rhamp_family_117)
names(families) <- c("114", "115", "116", "117")

#set colour pallet
cbbPalette <- c("#CCCCCC", "#444444")

#create empty list for storing plots
plot_list <- list()

# Loop and create bar plot for all families: mortality (dead versus alive) versus number of alt alleles. Assumes your data frame is named rhamp_family_114, rhamp_family_115, etc.
for (name in c("114", "115", "116", "117")) {
 
  mort_barplot.df <- get(paste0("rhamp_family_", name))
  
  # Convert num_alt and mortality to factors
  mort_barplot.df$num_alt <- as.factor(mort_barplot.df$num_alt)
  mort_barplot.df$mortality <- as.factor(mort_barplot.df$mortality)
  
  # Group and summarize number of died/survivors within each family
  count_data <- mort_barplot.df %>%
    group_by(num_alt, mortality) %>%
    summarise(count = n())
  
  # Convert Mortality back to character for the plot
  count_data$mortality <- ifelse(count_data$mortality == "1", "Died", "Survived")
  
  # Create bar plot
  p <- ggplot(count_data, aes(x = num_alt, y = count, fill = mortality)) +
        geom_bar(stat = "identity", position = position_dodge(width = 1.0), colour = "black") +
        labs(x = "Number of Alternate Alleles", y = "# Individuals", fill = "Mortality") +
        theme_classic() + theme(legend.position = "none") +
        scale_fill_manual(values = cbbPalette) + theme(axis.text = element_text(size= 12),axis.title = element_text(size = 14)) 
  #if(name!="115"){
    
    #p <- p + theme(legend.position = "none")
    
  #}
  
  # Add the plot to the list
  plot_list[[name]] <- p
}

# Combine barplots in list into one figure
final_bar <- ggarrange(plotlist = plot_list, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)

#Save as PDF
pdf(file = "03_results/barplot_#individuals_mortality_by_num_alt_alleles.pdf", width = 8, height = 5.5)
tiff(file = "03_results/barplot_#individuals_mortality_by_num_alt_alleles.tiff", width = 8, height = 5.5, units = "in", res = 300)
print(final_bar)
dev.off()
print(final_bar)


####04. Determine proportion of each type of genotype per family ####
#Make new combined data frame
rhamp.prop <- rbind(rhamp_family_114,rhamp_family_115,rhamp_family_116,rhamp_family_117)

#Make order of genotypes homo.ref, het, homo.alt
rhamp.prop$majority.geno <- factor(rhamp.prop$majority.geno, levels = c("homo.ref", "het", "homo.alt"))
                                 
#Calculate proportions of survivors versus morts
rhamp.prop <- rhamp.prop %>%
  group_by(family, majority.geno) %>%
  summarise(count = n()) %>%
  mutate(prop = count / sum(count))

#Make bar plot for percentage of each genotype in mapping families. Relabel genotypes
geno.prop <- ggplot(rhamp.prop, aes(x = family, y = prop, fill = majority.geno)) + 
  geom_bar(stat = "identity", position = "dodge", colour = "black") + 
  geom_text(aes(label = paste0(round(prop*100, 0), "%")), position = position_dodge(width = 0.9), vjust = -0.25) + 
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) + ylab("Percentage in Family") +
  xlab("Family") + labs(fill = "Genotype") + 
  theme_classic() + theme(axis.text = element_text(size= 12),axis.title = element_text(size = 14)) + scale_fill_manual(values = c("homo.ref" = "white", "het" = "#CCCCCC", "homo.alt" = "#444444"), labels = c("homo.ref" = "ref/ref", "het" = "ref/alt", "homo.alt" = "alt/alt"))    

geno.prop

#Save as PDF
pdf(file = "03_results/proportion_genotype_mapping_family.pdf", width = 8, height = 5.5)
tiff(file = "03_results/proportion_genotype_mapping_family.tiff", width = 8, height = 5.5, units = "in", res = 300)
print(geno.prop)
dev.off()


####05.Make bar plot of %Survival and % Mortality for all families including non-mapping + control#### 
mortalitybarplot <- read.csv("00_archive/mortalitybarplot.csv")
family_parental_genotypes <- read.csv("02_input_data/OSU_MBP_parental_crosses.csv")

#Convert % to proportions
mortalitybarplot$survival <- mortalitybarplot$survival / 100
mortalitybarplot$death <- mortalitybarplot$death / 100

#Switch F118 to Control
mortalitybarplot$family[mortalitybarplot$family == "118"] <- "Control"

#convert to long format
mort_long <- reshape2::melt(mortalitybarplot, id.vars = "family", variable.name = "status", value.name = "percent")

#Relabel family as chr
mort_long$family <- as.character(mort_long$family)

#Set Colour pallete
cbbPalette <- c("#CCCCCC", "#444444")

#Make plot
mort_bar <- ggplot(mort_long, aes(x = family, y = percent, fill = status)) + geom_bar(stat = "identity", colour = "black") + labs(x = "Family", y = "Proportion", fill = "Status") + theme_classic() +
  scale_fill_manual(values = cbbPalette) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none") 

#Save as PDF
pdf(file = "03_results/barplot_prop_individuals_mortality_per_family.pdf", width = 8, height = 5.5)
tiff(file = "03_results/barplot_prop_individuals_mortality_per_family.tiff", width = 8, height = 5.5, units = "in", res = 300)
print(mort_bar)
dev.off()
print(mort_bar)


####06.Make bar plot of %Survival and % Mortality grouped by family genotypes####

#merge mortality proportions and family parental crosses by genotype. Requires file "OSU_MBP_parental_crosses
mort_proportions_and_family_crosses <- merge(mortalitybarplot, family_parental_genotypes, by = "family")

#make so genotype cross orders consistently regardless of which parent.
  order_function <- function(parent_1_genotype, parent_2_genotype) {
    if (parent_1_genotype < parent_2_genotype) {
    return(paste(parent_1_genotype, parent_2_genotype, sep = " x "))
    } else {
    return(paste(parent_2_genotype, parent_1_genotype, sep = " x "))
    }
  }

#create genotype cross column
mort_proportions_and_family_crosses$parental_cross <- mapply(order_function, mort_proportions_and_family_crosses$parent_1_genotype, mort_proportions_and_family_crosses$parent_2_genotype)

#calculate average mortality per genotype cross
average_mortality <- mort_proportions_and_family_crosses %>%
  group_by(parental_cross) %>%
  summarise(average_death = mean(death), average_survival = mean(survival))

#Convert average mortality to long format
average_mortality <- reshape2::melt(average_mortality, id.vars = "parental_cross", variable.name = "status", value.name = "percent")

#convert to long format
mort_proportions_and_family_crosses <- reshape2::melt(mort_proportions_and_family_crosses, 
id.vars = c("family", "parental_cross"), variable.name = "status", value.name = "percent",measure.vars = c("survival", "death"))

# Reorder parental cross genotypes 
reorder_genotype <- c("(ref/ref) x (ref/ref)", "(alt/alt) x (ref/ref)", "(ref/alt) x (ref/alt)")
mort_proportions_and_family_crosses$parental_cross <- factor(mort_proportions_and_family_crosses$parental_cross, levels = reorder_genotype)

#Set colours                            
cbbPalette <- c("death" = "#CCCCCC", "survival" = "#444444")  

#Create bar plot grouped by genotype cross 
mort_bar_family_genotype <- ggplot(mort_proportions_and_family_crosses, aes(x = family, y = percent, fill = status)) + geom_bar(stat = "identity", position = "stack", colour = "black") + labs(x = "Family", y = "Proportion") + theme_classic() + scale_fill_manual(values = cbbPalette) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none") + facet_grid(~ parental_cross, scales = "free_x", space = "free_x")

pdf(file = "03_results/barplot_prop_individuals_mort_parental_cross_x_families.pdf", width = 8, height = 5.5)
tiff(file = "03_results/barplot_prop_individuals_mort_parental_cross_x_families.tiff", width = 8, height = 5.5, units = "in", res = 300)

print(mort_bar_family_genotype)
dev.off()
print(mort_bar_family_genotype)

#bar plot based on average mortality per parental cross 
cbbPalette <- c("average_survival" = "#CCCCCC", "average_death" = "#444444")

average_morts_parental_crosses_no_fam_barplot <- ggplot(average_mortality, aes(x = parental_cross, y = percent, fill = status)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(x = "Parental Cross", y = "Proportion") +
  theme_classic() +theme(legend.position = "none") + scale_fill_manual(values = cbbPalette)

pdf(file = "03_results/average_morts_parental_crosses_no_fam_barplot.pdf", width = 8, height = 5.5)
tiff(file = "03_results/average_morts_parental_crosses_no_fam_barplot.tiff", width = 8, height = 5.5, units = "in", res = 300)
print(average_morts_parental_crosses_no_fam_barplot)
dev.off()
print(average_morts_parental_crosses_no_fam_barplot)


####07. Basic violin and boxplot for genotype versus mortality####

#Looping to make box plot and violin plot for every family.

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
  tiff(paste0("boxplot_", f, ".tiff"), width =8.95275591, height = 10.102322, units = "in", res = 300) 
  
  
  boxplot <- ggplot(data = rhamp_family) + geom_boxplot(aes(x = majority.geno, y = day_of_death, fill = majority.geno)) + labs(x = "Genotype", y = "Day of Death") + theme_classic() + xlab("Genotype") + ylab("Day of Death") + labs(fill = "Genotype", colour = "Genotype") + scale_y_continuous(labels = c("3", "4", "5", "6", "Survivors")) + theme(axis.text = element_text(size = 28),axis.title = element_text(size = 22),legend.key.size = unit(4, 'cm'),legend.title = element_text(size = 22),legend.text = element_text(size = 22), legend.position = "none") + annotate("text", x = Inf, y = Inf, label = paste0("F", f), hjust = 1.4, vjust = 1, size = 15, fontface = "bold") 
  
  print(boxplot)
  graphics.off() 
  boxplots[[f]] <- boxplot
  
  # save and plot the violin plot for each family
  png(paste0("violin_", f, ".png"), width = 6, height = 4, units = "in", res = 300) 
  rhamp.violin <- ggplot(data= rhamp_family)+ geom_violin(aes(x=majority.geno, y= day_of_death, group = majority.geno, fill = majority.geno, colour = majority.geno)) + theme_classic() + xlab("Genotype") + ylab("Day of Death") + labs(fill = "Genotype", colour = "Genotype") + theme(axis.text = element_text(size = 18)) + scale_y_continuous(labels = c("3", "4", "5", "6", "Survivors")) + annotate("text", x = Inf, y = Inf, label = paste0("F", f), hjust = 1.4, vjust = 1, size = 15, fontface = "bold")
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
