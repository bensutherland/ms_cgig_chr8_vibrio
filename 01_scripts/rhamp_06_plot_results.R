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




head(rhamp_families)
#Reorder genotypes
rhamp_families$majority.geno <- factor(rhamp_families$majority.geno, levels = c ("homo.ref", "het", "homo.alt"))

rhamp_families <- rhamp_families %>%
  mutate(Mortality = ifelse(`day_of_death` %in% 3:6, 1, 0))

#Separate into indiivudal families:
rhamp_family_114 <-subset(rhamp_families, family == "114")
rhamp_family_115 <-subset(rhamp_families, family == "115")
rhamp_family_116 <-subset(rhamp_families, family == "116")
rhamp_family_117 <-subset(rhamp_families, family == "117")



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


#### Making a bar plot based on genotype and alive vs dead per family 


# List of all families
families <- list(rhamp_family_114, rhamp_family_115, rhamp_family_116, rhamp_family_117)
names(families) <- c("114", "115", "116", "117")

plot_list <- list()

# Loop through each family
for (name in names(families)) {
  # Convert majority.geno and Mortality to factor
  families[[name]]$majority.geno <- as.factor(families[[name]]$majority.geno)
  families[[name]]$Mortality <- as.factor(families[[name]]$Mortality)
  
  # Create count data
  count_data <- families[[name]] %>%
    group_by(majority.geno, Mortality) %>%
    summarise(Count = n())
  
  # Convert Mortality back to character for the plot
  count_data$Mortality <- ifelse(count_data$Mortality == "1", "Died", "Survived")
  
  # Create bar plot
  p <- ggplot(count_data, aes(x=majority.geno, y=Count, fill=Mortality)) +
    geom_bar(stat="identity", position=position_dodge()) +
    labs(x="Genotype", y="Count", fill="Mortality") +
    theme_classic() + ggtitle(paste("F", name, sep=""))
  
  # Add the plot to the list
  plot_list[[name]] <- p
}

# Combine all plots into one
combined_plot <- do.call(grid.arrange, c(plot_list, ncol=2))

# Save the combined plot as a TIFF file
ggsave("combined_plot.tiff", combined_plot, width = 11, height = 9)


#######Determine proportion of each type of genotype per family 
#Make new data frame
rhamp.prop <- as.data.frame(rhamp_families)
#Calculate proportions 
rhamp.prop <- rhamp.prop %>%
  group_by(family, majority.geno) %>%
  summarise(count = n()) %>%
  mutate(prop = count / sum(count))
#Make bar plot
geno.prop <- ggplot(rhamp.prop, aes(x = family, y = prop, fill = majority.geno)) + geom_bar(stat = "identity", position = "dodge") + geom_text(aes(label = paste0(round(prop*100, 0), "%")), position = position_dodge(width = 0.9), vjust = -0.25) + scale_y_continuous(labels = scales::percent, limits = c(0, 1)) + ylab("Percentage in Family") + xlab("Family") + labs(fill = "Genotype") + theme_classic() +  theme(axis.text = element_text(size = 28), axis.title = element_text(size = 22), legend.key.size = unit(2, 'cm'), legend.title = element_text(size = 22), legend.text = element_text(size = 22)) 
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


