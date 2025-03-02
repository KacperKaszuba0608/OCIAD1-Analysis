# INPUT: preprocessed data 'LipidGroups_processed_....csv'
# OUTPUT: boxplot of lipid classes

# Generating the boxplot showing different lipid classes

library(tidyverse)
library(plotly)
library(ggrepel)


#Loading data (specify which file to use, MITOS or TOTALS)
lipid.groups <- read_csv('LipidGroups_processed_....csv')

#Filtering the data to exclude the lipid classes with less than 25 data points
result <- filter_classes(lipid.groups, 25)
classes_filt <- result$filtered_classes
lipid.groups_filt <- result$filtered_df

#Filtering out unidentified lipids
lipid.groups_filt <- filter(lipid.groups_filt, `Lipid Class` != 'Unidentified')

#Palette describing the lipid category
palette <- c("Cholesteryl Ester" = "#48A0CF", "Fatty Acyl" = "#EC6B63", "Glycerolipid" = "#965DA6", "Phospholipid" = "#62C19B", "Sphingolipid" = "#FAC812", "Unidentified" = "#BEBEBE", "Other" = "#7FDFFD")

#Boxplot
box_lipidclass <- ggplot(lipid.groups_filt, aes(x = reorder(`Lipid Class`, FC, FUN = function(x) -median(x)), y = FC)) + 
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  stat_boxplot(geom ='errorbar') + 
  geom_boxplot(alpha = 0, col = '#4B4B4B', outlier.shape = NA) +
  geom_jitter(aes(color = `Lipid Category`, alpha = as.factor(Significant)), width = 0.2, size = 1.5, shape = 16) +
  ylab('Fold Change KO/WT') +
  theme_bw() +
  theme(axis.title.y = element_blank(), legend.position="none", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  guides(alpha = FALSE) +
  scale_alpha_manual(values = c(0.37, 0.8)) +
  scale_color_manual(values = palette) +
  scale_y_continuous(position="right") +
  coord_flip()

box_lipidclass