# INPUT: preprocessed data 'LipidGroups_processed_....csv'
# OUTPUT: boxplot containing information about fatty acids

# Generating the boxplot containing information about ether and even/odd lipids

library(tidyverse)
library(plotly)
library(ggrepel)


#Loading data (specify which file to use, MITOS or TOTALS)
lipid.groups <- read_csv('LipidGroups_processed_....csv')

#Combining information about ether and even/odd chains
lipid.groups$Interaction_Ether_EvenOdd <- factor(interaction(lipid.groups$Ether, lipid.groups$`Length Even`), 
                                                 labels = c('Odd\nNon-ether', 'Odd\nEther', 'Even\nNon-ether', 'Even\nEther'))

#Palette describing the lipid category
palette <- c("Cholesteryl Ester" = "#48A0CF", "Fatty Acyl" = "#EC6B63", "Glycerolipid" = "#965DA6", "Phospholipid" = "#62C19B", "Sphingolipid" = "#FAC812", "Unidentified" = "#BEBEBE", "Other" = "#7FDFFD")

#Boxplot (coded for only phospholipids, which can be modified)
box_fattyacids <- ggplot(filter(lipid.groups, `Lipid Category` == 'Phospholipid'), aes(x = reorder(`Interaction_Ether_EvenOdd_Sat`, FC, FUN = function(x) - median(x)), y = `FC`)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  stat_boxplot(geom ='errorbar', position = position_dodge(width = 0.9)) + 
  geom_boxplot(mapping = NULL, data = NULL, alpha = 0, col = '#4B4B4B', outlier.shape = NA, position = position_dodge(width = 0.9)) +
  geom_point(data = filter(lipid.groups, Significant == 0, `Lipid Category` == 'Phospholipid'), aes(color = `Lipid Category`), 
             alpha = 0.2, size = 2, shape = 16, position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.9)) +
  geom_point(data = filter(lipid.groups, Significant == 1, `Lipid Category` == 'Phospholipid'), aes(color = `Lipid Category`), 
             alpha = 0.8, size = 2, shape = 16, position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.9)) +
  ylab("Fold Change KO/WT") +
  scale_color_manual(values = palette) +
  scale_alpha_manual(values = c(0.37, 0.8)) +
  scale_y_continuous(position="right") +
  theme_bw() +
  theme(legend.position="none", axis.title.y = element_blank()) +
  scale_size(range = c(0.5, 7)) +
  coord_flip()

box_fattyacids