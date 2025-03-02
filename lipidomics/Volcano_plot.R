# INPUT: preprocessed data 'LipidGroups_processed_....csv'
# OUTPUT: volcano plot

# Generating the volcano plot

library(tidyverse)
library(plotly)
library(ggrepel)


#Loading data (specify which file to use, MITOS or TOTALS)
lipid.groups <- read_csv('LipidGroups_processed_....csv')

#Specifying which lipids to label on the plot based on their significance
FC_cutoff = 2.5
log10_pval_cutoff = 2

lipid.groups <- lipid.groups %>%
  mutate(Label = case_when(
    -log10(pvalue) > log10_pval_cutoff & abs(FC) > FC_cutoff & is.character(Identification) ~ as.character(Identification),
    TRUE ~ ''
  ))

#Palette describing the lipid category
palette <- c("Cholesteryl Ester" = "#48A0CF", "Fatty Acyl" = "#EC6B63", "Glycerolipid" = "#965DA6", "Phospholipid" = "#62C19B", "Sphingolipid" = "#FAC812", "Unidentified" = "#BEBEBE", "Other" = "#7FDFFD")

#Volcano plot
v1 <- ggplot(lipid.groups, aes(x = FC, y = -log10(pval), color = `Lipid Category`, shape = Significant&abs(FC)>1)) +
  geom_hline(yintercept = -log10(0.05), color = 'lightgrey', alpha = 0.6, linetype = 2) +
  geom_vline(xintercept = 1, color = 'lightgrey', alpha = 0.6, linetype = 2) +
  geom_vline(xintercept = -1, color = 'lightgrey', alpha = 0.6, linetype = 2) +
  geom_point(data = filter(lipid.groups, `Lipid Class` == 'Unidentified'), stroke = 1.2, alpha = 0.5) +
  geom_point(data = filter(lipid.groups, `Lipid Class` != 'Unidentified'), stroke = 1.2, alpha = 0.6) +
  xlab("Fold Change KO/WT") +
  ylab("-log10(p-value)") +
  theme_bw() +
  theme(legend.position = c(0.97, 0.03), legend.justification = c("right", "bottom"), legend.title=element_blank(), 
        legend.text = element_text(size=9), legend.background = element_blank(), legend.key = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(shape = 'none', alpha = 'none') +
  scale_shape_manual(values = c(1, 16)) +
  xlim(-max(abs(lipid.groups$FC)), max(abs(lipid.groups$FC))) +
  scale_color_manual(values = palette) +
  geom_text_repel(aes(label = Label), max.overlaps = Inf, verbose = TRUE, min.segment.length = 0, color = '#4B4B4B', box.padding = 1)

v1
