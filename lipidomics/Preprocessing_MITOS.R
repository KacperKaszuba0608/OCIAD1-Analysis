# INPUT: data from 'LipidGroups_intensities_MITOS.csv'
# OUTPUT: data frame with preprocessed data

# Preprocessing of the data: calculating the fold change, statistical significance, and adding all necessary labels

library(tidyverse)
source('./lipidomics/fxn.R')


#Loading data
lipid.groups <- read_csv('LipidGroups_intensities_MITOS.csv')

#Extracting just the numerical intensity data (removing the first three columns containing labels)
lipid.groups_mitos <- lipid.groups[-1:-3]

#Normalizing the data
lipid.groups_mitos <- normalise(lipid.groups_mitos)

#Calculating the fold change
lipid.groups_mitos_FC <- fc(lipid.groups_mitos, 3)

#Performing the t-test
lipid.groups_mitos_stats <- ttest(lipid.groups_mitos, 3)

#Combining the fold change and t-test data into one data frame
lipid.groups_mitos_processed <- add_column(lipid.groups_mitos_FC, 
                                 pvalue = lipid.groups_mitos_stats$pvalue, 
                                 Significant = lipid.groups_mitos_stats$Significant)

#Adding back the columns with labels
lipid.groups_mitos_processed <- lipid_groups_mitos_processed |> 
  add_column(Identification = lipid.groups$Identification, .before = 1) |>
  add_column(`Lipid Class` = lipid.groups$`Lipid Class`, .before = 1) |>
  add_column(`Unique ID` = lipid.groups$`Unique ID`, .before = 1)

#Categorizing the lipids 
lipid.groups_mitos_processed <- categorize_lipids(lipid.groups_mitos_processed)

#Export data
write.csv(lipid.groups_mitos_processed, file = './data/LipidGroups_processed_MITOS.csv',
          row.names = FALSE)