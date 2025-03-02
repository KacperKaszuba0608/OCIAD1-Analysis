# INPUT: preprocessed data 'LipidGroups_processed_TOTALS.csv'
# OUTPUT: data frame with fatty acid chain information

# Extracting fatty acid information about identified lipids

library(tidyverse)
source('./lipidomics/fxn.R')


#Loading data
lipid.groups <- read_csv('LipidGroups_processed_TOTALS.csv')

#Removing unidentified lipids
lipid.groups_totals_notna <- lipid.groups |>
  drop_na(Identification)

#Extracting fatty acid information
lipid.groups_totals_fattyacids <- fatty_acids(lipid.groups_totals_notna)

#Export data
write.csv(lipid.groups_totals_fattyacids, file = './data/LipidGroups_fattyacids_TOTALS.csv',
          row.names = FALSE)
