# INPUT: preprocessed data 'LipidGroups_processed_MITOS.csv'
# OUTPUT: data frame with fatty acid chain information

# Extracting fatty acid information about identified lipids

library(tidyverse)
source('./lipidomics/fxn.R')


#Loading data
lipid.groups <- read_csv('LipidGroups_processed_MITOS.csv')

#Removing unidentified lipids
lipid.groups_mitos_notna <- lipid.groups |>
  drop_na(Identification)

#Extracting fatty acid information
lipid.groups_mitos_fattyacids <- fatty_acids(lipid.groups_mitos_notna)

#Export data
write.csv(lipid.groups_mitos_fattyacids, file = './data/LipidGroups_fattyacids_MITOS.csv',
          row.names = FALSE)
