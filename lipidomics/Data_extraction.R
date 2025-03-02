# INPUT: raw data from `LipidGroups.csv`
# OUTPUT: data frames with the extracted intensities and lipid labels for MITOS and TOTALS

# Extracting lipid intensities from the raw data

library(tidyverse)
source('./lipidomics/fxn.R')


#Loading data
lipid.groups <- read_csv('LipidGroups.csv')

#Sorting the data alphabetically, flagging unidentified lipids, and creating a unique identifier for each measured lipid
lipid.groups <- lipid.groups |> 
  arrange(Identification) |> 
  mutate(`Lipid Class` = ifelse(is.na(`Lipid Class`), 'Unidentified', `Lipid Class`),
         `Unique ID` = paste0('ID_', 'RT', `Retention Time (min)`, '_', 'mz', `Quant Ion`, '_', `Polarity`)) %>%
  select(`Unique ID`, everything())


#Separately extracting all MITOS and TOTALS measurements
result <- get_sets(lipid.groups, 'MITO', 'TOT')
lipid.groups_mitos <- result$columns[[1]]
lipid.groups_totals <- result$columns[[2]]


#Grouping and rearranging the knockout and wild-type measurements in the MITOS data
result <- get_sets(lipid.groups_mitos, 'KO', 'WT')
lipid.groups_mitos_ko_names <- result$column_names[[1]]
lipid.groups_mitos_wt_names <- result$column_names[[2]]
lipid.groups_mitos <- select(lipid.groups_mitos, c(lipid.groups_mitos_ko_names, lipid.groups_mitos_wt_names))

#Grouping and rearranging the knockout and wild-type measurements in the TOTALS data
result <- get_sets(lipid.groups_totals, 'KO', 'WT')
lipid.groups_totals_ko_names <- result$column_names[[1]]
lipid.groups_totals_wt_names <- result$column_names[[2]]
lipid.groups_totals <- select(lipid.groups_totals, c(lipid.groups_totals_ko_names, lipid.groups_totals_wt_names))


#Adding back the columns with labels to MITOS data
lipid.groups_mitos <- lipid_groups_mitos |> 
  add_column(Identification = lipid.groups$Identification, .before = 1) |>
  add_column(`Lipid Class` = lipid.groups$`Lipid Class`, .before = 1) |>
  add_column(`Unique ID` = lipid.groups$`Unique ID`, .before = 1)

#Adding back the columns with labels to TOTALS data
lipid.groups_totals <- lipid_groups_totals |> 
  add_column(Identification = lipid.groups$Identification, .before = 1) |>
  add_column(`Lipid Class` = lipid.groups$`Lipid Class`, .before = 1) |>
  add_column(`Unique ID` = lipid.groups$`Unique ID`, .before = 1)

#Export data
write.csv(lipid.groups_mitos, file = './data/LipidGroups_intensities_MITOS.csv',
          row.names = FALSE)
write.csv(lipid.groups_totals, file = './data/LipidGroups_intensities_TOTALS.csv',
          row.names = FALSE)