# INPUT: data from 'LipidGroups_intensities_TOTALS.csv'
# OUTPUT: data frame with preprocessed data

# Preprocessing of the data: calculating the fold change, statistical significance, and adding all necessary labels

library(tidyverse)
source('./lipidomics/fxn.R')


#Loading data
lipid.groups <- read_csv('LipidGroups_intensities_TOTALS.csv')

#Extracting just the numerical intensity data (removing the first three columns containing labels)
lipid.groups_totals <- lipid.groups[-1:-3]

#Normalizing the data
lipid.groups_totals <- normalise(lipid.groups_totals)

#Calculating the fold change
lipid.groups_totals_FC <- fc(lipid.groups_totals, 3)

#Performing the t-test
lipid.groups_totals_stats <- ttest(lipid.groups_totals, 3)

#Combining the fold change and t-test data into one data frame
lipid.groups_totals_processed <- add_column(lipid.groups_totals_FC, 
                                            pvalue = lipid.groups_totals_stats$pvalue, 
                                            Significant = lipid.groups_totals_stats$Significant)

#Adding back the columns with labels
lipid.groups_totals_processed <- lipid_groups_totals_processed |> 
  add_column(Identification = lipid.groups$Identification, .before = 1) |>
  add_column(`Lipid Class` = lipid.groups$`Lipid Class`, .before = 1) |>
  add_column(`Unique ID` = lipid.groups$`Unique ID`, .before = 1)

#Categorizing the lipids 
lipid.groups_totals_processed <- categorize_lipids(lipid.groups_totals_processed)

#Export data
write.csv(lipid.groups_totals_processed, file = './data/LipidGroups_processed_TOTALS.csv',
          row.names = FALSE)
