# INPUT: 
# 1. raw data (LipidGroups.csv)
# 2. Two files after data extraction (MITOS and TOTALS)
# 3. Two files after preprocessing and fatty acid analysis (MITOS and TOTALS)

# OUTPUT: supplementary file in excel format

library(dplyr)

raw_data <- readr::read_csv('./data/LipidGroups.csv')

intensities_MITOS <- readr::read_csv('./data/LipidGroups_intensities_MITOS.csv')
intensities_TOTALS <- readr::read_csv('./data/LipidGroups_intensities_TOTALS.csv')

processed_MITOS <- readr::read_csv('./data/LipidGroups_fattyacids_MITOS.csv')
processed_TOTALS <- readr::read_csv('./data/LipidGroups_fattyacids_TOTALS.csv')


#### saving results ####
readme <- data.frame("COLUMN NAME" = c('Unique ID',
                                       'Lipid Category',
                                       'Lipid Class',
                                       'Identification',
                                       'VL24/23/22 OCIAD1 KO MITOS/TOTALS',
                                       'VL24/23/22 WT MITOS/TOTALS',
                                       'FC',
                                       'FC up/down',
                                       'pvalue',
                                       'Significant',
                                       'Chain Lengths',
                                       'Length Even',
                                       'Double Bonds',
                                       'Saturation',
                                       'Ether',
                                       ),
                     "DESCRIPTION" = c('Unique identifier for each row, containing information about whether the lipid was identified, its retention time, quant ion, and polarity',
                                       'Name of the wider category to which the lipid belongs',
                                       'Name of the class to which the lipid belongs',
                                       'Specific identification of the lipid, containing information about its head group and fatty acid chains',
                                       'Intensities of the lipids in the OCIAD1 knockout condition',
                                       'Intensities of the lipids in the wild-type condition',
                                       'Log2 fold change based on means of normalized intenisties',
                                       'Whether the change between the conditions (KO and WT) is positive or negative (1: positive change, 0: negative change)',
                                       'P-value from the t-test',
                                       'Whether the change between the conditions is statistically significant (p-value <0.05)',
                                       'Lengths of identified fatty acid chains (sometimes the lengths are combined into a single sum)',
                                       'Whether the length of at all one fatty acid chain is even (1: all fatty acid chains are even, 0: at least one fatty acid chain is odd)',
                                       'The number of double bonds in the identified fatty acid chains (sometimes the numbers are combined into a single sum)',
                                       'The type of saturation of the fatty acid chains (0: SFA (no double bonds), 1: MUFA (at least one chain has only one double bond), 2: PUFA (at least one chain has more than one double bond))',
                                       'Whether the lipid is an ether lipid'
                                       ))

list_of_sheets <- list('README' = readme,
                       'OCIAD1_lipid_raw_data' = raw_data,
                       'OCIAD1_lipid_intensities_MITOS' = intensities_MITOS,
                       'OCIAD1_lipid_intensities_TOTALS' = intensities_TOTALS,
                       'OCIAD1_lipid_processed_MITOS' = processed_MITOS,
                       'OCIAD1_lipid_processed_TOTALS' = processed_TOTALS
                       )

openxlsx::write.xlsx(list_of_sheets, file = "./supplementary_files/SuppFile_Lipidomics.xlsx",
                     keepNA=TRUE, na.string='NA')