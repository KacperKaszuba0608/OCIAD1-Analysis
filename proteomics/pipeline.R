dir.create('./data')
dir.create('./data/cleaned')
dir.create('./data/GO_enrichment')
dir.create('./supplementary_files')
dir.create('./plots')

# Totals processing
source('./proteomics/imputation_TOTALS.R')
source('./proteomics/normalization_TOTALS.R')
source('./proteomics/fold_change_calculation_TOTALS.R')

# Mitos processing
source('./proteomics/imputation_MITOS.R')
source('./proteomics/normalization_MITOS.R')
source('./proteomics/fold_change_calculation_MITOS.R')

source('./proteomics/GO_enrichment.R')
source('./proteomics/organelles_comparison.R')
source('./proteomics/supplementary_excel_file.R')
