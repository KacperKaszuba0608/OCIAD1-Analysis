suppressWarnings(dir.create('./data'))
suppressWarnings(dir.create('./data/cleaned'))
suppressWarnings(dir.create('./data/GO_enrichment'))
suppressWarnings(dir.create('./supplementary_files'))
suppressWarnings(dir.create('./plots'))

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

message('Done!')
