# NORMALIZATION with EigenMS for TOTALS

library(tidyverse)
source('./proteomics/code/EigenMS/EigenMS/EigenMS.R')

# Load imputed data
protein.groups <- readr::read_tsv('./data/cleaned/proteinGroupsImputed_TOTALS.txt',
                                  show_col_types = FALSE)

# Preparing data to normalization
lfq <- protein.groups |>
    select(starts_with('LFQ Intensity') 
           & (ends_with('22') | ends_with('23') | ends_with('24')) 
           & contains('TOTALS'))

colnames(lfq) <- gsub('_TOTALS_', '_', gsub('LFQ intensity ', '', colnames(lfq)))
lfq[lfq == 0] <- NA

# Performing normalization using EigenMS method
treatment = as.factor(c('KO', 'KO', 'KO', 'WT', 'WT', 'WT'))
prot.info <- data.frame(prot_ID = paste('prot_', 1:nrow(lfq), sep = '')) # add info about LFQ value
lfq.eig1 <- eig_norm1(lfq, 
                      treatment = treatment, 
                      prot.info = prot.info)

# Performing eig normalization
lfq.eig_norm <- eig_norm2(lfq.eig1)

lfq.eigen <- as.data.frame(lfq.eig_norm$norm_m)

# Export results
protein.groups.export <- protein.groups

protein.groups.export$`LFQ intensity KO_TOTALS_22` <- lfq.eigen$KO_22
protein.groups.export$`LFQ intensity KO_TOTALS_23` <- lfq.eigen$KO_23
protein.groups.export$`LFQ intensity KO_TOTALS_24` <- lfq.eigen$KO_24
protein.groups.export$`LFQ intensity WT_TOTALS_22` <- lfq.eigen$WT_22
protein.groups.export$`LFQ intensity WT_TOTALS_23` <- lfq.eigen$WT_23
protein.groups.export$`LFQ intensity WT_TOTALS_24` <- lfq.eigen$WT_24

write.table(protein.groups.export, file='./data/cleaned/proteinGroupsNormalized_TOTALS.txt',
            row.names = FALSE, sep = '\t')
