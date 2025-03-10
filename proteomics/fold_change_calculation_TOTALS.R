# INPUT: a normalized data from `proteinGroupsNormalized_TOTALS.txt`
# OUTPUT: a data frame with log2 fold change values for TOTALS

# Calculating log 2 Fold Change for TOTALS

library(tidyverse)
source('./proteomics/fxn.R')

# load necessarry data
mitocarta <- read.csv('./data/Human.MitoCarta3.0.csv')
protein.groups <- read_tsv('./data/cleaned/proteinGroupsNormalized_TOTALS.txt',
                           show_col_types = FALSE)

metadata <- merge(protein.groups, mitocarta, by.x = 'Gene names', by.y = 'Symbol', all.x = TRUE)
metadata <- metadata |>
    select(`Protein IDs`, `Gene names`, `Protein names`, MitoCarta3.0_SubMitoLocalization, 
           MitoCarta3.0_MitoPathways, missingness, contains('imputed'))

# setting cut offs
cutoff.p = 0.05 # alpha for t-test
cutoff.FC = 1 # fold change cut-off

# extracting LFQ values
lfq <- protein.groups |>
    select(starts_with('LFQ Intensity') 
           & (ends_with('22') | ends_with('23') | ends_with('24')) 
           & contains('TOTALS'))

colnames(lfq) <- gsub(' intensity ', '_', colnames(lfq))
lfq$protein.ids <- protein.groups$`Protein IDs`

# calculate p-value with unpaired t-test
pvalue <- lfq |>
    mutate("p_OCIAD1_TOTALS" = apply(lfq, 1, ttest,
                                    grp1=grep("KO", colnames(lfq)), 
                                    grp2=grep("WT", colnames(lfq)))) |>
    mutate(sig_TOTALS_OCIAD1 = ifelse(p_OCIAD1_TOTALS < cutoff.p, TRUE, FALSE)) |>
    relocate(protein.ids, .before = 1)

# Histogram of p-values
ggplot(data = pvalue, aes(x = p_OCIAD1_TOTALS, fill=sig_TOTALS_OCIAD1)) + 
    geom_histogram(bins = 100) + ggtitle("TOTALS - p.value")

# Preaparing data to calculations of Fold Change
lfq_long <- lfq |>
    pivot_longer(1:6, names_to = "sample", values_to = "LFQvalue")

# Add a new column with groups KO or WT
lfq_long$group <- "TOTALS_OCIAD1_KO"
lfq_long$group[grep("WT",lfq_long$sample)] <- "TOTALS_OCIAD1_WT"

# Calculation of means for groups KO and WT
means <- lfq_long |>
    group_by(group, protein.ids) |>
    summarise(mean=mean(LFQvalue, na.rm=TRUE)) |>
    ungroup()

# Calculate Fold Change
means <- pivot_wider(means, names_from = group, values_from = mean) |>
    mutate(FC_TOTALS_OCIAD1 = TOTALS_OCIAD1_KO - TOTALS_OCIAD1_WT)

# Preparing data for volcano plot
FC <- merge(means, lfq, by = "protein.ids")

final_table <- merge(FC, pvalue[,c("protein.ids", "p_OCIAD1_TOTALS", "sig_TOTALS_OCIAD1")], by = "protein.ids") |>
    select(-TOTALS_OCIAD1_KO, -TOTALS_OCIAD1_WT) |>
    mutate(sig_TOTALS_OCIAD1 = ifelse(p_OCIAD1_TOTALS < cutoff.p & abs(FC_TOTALS_OCIAD1) > cutoff.FC, TRUE, FALSE))

final_table <- merge(metadata, final_table, by.y = 'protein.ids', by.x = 'Protein IDs') |>
    select(`Protein IDs`, `Gene names`, `Protein names`, MitoCarta3.0_SubMitoLocalization, 
           MitoCarta3.0_MitoPathways, FC_TOTALS_OCIAD1, p_OCIAD1_TOTALS, sig_TOTALS_OCIAD1,
           LFQ_KO_TOTALS_22, LFQ_KO_TOTALS_23, LFQ_KO_TOTALS_24, LFQ_WT_TOTALS_22, 
           LFQ_WT_TOTALS_23, LFQ_WT_TOTALS_24, missingness, contains('imputed'))

write.csv(final_table, file = './data/cleaned/OCIAD1_proteomics_totals_process.csv',
          row.names = FALSE)