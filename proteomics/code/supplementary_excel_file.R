# INPUT: 
# 1. Two files after log2 fold change calculation (MITOS and TOTALS)
# 2. raw data (proteinGroups.txt)
# 3. Results of GO enrichment (4 files: 2 up and 2 down)

# OUTPUT: supplementary file in excel format

library(dplyr)

fc_output_M <- readr::read_csv('./data/cleaned/OCIAD1_proteomics_mitos_process.csv')
fc_output_T <- readr::read_csv('./data/cleaned/OCIAD1_proteomics_totals_process.csv')
raw_data <- readr::read_tsv('./data/proteinGroups.txt') |>
    select(`Protein IDs`, `Gene names`, contains('LFQ intensity') 
           & (ends_with('22') | ends_with('23') | ends_with('24')))

df_go_mitos_up <- readr::read_csv("./data/GO_enrichment/OCIAD1_proteomics_mitos_GOterms_up.csv", show_col_types = F) |>
    select(ONTOLOGY, ID, Description, GeneRatio, BgRatio, pvalue, 
           p.adjust, qvalue, geneID, Count)

df_go_mitos_down <- readr::read_csv("./data/GO_enrichment/OCIAD1_proteomics_mitos_GOterms_down.csv", show_col_types = F) |>
    select(ONTOLOGY, ID, Description, GeneRatio, BgRatio, pvalue, 
           p.adjust, qvalue, geneID, Count)

df_go_totals_up <- readr::read_csv("./data/GO_enrichment/OCIAD1_proteomics_totals_GOterms_up.csv", show_col_types = F) |>
    select(ONTOLOGY, ID, Description, GeneRatio, BgRatio, pvalue, 
           p.adjust, qvalue, geneID, Count)

df_go_totals_down <- readr::read_csv("./data/GO_enrichment/OCIAD1_proteomics_totals_GOterms_down.csv", show_col_types = F) |>
    select(ONTOLOGY, ID, Description, GeneRatio, BgRatio, pvalue, 
           p.adjust, qvalue, geneID, Count)

#### saving results ####
readme <- data.frame("COLUMN NAME" = c('Protein IDs',
                                       'Gene names',
                                       'Protein names',
                                       'MitoCarta3.0_SubMitoLocalization',
                                       'MitoCarta3.0_MitoPathways',
                                       'FC_MITOS/TOTALS_OCIAD1',
                                       'p_OCIAD1_MITOS/TOTALS',
                                       'sig_MITOS/TOTALS_OCIAD1',
                                       'LFQ_WT(KO)_MITOS(TOTALS)_(22,23,24)',
                                       'Ontology',
                                       'ID',
                                       'Description',
                                       'GeneRatio',
                                       'BgRatio',
                                       'pvalue',
                                       'p.adjust',
                                       'qvalue',
                                       'geneID',
                                       'Count'),
                     "DESCRIPTION" = c('Identifiers of proteins contained in the protein group.',
                                       'Names of genes this peptide is associated with.',
                                       'Names of proteins this peptide is associated with.',
                                       'Annotation to most likely compartment: matrix, inner membrane (MIM), intermembrane space (IMS), outer membrane (MOM), mitochondrial membrane, or unknown.',
                                       'Annotation to a hierarchy of 149 biological pathways, based on literature.',
                                       'Log2 Fold Change based on means of normalized LFQ intenisties.',
                                       'p-value from t-test',
                                       'Siginificance of peptide (p < 0.05 & log2(FC) > 1)',
                                       'Imputed and normalized LFQ intensity values for peptides.',
                                       'Category such as BP (biological process), CC (cellular component), MF (molecular function)',
                                       'Unique GO identifier',
                                       'Description of the GO term',
                                       'The ratio of genes in the input dataset associated with the term',
                                       'The ratio of genes in the entire genome associated with the term',
                                       'Statistical significance of the enrichment',
                                       'The corrected p-value',
                                       'The p-value corrected using the False Discovery Rate approach',
                                       'Specific genes in the input dataset associated with the term',
                                       'The number of the genes in the input dataset associated with the term'))

list_of_sheets <- list("README" = readme,
                       "OCIAD1_prot_mitos_process" = fc_output_M,
                       "OCIAD1_prot_totals_process" = fc_output_T,
                       "OCIAD1_prot_mitos_GOdown" = df_go_mitos_up,
                       "OCIAD1_prot_mitos_GOup" = df_go_mitos_down,
                       "OCIAD1_prot_totals_GOdown" = df_go_totals_up,
                       "OCIAD1_prot_totals_GOup" = df_go_totals_down,
                       "OCIAD1_prot_raw_data" = raw_data)

openxlsx::write.xlsx(list_of_sheets, file = "./supplementary_files/SuppFile2_Proteomics.xlsx",
                     keepNA=TRUE, na.string='NA')