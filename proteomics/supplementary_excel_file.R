# INPUT: 
# 1. Two files after log2 fold change calculation (MITOS and TOTALS)
# 2. raw data (proteinGroups.txt)
# 3. Results of GO enrichment (4 files: 2 up and 2 down)
# 4. Organelles comparison

# OUTPUT: supplementary file in excel format

library(dplyr)

organelles <- readr::read_csv('./data/organelles.csv')

fc_output_M <- readr::read_csv('./data/cleaned/OCIAD1_proteomics_mitos_process.csv') |> 
    select(-contains('MitoCarta'), -contains('imputed'))
fc_output_M <- merge(fc_output_M, organelles, by=c('Protein IDs', 'Gene names'), all.x = TRUE) |>
    distinct()

fc_output_T <- readr::read_csv('./data/cleaned/OCIAD1_proteomics_totals_process.csv') |> 
    select(-contains('MitoCarta'), -contains('imputed'))
fc_output_T <- merge(fc_output_T, organelles, by=c('Protein IDs', 'Gene names'), all.x = TRUE) |>
    distinct()

raw_data <- readr::read_tsv('./data/proteinGroups.txt') |>
    filter(is.na(`Only identified by site`),
           is.na(Reverse),
           is.na(`Potential contaminant`)) |>
    select(`Protein IDs`, `Gene names`, contains('LFQ intensity') 
           & (ends_with('22') | ends_with('23') | ends_with('24')),
           `Peptide sequences`)

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
                                       'FC_MITOS/TOTALS_OCIAD1',
                                       'p_OCIAD1_MITOS/TOTALS',
                                       'sig_MITOS/TOTALS_OCIAD1',
                                       'LFQ_WT(KO)_MITOS(TOTALS)_(22,23,24)',
                                       'missingness',
                                       'MitoCarta3.0',
                                       'Peroxisome DB',
                                       'One GO Mitochondrion',
                                       'One GO Peroxisome',
                                       'One GO ER',
                                       '###################',
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
                                       'Names of genes this protein is associated with.',
                                       'Names of proteins this protein is associated with.',
                                       'Log2 Fold Change based on means of normalized LFQ intenisties.',
                                       'p-value from t-test.',
                                       'Siginificance of protein (p < 0.05 & log2(FC) > 1).',
                                       'Imputed and normalized LFQ intensity values for proteins.',
                                       'Type of missigness: compelete (complete), missing at random (MAR), missing nor at random (MNAR), worse p-value (NA)',
                                       'Logical column containing if protein is in MitoCarta3.0 data base (TRUE/FALSE)',
                                       'Assigned organelle name based on the PeroxisomeDB. (TRUE/FALSE)',
                                       'Assigned mitochondrion using gene onthology for mitochondrion GO ID (TRUE/FALSE)',
                                       'Assigned mitochondrion using gene onthology for peroxisome GO ID (TRUE/FALSE)',
                                       'Assigned mitochondrion using gene onthology for endoplasmatic reticulum GO ID (TRUE/FALSE)',
                                       '###################',
                                       'Category such as BP (biological process), CC (cellular component), MF (molecular function).',
                                       'Unique GO identifier.',
                                       'Description of the GO term.',
                                       'The ratio of genes in the input dataset associated with the term.',
                                       'The ratio of genes in the entire genome associated with the term.',
                                       'Statistical significance of the enrichment.',
                                       'The corrected p-value.',
                                       'The p-value corrected using the False Discovery Rate approach.',
                                       'Specific genes in the input dataset associated with the term.',
                                       'The number of the genes in the input dataset associated with the term.')
                     )

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