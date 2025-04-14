# OCIAD1 KO Analysis :dna:

## Abstract

This repository contains the code of the statisticial and data analysis of the data from the article 
"Integrated proteome and lipidome analyses place OCIAD1 at mitochondria-peroxisome intersection balancing lipid metabolism" (https://doi.org/10.1242/jcs.263729).

## Repository content &#x1F333; &#x1F4C1;

```
│   README.md
│
├───lipidomics
│       Data_extraction.R
│       Fatty_acids_boxplot.R
│       Fatty_acid_analysis_MITOS.R
│       Fatty_acid_analysis_TOTALS.R
│       fxn.R
│       Lipid_class_boxplot.R
│       Preprocessing_MITOS.R
│       Preprocessing_TOTALS.R
│       Supplementary_excel_file.R
│       Volcano_plot.R
│
└───proteomics
        EigenMS.zip
        fold_change_calculation_MITOS.R
        fold_change_calculation_TOTALS.R
        fxn.R
        GO_enrichment.R
        imputation_MITOS.R
        imputation_TOTALS.R
        normalization_MITOS.R
        normalization_TOTALS.R
        organelles_comparison.R
        peroxisomeDB.R
        pipeline.R
        plotting.R
        split_violin_and_statistical_analysis.R
        supplementary_excel_file.R
```

To see the result you should clone this repository to your own local machine, 
extract the `EigenMS.zip` folder, create some folders like `data` and inside of 
it a `cleaned` and a `GO_enrichment` folder, set the working directory to `OCIAD1-Analysis` and then 
run the files. 

## Proteomics analysis

We can break down proteomics analysis into the following steps:
1. Data preparation and calculation of all necessary measures.
2. Visualisation of the results.

### 1st step (Preparing Data)

1. Run the `imputation_MITOS.R` and `imputation_TOTALS.R` files first. The imputation has been done using the Ludovic method from [`protti`](https://doi.org/10.1093/bioadv/vbab041) library. The input for both code files is the `OCIAD1_prot_raw_data` sheet from the supplementary file.
2. Next is the normalization, so run the files `normalization_MITOS.R` and `normalization_TOTALS.R`. The normalization was performed using the EigenMS method.
3. The last process is the fold change calculation. For this you have to run the files `fold_change_calculation_MITOS.R` and `fold_change_calculation_TOTALS.R`.
4. Now you can also analyse the data in terms of gene ontology (GO term analysis). To do this, run the `GO_enrichment.R` file.

### 2nd step (Visualizing the results)

1. The data is ready for plotting using the `plotting.R` file.
2. To plot the split violin and also perform statistical analysis, run the `split_violin_and_statistical_analysis.R` file.
3. To create a supplementary file with all the data, run the `supplementary_excel_file.R` file. Note: 
Before you running this file, you need to create a file with organelles comparison. To do this run `organelles_comparison.R`.

Another way to run the analysis is to run the `pipeline.R` file and most of the analysis will be done except plotting files.

Good luck! 😉

## Lipidomics analysis

Similarly, the lipidomics analysis can be broken up into into the following steps:
1. Data preparation and calculation of all necessary measures.
2. Visualisation of the results.

### 1st step (Preparing Data)

1. The `Data_extraction.R` file takes in the raw data and extracts the information we need and separates it for MITOS and TOTALS.
2. The `Preprocessing_MITOS.R` and `Preprocessing_TOTALS.R` files perform the normalisation, fold change calculation, and statistical testing of the data.
3. The `Fatty_acids_MITOS.R` and `Fatty_acids_TOTALS.R` files extract useful information about fatty acid chains, such as the chain lengths and number of double bonds, from the lipid identifiers.

### 2nd step (Visualizing the results)

1. The data is ready for plotting using the `Volcano_plot.R`, `Lipid_class_boxplot.R`, and `Fatty_acids_boxplot.R` files.
2. The `Lipid_class_boxplot.R` creates a boxplot and scatter plot combination of fold changes for each lipid class.
3. The `Fatty_acids_boxplot.R` creates a boxplot and scatter plot combination of fold changes for each combination of fatty acid chain characteristics, even/odd length of the chains and whether they are ether lipids.
4. A supplementary file with all the data can be generated using the `Supplementary_excel_file.R` file.

Enjoy!

Authors of code: [Vanessa Linke](https://github.com/vanilink), [Mateusz Chodkowski](https://github.com/matiich), [Kacper Kaszuba](https://github.com/KacperKaszuba0608)

## References

1) PMID: 19602524. "Normalization of peak intensities in bottom-up MS-based proteomics using singular value decomposition".
Karpievitch YV, Taverner T, Adkins JN, Callister SJ, Anderson GA, Smith RD, Dabney AR.
Bioinformatics 2009

2) "Metabolomics data normalization with EigenMS"
Karpievitch YV, Nikolic SB, Wilson R, Sharman JE, Edwards LM.
PLoS One 2014

3) Jan-Philipp Quast, Dina Schuster, Paola Picotti. protti: an R package for comprehensive 
data analysis of peptide- and protein-centric bottom-up proteomics data. Bioinformatics 
Advances, Volume 2, Issue 1, 2022, vbab041, https://doi.org/10.1093/bioadv/vbab041

4) Rath, S.,* Sharma, R.*, Gupta, R.*, ..., Calvo, S.E., Mootha, V.K.. MitoCarta3.0: 
an updated inventory of the mitochondrial proteome, now with sub-organelle localization
and pathway annotations (2020). Nucleic Acids Research [Pubmed: 33174596](https://pubmed.ncbi.nlm.nih.gov/33174596/)

5) Schluter, A., Real-Chicharro, A., Gabaldon, T., Sanchez-Jimenez, F. and Pujol, A. 
(2010) PeroxisomeDB 2.0: An integrative view of the global peroxisomal metabolome. 
Nucleic Acids Res, 38, D800-5. doi: https://doi.org/10.1093/nar/gkp935
