# OCIAD1 Analysis :dna:

Repository content &#x1F333; &#x1F4C1;

```
│   README.md
│
├───lipidomics
│       Lipid_analysis_final.Rmd
│
└───proteomics
    │   fold_change_calculation_MITOS.R
    │   fold_change_calculation_TOTALS.R
    │   fxn.R
    │   GO_enrichment.R
    │   imputation_MITOS.R
    │   imputation_TOTALS.R
    │   normalization_MITOS.R
    │   normalization_TOTALS.R
    │   organelles_comparison.R
    │   peroxisomeDB.R
    │   plotting.R
    │   split_violin_and_statistical_analysis.R
    │   supplementary_excel_file.R
    │
    └───EigenMS.zip
```

This repository contains the statistical analysis of the distributions of the mass spectrometry samples. 
To see the result you should clone this repository to your own local machine, extract the `EigenMS.zip` 
folder, set the working directory to `OCIAD1-Analysis` and then run the files.

## Proteomics analysis

We can break down proteomics analysis into the following steps:
1. Data preparation and calculation of all necessary measures.
2. Visualisation of the results.

### 1st step (Preparing Data)

1. Run the `imputation_MITOS.R` and `imputation_TOTALS.R` files first. The imputation has been done using the Ludovic method from [`protti`](https://doi.org/10.1093/bioadv/vbab041) library.
2. Next is the normalization, so run the files `normalization_MITOS.R` and `normalization_TOTALS.R`. The normalization was performed using the EigenMS method.
3. The last process is the fold change calculation. For this you have to run the files `fold_change_calculation_MITOS.R` and `fold_change_calculation_TOTALS.R`.
4. Now you can also analyse the data in terms of gene ontology (GO term analysis). To do this, run the `GO_enrichment.R` file.

### 2nd step (Visualizing the results)

1. The data is ready for plotting using the `plotting.R` file.
2. To plot the split violin and also perform statistical analysis, run the `split_violin_and_statistical_analysis.R` file.
3. To create a supplementary file with all the data, run the `supplementary_excel_file.R` file. Note: 
Before you running this file, you need to create a file with organelles comparison. To do this run `organelles_comparison.R`.

Good luck! 😉

## Lipidomics analysis

Open and run the `Lipid_analysis_final.Rmd` file.

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