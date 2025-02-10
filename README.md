# OCIAD1 Analysis :dna:

Repository content &#x1F333; &#x1F4C1;

```
│   GO_vs_peroxisomeDB.csv
│   README.md
│
├───data
│   │   Human.MitoCarta3.0.csv
│   │   peroxisomeDB.csv
│   │   proteinGroups.txt
│   │
│   ├───cleaned
│   └───GO_enrichment
│
├───lipidomics
│       Lipid_analysis_final.Rmd
│
├───proteomics
│   └───code
│       │   EigenMS.zip
│       │   fold_change_calculation_MITOS.R
│       │   fold_change_calculation_TOTALS.R
│       │   fxn.R - self-written functions
│       │   GO_enrichment.R
│       │   imputation_MITOS.R - imputation using the Ludovic method
│       │   imputation_TOTALS.R - imputation using the Ludovic method
│       │   normalization_MITOS.R - normalization using the EigenMS method
│       │   normalization_TOTALS.R - normalization using the EigenMS method
│       │   peroxisomeDB.R
│       │   plotting.R
│       │   supplementary_excel_file.R
│       │
│       └───EigenMS
│
└───supplementary_files
```

This repository contains the statistical analysis of the distributions of the mass spectrometry samples. 
To see the result you should clone this repository to your own local machine, extract the `EigenMS.zip` 
folder, set the working directory to `OCIAD1-Analysis` and then run the files.

## Proteomics analysis

1. Run the `imputation_MITOS.R` and `imputation_TOTALS.R` files first.
2. Next is the normalization, so run the files `normalization_MITOS.R` and `normalization_TOTALS.R`.
3. The last process is the fold change calculation. For this you have to run the files `fold_change_calculation_MITOS.R` and `fold_change_calculation_TOTALS.R`.
4. Now you can also analyse the data in terms of gene ontology (GO term analysis). To do this, run the `GO_enrichment.R` file.
5. The data is ready for plotting using the `plotting.R` file.
6. You can also create a supplementary file with all the data using the `supplementary_excel_file.R` file.

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
