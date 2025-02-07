# OCIAD1 Analysis :dna:

Repository content &#x1F333; &#x1F4C1;

```
│   README.md
│
├───data
│   │   Human.MitoCarta3.0.csv
│   │   proteinGroups.txt
│   │
│   └───cleaned
├───lipidomics
│       Lipid_analysis_final.Rmd
│
└───proteomics
    └───code
        │   EigenMS.zip
        │   fold_change_calculation_MITOS.R
        │   fold_change_calculation_TOTALS.R
        │   fxn.R
        │   imputation_MITOS.R
        │   imputation_TOTALS.R
        │   normalization_MITOS.R
        │   normalization_TOTALS.R
        │
        └───EigenMS
```

This repository contains the statistical analysis of the distributions of the mass spectrometry samples. 
To see the result you should clone this repository to your own local machine, extract the EigenMS.zip 
folder, set the working directory to OCIAD1-Analysis and then run the files.

## Proteomics analysis

1. Run the `imputation_MITOS.R` and `imputation_TOTALS.R` files first.
2. Next is the normalization, so run the files `normalization_MITOS.R` and `normalization_TOTALS.R`.
3. The last process is the fold change calculation. For this you have to run the files `fold_change_calculation_MITOS.R` and `fold_change_calculation_TOTALS.R`.
4. The data are ready for plotting with the file `plotting.R`.

Good luck! 😉

## Lipidomics analysis

Open and run the `Lipid_analysis_final.Rmd` file.

Enjoy!

Authors of code: [Vanessa Linke](https://github.com/vanilink), [Mateusz Chodkowski](https://github.com/matiich), [Kacper Kaszuba](https://github.com/KacperKaszuba0608)
