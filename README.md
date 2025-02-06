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

In this repository you will find the statistical analysis of the distributions of the
mass spectrometry samples. To see the result you should clone this repository
to your own local machine, extract the `EigenMS.zip` folder and then run files.

## Proteomics analysis

1. Firstly run files `imputation_MITOS.R` and `imputation_TOTALS.R`.
2. Next is normalization, so run files `normalization_MITOS.R` and `normalization_TOTALS.R`
3. The last process is fold change calculation. To achive that you have to run files `fold_change_calculation_MITOS.R` and `fold_change_calculation_TOTALS.R`
4. The data are ready to plotting with file `plotting.R`.

Good luck! 😉

## Lipidomics analysis

Open and run the file `Lipid_analysis_final.Rmd`.

Enjoy!

Authors of code: [Vanessa Linke](https://github.com/vanilink), [Mateusz Chodkowski](https://github.com/matiich), [Kacper Kaszuba](https://github.com/KacperKaszuba0608)