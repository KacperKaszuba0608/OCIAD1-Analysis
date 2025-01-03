# OCIAD1 Analysis :dna:

In this repository you will find the statistical analysis of the distributions of the
mass spectrometry samples. To see the result you should clone this repository
to your own local machine, extract the `EigenMS.zip` folder and then run one of the following files:
1. `MITOS_norm.Rmd` for an analysis of an isolated mitochondria.
2. `TOTALS_norm.Rmd` for an analysis of a total cells.

When you ran this file the report will be generated and you will find there:
* ditributions of raw data and imputed data,
* tables with skewness for different transformations of data,
* distributions of normalized data with EigenMS method,
* statistics metrics like:
    * **CV** - Coefficient of Variation;
    * **ICC** - Intraclass Correlation Coefficient;
    * **unpaired t-test**,
* in the section `Raw VS Normalized` you will find the interactive plots

Enjoy! 😉
