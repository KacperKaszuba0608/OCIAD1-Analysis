Zip file contains the following files:

EigenMS.Rtest_EigenMS.pdftest_EigenMS.Rmdtest_peptides.txt
EigenMS_test1.R
README.txt


October 2017
Removed an error in normalizion function dealing with a single treatment group. 
Retested single and multiple group normalization after correction. 
Added test_EigenMS.Rmd which produces the “test_EigenMS.pdf” file that provides 
a few simple examples of running EigenMS. 


March 2015
This version removed the rescaling (addition of white noise) component of 
EigenMS. Rescaling in some cases produced rows of NAs under some 
missingness patterns. Some argue that rescaling is unnecessary as we remove the systematic bias which was added due to some technical artefact. Thus there is not need to add any additional white noise… Metabolomics data has very high variation and in our experiments no rescaling needed to be performed after EigenMS normalisation. In proteomic experiments there was also variation present after normalisation with the EigenMS-estimated number of bias trends. Users should be careful when changing the automatically estimated number of bias trends to be eliminated. Increasing the number will increase the chance of over-normalizing and producing identical values in all samples for a treatment group. 


