# This file contains 2 examples of running EigenMS 
# 
# test_peptides.txt dataset provides 4 samples with 2000 peptdies (can be metabolites) each
# Test 1 shows how to run EigenMS on data with 1 factor with 2 levels (can have more levels)
# Test 2 expands the same dataset by duplicating the 4 columns and creates another treatment
#        groupping. This is a fake dataset. 
# 
# There are 2 functions that need to be run to complete the normalization process:
#       eig_norm1 & eig_norm2
# Variables retuned by the 2 functions in EigenMS contain many subvariables.
# Normalized dta is stored in the matrix 'norm_m' returend by eig_norm2


# EigenMS example, data is already on LOG-SCALE, 
# otherwise replace 0's with NAs and take log base 2
ddata = read.table("test_peptides.txt", header=TRUE) 
dim(ddata) # 1946 peptides x 4 samples
m_logInts = ddata[,3:6]

# if not on log scale and if has 0's for missing values, uncomment the following 3 lines:
# m_logInts[m_logInts==0] = NA  #  3.4% missing values, remove 0's, replace with NA's
# m_logInts = log2(m_logInts)
m_nummiss = sum(is.na(m_logInts)) # 1946 total values, 681 missing values

# plot boxplots for each sample
par(mar=c(10,3,3,3)) # allows to have nice vertical labels!!! this provides room for them
par(mfcol=c(1,1))
boxplot(m_logInts, las=2) 

m_numtot = dim(m_logInts)[1] * dim(m_logInts)[2] # 8000 total observations
m_percmiss = m_nummiss/m_numtot  # 8.75% percent missing observations
# plot numbr of missing values for each sample
par(mfcol=c(1,1))
barplot(colSums(is.na(m_logInts)), main="Numbers of missing values in samples (grp order)") #570 x 600

# source in the EigenMS and correlation heatmap functions
source("EigenMS.R")

# define parameter prot.info, 2 column data frame with IDs for metabolites or peptides
# in case of matabolites the 2 columns are identical. 
# For peptides 1st column must contain unique peptide ID (usually sequences)
# 2nd column can contain protein IDs
m_prot.info = ddata[,1:2] # all unique metabolite/peptide IDs, otherwise not possible to have as row names


# Example part 1 - single factor
# set up treatent group information
grps = as.factor(c(1,2,1,2)) # 1 = control; 2 = treatment
m_ints_eig1 = eig_norm1(m=m_logInts,treatment=grps,prot.info=m_prot.info)
m_ints_eig1$h.c
m_ints_norm1 = eig_norm2(rv=m_ints_eig1) 
par(mfcol=c(1,1))
boxplot(m_ints_norm1$norm_m, las=2)


# TEST part 2 - provide 1 treatment grou for the same dataset
# Note that we will be removing most of the biological variation from the data. 
# Thus this is just a test to see if code to normalize 1 treatment groups works. 
grps2 = as.factor(c(1,1,1,1))
m_ints_eig1_v2 = eig_norm1(m=m_logInts,treatment=grps2,prot.info=m_prot.info)
m_ints_norm1_v2 = eig_norm2(rv=m_ints_eig1_v2) 
par(mfcol=c(1,1)) 
boxplot(m_ints_norm1_v2$norm_m, las=2) # less variation for S3
boxplot(m_ints_norm1$norm_m, las=2) 




# TEST part 3 - nested treatment groups
# make the matrix of data about 2 times larger by copying columns and adding a bit of variation
m_logInts2 = cbind(m_logInts, m_logInts+runif(2000, 0, .3))
gr1 = as.factor(c(1,2,1,2,1,2,1,2)) # similar to S1, S2, ... in biological data
gr2 = as.factor(c(1,1,1,1,2,2,2,2))
grps2 = data.frame(gr1,gr2)
par(mar=c(10,4,4,4)) # important allows to have nice vertical labels
par(mfcol=c(1,1))
boxplot(m_logInts2, las=2) # 700 x 700


m_ints_eig2 = eig_norm1(m=m_logInts2,treatment=grps2,prot.info=m_prot.info)
m_ints_norm2 = eig_norm2(rv=m_ints_eig2) 
par(mfcol=c(1,1))
boxplot(m_ints_norm2$norm_m, las=2)  # normalized 



