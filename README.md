# fracpred
This R package is to predict the first dimension fractions containing target peptides in a 2D LC-MS/MS proteomics analysis. 

To use the package, you first need to install the devtools package. 
install.packages("devtools")

Load the devtools package.
library(devtools)

Install the "fracpred" package from GitHub
install_github("zhuhaojie/fracpred")

Load the fracpred package.
library(fracpred)

Use the "fracpred" function to predict the rentention times of the peptides of interest
fracpred()

Please note that the CSV files should at least contain "Peptide Sequence" and "Peptide Retention Time" columns. 
The files can be exported from the Skyline software. 
