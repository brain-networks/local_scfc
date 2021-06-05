# local_scfc
Calculate structure-function coupling at a local level as in Esfahlani et al (2021).

# dependencies and requirements
This script uses MATLAB. It was tested using MATLAB 2020a on macOS Mojave 10.14.6. It requires no additional software packages.

# installation and usage
Download and unzip the directory, open MATLAB, and then navigate to wherever you downloaded the files. You can run the script from within that directory.   

# what is here
* Group-representative sc/fc/euclidean distance in 400-node parcellation.
* Directory containing functions for generating predictors based on sc matrix.
* Example script to implement the procedure.

# what does the script do?
1. Reads in sc/fc/euclidean distance data
2. Generates a set of predictors
3. Calculates local (regional) correlations of predictors with FC.

# if you want to apply this method to your data
* Replace variables sc, fc, d with your own structural connectivity, functional connectivity, and inter-regional Euclidean distance matrix.

Runtime for entire script is < 10 s.

If you use this code, please cite:
Esfahlani, F. Z., Faskowitz, J., Slack, J., Misic, B., & Betzel, R. (2021). Local structure-function relationships in human brain networks across the human lifespan. bioRxiv. [(link to paper)](https://www.biorxiv.org/content/10.1101/2021.05.23.445128v1.abstract)
