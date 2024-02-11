# Cosmic-distance-ladder-fits
This repository contains the files needed to make least-squares fit using information about the cosmic distance ladder method.

The data used are in the form of vectors in the folder "data-vectors", the theoretical equations of each model are in the form of matrices in the folder "equations-matrices" and the information about the covariance matrices for each case are in the folder "covariance-matrices". 

An example: To fit a data sample with non-standardized apparent magnitudes of supernovae and considering the standard color correction model, it is necessary use "Y_r22-pantheon+_unc.fits" (folder "data-vectors"), with "L_I.fits" (folder "equations-matrices"), and "C_r22-pantheon+_unc.fits" ("covariance-matrices") and "fit_I.py" (folder "fit-files").
