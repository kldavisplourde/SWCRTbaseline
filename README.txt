The ZIP folder contains files for implementing the sample size and power calculations under a bivariate, change score, and ANCOVA linear mixed effects models as shown in "Designing stepped wedge cluster randomized trials with a baseline measurement of the outcome" by Davis-Plourde, Goldfeld, Allore, Taljaard, and Li.

List of Files:
1) Application_TSOS.r = R file for generating the Application study.
2) calPower.r = R file for generating the power for the bivariate, change score, and ANCOVA linear mixed effects models under a SW-CRT design.
3) copri_test_HH.r = R file for conducting the Simulation study.
4) EM_uncorrected_BV.r = R code for EM Algorithm (used to fit bivariate linear mixed model).

NOTES:  1) Simulation study additionally requires the installation of the nlme, doMC, doRNG, and lmeInfo packages. The EM algorithm requires the installation of the mvtnorm, numDeriv, and sjlabelled packages.
	2) You will need to change path names before running the programs. 
	3) Latest version of all files are available on GitHub: https://github.com/kldavisplourde/SWCRTbaseline