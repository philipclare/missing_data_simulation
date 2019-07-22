# Comparison of methods of adjusting for time-varying confounding under misspecification – A Monte-Carlo simulation study
## Stata and R Analysis Code

This repository contains the Stata and R code used in the missing data simulation by Clare et al. 2019

The Stata code creates a series of quasi-random datasets (3 different datasets were used in the simulation) using a pre-specified data structure.
Analysis code runs all analyses on those datasets, and saves the results. Note that the code is written to run on the UNSW Katana cluster computer, which uses a scheduler to sequentially call the R script and pass it the particular iterations of the data to be processed in each step. To run the code on a standard computer, the code can be edited so the parameters passed by the Katana scheduler are defined internally.

| Description | Code |
| --- | --- |
| S1 - Data creation of Dataset 1 - Stata Code | [Data creation code](Code/S1_data_creation_dataset1.do) |
| S2 - Data creation of Dataset 2 - Stata Code | [Data creation code](Code/S2_data_creation_dataset2.do) |
| S3 - Data creation of Dataset 3 - Stata Code | [Data creation code](Code/S3_data_creation_dataset3.do) |
| S4 - Analysis - R Code | [Analysis code](Code/S4_analysis_code.R) |



