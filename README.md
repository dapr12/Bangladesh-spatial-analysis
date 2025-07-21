# Bayesian Spatial Analysis of Unemployment in Bangladesh (2013)

## Overview

This repository contains the R code and analysis for a geospatial study of unemployment patterns across the districts of Bangladesh. Using the **2013 Bangladesh Labour Force and Child Labour Survey (LFS)** microdata, this project implements a **Bayesian hierarchical model (BYM2)** to:

1.  Model the district-level unemployment rate while accounting for spatial autocorrelation.
2.  Identify the effect of urbanization on unemployment.
3.  Visualize and identify geographic hotspots and coldspots of unemployment risk that are not explained by covariates alone.

The primary tool used for this analysis is the `R-INLA` package, valued for its efficiency in fitting complex Bayesian models.

## Key Features

- **Data Aggregation:** Scripts for processing and aggregating the raw LFS 2013 survey data to the district level, applying survey weights for representative estimates.
- **Spatial Modeling:** Implementation of a Besag-York-Mollié (BYM2) model to separate structured (spatially correlated) and unstructured random effects.
- **Visualization:** Generation of high-quality choropleth maps using `ggplot2` and `sf` to display raw rates, model-smoothed rates, and latent spatial risk patterns.
- **Reproducibility:** The code is structured to be run from start to finish, ensuring the analysis is fully reproducible.

## Methodology

The model analyzes the count of unemployed persons in each district as a binomial outcome, with the total labor force as the number of trials. The key predictor variable is the proportion of the district's population that is urban. The BYM2 framework allows us to robustly estimate these effects by "borrowing strength" from neighboring districts, providing a clearer picture of the true underlying spatial patterns.

## How to Use

1.  Clone the repository: `git clone https://github.com/your-username/bangladesh-unemployment-spatial-analysis.git`
2.  Open the R project or script file in RStudio.
3.  Ensure you have the required packages installed (the script includes an installer).
4.  Place the `LFS-2013-By Quarter.dta` file in the appropriate directory (e.g., a `/data` subfolder).
5.  Update the file path in the main R script and run it.

## Data Source

The microdata used in this project is from the **Bangladesh - Labour Force and Child Labour Survey 2013**, conducted by the Bangladesh Bureau of Statistics (BBS) and made available via the ILO Microdata Repository. Due to data use agreements, the raw `.dta` file is not included in this repository.
