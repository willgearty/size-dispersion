This repository contains the required R code and data files to produce all of the analyses and plots in:  
‘Investigating the biotic and abiotic drivers of body size disparity in communities of non-volant terrestrial mammals’  
William Gearty, Lawrence H. Uricchio, and S. Kathleen Lyons

## Before running the code
Before running this R script, we suggest that you download this folder and set it as your R working directory. Required data files should then load when called in each script, and plot files will save within the same folder. Note that some data must be downloaded separately before the analyses can be performed ([see below](#data)).

## Running the code
Run `mammal_variance.R` to perform all analyses and produce all plots from the manuscript. Note that the analyses can take a few days to run completely.

## Required R packages
Unless otherwise stated, the following packages can all be installed from CRAN:

deeptime  
dggridR  
ggforce  
ggtern  
grid  
maptools (was retired in 2023; still available on R-Forge)  
moments  
MuMIn  
parallel  
patchwork  
quantreg  
relaimpo  
rgbif  
sf  
spData  
spDep  
stringr  
terra  
tidyverse  
viridis  

  
These scripts have been tested using R version 4.3.3 - 
Copyright (C) 2024 The R Foundation for Statistical Computing.

## Data
The IUCN range maps for terrestrial mammals are too large (~1GB total) to include in this repository. They should be downloaded from the IUCN here: https://www.iucnredlist.org/resources/spatial-data-download. Once downloaded, the files should be unzipped and placed into a into a folder named "MAMMALS_TERRESTRIAL_ONLY".

Taxonomy (taxonomy.csv) and synonym (synonyms.csv) files are included for IUCN terrestrial mammals (accessed 02/11/2024), but these may need to be updated from the IUCN website: https://www.iucnredlist.org/search?dl=true&permalink=6da36cdc-5f19-46e7-b980-aeb06b2b6208.

The EltonTraits database contains species-level foraging attributes for mammals (and birds) and is copied from figshare here: https://doi.org/10.6084/m9.figshare.c.3306933.v1. The publication describing this database is available here: https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1890/13-1917.1.

The Human Footprint Index (for 1993 and 2009) data are copied from the NASA Socioeconomic Data and Applications Center: https://sedac.ciesin.columbia.edu/data/set/wildareas-v3-2009-human-footprint.

The Millenium Assessment's deforestation hotspot data are also copied from the NASA Socioeconomic Data and Applications Center: https://sedac.ciesin.columbia.edu/data/set/ma-rapid-land-cover-change.

The WorldClim (v. 2.1) bioclimatic variable data are too large to include in this repository (~10GB total). They can be downloaded (30 second resolution) from the WorldClim website here: https://www.worldclim.org/data/worldclim21.html. Once downloaded, the files should be unzipped and placed into a folder named "WorldClim".

The elevation data are coped from the USGS Global multi-resolution terrain elevation data 2010 (GMTED2010): https://www.usgs.gov/coastal-changes-and-impacts/gmted2010. The data (mean statistic, 30 second resolution) can be downloaded here: https://topotools.cr.usgs.gov/gmted_viewer/gmted2010_global_grids.php.
