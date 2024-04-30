# Broughton biophysical classification for kelp at three spatial extents

__Created:__      2024/02/12 Forgive me father, its been over 2 years since my last R project.
__Main author:__  Edward Gregr  
__Affiliation:__  SciTech Environmental Consulting   
__Location:__     Vancouver, BC   
__Contact:__      e-mail: ed@scitechconsulting.com | tel: 604-612-8324
__Last Update:__  2024/02/12   
__Version:__      R version 4.3.2 x64

- [Objective](#objective)
- [Summary](#summary)
- [Status](#status)
- [Contents and methods](#contents and methods)
  + [Subsections within contents](#subsections-within-contents)
- [Data management](#data management)
- [Requirements](#requirements)
- [Caveats](#caveats)
- [Uncertainty](#uncertainty)
- [Acknowledgements](#acknowledgements)
- [References](#references)

## Objective
Following the methods of Mora-Soto et al. 2024, create a k-means clustering of biophysical layers to support the representation of kelp habitat in British Columbia, with a particular focus on the Broughton region.

From the DFO statement of work: This project will provide an R script to generate clusters representing areas of common environmental conditions from existing spatial layers, with examples of outputs generated across two spatial extents (i.e., Broughton Archipelago and the Queen Charlotte Strait region).

## Summary   

## Status
2024/03/22: Developmental milestone with data loaded, scaled, and displayed prior to clustering. Good learnings around rasters.

Next steps: cross-correlation table; data transformations; update with latest data; rubric for finalizing predictors, cluster number, cluster stability methods. 



## Contents and methods


### Data loading
Source data include rasters of regional 20 m predictors, and independent observations of kelp communities from the DFO field programs.

Data structures resulting from the data load include: ... These processed observations are used for the analysis and are included in the distributed code package.  

### Data analysis
Analytical steps include: ...
Data structures ... 
Using RMD to produce outputs ... 

## Data management  

## Requirements
Working directories are required to re-build the data. At a minimum, a source and output directory are needed. These are located near the top of the substrate_function.R script.

## Caveats
Versioning of tidyverse packages has been an issue during development, as functions continue to be depreciated. 

## Acknowledgements

## References
Mora-Soto et al. 2024.

