# Broughton physical classification for kelp version 2

__Created:__      2024/05/29 Cloned from Broughton version 1 to extend to LSSM work.
__Main author:__  Edward Gregr  
__Affiliation:__  SciTech Environmental Consulting   
__Location:__     Vancouver, BC   
__Contact:__      e-mail: ed@scitechconsulting.com | tel: 604-612-8324
__Last Update:__  2024/05/29   
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
Building on the k-means clustering of physical layers done for DFO/MSEA, this project focuses on building clusters for the Broughton region to support the creation of relevant oceanographic contexts. 

## Summary   

## Status
2024/05/29: Cloned from DFO work to allow clean delivery and further development. 

Next steps: Explore data from the Bianucci FVCOMM model and other data prepared by Romina Barossa. 

## Contents and methods


### Data loading
Source data include ... 

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

