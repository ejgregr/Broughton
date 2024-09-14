---
title: "Cluster analysis of nearshore physical predictors"
author: "Edward Gregr"
date: "September 12, 2024"
output:
  html_document: 
    keep_md: yes
    theme: yeti
    highlight: haddock
    toc: yes
    toc_float: yes
---
<style type="text/css">
#TOC {
  color: black;
}
</style>


This work is Part One of a classification exercise intended to partition nearshore coastal waters into ecologigally meaningful clusters. The objective here is available kelp habitat - floating plants attached to hard substrate. Part Two of the work (the Broughton LSSM) explores more general(?) classifications.

This script describes the classification framework delivered to DFO, defined using a substantial body of existing data, to support the organisation and analysis of relevant macroalgal data for the Pacific Regions. 

# Objective
This project provides an R script to generate clusters representing areas of common environmental conditions from existing spatial layers. Examples of outputs are generated across three spatial extents: The full QCS 'data' region, a smaller, region with higher data quality, and a small, self-similar area within the Broughton Archipelago.

## Background   
Following the methods of Mora-Soto et al. 2024, create a k-means clustering of biophysical layers to support the representation of kelp habitat in British Columbia, with a particular focus on the Broughton region.

K-means is described as sensitive to the 'shape' of the data (e.g., range, distribution). Much attention was therefore paid to outliers and skewness. All predictors were centered, then transformed if necessary, and finally scaled so all predictors contribute equally.

## Data summary

| Process         | Predictor          | Description and rationale  |
|-----------------|-----------------|-----------------|
| Light, Energy   | bathymetry      | Kelps are light restricted and thus depth is essential for creating clusters limited to the photic zone. Following Barbosa et al., we used the bathymetry to mask out unsuitable depths rather than use it as a classifying variable.                                   |
| Circulation     | circ_mean_summer     | intended proxy for water mass movements     |
|                 | tidal_mean_summer   | A representation of diurnal and fortnight tidal ranges |
| Freshwater input | fresh_index        | Kelps do better in higher salinity waters, as salinity tends to be correlated with both cooler temperatures and increased nutrients. We used a simple diffusion model developed by MSEA which shows potentially low salinity areas based on point estimates of riverine input (from the BC Freshwater Atlas).|
|Relative exposure index | REI | In addition to influencing substrate, exposure is also an indicator of mixing and can also influence zoospore settlement. This updated exposure layer from MSEA includes depth attenuated, relative exposure index combines fetch with wind and depth. |
|Salinity	| salinity_mean_summer    |         |
|         |  salinity_range           |       |
|Temperature 	| sentinel_max          |       |
|             | sentinel_mean         |       |
|             | sentinel_sd           |       |
|             | temp_mean_summer      |       |
|             |temp_range             |       |




# Source data

Dimensions of the loaded Raster stack.

```
## [1] 4937 9684    6
```

<div class="figure" style="text-align: center">
<img src="C:/Data/Git/Broughton/Broughton_DFO_files/figure-html/unnamed-chunk-3-1.png" alt="__*Figure 1: Maps of the selected, unmodified predictors.*__" width="100%" />
<p class="caption">__*Figure 1: Maps of the selected, unmodified predictors.*__</p>
</div>

<div class="figure" style="text-align: center">
<img src="C:/Data/Git/Broughton/Broughton_DFO_files/figure-html/unnamed-chunk-4-1.png" alt="__*Figure 2: Histograms of the selected, unmodified predictors.*__" width="100%" />
<p class="caption">__*Figure 2: Histograms of the selected, unmodified predictors.*__</p>
</div>

Examine correlations across predictors, and identify those beyond a threshold of 0.6.

Table: Correlation matrix for assessing predictor cross-correlations

|                         | bathymetry| circ_mean_summer| qcs_freshwater_index|  SUBSTRATE|    rei_qcs| salinity_mean_summer| salinity_range| standard_deviation_slope| temp_mean_summer| temp_range| tidal_mean_summer| sst_sentinel_20m_bi_max| sst_sentinel_20m_bi_mean| sst_sentinel_20m_bi_sd|
|:------------------------|----------:|----------------:|--------------------:|----------:|----------:|--------------------:|--------------:|------------------------:|----------------:|----------:|-----------------:|-----------------------:|------------------------:|----------------------:|
|bathymetry               |         NA|       -0.0596522|           -0.0792392|  0.1228241| -0.1211747|            0.1970434|      0.0063756|               -0.1986177|       -0.1713716| -0.1092642|        -0.0764030|              -0.1283650|               -0.1285915|             -0.1936621|
|circ_mean_summer         |         NA|               NA|           -0.0313665| -0.2496822| -0.2514569|           -0.1021889|     -0.1130768|                0.0589342|        0.0203808| -0.1049203|         0.6832543|              -0.0017839|               -0.1010601|              0.0842129|
|qcs_freshwater_index     |         NA|               NA|                   NA| -0.0149079| -0.0347575|           -0.0534843|     -0.0301344|                0.0630604|        0.0351385|  0.0145500|        -0.0254697|               0.0046085|                0.0193847|              0.0433353|
|SUBSTRATE                |         NA|               NA|                   NA|         NA|  0.3002426|            0.1367176|      0.1288244|               -0.3615454|        0.0039791|  0.0809443|        -0.2065128|              -0.0360729|               -0.0036349|             -0.1583280|
|rei_qcs                  |         NA|               NA|                   NA|         NA|         NA|            0.2796995|      0.1798324|               -0.2156858|       -0.1569911| -0.0015175|        -0.2527533|              -0.1920590|               -0.0930388|             -0.2813880|
|salinity_mean_summer     |         NA|               NA|                   NA|         NA|         NA|                   NA|     -0.4908089|               -0.1970693|       -0.9481338| -0.8141805|        -0.3309511|              -0.5438599|               -0.5547659|             -0.5534339|
|salinity_range           |         NA|               NA|                   NA|         NA|         NA|                   NA|             NA|               -0.0273952|        0.5574500|  0.7620925|         0.0045803|               0.1685611|                0.2325462|              0.1004962|
|standard_deviation_slope |         NA|               NA|                   NA|         NA|         NA|                   NA|             NA|                       NA|        0.0893897|  0.0278378|         0.0265606|               0.1404021|                0.1714052|              0.2131796|
|temp_mean_summer         |         NA|               NA|                   NA|         NA|         NA|                   NA|             NA|                       NA|               NA|  0.9089467|         0.3212514|               0.5181216|                0.5111020|              0.5274607|
|temp_range               |         NA|               NA|                   NA|         NA|         NA|                   NA|             NA|                       NA|               NA|         NA|         0.2048309|               0.4232786|                0.4586944|              0.4073992|
|tidal_mean_summer        |         NA|               NA|                   NA|         NA|         NA|                   NA|             NA|                       NA|               NA|         NA|                NA|               0.1722931|               -0.0281543|              0.2882165|
|sst_sentinel_20m_bi_max  |         NA|               NA|                   NA|         NA|         NA|                   NA|             NA|                       NA|               NA|         NA|                NA|                      NA|                0.7799091|              0.7998959|
|sst_sentinel_20m_bi_mean |         NA|               NA|                   NA|         NA|         NA|                   NA|             NA|                       NA|               NA|         NA|                NA|                      NA|                       NA|              0.4994191|
|sst_sentinel_20m_bi_sd   |         NA|               NA|                   NA|         NA|         NA|                   NA|             NA|                       NA|               NA|         NA|                NA|                      NA|                       NA|                     NA|



Table: Predictor variables that exceed 0.6 threshold

|                        | bathymetry| circ_mean_summer| qcs_freshwater_index|  SUBSTRATE|    rei_qcs| salinity_mean_summer| salinity_range| standard_deviation_slope| temp_mean_summer| temp_range| tidal_mean_summer| sst_sentinel_20m_bi_max| sst_sentinel_20m_bi_mean| sst_sentinel_20m_bi_sd|
|:-----------------------|----------:|----------------:|--------------------:|----------:|----------:|--------------------:|--------------:|------------------------:|----------------:|----------:|-----------------:|-----------------------:|------------------------:|----------------------:|
|circ_mean_summer        |         NA|               NA|           -0.0313665| -0.2496822| -0.2514569|           -0.1021889|     -0.1130768|                0.0589342|        0.0203808| -0.1049203|         0.6832543|              -0.0017839|               -0.1010601|              0.0842129|
|salinity_range          |         NA|               NA|                   NA|         NA|         NA|                   NA|             NA|               -0.0273952|        0.5574500|  0.7620925|         0.0045803|               0.1685611|                0.2325462|              0.1004962|
|temp_mean_summer        |         NA|               NA|                   NA|         NA|         NA|                   NA|             NA|                       NA|               NA|  0.9089467|         0.3212514|               0.5181216|                0.5111020|              0.5274607|
|sst_sentinel_20m_bi_max |         NA|               NA|                   NA|         NA|         NA|                   NA|             NA|                       NA|               NA|         NA|                NA|                      NA|                0.7799091|              0.7998959|

__*Table 1: Transforms, how they were done and the results.*__

```
## [1] "Hi World!"
```



<div class="figure" style="text-align: center">
<img src="C:/Data/Git/Broughton/Broughton_DFO_files/figure-html/FixedHists-1.png" alt="__*Figure 3: Histograms of selected, transformed and scaled predictor data.*__" width="100%" />
<p class="caption">__*Figure 3: Histograms of selected, transformed and scaled predictor data.*__</p>
</div>

# Results
## Part 1 - Cluster number
<div class="figure" style="text-align: center">
<img src="C:/Data/Git/Broughton/Broughton_DFO_files/figure-html/ScreePlot-1.png" alt="__*Figure 4: Scree plot showing the total within sum-of-squares across a range of cluster numbers.*__" width="100%" />
<p class="caption">__*Figure 4: Scree plot showing the total within sum-of-squares across a range of cluster numbers.*__</p>
</div>

<div class="figure" style="text-align: center">
<img src="C:/Data/Git/Broughton/Broughton_DFO_files/figure-html/HeatMap-1.png" alt="__*Figure 5: Heat map showing within cluster standard deviation of the predictors.*__" width="100%" />
<p class="caption">__*Figure 5: Heat map showing within cluster standard deviation of the predictors.*__</p>
</div>

<div class="figure" style="text-align: center">
<img src="C:/Data/Git/Broughton/Broughton_DFO_files/figure-html/SilhouettePlot-1.png" alt="__*Figure 6: Silhouette plot showing pixel membership in each cluster and silhouette widths.*__" width="100%" />
<p class="caption">__*Figure 6: Silhouette plot showing pixel membership in each cluster and silhouette widths.*__</p>
</div>


## Part 2: Clusters and predictor loadings 


```
## [[1]]
## Standard deviations (1, .., p=6):
## [1] 1.5781901 1.1682439 0.9413997 0.7920353 0.7157889 0.3444054
## 
## Rotation (n x k) = (6 x 6):
##                                 PC1        PC2         PC3        PC4        PC5          PC6
## rei_qcs                  -0.3180046 -0.5807007  0.10703319  0.3963809 -0.5802245  0.237540329
## salinity_mean_summer     -0.5824987  0.1478529 -0.08785136 -0.1770014 -0.2639274 -0.728100165
## standard_deviation_slope  0.1964921  0.5933617  0.47347180  0.5321103 -0.2978589 -0.115221023
## temp_range                0.4640925 -0.4899460  0.03339078  0.3248114  0.2066680 -0.628682758
## tidal_mean_summer         0.3226027  0.1914152 -0.81394538  0.1576624 -0.4145369 -0.009074291
## sst_sentinel_20m_bi_max   0.4488929 -0.1102783  0.30500567 -0.6309132 -0.5388966 -0.069602472
## 
## [[2]]
```

<div class="figure" style="text-align: center">
<img src="C:/Data/Git/Broughton/Broughton_DFO_files/figure-html/PCAPlots-1.png" alt="__*Figure 7: PCA Plots for a) dimensions 1 and 2, and b) dimensions 3 and 4.*__" width="100%" />
<p class="caption">__*Figure 7: PCA Plots for a) dimensions 1 and 2, and b) dimensions 3 and 4.*__</p>
</div>

```
## 
## [[3]]
```

<div class="figure" style="text-align: center">
<img src="C:/Data/Git/Broughton/Broughton_DFO_files/figure-html/PCAPlots-2.png" alt="__*Figure 7: PCA Plots for a) dimensions 1 and 2, and b) dimensions 3 and 4.*__" width="100%" />
<p class="caption">__*Figure 7: PCA Plots for a) dimensions 1 and 2, and b) dimensions 3 and 4.*__</p>
</div>

<div class="figure" style="text-align: center">
<img src="C:/Data/Git/Broughton/Broughton_DFO_files/figure-html/ViolinPlot-1.png" alt="__*Figure 8: Violing plots showing distribiton of predictors in each of the k-means clusters.*__" width="100%" />
<p class="caption">__*Figure 8: Violing plots showing distribiton of predictors in each of the k-means clusters.*__</p>
</div>


__*Table XX: Need some tablulated PCA loadings ... .*__




# Addendums

## R packages

dplyr: A fast, consistent set of tools for working with data frame like objects, both in memory and out of memory.  

ggplot2: A system for 'declaratively' creating graphics based on ``The Grammar of Graphics''. You provide the data, tell 'ggplot2' how to map variables to aesthetics, what graphical primitives to use, and it takes care of the details. 
Hmisc: Contains many functions useful for data analysis, high-level graphics, utility operations, functions for computing sample size and power, simulation, importing and annotating datasets,
imputing missing values, advanced table making, variable clustering, character string manipulation, conversion of R objects to LaTeX and html code, and recoding variables. Used here to add standard deviation to as lines to ggplot.

knitr: Provides a general-purpose tool for dynamic report generation in R using Literate Programming techniques. Used here for Markdown formatting as html or pdf.

lubridate : Functions to work with date-times and time-spans: fast and user friendly parsing of date-time data, extraction and updating of components of a date-time (years, months, days, hours, minutes, and seconds), algebraic manipulation on date-time and time-span objects. Used here to pool dates to month.

markdown: R Markdown allows the use of knitr and Pandoc in the R environment. This package translates R Markdown to standard Markdown that knitr and Pandoc render to the desired output format (e.g., PDF, HTML, Word, etc).

purrr: A complete and consistent functional programming toolkit for data manupulation in R.  

readxl: Imports excel files into R.  

reshape2: Flexibly restructure and aggregate data using just two functions: melt and 'dcast'. Used here for ggplot() support. 

tibble: Functions for manipulating the tibble data format in R. Used here for the deframe() function.  

vegan: Ordination methods, diversity analysis and other functions for community and vegetation ecologists.  

