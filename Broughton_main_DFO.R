#################################################################################
# Script:  Broughton_main_DFO.R - DFO CLASSIFICATION VERSION
# Created: February 2024. EJG
# 
# This script sources the necessary libraries and functions, coordinates the 
# analysis, and creates data structures that are then 'knitted' together (I think).
# So the idea is to run bits in here, and then 'Render' the RMD script. 
# Seems straightforward. :)
#
# Updates: 
# 2024/04/29: Steady development past few weeks; alll necessary pieces now developed. 
# 2024/04/29: Git repo created and working version pushed prior to RMarkdown development.
# 2024/05/02: Completed smoothing pass thru code; sorted raster plotting. Pushed.
# 2024/05/07: Another pass thru, adding some controls. Ready for RMD work.Pushed.
# 2024/05/29: After some exploratory work with Romina's tifs from SPECTRAL lab, forked the 
#   code so this version focuses on DFO data/objectives and Broughton2 attends to the LSSM needs.
# 2024/05/29: Back after summer. Nothing this WILL NOT use Romina's data.
# 2024/08/30: Back on this. Renamed the 2 projects for clarity, and updated here with LSSM progress. 
# 2024/09/10: Working now. Have spent days looking at distributions, outliers, and skew. 
#   Substrate has now joined bathy as a necessary characteristic. Almost ready to to start running 
#   some RMD reports and comparing results.
# 2024/09/11: A few minor(ish) changes: consolidate all changes to data (ie, transforming, centering, 
#   scaling) in one place. Bathy and Substrate applied as exclusions. 
#################################################################################

print('Starting Broughton - DFO Version ...')
rm(list=ls(all=T))  # Erase environment.

# Load necessary packages and functions ... 
source( "broughton_functions.R" )
# source( "Plot_Functions.R" )

# Directories ...
#-- Source and output directories. Will be created if doesn't exist, overwritten if it does.
raster_dir <- 'C:/Data/SpaceData/Classification/MSEA'
data_dir   <- 'C:/Data/Git/Broughton/Data'
rmd_dir    <- 'C:/Data/Git/Broughton' 

# Processing FLAGS...
loadtifs <- T # If true the data will be re-loaded from TIFs, else it will be loaded from rData.
clipdata <- F # If true a spatial subset of the data will be taken based on a polygon shape file. 
scaledat <- T # If true, imported data will be scaled using raster::scale().
reclust  <- T # If true, re-cluster full data set prior to mapping, else predict to unclassified pixels.
addKmDat <- T

#---- Part 1 of 3: Load, clean, and prepare predictor data.  ----
# If loadtifs == TRUE then run all this else load the processed data.

tif_stack <- stack()
today <- format(Sys.Date(), "%Y-%m-%d")

if (loadtifs) {
  print( "Loading predictors ... ")
  src_stack <- LoadPredictors( raster_dir, addKmDat )
  print( "Data loaded.")

  tif_stack <- src_stack
  
  # The landside on the MSEA bathymetry is not of interest to the classification.
  # Deeper areas are also unsuitable for kelps. This removes those pixels from all rasters
  tif_stack <- DropNonHabitat( tif_stack, -5, 40 )
  print('Unsuitable elevations removed.')

  if (clipdata) {
    print( "clipping TIFs ... ")
    amask <- shapefile("C:\\Data\\SpaceData\\Broughton\\broughton_region.shp")
    #amask <- shapefile("C:\\Data\\SpaceData\\Broughton\\broughton_smaller.shp")
    tif_stack <- ClipPredictors( tif_stack, amask )
    print('Rasters clipped.')
  }
  
  if (scaledat) {
    print( "Centering  TIFs ... ")
    tmp_stack <- scale( tif_stack, scale=F )
    tif_stack <- tmp_stack
    print('Rasters centered.')
  }

  save( src_stack, file = paste0( data_dir, '/src_stack_', today, '.rData' ))
  save( tif_stack, file = paste0( data_dir, '/tif_stack_', today, '.rData' ))
  save( tif_stack, file = paste0( data_dir, '/tifs_DFO_scaled_QCS_', today, '.rData' ))
  
} else {
  print( 'Loading project data ... ')
  # Ideally meaningfully named and tested so no source is required.
#  load( paste0( data_dir, '/tifs_DFO_scaled_QCS_2024-09-05.rData' ))
#  load( paste0( data_dir, '/tifs_DFO_centred_QCS_2024-09-05.rData' ))
  load( paste0( data_dir, '/tif_stack_2024-09-05.rData' ))
}


#---- Final data preparation ----
#NOTE: Everything from here down is in scaled units. 

#-- Intergerize scaled data to reduce data size and improve performance. 
# NOTE: 2024/06/03: Integerize throws a warning (layers with no data) on the MSEA data set.
#   but no problem evident in results. Less useful with trimmed data, but may still have 
#   value for larger data sets.

#prepped_layers <- Integerize( tif_stack )
prepped_layers <- tif_stack

#-- Quick correlation across data layers
x <- getValues( prepped_layers )
x_clean <- x[ complete.cases(x), ]
cor_table <- cor( x_clean )
cor_table[lower.tri(cor_table)] <- NA
(cor_table >= 0.6) & (cor_table != 1)

#-- Remove correlated layers (see RMD document) 
prepped_layers <- dropLayer( prepped_layers, c("bathymetry", "SUBSTRATE", "salinity_range", "temp_mean_summer", "circ_mean_summer",
                                    "sst_sentinel_20m_bi_mean", "sst_sentinel_20m_bi_sd") )
names( prepped_layers )

# RENAME variables - after selection.

#-- Visualize source data
dim( prepped_layers )
names( prepped_layers )
plot( prepped_layers )
raster::hist(prepped_layers, nclass=50)
par(mfrow = c(1, 1))

#-- Finalize the data for exploring and creating the clusters.
  # Extract the data for the cluster analyses
stack_data <- getValues( prepped_layers )

# REMOVE fix some (hard-coded) distributions by adding ceilings and root transforms.
t_stack_data <- MakeMoreNormal( stack_data )

# remove any rows with an NA
# NB: This decouples the data from the RasterStack and requires re-assembly for mapping
# THESE are the two key data structures used in subsequent steps
clean_idx <- complete.cases(stack_data)
stack_data_clean <- stack_data[ clean_idx, ]

dim( stack_data )
dim( stack_data_clean )


#---- Part 2 of 3: Cluster number selection ----

set.seed <- 42 # Seed for reproducibility
randomz  <- 20 # the number of randomizations for kmeans to do.
imax     <- 25 # maximum iterations to try for convergence

#---- Part 2a: Explore number of clusters using Within-sum-of-squares scree plot ----
# Runs kmeans with increasing number of clusters

nclust   <- 18 # number of clusters for scree plot
nsample  <- 50000 # Scree plot needs a subsample to run reasonably. 
plotme <- MakeScreePlot( stack_data_clean, nclust, randomz, imax, nsample )
plotme

#---- Create a working set of N clusters (N based on scree plot) to further assess cluster number. ----

nclust  <- 7 # the number of clusters based on scree plot, above.
nsample <- 500000 # a larger sample for more robust classification

sidx <- sample( 1:length( stack_data_clean[ , 1] ), nsample )
samp <- stack_data_clean[ sidx, ]
cluster_result <- kmeans(samp, centers=nclust, nstart=randomz, iter.max=imax) 

#---- Part 2b: Create heat map of within-cluster standard deviations ----

# Define color palette
pal_heat <- rev( brewer.pal(n = nclust, name = "RdYlBu")) # heat map palette

profile_data <- as.data.frame( cbind(cluster = cluster_result$cluster, samp ) )

cluster_sd <- profile_data %>%
  group_by(cluster) %>%
  summarise_all(sd)

x <- as.data.frame( cluster_sd )
head(x)
xm <- melt( x, id.var = "cluster" )

z_heat <- ggplot(xm, aes(x=cluster, y=variable, fill=value) ) +
  geom_tile() +
  scale_fill_gradientn(colours = pal_heat) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(title = "Within-cluster Standard Deviation", x = "Clusters", y = "Attributes", fill = "Value")
z_heat

#---- Part 2c: Examine silhouette plot of the WORKING clusters  ----
# Uses the predictor values and the corresponding assigned cluster
# Need to subsample from the cluster result above as distance matrix take long time.

# Take a subsample of the clusters and the predictors for the silhouette plot. 
sil_n <- 10000
silx <- sample( 1:length( samp[ , 1] ), sil_n )

cs <- cluster_result$cluster[ silx ]
ss <- samp[ silx, ]

# Calculate a distance matrix between the predictors and plot assigned to each cluster.
# both these steps are time consuming, hence a smaller sample.
c_dist <- dist(ss)
sk <- silhouette(cs, c_dist)
#mean( sk[,"sil_width"] )

par(mfrow = c(1, 1))
plot(sk, col = 1:nclust, border=NA, main = "Hi World" )


#---- Part 3 of 3: Detailed examination of N clusters  ----
#---- Part 3a: Show cluster groupings using PCA ----

#-- Can take some time so it makes its own cluster 
pca_n <- 25000

pca_plots <- ClusterPCA( pca_n, nclust ) # uses global variable stack_data_clean
pca_plots

#Percentage of variance explained by dimensions
#eigenvalue <- round(get_eigenvalue(res_pca), 1)
#var_percent <- eigenvalue$variance.percent


#---- Part 3b: Violins of predictor contributions to WORKING clusters ----

x <- as.data.frame( samp )
x$cluster <- as.factor( cluster_result$cluster )

y <- x %>%
  pivot_longer(cols = -cluster, names_to = "predictor", values_to = "value")

# Create violin plot
ggplot(y, aes(x = cluster, y = value, fill = cluster)) +
  geom_violin(trim = FALSE) +
  facet_wrap(~ predictor, scales = "free_y") +
  theme_minimal() +
  labs(title = "Violin Plots of Predictors Across k-means Clusters",
       x = "Cluster",
       y = "Value")


#---- Part 3c: Spatialize the WORKING clusters ----
# NB: To show a comprehensive map, can either:
#     a) re-cluster the entire data set (using imax and randomz from above) or
#     b) Predict to the unsampled portion of the raster. 

# initialize target data structure 
cluster_raster <- prepped_layers[[1]]
dataType(cluster_raster) <- "INT1U"
cluster_raster[] <- NA

if (reclust == T) {
  # Re-cluster using all the clean data  ... 
  # less than 1 min with iter.max = 20, nstart = 20 for smallest region
  cluster_result <- kmeans(stack_data_clean, centers = nclust, nstart = randomz, iter.max = imax)
  
  # Assign the clustered values ... 
  # extract values from the target cluster
  new_values <- values( cluster_raster )
  # replace non-NA values with the cluster results
  new_values[ clean_idx ] <- cluster_result$cluster
  # put the updated values back on the target cluster
  values( cluster_raster ) <- new_values  
} else {
  # Predict values for unclustered cells. Can be more time-consuming than re-classifying everything. 
  values( cluster_raster ) <- transferCluster( values(cluster_raster), cluster_result )
}

#--- Display the results, first as histogram then as map.  
raster::hist( values(cluster_raster ))

# Define color palette
pal_clust <- brewer.pal(n = nclust, "Accent") # Max for Accent is 8

ckey <- list( at=0:nclust, 
              title="Clusters", 
              col = pal_clust)
myTheme <- rasterTheme( region = pal_clust )
z_map <- levelplot( cluster_raster, margin = F, 
           colorkey = ckey,
           par.settings = myTheme,
           main = "K-means Clustering Result - Local extents" )
z_map

writeRaster( cluster_raster, paste0( data_dir, "/SPEC_7clusterb.tif"), overwrite=TRUE)


#---- 


#---- Knit and render Markdown file -----
### Process file 
# To HTML ... 
rmarkdown::render( "Broughton_DFO.Rmd",   
                   output_format = 'html_document',
                   output_dir = rmd_dir )  

# 2024/04/29: It looks like this has gotten easier in the last 2 years ... version up!

# To PDF:
# First had to install the library tinytex.
# then run >tinytex::install_tinytex()
# ... and done. 
rmarkdown::render( "Broughton_DFO.Rmd",   
                   output_format = 'pdf_document',
                   output_dir = rmd_dir )  


#---- Some details on correlation analysis ... ----
#---- Correlation across UN-scaled datalayers ----
# foo <- getValues( scaled_layers )
# foo_clean <- na.omit(stack_data)
# pick <- sample( 1:length( foo_clean[ , 1] ), 10000 )
# plot( foo_clean[ pick,'rugosity'] ~ foo_clean[ pick,'standard_deviation_slope'] )

# RUGOSITY is a bit of a problem distribution
#cellStats( data_layers$rugosity, stat="range" )
#raster::hist(log( data_layers$rugosity+10 ), nclass=50)
# Look at bottom roughness relationships  
#pick <- sample( 1:length( stack_data_clean[ , 1] ), 10000 )
#plot( stack_data_clean[ pick,'rugosity'] ~ stack_data_clean[ pick,'standard_deviation_slope'] )

#---- 

#---- Plot any classified outliers in classified space --
# create a df w the centers for each cluster to facilitate distance calculation.
centers <- cluster_result$centers[cluster_result$cluster, ] 
head(centers)
# calculate distance
distances <- sqrt(rowSums((stack_data_clean - centers )^2))
head(x)
outlier_idx <- order(distances, decreasing=T)[1:1000000]
# a subsample to reduce the plot size ... 
ssidx <- sample( 1:length( stack_data_clean[ , 1] ), 10000 )

plot(stack_data_clean[ ssidx, c("rei_qcs", "qcs_freshwater_index") ], pch=19, col=cluster_result$cluster[ ssidx ], cex=1)
points(cluster_result$centers[ ,c("rei_qcs", "qcs_freshwater_index") ], col=1:3, pch=15, cex=2)
#----


# FIN.






