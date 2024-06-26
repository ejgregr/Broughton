#----------------------------------------------------------------------------
# Script:  broughton_functions.R
# Created: February 2024. EJG
#
# Purpose: Support building and evaluation of Seaweed Clusters for Coastal BC.
# Illustrated at 3 spatial extents using the Broughton Region.
#
# Notes:
#  - 2024/02/12: Created from Substrate code streamlined for the paper submission.


#================================== Load require packages =================================

# check for any required packages that aren't installed and install them
required.packages <- c( "ggplot2", "reshape2", "tidyr","dplyr", "raster", "stringr", "rasterVis",
                        "RColorBrewer", "factoextra", "ggpubr", "cluster", "rmarkdown","lubridate" )

# "diffeR", "vegan", "ranger", "rgdal", "e1071", "forcats", "measures", "caret", "PresenceAbsence"
# "randomForest", "spatialEco", "xlsx", "robustbase", "biomod2", "sp", "magrittr", "tinytex", "rmarkdown", "binr", 'gwxtab'

uninstalled.packages <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]

# install any packages that are required and not currently installed
if(length(uninstalled.packages)) install.packages(uninstalled.packages)

# require all necessary packages
lapply(required.packages, require, character.only = TRUE)
#lapply(required.packages, library, character.only = TRUE)


#=========================== Data sources and constants =====================================

#-- Set source and output directories. Directory will be created if doesn't exist; file will be overwritten if it does.
#raster.dir  <- 'C:/Data/SpaceData/Substrate2019/Predictors/QCS'
raster_dir <- 'C:/Data/Git/Broughton_descriptors/MSEA'
#raster_dir <- 'C:/Data/Git/Broughton/Data/Predictors'
data_dir   <- 'C:/Data/Git/Broughton/Data'
rmd_dir    <- 'C:/Data/Git/Broughton' 

# proj4 string for albers projection with NAD83 datum
spat_ref <- '+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'


#===================================== Functions =========================================

#---- Loads predictors from specified subdirectory -----
LoadPredictors <- function( pred_dir ) {
  
  print( list.files(path = pred_dir, pattern = '\\.tif$', full.names = FALSE) )

  # make a list of predictor rasters
  raster.list <- list.files(path = pred_dir, pattern = '\\.tif$', full.names = TRUE)
  
  
  # make list of raster names (without extension)
  raster.names <- lapply(raster.list, FUN = function(raster.layer){
    substr(basename(raster.layer), 1, nchar(basename(raster.layer)) - 4)
  } )
  
  # create a raster stack from raster_list
  raster.stack <- raster::stack(x = raster.list)
  return(raster.stack)
}


#---- ClipPredictors: masks away deeper water and limits extents based on a polygon mask
ClipPredictors <- function( stack_in, the_mask){
  z <- stack()
  for (i in 1:dim( stack_in)[[3]] ) {
    x <- raster::crop(stack_in[[i]], the_mask )
    y <- raster::mask(x, the_mask)
    z <- stack( z, y )
    print( paste0('layer ', i, ' clipped.'))
  }
  return( z )
}


ScalePredictors <- function( the_stack ){
# A little shimmy to avoid scaling substrate and name it consistently
#  scaled_stack <- scale( dropLayer( the_stack, "SUBSTRATE") )
#  scaled_stack <- stack( scaled_stack, trim_layers$SUBSTRATE )
# safe to assume its the last layer in the stack
#  names( scaled_stack[[ dim(scaled_stack)[[3]] ]] ) <- "substrate"
  
} 


#---- Returns a stack of integerized rasters from a raster stack ----
Integerize <- function( in_layers, sig = 1000 ) {
  int_layers <- brick()
  for (i in 1:nlayers(in_layers)) {
    raster_to_scale <- in_layers[[i]]
    int_raster <- as.integer( raster_to_scale * sig )
    int_layers <- stack( int_layers, int_raster )
  }
  
  names( int_layers ) <- names( in_layers )
  return( int_layers )
}

 
#---- Set unsuitable elevations to NA ----
DropNonHabitat <- function( data_in ){
  
  # To avoid spurious scaling and classification, and outliers, this 
  # removes all data above the HHWL (assumed to be 5 m) AND below 40 m depth. 
  # NOTE: Hard-coded for three MSEA layers.
  
  trim_layers <- data_in
  
  trim_data <- getValues( trim_layers$bathymetry )
  trim_idx <- (trim_data < -5 | trim_data > 40 )
  
  trim_data[ trim_idx ] <- NA
  trim_layers$bathymetry <- setValues( trim_layers$bathymetry, as.integer( trim_data ))
  trim_layers$bathymetry <- setMinMax( trim_layers$bathymetry )
  
  trim_data <- getValues( trim_layers$rugosity )
  trim_data[ trim_idx ] <- NA
  trim_layers$rugosity <- setValues( trim_layers$rugosity, trim_data )
  trim_layers$rugosity <- setMinMax( trim_layers$rugosity )
  
  trim_data <- getValues( trim_layers$standard_deviation_slope )
  trim_data[ trim_idx ] <- NA
  trim_layers$standard_deviation_slope <- setValues( trim_layers$standard_deviation_slope, trim_data )
  trim_layers$standard_deviation_slope <- setMinMax( trim_layers$standard_deviation_slope )
  print( "Nonhabitat removed.")

  return( trim_layers )
}


#---- MakeScreePlot: returns a ggplot. ----
# samp is optional, uses all dat if omitted.
MakeScreePlot <- function( indat, nclust, nrand, maxi, sampsize = 0 ){
  #initialize list for results
  wss <- numeric(nclust) 
  
  #subsample as requested
  if (sampsize > 0) {
    samp <- sample( 1:length( indat[ , 1] ), sampsize )
    dat <- indat[ samp, ]
  } else dat <- indat
  
  for (i in 1:nclust) {
    # Fit the model: km.out
    print( paste0( "centers ",i))
    km.out <- kmeans(dat, centers = i, nstart = nrand, iter.max = maxi)
    # Save the within cluster sum of squares
    wss[i] <- km.out$tot.withinss
    
    # calculate the silhouete width
  }
  
  # Produce the scree plot ... using a tibble, I guess. 
  wss_df <- tibble(clusters = 1:nclust, wss = wss)
  scree_plot <- ggplot(wss_df, aes(x = clusters, y = wss, group = 1)) +
    geom_point(size = 4)+
    geom_line() +
    scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
    xlab('Number of clusters') +
    ylab('Total within-cluster sum of squares') +
    geom_hline(
      yintercept = wss, 
      linetype = 'dashed')
  
  return( scree_plot)
}


#---- ClusterPCA: Returns a pair of PCA plots showing the separation of the clusters.
# NOTE: Uses a relatively small subset of the overall data so these will change with a different sample.
ClusterPCA <- function( n_samp, clustnum ) {
  
  ssidx <- sample( 1:length( stack_data_clean[ , 1] ), n_samp )
  ssamp <- stack_data_clean[ ssidx, ]
  
  # re-run cluster for smaller sample.
  cluster_result <- kmeans(ssamp, centers = clustnum, nstart = randomz) # less than 10 seconds
  csamp <- cluster_result$cluster
  
  # Create the data structure for PCA profiling (combining predictor data with the clusters). 
  # Put cluster # first so easy to find for PCA
  p_data <- as.data.frame( cbind(cluster = csamp, ssamp ) )
  
  res_pca <- prcomp(p_data[,-1],  scale = TRUE)
  # PC coordinates of individual raster cells
  ind_coord <- as.data.frame(get_pca_ind(res_pca)$coord)
  # Add clusters from the classification
  ind_coord$cluster <- factor( p_data$cluster )
  # Data inspection
  #head(ind_coord)
  
  # Percentage of variance explained by dimensions
  eigenvalue <- round(get_eigenvalue(res_pca), 1)
  var_percent <- eigenvalue$variance.percent
  
  # Look at the clusters for the first 4 PCs
  a <- ggscatter(
    ind_coord, x = "Dim.1", y = "Dim.2", 
    color = "cluster", palette = "simpsons", ellipse = TRUE, ellipse.type = "convex",
    size = 1.5,  legend = "right", ggtheme = theme_bw(),
    xlab = paste0("Dim 1 (", var_percent[1], "% )" ),
    ylab = paste0("Dim 2 (", var_percent[2], "% )" )
  ) +
    stat_mean(aes(color = cluster), size = 4)
  
  b <- ggscatter(
    ind_coord, x = "Dim.3", y = "Dim.4", 
    color = "cluster", palette = "simpsons", ellipse = TRUE, ellipse.type = "convex",
    size = 1.5,  legend = "right", ggtheme = theme_bw(),
    xlab = paste0("Dim 3 (", var_percent[3], "% )" ),
    ylab = paste0("Dim 4 (", var_percent[4], "% )" )
  ) +
    stat_mean(aes(color = cluster), size = 4)
  
  return( list(res_pca, a, b))
}


#---- PredictClusters: returns cluster assignments for un-clustered pictures. ----
# Assembly required with the classified pixels before plotting.
# Function by ChatGPT.
PredictClusters <- function(newdata, kmeans_model) {
  centers <- kmeans_model$centers
  dists <- as.matrix(dist(rbind(centers, newdata)))
  dists <- dists[-(1:nrow(centers)), 1:nrow(centers)]
  cluster_assignments <- apply(dists, 1, which.min)
  return(cluster_assignments)
}



 
####------------Depreciated or replaced functions----------------------------

#---- TrimStack: limit extents based on a polygon mask (vestigial)
# NOTE: This function now superseded by the crop/mask approach in ClipPredictors.
TrimStack <- function( stack_in, padsize ) {
  # Assigning extents here worked but results odd - plots appeared but looked empty.
  y <- stack()
  # estimate extents from first raster ... 
  x <- trim(stack_in[[1]], padsize )
  x_ext <- extent(x)
  x_ext <- ceiling( x_ext )
  extent( x ) <- x_ext
  y <- stack(y, x)
  print( paste0('layer 1 trimmed.'))
  
  # use these extents for the remaining layers ...
  for (i in 2:dim( stack_in)[[3]] ) {
    x <- setExtent( stack_in[[i]], x_ext, keepres=T, snap=T )
    y <- stack( y, x )
    print( paste0('layer ', i, ' trimmed.'))
  }
  return( y )  
}


# if (spacesub) {
#   
#   sMask <- shapefile("C:\\Data\\SpaceData\\Broughton\\broughton_small.shp")
#   smMask <- shapefile("C:\\Data\\SpaceData\\Broughton\\broughton_smaller.shp")
#   
#   # Data extents No mask and polygons
#   x  <- prepped_layers
#   
#   # Use only one of the loaded shapefiles
#   x  <- raster::mask(x, sMask )
#   # x  <- raster::mask(x, smMask )
#   
#   # Plot first layer in the resulting stack for inspection
#   plot(trim(x[[1]]), main = 'Local')
#   
#   # Plot masks over raster. NB: will only see if raster extends beyond the masks. 
#   sp::plot( sMask, border = "darkgreen", lwd = 2, add=TRUE)
#   sp::plot( smMask, border = "red", lwd = 2, add=TRUE)
#   
#   prepped_layers <- x
#   rm('x')
# }


# Apply scaling to standardize rasters. Does either min/max or z-score.
# Returns a scaled raster.
Scale.Raster <- function( in.raster, scale.how = "mm", intshift = 'T' ) {
  
  if (scale.how == "z") {
    mean.value <- cellStats( in.raster, stat='mean', na.rm=TRUE )
    sd.value   <- cellStats( in.raster, stat='sd', na.rm=TRUE )  
  } else {
    old.min    <- cellStats( in.raster, stat='min', na.rm=TRUE )  
    old.max    <- cellStats( in.raster, stat='max', na.rm=TRUE )  
  }
  
  if (scale.how == "z") {
    # Perform z-score normalization
    scaled.raster <- (in.raster - mean.value) / sd.value
  } else {
    # Perform min/max scaling
    scaled.raster <- (in.raster - old.min) / (old.max - old.min)
  }
  
  if (intshift) {
    scaled.raster <- as.integer( scaled.raster * 1000 )
  }
  
  return( scaled.raster )
}


# fin
