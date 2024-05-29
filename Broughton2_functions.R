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

# "diffeR", "vegan", "ranger", "e1071", "forcats", "measures", "caret", "PresenceAbsence"
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
#raster_dir <- 'C:/Data/Git/Broughton/Data/OceanOnly'
raster_dir <- 'C:/Data/Git/Broughton/Data/Predictors'
data_dir   <- 'C:/Data/Git/Broughton/Data'
rmd_dir    <- 'C:/Data/Git/Broughton' 

# proj4 string for albers projection with NAD83 datum
spat_ref <- '+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'

# Define depth classes, in case needed
#z.breaks <- c( -5000, -50, -20, -10, -5, 0, 1000)
z_breaks <- c( -1000, 0, 5, 10, 20, 50, 100, 200, 5000)
z_ribs   <- c('ITD', '0-5', '5-10', '10-20', '20-50', '50-100', '100-200', '200+')


#===================================== Functions =========================================


#---- The full set of load, trim, and scale operations. ----
# Creates global raster stack for analysis. 
prepareData <- function() {
  
}
  

#---- MakeScreePlot: retuns a ggplot. ----
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


#---- Loads predictors produced by SPECTRAL lab from specified subdirectory -----
LoadPredictors <- function( pred_dir ) {
  
  print( list.files(path = pred_dir, pattern = '\\.tif$', full.names = FALSE) )

  # make a list of predictor rasters
  raster.list <- list.files(path = pred_dir, pattern = '\\.tif$', full.names = TRUE)
  
  
  # make list of raster names (without extension)
  raster.names <- lapply(raster.list, FUN = function(raster.layer){
    substr(basename(raster.layer), 1, nchar(basename(raster.layer)) - 4)
  } )
  
  # create a raster stack from raster_list
  # raster.stack <- raster::stack(x = raster.list)
  
  # Go thru one at a time so that we can set NA values ... needed for Romina's rasters
  raster.stack <- stack()
  for (i in raster.list) {
    new.layer <- raster( i )
    if (cellStats(new.layer, min) < -5) {
      values(new.layer)[ values(new.layer) < 0 ] <- NA }
    raster.stack <- addLayer( raster.stack, new.layer)
  }
  
  return(raster.stack)
}


ScalePredictors <- function( the_stack ){
# A little shimmy to avoid scaling substrate and name it consistently
#  scaled_stack <- scale( dropLayer( the_stack, "SUBSTRATE") )
#  scaled_stack <- stack( scaled_stack, trim_layers$SUBSTRATE )
# safe to assume its the last layer in the stack
#  names( scaled_stack[[ dim(scaled_stack)[[3]] ]] ) <- "substrate"
  
  x <- scale( the_stack )
  print('saving ...')
  today <- format(Sys.Date(), "%Y-%m-%d")
  save( x, file = paste0( data_dir, '/scaled_rasters_', today, '.rData' ))
  
  return (x)
} 


ClipPredictors <- function( data_layers, in_mask ){
# Restrict extent of predictors to specified mask
  
  x  <- raster::mask(data_layers, in_mask )
  today <- format(Sys.Date(), "%Y-%m-%d")
  save( x, file = paste0( data_dir, '/clipped_rasters_', today, '.rData' ))
  print('Clip rasters saved.')
  
  return(x)  
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



#-------------------------------------------------------------
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




#--------- VESTIGIAL BELOW -------------



# PlotMasksOnRaster <- function () {
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
# 
# # Obsolete for Romina's layers as they are all trimmed to 40 m depth, and land is NA. 
# # Resurect only if need to remove land and/or deep from bathy and derivates before scaling.
# 

# #---- Trim Raster - not quite working ... ----
# Trim.Raster <- function( a_raster, dType = 'NA' ) {
#   trimmed_data <- getValues( a_raster )
#   trimmed_data[ trimmed_data < -5 ] <- NA
#   if (dType == 'INT') {
#     a_raster <- setValues( a_raster, as.integer( trimmed_data )) }
#   else {
#     a_raster <- setValues( a_raster, trimmed_data ) }
#   return( a_raster )
# }
# BathyTrim <- function () {
#   # Removing all data above the HHWL, assumed to be 5 m.
#   # This is to avoid spurious outliers, and for visuallzation.
#   trim_layers <- data_layers
#   trim_data <- getValues( trim_layers$bathymetry )
#   trim_idx <- trim_data < -5
#   trim_data[ trim_idx ] <- NA
#   trim_layers$bathymetry <- setValues( trim_layers$bathymetry, as.integer( trim_data ))
#   trim_layers$bathymetry <- setMinMax( trim_layers$bathymetry )
#   
#   trim_data <- getValues( trim_layers$rugosity )
#   trim_data[ trim_idx ] <- NA
#   trim_layers$rugosity <- setValues( trim_layers$rugosity, trim_data )
#   trim_layers$rugosity <- setMinMax( trim_layers$rugosity )
#   
#   trim_data <- getValues( trim_layers$standard_deviation_slope )
#   trim_data[ trim_idx ] <- NA
#   trim_layers$standard_deviation_slope <- setValues( trim_layers$standard_deviation_slope, trim_data )
#   trim_layers$standard_deviation_slope <- setMinMax( trim_layers$standard_deviation_slope )
#   print( "Data trimmed ... ")
#   
#   today <- format(Sys.Date(), "%Y-%m-%d")
#   save( trim_layers, file = paste0( data_dir, '/trimmed_rasters_', today, '.rData' ))
#   
# } else {
#   print('Loading trimmed data ...')
#   load( paste0( data_dir, '/trimmed_rasters_2024-05-01.rData' ))
#   return( trimmed_layers )
# }  
# }





#-- 2024/02/12: Functions below not used ... 


# 

 
# fin.
