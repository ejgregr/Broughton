raster.list <- list.files(path = raster_dir, pattern = '\\.tif$', full.names = TRUE)

x <- raster( raster.list[ 1] )
str(x)
dim(x)

if (cellStats(x, min) < -5) {
  values(x)[ values(x) < 0 ] <- NA }

raster::plot(x, maxpixels = 2000000 )

writeRaster( x, paste0( data_dir, "/foo.tif"), overwrite=TRUE)





# make a list of predictor rasters
raster.list <- list.files(path = raster_dir, pattern = '\\.tif$', full.names = TRUE)
  
  # make list of raster names (without extension)
  raster.names <- lapply(raster.list, FUN = function(raster.layer){
    substr(basename(raster.layer), 1, nchar(basename(raster.layer)) - 4)
  } )
  
  # create a raster stack from raster_list

raster.stack <- raster::stack(x = raster.list)

names(raster.stack)
foo <- getValues(raster.stack$bathymetry)
min( na.omit(foo) )
max( na.omit(foo) )


names(scaled_layers)

str(scaled_layers$bathymetry)

x <- getValues( scaled_layers[[1]] )


names(data_layers)






# 
# in.raster <- data.layers[[3]]
# foo <- getValues( in.raster )
# 
# old.min <- min( na.omit(foo) )
# old.max <- max( na.omit(foo) )
# scaled.foo <- (foo - old.min) / (old.max - old.min)
# #scaled.foo <- as.integer( scaled.foo * 1000 )
# 
# 
# bar <- as.integer( in.raster )
# bar <- setValues( bar, scaled.foo )
# 
# 
# cellStats(bar, min, na.rm = TRUE)
# cellStats(bar, max, na.rm = TRUE)
# 
# plot(bar)
# raster::hist(bar, nclass=50, main = name(bar))
# 




if (scale.how == "z") {
  # Perform z-score normalization
  mean.value <- mean( na.omit( foo ))
  sd.value   <- sd( na.omit( foo ))
  scaled.foo <- (foo - mean.value) / sd.value }
else {
  # Perform min/max scaling
}

if (intshift) {
  scaled.foo <- as.integer( scaled.foo * 1000 )
}