
names(data.layers)
names(scaled.layers)







names( scaled.layers ) <- names(data.layers)



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