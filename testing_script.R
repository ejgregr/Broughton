load( paste0( data_dir, '/tif_stack_2024-09-05.rData' ))
prepped_layers <- dropLayer( tif_stack, c("bathymetry", "SUBSTRATE", "salinity_range", "temp_mean_summer", "circ_mean_summer",
                                               "sst_sentinel_20m_bi_mean", "sst_sentinel_20m_bi_sd") )

stack_data <- getValues(prepped_layers)

#---- Outlier Detection ----
skews <- vector()
for (i in 1:dim(stack_data)[2]) {
  j <- skewness( stack_data[,i], na.rm = T )
  skews <- c( skews, j)
}
names(skews) <-  colnames(stack_data)

skews <- as.data.frame(skews)
skews <- cbind( skews, "lambda" = 0)
skews



j <- stack_data[ !is.na(stack_data[,1]),1]
skewness(j)

if (abs( skewness( j )) > 1) {
  j <- j+abs(min(j))
  bc_result <- boxcox( lm(j+1~1) )
  lambda <- bc_result$x[ which.max( bc_result$y )]
  skews[i, "lambda"] <- lambda
  t_data <- (j^lambda - 1) / lambda
} else {
  t_data <-i }

skewness(t_data[1:2000], na.rm=T) #doesnt work with above t_data. only for some subsets



j <- j+abs(min(j))
t_data <- car::powerTransform( j+1 )

skewness(t_data, na.rm=T)

sum(j==0)/length(j)
j_trans <- j^(1/5) #not asn()
histogram(j_trans)
skewness(j_trans, na.rm=T)


# Drop all the zeros?
range(j)
j <- stack_data[ !is.na(stack_data[,1]),1]
j <- j[ !(j==0)]



foo <- fixDistributions( stack_data_clean, skews)

#Input is a matrix with columns as predictors, and the skew vector
fixDistributions <- function( a_stack, skews ){
  t_stack <- stack
  
  for (i in 1:dim(a_stack)[2]) {
    print(i)
    j <- a_stack[,i] # ith column
    if (abs( skewness( j )) > 1) {
        
      
      skews[i, "lambda"] <- lambda
      t_data <- (j^lambda - 1) / lambda
      t_stack <- stack( t_stack, t_data)
    } else {
      t_stack <- stack(t_stack, i) }
  }
  return(t_stack)
}















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