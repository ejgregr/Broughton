#################################################################################
# Script:  Broughton_main.R
# Created: February 2024. EJG
# 
# This script sources the necessary libraries and functions, coordinates the 
# analysis, and creates data structures that are then 'knitted' together (I think).
# So the idea is to run bits in here, and then 'Render' the RMD script. 
# Seems straightforward. :)
#
#
# Updates: 
# 2024/04/29: Steady development past few weeks; alll necessary pieces now developed. 
# 2024/04/29: Git repo created and working version pushed prior to RMarkdown development.
# 2024/05/02: Completed smoothing pass thru code; sorted raster plotting. Pushed.
# 2024/05/07: Another pass thru, adding some controls. Ready for RMD work.Pushed.

# TO DO: 
#  Design RMD report
#  Add config section to allow an RMD report to be built for selected extents.
#################################################################################

print('Starting Broughton ...')
rm(list=ls(all=T))  # Erase environment.

# Load necessary packages and functions ... 
source( "broughton_functions.R" )
# source( "Plot_Functions.R" )

# Processing FLAGS. When set to true, the data structure will be re-built from imported data. 
# Otherwise, it will be loaded from rData. Watch for dependencies.

loadtifs <- F
clipdata <- T # True if a spatial subset of the data is desired. Requires a polygon shape file. 
scaledat <- T # only needed if new layers loaded 


#---- Part 1 of 8: Load, clean, and display rasters for processing  ----
# If loadtifs == TRUE then run all this else load the processed data.

tif_stack <- stack()
dat_brick <- brick()

if (loadtifs) {
  print( "Loading predictors ... ")
  tif_stack <- LoadPredictors( raster_dir )

  today <- format(Sys.Date(), "%Y-%m-%d")
  save( tif_stack, file = paste0( data_dir, '/source_rasters_', today, '.rData' ))
  print( "Data saved. ")

  dat_brick <- tif_stack
  
  if (clipdata) {
#    maskme <- shapefile("C:\\Data\\SpaceData\\Broughton\\broughton_small.shp")
    maskme <- shapefile("C:\\Data\\SpaceData\\Broughton\\broughton_smaller.shp")
    tmp_stack <-  ClipPredictors( tif_stack, maskme )
    dat_brick <- tmp_stack
  } else {
    dat_brick <- tif_stack
  }
  
  if (scaledat){
    print( "Scaling TIFs ... ")
    tmp_stack <- ScalePredictors( dat_brick )
    dat_brick <- tmp_stack
  }

  } else {
  print( 'Loading project data ... ')
  load( paste0( data_dir, '/source_rasters_2024-05-27.rData' ))
}
 
# TIF and DAT stacks created.

### NOTE: All layers in the brick have the same geometry!
### NOTE: There is an interaction here with scaling. Perhaps integerization should be first?

names(tif_stack)
names(dat_brick)


# Intergerize scaled data to reduce data size and improve performance. 
prepped_layers <- Integerize( dat_brick )
save_brick <- brick( prepped_layers )
today <- format(Sys.Date(), "%Y-%m-%d")
save( save_brick, file = paste0( data_dir, '/prepped_brick_MSEA_QCS', today, '.rData' ))


prepped_layers <- stack( save_brick )


dim(prepped_layers)
names(prepped_layers)

plot( prepped_layers )
raster::hist(prepped_layers, nclass=50)
par(mfrow = c(1, 1))

# Quick correlation across data layers
x <- getValues( prepped_layers )
x_clean <- x[ complete.cases(x), ]
y <- cor( x_clean )
y

# Remove highly correlated layers.
  # Drop rugosity for now ... 
prepped_layers <- dropLayer( prepped_layers, c("rugosity" ))

prepped_layers <- dropLayer( prepped_layers, c("julBSpd_ave", "julBT_ave", "julSS_ave",
                        "julSS_min", "julSSpd_ave", "julST_ave", "julST_max", "northness", "tidal_cur") ) 

# All prepped_layers have same length of values 


#---- Part 2 of 8: Final data preparation ----


# Final data adjustments 

# Extract the data for the cluster analyses
stack_data   <- getValues( prepped_layers )

# remove any rows with an NA
# NB: This decouples the data from the RasterStack and 
#   needs to be put back together for some analyses below.
# ALSO: Blows up with a single layer in the stack 


clean_idx <- complete.cases(stack_data)
stack_data_clean <- stack_data[ clean_idx, ]

dim(stack_data)

#na_positions <- apply(stack_data, 1, anyNA) # takes a long time ... 
na_positions <- !clean_idx 

str(stack_data)
length(stack_data)
length(clean_idx)
length(na_positions)

SO ... complete.cases() appears to be working ... 




#---- Part 3 of 8: Explore number of clusters using Within-sum-of-squares plot ----
# Run multiple times with different number of clusters
# Initialize total within sum of squares errors: wss

set.seed <- 42 # Seed for reproducibility
randomz  <- 20 # the number of randomizations for kmeans to do.
imax     <- 25 # maximum iterations to try for convergence
nclust   <- 18 # number of clusters for scree plot

# Needs a subsample to run reasonably. Running all data is minutes per iteration ... 
#*** Deal with the warnings from kmeans
plotme <- MakeScreePlot( stack_data_clean, nclust, randomz, imax, 25000 )
plotme


#---- Part 4 of 8: Create and examine the clusters. ----
nclust <- 4 # the number of clusters based on Part 2, above.

# Perform k-means clustering. 500k good for exploration. 
# The full data set takes a few minutes with above randomz. 
sidx <- sample( 1:length( stack_data_clean[ , 1] ), 500000 )
samp <- stack_data_clean[ sidx, ]
cluster_result <- kmeans(samp, centers=nclust, nstart=randomz, iter.max=imax) 

# Summary of centres and relationship to predictors 
cluster_result$centers
x <- as.data.frame( cluster_result$centers )
for (i in names(x)){
  print( paste0( i, ":  ", round( range(x[i])[1], 2), " to ", round( range(x[i])[2], 2 ),
                "   Extent = ", round( abs( range(x[i])[1] - range(x[i])[2]), 2 )
        ))
}

# Expose and examine the cluster results
# NB: there are ~103M raster cells in the Broughton region. 
length( cluster_result$cluster )

#---- Part 5 of 8: Examine silhouette plot of the current cluster result ----
# Requires the data (i.e., cell attributes) and the corresponding assigned cluster
# Need to subsample from the cluster result above as distance matrix take long time.

# Calculating dissimilarity matrix (dist() below) takes a long time ... 
# 100000 way to much ... took 30ish min. 50k is ~1 min to completion. 
silx <- sample( 1:length( samp[ , 1] ), 10000 )

cs <- cluster_result$cluster[ silx ]
ss <- samp[ silx, ]

c_dist <- dist(ss)
sk <- silhouette(cs, c_dist) #also timeconsuming ... 

#plot(sk, border=NA )
par(mfrow = c(1, 1))
plot(sk, col = 1:nclust, border=NA )

#---- Part 6 of 8: Heat map of within-cluster standard deviations ----

# Define color palette
pal_heat <- rev( brewer.pal(n = nclust, name = "RdYlBu")) # heat map palette

# Create a smaller sample for some of the time-intensive subsequent steps  ... 
sidx <- sample( 1:length( stack_data_clean[ , 1] ), 10000 )
samp <- stack_data_clean[ sidx, ]
# re-run cluster for smaller sample.
cluster_result <- kmeans(samp, centers = nclust, nstart = randomz) # less than 10 seconds
csamp <- cluster_result$cluster

# Put the pieces together for the PCA by combining the data and the cluster. 
# Put cluster # first so easy to find for PCA below
profile_data <- as.data.frame( cbind(cluster = csamp, samp ) )

dim( profile_data )
names( profile_data )
sort( unique( profile_data$cluster ))

cluster_sd <- profile_data %>%
  group_by(cluster) %>%
  summarise_all(sd)

x <- as.data.frame( cluster_sd )
head(x)
xm <- melt( x, id.var = "cluster" )

z <- ggplot(xm, aes(x=cluster, y=variable, fill=value) ) +
  geom_tile() +
  scale_fill_gradientn(colours = pal_heat) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(title = "Within-cluster Standard Deviation", x = "Clusters", y = "Attributes", fill = "Value")
z


#---- Part 7 of 8: Show cluster groupings using PCA ----
# NB: This takes some time on the full data set. Uses profile_data from Part 5.
# Use profile data created from samples above. 

res.pca <- prcomp(profile_data[,-1],  scale = TRUE)
  # PC coordinates of individual raster cells
ind.coord <- as.data.frame(get_pca_ind(res.pca)$coord)
  # Add clusters from the classification
ind.coord$cluster <- factor( csamp )
# Data inspection
head(ind.coord)

# Percentage of variance explained by dimensions
eigenvalue <- round(get_eigenvalue(res.pca), 1)
variance.percent <- eigenvalue$variance.percent
eigenvalue

# Look at the clusters for the first 4 PCs
ggscatter(
  ind.coord, x = "Dim.1", y = "Dim.2", 
  color = "cluster", palette = "simpsons", ellipse = TRUE, ellipse.type = "convex",
  size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", variance.percent[1], "% )" ),
  ylab = paste0("Dim 2 (", variance.percent[2], "% )" )
) +
  stat_mean(aes(color = cluster), size = 4)

ggscatter(
  ind.coord, x = "Dim.3", y = "Dim.4", 
  color = "cluster", palette = "simpsons", ellipse = TRUE, ellipse.type = "convex",
  size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 3 (", variance.percent[3], "% )" ),
  ylab = paste0("Dim 4 (", variance.percent[4], "% )" )
) +
  stat_mean(aes(color = cluster), size = 4)


#---- Part 8 of 8: Spatialize the cluster results ----
# NB: Here have the option to re-cluster the entire data set for mapping.
#   Uses iter.max and nstart from above.

# Prepare a new raster for the cluster visualization
cluster_raster <- prepped_layers[[1]]
# dataType(cluster_raster) <- "INT1U" seems to fuck things up somehow ... 

length(cluster_result$cluster)

# Flag for re-clustering using all data prior to rebuilding.
allDat <- T
if (allDat) {
  # Cluster the entire data set for mapping ... 
    # less than 1 min with iter.max = 20, nstart = 20 for smallest region
  
  
  clean_idx <- complete.cases(stack_data)
  stack_data_clean <- stack_data[ clean_idx, ]
  na_positions <- !clean_idx
  
  cluster_result <- kmeans(stack_data_clean, centers = nclust, nstart = randomz, iter.max = imax)
  
  # Assign the clustered values ... 
  new_values <- values( cluster_raster )
  new_values[ clean_idx ] <- cluster_result$cluster
  values( cluster_raster ) <- new_values  

  raster::hist(cluster_result$cluster)
} else {
  # This part is a bit harder as need to deal with the sub-sampling ... 
  
  
}

range(cluster_result$cluster)
cellStats(cluster_raster, range)

# ???
  


# Plot the cluster raster ... 
# Define color palette
pal_clust <- brewer.pal(n = nclust, "Accent")

# trim and check the raster before plotting
a <- trim(cluster_raster)

# plot( a , col = pal_clust,
#      main = "K-means Clustering Result" )

ckey <- list( at=0:nclust, 
              title="Clusters", 
              col = pal_clust)
myTheme <- rasterTheme( region = pal_clust )
levelplot( a, margin = F, 
           colorkey = ckey,
           par.settings = myTheme,
           main = "K-means Clustering Result - Local extents" )


writeRaster( a, paste0( data_dir, "/foo.tif"), overwrite=TRUE)



#---- ----
#---- Outlier Detection ----
# Work with the full data set. 
cluster_result$centers

  # create a df w the centers for each cluster to facilitate distance calculation.
centers <- cluster_result$centers[cluster_result$cluster, ] 
head(centers)
  # calculate distance
distances <- sqrt(rowSums((stack_data_clean - centers )^2))
head(x)
outlier_idx <- order(distances, decreasing=T)[1:1000000]

#print(stack_data_clean[outlier_idx,])

# Redo the cluster work with the outliers removed ... 
dim( stack_data_clean )
stack_data_clean <- stack_data_clean[ -outlier_idx,]
dim( stack_data_clean )

# Plot the outliers.
  # Note: These outliers are in cluster space, but the plots below are in the original 
  # space of the scaled variables. And while there are clearly a few outliers in the 
  # scaled variables (for which we should do a page of box plots), its the combination 
  # that creates ouitliers in the cluster analysis. Interesting.

  # NB: Testing clustering with outliers removed (up to 1M) creates a bit less dispersion 
  # but the shapes remain. 

# a subsample to reduce the plot size ... 
ssidx <- sample( 1:length( stack_data_clean[ , 1] ), 1000 )

plot(stack_data_clean[ ssidx, c("bathymetry", "fetch") ], pch=19, col=cluster_result$cluster[ ssidx ], cex=1)
points(cluster_result$centers[ ,c("bathymetry", "fetch") ], col=1:3, pch=15, cex=2)
#points(stack_data_clean[ outlier_idx, c("bathymetry", "fetch") ], pch="+", col=4, cex=3)

p_axes <- c("rugosity", "standard_deviation_slope") 
plot(stack_data_clean[ ssidx, p_axes ], pch=19, col=cluster_result$cluster[ ssidx ], cex=1)
points(cluster_result$centers[ , p_axes ], col=1:3, pch=15, cex=2)
points(stack_data_clean[ outlier_idx, p_axes ], pch="+", col=4, cex=3)

p_axes <- c("rugosity", "tidal") 
plot(stack_data_clean[ ssidx, p_axes ], pch=19, col=cluster_result$cluster[ ssidx ], cex=1)
points(cluster_result$centers[ , p_axes ], col=1:3, pch=15, cex=2)
points(stack_data_clean[ outlier_idx, p_axes ], pch="+", col=4, cex=3)




#---- Knit and render Markdown file -----
### Process file 
# To HTML ... 
rmarkdown::render( "Broughton.Rmd",   
                   output_format = 'html_document',
                   output_dir = rmd_dir )  

# 2024/04/29: It looks like this has gotten easier in the last 2 years ... version up!

# To PDF:
# First had to install the library tinytex.
# then run >tinytex::install_tinytex()
# ... and done. 
rmarkdown::render( "Broughton.Rmd",   
                   output_format = 'pdf_document',
                   output_dir = rmd_dir )  


#---- Some details on correlation analysis ... ----
#-- Correlation across UN-scaled data layers ... 
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

#---- ----

# FIN.






