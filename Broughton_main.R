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
# 2024/05/29: After some exploratory work with Romina's tifs from SPECTRAL lab, forked the 
#   code so this version focuses on DFO data/objectives and Broughton2 attends to the LSSM needs.

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

loadtifs <- T
clipdata <- T # Restrict the spatial extents of the tifs using a supplied polygon shape file mask
scaledat <- T # only needed if new layers loaded 

#---- Part 1 of 8: Load, clean, and display rasters for processing  ----
# If loadtifs == TRUE then run all this else load the processed data.

tif_stack <- stack()
today <- format(Sys.Date(), "%Y-%m-%d")

if (loadtifs) {
  print( "Loading predictors ... ")
  src_stack <- LoadPredictors( raster_dir )
  print( "Data loaded.")

  tif_stack <- src_stack
  
  if (clipdata) {
    print( "clipping TIFs ... ")
    #amask <- shapefile("C:\\Data\\SpaceData\\Broughton\\broughton_small.shp")
    amask <- shapefile("C:\\Data\\SpaceData\\Broughton\\broughton_smaller.shp")
    tif_stack <- ClipPredictors( tif_stack, amask )
    print('Rasters clipped.')
  }

  # The landside on the MSEA bathymetry is not of interest to the classification.
  # Deeper areas are also unsuitable for kelps. This removes those  pixels from 
  # the affected rasters. HARD-CODED! 
  tif_stack <- DropNonHabitat( tif_stack )
  print('Unsuitable elevations removed.')

  if (scaledat) {
    print( "Scaling TIFs ... ")
    tmp_stack <- scale( tif_stack )
    tif_stack <- tmp_stack
    print('Rasters Scaled.')
  }

  save( src_stack, file = paste0( data_dir, '/src_stack_', today, '.rData' ))
  save( tif_stack, file = paste0( data_dir, '/tif_stack_', today, '.rData' ))
  save( tif_stack, file = paste0( data_dir, '/tifs_MSEA_scaled_smMask_', today, '.rData' ))
  
} else {
  print( 'Loading project data ... ')
  # Ideally meaningfully named and tested so no source is required.
  load( paste0( data_dir, '/tifs_MSEA_scaled_smMask_2024-05-29.rData' ))
}
 

plot( tif_stack )

a <- tif_stack[[1]]
writeRaster( a, paste0( data_dir, "/foo.tif"), overwrite=TRUE)



#-- Intergerize scaled data to reduce data size and improve performance. 
  # NOTE: Throws an error on the lage MSEA data set ... no problem evident in the results. 
prepped_layers <- Integerize( tif_stack )

dim(prepped_layers)
names(prepped_layers)


#-- Visualize source data
plot( prepped_layers )
raster::hist(prepped_layers, nclass=50)
par(mfrow = c(1, 1))

#-- Quick correlation across data layers
x <- getValues( prepped_layers )
x_clean <- x[ complete.cases(x), ]
y <- cor( x_clean )
y

#-- Remove any unwanted layer.
  # Drop rugosity as correlated with SDSlope 
prepped_layers <- dropLayer( prepped_layers, "rugosity" )
names( prepped_layers )


#---- Part 2 of 8: Final data preparation ----

# Extract the data for the cluster analyses
stack_data <- getValues( prepped_layers )

# remove any rows with an NA
# NB: This decouples the data from the RasterStack and requires re-assembly for mapping
clean_idx <- complete.cases(stack_data)
stack_data_clean <- stack_data[ clean_idx, ]

dim( stack_data )
dim( stack_data_clean )

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
nclust <- 10 # the number of clusters based on Part 2, above.

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
cluster_result <- kmeans(samp, centers = nclust, nstart = radomz) # less than 10 seconds
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

allDat <- F

if (allDat) {
# Prepare a zeroed raster for the cluster assignment
  cluster_raster <- prepped_layers$bathymetry
  dataType(cluster_raster) <- "INT1U"
  cluster_raster[] <- NA
  
  # Cluster the entire data set for mapping ... 
    # less than 1 min with iter.max = 20, nstart = 20 for smallest region
  cluster_result <- kmeans(stack_data_clean, centers = nclust, nstart = randomz, iter.max = imax)
  
  # Assign the clustered values ... 
  cluster_raster[ complete.cases( stack_data ), ] <- cluster_result$cluster
  raster::hist(cluster_result$cluster)
} else {
  
  cluster_raster <- prepped_layers$bathymetry
  dataType(cluster_raster) <- "INT1U"
  cluster_raster[] <- NA
  
  # Assign the clustered values ... 
  cluster_raster[ complete.cases( stack_data ), ] <- cluster_result$cluster
  #raster::hist(cluster_result$cluster)
  
}
# Define color palette
pal_clust <- brewer.pal(n = nclust, "Accent")

# trim and check the raster before plotting
a <- trim(cluster_raster)
unique( values( a$bathymetry ))

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

#---- ----


# FIN.






