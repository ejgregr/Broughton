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

loadtifs <- T # FROM the specified raster data directory
clipdata <- F # Restrict the spatial extents of the tifs using a supplied polygon shape file mask
scaledat <- T # only needed if new layers loaded 

#---- Part 1 of 3: Load, clean, and prepare predictor data.  ----
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
  save( tif_stack, file = paste0( data_dir, '/tifs_MSEA_scaled_QCS_', today, '.rData' ))
  
} else {
  print( 'Loading project data ... ')
  # Ideally meaningfully named and tested so no source is required.
  load( paste0( data_dir, '/tifs_MSEA_scaled_smMask_2024-05-30.rData' ))
}
 

#---- Final data preparation ----
#NOTE: Everything from here down is in scaled units. 

#-- Intergerize scaled data to reduce data size and improve performance. 
# NOTE: 2024/06/03: Integerize throws a warning (layers with no data) on the MSEA data set.
#   but no problem evident in results. Less useful with trimmed data, but may still have 
#   value for larger data sets.
#prepped_layers <- Integerize( tif_stack )
prepped_layers <- tif_stack

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

#-- Remove any unwanted layers.
  # Drop rugosity as correlated with SDSlope 
prepped_layers <- dropLayer( prepped_layers, "rugosity" )
names( prepped_layers )


#-- Finalize the data for exploring and creating the clusters.
  # Extract the data for the cluster analyses
stack_data <- getValues( prepped_layers )
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
nsample  <- 25000 # Scree plot needs a subsample to run reasonably. 
plotme <- MakeScreePlot( stack_data_clean, nclust, randomz, imax, nsample )
plotme

#---- Create a working set of N clusters (N based on scree plot) to further assess cluster number. ----

nclust  <- 8 # the number of clusters based on scree plot, above.
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

z <- ggplot(xm, aes(x=cluster, y=variable, fill=value) ) +
  geom_tile() +
  scale_fill_gradientn(colours = pal_heat) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(title = "Within-cluster Standard Deviation", x = "Clusters", y = "Attributes", fill = "Value")
z

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

#plot(sk, border=NA )
par(mfrow = c(1, 1))
plot(sk, col = 1:nclust, border=NA )


#---- Part 3: Detailed examination of N clusters  ----
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

reclust = F

# initialize target data structure 
cluster_raster <- prepped_layers[[1]]
dataType(cluster_raster) <- "INT1U"
cluster_raster[] <- NA

if (reclust == F) {

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
  
  #--- Two steps here: First assign clusters to the rest of the clean stack data, 
  #     THEN put the clean data back in the raster.
  # Uses samp and sidx from lines ~150 above.
  # Find the things not in the sample (using sidx from line ~150 above).

  not_sidx <- setdiff( 1:dim(stack_data_clean)[[1]], sidx )
  
  # the data to predict clusters for
  new_dat <- stack_data_clean[ not_sidx, ]
  # Process the new data in chunks of 10k to ensure performance
  chunk_size <- 10000
  n_chunks <- ceiling(nrow(new_dat) / chunk_size)
  
  print( paste0( "predicting clusters for ", dim(new_dat)[[1]], " pixels using ", 
                 n_chunks, " chunks of ", chunk_size, ". Stand by ... ") )
  
  # Where we are putting our new chunks
  predicted_clusters <- vector("integer", nrow(new_dat))
  
  for ( i in seq_len(n_chunks) ) {
    start_index <- (i - 1) * chunk_size + 1
    end_index <- min(i * chunk_size, nrow(new_dat))
    
    chunk <- new_dat[start_index:end_index, ]
    predicted_clusters[start_index:end_index] <- PredictClusters(chunk, cluster_result)
    cat(i, " ")
  }
  
  length(predicted_clusters)
  length(cluster_result$cluster)
  # Combine the results for the cleaned data ... 
  clean_clusts <- array( 1:dim(stack_data_clean)[1] )
  clean_clusts[ sidx ] <- cluster_result$cluster
  clean_clusts[ not_sidx ] <- predicted_clusters
  
  # extract values from the target cluster
  new_values <- values( cluster_raster )
  # replace clustered values with the cluster results
  new_values[ clean_idx ] <- clean_clusts
  # reassign to the plotting raster
  values( cluster_raster ) <- new_values 
  
}


#--- Display the results, first as histogram then as map.  
raster::hist( values(cluster_raster ))

# Define color palette
pal_clust <- brewer.pal(n = nclust, "Accent")

ckey <- list( at=0:nclust, 
              title="Clusters", 
              col = pal_clust)
myTheme <- rasterTheme( region = pal_clust )
levelplot( cluster_raster, margin = F, 
           colorkey = ckey,
           par.settings = myTheme,
           main = "K-means Clustering Result - Local extents" )


writeRaster( cluster_raster, paste0( data_dir, "/QCS_8cluster.tif"), overwrite=TRUE)


#---- 

#----Outlier Detection ----
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






