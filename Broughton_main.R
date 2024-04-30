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
# Update: 2024/04/30 EJG
# Steady development has led to all the pieces being largely develped. 
#################################################################################

print('Starting Broughton ...')
rm(list=ls(all=T))  # Erase environment.

# Load necessary packages and functions ... 
source( "broughton_functions.R" )
# source( "Plot_Functions.R" )


#----------------------------
### Process file 
# To HTML ... 
rmarkdown::render( "Broughton.Rmd",   
                   output_format = 'html_document',
                   output_dir = results.dir )  

# 2024/04/29: It looks like this has gotten easier in the last 2 years ... version up!

# To PDF:
# First had to install the library tinytex.
# then run >tinytex::install_tinytex()
# ... and done. 
rmarkdown::render( "Broughton.Rmd",   
                   output_format = 'pdf_document',
                   output_dir = results.dir )  
#-------------------------------------------------





#---- Part 1 of 4: Load and prepare rasters for classification ----

# Processing FLAGS. When set to true, the data structure will be re-built from imported data. 
# Otherwise, it will be loaded from rData. Watch for dependencies.

loadtifs  <- F
scaleData <- F
# trimLand  <- F ## Not sure we need this ...

  # Load data ... 
if (loadtifs == T){

  print( "Loading predictors ... ")
  print( list.files(path = raster_dir, pattern = '\\.tif$', full.names = FALSE) )
  
  data_layers <- Load.Predictors( raster_dir )
  print( "Data loaded ... ")
  
  today <- format(Sys.Date(), "%Y-%m-%d")
  save( data_layers, file = paste0( data_dir, '/source_rasters_', today, '.rData' ))
  print( "Data saved ... ")
  
} else {
    load( paste0( data_dir, '/source_rasters_2024-04-12.rData' ))}


# Trim the land part away. Seems desirable for most applications ... 
# BUT WHY?? Can't be for the clustering cuz that removes NAs ... 
# THINK!!

names(data_layers)

# if (trimland == T){
#   trim_data <- getValues( data_layers$bathymetry )
#   trim_idx <- trim_data < -5 #for rugos and stdevslope
#   trim_data[ trim_idx ] <- NA
#   data_layers$bathymetry <- setValues( data_layers$bathymetry, as.integer( trim_data ))
#   data_layers$bathymetry <- setMinMax( data_layers$bathymetry )
#   
#   trim_data <- getValues( data_layers$rugosity )
#   trim_data[ trim_idx ] <- NA
#   data_layers$rugosity <- setValues( data_layers$rugosity, trim_data )
#   data_layers$rugosity <- setMinMax( data_layers$rugosity )
#   
#   trim_data <- getValues( data_layers$standard_deviation_slope )
#   trim_data[ trim_idx ] <- NA
#   data_layers$standard_deviation_slope <- setValues( data_layers$standard_deviation_slope, trim_data )
#   data_layers$standard_deviation_slope <- setMinMax( data_layers$standard_deviation_slope )
#   print( "Data trimmed ... ")
#   
#   today <- format(Sys.Date(), "%Y-%m-%d")
#   save( data_layers, file = paste0( data_dir, '/trimmed_rasters_', today, '.rData' ))
#   print( "Data saved ... ")
#   
# } else {
#   load( paste0( data_dir, '/trimmed_rasters_2024-04-12.rData' ))}

  # Now scale the data for analysis and save.
if (scaleData == T){

  print('scaling ...')
  scaled_layers <- scale( data_layers )
  print('saving ...')
  today <- format(Sys.Date(), "%Y-%m-%d")
  names( scaled_layers ) <- names( data_layers )
  save( scaled_layers, file = paste0( data_dir, '/scaled_rasters_', today, '.rData' ))

} else {
  load( paste0( data_dir, '/scaled_rasters_2024-04-26.rData' )) }

# Intergerize scaled data to see if performance improves ... 
# Also turns the scaled brick back into a RasterStack.
scaled_ilayers <- Integerize( scaled_layers )

# Data inspection - NB can't plot a RasterBrick, apparently
names(data_layers)
plot( data_layers )
raster::hist(data_layers, nclass=50)

plot( scaled_ilayers )
raster::hist(scaled_ilayers, nclass=50)
par(mfrow = c(1, 1))

  # Quick correlation across data layers
x <- getValues( scaled_ilayers )
x_clean <- x[ complete.cases(x), ]
y <- cor( x_clean )
y

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


  # Final data prep

# Select StdDevSlope  over Rugosity for now ... 
layers_to_use <- c("bathymetry", "fetch", "standard_deviation_slope", "tidal")
#layers_to_use <- c("bathymetry", "fetch", "rugosity", "tidal")
scaled_ilayers <- subset(scaled_ilayers, layers_to_use )
names(scaled_ilayers)


#--------------------------
### Spatial subsetting here

sMask <- shapefile("C:\\Data\\SpaceData\\Broughton\\broughton_small.shp")
smMask <- shapefile("C:\\Data\\SpaceData\\Broughton\\broughton_smaller.shp")
#aMask   <- boundaries[[1]]
str(sMask)

# Data extents No mask and polygons
x  <- scaled_ilayers[[1]]
#Regional
x  <- raster::mask(scaled_ilayers[[1]], sMask )
# Local
x  <- raster::mask(scaled_ilayers[[1]], smMask )

plot(trim(x), main = 'Local')
sp::plot( sMask, border = "darkgreen", lwd = 2, add=TRUE)
sp::plot( smMask, border = "red", lwd = 2, add=TRUE)



# Extract the data for the cluster analyses
x <- raster::mask(scaled_ilayers, smMask )
stack_data <- getValues( x )

# remove any rows with an NA
clean_idx <- complete.cases(stack_data)
stack_data_clean <- stack_data[ clean_idx, ]





#---- Part 2 of 5: Explore number of clusters. ----
#---- Within sum of squares plot - cluster number selection. ----
# Run multiple times with different number of clusters
# Initialize total within sum of squares errors: wss

# Needs a subsample to run reasonably ... 
sidx <- sample( 1:length( stack_data_clean[ , 1] ), 50000 )
stack_sample <- stack_data_clean[ sidx, ]
clust_sample <- cluster_result$cluster[ sidx ]

n <- 18
wss <- numeric(n)

# Look over 1 to n possible clusters
for (i in 1:n) {
  
  ss <- clust_sample #stack_data_clean
  
  # Fit the model: km.out
  km.out <- kmeans(ss, centers = i, nstart = 20)
  
  # Save the within cluster sum of squares
  wss[i] <- km.out$tot.withinss
}

# Produce a scree plot
wss_df <- tibble(clusters = 1:n, wss = wss)

scree_plot <- ggplot(wss_df, aes(x = clusters, y = wss, group = 1)) +
  geom_point(size = 4)+
  geom_line() +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
  xlab('Number of clusters') +
  ylab('Total within-cluster sum of squares') +
  geom_hline(
    yintercept = wss, 
    linetype = 'dashed')

scree_plot




#---- Part 3 of 5: Create the k-means clusters. ----

 # Number of clusters or from above. 
k <- 8 # Based on scree plot of total sum of squares

# Perform k-means clustering. Number of clusters in supporting script.
# A few randomizations in the start makes a big difference. Note seed. 
cluster_result <- kmeans(stack_data_clean, centers = k, nstart = 5) # less than 10 seconds

#names(cluster_result)
cluster_result$centers
#format(cluster_result$tot.withinss, scientific = TRUE)

# Expose and examine the cluster results
# NB: there are ~103M raster cells in the Broughton region. 
length( cluster_result$cluster )
range( cluster_result$cluster )

#--- Spatialize the cluster results ... 
  # Prepare a zeroed raster for the cluster assignment
cluster_raster <- scaled_ilayers$bathymetry
dataType(cluster_raster) <- "INT1U"
cluster_raster[] <- NA

  # Assign the clustered values ... 
cluster_raster[ complete.cases( stack_data ), ] <- cluster_result$cluster
#raster::hist(cluster_result$cluster)

### 2024/04/02: Pretty sure I have a good cluster raster at this point. 

#---- Unsatisfying original rasterVis plots ----

plot(cluster_raster,
     main = "K-means Clustering Result" )

rasterVis::levelplot(cluster_raster,
                     main = "K-means Clustering Result" )


#----  Cluster profiling ----

# Create a sample for some of the time-intensive subsequent steps  ... 
sidx <- sample( 1:length( stack_data_clean[ , 1] ), 10000 )
stack_sample <- stack_data_clean[ sidx, ]
clust_sample <- cluster_result$cluster[ sidx ]
head(clust_sample)

# Combine the data and the cluster - complete cases only. 
# Put cluster # first so easy to remove for PCA below
profile_data <- as.data.frame( cbind(cluster = clust_sample, stack_sample ) )
dim( profile_data )
names( profile_data )
head(profile_data)

cluster_sd <- profile_data %>%
  group_by(cluster) %>%
  summarise_all(sd)

x <- as.data.frame( cluster_sd )
x <- x[, -4] # leave out rugosity ... 
head(x)
xm <- melt( x, id.var = "cluster" )

z <- ggplot(xm, aes(x=cluster, y=variable, fill=value) ) +
  geom_tile() +
  scale_fill_gradientn(colours = pal_heat) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(title = "Within-cluster Standard Deviation", x = "Clusters", y = "Attributes", fill = "Value")
z

# Show cluster groupings
  # NB: This takes some time on the full data set. 
  # Use profile data created from samples above. 

  # Dimension reduction using PCA
res.pca <- prcomp(profile_data[,-1],  scale = TRUE)
  # PC coordinates of individual raster cells
ind.coord <- as.data.frame(get_pca_ind(res.pca)$coord)
  # Add clusters from the classification
ind.coord$cluster <- factor( clust_sample )
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




# Silhouette plot ... 

# Use a smaller sample ... 
sx <- sample( 1:length( stack_sample[ , 1] ), 1000 )

cs <- clust_sample[ sx ]
ss <- stack_sample[ sx, ]

# Need a dissimilarity matrix - this takes a long time ... 
c_dist <- dist(ss)

sk <- silhouette(cs, c_dist) 
str(sk)
str( sk[ 1:100, ])

plot(sk)

plot(sk, nmax.lab = 40, max.strlen = 5,
     main = NULL, sub = NULL, xlab = expression("Silhouette width "* s[i]),
     col = 2:(k + 1), do.col.sort = length(col) > 1, border = 0,
#     col = "blue", do.col.sort = length(col) > 1, border = 0,
     cex.names = par("cex.axis"), do.n.k = TRUE, do.clus.stat = TRUE)


fviz_silhouette(sk)


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





# FIN.


