# Exploring various plot options. 
# This is not ready for prime time code. 


#---- Alternatives for plotting cluster maps

legend_labels <- c("Ca1", "C2", "C3", "C4", "C5", "C6")

# Define custom legend colors
#k = 6
legend_colors <- c("blue", "yellow", "purple", "red", "pink", "green")

# Plot the cluster raster with custom legend

rasterVis::levelplot(cluster_raster, col.regions = rainbow(k),
                     main = "K-means Clustering Result")

rasterVis::levelplot(cluster_raster, col.regions = legend_colors,
                     main = "K-means Clustering Result",
                     #                     at = 1:length(legend_labels),
                     at = 1:6,
                     colorkey = list(labels = list(at = 1:length(legend_labels), labels = legend_labels)))

plot( cluster_raster )
legend("topright", legend = sort( unique(cluster.result$cluster )), fill = 1:length(sort(unique(cluster.result$cluster))) )





#---- Plot region map showing centroids ...
# e.g.,

cluster_data <- data.frame(
  ClusterID = c(1, 2, 3),
  Longitude = c(-122, -121, -123),
  Latitude = c(37, 38, 36)
)

# Convert data to sf object
cluster_sf <- st_as_sf(cluster_data, coords = c("Longitude", "Latitude"))

# Plot cluster centroid map
ggplot() +
  geom_sf(data = cluster_sf, aes(color = as.factor(ClusterID))) +
  geom_point(data = cluster_data, aes(x = Longitude, y = Latitude, color = as.factor(ClusterID)), size = 3) +
  labs(title = "Cluster Centroid Map", x = "Longitude", y = "Latitude") +
  scale_color_discrete(name = "Cluster ID") +
  theme_minimal()


stack_data_matrix <- as.matrix( stack_data_clean )
stack_sample <- stack_data_matrix[1:1000, ]

str(stack_data_matrix)

hcluster <- hclust( dist( stack_sample ), method = "centroid", members = NULL )
dendro <- as.dendrogram( hcluster )
plot( dendro )



# Plotting as rasters ... 
rmat <- as.matrix( cluster_stack )
# Remove rows with NA values (assuming NA values represent no data)
cleaned_matrix <- rmat[complete.cases( rmat), ]
# Convert matrix back to raster
cleaned_raster <- raster(cleaned_matrix, xmn=min(r), xmx=max(r), ymn=min(r), ymx=max(r), crs=crs(r))
# Plot the cleaned raster
plot(cleaned_raster)
