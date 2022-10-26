####random forest species distribution models
library(tidyr); library(dplyr); library(spatialsample); library(sf)

#------------------------------------------------------------------------#
#read in data & merge by row id  (or merge from output of above code)    #
#------------------------------------------------------------------------#

lulc <- read.csv("data/b_tridactylus_ter_mammals_lulc_cleaned_Oct2022.csv")
amazon_basin_pnts <-  read.csv("data/b_tridactylus_ter_mammals_amazon_thinned_Oct22.csv")

data0 <- left_join(amazon_basin_pnts, lulc, by = "row_code")

#------------------------------------#
#get fold id by k-means clustering   #
#------------------------------------#

#convert analysis_data to a spatial object
data0_sf <- st_as_sf(x = data0, 
                             coords = c("lon", "lat"),
                             crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

#identify groups of 5 clusters using the spatialsample package
clusters <- spatial_clustering_cv(data0_sf, v = 5) #k-means clustering to identify 5 cross-validation folds

#for loop to create a dataframe that assigns a fold number to each data point
splits_df <- c()
for(i in 1:5){
  new_df <- assessment(clusters$splits[[i]]) #extract points in fold number i
  new_df$fold <- i
  new_df <- new_df[,c("row_code", "fold")]
  splits_df <- rbind(splits_df, new_df) #bind all points x fold id together
}

splits_df <- st_drop_geometry(splits_df) #drop shapefiles

#final data - merge cluster id to final dataset for analysis
analysis_data <- merge(data0, splits_df, by = "row_code")

#sanity check: check how many data points are in each fold
table(analysis_data$fold)

write.csv(analysis_data, "data/b_tridactylus_ter_mammals_finalData_Oct22.csv")
