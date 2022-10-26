###cleaning data for SDM

#load libraries
library(tidyr); library(dplyr); library(PerformanceAnalytics)

#-----------------------------------#
#read in datasets                   #
#-----------------------------------#

occ_data <- read.csv("data/b_tridactylus_ter_mammals_amazon_thinned_Oct22.csv")
mapbiomas <- read.csv("data/b_tridactylus_ter_mammals_lulc_Oct2022.csv")
#human_population <- read.csv("...") skipping this for now

#-----------------------------------#
#update label MAPBIOMAS classes     #
#-----------------------------------#

#you can look at the classes included in the data using:
#unique(mapbiomas$class)

#relabel each class to make it easier to see results
mapbiomas$class[mapbiomas$class == 3] <- "forest_formation"
mapbiomas$class[mapbiomas$class == 4] <- "savannah_formation"
mapbiomas$class[mapbiomas$class == 5] <- "mangrove"
mapbiomas$class[mapbiomas$class == 6] <- "flooded_forest"
mapbiomas$class[mapbiomas$class == 11] <- "wetland"
mapbiomas$class[mapbiomas$class == 12] <- "grassland"
mapbiomas$class[mapbiomas$class == 13] <- "non_forested_natural"
mapbiomas$class[mapbiomas$class == 14] <- "farming"
mapbiomas$class[mapbiomas$class == 24] <- "urban"
mapbiomas$class[mapbiomas$class == 25] <- "other_non_vegetated"
mapbiomas$class[mapbiomas$class == 27] <- "not_observed"
mapbiomas$class[mapbiomas$class == 29] <- "rocky_outcrop"
mapbiomas$class[mapbiomas$class == 30] <- "mining"
mapbiomas$class[mapbiomas$class == 33] <- "river_lake_ocean"

#----------------------------------------------------------#
#summarize average area per class per point across years   #
#----------------------------------------------------------#

mapbiomas_mean <- mapbiomas %>%
                  group_by(row_code, class) %>%
                  summarise(mean_area = mean(area))

#----------------------------------------------------------#
#go from wide to long so each class is a unique column     #
#----------------------------------------------------------#

mapbiomas_mean_wide <- mapbiomas_mean %>% 
  pivot_wider(names_from = class, values_from = mean_area) 

#change NAs to zero as NA means the landclass is not present
mapbiomas_mean_wide[is.na(mapbiomas_mean_wide)] <- 0


#----------------------------------------------------------#
#summarize average population per point across years       #
#----------------------------------------------------------#

##skipping this for now

#----------------------------------------------------------#
#pearsons correlation analysis                             #
#----------------------------------------------------------#

# correlation in absolute terms
corr <- abs(cor(mapbiomas_mean_wide[2:ncol(mapbiomas_mean_wide)])) 

#correlation plot with points on lower left of matrix, correlation coeffecients on upper right, and distributions on the diagonal
chart.Correlation(mapbiomas_mean_wide[2:ncol(mapbiomas_mean_wide)], 
                  histogram = TRUE, method = "pearson")


write.csv(mapbiomas_mean_wide, "data/b_tridactylus_ter_mammals_lulc_cleaned_Oct2022.csv")





