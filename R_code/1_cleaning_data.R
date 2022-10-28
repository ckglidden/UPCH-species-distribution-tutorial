###cleaning mapbiomas data for SDM

#load libraries
library(tidyr); library(dplyr); library(PerformanceAnalytics)

#-----------------------------------#
#read in datasets                   #
#-----------------------------------#

mapbiomas <- read.csv("data/a_chamek_ter_mammals_lulc_Oct2022.csv")

#-----------------------------------#
#update label MAPBIOMAS classes     #
#-----------------------------------#

#you can look at the classes included in the data using:
unique(mapbiomas$class)

#relabel each class to make it easier to see results
mapbiomas$class[mapbiomas$class == 3] <- "forest_formation"
#mapbiomas$class[mapbiomas$class == 4] <- "savannah_formation"
mapbiomas$class[mapbiomas$class == 6] <- "flooded_forest"
mapbiomas$class[mapbiomas$class == 11] <- "wetland"
#mapbiomas$class[mapbiomas$class == 12] <- "grassland"
#mapbiomas$class[mapbiomas$class == 13] <- "non_forested_natural"
mapbiomas$class[mapbiomas$class == 14] <- "farming"
mapbiomas$class[mapbiomas$class == 24] <- "urban"
#mapbiomas$class[mapbiomas$class == 25] <- "other_non_vegetated"
#mapbiomas$class[mapbiomas$class == 27] <- "not_observed"
#mapbiomas$class[mapbiomas$class == 29] <- "rocky_outcrop"
mapbiomas$class[mapbiomas$class == 30] <- "mining"
mapbiomas$class[mapbiomas$class == 33] <- "river_lake_ocean"
#mapbiomas$class[mapbiomas$class == 34] <- "glacier"

#any lulc that was not defined by classIDs in you MAPBIOMAS code will be a zero here, remove these datapoints
mapbiomas <- mapbiomas[mapbiomas$class != 0, ]

#----------------------------------------------------------#
#go from wide to long so each class is a unique column     #
#----------------------------------------------------------#

mapbiomas_wide <- mapbiomas %>% 
  pivot_wider(names_from = class, values_from = area) 

#change NAs to zero as NA means the landclass is not present
mapbiomas_wide[is.na(mapbiomas_wide)] <- 0


#---------------------------------------------------------#
#summarize average area per class per point across years  #
#and difference in area over the study period             #
#---------------------------------------------------------#

#there prob is a more effecient way to loop through columns...
mapbiomas_mean_diff <- mapbiomas_wide %>%
                  group_by(row_code) %>%
                  summarise(#mean per class
                            mean_forest = mean(forest_formation),
                            mean_river_lake_ocean = mean(river_lake_ocean),
                            mean_farming = mean(farming),
                            mean_flooded_forest = mean(flooded_forest),
                            mean_wetland = mean(wetland),
                            mean_urban = mean(urban),
                            mean_mining = mean(mining),
                            #difference over study period per class
                            diff_forest_formation = forest_formation[match(2020, year)] - forest_formation[match(2001, year)],
                            diff_farming = farming[match(2020, year)] - farming[match(2001, year)],
                            diff_flooded_forest = flooded_forest[match(2020, year)] - flooded_forest[match(2001, year)],
                            diff_wetland = wetland[match(2020, year)] - wetland[match(2001, year)],
                            diff_urban = urban[match(2020, year)] - urban[match(2001, year)])

#----------------------------------------------------------#
#pearsons correlation analysis                             #
#----------------------------------------------------------#

# correlation in absolute terms
corr <- abs(cor(mapbiomas_mean_diff[2:ncol(mapbiomas_mean_diff)])) 

#correlation plot with points on lower left of matrix, correlation coeffecients on upper right, and distributions on the diagonal
chart.Correlation(mapbiomas_mean_diff[2:ncol(mapbiomas_mean_diff)], 
                  histogram = TRUE, method = "pearson")

#nothing with correlation > 0.7

write.csv(mapbiomas_mean_diff, "data/a_chamek_ter_mammals_lulc_cleaned_Oct2022.csv")





