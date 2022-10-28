###cleaning mapbiomas data for SDM

#load libraries
library(tidyr); library(dplyr); library(PerformanceAnalytics)

#-----------------------------------#
#read in datasets                   #
#-----------------------------------#

mapbiomas <- read.csv("data/a_chamek_ter_mammals_lulc_Oct2022.csv")
climate <- read.csv("data/a_chamek_ter_mammals_climate_Oct2022.csv"); climate <- climate[  c("row_code", "bio13_precip_wettest_month", "cmi_min")]
amazon_basin_pnts <-  read.csv("data/a_chamek_ter_mammals_amazon_thinned_Oct22.csv")

#---------------------------------------#
#update label for MAPBIOMAS classes     #
#---------------------------------------#

#you can look at the classes included in the data using:
unique(mapbiomas$class)

#relabel each class to make it easier to see results
mapbiomas$class[mapbiomas$class == 3] <- "forest_formation"
mapbiomas$class[mapbiomas$class == 14] <- "farming"

#any lulc that was not defined by classIDs in you MAPBIOMAS code will be a zero here, remove these rows
mapbiomas <- mapbiomas[mapbiomas$class != 0, ]

#----------------------------------------------------------#
#go from wide to long so each class is a unique column     #
#----------------------------------------------------------#

mapbiomas_wide <- mapbiomas %>% 
  pivot_wider(names_from = class, values_from = area) 

#change NAs to zero as NA means the land class is not present (so there is 0 area)
mapbiomas_wide[is.na(mapbiomas_wide)] <- 0


#---------------------------------------------------------#
#summarize average area per class per point across years  #
#and difference in area over the study period             #
#---------------------------------------------------------#

mapbiomas_mean_diff <- mapbiomas_wide %>%
                  group_by(row_code) %>%
                  summarise(
                            #mean per class
                            mean_forest = mean(forest_formation),
                            mean_farming = mean(farming),
                            
                            #difference over study period per class
                            diff_forest_formation = forest_formation[match(2020, year)] - forest_formation[match(2001, year)])

#save dataframe
write.csv(mapbiomas_mean_diff, "data/a_chamek_ter_mammals_lulc_cleaned_Oct2022.csv")

#--------------------------------------------------------------------------------------------------------#
#merge data frames with response (presence / absence labels), climate covariates, and lulc covariates    #
#--------------------------------------------------------------------------------------------------------#

covariates <- left_join(climate, mapbiomas_mean_diff, by = "row_code")
data0 <- left_join(amazon_basin_pnts, covariates, by = "row_code")

#----------------------------------------------------------#
#pearsons correlation analysis                             #
#----------------------------------------------------------#

# correlation in absolute terms
corr <- abs(cor(data0[7:ncol(data0)])) 

#correlation plot with points on lower left of matrix, correlation coeffecients on upper right, and distributions on the diagonal
chart.Correlation(data0[7:ncol(data0)], 
                  histogram = TRUE, method = "pearson")

#nothing with correlation > 0.7

#--------------------------------------------#
#get  fold  id  by  block  clustering      #
#--------------------------------------------#

#convert  analysis_data  to  a  spatial  object
data0_sf  <-  st_as_sf(x  =  data0,
                       coords  =  c("lon", "lat"),
                       crs  =  "+proj=longlat +datum=WGS84 +ellps=WGS84")

#identify  groups  of  5  clusters  using  the  spatialsample  package
set.seed(99)  #set  seed  to  get  same  split  each  time
clusters  <-  spatial_block_cv(data0_sf, 
                               method = "random", n = 30, #method for how blocks are oriented in space & number of blocks
                               relevant_only = TRUE,  v = 3)  #k-means  clustering  to  identify  cross-validation  folds  (3  is  too  few  to  be  robust  but  using  here  to  save  time)

#for  loop  to  create  a  dataframe  that  assigns  a  fold  number  to  each  data  point
splits_df  <-  c()
for(i  in  1:3){
  new_df  <-  assessment(clusters$splits[[i]])  #extract  points  in  fold  number  i
  new_df$fold  <-  i
  new_df  <-  new_df[  c("row_code", "fold")]
  splits_df  <-  rbind(splits_df, new_df)  #bind  all  points  x  fold  id  together
}

splits_df  <-  st_drop_geometry(splits_df)  #drop  shapefiles

#final  data  -  merge  cluster  id  to  final  dataset  for  analysis
analysis_data  <-  merge(data0, splits_df, by  =  "row_code")

#sanity  check:  check  how  many  data  points  are  in  each  fold
table(analysis_data$fold)


#write  df  to  save  for  later
write.csv(analysis_data, "data/a_chamek_ter_mammals_finalData_Oct22.csv")



