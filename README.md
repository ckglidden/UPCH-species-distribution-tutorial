# UPCH geospatial analysis workshop

## Species distribution modeling tutorial
> This github page includes a tutorial of how to model the distribution of _Bradypus tridactylus_ (the pale-throated sloth) throughout the Amazon Basin using MAPBIOMAS land-use/land-cover data.

1. [Presence-background species distribution models](https://github.com/ckglidden/UPCH-species-distribution-tutorial/edit/main/README.md#1-presence-background-species-distribution-models)

2. [Environmental covariate data - MAPBIOMAS LULC](https://github.com/ckglidden/UPCH-species-distribution-tutorial/edit/main/README.md#2-environmental-covariate-data)

3. [Machine-learning based SDMs](https://github.com/ckglidden/UPCH-species-distribution-tutorial/edit/main/README.md#3-machine-learning-based-sdms)


### 1. Presence-background species distribution models
> For this model we are using _B. tridactylus_ as a focal species and all other terrestrial mammals as background species. The background species helps to  us to understand the difference between the focal species and the average landscape over which mammals are sampled (thus accounting for sampling bias in the occurrence points). Point data was downloaded from GBIF using the ["0_download_gbif_points.R" code](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/R_code/0_download_gbif_points.R). 

![Figure 1. Distribution of points](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/b_tridactylus_sdm_point_distribution.png)
**Figure 1.** Distribution of the focal species (_B. tridactylus_) and background species thinned to 1km grid cells. The number of background points was further reduced by using a background point probability mask and sampling 1.5 * the no. of occurrence points (189 _B. tridactylus_ points; 283 bkg points). The code used to create the probability surface and sample the background points is in the ["0_download_gbif_points.R" code](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/R_code/0_download_gbif_points.R). An example of using a background mask to sample background points can be found in [Moyes et al. 2016](https://parasitesandvectors.biomedcentral.com/articles/10.1186/s13071-016-1527-0).

#### Uploading point data to GEE :mosquito:
>Once you have your occurrence data downloaded, upload it to GEE as a csv. Make sure the dataframe has numerical longitude (x) and latitude (y) coordinates (in decimal degrees) as separate columns. Make sure there is a row identifier to match the points to bind multiple datasets after downloading geospatial data. 

_Step 1._ Download occurrence dataset from the data folder: [b_tridactylus_ter_mammals_amazon_thinned_Oct22.csv](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/data/b_tridactylus_ter_mammals_amazon_thinned_Oct22.csv)

_Step 2._ Upload occurrence dataset as a GEE feature collections.</br>
<img src=https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/GEE_csv_asset.png width="250" height="380"></br>
**Figure 2.** Navigation for uploading csv as a feature collection.

### 2. Environmental covariate data
> Species distribution models model the probability a species occurs in pixel _x_ given the environmental conditions (covariates) in pixel _x_. Here, we will download land-use / land-cover data to use as environmental covariates in our model. 

#### Downloading MAPBIOMAS data :earth_americas: </br>

<img src=https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/mapbiomas_amazonia_example.png width="900" height="550"></br>
**Figure 3.** The MAPBIOMAS Panamazonia platform.

_Step 3._ Explore the [MAPBIOMAS website](https://plataforma.panamazonia.mapbiomas.org/)

_Step 4._ Identify LULC categories / datasets of interest using the [MAPBIOMAS legend](https://s3.amazonaws.com/amazonia.mapbiomas.org/leyenda/C%C3%B3digo_de_la_Leyenda_-_colecci%C3%B3n_3.pdf) (make sure it is the correct legend for the MAPBIOMAS assett)

_Step 5._ Explore [MAPBIOMAS GEE code](https://github.com/mapbiomas-brazil/user-toolkit). Many of the functions you need to work with MAPBIOMAS data has pre-written code, you just need to adapt it to your dataset.

_Step 6._ We will build our species distribution model at a 1km^2 resolution. MAPBIOMAS is available at a finer scale resolution (30m) so it is possible to create the model at a finer spatial resolution. Note: the resolution of all covariates should be equal to the resolution of the coarsest variable (given that the precision of your occurrence point is as or more precise than that resolution). The data we fed into GEE are points. We will create a 1km buffer around each point and then use the remainder of the code to calculate area of each land class per grid cell.


```
//read in csv of points (lat, long, & a row identifier)
var amazon_mammals = ee.FeatureCollection('users/cglidden/b_tridactylus_ter_mammals_amazon_thinned_Oct22'); 

/////functions to create buffer zone around each occurence point
function bufferPoints(radius, bounds) {
  return function(pt) {
    pt = ee.Feature(pt);
    return bounds ? pt.buffer(radius).bounds() : pt.buffer(radius);
  };
}

/////Implement function to make each point is a 1km^2 polygon
var pointBuffers = amazon_mammals.map(bufferPoints(500, true)); //true = square pixel

```

[Caroline's GEE code](https://code.earthengine.google.com/4e4104e7bb0bac0ab475e34d7681295e)  -- **break this into more steps and include code chunks below**
_Step 8._ Calculating area...


_Step 9._ Exporting to a table...


_Step 10._ Skip actually running the GEE code for now and download the ["ter_mammals_lulc_Oct22.csv" file](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/data/ter_mammals_lulc_Oct22.csv.csv) of pre-downloaded data from the data folder.


#### Clean MAPBIOMAS data :broom:

R code for the following section can be found in the ["1_cleaning_data.R"](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/R_code/1_cleaning_data.R)

_Step 11._ Using the data downloaded in step 7 and the code below, we will relabel MAPBIOMAS classes to make it easier to view results. We'll then aggregate LULC data by taking the mean across the study period. Note: Given the pace of LULC change, this is a really coarse way of aggregating the data and we likely loose a lot of signal.

```
#load libraries
library(tidyr); library(dplyr)

#-----------------------------------#
#read in dataset.                   #
#-----------------------------------#

mapbiomas <- read.csv("data/ter_mammals_lulc_Oct22.csv.csv")

#---------------------------------------#
#update label for MAPBIOMAS classes     #
#---------------------------------------#
mapbiomas$class[mapbiomas$class == 3] <- "forest_formation"
mapbiomas$class[mapbiomas$class == 6] <- "flooded_forest"
mapbiomas$class[mapbiomas$class == 11] <- "wetland"
mapbiomas$class[mapbiomas$class == 12] <- "grassland"
mapbiomas$class[mapbiomas$class == 14] <- "farming"
mapbiomas$class[mapbiomas$class == 24] <- "urban"
mapbiomas$class[mapbiomas$class == 25] <- "other_non_vegetated"
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

```

_Step 9._ Using the code below, we will now clean our LULC data a bit more by removing highly colinear variables. While machine learning can handle multicolinearity when making predictions, removing colinear variables can still be helpful for model interpretation. The correlation value depends on your questions and dataset, but we will use a 0.7 correlation cutoff in the code below. We will use a pair-wise analysis but another option is a variable inflation analysis (or you can use both).

```
#load libraries
library(PerformanceAnalytics)

#calculate a correlation matrix
corr <- abs(cor(mapbiomas_mean_wide[2:ncol(mapbiomas_mean_wide)])) 

#plot correlation and distribution of variables
chart.Correlation(mapbiomas_mean_wide[2:ncol(mapbiomas_mean_wide)], 
                  histogram = TRUE, method = "pearson")
                  
 #no variables are super correlated so the last step for this section is to merge "mapbiomas_mean_wide" by the occurrence set
```
<img src= https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/lulc_correlation.png></br>
**Figure 4.** Pearson's correlation coeffecients of variables.

### 3. Machine-learning based SDMs 
> short blurb on different algorithms that have been used and why we are using one specific one for the workshop (RF?)


#### Model tuning & testing using spatial cross-validation :white_check_mark:
_Step 10._  Here we will split our data in 5 folds (5 subsets) for 5-fold cross validation. It is important to test the perfomance of your model using a hold-out test set. This allows you to evaluate if your model is predicting generazliable patterns, or if it only learning the traing data (and thus "overfitting"). One way to test out-of-sample model performance is using k-fold cross validation. K-fold cross validation splits the data into k folds, it then trains and tests the model k times (where, for each iteration, one fold is a hold out fold and the remaning folds (k-1) are used for training the model). K-fold cross validation helps to test model performance across different subsets of data where the subsets are sampled without replacement. For many applications of species distribution modeling, it is ideal to use spatial cross-validation where folds are separated in space so to avoid issues of autocorrelation that arise from test and training points being very close to each other. See the figure 5. for a visual explanation. Here we will use the R package blockCV. Methods for splitting folds can be dependent on your data and study questions. View the [blockCV paper (Valavi et al. 2021)](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13107) and the [_spatialsample_](https://spatialsample.tidymodels.org/) rpackage to learn of different ways to split data. Below we will use k-means clustering because it is quick to implement.</br>

<img src= https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/spatialcv_visualization.png></br>
**Figure 5.** Visualization of spatial partitioning versus random test versus train set. Figure from towards data science ["Spatial CV Using Sickit-learn"](https://towardsdatascience.com/spatial-cross-validation-using-scikit-learn-74cb8ffe0ab9).

```####random forest species distribution models
library(spatialsample); library(sf)

#--------------------------------------------------------------------------------------------#
#read in occ data & covariate data, merge by row id  (or merge from output of above code)    #
#--------------------------------------------------------------------------------------------#

lulc <- read.csv("data/.  .csv")
amazon_basin_pnts <-  read.csv("data/b_tridactylus_ter_mammals_amazon_thinned_Oct22.csv.csv")

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
```

_Step 11._ Now we train the model on each set of k-1 folds and test it on the holdout fold. For each iteration, we tune the randomForest model to optimize model performance. The tuning step can also be used to prevent over-fitting, depending on your dataset and the parameter values you search over. There are different methods for tuning a machine-learning model. Below we use a [hypergrid search](https://afit-r.github.io/random_forests#tune) and select the final parameters based on the combination that yields the best model performance.

```
library(ranger)

#------------------------------------#
#define the grid to search over      #
#------------------------------------#
# hyperparameter grid search - keeping it small to save time

hyper_grid <- expand.grid(
  mtry       = seq(20, 30, by = 5), #the number of variables to randomly sample as candidates at each split
  node_size  = seq(3, 9, by = 3), #minimum number of samples within the terminal nodes
  sampe_size = c(.6, .70, .80), #the number of samples to train on
  num.trees  = c(500, 1000), #number of trees
  OOB_RMSE   = 0
)

#------------------------------------------------------------#
#tune, train model, calculate out-of-sample performance      #
#------------------------------------------------------------#

for(i in 1:5){

#tune model


#train model



#save model performance results



}

#------------------------------------------------------------#
#calculate average out of sample performance                 #
#------------------------------------------------------------#







```

#### Training of final model :woman_technologist:

```


``` 

#### Model interpretation :bar_chart: :chart_with_upwards_trend:

##### Variable importance
> blurb about different variable importance (permutation importance, shapley, etc, etc), choose on for workshop

include variable importance plot?

```
code to quantify variable importance and plot it
```

##### PDPs
> blurb defining pdps

PDP plot?

```
code for generating pdp plots
```

##### Model predictions
> blurb about how to create predictions and interpretation of predictions

distribution map?

```
code for generating distribution map based on geoTIFFs of prediction variables? might need to create it for a small area, otherwise files might be too big for github
```
