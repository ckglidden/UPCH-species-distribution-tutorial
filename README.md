# UPCH geospatial analysis workshop

## Species distribution modeling tutorial
> This github page includes a tutorial of how to model the distribution of _Corythopis torquatus_ (the Ringed antpipit) using MAPBIOMAS land-use/land-cover data.

1. [Presence-background species distribution models](https://github.com/ckglidden/UPCH-species-distribution-tutorial/edit/main/README.md#1-presence-background-species-distribution-models)

2. [Environmental covariate data - MAPBIOMAS LULC](https://github.com/ckglidden/UPCH-species-distribution-tutorial/edit/main/README.md#2-environmental-covariate-data)

3. [Machine-learning based SDMs](https://github.com/ckglidden/UPCH-species-distribution-tutorial/edit/main/README.md#3-machine-learning-based-sdms)


### 1. Presence-background species distribution models
> For this model we are using _C. torquatus_ as a focal species, and all other passerine birds as background species. The background species helps to  us to understand the difference between the focal species and the average landscape over which songbirds are sampled (thus accounting for sampling bias in the occurrence points). Point data was downloaded from GBIF using the ["download_gbif_points.R" code](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/R_code/download_gbif_points.R). 

![Figure 1. Distribution of points](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/c_torquatus_sdm_point_distribution.png)
**Figure 1.** Distribution of the focal species (_C. torquatus_) and background species thinned to 100m grid cells. The number of background points was further reduced by using a background point probability mask and sampling 2 * the no. of occurrence points (65 _C. torquatus_ points; 130 bkg points). The code used to create the probability surface and sample the background points is in the ["download_gbif_points.R" code](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/R_code/download_gbif_points.R). An example of using a background mask to sample background points can be found in [Moyes et al. 2016](https://parasitesandvectors.biomedcentral.com/articles/10.1186/s13071-016-1527-0).

#### Uploading point data to GEE :mosquito:
>Once you have your occurrence data downloaded, you should upload it to GEE so that you can download geospatial covariates per point. Make sure the dataframe has numerical longitude (x) and latitude (y) coordinates (in decimal degrees) as separate columns. I also add a row identifier in case I need to match the points to other datasets or bind multiple datasets after downloading geospatial data. 

_Step 1._ Download occurrence dataset from the data folder: [final_passerine_data_Oct20_2022.csv](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/data/final_passerine_dataset_Oct20_2022.csv)

_Step 2._ Upload occurrence dataset as a GEE feature collections.</br>
<img src=https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/GEE_csv_asset.png width="250" height="380"></br>
**Figure 2.** Navigation for uploading csv as a feature collection.

### 2. Environmental covariate data
> short blurb on how environmental data is used to predict the distribution of focal species

#### Downloading MAPBIOMAS data :earth_americas: </br>

<img src=https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/mapbiomas_amazonia_example.png></br>
**Figure 3.** The MAPBIOMAS Panamazonia platform.

_Step 3._ Explore the [MAPBIOMAS website](https://plataforma.panamazonia.mapbiomas.org/)

_Step 4._ Identify LULC categories / datasets of interest using the [MAPBIOMAS legend](https://s3.amazonaws.com/amazonia.mapbiomas.org/leyenda/C%C3%B3digo_de_la_Leyenda_-_colecci%C3%B3n_3.pdf) (make sure it is the correct legend for the MAPBIOMAS assett)

_Step 5._ Explore [MAPBIOMAS GEE code](https://github.com/mapbiomas-brazil/user-toolkit). Many of the functions you need to work with MAPBIOMAS data has pre-written code, you just need to adapt it to your dataset.

_Step 6._ Navigate to [Caroline's GEE code](https://code.earthengine.google.com/4e4104e7bb0bac0ab475e34d7681295e)  -- we will walk through this code as a group, it will let you download LULC area per year for each point in the feature collection.

_Step 7._ Skip actually running the GEE code for now and download the XXX csv of pre-downloaded data (include link) from the data folder.


#### Clean MAPBIOMAS data :broom:

_Step 8._ Using the data downloaded in step 7 and the code below, we will relabel MAPBIOMAS classes to make it easier to view results and then aggregate LULC data by taking the mean across the study period. Note: Given the pace of LULC change, this is a really coarse way of aggregating the data and we likely loose a lot of signal.

```
#load libraries
library(tidyr); library(dplyr)

#-----------------------------------#
#read in datasets                   #
#-----------------------------------#

occ_data <- read.csv("data/final_passerine_dataset_Oct20_2022.csv")
mapbiomas <- read.csv("data/passerine_lulc_Oct2022.csv")

#-----------------------------------#
#update label MAPBIOMAS classes     #
#-----------------------------------#
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

_Step 9._ Using the code below, we will now clean our LULC data a bit by removing highly colinear variables. While machine learning can handle multicolinearity when making predictions, removing colinear variables can still be helpful for model interpretation to remove highly correlated variables. The correlation value depends on your questions and dataset, but we will use a 0.7 correlation cutoff in the code below. We will use a pair-wise analysis below but another option is a variable inflation analysis (or you can use both).

```
#load libraries
library(PerformanceAnalytics)

#correlation matrix
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

```
code chunk to demonstrate how points are clustered for spatial CV
```

couple sentences about what needs to be tuned for RF, how to account for unbalanced data, and different metrics for testing the data? have a ready dataset with splits...

```
for loop with model tuning & testing for each iteration??
```

#### Training of final model :woman_technologist:

```
code chunk for final model w/ all of the data - should we incorporate uncertainty around variable importance at this point?
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
