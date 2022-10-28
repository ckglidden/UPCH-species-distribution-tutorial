# UPCH geospatial analysis workshop

## Species distribution modeling tutorial
> This github page includes a tutorial of how to model the distribution of _Marmosops noctivagus_ (the white-bellied slender opposum) throughout the Amazon Basin using MAPBIOMAS land-use/land-cover data.

1. [Presence-background species distribution models](https://github.com/ckglidden/UPCH-species-distribution-tutorial/edit/main/README.md#1-presence-background-species-distribution-models)

2. [Environmental covariate data - MAPBIOMAS LULC](https://github.com/ckglidden/UPCH-species-distribution-tutorial/edit/main/README.md#2-environmental-covariate-data)

3. [Machine-learning based SDMs](https://github.com/ckglidden/UPCH-species-distribution-tutorial/edit/main/README.md#3-machine-learning-based-sdms)


### 1. Presence-background species distribution models
> For this model we are using _Ateles chamek_ (the Peruvian spider moneky) as a focal species and all other terrestrial mammals as background species. The background species helps to  us to understand the difference between the focal species and the average landscape over which mammals are sampled (thus accounting for sampling bias in the occurrence points). Point data was downloaded from GBIF using the ["0_download_gbif_points.R" code](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/R_code/0_download_gbif_points.R). 

![Figure 1. Distribution of points](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/a_chamek_sdm_point_distribution.png)
**Figure 1.** Distribution of points for the focal species (_A. chamek_) (blue) and background species (grey) thinned to 1km grid cells across the Amazon Basin. The number of background points was further reduced by using a background point probability mask and sampling 3 * the no. of occurrence points. The code used to create the probability surface and sample the background points is in the ["0_download_gbif_points.R" code](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/R_code/0_download_gbif_points.R). An example of using a background mask to sample background points can be found in [Moyes et al. 2016](https://parasitesandvectors.biomedcentral.com/articles/10.1186/s13071-016-1527-0).

#### Uploading point data to GEE :mosquito:
>Once you have your occurrence data downloaded, upload it to GEE as a csv. Make sure the dataframe has numerical longitude (x) and latitude (y) coordinates (in decimal degrees) as separate columns. Make sure there is a row identifier to match the points to bind multiple datasets after downloading geospatial data. 

_Step 1._ Download occurrence dataset from the data folder: [a_chamek_ter_mammals_amazon_thinned_Oct22.csv](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/data/a_chamek_ter_mammals_amazon_thinned_Oct22.csv)

_Step 2._ Upload occurrence dataset as a GEE feature collections.</br>
<img src=https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/GEE_csv_asset.png width="250" height="380"></br>
**Figure 2.** Navigation for uploading csv as a feature collection.

### 2. Environmental covariate data
> Species distribution models model the probability a species occurs in pixel _x_ given the environmental conditions (covariates) in pixel _x_. Here, we will download land-use / land-cover and climate data to use as environmental covariates in our model. 

#### Downloading MAPBIOMAS data :earth_americas: </br>

<img src=https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/mapbiomas_amazonia_example.png width="900" height="550"></br>
**Figure 3.** The MAPBIOMAS Panamazonia platform.

_Step 3._ Explore the [MAPBIOMAS website](https://plataforma.panamazonia.mapbiomas.org/)

_Step 4._ Identify LULC categories / datasets of interest using the [MAPBIOMAS legend](https://s3.amazonaws.com/amazonia.mapbiomas.org/leyenda/C%C3%B3digo_de_la_Leyenda_-_colecci%C3%B3n_3.pdf) (make sure it is the correct legend for the MAPBIOMAS assett)

_Step 5._ Explore [MAPBIOMAS GEE code](https://github.com/mapbiomas-brazil/user-toolkit). Many of the functions you need to work with MAPBIOMAS data has pre-written code, you just need to adapt it to your dataset.

_Step 6._ We will build our species distribution model at a 1km^2 resolution. MAPBIOMAS is available at a finer scale resolution (30m) so it is possible to create the model at a finer spatial resolution. Note: the resolution of all covariates should be equal to the resolution of the coarsest variable (given that the precision of your occurrence point is as or more precise than that resolution). The data we fed into GEE are points. We will create a 1km buffer around each point and then use the remainder of the code to calculate area of each land class per grid cell.


```
//read in csv of points (make sure it has lat, long, & a row identifier)
var amazon_mammals = ee.FeatureCollection('users/cglidden/a_chamek_ter_mammals_amazon_thinned_Oct22'); 

/////functions to create buffer zone around each occurence point
function bufferPoints(radius, bounds) {
  return function(pt) {
    pt = ee.Feature(pt);
    return bounds ? pt.buffer(radius).bounds() : pt.buffer(radius);
  };
}

/////Implement function to make each point is a 1km^2 polygon
var pointBuffers = amazon_mammals.map(bufferPoints(500, true)); //true = square pixel


// Paint FeatureCollection to GEE map for visualization.
var fcVis = pointBuffers.draw({color: '800080', pointRadius: 10});
Map.setCenter(-69.60, -12.39, 10) //coordinates & degree to zoom in
Map.addLayer(fcVis);

```
<img src=https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/gee_feature_collection.png></br>
**Figure 4**. The feature collection (1km square buffers around each point) visualized using the _Map.addLayer(fcVis)_ function from the code above. </br>
_Step 8._ Calculating area...


_Step 9._ Exporting to a table...


_Step 10._ Skip actually running the GEE code for now and download the ["m_noctivagus_ter_mammals_lulc_Oct22.csv" file](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/data/a_chamek_ter_mammals_lulc_Oct22.csv.csv) of pre-downloaded data from the data folder.


#### Clean MAPBIOMAS data :broom:

R code for the following section can be found in the ["1_cleaning_data.R"](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/R_code/1_cleaning_data.R)

_Step 11._ Using the data downloaded in step 10 and the code below, we will relabel MAPBIOMAS classes to make it easier to view results. We'll then aggregate LULC data by taking the mean LULC from 2010-2020. Note: Given the pace of LULC change, this is a really coarse way of aggregating the data and we likely loose a lot of signal.  

```
#load libraries
library(tidyr); library(dplyr)

#-----------------------------------#
#read in datasets                   #
#-----------------------------------#

occ_data <- read.csv("data/m_noctivagus_ter_mammals_amazon_thinned_Oct22.csv")
mapbiomas <- read.csv("data/m_noctivagus_ter_mammals_lulc_Oct2022.csv")

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
mapbiomas$class[mapbiomas$class == 30] <- "mining"
mapbiomas$class[mapbiomas$class == 33] <- "river_lake_ocean"
mapbiomas$class[mapbiomas$class == 34] <- "glacier"

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

_Step 9._ Using the code below, we will now clean our LULC data a bit more by removing highly colinear variables. While machine learning can handle multicolinearity when making predictions, removing colinear variables can still be helpful for model interpretation. The correlation value depends on your questions and dataset, but we will use a 0.75 correlation cutoff in the code below. We will use a pair-wise analysis but another option is a variable inflation analysis (or you can use both).

```
#load libraries
library(PerformanceAnalytics)

#----------------------------------------------------------#
#pearsons correlation analysis                             #
#----------------------------------------------------------#

# correlation in absolute terms
corr <- abs(cor(mapbiomas_mean_wide[2:ncol(mapbiomas_mean_wide)])) 

#correlation plot with points on lower left of matrix, correlation coeffecients on upper right, and distributions on the diagonal
chart.Correlation(mapbiomas_mean_wide[2:ncol(mapbiomas_mean_wide)], 
                  histogram = TRUE, method = "pearson")

#no variables are super correlated so we will keep them for now and write the cleaned datatable to a csv for later use
write.csv(mapbiomas_mean_wide, "data/m_noctivagus_ter_mammals_lulc_cleaned_Oct2022.csv")

                 
```
<img src= https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/lulc_correlation.png></br>
**Figure 5.** Pearson's correlation coeffecients of variables.

### 3. Machine-learning based SDMs 
> short blurb on different algorithms that have been used and why we are using one specific one for the workshop (RF?)


#### Model tuning & testing using spatial cross-validation :white_check_mark:
_Step 10._  Here we will split our data in 5 folds (5 subsets) for 5-fold cross validation. It is important to test the perfomance of your model using a hold-out test set. This allows you to evaluate if your model is predicting generazliable patterns, or if it only learning the traing data (and thus "overfitting"). One way to test out-of-sample model performance is using k-fold cross validation. K-fold cross validation splits the data into k folds, it then trains and tests the model k times (where, for each iteration, one fold is a hold out fold and the remaning folds (k-1) are used for training the model). K-fold cross validation helps to test model performance across different subsets of data where the subsets are sampled without replacement. For many applications of species distribution modeling, it is ideal to use spatial cross-validation where folds are separated in space so to avoid issues of autocorrelation that arise from test and training points being very close to each other. See the figure 5. for a visual explanation. Here we will use the R package _spatialsample_. Methods for splitting folds can be dependent on your data and study questions. View the [blockCV paper (Valavi et al. 2021)](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13107) and the [_spatialsample_](https://spatialsample.tidymodels.org/) rpackage to learn of different ways to split data. Below we will use block clustering because it is quick to implement.</br>

<img src= https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/spatialcv_visualization.png></br>
**Figure 5.** Visualization of spatial partitioning versus random test versus train set. Figure from towards data science ["Spatial CV Using Sickit-learn"](https://towardsdatascience.com/spatial-cross-validation-using-scikit-learn-74cb8ffe0ab9).

```####random forest species distribution models
library(spatialsample); library(sf)

#--------------------------------------------------------------------------------------------#
#read in occ data & covariate data, merge by row id  (or merge from output of above code)    #
#--------------------------------------------------------------------------------------------#

lulc <- read.csv("data/m_noctivagus_ter_mammals_lulc_cleaned_Oct2022.csv")
amazon_basin_pnts <-  read.csv("data/m_noctivagus_ter_mammals_amazon_thinned_Oct22.csv")

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

splits_df <- st_drop_geometry(splits_df) #drop geometry (alternatively, you could keep geometry and save this as a shapefile

#final data - merge cluster id to final dataset for analysis
analysis_data <- merge(data0, splits_df, by = "row_code")

#sanity check: check how many data points are in each fold
table(analysis_data$fold)
```

_Step 11._ Now we train the model on each set of k-1 folds and test it on the holdout fold. For each iteration, we tune the randomForest model to optimize model performance. The tuning step can also be used to prevent over-fitting, depending on your dataset and the parameter values you search over. There are different methods for tuning a machine-learning model. Below we use a [hypergrid search](https://afit-r.github.io/random_forests#tune) and select the final parameters based on the combination that yields the best model performance. If you had trouble running _Step 10._ you can download the ["m_noctivagus_ter_mammals_finalData_Oct22.csv"](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/data/b_tridactylus_ter_mammals_finalData_Oct22.csv) to use for the next few steps.

```
library(ranger)

#------------------------------------------------#
#run spatial cv to evaluate model performance    #
#------------------------------------------------#



#------------------------------------#
#tune, train model,                  #
#------------------------------------#

#create empty dataframe to for loop to store results, one row for each fold
rf_performance <- data.frame(model = rep("RF", 5),
                             fold_id = 1:5,
                             auc = rep(NA, 5),
                             sensitivity = rep(NA, 5),
                             specificity = rep(NA, 5),
                             oob_error = rep(NA, 5),
                             presence = rep(NA, 5), #number of presence points in the fold
                             background = rep(NA, 5)) #number of bkg points in the fold

for(i in 1:5){
  
  train <- analysis_data[analysis_data$fold != i, ]
  test <- analysis_data[analysis_data$fold == i, ]
  
  #remove any rows with NAs bc RF can't handle missing data
  train_complete <- train[complete.cases(train),]
  test_complete <- test[complete.cases(test),]
  
  
  #------------------------------------#
  #define the grid to search over      #
  #------------------------------------#
  # hyperparameter grid search - keeping it small to save time
  # the function below creates a grid with all combinations of parameters
  
  hyper_grid <- expand.grid(
    mtry       = seq(1, 4, by = 1), #the number of variables to randomly sample as candidates at each split
    node_size  = seq(3, 9, by = 3), #minimum number of samples within the terminal nodes
    sampe_size = c(.6, .70, .80), #the number of samples to train on
    num.trees  = c(500, 1000), #number of trees
    OOB_RMSE   = 0
  )
  
  #tune model
  for(j in 1:nrow(hyper_grid)) {
    
    # train model
    model <- ranger(
      formula = presence ~ farming + urban + flooded_forest + forest_formation + river_lake_ocean, 
      data = train_complete, 
      num.trees = hyper_grid$num.trees[j],
      mtry = hyper_grid$mtry[j],
      min.node.size = hyper_grid$node_size[j],
      sample.fraction = hyper_grid$sampe_size[j],
      classification = TRUE,
      seed = 123
    )
    
    # add OOB error to grid
    hyper_grid$OOB_RMSE[i] <- sqrt(model$prediction.error)
  }
  
  #arrange the hypergrid so the lowest out-of-bag error (best performing set of parameters) is in the first row
  hyper_grid2 <- hyper_grid %>% 
    dplyr::arrange(OOB_RMSE)
  
  #train model
  train_model <- ranger(
    formula = presence ~ farming + urban + flooded_forest + forest_formation + river_lake_ocean, 
    data = train_complete, 
    num.trees = hyper_grid2$num.trees[1],
    mtry = hyper_grid2$mtry[1],
    min.node.size = hyper_grid2$node_size[1],
    sample.fraction = hyper_grid2$sampe_size[1],
    classification = TRUE,
    seed = 123)
  
  #save model performance results
  pred0 <- predict(train_model, data=test_complete); pred <- pred0$predictions
  auc <- pROC::roc(response=test_complete[,"presence"], predictor=pred, levels=c(0,1), auc = TRUE)
  rf_performance[i, "auc"] <- auc$auc
  best.threshold <- pROC::coords(auc, "best", ret = "threshold")
  metrica.format <- data.frame(cbind(ifelse(test_complete[,"presence"]==1,1,0)),ifelse(pred >= best.threshold[1,1],1,0)); colnames(metrica.format) <- c("labels","predictions"); rownames(metrica.format) <- 1:dim(metrica.format)[1]
  sensitivity <- metrica::recall(data = metrica.format, obs = labels, pred = predictions)$recall 
  rf_performance[i, "sensitivity"] <- sensitivity 
  specificity <- metrica::specificity(data = metrica.format, obs = labels, pred = predictions)$spec
  rf_performance[i, "specificity"] <- specificity
  rf_performance[i, "oob_error"] <- train_model$prediction.error
  rf_performance[i, "presence"] <- nrow(subset(test, presence == 1))
  rf_performance[i, "background"] <- nrow(subset(test, presence == 0))
  
  
}

#------------------------------------------------------------#
#calculate average out of sample performance                 #
#------------------------------------------------------------#

model_performance <- data.frame(metric = names(rf_performance)[2:ncol(rf_performance)],
                                mean_metric = colMeans(rf_performance[2:ncol(rf_performance)]))




```

#### Training of final model :woman_technologist:

```
#------------------------------------------------------------#
#train final model                                           #
#------------------------------------------------------------#

final_model <- ranger(
  formula = presence ~ farming + urban + flooded_forest + forest_formation + river_lake_ocean, 
  data = analysis_data[complete.cases(analysis_data),], 
  num.trees = ...,
  mtry = ...,
  min.node.size = ...,
  sample.fraction =...,
  classification = TRUE,
  importance = 'permutation', #specify this to get variable importance in the next step
  seed = 123)

``` 

#### Model interpretation :bar_chart: :chart_with_upwards_trend:

##### Variable importance
_Step 12._ There are many ways to calculate variable importance. Here, we will use a intuitive and model agnostic measure of variable importance. In sum, we will calculate change in model performance when a focal variable is randomly permuted, which will tell us the degree to which the variable contributes to model performance and thus accuracy of model predictions. **need to define y axis better**


```
#------------------------------------------------------------#
#variable importance                                         #
#------------------------------------------------------------#

#extract model results to get permutation importance
permutation_importance <- data.frame(variable = rownames(as.data.frame(final_model$variable.importance)),
                                     importance = as.vector(final_model$variable.importance)) 
  
#plot importance
ggplot(permutation_importance, aes(x = variable, y = importance)) +
  geom_bar(stat="identity") +
  ggtitle("permutation importance") +
  theme_classic()

```

include variable importance plot?

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

### Challenge questions
1. How would you add uncertainty to variable importance, functional relationships (pdps), or prediction maps?
2. Try using GEE to add additional climate or LULC covariates. How does you model perfromance and interpretation change?
3. Pick a new focal species to create an SDM. You will have to use the "notThinned" datasets to select a focal species and create background points.
