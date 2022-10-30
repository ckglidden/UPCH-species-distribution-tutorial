# UPCH geospatial analysis workshop

## Species distribution modeling tutorial
> This github page includes a tutorial of how to model the distribution of _Ateles chamek_ (the Peruvian spider monkey) throughout the Amazon Basin using MAPBIOMAS land-use/land-cover data and CHELSA climate data.

1. [Presence-background species distribution models](https://github.com/ckglidden/UPCH-species-distribution-tutorial/edit/main/README.md#1-presence-background-species-distribution-models)

2. [Environmental covariate data - MAPBIOMAS LULC](https://github.com/ckglidden/UPCH-species-distribution-tutorial/edit/main/README.md#2-environmental-covariate-data)

3. [Machine-learning based SDMs](https://github.com/ckglidden/UPCH-species-distribution-tutorial/edit/main/README.md#3-machine-learning-based-sdms)


&nbsp;  
&nbsp;  


### 1. Presence-background species distribution models
> For this model we are using _Ateles chamek_ (the Peruvian spider moneky) as a focal species and all other terrestrial mammals as background species. The background species helps to  us to understand the difference between the focal species and the average landscape over which mammals are sampled (thus accounting for sampling bias in the occurrence points). Point data was downloaded from GBIF using the ["0_download_gbif_points.R" code](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/R_code/0_download_gbif_points.R). 

![Figure 1. Distribution of points](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/a_chamek_sdm_point_distribution.png)
**Figure 1.** Distribution of points for the focal species (_A. chamek_) (blue) and background species (grey) thinned to 1km grid cells across the Amazon Basin. The number of background points was further reduced by using a background point probability mask and sampling 3 * the no. of occurrence points. The code used to create the probability surface and sample the background points is in the ["0_download_gbif_points.R" code](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/R_code/0_download_gbif_points.R). An example of using a background mask to sample background points can be found in [Moyes et al. 2016](https://parasitesandvectors.biomedcentral.com/articles/10.1186/s13071-016-1527-0).

&nbsp;  

#### Uploading point data to GEE :mosquito:
>Once you have your occurrence data downloaded, upload it to GEE as a csv. Make sure the dataframe has numerical longitude (x) and latitude (y) coordinates (in decimal degrees) as separate columns. Make sure there is a row identifier to match the points to bind multiple datasets after downloading geospatial data. 

_Step 1._ Download occurrence dataset from the data folder: [a_chamek_ter_mammals_amazon_thinned_Oct22.csv](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/data/a_chamek_ter_mammals_amazon_thinned_Oct22.csv)

_Step 2._ Upload occurrence dataset as a GEE feature collections.</br>
<img src=https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/GEE_csv_asset.png width="250" height="380"></br>
**Figure 2.** Navigation for uploading csv as a feature collection.

&nbsp;  
&nbsp;  

### 2. Environmental covariate data
> Species distribution models model the probability a species occurs in pixel _x_ given the environmental conditions (covariates) in pixel _x_. Here, we will download land-use / land-cover and climate data to use as environmental covariates in our model. For the purposes of this tutorial, we will run the initial model with a small number of covariates. 

&nbsp;  

#### Downloading MAPBIOMAS data :earth_americas: </br>

<img src=https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/mapbiomas_amazonia_example.png width="900" height="550"></br>
**Figure 3.** The MAPBIOMAS Panamazonia platform.

&nbsp;  

_Step 3._ Explore the [MAPBIOMAS website](https://plataforma.panamazonia.mapbiomas.org/)

_Step 4._ Identify LULC categories / datasets of interest using the [MAPBIOMAS legend](https://s3.amazonaws.com/amazonia.mapbiomas.org/leyenda/C%C3%B3digo_de_la_Leyenda_-_colecci%C3%B3n_3.pdf) (make sure it is the correct legend for the MAPBIOMAS assett)

_Step 5._ Explore [MAPBIOMAS GEE code](https://github.com/mapbiomas-brazil/user-toolkit). Many of the functions you need to work with MAPBIOMAS data has pre-written code, you just need to adapt it to your dataset.

_Step 6._ We will build our species distribution model at a 1km^2 resolution. MAPBIOMAS is available at a finer scale resolution (30m) so it is possible to create the model at a finer spatial resolution. Note: the resolution of all covariates should be equal to the resolution of the coarsest variable (given that the precision of your occurrence point is as or more precise than that resolution). The data we fed into GEE are points. We will create a 1km buffer around each point (step 6) and then calculate area of each land class per grid cell (step 7-8). A GEE script for steps 6-8 can be found [here](https://code.earthengine.google.com/8fc7e788d5971db28c29447100bdb778).


```
//call csv of points (make sure it has lat, long, & a row identifier)
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
**Figure 4**. Part of the feature collection (1km square buffers around each point) visualized using the _Map.addLayer(fcVis)_ function from the code above. </br>

&nbsp;  

_Step 7._ We will now calculate area per each feature (point + buffer) per year in the study period (2001-2020) in the feature collection. _A. chamek_ typically live in lowland forests but are listed as endangered due to habitat loss. For our initial model we will include lulc variables related to forest cover and farming.


_Step 8._ We will export a table with land class area per year per point as a csv to our GEE folder. You can also export the table as a shape file. 

&nbsp;  

#### Downloading climate data :thermometer: </br>
_Step 9._ Species distributions are usually, at least in part, dictated by interations between land cover and climate. We will use the [CHELSA dataset](https://chelsa-climate.org/) to add climatologies to our model. The CHELSA data is at a 1km^2 resolution, which is why we scaled the MAPBIOMAS data up to 1km. We will use the same feature collection from steps 6-8. A GEE script can be found [here](https://code.earthengine.google.com/14b1a32976d3097b5eca6be97cf84559).

```
/////CHELSA climatologies
var bio13_precip_wettest_month = ee.Image('users/cglidden/CHELSA_bio13_1981-2010_V21').select('b1')
            .resample("bilinear")
            .reproject("EPSG:4326", null, 1000) //make sure variable is at 1km res
            .reduce(ee.Reducer.mean()) //get mean climatology per 1km pixel
            .rename('bio13_precip_wettest_month'); //rename band - this will be the column name in your csv

var cmi_min = ee.Image('users/cglidden/CHELSA_1981_2010/CHELSA_cmi_min_1981-2010_V21').select('b1')
            .resample("bilinear")
            .reproject("EPSG:4326", null, 1000)
            .reduce(ee.Reducer.mean())
            .rename('cmi_min');

///compile image to export both variables in one csv
var final_image = bio13_precip_wettest_month.addBands([cmi_min])

  
////export table to drive
Export.table.toDrive({
  collection: final_image.reduceRegions({
    collection: pointBuffers,
    reducer: ee.Reducer.mean(), //now get mean climatology per pointBuffer and export to table
    scale: 1000, //1km resolution
    tileScale: 16
    }),
  description: 'a_chamek_ter_mammals_climate_Oct2022',
  folder: 'GEEexports',
  fileFormat: 'csv',
});
```

&nbsp;  

_Step 10._ Skip actually running the GEE code for now and download the ["a_chamek_ter_mammals_lulc_Oct22.csv" file](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/data/a_chamek_ter_mammals_lulc_Oct22.csv) and th ["a_chamek_ter_mammals_climate_Oct22.csv" file](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/data/a_chamek_ter_mammals_climate_Oct22.csv) pre-downloaded data from the data folder.

&nbsp;  

#### Clean covariate data :broom:

R code for the following sections can be found in the ["1_cleaning_data.R"](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/R_code/1_cleaning_data.R)

&nbsp;  

_Step 11._ Using the data downloaded in step 10 and the code below, we will relabel MAPBIOMAS classes to make it easier to view results. We'll then aggregate LULC data by taking the mean LULC from 2001-2020. Since _A. chamek_ is affected by deforestation we will also look at difference in forest cover between 2001-2020. We will then merge the lulc data with the occurrence data and climate data. The climate data is exported from GEE in an immediately useable form. Note: Given the pace of LULC change, this is a really coarse way of aggregating the data and we likely loose a lot of signal.  

```
#load libraries
library(tidyr); library(dplyr)

#-----------------------------------#
#read in datasets                   #
#-----------------------------------#

mapbiomas <- read.csv("data/a_chamek_ter_mammals_lulc_Oct2022.csv")

climate <- read.csv("data/a_chamek_ter_mammals_climate_Oct2022.csv")
#subset climate data to variables of interest (not necessary but makes the final dataset cleaner
climate <- climate[  c("row_code", "bio13_precip_wettest_month", "cmi_min")]

amazon_basin_pnts <-  read.csv("data/a_chamek_ter_mammals_amazon_thinned_Oct22.csv")

#---------------------------------------#
#update label for MAPBIOMAS classes     #
#---------------------------------------#

#you can determine the classes included in the data using:
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
                            
                            #difference over study period
                            diff_forest_formation = forest_formation[match(2020, year)] - forest_formation[match(2001, year)])

#save dataframe
write.csv(mapbiomas_mean_diff, "data/a_chamek_ter_mammals_lulc_cleaned_Oct2022.csv")

#--------------------------------------------------------------------------------------------------------#
#merge data frames with response (presence / absence labels), climate covariates, and lulc covariates    #
#--------------------------------------------------------------------------------------------------------#

covariates <- left_join(climate, mapbiomas_mean_diff, by = "row_code")
data0 <- left_join(amazon_basin_pnts, covariates, by = "row_code")

```

&nbsp;  

_Step 9._ Using the code below, we will now clean our covariate data a bit more by removing highly colinear variables. While machine learning can handle multicolinearity when making predictions, removing colinear variables can still be helpful for model interpretation. The correlation value depends on your questions and dataset but we will use a 0.7 correlation cutoff in the code below. We will use a pair-wise analysis but another option is a variable inflation analysis (or you can use both).

```
#load libraries
library(PerformanceAnalytics)

# correlation in absolute terms
corr <- abs(cor(data0[7:ncol(data0)])) 

#correlation plot with points on lower left of matrix, correlation coeffecients on upper right, and distributions on the diagonal
chart.Correlation(data0[7:ncol(data0)], 
                  histogram = TRUE, method = "pearson")

#nothing with correlation > 0.7 so we keep all variabes for the analysis
                 
```
<img src= https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/lulc_correlation.png></br>
**Figure 5.** Point plot of pairwise realtionships, covariate distribution, and Pearson's correlation coeffecients of variables.

&nbsp;  

### Challenge Questions
1. Try using GEE to add additional climate or LULC covariates to the final analysis dataset.

&nbsp;  
&nbsp;  

### 3. Machine-learning based SDMs 
> short blurb on different algorithms that have been used and why we are using one specific one for the workshop (RF?)

&nbsp;  

#### Model tuning & testing using spatial cross-validation :white_check_mark:

_Step 10._  Next we will split our data in 3 folds (3 subsets) for 3-fold cross validation. It is important to test the perfomance of your model using a hold-out test set. This allows you to evaluate if your model is predicting generazliable patterns, or if it only learning the traing data (and thus "overfitting"). One way to test out-of-sample model performance is using k-fold cross validation. K-fold cross validation splits the data into k folds, it then trains and tests the model k times (where, for each iteration, one fold is a hold out fold and the remaning folds are used for training the model). K-fold cross validation helps to test model performance across different subsets of data where the subsets are sampled without replacement. For many applications of species distribution modeling, it is ideal to use spatial cross-validation where folds are separated in space so to avoid issues of autocorrelation that arise from test and training points being very close to each other. See the figure 5. for a visual explanation. Here we will use the R package _spatialsample_. Methods for splitting folds can be dependent on your data and study questions. View the [blockCV paper (Valavi et al. 2021)](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13107) and the [_spatialsample_](https://spatialsample.tidymodels.org/) rpackage to learn of different ways to split data. Below we will use block clustering because it is quick to implement.</br>

<img src= https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/spatialcv_visualization.png></br>
**Figure 5.** Visualization of spatial partitioning versus random test versus train set. Figure from towards data science ["Spatial CV Using Sickit-learn"](https://towardsdatascience.com/spatial-cross-validation-using-scikit-learn-74cb8ffe0ab9).

```
library(spatialsample); library(sf)

#--------------------------------------------#
#get  fold  id  by  block  clustering      #
#--------------------------------------------#

#convert  analysis_data  to  a  spatial  object
data0_sf  <-  st_as_sf(x  =  data0,
                       coords  =  c("lon", "lat"),
                       crs  =  "+proj=longlat +datum=WGS84 +ellps=WGS84")

#identify  groups  of  3  clusters  using  the  spatialsample  package
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

splits_df  <-  st_drop_geometry(splits_df)  #drop  geometry

#final  data  -  merge  cluster  id  to  final  dataset  for  analysis
analysis_data  <-  merge(data0, splits_df, by  =  "row_code")

#sanity  check:  check  how  many  data  points  are  in  each  fold
table(analysis_data$fold)


#write  df  to  save  for  later
write.csv(analysis_data, "data/a_chamek_ter_mammals_finalData_Oct22.csv")

```

&nbsp;  

_Step 11._ Now we train the model on each set of folds and test it on the holdout fold. For each iteration, we tune the randomForest model to optimize model performance. The tuning step can also be used to prevent over-fitting, depending on your dataset and the parameter values you search over. There are different methods for tuning a machine-learning model. Below we use a [hypergrid search](https://afit-r.github.io/random_forests#tune) and select the final parameters based on the combination that yields the best model performance. If you had trouble running _Step 10._ you can download the ["a_chamek_ter_mammals_finalData_Oct22.csv"](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/data/a_chamek_ter_mammals_finalData_Oct22.csv) to use for the next few steps. R code for the following sections can be found in the ["2_random_forest_sdm.R"](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/R_code/2_random_forest_sdm.R) </br>

```
library(ranger)

#------------------------------------------------#
#reduce dataset to columns of interest           # this is stylistic as you could also specify the covariates in the forumal call below
#------------------------------------------------#
analysis_data_v2  <-  data0[  c("presence", "fold", "bio13_precip_wettest_month", "cmi_min",
                                        "mean_forest",  "mean_farming", "diff_forest_formation")]

#for many ML algorithms  you should make sure your response variables are not grouped in the data-frame (e.g.  all 1s and 0s next to each other)
analysis_data_v2 <- analysis_data_v2[sample(1:nrow(analysis_data_v2)), ]


#------------------------------------#
#tune, train, model                  #
#------------------------------------#

#create  empty  dataframe  to  for  loop  to  store  results    one  row  for  each  fold
rf_performance  <-  data.frame(model  =  rep("RF", 3),
                               fold_id  =  1:3,
                               auc  =  rep(NA, 3),
                               sensitivity  =  rep(NA, 3),
                               specificity  =  rep(NA, 3),
                               oob_error  =  rep(NA, 3),
                               presence  =  rep(NA, 3),    #number  of  presence  points  in  the  fold
                               background  =  rep(NA, 3))  #number  of  bkg  points  in  the  fold

#create  empty  dataframe  to  store  parameters  used  to  train  each  model
hypergrid_final  <-  data.frame(mtry  =  rep(NA,3),    #the  number  of  variables  to  randomly  sample  as  candidates  at  each  split
                                                            node_size    =  rep(NA, 3),    #minimum  number  of  samples  within  the  terminal  nodes
                                                            sampe_size  =  rep(NA, 3))    #the  number  of  samples  to  train  on

    
    
for(i  in  1:3){  #  run  one  iteration  per  fold
    
    train  <-  analysis_data_v2[analysis_data_v2$fold !=  i, ];  train  <-  train[-2]
    test  <-  analysis_data_v2[analysis_data_v2$fold ==  i, ];  test  <-  test[-2]
    
    #remove  any  rows  with  NAs  bc  RF  can't  handle  missing  data
    train_complete  <-  train[complete.cases(train), ]
    test_complete  <-  test[complete.cases(test), ]
    
    
    #-----------------------------------------------#
    #define  the  grid  to  search  over            #
    #-----------------------------------------------#
    #  the  function  below  creates  a  grid  with  all  combinations  of  parameters
    
    hyper_grid  <-  expand.grid(
        mtry =  seq(1, 3, by  =  1),    #the  number  of  variables  to  randomly  sample  as  candidates  at  each  split
        node_size =  seq(1,4, by  =  1),    #shallow trees
        sampe_size  =  c(.6, .70, .80),    #the  number  of  samples  to  train  on
        OOB_RMSE =  0
    )
    
    #tune  model
    for(j  in  1:nrow(hyper_grid)){
        
        #  train  model
        model  <-  ranger(
            formula  =  presence  ~  .,    
            data  =  train_complete,    
            num.trees  =  2000,  
            mtry  =  hyper_grid$mtry[j],
            min.node.size  =  hyper_grid$node_size[j],  
            sample.fraction  =  hyper_grid$sampe_size[j], 
            probability  =  TRUE, 
            replace = TRUE,
            splitrule = "hellinger",
            seed  =  123
            )
        
        #  add  OOB  error  to  grid
        hyper_grid$OOB_RMSE[i]  <-  sqrt(model$prediction.error)
    }
    
    #arrange  the  hypergrid  so  the  lowest  out-of-bag  error  (best  performing  set  of  parameters)  is  in  the  first  row
    hyper_grid2  <-  hyper_grid  %>%  
        dplyr::arrange(OOB_RMSE)
    
    #train  model
    train_model  <-  ranger(
        formula  =  presence  ~  .,    
        data  =  train_complete,   
        #use  the  first  row  of  the  grid  as  model  parameters
        num.trees  =  2000,  
        mtry  =  hyper_grid2$mtry[1], 
        min.node.size  =  hyper_grid2$node_size[1], 
        sample.fraction  =  hyper_grid2$sampe_size[1],
        probability  =  TRUE, 
        replace = TRUE,
        splitrule = "hellinger",
        seed  =  123)
    
    #save  model  performance  results
    pred0  <-  predict(train_model, data=test_complete);  pred  <-  pred0$predictions[,1]
    auc  <-  pROC::roc(response=test_complete[,"presence"], predictor=pred, levels=c(0, 1), auc  =  TRUE)
    rf_performance[i, "auc"]  <-  auc$auc
    best.threshold  <-  pROC::coords(auc, "best", ret  =  "threshold")
    metrica.format  <-  data.frame(cbind(ifelse(test_complete[,"presence"]==1, 1, 0)),  ifelse(pred  >=  best.threshold[1, 1], 1, 0));  colnames(metrica.format)  <-  c("labels", "predictions");  rownames(metrica.format)  <-  1:dim(metrica.format)[1]
    sensitivity  <-  metrica::recall(data  =  metrica.format, obs  =  labels, pred  =  predictions)$recall  
    rf_performance[i, "sensitivity"]  <-  sensitivity  
    specificity  <-  metrica::specificity(data  =  metrica.format, obs  =  labels, pred  =  predictions)$spec
    rf_performance[i, "specificity"]  <-  specificity
    rf_performance[i, "oob_error"]  <-  train_model$prediction.error
    rf_performance[i, "presence"]  <-  nrow(subset(test, presence  ==  1))
    rf_performance[i, "background"]  <-  nrow(subset(test, presence  ==  0))
    
    #save hypergrid results to use for final model
    hypergrid_final[i, "mtry"]  <-  hyper_grid2$mtry[1]
    hypergrid_final[i, "node_size"]  <-  hyper_grid2$node_size[1]
    hypergrid_final[i, "sampe_size"]  <-  hyper_grid2$sampe_size[1]

}

```
**include model performance table**

&nbsp;  

#### Training of final model :woman_technologist:
_Step 12._ Now train the final model that will be used for interpretation and prediction. Below, I used the average hyperpameters from the tuning steps.

```
#------------------------------------------------------------#
#train  final  model                                                                                      #
#------------------------------------------------------------#

final_model  <-  ranger(
    formula  =  presence  ~  .,    
    data  =  analysis_data_v2[complete.cases(analysis_data_v2), -2],    #complete case dataset without fold column
    #parameters used here are the averages from hypergrid_final
    num.trees  =  2000 ,
    mtry  =  2,  
    min.node.size  =  1,  
    sample.fraction  =  0.6,
    probability  =  TRUE,  
    replace = TRUE,
    splitrule = "hellinger",
    importance  =  'permutation',    #specify  this  to  get  variable  importance  in  the  next  step
    seed  =  123)

``` 

&nbsp;  

#### Model interpretation :bar_chart: :chart_with_upwards_trend:

##### Variable importance
_Step 12._ There are many ways to calculate variable importance. Here, we will use a intuitive and model agnostic measure of variable importance. In sum, we will calculate change in model performance when a focal variable is randomly permuted, which will tell us the degree to which the variable contributes to model performance and thus accuracy of model predictions. **need to define y axis better**


```
#------------------------------------------------------------#
#variable  importance                                                                                  #
#------------------------------------------------------------#

#extract  model  results  to  get  permutation  importance
permutation_importance  <-  data.frame(variable  =  rownames(as.data.frame(final_model$variable.importance)),  
                                                                          importance  =  as.vector(final_model$variable.importance))
    
#plot  importance
ggplot(permutation_importance, aes(x  =  variable, y  =  importance))  +
    geom_bar(stat="identity")  +
    ggtitle("permutation  importance")  +
    coord_flip()  +
    theme_classic()


```

<img src= https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/variable_importance_plot.png>

**Figure x** variable importance

&nbsp;  

##### Partial dependence plots
Partial dependence plots (PDPs) depict the relationship between the probability of species occurrence and the variable of interest across a range of values for that variable. At each value of the variable, the model is evaluated for all values of the other covariates. The final model output is the average predicted probability across all model inputs. ![image](https://user-images.githubusercontent.com/44818401/198905725-48ddf02b-2b0a-4432-a1f3-5501c97631a9.png)



```
#------------------------------------------------------------#
#pdps                                                                                                                #
#------------------------------------------------------------#

#list of covariate names to generate pdps for and loop through
var_names  <-  names(analysis_data_v2[complete.cases(analysis_data_v2),  -c(1, 2)]) #analysis data set exlucing 'presence' and 'fold' column

#dataframe  to  make  partial  dependence  plots
pd_df  =  data.frame(matrix(vector(), 0, 3, dimnames=list(c(), c('variable', 'value', 'yhat'))),  
                                      row.names  =  NULL, stringsAsFactors=F)

#loop  through  each  variable
for  (j  in  1:length(var_names))  { 
    
    output  <-  as.data.frame(pdp::partial(final_model, pred.var  =  var_names[j], prob  =  TRUE, train  =  analysis_data_v2[complete.cases(analysis_data_v2), -2]))
    
    loop_df  <-  data.frame(variable = rep(var_names[j], length(output[[1]])),
                            value  =  output[[1]],
                            yhat  =  output[[2]])
    
    pd_df  <-  rbind(pd_df, loop_df)  
}


#plot  pdps  for  each  variable
ggplot(pd_df, aes(x  =  value, y=  yhat))  +
    geom_smooth()  +
    ylab("probability")  +
    facet_wrap(~variable, scales  =  "free")  +
    theme_bw()
```
<img src= https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/pdp_plot.png>

**Figure x** pdps...

&nbsp;  

##### Model predictions
> blurb about how to create predictions and interpretation of predictions

distribution map?

```
code for generating distribution map based on geoTIFFs of prediction variables? might need to create it for a small area, otherwise files might be too big for github
```

&nbsp; 

### Challenge questions
1. How would you add uncertainty to variable importance, functional relationships (pdps), or prediction maps?
2. How does you model perfromance and interpretation change when you add / remove covariates?
3. Pick a new focal species to create an SDM. You will have to use the "notThinned" datasets to select a focal species and create background points.
