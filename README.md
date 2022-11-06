# UPCH geospatial analysis workshop

## Species distribution modeling tutorial
> This github page includes a tutorial of how to model the distribution of _Ateles chamek_ (the Peruvian spider monkey) throughout the Amazon Basin using MAPBIOMAS land-use/land-cover data and CHELSA climate data.

1. Species presence and background data

2. Environmental data - MAPBIOMAS LULC & CHELSA Climate

3. Machine-learning based SDMs - Random Forest


## 1. Species presence and background data :monkey:
> For this model we are using _Ateles chamek_ (the Peruvian spider moneky) as a focal species and all other terrestrial mammals as background species. The background species help us to understand the difference between the focal species and the average landscape over which mammals are sampled (thus accounting for sampling bias in the occurrence points). Point data was downloaded from GBIF using the ["0_download_gbif_points.R" code](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/R_code/0_download_gbif_points.R). To save time we will jump in using the already downloaded dataset as referenced below, but the code will be here for your future reference. 

![Figure 1. Distribution of points](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/a_chamek_sdm_point_distribution.png)
**Figure 1.** Distribution of points for the focal species (_A. chamek_) (blue) and background species (grey) thinned to 1km grid cells across the Amazon Basin. The number of background points was further reduced by using a background point probability mask and sampling 3 * the no. of occurrence points. The code used to create the probability surface and sample the background points is in the ["0_download_gbif_points.R" code](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/R_code/0_download_gbif_points.R). An example of using a background mask to sample background points can be found in [Moyes et al. 2016](https://parasitesandvectors.biomedcentral.com/articles/10.1186/s13071-016-1527-0).

_Step 1._ Download occurrence dataset from the data folder: [a_chamek_ter_mammals_amazon_thinned_Oct22.csv](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/data/a_chamek_ter_mammals_amazon_thinned_Oct22.csv)

#### Uploading point data to GEE
>Once you have your occurrence data downloaded, upload it to GEE as a csv. Make sure the dataframe has numerical longitude (x) and latitude (y) coordinates (in decimal degrees) as separate columns. Make sure there is a row identifier to match the points to bind multiple datasets after downloading geospatial data. 

_Step 2._ Upload occurrence dataset as a GEE feature collection. Make sure to specify which are the "lat" and "lon" columns so GEE knows to take them in as numerics.</br>

&nbsp;   

<img src=https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/GEE_csv_asset.png width="250" height="380"></br>

&nbsp;   

**Figure 2.** Navigation for uploading csv as a feature collection.
 
&nbsp;  

## 2. Environmental data :deciduous_tree:
> Species distribution models model the probability a species occurs in pixel _x_ given the environmental conditions (covariates) in pixel _x_. Here, we will download land-use / land-cover and climate data to use as environmental covariates in our model. For the purposes of this tutorial, we will run the initial model with a small number of covariates. 

&nbsp;  

#### Downloading MAPBIOMAS data :earth_americas: </br>

<img src=https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/mapbiomas_amazonia_example.png width="900" height="550"></br>
**Figure 3.** The MAPBIOMAS Panamazonia platform.

&nbsp;  

_Step 3._ Explore the [MAPBIOMAS website](https://plataforma.panamazonia.mapbiomas.org/)

_Step 4._ Identify LULC categories / datasets of interest using the [MAPBIOMAS legend](https://s3.amazonaws.com/amazonia.mapbiomas.org/leyenda/C%C3%B3digo_de_la_Leyenda_-_colecci%C3%B3n_3.pdf) (make sure it is the correct legend for the MAPBIOMAS assett)

_Step 5._ Explore [MAPBIOMAS GEE code](https://github.com/mapbiomas-brazil/user-toolkit). MAPBIOMAS has pre-written code for most of the functions you need for SDMs.

_Step 6._ We will build our species distribution model at a 1km^2 resolution. MAPBIOMAS is available at a finer scale resolution (30m) so it is possible to create the model at a finer spatial resolution. Note: the resolution of all covariates should be equal to the resolution of the coarsest variable (given that the precision of your occurrence point is as or more precise than that resolution). The data we fed into GEE are points. We will create a 1km buffer around each point (step 6) and then calculate area of each land class per grid cell (step 7-10). A GEE script for steps 6-10 can be found [here](https://code.earthengine.google.com/968c3eb7d67d49c47ca18fe74ce26043).


```
//call csv of points (make sure it has lat, long, & a row identifier)
var amazon_mammals = ee.FeatureCollection('users/cglidden/a_chamek_ter_mammals_amazon_thinned_Oct22'); // replace "cglidden" with your username

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

_Step 7._ We will now calculate area per each feature (point + buffer) per year in the study period (2001-2020) in the feature collection. _A. chamek_ typically live in lowland forests but are listed as endangered due to habitat loss. For our initial model we will include lulc variables related to forest cover and farming. The MAPBIOMAS code first sets the calculations up by setting a number of parameters / variables including, among other parameters / variables, setting the identifying property of each point as 'attribute', defining territories (the features), the image (MAPBIOMAS collection), classIds (the land-class you want to calculate area for), the name of the export file, and the file to export to. You can follow along with the GEE code below or in the file on [Caroline' GEE code editor, linked here and above](https://code.earthengine.google.com/968c3eb7d67d49c47ca18fe74ce26043). The code below also includes code to print the MAPBIOMAS image to view the structure and to map the first year of MAPBIOMAS data (the first band of the image). 

<img src= https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/mapbiomas_structure.png width="900" height="450"></br>
**Figure 5.** The print of the MAPBIOMAS image. You will see that each year of LULC data is a band in the image. This is important infomation to know when it comes to deciding how to extract the data you want. 

```
// Asset = mapbiomas collection that you would like to use

var asset = 'projects/mapbiomas-raisg/public/collection3/mapbiomas_raisg_panamazonia_collection3_integration_v2';

// Numeric attribute to index the point
var attribute = "row_code";

// Output csv name
var outputName = 'a_chamek_ter_mammals_lulc_Oct2022';

// Change the scale if you need.
var scale = 30;

// Define a list of years to export
var years = ['2001','2002','2003','2004','2005','2006','2007','2008','2009',
              '2010', '2011', '2012', '2013', '2014', '2015', '2016',
              '2017', '2018', '2019', '2020'];

// A list of class ids you are interested
var classIds = [
 // 1, //Bosque
    3, // Formação Florestal
//  4, // Formação Savânica
//  5, // Mangue
//  6, // Bosque inundable
//  10, // Formación Natural no Forestal
//  11, // Formación Natural no Forestal Inundable
//  12, // Formación Campestre o herbazal
//  29, // Afloramento Rochoso
//  13, // Outra Formação não Florestal
    14, //Uso Agropecuario
//  22, // Área não Vegetada
//  24, // Infraestrutura Urbana
//  30, // Mineração
//  25, // Outra Área não Vegetada
//  26, // Cuerpo de agua
//  33, // Rio, Lago e Oceano
//  34, //Glaciar
//  27, //No observado
];    
    

// Define a Google Drive output folder 
//if this folder does not already exist, this code will create the folder in your google drive
var driverFolder = 'GEEexports';

/**
 * 
 */
/**
 * 
 */
// Territory
var territory = ee.FeatureCollection(pointBuffers);

// LULC mapbiomas image - 36 bands, one image per band, band named "classification_year"
var mapbiomas = ee.Image(asset).selfMask();
print(mapbiomas) print to view image collection

var singleBandVis = { //set mapping parameters
  min: 1, //min classId
  max: 35, //max classId
  //define color ramp that will be stretched from 1 to 35
  palette: ['08A11A', '71FDD9', 'FAD7A0', 'EBB3FC', '3498DB', '1F618D']};
Map.setCenter(-69.60, -12.39, 4) //coordinates & degree to zoom in
Map.addLayer(mapbiomas.select("classification_1985"), singleBandVis, 'MAPBIOMAS')

// Image area in km2, one band named "area"
var pixelArea = ee.Image.pixelArea().divide(1000000);

// Geometry to export
var geometry = mapbiomas.geometry();


```
<img src= https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/mapbiomas_example_1985.png width="800" height="550"></br>
**Figure 6.** A map of the MAPBIOMAS data from 1985 (produced from the code above). 

_Step 8._ Now the MAPBIOMAS code defines functions that will be used to calculate area. This code is complex and hard to understand, so don't feel overwhelemed if it does not make sense the first time we go through it. The second function takes in a band of a MAPBIOMAS image, a feature, and the feature geometry and returns area per land-class per feature (data point) as a table. The second feature initially produces a complex object because it converts the feature to an image, and then maps over the MAPBIOMAS image and the feature image to produce the sum of class pixels per feature. The first function is used within the second function to convert that complex object into a table that we can export. 

```
/**
 * Convert a complex ob to feature collection
 * @param obj 
 */
var convert2table = function (obj) {

    obj = ee.Dictionary(obj); //create dictionary 

    var territory = obj.get('territory'); //extract a property from a feature

    var classesAndAreas = ee.List(obj.get('groups')); //extract a property from a feature and create a list

    var tableRows = classesAndAreas.map( //apply to every image in image collection
        function (classAndArea) {
            classAndArea = ee.Dictionary(classAndArea);

            var classId = classAndArea.get('class'); //extract class
            var area = classAndArea.get('sum'); //extract sum of pixels

            var tableColumns = ee.Feature(null)
                .set(attribute, territory) //set row_code
                .set('class', classId) //set id to class column
                .set('area', area); //set area to area column

            return tableColumns; //return table as funtion output
        }
    );

    return ee.FeatureCollection(ee.List(tableRows));
};

/**
 * Calculate area crossing a cover map (deforestation, mapbiomas)
 * and a region map (states, biomes, municipalites)
 * @param image 
 * @param territory 
 * @param geometry
 */
 
 
 
 //the following code creates a function that takes in a band of a MAPBIOMAS image, a feature, and the feature geometry and returns areas as a table
var calculateArea = function (image, territory, geometry) {

    //create reducer to get sum per class per territory
    var reducer = ee.Reducer.sum().group(1, 'class').group(1, 'territory');

  // now create a territories object that for each 1km pixel, add info from each feature, add the image, and then sums per group combination (feature x land-class)
  //each pixel using the reducer above
    var territoriesData = pixelArea.addBands(territory).addBands(image)
        .reduceRegion({
            reducer: reducer,
            geometry: geometry,
            scale: scale,
            maxPixels: 1e12
        });

    territoriesData = ee.List(territoriesData.get('groups')); //list for each group

    var areas = territoriesData.map(convert2table); //use code above to turn territories object into a table, obj = territoriesData

    areas = ee.FeatureCollection(areas).flatten(); //flattens collection of collections

    return areas;
};
```

_Step 9._ Now we apply the above functions to the MAPBIOMAS image and our feature collection, and map over the MAPBIOMAS image so that we create a table for each year (i.e., band of the MAPBIOMAS image that we want, defined in the year list in the code above). 

```
//for each year, calculate area using functions above
var areas = years.map(
    function (year) {
        var image = mapbiomas.select('classification_' + year); //select the band for a year

        var areas = territory.map( //for each feature (territory) run the calculateArea function
            function (feature) {
                return calculateArea( //function from above where you need to define (image, feature, geometry)
                  //image: the MAPBIOMAS band defined above, where land classes that are not defined in class id list are remapped to 0
                    image.remap(classIds, classIds, 0),
                    //feature: the feature converted to an image with row_code as a property
                    ee.Image().int64().paint({
                        'featureCollection': ee.FeatureCollection(feature),
                        'color': attribute
                    }),
                    feature.geometry() //geometry: feature geometry
                );
            }
        );

        areas = areas.flatten(); //flatten collection

        // set additional properties
        areas = areas.map(
            function (feature) {
                return feature.set('year', year);
            }
        );

        return areas;
    }
);

areas = ee.FeatureCollection(areas).flatten(); //flatten collection of collections
```

_Step 10._ Now export the new feature collection (a table with area per land-class per data point) to the Google Drive folder that you named in the code above.

```
Export.table.toDrive({
    collection: areas,
    description: outputName,
    folder: driverFolder,
    fileNamePrefix: outputName,
    fileFormat: 'CSV'
});
```

&nbsp;  

#### Downloading historical climate data :thermometer: </br>
_Step 11._ Species distributions are usually, at least in part, dictated by interations between land cover and climate. We will use the [CHELSA dataset](https://chelsa-climate.org/) to add climatologies to our model. The CHELSA climatologies are averages from 1985-2010 a 1km^2 resolution, which is why we scaled the MAPBIOMAS data up to 1km. We will use the same feature collection from steps 6-10. A GEE script can be found [here](https://code.earthengine.google.com/14b1a32976d3097b5eca6be97cf84559).

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

_Step 12._ Skip actually running the GEE code for now and download the ["a_chamek_ter_mammals_lulc_Oct22.csv" file](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/data/a_chamek_ter_mammals_lulc_Oct22.csv) and the ["a_chamek_ter_mammals_climate_Oct22.csv" file](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/data/a_chamek_ter_mammals_climate_Oct22.csv) pre-downloaded data from the data folder.

&nbsp;  

#### Clean covariate data :broom:

R code for the following sections can be found in the ["1_cleaning_data.R"](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/R_code/1_cleaning_data.R)

&nbsp;  

_Step 13._ Using the data downloaded in step 12 and the code below, we will relabel MAPBIOMAS classes to make it easier to view results. We'll then aggregate LULC data by taking the mean LULC from 2001-2020. Since _A. chamek_ is affected by deforestation we will also look at difference in forest cover between 2001-2020. We will then merge the lulc data with the occurrence data and climate data. The climate data is exported from GEE in an immediately useable form. Note: Given the pace of LULC change, this is a really coarse way of aggregating the data and we likely loose a lot of signal.  

```
#download the github folder and set it to your working directory
rm(list=ls())
setwd("~/Desktop/UPCH-species-distribution-tutorial-main/")

#load libraries
library(tidyr); library(dplyr)

#-----------------------------------#
# read in datasets                  #
#-----------------------------------#

mapbiomas <- read.csv("data/a_chamek_ter_mammals_lulc_Oct2022.csv")

climate <- read.csv("data/a_chamek_ter_mammals_climate_Oct2022.csv")
#subset climate data to variables of interest (not necessary but makes the final dataset cleaner
climate <- climate[  c("row_code", "bio13_precip_wettest_month", "cmi_min")]

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

mapbiomas_wide <- mapbiomas[2:5] %>% 
  pivot_wider(names_from = class, 
              values_from = area,
              values_fn = list(area = sum),
              values_fill = list(area = 0)) 


#---------------------------------------------------------#
# summarize average area per class per point across years #
# and difference in area over the study period            #
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
# merge data frames with response (presence / absence labels), climate covariates, and lulc covariates   #
#--------------------------------------------------------------------------------------------------------#

covariates <- left_join(climate, mapbiomas_mean_diff, by = "row_code")
data0 <- left_join(amazon_basin_pnts, covariates, by = "row_code")

```

&nbsp;  

## 3. Multicollinearity analysis :raised_eyebrow: :chart_with_upwards_trend:
_Step 14._ Using the code below, we will now clean our covariate data a bit more by removing highly colinear variables. While machine learning can handle multicolinearity when making predictions, removing colinear variables can still be helpful for model interpretation. The correlation value depends on your questions and dataset but we will use a 0.7 correlation cutoff in the code below. We will use a pair-wise analysis but another option is a variable inflation analysis (or you can use both).

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
**Figure 7.** Point plot of pairwise realtionships, covariate distribution, and Pearson's correlation coeffecients of variables.

&nbsp;  

### Challenge Questions
1. Try using GEE to add additional climate or LULC covariates to the final analysis dataset.

&nbsp;  
&nbsp;  


## 4. Spatial cross-validation :white_check_mark:

_Step 15._  Next we will split our data into 3 folds (3 subsets) for 3-fold cross validation. It is important to test the perfomance of your model using a hold-out test set. This allows you to evaluate if your model is predicting generazliable patterns, or if it is only learning the traing data (and thus "overfitting"). One way to test out-of-sample model performance is using k-fold cross-validation. K-fold cross-validation splits the data into k folds, it then trains and tests the model k times (where, for each iteration, one fold is a hold out fold and the remaning folds are used for training the model). K-fold cross-validation helps to test model performance across different subsets of data where the subsets are sampled without replacement. For many applications of species distribution modeling, it is ideal to use spatial cross-validation where folds are separated in space so to avoid issues of autocorrelation that arise from test and training points being very close to each other. See figure 5. for a visual explanation. Here we will use the R package _spatialsample_. Methods for splitting folds can be dependent on your data and study questions. View the [blockCV paper (Valavi et al. 2021)](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13107) and the [_spatialsample_](https://spatialsample.tidymodels.org/) rpackage to learn of different ways to split data. Below we will use block clustering because it is quick to implement.</br>

<img src= https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/spatialcv_visualization.png width="700" height="350"></br>
**Figure 8.** Visualization of spatial partitioning versus random test versus train set. Figure from towards data science ["Spatial CV Using Sickit-learn"](https://towardsdatascience.com/spatial-cross-validation-using-scikit-learn-74cb8ffe0ab9).

```
library(spatialsample); library(sf)

#--------------------------------------------#
# get  fold  id  by  block  clustering       #
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

#sanity  check:  check  how  many  data  points  are  in  each  fold (make sure folds have adequate data in them)
table(analysis_data$fold)


#write  dataframe  to  save  for  later
write.csv(analysis_data, "data/a_chamek_ter_mammals_finalData_Oct22.csv")

```

&nbsp;  

## 5. Choice of statistical model -- Random Forest :computer:
> There are multiple statistical models (e.g., generalized additive models) and machine learning models (e.g., maxent, random forest, boosted regression tree) that can be used for SDMs. You can read more about different SDM modeling methods in [Valavi et al. 2022](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecm.1486). Here, we will use a random forest model as they typically have high accuracy (as they deal with non-linearity and higher-order interactions well) and are fast to implement. Random forest is a classification algorithm that takes the average of multiple decision trees. You can learn more about random forests by watching [this video](https://www.youtube.com/watch?v=v6VJ2RO66Ag). 

&nbsp;  

_Step 16._ Now we train the model on each set of k-1 folds and test it on the holdout fold. For each iteration, we tune the randomForest model to optimize model performance. Depending on your dataset and the parameter values you search over, the tuning step can also be used to prevent over-fitting. There are different methods for tuning a machine-learning model. Below we use a [hypergrid search](https://afit-r.github.io/random_forests#tune) and select the final parameters based on the combination that yields the best model performance. If you had trouble running _Step 10._ you can download the ["a_chamek_ter_mammals_finalData_Oct22.csv"](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/data/a_chamek_ter_mammals_finalData_Oct22.csv) to use for the next few steps. R code for the following sections can be found in the ["2_random_forest_sdm.R"](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/R_code/2_random_forest_sdm.R) </br>

```
library(ranger)

#------------------------------------#
#tune, train, model                  #
#------------------------------------#

#first  reduce  data  down  to  covariates  of  interest  (or  you  could  specify  it  in  the  formula  below)

analysis_data_v2  <-  analysis_data[  c("presence", "fold", "bio13_precip_wettest_month", "cmi_min",
                                        "mean_forest",  "mean_farming", "diff_forest_formation")]

#for many ML algorithms  you should make sure your response variables are not grouped in the data-frame (i.e.,  all 1s and 0s next to each other)
analysis_data_v2 <- analysis_data_v2[sample(1:nrow(analysis_data_v2)), ]

#create  empty  dataframe  to  for  loop  to  store  results    one  row  for  each  fold
rf_performance  <-  data.frame(model  =  rep("RF", 3),
                               fold_id  =  1:3,
                               auc  =  rep(NA, 3),
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
        hyper_grid$OOB_RMSE[j]  <-  sqrt(model$prediction.error)
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
    rf_performance[i, "presence"]  <-  nrow(subset(test, presence  ==  1))
    rf_performance[i, "background"]  <-  nrow(subset(test, presence  ==  0))
    
    #save hypergrid results to use for final model
    hypergrid_final[i, "mtry"]  <-  hyper_grid2$mtry[1]
    hypergrid_final[i, "node_size"]  <-  hyper_grid2$node_size[1]
    hypergrid_final[i, "sampe_size"]  <-  hyper_grid2$sampe_size[1]

}

#mean auc = 0.71
mean(rf_performance$auc)


```

&nbsp;  

#### Training of final model :woman_technologist:
_Step 17._ Now train the final model that will be used for interpretation and prediction. Below, I used the average hyperpameters from the tuning steps.

```
#------------------------------------------------------------#
# train  final  model                                        #
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

## 6. Model interpretation and evaluation :bar_chart: :chart_with_downwards_trend:

##### Variable importance
_Step 18._ There are many ways to calculate variable importance. Here, we will use an intuitive and model agnostic measure of variable importance. In sum, we will calculate change in model performance when a focal variable is randomly permuted, which will tell us the degree to which the variable contributes to model performance and thus accuracy of model predictions.


```
#------------------------------------------------------------#
# variable  importance                                       #
#------------------------------------------------------------#

#extract  model  results  to  get  permutation  importance
permutation_importance  <-  data.frame(variable  =  rownames(as.data.frame(final_model$variable.importance)),  
                                                                          importance  =  as.vector(final_model$variable.importance))
    
#plot  importance
ggplot(permutation_importance, aes(x  =  variable, y  =  importance))  +
    geom_bar(stat="identity")  +
    ggtitle("permutation  importance")  +
    ylab("importance (change in model error)") + 
    coord_flip()  +
    theme_classic()


```

<img src= https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/variable_importance_plot.png width="500" height="600">

**Figure 9.** Permutation variable importance for each covariate in the model. The y-axis is the change in model error when the variable is permuted.

&nbsp;  

##### Partial dependence plots
_Step 19._ Now you will visualize your model results using partial dependence plots. Partial dependence plots (PDPs) depict the relationship between the probability of species occurrence and the variable of interest across a range of values for that variable. At each value of the variable, the model is evaluated for all values of the other covariates. The final model output is the average predicted probability across all combinations of the other covariates. 



```
#------------------------------------------------------------#
# pdps                                                       #
#------------------------------------------------------------#

#try plotting a PDP for just one variable
pdp::partial(final_model, pred.var  =  "mean_forest", prob  =  TRUE, plot = TRUE, train  =  analysis_data_v2[complete.cases(analysis_data_v2), -2])
#train = data without NAs & without "fold" column

###now run a for loop to get plotting data for all variables in the model (or in the 'var_names' list
#list of covariate names to generate pdps for and loop through
var_names  <-  names(analysis_data_v2[complete.cases(analysis_data_v2),  -c(1, 2)]) #analysis dataset exlucing 'presence' and 'fold' column

#dataframe  to  make  partial  dependence  plots
pd_df  =  data.frame(matrix(vector(), 0, 3, dimnames=list(c(), c('variable', 'value', 'yhat'))),  
                     row.names  =  NULL, stringsAsFactors=F)

#loop  through  each  variable
for  (j  in  1:length(var_names))  { 
    
    output  <-  as.data.frame(pdp::partial(final_model, pred.var  =  var_names[j], prob  =  TRUE, train  = analysis_data_v2[complete.cases(analysis_data_v2), -2]))
    
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
    theme_bw(base_size = 14)

```
<img src= https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/pdp_plot.png>

**Figure 10.** Partial dependence plots for each covariate in the model.

&nbsp;  

##### Model predictions
_Step 20._ Using the code below, you will use the final model to map the distribution of _A. chamek_ across the Amazon Basin. Once you validate your model and are happy with the biological interpration gained from variable importance, you can use the model to map the distribution of the species. To do this, you create a grid of your area of interest (the area that the occurrence points were distributed over) and determine the environmental covariates in each grid cell. You then use your model to predict the probability of vector occurrence within each grid cell and display this data using a map.

```
code for generating distribution map based on geoTIFFs of prediction variables? might need to create it for a small area, otherwise files might be too big for github
```

&nbsp; 

### Challenge questions
1. How would you add uncertainty to variable importance, functional relationships (pdps), or prediction maps?
2. How does you model perfromance and interpretation change when you add / remove covariates?
3. Pick a new focal species to create a SDM. You will have to use the "notThinned" datasets to select a focal species and create background points.
