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

Step 1. Download occurrence dataset from the data folder: [final_passerine_data_Oct20_2022.csv](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/data/final_passerine_dataset_Oct20_2022.csv)

Step 2. Upload occurrence dataset as a GEE feature collections


### 2. Environmental covariate data
> short blurb on how environmental data is used to predict the distribution of focal species

#### Downloading MAPBIOMAS data :earth_americas:

Maybe include a screengrab of MAPBIOMAS?

Step 3. Explore the [MAPBIOMAS website](https://plataforma.panamazonia.mapbiomas.org/)

Step 4. Identify LULC categories of interest using the [MAPBIOMAS legend](https://mapbiomas.org/en/legend-codes) (make sure it is the correct legend for the MAPBIOMAS assett)

Step 5. navigate to Caroline's GEE code -- we will walk through this code as a group, it will let you download LULC data per year for each point in the feature collection.

Step 6. Skip actually running the GEE code for now (or try running it in the background) and download this csv of pre-downloaded data (include link) from the data folder.


#### Clean MAPBIOMAS data :broom:

```

cleaning code chunk - maybe take average of LULC over a ten year period for simplicty sake?

```
```
vif or correlation analysis to remove highly correlated variables
```

### 3. Machine-learning based SDMs 
> short blurb on different algorithms that have been used and why we are using one specific one for the workshop (RF?)


#### Model tuning & testing using spatial cross-validation :white_check_mark:

```
code chunk to demonstrate how points are clustered for spatial CV
```

couple sentences about what needs to be tuned for RF, how to account for unbalanced data, and different metrics for testing the data?

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
