# UPCH geospatial analysis workshop

## Species distribution modeling tutorial
> This github page includes a tutorial of how to model the distribution of _Corythopis torquatus_ (the Ringed antpipit) using MAPBIOMAS land-use/land-cover data.

1. [Presence-background species distribution models](https://github.com/ckglidden/UPCH-species-distribution-tutorial/edit/main/README.md#1-presence-background-species-distribution-models)

2. [Environmental covariate data - MAPBIOMAS LULC](https://github.com/ckglidden/UPCH-species-distribution-tutorial/edit/main/README.md#2-environmental-covariate-data)

3. [Machine-learning based SDMs](https://github.com/ckglidden/UPCH-species-distribution-tutorial/edit/main/README.md#3-machine-learning-based-sdms)


### 1. Presence-background species distribution models
> For this model we are using _C. torquatus_ as a focal species, and all other passerine birds as background species. The background species helps to  us to understand the difference between the focal species and the average landscape over which songbirds are sampled. Point data was downloaded from GBIF using the "download_gbif_points.R" code. Note: there are many ways that SDMs can be modeled when lacking true absences of the focal species.

![Figure 1. Distribution of points](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/final_figures/c_torquatus_sdm_point_distribution.png)

#### Uploading point data to GEE :mosquito:

a. download occurrence dataset: [final_passerine_data_Oct20_2022.csv](https://github.com/ckglidden/UPCH-species-distribution-tutorial/blob/main/data/final_passerine_dataset_Oct20_2022.csv)

b. upload occurrence dataset as an assett


### 2. Environmental covariate data
> short blurb on how environmental data is used to predict the distribution of focal species

#### Downloading MAPBIOMAS data :earth_americas:

Maybe include a screengrab of MAPBIOMAS?

1. Explore the [MAPBIOMAS website](https://plataforma.panamazonia.mapbiomas.org/)

2. Identify LULC categories of interest using the [MAPBIOMAS legend](https://mapbiomas.org/en/legend-codes) (make sure it is the correct legend for the MAPBIOMAS assett)

3. navigate to Caroline's GEE code - in this code you will upload the assett, specifiy the output file & directory, and choose the LULC and years you want to download data for


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
