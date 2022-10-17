# UPCH geospatial analysis workshop
## Species distribution modeling tutorial
> This github page includes a tutorial of how to model the distribution of XX species using MAPBIOMAS land-use/land-cover data.


### 1. Uploading point data to GEE :mosquito:

a. download occurrence dataset (XYZ - csv file available on this github page)

b. upload occurrence dataset as an assett


### 2. Downloading MAPBIOMAS data :earth_americas:

Maybe include a screengrab of MAPBIOMAS?

1. Explore the [MAPBIOMAS website](https://plataforma.brasil.mapbiomas.org/?activeBaseMap=9&layersOpacity=100&activeModule=coverage&activeModuleContent=coverage%3Acoverage_main&activeYear=2021&mapPosition=-15.072124%2C-51.416016%2C4&timelineLimitsRange=1985%2C2021&baseParams[territoryType]=1&baseParams[territories]=1%3BBrasil%3B1%3BPa%C3%ADs%3B-33.751177993999946%3B-73.9904499689999%3B5.271841077000019%3B-28.847639913999956&baseParams[activeClassesLevelsListItems]=1%2C7%2C8%2C9%2C10%2C2%2C11%2C12%2C13%2C14%2C15%2C16%2C3%2C17%2C18%2C27%2C37%2C38%2C39%2C40%2C41%2C28%2C42%2C43%2C44%2C19%2C20%2C4%2C21%2C22%2C23%2C24%2C5%2C25%2C26%2C6)

2. Identify LULC categories of interest using the [MAPBIOMAS legend](https://mapbiomas.org/en/legend-codes) (make sure it is the correct legend for the MAPBIOMAS assett)

3. navigate to Caroline's GEE code - in this code you will upload the assett, specifiy the output file & directory, and choose the LULC and years you want to download data for


### 3. Clean MAPBIOMAS data :broom:

```

cleaning code chunk - maybe take average of LULC over a ten year period for simplicty sake?

```

