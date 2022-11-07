#set your working directory
setwd("~/Desktop/UPCH-species-distribution-tutorial-main/")

##note that folder structure of the github may have changed since writing this code
#install and load packages
#install.packages(c("sf","rbif", "dplyr", "raster", "ggplot2"))
library(sf); library(rgbif); library(dplyr); library(raster); library(ggplot2)

#-------------------------------------------------------------#
#download gbif data                                           #
#-------------------------------------------------------------#

#read in IUCN species names for songbird in Peru
iucn_species <- read_sf("/Users/carolineglidden/Documents/GitHub/iucn_MDD_species/neotropics_terrestrial_mammals") #this shape file was downloaded from ICUN redlist site
names(iucn_species)[3] <- "scientificName" #needs to match gbif to work
unique_names <- unique(iucn_species$scientificName) #names to loop through, can also just use this function in the loop

#create a dataframe with the scientific name & all occ of songbirds
latlong <- data.frame()

#for loop to get lat long of each species from GBIF
for (i in 1:nrow(iucn_species)){
  tryCatch({
    b <- data.frame(occ_search(scientificName = unique_names[i])$data)
    #filter so that only the recorded occurrences are in the dataframe
    if("decimalLatitude" %in% colnames(b)){
      c <- subset(b, select=c("scientificName", "decimalLatitude", "decimalLongitude", "year"))
      #some of the occurrences have weird names, so we just rename that column to their scientific name
      c[1:nrow(b), 1] <- unique_names[i]
      latlong <- rbind(latlong, c)
    }
  }, error=function(e){})
}

#get rid of all the occurrences with no lat longs and create a csv
g <- latlong[complete.cases(latlong), ]

######For Amazon specific occurrence points
ab0 <- read_sf("Amazon_Basin")
#sf_use_s2(FALSE) # may need to switch off spherical geometry
ab <- st_union(ab0); plot(ab) # create one figure and plot it to make sure it looks correct

pnts_sf_end <- st_as_sf(g, coords = c('decimalLongitude', 'decimalLatitude'), crs = st_crs(ab)) #make sure crs is same between points & shape file

pnts_sf_end$indicator <- st_within(pnts_sf_end, ab) %>% lengths > 0 #indicator = 1 if points falls in MDD

pnts_ab_end <- subset(pnts_sf_end, indicator == TRUE)

pnts_ab_end <- pnts_ab_end %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) # add lat/lon as unique columns, and drop geometry below (converts data to a regular dataframe)

pnts_ab_end <- st_drop_geometry(pnts_ab_end) #alternatively, you could save this as a shape file and keep the geometry

#######For Madre de Dios specific occurrence points
#clip occurrence points to only ones that fall in MDD
# mdd0 <- read_sf("Madre_de_Dios"); 
# mdd <- st_union(mdd0)
# # 
# pnts_sf <- st_as_sf(g, coords = c('decimalLongitude', 'decimalLatitude'), crs = st_crs(mdd))
# pnts_sf$indicator <- st_within(pnts_sf, mdd) %>% lengths > 0 #indicator = if points falls in MDD
#   
# pnts_mdd <- subset(pnts_sf, indicator == TRUE)
# 
# pnts_mdd <- pnts_mdd %>%
#   dplyr::mutate(lon = sf::st_coordinates(.)[,1],
#                 lat = sf::st_coordinates(.)[,2]) # add lat/lon as unique columns, and drop geometry below (converts data to a regular dataframe)
# 
# pnts_mdd <- st_drop_geometry(pnts_mdd) 
# 
# write.csv(pnts_mdd,"data/XXX_mdd.csv", row.names = TRUE) #change XXX with identifiers of species

#-------------------------------------------------------------#
#add invasive species                                         #
#-------------------------------------------------------------#

#the above list of IUCN species only included endemic species in the Amazon, we can use thist paper to add points from invasive species
#https://esajournals.onlinelibrary.wiley.com/doi/epdf/10.1002/ecy.3115
#you could also use the paper to generate a list of species and download points off of gbif

alien_species <- read.csv("data/NEOTROPICAL_ALIEN_MAMMALS_OCCURENCE_v1_0.csv")
alien_species <- alien_species[,c("LONG_X", "LAT_Y", "SPECIES", "RECORD_YEAR")] #subset to columns of interest
names(alien_species) <- c("lon", "lat", "scientificName", "year") #name to match endemic data
alien_species <- alien_species[complete.cases(alien_species), ] #only include pnts with all data

pnts_sf_inv <- st_as_sf(alien_species, coords = c('lon', 'lat'), crs = st_crs(ab)) #make sure crs is same between points & shape file

pnts_sf_inv$indicator <- st_within(pnts_sf_inv, ab) %>% lengths > 0 #indicator = 1 if points falls in MDD

pnts_ab_inv <- subset(pnts_sf_inv, indicator == TRUE)

pnts_ab_inv <- pnts_ab_inv %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) # add lat/lon as unique columns, and drop geometry below (converts data to a regular dataframe)

pnts_ab_inv <- st_drop_geometry(pnts_ab_inv) #alternatively, you could save this as a shape file and keep the geometry

pnts_ab <- rbind(pnts_ab_end, pnts_ab_inv) #could add another indicator here distinguishing invasive from endemic
write.csv(pnts_ab,"data/ter_mammals_amazon_notThinned_Oct22.csv", row.names = TRUE) #save output with all datapoints (pre-thinning)


#--------------------------------------------------------------------------------------------------------------#
#visually check distribution of each species to decide focal species / look at background points               #
#--------------------------------------------------------------------------------------------------------------#
#usually you know your focal species, so this step is just for the purpose of trying out different models for different focal species
#let's first subset the data to between 2010 - 2020, so the points occur within the dates that MAPBIOMAS has available

#also, we will take average land-use / land-cover, a long-term average is coarser resolution than short term -- the short (er) time period might retain more signal

pnts_ab <- subset(pnts_ab, year > 2000 & year < 2021)

# Remove plot axis
library(ggplot2)

#subset species with a high number of occ points, then look at distribution
pts_per_species <- as.data.frame(table(pnts_ab$scientificName))
quantile(pts_per_species$Freq, prob = 0.75) #47 points
plot_species <- subset(pts_per_species, Freq > 47); names <- as.vector(plot_species$Var1)

# Plot each species - slow process but still useful to visualize data
ggplot() +
  geom_sf(data=ab, color="#2D3E50", fill="lightgrey", size=.15, alpha = 0.5, show.legend = FALSE) +
  geom_point(data =  subset(pnts_ab, scientificName != names[1]),
             aes(x = lon, y = lat), color = "grey", alpha = 0.5) + #plot potential bkg points
  geom_point(data = subset(pnts_ab, scientificName == names[1]),
             aes(x = lon, y = lat), color = "blue", alpha = 0.5) + #plot focal species
  theme_minimal()


#candidate species for amazon terrestrial mammals (enough presence points that do not completely overlap with background)
#Ateles chamek - 87 unique points (endangered, Peruvian spider monkey)


#-------------------------------------------------------------#
#double check how many unique points in a 1000m raster         #
#-------------------------------------------------------------#

#make raster with 1000m sq grid cells
st_bbox(ab) #get limits to put in raster
r <- raster(xmn = -79.699771, xmx = -44.491086, ymn = -20.493752, ymx = 8.663513, res = 0.0083)

#one point per grid cell
s <- dismo::gridSample(pnts_ab[pnts_ab$scientificName == "Ateles chamek", c("lon", "lat")], r, n=1) #87 obs for focal species

s0 <- dismo::gridSample(pnts_ab[pnts_ab$scientificName != "Ateles chamek", c("lon", "lat")], r, n=1) #6914 obs for background species

#-------------------------------------------------------------#
#Create background mask using probability sampling            #
#-------------------------------------------------------------#

background <- pnts_ab[pnts_ab$scientificName != "Ateles chamek", ] #exclude focal species

bg_species_list <- unique(background$scientificName)


#-----------------------------------------------------#
# Extract number of per grid cell                     #
#-----------------------------------------------------#

bg_points <- background %>% dplyr::select(c(lon, lat)) %>%
  as.matrix()

bg_longlat <- cellFromXY(r, bg_points) %>% as.data.frame() %>%
  cbind(background$year, background$scientificName) %>%
  mutate(count = 1) %>% setNames(c("cell","year", "scientificName","count")) %>%
  group_by(cell) %>% dplyr::summarize(count = sum(count),
                                      scientificName = scientificName,
                                      max_year = max(year),
                                      avg_year = mean(year)) %>%
  arrange(desc(count)) %>%
  mutate(lon = xFromCell(r, cell),  # Acquire longitude (x) and latitude (y) from cell centroids
         lat = yFromCell(r, cell)) %>%
  dplyr::select(-cell) %>% # Cell number is now obsolete, since will be working from (x,y) as an sf object
  filter(!is.na(lon) & !is.na(lat)) # Remove the NA locations

bg_mask_sf <- st_as_sf(bg_longlat, coords = c("lon","lat"),
                       agr = "constant", remove = FALSE, crs = 4326)

# Random sample bg without replacement from weighted bias mask at (1.5x occ) multiplier
set.seed(99)
multiplier <- 3

bg_mask_weights <- bg_mask_sf %>%
  mutate(weight = count/sum(count))

bg_mask_df <- bg_mask_sf[sample(nrow(bg_mask_weights),
                                size = multiplier * nrow(s), # s is from earlier code chunk, # of occ points
                                replace = FALSE,
                                prob = bg_mask_weights$weight),]

#make bkg dataset match presence
bg_mask_df <- st_drop_geometry(bg_mask_df)
names(bg_mask_df)[c(4)] <- c("year"); bg_mask_df <- bg_mask_df[, c("scientificName", "year", "lon", "lat")]
#make sure bkg points are labeled
bg_mask_df$presence <- 0

#subset presence points to thin set, make sure it is lableled
occ_points <- pnts_ab[row.names(s), c("scientificName", "year", "lon", "lat")]; occ_points$presence <- 1 

#final passerine occ set
final_pass <- rbind(occ_points, bg_mask_df)

#add in row identifier for GEE
final_pass$row_code <- seq(1, nrow(final_pass), by = 1)
write.csv(final_pass, "data/a_chamek_ter_mammals_amazon_thinned_Oct22.csv")

#-------------------------------------------------------------#
#final figure to visualize distribution of points             #
#-------------------------------------------------------------#


# Plot each species
point_distribution <- ggplot() +
  geom_sf(data=ab, color="#2D3E50", fill="lightgrey", alpha = 0.5, size=.15, show.legend = FALSE) +
  geom_jitter(data =  subset(final_pass, scientificName != "Ateles chamek"),
             aes(x = lon, y = lat), color = "darkgrey", alpha = 0.75) + #plot potential bkg points
  geom_jitter(data = subset(final_pass, scientificName == "Ateles chamek"),
             aes(x = lon, y = lat), color = "blue", alpha = 0.5) + #plot focal species
  theme_minimal()


ggsave("final_figures/a_chamek_sdm_point_distribution.png", point_distribution, dpi = 300)

