#==============================================================================
# Species Distribution Modeling: Data Preparation
#==============================================================================
# Purpose: Download GBIF occurrence data for Peruvian mammals, spatially filter
#          to Amazon Basin, thin occurrence points to reduce spatial bias, and
#          create presence/background dataset for SDM modeling

#set your working directory
setwd("~/Desktop/UPCH-species-distribution-tutorial-main/")

##note that folder structure of the github may have changed since writing this code
#install and load packages
#install.packages(c("sf","rbif", "dplyr", "raster", "ggplot2"))
library(sf); library(rgbif); library(dplyr); library(purrr); library(ggplot2)

#-------------------------------------------------------------#
#download gbif data                                           #
#-------------------------------------------------------------#

#read in IUCN species names for mammals in Peru
iucn_species <- read.csv('data/peru_mammal_list.csv') %>% pull(scientificName)

final_species <- c("Ateles chamek",  # randomly sampling species + focal species for this example bc download is slow
                   sample(setdiff(iucn_species, "Ateles chamek"), 999))

# Download with rate limiting and progress bar
latlong <- map_dfr(final_species[1:10], function(species) {
  Sys.sleep(0.3)  # 300ms delay = ~200 requests/minute (GBIF limit is 300/min)
  
  tryCatch({
    result <- occ_search(scientificName = species, 
                         fields = c("scientificName", "decimalLatitude", 
                                    "decimalLongitude", "year"),
                         limit = 1000)
    
    if(!is.null(result$data)) {
      return(result$data)
    } else {
      return(NULL)
    }
    
  }, error = function(e) {
    message("Failed: ", species)
    return(NULL)
  })
}, .progress = TRUE)  # shows progress bar

# Remove rows with missing coordinates
g <- latlong %>% 
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude))

######For Amazon specific occurrence points
ab0 <- read_sf("Amazon_Basin")
sf_use_s2(FALSE)
pnts_sf_end <- st_as_sf(g, coords = c('decimalLongitude', 'decimalLatitude'),
                        crs = st_crs(ab0))

# Just get points inside the basin
points_inside <- st_filter(pnts_sf_end, ab0)

pnts_ab_end <- points_inside %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) # add lat/lon as unique columns, and drop geometry below (converts data to a regular dataframe)

pnts_ab <- as.data.frame(st_drop_geometry(pnts_ab_end)) # alternatively, you could save this as a shape file and keep the geometry
# write.csv(pnts_mdd,"data/XXX_ab.csv", row.names = TRUE) #change XXX with identifiers of species

#######For Madre de Dios specific occurrence points
# Clip occurrence points to only ones that fall in MDD
# mdd0 <- read_sf("Madre_de_Dios")
# sf_use_s2(FALSE)
# 
# pnts_sf <- st_as_sf(g, coords = c('decimalLongitude', 'decimalLatitude'), 
#                     crs = st_crs(mdd0))
# 
# pnts_mdd <- st_filter(pnts_sf, mdd0)
# 
# pnts_mdd <- pnts_mdd %>%
#   dplyr::mutate(lon = sf::st_coordinates(.)[,1],
#                 lat = sf::st_coordinates(.)[,2]) # add lat/lon as unique columns, and drop geometry below (converts data to a regular dataframe)
# 
# pnts_mdd <- st_drop_geometry(pnts_mdd)
# 
# write.csv(pnts_mdd,"data/XXX_mdd.csv", row.names = TRUE) #change XXX with identifiers of species

#----------------------------------------------------------------------------------------------------------------------#
# Visually check distribution of each species to choose focal species / inspect background points   - can skip this    #
#----------------------------------------------------------------------------------------------------------------------#
# Usually you already know your focal species; this step is mainly for exploring candidate species
# and understanding how presence points compare to the available background.
# First, subset points to 2001–2020 to align with MAPBIOMAS availability.

pnts_ab <- read.csv("data/ter_mammals_amazon_notThinned_Oct22.csv") ## READ IN DATA TO OUTLINE THINNING STEPS

pnts_ab <- subset(pnts_ab, year > 2000 & year < 2021)

# Subset species with a high number of occurrence points, then inspect distributions
pts_per_species <- as.data.frame(table(pnts_ab$scientificName))
quantile(pts_per_species$Freq, prob = 0.75) # 47 points
plot_species <- subset(pts_per_species, Freq > 47)
names <- as.vector(plot_species$Var1)

# Plot each species (slow but useful for visual checks)
ggplot() +
  geom_sf(data = ab0, color = "#2D3E50", fill = "lightgrey", size = .15, alpha = 0.5, show.legend = FALSE) +
  geom_point(
    data = subset(pnts_ab, scientificName != names[1]),
    aes(x = lon, y = lat), color = "grey", alpha = 0.5
  ) + # potential background points
  geom_point(
    data = subset(pnts_ab, scientificName == names[1]),
    aes(x = lon, y = lat), color = "blue", alpha = 0.5
  ) + # focal species
  theme_minimal()

# Candidate species for Amazon terrestrial mammals (enough presence points that do not completely overlap with background)
# Ateles chamek - 87 unique points (endangered, Peruvian spider monkey)

#-------------------------------------------------------------#
# Double check how many unique points per 1000 m raster cell   #
#-------------------------------------------------------------#

st_bbox(ab0) # get limits to set raster extent
r <- raster(
  xmn = -79.699771, xmx = -44.491086,
  ymn = -20.493752, ymx = 8.663513,
  res = 0.0083
)

# One point per grid cell
s  <- dismo::gridSample(
  pnts_ab[pnts_ab$scientificName == "Ateles chamek (Humboldt, 1812)", c("lon", "lat")],
  r, n = 1
) # 87 obs for focal species

s0 <- dismo::gridSample(
  pnts_ab[pnts_ab$scientificName != "Ateles chamek (Humboldt, 1812)", c("lon", "lat")],
  r, n = 1
) # 6914 obs for background species

#-------------------------------------------------------------#
# Create background mask using probability sampling            #
#-------------------------------------------------------------#

background <- pnts_ab[pnts_ab$scientificName != "Ateles chamek", ] # exclude focal species

#-----------------------------------------------------#
# Extract number per grid cell                         #
#-----------------------------------------------------#

bg_points <- background %>%
  dplyr::select(lon, lat) %>%
  as.matrix()

bg_longlat <- cellFromXY(r, bg_points) %>%
  as.data.frame() %>%
  cbind(background$year, background$scientificName) %>%
  mutate(count = 1) %>%
  setNames(c("cell", "year", "scientificName", "count")) %>%
  group_by(cell) %>%
  dplyr::summarize(
    count = sum(count),
    scientificName = scientificName,
    max_year = max(year),
    avg_year = mean(year)
  ) %>%
  arrange(desc(count)) %>%
  mutate(
    lon = xFromCell(r, cell),
    lat = yFromCell(r, cell)
  ) %>%
  dplyr::select(-cell) %>%
  filter(!is.na(lon) & !is.na(lat))

bg_mask_sf <- st_as_sf(
  bg_longlat,
  coords = c("lon", "lat"),
  agr = "constant",
  remove = FALSE,
  crs = 4326
)

# Randomly sample background points without replacement, weighted by sampling intensity
set.seed(99)
multiplier <- 2

bg_mask_weights <- bg_mask_sf %>%
  mutate(weight = count / sum(count))

bg_mask_df <- bg_mask_sf[
  sample(
    nrow(bg_mask_weights),
    size = multiplier * nrow(s),
    replace = FALSE,
    prob = bg_mask_weights$weight
  ),
]

# Format background dataset
bg_mask_df <- st_drop_geometry(bg_mask_df)
names(bg_mask_df)[c(4)] <- c("year")
bg_mask_df <- bg_mask_df[, c("scientificName", "year", "lon", "lat")]
bg_mask_df$presence <- 0

# Format presence points
occ_points <- pnts_ab[row.names(s), c("scientificName", "year", "lon", "lat")]
occ_points$presence <- 1

# Final presence/background dataset
final_pass <- rbind(occ_points, bg_mask_df)

# Add row identifier for GEE
final_pass$row_code <- seq(1, nrow(final_pass), by = 1)
write.csv(final_pass, "data/a_chamek_ter_mammals_amazon_thinned_Oct22.csv")

#-------------------------------------------------------------#
# Final figure to visualize distribution of points             #
#-------------------------------------------------------------#

point_distribution <- ggplot() +
  geom_sf(data = ab, color = "#2D3E50", fill = "lightgrey", alpha = 0.5, size = .15, show.legend = FALSE) +
  geom_jitter(
    data = subset(final_pass, scientificName != "Ateles chamek"),
    aes(x = lon, y = lat), color = "darkgrey", alpha = 0.75
  ) +
  geom_jitter(
    data = subset(final_pass, scientificName == "Ateles chamek"),
    aes(x = lon, y = lat), color = "blue", alpha = 0.5
  ) +
  theme_minimal()

ggsave("final_figures/a_chamek_sdm_point_distribution.png",
       point_distribution, dpi = 300)
