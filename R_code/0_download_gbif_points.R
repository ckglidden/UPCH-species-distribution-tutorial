setwd("~/Documents/GitHub/UPCH-species-distribution-tutorial")

##note that folder structure of the github may have changed since writing this code

library(sf); library(rgbif); library(dplyr); library(raster); library(ggplot2)

#-------------------------------------------------------------#
#download gbif data                                           #
#-------------------------------------------------------------#

#read in IUCN species names for songbird in Peru
iucn_species <- read_sf("iucn_peruvian_passerines") #this shape file was downloaded from ICUN redlist site
names(iucn_species)[3] <- "scientificName" #needs to match gbif to work

#create a dataframe with the scientific name & all occ of songbirds
latlong <- data.frame()

#for loop to get lat long of each species from GBIF
for (i in 1:nrow(iucn_species)){
  tryCatch({
    b <- data.frame(occ_search(scientificName = iucn_species$scientificName[i])$data)
    #filter so that only the recorded occurrences are in the dataframe
    if("decimalLatitude" %in% colnames(b)){
      c <- subset(b, select=c("scientificName", "decimalLatitude", "decimalLongitude", "year"))
      #some of the occurrences have weird names, so we just rename that column to their scientific name
      c[1:nrow(b), 1] <- iucn_species$scientificName[i]
      latlong <- rbind(latlong, c)
    }
  }, error=function(e){})
}

#get rid of all the occurrences with no lat longs and create a csv
g <- latlong[complete.cases(latlong), ]

#clip occurrence points to only ones that fall in MDD
mdd0 <- read_sf("Madre_de_Dios"); 
mdd <- st_union(mdd0)

pnts_sf <- st_as_sf(g, coords = c('decimalLongitude', 'decimalLatitude'), crs = st_crs(mdd))

pnts_sf$indicator <- st_within(pnts_sf, mdd) %>% lengths > 0 #indicator = if points falls in MDD
  
pnts_mdd <- subset(pnts_sf, indicator == TRUE)

pnts_mdd <- pnts_mdd %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) # add lat/lon as unique columns, and drop geometry below (converts data to a regular dataframe)

pnts_mdd <- st_drop_geometry(pnts_mdd) 

write.csv(pnts_mdd,"passerine_occ_pts_mdd.csv", row.names = TRUE)

#-------------------------------------------------------------#
#visually check distribution of each species                  #
#-------------------------------------------------------------#

# Remove plot axis
library(ggplot2)

no_axis <- theme(axis.title=element_blank(),
                 axis.text=element_blank(),
                 axis.ticks=element_blank())

names <- unique(pnts_mdd$scientificName) #look at distribution by species
#look by genus or family?

# Plot each species
ggplot() +
  geom_sf(data=mdd, color="#2D3E50", fill="lightgrey", size=.15, show.legend = FALSE) +
  geom_point(data = subset(pnts_mdd, scientificName == names[15]), 
             aes(x = lon, y = lat), alpha = 0.5) +
  theme_minimal() +
  no_axis

#candidate species
#"Poecilotriccus albifacies" - 478, 59 unique points
#"Corythopis torquatus" - 244, 65 unique points

#-------------------------------------------------------------#
#double check how many unique points in a 100m raster         #
#-------------------------------------------------------------#

#make raster with 100m grid cells
st_bbox(mdd)
r <- raster(xmn = -72.428739, xmx = -68.652279, ymn = -13.341717, ymx = -9.873393, res = 0.001)

#one point per grid cell -- 245 unique grid points distributed throughout BR
s <- dismo::gridSample(pnts_mdd[pnts_mdd$scientificName == "Corythopis torquatus", c("lon", "lat")], r, n=1) #65 obs for focal species

s0 <- dismo::gridSample(pnts_mdd[pnts_mdd$scientificName != "Corythopis torquatus", c("lon", "lat")], r, n=1) #178 obs for background species

#-------------------------------------------------------------#
#Create background mask using probability sampling            #
#-------------------------------------------------------------#

background <- pnts_mdd[pnts_mdd$scientificName != "Corythopis torquatus", ]

bg_species_list <- unique(background$scientificName)


#------------------------------------------------------
# Extract number of songbirds per grid cell
#------------------------------------------------------

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
                       agr = "constant", crs = 4326)

# Random sample bg without replacement from weighted bias mask at (2x occ) multiplier
set.seed(9)
multiplier <- 2

bg_mask_weights <- bg_mask_sf %>%
  mutate(weight = count/sum(count))

bg_mask_df <- bg_mask_sf[sample(nrow(bg_mask_weights),
                                size = multiplier * nrow(s), # s is from earlier code chunk, # of occ points
                                replace = FALSE,
                                prob = bg_mask_weights$weight),]

#make bkg dataset match presence
bg_mask_df <- bg_mask_df %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]); bg_mask_df <- st_drop_geometry(bg_mask_df)
names(bg_mask_df)[c(4)] <- c("year"); bg_mask_df <- bg_mask_df[, c("scientificName", "year", "lon", "lat")]
#make sure bkg points are labeled
bg_mask_df$presence <- 0

#subset presence points to thin set, make sure it is lableled
occ_points <- pnts_mdd[row.names(s), c("scientificName", "year", "lon", "lat")]; occ_points$presence <- 1 

#final passerine occ set
final_pass <- rbind(occ_points, bg_mask_df)

#add in row identifier for GEE
final_pass$row_code <- seq(1, nrow(final_pass), by = 1)
write.csv(final_pass, "final_passerine_dataset_Oct20_2022.csv")

#-------------------------------------------------------------#
#final figure to visualize distribution of points             #
#-------------------------------------------------------------#

# Remove plot axis
no_axis <- theme(axis.title=element_blank(),
                 axis.text=element_blank(),
                 axis.ticks=element_blank())

##for mapping purposes
final_pass$species <- ifelse(final_pass$scientificName == "Corythopis torquatus", 
                             "C_torquatus", "bkg_species")

# Plot each species
point_distribution <- ggplot() +
  geom_sf(data=mdd, color="#2D3E50", fill="lightgrey", size=.15, show.legend = FALSE) +
  geom_point(data = final_pass, 
             aes(x = lon, y = lat, color = species), size = 1, alpha = 0.5) +
  theme_minimal() +
  no_axis

ggsave("c_torquatus_sdm_point_distribution.tiff", point_distribution, dpi = 300)

