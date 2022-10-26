setwd("~/Documents/GitHub/UPCH-species-distribution-tutorial")

##note that folder structure of the github may have changed since writing this code

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
ab <- st_union(ab0)

pnts_sf <- st_as_sf(g, coords = c('decimalLongitude', 'decimalLatitude'), crs = st_crs(ab)) #make sure crs is same between points & shape file

pnts_sf$indicator <- st_within(pnts_sf, ab) %>% lengths > 0 #indicator = 1 if points falls in MDD

pnts_ab <- subset(pnts_sf, indicator == TRUE)

pnts_ab <- pnts_ab %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) # add lat/lon as unique columns, and drop geometry below (converts data to a regular dataframe)

pnts_ab <- st_drop_geometry(pnts_ab) #alternatively, you could save this as a shape file and keep the geometry

write.csv(pnts_ab,"data/ter_mammals_amazon_notThinned_Oct22.csv", row.names = TRUE) #save output with all datapoints (pre-thinning)


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
#visually check distribution of each species                  #
#-------------------------------------------------------------#

# Remove plot axis
library(ggplot2)

no_axis <- theme(axis.title=element_blank(),
                 axis.text=element_blank(),
                 axis.ticks=element_blank())

#subset species with a high number of occ points, then look at distribution
pts_per_species <- as.data.frame(table(pnts_ab$scientificName))
quantile(pts_per_species$Freq, prob = 0.75) #83 points
plot_species <- subset(pts_per_species, Freq > 83); names <- as.vector(plot_species$Var1)

# Plot each species
ggplot() +
  geom_sf(data=ab, color="#2D3E50", fill="lightgrey", size=.15, show.legend = FALSE) +
  geom_point(data = subset(pnts_ab, scientificName == names[1]), 
             aes(x = lon, y = lat), alpha = 0.5) +
  theme_minimal() +
  no_axis

#candidate species for amazon terrestrial mammals 
#Akodon dayi - 45 unique obs
#Akodon subfuscus - 49 unique obs
#Alouatta macconnelli - 87 unique obs
#Ametrida centurio - 60 unique obs
#Anoura peruana - 45 unique obs
#Aotus vociferans - 75 unique obs
#Artibeus concolor - 80 unique obs
#Artibeus obscurus - 89 unique obs
#Atelocynus microtis - 64 unique obs (short-eared dog) #2
#Bradypus tridactylus - 189 unique obs (pale throated sloth)
#Cormura brevirostris - 111 unique obs
#Euryoryzomys nitidus - 87 unique obs
#Marmosops noctivagus - 86 unique obs (white-bellied slender opossum) #1
#Philander andersoni - 63 unique obs
#Rhinophylla fischerae - 85 unique obs
#Sturnira magna - 67 unique obs #3
#Uroderma magnirostrum - 79 obs


#-------------------------------------------------------------#
#double check how many unique points in a 100m raster         #
#-------------------------------------------------------------#

#make raster with 1000m grid cells
st_bbox(ab) #get limits to put in raster
r <- raster(xmn = -79.699771, xmx = -44.491086, ymn = -20.493752, ymx = 8.663513, res = 0.0083)

#one point per grid cell
s <- dismo::gridSample(pnts_ab[pnts_ab$scientificName == "Atelocynus microtis", c("lon", "lat")], r, n=1) #64 obs for focal species

s0 <- dismo::gridSample(pnts_ab[pnts_ab$scientificName != "Atelocynus microtis", c("lon", "lat")], r, n=1) #6914 obs for background species

#-------------------------------------------------------------#
#Create background mask using probability sampling            #
#-------------------------------------------------------------#

background <- pnts_ab[pnts_ab$scientificName != "Atelocynus microtis", ] #exclude focal species

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
set.seed(909)
multiplier <- 2

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
write.csv(final_pass, "data/a_microtis_ter_mammals_amazon_thinned_Oct22.csv")

#-------------------------------------------------------------#
#final figure to visualize distribution of points             #
#-------------------------------------------------------------#

# Remove plot axis
no_axis <- theme(axis.title=element_blank(),
                 axis.text=element_blank(),
                 axis.ticks=element_blank())

##for mapping purposes
final_pass$species <- ifelse(final_pass$scientificName == "Atelocynus microtis", 
                             "Atelocynus microtis", "background species")

# Plot each species
point_distribution <- ggplot() +
  geom_sf(data=ab, color="#2D3E50", fill="lightgrey", size=0.15, show.legend = FALSE) +
  geom_jitter(data = final_pass, 
             aes(x = lon, y = lat, color = species), size = 0.5, alpha = 0.75) +
  theme_minimal() +
  no_axis

ggsave("final_figures/a_microtis_sdm_point_distribution.png", point_distribution, dpi = 300)

