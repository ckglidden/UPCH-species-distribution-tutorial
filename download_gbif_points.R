setwd("~/Documents/GitHub/UPCH-species-distribution-tutorial")

library(sf); library(rgbif); library(dplyr); library(raster)

#read in IUCN species names for songbird in Peru
iucn_species <- read_sf("iucn_peruvian_passerines") #this shape file was downloaded from ICUN redlist site
names(iucn_species)[3] <- "scientificName" #needs to match gbif to work

#create a dataframe with the scientific name, and all the occurrences of the sandfly
latlong <- data.frame()

#for loop to get lat long
for (i in 1:1379){
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

pnts_sf$indicator <- st_within(pnts_sf, mdd) %>% lengths > 0
  
pnts_mdd <- subset(pnts_sf, indicator == TRUE)

pnts_mdd <- pnts_mdd %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2])

pnts_mdd <- st_drop_geometry(pnts_mdd) 

write.csv(pnts_mdd,"passerine_occ_pts_mdd.csv", row.names = TRUE)

#-------------------------------------------------------------#
#visually check resolution of species                         #
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
s <- dismo::gridSample(pnts_mdd[pnts_mdd$scientificName == "Corythopis torquatus", 4:5], r, n=1) #65 obs for focal species

s0 <- dismo::gridSample(pnts_mdd[pnts_mdd$scientificName != "Corythopis torquatus", 4:5], r, n=1) #178 obs for background species

#now create final thinned dataset




