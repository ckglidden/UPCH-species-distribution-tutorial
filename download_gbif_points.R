setwd("~/Documents/GitHub/UPCH-species-distribution-tutorial")

library(sf); library(rgbif); library(dplyr)

#read in IUCN species names for bats in Peru
iucn_species <- read_sf("iucn_peruvian_bats") #this shape file was downloaded from ICUN redlist site
names(iucn_species)[3] <- "scientificName" #needs to match gbif to work

#create a dataframe with the scientific name, and all the occurrences of the sandfly
latlong <- data.frame()

#for loop to get lat long
for (i in 1:299){
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

write.csv(pnts_mdd,"bat_occ_pts_mdd.csv", row.names = TRUE)

###plot test
# Remove plot axis
library(ggplot2)

no_axis <- theme(axis.title=element_blank(),
                 axis.text=element_blank(),
                 axis.ticks=element_blank())

names <- unique(pnts_mdd$scientificName) #look at distribution by species
#look by genus or family?

# Plot "Chiroderma villosum"
ggplot() +
  geom_sf(data=mdd, color="#2D3E50", fill="lightgrey", size=.15, show.legend = FALSE) +
  geom_point(data = subset(pnts_mdd, scientificName == names[79]), 
             aes(x = lon, y = lat), alpha = 0.5) +
  theme_minimal() +
  no_axis

#list to re-check
#"Platyrrhinus brachycephalus"
#"Mesophylla macconnelli"
#"Lonchophylla thomasi"
#"Uroderma magnirostrum"

ggsave("glabrata_all_aquatic_pts_map.tiff", glabrata_points, dpi = 600)



#double check how many unique points in a 1km raster grid...or 30m?