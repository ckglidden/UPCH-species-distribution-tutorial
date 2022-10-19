setwd("~/Documents/GitHub/UPCH-species-distribution-tutorial")

library(sf); library(rgbif)

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

#clip occurrence points to only ones that fall in MDD