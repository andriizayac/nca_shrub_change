### import libraries
pckgs <- c("raster", "rgdal", "sf", "gplyr", "rgeos")
sapply(pckgs, require, character.only = T)

# import NCA boundary and create a buffer
nca <- readOGR(paste0(path,"raw_data/BDY_NOC_NLCSNMNCA_PUB_24K_POLY/BDY_NOC_NLCSNMNCA_PUB_24K_POLY.shp"))[2,] %>%
  gBuffer(width = 13000) %>%
  spTransform(CRSobj = CRS("+init=epsg:26911")) %>%
  spTransform(CRSobj = CRS("+init=epsg:4326")) %>%
  as("SpatialPolygonsDataFrame")

# fires <- readOGR("raw_data/fireshape/US_Wildfires_1878_2019.shp") %>% 
#   filter(FireYear > 1960) %>%
#   spTransform(CRSobj = CRS("+init=epsg:4326")) %>% 
#   buffer(width = 0, dissolve = FALSE) %>% 
#   crop(nca)

# writeOGR(fires, "data/", "firedata_1961_2019_epsg4326", driver = "ESRI Shapefile")
# writeOGR(nca, "data/", "nca_boundary_epsg4326", driver = "ESRI Shapefile")

