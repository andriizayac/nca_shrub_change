# the raw data represent a crop + mask subset of the NLCD data set 
# (https://data.usgs.gov/datacatalog/data/USGS:5ed816eb82ce7e579c67004a) 
# using the outlines of the Morley Nelson NCA + 13,000 m buffer 
# (accounting for collard raptors flying outside the NCA boundaries)

### import libraries
pckgs <- c("raster", "dplyr", "sf", "lme4", "ggplot2")
sapply(pckgs, require, character.only = T)

### helper functions
vec <- function(x) {return( t(t(as.vector(x))) )} # from matrixcalc

# path to data
path <- "~/../../Volumes/JenCruz/Common/QCLData/Habitat/NCA/Shrub_model/"


# The following chunk only had to run once
### import raw data 
# sage <- brick(list.files("raw_data", full.names = T, pattern = "sage")) 
# shrub <- brick(list.files("raw_data", full.names = T, pattern = "shrub"))

# replace '101' values with NA; reclassify is a memory-friendly function
# sage1 <- reclassify(sage, cbind(100, Inf, NA), right = TRUE)
# shrub1 <- reclassify(shrub, cbind(100, Inf, NA), right = TRUE)

### aggregate and export
# aggregate(shrub1, fact = 5) %>% # 150m or 300m aggregation
#   writeRaster(filename = "data/NCA_1985_2018_buff_shrub.tif", overwrite = T) 
# 

### import aggregated data
nca_150 <- brick(paste0(path, "data/NCA_1985_2018_buff_sage_150m.tif"))
nca_300 <- brick(paste0(path, "data/NCA_1985_2018_buff_sage_300m.tif"))
# --- subset
nca <- brick(paste0(path, "data/NCA_1985_2018_buff_shrub.tif")) %>% 
  crop(extent(.)) 

fires <- st_read(paste0(path, "data/firedata_1961_2019_epsg4326.shp")) %>% 
  st_transform(nca@crs)
ncapol <- st_read(paste0(path, "data/nca_boundary_epsg4326.shp")) %>% 
  st_transform(nca@crs)

nca %>% 
  subset(subset = c(1:5)) %>% 
  mask(fires[fires$FireYear == 1986,]) %>% 
  as.data.frame(na.rm = T, xy = FALSE) %>% 
  data.matrix() %>%
  t() %>% 
  matplot(type = "l", col = rgb(0,0,0,.05))

# === reclassify data into dummy variables: 
# --- the dummy vars are: recovery=101, stasis=102, burn=103
# --- never burnt
stas <- nca %>%
  #mask(fires[fires$FireYear > 1984,], inverse = T) %>% 
  reclassify(rcl = t( matrix(c(-1, 100, 102), nr = 3)) )

recov <- nca %>%
  mask(fires[fires$FireYear > 1984,]) %>% 
  reclassify(rcl = t( matrix(c(-1, 100, 101), nr = 3)) )

sr = mosaic(stas, recov, fun = "max")

# --- first year as NA because we don't know what preceded it
l1 <- nca %>% 
  subset(1) %>% 
  mask(fires[fires$FireYear == 1985, ]) %>% 
  reclassify(rcl = t( matrix(c(-1, 100, 103), nr = 3)) )

# --- reclassify the rest of the layers
for(i in 2:dim(nca)[3] ) {
  l2 <- nca %>% 
    subset(i) %>% 
    mask(fires[fires$FireYear == 1985 + (i-1), ]) %>% 
    reclassify(rcl = t( matrix(c(-1, 100, 103), nr = 3)) )
    
  l1 <- stack(l1, l2)
}

# --- mosaic treatments together
# this output contains marks for stasis and burn events. needs recovery marks that are implemented in the next step
# 101 - recovery; 102 - stasis; 103 - burn
finr <- mosaic(stas, l1, fun = "max")

# === convert cover and treatment rasters into data frames
a = as.data.frame(nca, xy = FALSE, na.rm = TRUE)
trtm <- as.data.frame(finr, xy = TRUE, na.rm = TRUE)
rxy <- select(trtm, x, y)
trtm <- select(trtm, -x, -y)

# this loop fill in recovery indices after disturbance
for(i in 2:ncol(a)) {
  y <- which(trtm[,i-1] == 103)
  trtm[y, i] <- 101
  z <- which(trtm[, i-1] == 101 & trtm[, i] != 103)
  trtm[z, i] <- 101
}
# last, create time since disturbance variable
tsd <- matrix(NA, nr = nrow(a), nc = ncol(a))
for(i in 2:32) {
  y <- which(trtm[, i-1] == 103)
  tsd[y, i] <- 0
  z <- which(is.finite(tsd[, i]))
  tsd[z, i+1] <- tsd[z, i] + 1
}

id <- rep(rownames(a), ncol(a)-1)
Nt0 <- vec( data.matrix(a[, -ncol(a)]) )
Nt1 <- vec( data.matrix(a[, -1]) ) 
tr <- vec( data.matrix(trtm[, -ncol(trtm)]) ) %>% 
  dplyr::recode(`101` = "recovery", `102` = "stasis", `103` = "burn") %>% 
  factor()
timerec <- vec(tsd[, -1])
x <- rep(rxy[,"x"], ncol(a)-1)
y <- rep(rxy[,"y"], ncol(a)-1)

df <- data.frame(ID = id, nt0 = Nt0, nt1 = Nt1, 
                 treatment = tr, timerec = timerec) %>% 
  bind_cols(x = x, y = y)
  # group_by(treatment) %>%
  # filter(nt0 != 0 & nt1 != 0) %>%
  # # sample_frac(.0051) #%>%
  # sample_n(500) %>%
  # ungroup()

# write.csv(df, file = paste0(path, "data/sage_dat.csv"), row.names = FALSE)

df <- read.csv(paste0(path, "data/sage_dat.csv"))  

# grid sample

library(dismo)
xy <- df %>% 
  select(ID, x, y) %>% 
  group_by(ID) %>% 
  summarize_all(mean) %>% as.data.frame()

rownames(xy) <- xy$ID
res <- 1800
r <- raster(extent(range(xy$x), range(xy$y)) + res)
res(r) <- res
s <- gridSample(xy[,2:3], r, n = 25)

s$ID <- as.numeric(rownames(s))

dffin <- filter(df, ID %in% s$ID)

write.csv(dffin, file = paste0(path, "shrub_dat_subset.csv"), row.names = FALSE)
