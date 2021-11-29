pkgs <- c("lme4", "dplyr", "raster")
sapply(pkgs, require, character.only = TRUE)

path <- "~/../../Volumes/JenCruz/Common/QCLData/Habitat/NCA/Shrub_model/"

# Code example from Trevor:
# === convert time since fire from numeric to categories of time since fire
years_since_fire <- c(-1:31)
fire_list <- vector("list", 15)
for(i in 1:15) {
  
  char_vec <- rep(NA, times = length(years_since_fire))
  
  for(j in 1:length(years_since_fire)){ 
    char_vec[j] <- ifelse(years_since_fire[j] >= i, 
                          paste(">", i, "years post fire"), 
                          as.character(j))
  }
  fire_list[[i]] <- char_vec
}


# === Sagebrush models
dfs <- read.csv(paste0(path, "data/sage_dat_subset.csv")) %>% 
  mutate(timerec = ifelse(is.na(timerec), 50, timerec))
# === create categorical predictors to the df
K <- 15
catmat <- as.data.frame(matrix(NA, nr = nrow(dfs), nc = K))
for(i in 1:K) {
  catmat[, i] <- ifelse(dfs$timerec > i, 
                        paste0(">", i, "years"),
                        as.character(dfs$timerec))
}

names(catmat) <- paste0("timerec_", 1:K)
dfout <- dfs %>% 
  bind_cols(catmat) %>% 
  mutate(across(starts_with("timerec_"), as.factor))

# write.csv(dfout, file = paste0(path, "data/sage_dat_subset_tsd.csv"), row.names = FALSE)

# === run model loop and store AIC
dfout <- read.csv(paste0(path, "data/sage_dat_subset_tsd.csv"))

mouts <- list()
for(i in 1:K){
  m <- dfout %>% 
    mutate(tsd = dfout[,7+i]) %>% 
    filter(nt0 > 0) %>% 
    glm(nt1 ~ 1 + log(nt0) + tsd, offset = log(nt0), family = poisson(), data = .)
  
  mouts[[i]] <- m
}

# generate AIC table for model comparison
AICcmodavg::aictab(mouts)
# save output
saveRDS(mouts, file = paste0(path, "sage_models.rds"))
# ===

# === Sagebrush models
dfs <- read.csv(paste0(path, "data/shrub_dat_subset.csv")) %>% 
  mutate(timerec = ifelse(is.na(timerec), 50, timerec))


