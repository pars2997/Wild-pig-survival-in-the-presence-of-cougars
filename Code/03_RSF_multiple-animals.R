##################################################################
###### Cougar feeding site RSF and risk layer creation ###########
############ Mitchell Parsons ####################################
############ February 22, 2024 ###################################

#### Load packages ####
library(sf)
library(nngeo)
library(terra)
library(amt)
library(tidyverse)
library(broom)
library(lme4)
library(AICcmodavg)

#### Read in data ####

# clusters is each investigated GPS cluster and what we found there
# Read in and filter to only current clusters where a prey item was found
# current excludes very old clusters that were investigated, dens, and day beds
clusters <- read.csv("../RawData/cougar_cluster_investigations.csv")
clusters <- clusters[clusters$type == "current",]
clusters <- clusters[clusters$carcass_found == "Y",]
# correct longitudes that were entered without the negative
clusters[clusters$carcass_lon > 0,"carcass_lon"] <- clusters[clusters$carcass_lon > 0,"carcass_lon"] * -1
# filter to the three cougars we have many kills for
clusters <- clusters[clusters$cougar_id %in% c("SF09","SF01","SM06"),]


#### Data processing ####

# create a simple feature object for plotting, creating home ranges, and 
# extracting covariates. crs 4326 is WGS 84 lat long, 32610 is UTM zone 10N

clusters <- st_as_sf(clusters, coords = c("carcass_lon", "carcass_lat"), crs = 4326)
clusters <- st_transform(clusters, 32610)

# lump prey into broad categories
clusters <- clusters %>% 
  mutate(prey_cat = case_when(
    species == "deer" ~ "Deer",
    species == "pig" ~ "Pig",
    species %in% c("coyote","bobcat","opossum") ~ "Carnivore",
    TRUE ~ "Other"
  ))

# Order prey by broad categories for plotting
clusters$prey_cat <- ordered(clusters$prey_cat,levels = c("Deer","Pig","Carnivore","Other"))

# extract x and y coordiantes as their own values for creating tracks
clusters$x <- st_coordinates(clusters$geometry)[,1]
clusters$y <- st_coordinates(clusters$geometry)[,2]
clusters$month <- month(mdy(clusters$form_date))
clusters <- clusters %>% 
  mutate(season = case_when(
    month %in% c(4,5,6,7) ~ "summer",
    month %in% c(10,11,12,1,2) ~ "winter"
  ))


# create tracks out of the GPS data
lion_track <- track(clusters$x, clusters$y,id = clusters$cougar_id,season = clusters$season)

# Group the track by individual cougars
lion_track_group <- lion_track %>% 
  nest(data = -"id")

#### Home ranges and random points ####

# Read in MCP files for each individual cougar
SF9_hr <- read_sf("../ProcessedData/SF9_HomeRange_MCP.shp")
SF1_hr <- read_sf("../ProcessedData/SF1_HomeRange_MCP.shp")
SM6_hr <- read_sf("../ProcessedData/SM6_HomeRange_MCP.shp")

# Add a 500m buffer around each home range. This is the area from which
# we will draw random points.
SF9_hr <- SF9_hr %>%  
  sf::st_buffer(dist = 500)
SF1_hr <- SF1_hr %>%  
  sf::st_buffer(dist = 500)
SM6_hr <- SM6_hr %>%  
  sf::st_buffer(dist = 500)

# Add homes ranges to our nested data frame
lion_track_group$hrs <- list(SF9_hr,SF1_hr,SM6_hr)

# Create random points for each individual cougar
# creating 1000 random points for each actual kill
# making points within the MCP+500 of each cougar

rand <- list()
for(i in 1:nrow(lion_track_group)){
  rand[[i]] <- random_points(lion_track_group$hrs[[i]], 
                             n = nrow(lion_track_group$data[[i]])*1000,
                             presence = lion_track_group$data[[i]])
}
# add random points to the nested dataframe
lion_track_group$rand <- rand



#### Read in GIS covariates that will be used in the RSF ####
# is.land is used to mask cells that occur in the ocean to avoid unrealistic
# prediction values

is.land <- rast("../ProcessedData/is_land.tiff")
forest_dist <- rast("../ProcessedData/log_forest_dist.tiff")*is.land
grass_dist <- rast("../ProcessedData/log_grass_dist.tiff")*is.land
shrub_dist <- rast("../ProcessedData/log_shrub_dist.tiff")*is.land
riparian_dist <- rast("../ProcessedData/log_riparian_dist.tiff")*is.land
TRI <- rast("../ProcessedData/resamp_TRI.tiff")*is.land
TPI <- rast("../ProcessedData/resamp_TPI.tiff")*is.land
DEM <- rast("../ProcessedData/resamp_DEM.tiff")*is.land


windowsize <- 5
window <- matrix(rep(1,windowsize^2),ncol = windowsize, nrow = windowsize)
TRI <- focal(x = TRI, w = window, fun = mean)

# create aspect from the DEM
# transform aspect into North-South and East-West components
aspect <- terra::terrain(DEM,"aspect")
aspect_rad <- (aspect * pi)/180
aspectNS <- cos(aspect_rad)
aspectEW <- sin(aspect_rad)

# create template raster that covers the needed area
xmin <- as.numeric(round(st_bbox(SM6_hr))[1])
ymin <- as.numeric(round(st_bbox(SM6_hr))[2])-10000
xmax <- as.numeric(round(st_bbox(SM6_hr))[3])+5000
ymax <- as.numeric(round(st_bbox(SM6_hr))[4])+5000

r <- rast(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,nrow = 30, ncol = 30)

# crop all covariates to the needed area to decrease storage size

forest_dist <- crop(forest_dist,r)
grass_dist <- crop(grass_dist,r)
shrub_dist <- crop(shrub_dist,r)
riparian_dist <- crop(riparian_dist,r)
TRI <- crop(TRI,r)
TPI <- crop(TPI,r)
DEM <- crop(DEM,r)
aspectNS <- crop(aspectNS,r)
aspectEW <- crop(aspectEW,r)

# And finally scale all covariates to facilitate model convergence and 
# interpretation of parameters

DEM <- terra::scale(DEM)

# create a raster stack to extract covariates from and rename to short
# yet effective names
allcovs <- c(forest_dist,grass_dist,shrub_dist,
             riparian_dist,TRI,TPI,DEM,aspectNS,aspectEW)
names(allcovs) <- c("forest","grass","shrub","riparian","TRI","TPI",
                    "elev","aspectNS","aspectEW")

# Finally, extract covariate values for each used and available location
lion_track_group <- lion_track_group %>% 
  mutate(covs = lapply(rand,extract_covariates,allcovs))

cor(cbind(lion_track_group$covs[[1]]$forest,lion_track_group$covs[[1]]$shrub,
    lion_track_group$covs[[1]]$grass,lion_track_group$covs[[1]]$riparian))


#### Fit RSF for moutnain lion kill sites
# this fits an individual model for each individual cougar
lion_track_group <- lion_track_group %>% 
  mutate(rsf = lapply(covs, fit_rsf, case_ ~ forest + 
                        riparian + 
                        shrub +
                        TRI +
                        aspectNS*aspectEW))

# extract AICc values for the fit models
aic.models <- lion_track_group %>% 
  # This will create a list column
  mutate(aics = lapply(rsf, function(each_rsf) {
    res <- AICc(each_rsf$model)
    return(res)
  })) %>% 
  dplyr::select(aics)

unlist(aic.models)

# extract parameter estimates for each cougar's rsf
uni_beta <- lion_track_group %>% 
  # This will create a list column
  mutate(betas = lapply(rsf, function(each_rsf) {
    res <- tidy(each_rsf$model)
    return(res)
  })) %>% 
  # Now we want to drop the other columns and unnest
  dplyr::select(id, betas) %>% 
  unnest(cols = betas)
# look at unweighted means and SEs for each parameter estimate
uni_beta %>% 
  group_by(term) %>% 
  summarize(mean = mean(estimate),
            se = sd(estimate)/sqrt(3),
            t = mean/se) 

# add inverse variance for weighting parameters
uni_beta <- uni_beta %>% 
  mutate(inv_var = 1/(std.error^2))

# Loop through all model terms, run a linear model on parameter estimates,
# and calculate weighted means, ses, and p values for each parameter
terms <- unique(uni_beta$term)
param.est <- vector()
se.est <- vector()
t.est <- vector()
p.est <- vector()
ci.95.L <- vector()
ci.95.U <- vector()

for(i in 1:length(terms)){
  temp <- lm(estimate ~ 1,
             data = uni_beta,
             subset = term == terms[i],
             weights = inv_var)
  param.est[i] <- round(summary(temp)$coefficients[1],3)
  se.est[i] <- round(summary(temp)$coefficients[2],3)
  t.est[i] <- round(summary(temp)$coefficients[3],3)
  p.est[i] <- round(summary(temp)$coefficients[4],3)
  ci.95.L[i] <- round(confint(temp)[1],3)
  ci.95.U[i] <- round(confint(temp)[2],3)
  
}

# combine weighted estimates and confidence intervals into a data frame 
# and label columns
weighted.est <- data.frame(terms,param.est,se.est,t.est,p.est,ci.95.L,ci.95.U)
colnames(weighted.est) <- c("Parameter","Estimate","SE","t","p","Lower 95 CI","Upper 95 CI")
weighted.est


# write.csv(weighted.est,file = "lionkillrsf.csv",row.names = F)
# fullaspect <- unlist(aic.models)

#### predict cougar risk based on top model ####
# use parameter estiamtes and raster layers to create a raster of cougar risk
logit.lion_pred_risk <- weighted.est$Estimate[1] + 
  weighted.est$Estimate[2]*forest_dist +
  weighted.est$Estimate[3]*riparian_dist +
  weighted.est$Estimate[4]*shrub_dist +
  weighted.est$Estimate[5]*TPI + 
  weighted.est$Estimate[6]*aspectNS + 
  weighted.est$Estimate[7]*aspectEW + 
  weighted.est$Estimate[8]*aspectNS*aspectEW

#### Calculate backtransformed cougar predation risk ####
lion_pred_risk <- exp(logit.lion_pred_risk)/(1 + exp(logit.lion_pred_risk))
plot(lion_pred_risk)
points(clusters[clusters$cougar_id == "SF09",],col = "darkblue")
points(clusters[clusters$cougar_id == "SF01",],col = "grey20")
points(clusters[clusters$cougar_id == "SM06",],col = "maroon")
lines(SF9_hr, col = "darkblue",lty = 2,lwd = 2)
lines(SF1_hr, col = "grey20",lty = 2,lwd = 2)
lines(SM6_hr, col = "maroon",lty = 2,lwd = 2)

writeRaster(lion_pred_risk, "../ProcessedData/LionPredRisk_logdist.tiff", overwrite = T)