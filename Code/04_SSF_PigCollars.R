#############################################
#### Pig Step Selection Functions ###########
########## Mitchell Parsons##################
############ 10/6/2023#######################

#### load packages ####
library(amt)
library(ggplot2)
library(lubridate)
library(tidyverse)
library(ctmm)
library(terra)
library(sf)
library(suncalc)
library(tidyterra)
library(ggpubr)
library(adehabitatHR)
library(AICcmodavg)
library(broom)

#### Read in location data and do inital processing ####
# Seperate csv file for each pig which will be formatted and combined
FA05 <- read.csv("../RawData/PigLocations/collar20147_FA05_03-31-2021_09-17-2021.csv",header = F)
FA07 <- read.csv("../RawData/PigLocations/collar20147_FA07_2021_12_20-2022_01_08.csv",header = F)
FA01<- read.csv("../RawData/PigLocations/collar20149_FA01_AllDat.csv",header = F)
FA02 <- read.csv("../RawData/PigLocations/collar20150_FA02_2021-02-23_2023-02-21.csv",header = F)
FA09 <- read.csv("../RawData/PigLocations/collar26212_FA09_2022_04_15-2022_06_22.csv",header = F)
FA04 <- read.csv("../RawData/PigLocations/collar26215_FA04_2021_03_31-2021_05_07.csv",header = F)
MS09 <- read.csv("../RawData/PigLocations/collar33346_MS09_2023_06_02-2023_07_17.csv",header = F)
FA08 <- read.csv("../RawData/PigLocations/collar33347_FA08_2021-12-20_2024-01-01.csv",header = F)
FS04 <- read.csv("../RawData/PigLocations/collar33348_FS04_2021_03_31-2021_04_27.csv",header = F)
FS11 <- read.csv("../RawData/PigLocations/collar33349_FS11_2023-05-31_2024-01-01.csv",header = F)
FA10 <- read.csv("../RawData/PigLocations/collar33349_FA10_2023_01_26-2023_03_25.csv",header = F)
MA02 <- read.csv("../RawData/PigLocations/collar33351_MA02_2021_12_20-2022_07_04.csv",header = F)
FA11 <- read.csv("../RawData/PigLocations/collar33352_FA011_2023_01_26-2023_06_23.csv",header = F)
FA06_1 <- read.csv("../RawData/PigLocations/collar33353_FA06_2022-04-15_2023-05-11.csv",header = F)
MA04 <- read.csv("../RawData/PigLocations/collar33354_MA04_2023-01-27_2024-01-01.csv",header = F)
MA01 <- read.csv("../RawData/PigLocations/collar33356_MA01_02-23-2021_11-25-2023.csv",header = F)
FA06_2 <- read.csv("../RawData/PigLocations/collar33358_FA06_2021_03_31-2021_08_20.csv",header = F)
MA03 <- read.csv("../RawData/PigLocations/collar46254_MA03_2023-01-27_2023-04-10.csv",header = F)
FS10 <- read.csv("../RawData/PigLocations/collar46254_FS10_2023-05-31_2024-01-01.csv",header = F)
FA03 <- read.csv("../RawData/PigLocations/collar46255_FA03_2021_03_05-2021_05_30.csv",header = F)

# create a list of all pig collar data
pigcollars <- list(FA01, FA02, FA03,FA04,FA05,FA06_1,FA06_2,FA07,FA08,FA09,FA10,FA11,
                   FS04,FS10,FS11,
                   MA01,MA02,MA03,MA04,
                   MS09)

# Create a function that completes needed data processing
# removes column names from data, selects needed columns, and renames.
# The ensures numeric formating for needed columns
collardat.org <- function(x){
  x <- x[-1,c(1,2,3,4,13,14,15,16,49,50)]
  colnames(x) <- c("LocID","CollarID","UTC_Date","UTC_Time","Latitude","Longitude",
                       "Height","DOP","Easting","Northing")
  x <- x %>% 
    mutate(Latitude = as.numeric(Latitude),
           Longitude = as.numeric(Longitude),
           Height = as.numeric(Height),
           DOP = as.numeric(DOP),
           Easting = as.numeric(Easting)-10000000,
           Northing = as.numeric(Northing))
  return(x)
}
# apply the data processing function to each individual pig
pigcollars.org <- lapply(pigcollars, collardat.org)


# Add animal ID column to each pig's data
pigcollars.org[[1]]$AID <- "FA01"
pigcollars.org[[2]]$AID <- "FA02"
pigcollars.org[[3]]$AID <- "FA03"
pigcollars.org[[4]]$AID <- "FA04"
pigcollars.org[[5]]$AID <- "FA05"
pigcollars.org[[6]]$AID <- "FA06"
pigcollars.org[[7]]$AID <- "FA06"
pigcollars.org[[8]]$AID <- "FA07"
pigcollars.org[[9]]$AID <- "FA08"
pigcollars.org[[10]]$AID <- "FA09"
pigcollars.org[[11]]$AID <- "FA10"
pigcollars.org[[12]]$AID <- "FA11"
pigcollars.org[[13]]$AID <- "FS04"
pigcollars.org[[14]]$AID <- "FS10"
pigcollars.org[[15]]$AID <- "FS11"
pigcollars.org[[16]]$AID <- "MA01"
pigcollars.org[[17]]$AID <- "MA02"
pigcollars.org[[18]]$AID <- "MA03"
pigcollars.org[[19]]$AID <- "MA04"
pigcollars.org[[20]]$AID <- "MS09"


# Read in collar deployment data
# This will be used to filter data to only the times the collar was actively deployed
deployments <- read.csv("../RawData/pig_collar_deployments.csv")

head(deployments)
# rename columns and format dates
colnames(deployments) <- c("deployment_id","CollarID","AID","start_date",
                           "start_day","start_month","start_year",
                           "end_deployment","end_day","end_month","end_year",
                           "end_cause")
deployments$start_date <- mdy(deployments$start_date)
deployments$end_deployment <- mdy(deployments$end_deployment)

#### Combine pig data ####

# loop through each pig's data
# format dates and times
# Match with correct collar deployment data
# filter location data to the dates between collar deployment and 
# the end of the deployment
# Removes locations during collar testing and after recovery
piglocs.filtered <- data.frame()
for(i in 1:length(pigcollars.org)){
  temp <- pigcollars.org[[i]]
  temp$Date <- mdy(temp$UTC_Date)
  temp$Time <- hms(temp$UTC_Time)
  temp$dt <- mdy_hms(paste(temp$UTC_Date,temp$UTC_Time)) - hours(7)
  
  temp2 <- deployments[deployments$CollarID == temp$CollarID[1] & deployments$AID == temp$AID[1],]
  temp <- temp[temp$dt > temp2$start_date + hours(24) &
                 temp$dt < temp2$end_deployment - hours(24),]
  piglocs.filtered <- rbind(piglocs.filtered,temp)

}
# Removed missed fixes
piglocs.filtered <- piglocs.filtered[!is.na(piglocs.filtered$Latitude),]

# Remove locations with DOP >5
piglocs.filtered <- piglocs.filtered[piglocs.filtered$DOP <5,]
# how many locations for each individual?
table(piglocs.filtered$AID)

#### Correct dates and times and add night/day data column ####
# force timezone to pacific for all collars
# add add corrected date and time based on time zone
piglocs.filtered$dt <- force_tz(piglocs.filtered$dt,tzone = "America/Los_Angeles")

piglocs.filtered[is.na(piglocs.filtered$dt),"dt"] <- mdy_hms(paste(piglocs.filtered$UTC_Date[is.na(piglocs.filtered$dt)],
                                                                   piglocs.filtered$UTC_Time[is.na(piglocs.filtered$dt)])) - hours(6)


piglocs.filtered$Date <- date(piglocs.filtered$dt)

# Add sunrise and sunset times based on date and lat long
piglocs.filtered$sunrise <- getSunlightTimes(data = data.frame(date = piglocs.filtered$Date,
                                                               lat = piglocs.filtered$Latitude,
                                                               lon = piglocs.filtered$Longitude),
                                             tz = "America/Los_Angeles")$sunrise

piglocs.filtered$sunset <- getSunlightTimes(data = data.frame(date = piglocs.filtered$Date,
                                                               lat = piglocs.filtered$Latitude,
                                                               lon = piglocs.filtered$Longitude),
                                             tz = "America/Los_Angeles")$sunset

# add column for each location specifying if it occured at night
piglocs.filtered <- piglocs.filtered %>% 
  mutate(night = case_when(
    dt < (sunrise + hours(1)) | dt > (sunset - hours(1)) ~ 1,
    TRUE ~ 0
  ))

head(piglocs.filtered)

# Remove individual data files to clear space
rm(FA01,FA02,FA03,FA04,FA05,FA06_1,FA06_2,FA07,FA08,FA09,FA10,FA11,
   MA01,MA02,MA03,MA04,
   MS09,FS04,FS10,FS11,pigcollars.org,pigcollars,deployments)
gc()

# Jolon Road Polygon
# Check how many collared pigs crossed the road

# Create polygon of area across the road that could be used by pigs
pt1 <- c(35.969149, -121.176401)
pt2 <- c(36.023706, -121.177895)
pt3 <- c(36.017310, -121.048566)
pt4 <- c(35.919891, -121.033356)
pt5 <- c(35.936594, -121.064463)
pt6 <- c(35.939770, -121.128203)
pt7 <- pt1
polygon.df <- as.data.frame(rbind(pt1,pt2,pt3,pt4,pt5,pt6,pt7))

pol = st_polygon(
  list(
    cbind(
      polygon.df$V2[c(1,2,3,4,5,6,7)], 
      polygon.df$V1[c(1,2,3,4,5,6,7)])
  )
)

# Add CRS to the polygon
polc = st_sfc(pol, crs=4326)

# Filter GPS collar locations to only those within the polygon
temp <- st_as_sf(piglocs.filtered,coords = c("Longitude","Latitude"),crs = 4326)
temp2 <- st_filter(temp,polc)
# Check how many pigs used the area across the road and how many points they had there.
length(table(temp2$AID))
length(table(piglocs.filtered$AID))


#### Read in GIS Covariates ####
# is.land filters masks the ocean to prevent weird preditction values
is.land <- rast("../ProcessedData/is_land.tiff")
DEM <- rast("../ProcessedData/resamp_DEM.tiff")*is.land
propshrub <- rast("../ProcessedData/propshrub_25cell.tif")*is.land
propgrass <- rast("../ProcessedData/propgrass_25cell.tif")*is.land
propforest <- rast("../ProcessedData/propforested_25cell.tif")*is.land
propriparian <- rast("../ProcessedData/propriparian_25cell.tif")*is.land
forest_dist <- rast("../ProcessedData/log_forest_dist.tiff")*is.land
grass_dist <- rast("../ProcessedData/log_grass_dist.tiff")*is.land
shrub_dist <- rast("../ProcessedData/log_shrub_dist.tiff")*is.land
riparian_dist <- rast("../ProcessedData/log_riparian_dist.tiff")*is.land
TPI <- rast("../ProcessedData/resamp_TPI.tiff")*is.land
TRI <- rast("../ProcessedData/resamp_TRI.tiff")*is.land
LionRisk <- rast("../ProcessedData/LionPredRisk_logdist.tiff")

# resample all variables to be same extent as LionRisk
DEM <- resample(DEM,LionRisk)
propshrub <- resample(propshrub,LionRisk)
propgrass <- resample(propgrass,LionRisk)
propforest <- resample(propforest,LionRisk)
propriparian <- resample(propriparian,LionRisk)
forest_dist <- resample(forest_dist,LionRisk)
grass_dist <- resample(forest_dist,LionRisk)
shrub_dist <- resample(shrub_dist,LionRisk)
riparian_dist <- resample(riparian_dist,LionRisk)
TPI <- resample(TPI, LionRisk)
TRI <- resample(TRI, LionRisk)

# scale all variables to facilitate model convergence and parameter interp
DEM <- terra::scale(DEM)
LionRisk <- terra::scale(LionRisk)

# Use a moving window to estimate lion predation risk over an area
# that is windowsize x windowsize cells
windowsize <- 5
window <- matrix(rep(1,windowsize^2),ncol = windowsize, nrow = windowsize)


LionRisk_25 <- focal(x = LionRisk, w = window, fun = mean)
# writeRaster(LionRisk_25, filename = "../ProcessedData/LionRisk_25cell.tif")


#### Plots of pig locations/KDE's plotted over relative predation risk ####

# create sf object from pig location data
piglocs <- st_as_sf(piglocs.filtered,coords = c("Easting","Northing"),crs = 32610)

# create objects for two individual pigs, FA01 and FA08
temp <- piglocs[piglocs$AID == "FA08",]
# temp.day <- piglocs[piglocs$AID == "FA08" & piglocs$night == 0,]
# temp.night <- piglocs[piglocs$AID == "FA08" & piglocs$night == 1,]

temp2 <- piglocs[piglocs$AID == "FA01",]
# temp.day2 <- piglocs[piglocs$AID == "FA01" & piglocs$night == 0,]
# temp.night2 <- piglocs[piglocs$AID == "FA01" & piglocs$night == 1,]

# set min and max extent needed for each pig
xmin <- as.numeric(660250)
ymin <- as.numeric(3977250)
xmax <- as.numeric(667250)
ymax <- as.numeric(3984250)

xmin2 <- as.numeric(652000)
ymin2 <- as.numeric(3978500)
xmax2 <- as.numeric(658000)
ymax2 <- as.numeric(3984500)

# Create single column identifying if a location occurred
# during the day or at night
temp <- temp %>%
  mutate(Period = case_when(
    night == 0 ~ "Day",
    night == 1 ~ "Night"
  ))

temp2 <- temp2 %>%
  mutate(Period = case_when(
    night == 0 ~ "Day",
    night == 1 ~ "Night"
  ))

# add columns for x and y coordinates
temp$xcord <- st_coordinates(temp)[,1]
temp$ycord <- st_coordinates(temp)[,2]

# drop the geometry
FA08 <- st_drop_geometry(temp)
# split the dataframe into day and night compoents
# these are used to estimate KDE's for each period
FA08day <- FA08[FA08$night == 0,]
FA08night <- FA08[FA08$night == 1,]
# create spatial points data frame for kernelUD function
FA08_day_spdf <- SpatialPoints(coords = FA08day[,17:18])
FA08_night_spdf <- SpatialPoints(coords = FA08night[,17:18])
# calculate UDs for individual pig
FA08kud_day <- kernelUD(FA08_day_spdf,grid = 500)
FA08kud_night <- kernelUD(FA08_night_spdf,grid = 500)
# extract polygon for 50% UD (i.e., core area)
FA08_50_day <- getverticeshr(FA08kud_day,50)
FA08_50_night <- getverticeshr(FA08kud_night,50)

# Repeat previous steps for second pig
temp2$xcord <- st_coordinates(temp2)[,1]
temp2$ycord <- st_coordinates(temp2)[,2]

FA01 <- st_drop_geometry(temp2)
FA01day <- FA01[FA01$night == 0,]
FA01night <- FA01[FA01$night == 1,]
FA01_day_spdf <- SpatialPoints(coords = FA01day[,17:18])
FA01_night_spdf <- SpatialPoints(coords = FA01night[,17:18])
FA01kud_day <- kernelUD(FA01_day_spdf,grid = 500)
FA01kud_night <- kernelUD(FA01_night_spdf,grid = 500)
FA01_50_day <- getverticeshr(FA01kud_day,50)
FA01_50_night <- getverticeshr(FA01kud_night,50)

# convert UDs from lat long crs to UTM crs to match lion risk layer
FA01_50_day.2 <- st_as_sf(FA01_50_day,crs = 32610)
st_crs(FA01_50_day.2) <- 32610
FA01_50_night.2 <- st_as_sf(FA01_50_night,crs = 32610)
st_crs(FA01_50_night.2) <- 32610

# set colors for plotting
cols.fill = c("Day" = "grey70", "Night" = "grey20")

# plot FA01's uds over the lion risk layer
FA01_UDplot <- ggplot()+
  geom_spatraster(data = LionRisk_25)+
  xlim(xmin2,xmax2)+
  ylim(ymin2,ymax2)+
  scale_fill_viridis_c(option = "magma",direction = -1)+
  labs(fill = "Relative \npredation risk")+
  geom_sf(data = FA01_50_day.2,color = "cyan3",fill = NA,lwd = 0.8)+
  geom_sf(data = FA01_50_night.2, col = "grey20",fill = NA,lwd = 0.8)+
  # scale_color_manual(name = "Time of day",values = cols.fill)+
  theme_bw()+
  theme(text = element_text(size = 18),axis.text.x = element_text(angle=45, hjust=1))+
  theme(legend.position = "none")+
  ggplot2::annotate("text", x = 652400, y = 3984100, label = "A", size = 6)+
  xlab("Latitude")+
  ylab("Longitude")


FA01_UDplot


# convert FA08's uds from lat long crs to utm crs
FA08_50_day.2 <- st_as_sf(FA08_50_day,crs = 32610)
st_crs(FA08_50_day.2) <- 32610
FA08_50_night.2 <- st_as_sf(FA08_50_night,crs = 32610)
st_crs(FA08_50_night.2) <- 32610

# plot FA08's UDs over the lion risk layer
FA08_UDplot <- ggplot()+
  geom_spatraster(data = LionRisk_25)+
  xlim(xmin,xmax)+
  ylim(ymin,ymax)+
  scale_fill_viridis_c(option = "magma",direction = -1)+
  labs(fill = "Relative \npredation risk")+
  geom_sf(data = FA08_50_day.2,aes(color = "Day"),fill = NA,lwd = 0.8,show.legend = "line")+
  geom_sf(data = FA08_50_night.2, aes(color = "Night"),fill = NA,lwd = 0.8,show.legend = "line")+
  scale_color_manual(name = "Time of day",values = cols.fill)+
  theme_bw()+
  theme(text = element_text(size = 18),axis.text.x = element_text(angle=45, hjust=1))+
  ggplot2::annotate("text", x = 660600, y = 3983850, label = "B", size = 6)+
  xlab("Latitude")+
  ylab("")

FA08_UDplot

# plot both UD layers together
bothUDs <- ggarrange(plotlist = list(FA01_UDplot,FA08_UDplot),ncol = 2,widths = c(1,1.45))

# add save plot to a png for publication
png(file = "../Figures/Pig_UDplots_25cell.png",height = 6, width = 10,
    res = 300, units = "in",type = "cairo")
bothUDs
dev.off()


Both_UDplot <- ggplot()+
  geom_spatraster(data = LionRisk_25)+
  scale_fill_viridis_c(option = "mako",direction = -1)+
  labs(fill = "Relative probability of \ncougar kill occurrence\n")+
  xlim(xmin2,xmax)+
  ylim(ymin,ymax2)+
  # geom_sf(data = FA01_50_day.2,aes(color = "Day"),fill = NA,lwd = 1.5)+
  # geom_sf(data = FA01_50_night.2, aes(color = "Night"),fill = NA,lwd = 1.5)+
  # geom_sf(data = FA08_50_day.2,aes(color = "Day"),fill = NA,lwd = 1.5,show.legend = "line")+
  # geom_sf(data = FA08_50_night.2, aes(color = "Night"),fill = NA,lwd = 1.5,show.legend = "line")+
  coord_sf(datum = st_crs(FA01_50_day.2))+
  theme_bw()+
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle=45, hjust=1),
        legend.title = element_text(size = 15))+
  # ggplot2::annotate("text", x = 652400, y = 3984100, label = "A", size = 6)+
  xlab("Easting")+
  ylab("Northing")+
  scale_color_manual(name = "Time of day",values = cols.fill)

tiff(filename = "../Figures/lionriskbare.tiff",width = 9, height = 5,
     res = 300, compression = "lzw", units = "in")
Both_UDplot
dev.off()


# 
# 
# 
# # the next two plot plot the individual points over the lion risk layer
# # but I decided it was too messy and difficult to interpret
# FA08_plot <- ggplot()+
#   geom_spatraster(data = LionRisk)+
#   xlim(xmin,xmax)+
#   ylim(ymin,ymax)+
#   scale_fill_viridis_c(option = "magma",direction = -1)+
#   labs(fill = "Relative \npredation risk")+
#   geom_sf(data = temp,aes(col = factor(Period)),
#           alpha = 0.6,show.legend = FALSE,size = 1)+
#   scale_color_manual(values = c("cyan3","grey20"))+
#   theme_bw()+
#   facet_wrap(~Period)+
#   theme(text = element_text(size = 18),axis.text.x = element_text(angle=45, hjust=1))
# 
# 
# 
# FA01_plot <- ggplot()+
#   geom_spatraster(data = LionRisk)+
#   xlim(xmin2,xmax2)+
#   ylim(ymin2,ymax2)+
#   scale_fill_viridis_c(option = "magma",direction = -1)+
#   labs(fill = "Relative \npredation risk")+
#   geom_sf(data = temp2,aes(col = factor(Period)),
#           alpha = 0.4,show.legend = FALSE,size = 1)+
#   scale_color_manual(values = c("cyan3","grey20"))+
#   theme_bw()+
#   facet_wrap(~Period)+
#   theme(text = element_text(size = 18),axis.text.x = element_text(angle=45, hjust=1))
# 
# 
# 
# png(file = "../Figures/FA01_UD_plot.png",height = 6, width = 10,
#      res = 300, units = "in",type = "cairo")
# FA01_UDplot
# dev.off()
# 
# png(file = "../Figures/FA08_UD_plot.png",height = 6, width = 10,
#     res = 300, units = "in",type = "cairo")
# FA08_UDplot
# dev.off()

#### Fit iSSA models for pig steps ####

# combine all rasters into a raster stack and name accordingly
allcovs <- c(DEM,propshrub,propgrass,propforest,propriparian,
             forest_dist,grass_dist,shrub_dist,riparian_dist,
             TPI,TRI,LionRisk_25)
names(allcovs) <- c("Elev","PropShrub","PropGrass","PropForest","PropRiparian",
                    "Forest_Dist","Grass_Dist","Shrub_Dist","Riparian_Dist",
                    "TPI","TRI","Lion_Risk")
# remove individual raster layers to save space
rm(is.land,DEM,propshrub,propgrass,propforest,propriparian,
   forest_dist,grass_dist,shrub_dist,riparian_dist,
   TPI,TRI,LionRisk)

gc()

# check correlations of veg cover
# shrub.vals <- values(propshrub,na.rm = T)
# forest.vals <- values(propforest,na.rm = T)
# grass.vals <- values(propgrass,na.rm = T)
# riparian.vals <- values(propriparian,na.rm = T)
# 
# cor(shrub.vals,forest.vals)
# cor(shrub.vals,grass.vals)
# cor(shrub.vals,riparian.vals)
# cor(forest.vals,grass.vals)
# cor(forest.vals,riparian.vals)
# cor(grass.vals,riparian.vals)
# 
# no shrub and grass, no forest and grass

# Run iSSA
# This series of pigs
# 1 - creates tracks and nests by individual pig
# 2 - calcualte pigs step lengths and turn angles based on traks
# I made 0 length steps 0.5m because other wise it causes errors with turn angle
# 3 - creates random steps to pair with each observed step (15 random per obs)
# 4 - extracts covariates for both the start and end of each step
# 5 - runs an iSSA without the lion risk layer
# 6 - runs an iSSA with the lion risk layer


# prop shrub and grass converge but are correlated
# prop grass and riparian work
# prop shrub and riparian work
# prop forest doesn't work

pig_steps <- piglocs.filtered %>% 
  # nest data by pig ID
  nest(data = -"AID") %>% 
  # create tracks for each pig
  mutate(tracks = lapply(data, function(d){
    res <- make_track(d,.x = Easting, .y = Northing, .t = dt, crs = 32610)
    return(res)
  })) %>% 
  # calculate steps and correct for 0-length steps
  mutate(pig.stps = lapply(tracks, function(d){
    res <- steps(d)
    res$ta_[is.na(res$ta_)] <- 0.001
    res$sl_[res$sl_ == 0] <- 0.5
    res$direction_p[is.na(res$direction_p)] <- 0.001
    return(res)
  })) %>%
  # add 15 random steps for each observed step from observed distribution
  mutate(rand_stps = lapply(pig.stps,function(d){
    res <- random_steps(d, n_control = 15,
                   sl_distr = fit_distr((d$sl_ + 0.5), "gamma"),
                   ta_distr = fit_distr(d$ta_[!is.na(d$ta_)],"vonmises"))
    res <- extract_covariates(x = res,covariates = allcovs, where = "both")%>% 
      mutate(log_sl_ = log(sl_),
             cos_ta_ = cos(ta_),
             sunrise = getSunlightTimes(date = date(t1_),
                                        lat = piglocs.filtered$Latitude[1],
                                        lon = piglocs.filtered$Longitude[1],
                                         tz = "America/Los_Angeles")$sunrise,
             
             sunset = getSunlightTimes(date = date(t1_),
                                       lat = piglocs.filtered$Latitude[1],
                                       lon = piglocs.filtered$Longitude[1],
                                       tz = "America/Los_Angeles")$sunset,
             night = case_when(
               t1_ < (sunrise + hours(1)) | t1_ > (sunset - hours(1)) ~ 1,
               TRUE ~ 0)
             )
    return(res)
  })) %>% 
  # fit issa without lion risk
  mutate(issf_norisk = lapply(rand_stps, function(each_ind) {
    res <- each_ind %>% 
      fit_clogit(case_ ~ PropRiparian_end +
                   PropShrub_end +
                   Riparian_Dist_end + 
                   PropRiparian_end:night +
                   PropShrub_end:night +
                   Riparian_Dist_end:night + 
                   sl_ + log_sl_ + cos_ta_ +
                   log_sl_:PropRiparian_start +
                   log_sl_:PropShrub_start +
                   log_sl_:Riparian_Dist_start +
                   log_sl_:night +
                   cos_ta_:PropRiparian_start +
                   cos_ta_:PropShrub_start +
                   cos_ta_:Riparian_Dist_start +
                   cos_ta_:night +
                   strata(step_id_), 
                 model = TRUE)
    return(res)
  })) %>% 
  # fit issa with lion risk
  mutate(issf_risk = lapply(rand_stps, function(each_ind) {
    res <- each_ind %>% 
      fit_clogit(case_ ~ PropRiparian_end +
                   PropShrub_end +
                   Riparian_Dist_end + 
                   Lion_Risk_end +
                   PropRiparian_end:night +
                   PropShrub_end:night +
                   Riparian_Dist_end:night + 
                   Lion_Risk_end:night +
                   sl_ + log_sl_ + cos_ta_ +
                   log_sl_:PropRiparian_start +
                   log_sl_:PropShrub_start +
                   log_sl_:Riparian_Dist_start +
                   log_sl_:Lion_Risk_start +
                   log_sl_:night +
                   cos_ta_:PropRiparian_start +
                   cos_ta_:PropShrub_start +
                   cos_ta_:Riparian_Dist_start +
                   cos_ta_:Lion_Risk_start +
                   cos_ta_:night +
                   strata(step_id_), 
                 model = TRUE)
    return(res)
  }))

#### Review iSSA results ####

# extract parameter estimates for each individual pig's issa
issa_df_risk <- lapply(pig_steps$issf_risk, function(x){
  broom::tidy(x$model)
}) %>% 
  bind_rows()

# We could easily take the mean of each beta now
issa_df_risk %>% 
  group_by(term) %>% 
  summarize(beta = mean(estimate))

# for each term, calculate the weighted mean and standard error
# using inverse variance weighted regression
terms <- unique(issa_df_risk$term)
coefs <- vector()
SEs <- vector()
t <- vector()
p <- vector()
CI.low <- vector()
CI.high <- vector()

for(i in 1:length(terms)){
  temp <- lm(estimate ~ 1, 
             data = issa_df_risk,
             subset = term == terms[i],
             weights = 1/(std.error^2))
  coefs[i] <- summary(temp)$coefficients[1]
  SEs[i] <- summary(temp)$coefficients[2]
  t[i] <- summary(temp)$coefficients[3]
  p[i] <- summary(temp)$coefficients[4]
  CI.low[i] <- confint(temp)[1]
  CI.high[i] <- confint(temp)[2]
}
# combine estimates into a dataframe
mean_estimates_risk <- data.frame(term = terms,
                             coef = round(coefs,3),
                             SE = round(SEs,3),
                             t = round(t,3),
                             p = round(p,3),
                             ci.low = round(CI.low,3),
                             ci.high = round(CI.high,3))
mean_estimates_risk
# write.csv(mean_estimates_risk,file="pigssf_risk.csv", row.names = F)

# calculate AICc values for models iwth lion risk
AIC.risk <- vector()
for(i in 1:nrow(pig_steps)){
  temp <- pig_steps$issf_risk[[i]]$model
  AIC.risk[i] <- AICc(temp)
}

# write.csv(AIC.risk, file = "../ProcessedData/AICcRiskIncluded.csv")

# repeat the above process for the model without lion risk
issa_df_norisk <- lapply(pig_steps$issf_norisk, function(x){
  broom::tidy(x$model)
}) %>% 
  bind_rows()

# We could easily take the mean of each beta now
issa_df_norisk %>% 
  group_by(term) %>% 
  summarize(beta = mean(estimate))

terms <- unique(issa_df_norisk$term)
coefs <- vector()
SEs <- vector()
t <- vector()
p <- vector()
CI.low <- vector()
CI.high <- vector()

for(i in 1:length(terms)){
  temp <- lm(estimate ~ 1, 
             data = issa_df_norisk,
             subset = term == terms[i],
             weights = 1/(std.error^2))
  coefs[i] <- summary(temp)$coefficients[1]
  SEs[i] <- summary(temp)$coefficients[2]
  t[i] <- summary(temp)$coefficients[3]
  p[i] <- summary(temp)$coefficients[4]
  CI.low[i] <- confint(temp)[1]
  CI.high[i] <- confint(temp)[2]
}

mean_estimates_norisk <- data.frame(term = terms,
                                  coef = round(coefs,3),
                                  SE = round(SEs,3),
                                  t = round(t,3),
                                  p = round(p,3))
mean_estimates_norisk
# write.csv(mean_estimates_norisk,file="pigssf_norisk.csv", row.names = F)

# calculate AICc for each pig without risk
AIC.norisk <- vector()
for(i in 1:nrow(pig_steps)){
  temp <- pig_steps$issf_norisk[[i]]$model
  AIC.norisk[i] <- AICc(temp)
}

# write.csv(AIC.norisk, file = "../ProcessedData/AICcRiskExcluded.csv")

# See how AICc scores compare for models with and without risk
sum(AIC.norisk - AIC.risk > 0)
# adding risk improved AICc scores for 18 of 19 pigs
mean(AIC.norisk - AIC.risk)
# average deltaAIC was 39.4 for all individual pigs
AIC.norisk - AIC.risk
# the one that was not improved was only minorly worse with risk included
mean_estimates_norisk
mean_estimates_risk



#### To look at functional response: haven't dont this yet ####

uni_beta2 <- pig_steps %>% 
  # This will create a list column
  mutate(betas = lapply(issf_risk, function(each_rsf) {
    res <- tidy(each_rsf$model)
    return(res)
  })) %>% 
  # Also get the availability
  mutate(avail = lapply(rand_stps, function(each_df) {
    res <- each_df %>% 
      # Keep only the available points
      filter(case_ == FALSE) %>% 
      # Take the average of elevation
      summarize(ripdist_mean = mean(Riparian_Dist_end,na.rm=T),
                ripprop_mean = mean(PropRiparian_end,na.rm=T),
                shrub_mean = mean(PropShrub_end,na.rm=T),
                lionrisk_mean = mean(Lion_Risk_end),na.rm=T)
    return(res)
  })) %>% 
  # Now we want to drop the other columns and unnest
  dplyr::select(AID, betas, avail) %>% 
  unnest(cols = betas) %>% 
  unnest(cols = avail)

# What did we get?
uni_beta2
uni_beta2$inv_var <- 1/(uni_beta2$std.error^2)

# Let's get rid of the other terms to make this cleaner
risk_betas_night <- uni_beta2 %>% 
  filter(term == "Lion_Risk_end:night")
risk_betas_day <- uni_beta2 %>% 
  filter(term == "Lion_Risk_end")
# Quick exploratory plot. Does the estimate vary with availability?

risk_betas_night$night <- "Night"
risk_betas_day$night <- "Day"
risk_betas_all <- rbind(risk_betas,risk_betas_day)

sum(risk_betas_night$estimate >0)
sum(risk_betas_night$estimate <0)
sum(risk_betas_day$estimate >0)
sum(risk_betas_day$estimate <0)

rectdims <- data.frame(xmi = xmin,
                       xma = xmax,
                       ymi = ymin,
                       yma = ymax)
xmin = c(0,2.5,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5)
xmax = c(1.5,3.5,5.5,7.5,9.5,11.5,13.5,15.5,17.5,19.5)
ymin = -10
ymax = 10

risk_betas_all <-  risk_betas_all[order(risk_betas_all$night, risk_betas_all$estimate), ]

ggplot(risk_betas_all, aes(x = rep(1:19,times = 2), y = estimate, col = night)) + 
  geom_point(size = 2.5,position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),width = 0.25,
                position = position_dodge(width = 0.5))+
  scale_color_manual(name = "Time of day", values = c("cyan3","grey20"))+
  xlab("Index") +
  ylab(expression("Estimated " * beta)) +
  theme_bw()+
  geom_hline(yintercept = 0)+
  # facet_wrap(vars(night))+
  # theme(strip.background = element_blank(), strip.text = element_blank())+
  theme(axis.title = element_text(size = 20))+
  theme(axis.text = element_text(size = 16))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 18))+
  scale_x_continuous(breaks = 1:19)
  # theme(axis.text.x = element_blank())
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  # geom_hline(yintercept = mean_estimates_risk[mean_estimates_risk$term == "Lion_Risk_end","coef"],
  #            col = "grey",lty = 2) +
  # geom_hline(yintercept = mean_estimates_risk[mean_estimates_risk$term == "Lion_Risk_end:night","coef"],
  #            col = "grey",lty = 2)+

# Looks like there could be a relationship there. Let's fit a model using
# inverse-variance weighting.

# Calculate inverse variance
rip_betas$inv_var <- 1/(rip_betas$std.error^2)

# Fit
m_rip_avail <- lm(estimate ~ lionrisk_mean + ripdist_mean + ripprop_mean + shrub_mean,
                   data = uni_beta2[uni_beta2$term == "Lion_Risk_end:night",],
                   weights = inv_var)

summary(m_rip_avail)
uni_beta2
cor(data.frame(proprip_day = uni_beta2$estimate[uni_beta2$term == "PropRiparian_end"],
               propshrub_day = uni_beta2$estimate[uni_beta2$term == "PropShrub_end"],
               ripdist_day = uni_beta2$estimate[uni_beta2$term == "Riparian_Dist_end"],
               risk_day = uni_beta2$estimate[uni_beta2$term == "Lion_Risk_end"],
               proprip_night = uni_beta2$estimate[uni_beta2$term == "PropRiparian_end:night"],
               propshrub_night = uni_beta2$estimate[uni_beta2$term == "PropShrub_end:night"],
               ripdist_night = uni_beta2$estimate[uni_beta2$term == "Riparian_Dist_end:night"],
               risk_night = uni_beta2$estimate[uni_beta2$term == "Lion_Risk_end:night"]))
cor(pig_steps$rand_stps[[1]][25:35])

uni_beta2 <- pig_steps %>% 
  # This will create a list column
  mutate(betas = lapply(issf_risk, function(each_rsf) {
    res <- tidy(each_rsf$model, conf.int = T)
    return(res)
  })) %>% 
  # Also get the availability
  mutate(avail = lapply(rand_stps, function(each_df) {
    res <- each_df %>% 
      # Keep only the available points
      filter(case_ == FALSE) %>% 
      # Take the average of elevation
      summarize(ripdist_mean = mean(Riparian_Dist_end),
                ripprop_mean = mean(PropRiparian_end),
                shrub_mean = mean(PropShrub_end),
                lionrisk_mean = mean(Lion_Risk_end))
    return(res)
  })) %>% 
  # Now we want to drop the other columns and unnest
  dplyr::select(AID, betas, avail) %>% 
  unnest(cols = betas) %>% 
  unnest(cols = avail)
