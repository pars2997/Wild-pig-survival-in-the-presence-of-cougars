#### Lion Home Range Estimation ####
## Mitchell Parsons
## 10/6/2023


#### install packages ####

library(amt)
library(ggplot2)
library(lubridate)
library(tidyverse)
library(ctmm)
library(terra)
library(sf)

#### Read in location data ####
# separate file for each collar from SF9 and SM6

sf9_1 <- read.csv("../RawData/CougarLocations/collar44744_SF09_02-10-2021_12-21-2022.csv",header = F)
sf9_2 <- read.csv("../RawData/CougarLocations/collar44746_SF09_12-21-2022_12-31-2023.csv",header = F)
sm6_1 <- read.csv("../RawData/CougarLocations/collar44747_SM06_02-09-2021_12-18-2022.csv",header = F)
sm6_2 <- read.csv("../RawData/CougarLocations/collar45771_SM06_12-18-2022_12-31-2023.csv",header = F)

#### Select only needed columns, rename, and change data types as needed ####
# SF9 - collar 44744
head(sf9_1)
sf9_1 <- sf9_1[-1,c(1,2,3,4,13,14,15,16,49,50)]
colnames(sf9_1) <- c("LocID","CollarID","UTC_Date","UTC_Time","Latitude","Longitude",
                     "Height","DOP","Easting","Northing")
# set needed columns as numeric
sf9_1 <- sf9_1 %>% 
  mutate(Latitude = as.numeric(Latitude),
         Longitude = as.numeric(Longitude),
         Height = as.numeric(Height),
         DOP = as.numeric(DOP),
         Easting = as.numeric(Easting)-10000000,
         Northing = as.numeric(Northing))



# SF9 - collar 44746
head(sf9_2)
sf9_2 <- sf9_2[-1,c(1,2,3,4,13,14,15,16,49,50)]
colnames(sf9_2) <- c("LocID","CollarID","UTC_Date","UTC_Time","Latitude","Longitude",
                     "Height","DOP","Easting","Northing")
sf9_2 <- sf9_2 %>% 
  mutate(Latitude = as.numeric(Latitude),
         Longitude = as.numeric(Longitude),
         Height = as.numeric(Height),
         DOP = as.numeric(DOP),
         Easting = as.numeric(Easting)-10000000,
         Northing = as.numeric(Northing))

# SM6 - collar 44747
head(sm6_1)
sm6_1 <- sm6_1[-1,c(1,2,3,4,13,14,15,16,49,50)]
colnames(sm6_1) <- c("LocID","CollarID","UTC_Date","UTC_Time","Latitude","Longitude",
                     "Height","DOP","Easting","Northing")
sm6_1 <- sm6_1 %>% 
  mutate(Latitude = as.numeric(Latitude),
         Longitude = as.numeric(Longitude),
         Height = as.numeric(Height),
         DOP = as.numeric(DOP),
         Easting = as.numeric(Easting)-10000000,
         Northing = as.numeric(Northing))

# SM6 - collar 45771
head(sm6_2)
sm6_2 <- sm6_2[-1,c(1,2,3,4,13,14,15,16,49,50)]
colnames(sm6_2) <- c("LocID","CollarID","UTC_Date","UTC_Time","Latitude","Longitude",
                     "Height","DOP","Easting","Northing")
sm6_2 <- sm6_2 %>% 
  mutate(Latitude = as.numeric(Latitude),
         Longitude = as.numeric(Longitude),
         Height = as.numeric(Height),
         DOP = as.numeric(DOP),
         Easting = as.numeric(Easting)-10000000,
         Northing = as.numeric(Northing))

#rbind sf9 and sm6 collars into single data frame for each animal
sf9_dat <- rbind(sf9_1,sf9_2)
sm6_dat <- rbind(sm6_1,sm6_2)
rm(sf9_1,sf9_2,sm6_1,sm6_2)

# repeat for sf1 who only had 1 collar
# SF1 - collar 44748
sf1_dat <- read.csv("../RawData/CougarLocations/collar44748_SF1_02-05-2021_11-26-2021.csv",header = F)
sf1_dat <- sf1_dat[-1,c(1,2,3,4,13,14,15,16,49,50)]
colnames(sf1_dat) <- c("LocID","CollarID","UTC_Date","UTC_Time","Latitude","Longitude",
                     "Height","DOP","Easting","Northing")
sf1_dat <- sf1_dat %>% 
  mutate(Latitude = as.numeric(Latitude),
         Longitude = as.numeric(Longitude),
         Height = as.numeric(Height),
         DOP = as.numeric(DOP),
         Easting = as.numeric(Easting)-10000000,
         Northing = as.numeric(Northing))

# add animal ID column to each data frame
sf1_dat$AID <- "SF1"
sf9_dat$AID <- "SF9"
sm6_dat$AID <- "SM6"

# combine all locations into a single data frame
locs <- rbind(sf1_dat,sf9_dat,sm6_dat)


# set date and time and create time stamp
locs$Date <- mdy(locs$UTC_Date)
locs$Time <- hms(locs$UTC_Time)
head(locs)
locs$dt <- mdy_hms(paste(locs$UTC_Date,locs$UTC_Time)) - hours(7)

head(locs)

# reseparate into individual columns because I forgot to do something before
# filter to only those locations that occurred while the collar was deployed
# do this for each collar, then recombine
collar44748 <- locs[locs$CollarID == 44748,]
collar44748 <- collar44748[collar44748$dt > mdy("02/05/2021") & 
                             collar44748$dt < mdy("11/25/2021"),]

collar44744 <- locs[locs$CollarID == 44744,]
collar44744 <- collar44744[collar44744$dt > mdy("2/10/2021") & 
                             collar44744$dt < mdy("12/21/2022"),]

collar44746 <- locs[locs$CollarID == 44746,]
collar44746 <- collar44746[collar44746$dt > mdy("12/21/2022") & 
                             collar44746$dt < mdy("1/1/2024"),]

collar44747 <- locs[locs$CollarID == 44747,]
collar44747 <- collar44747[collar44747$dt > mdy("2/9/2021") & 
                             collar44747$dt < mdy("12/18/2022"),]

collar45771 <- locs[locs$CollarID == 45771,]
collar45771 <- collar45771[collar45771$dt > mdy("12/18/2022") & 
                             collar45771$dt < mdy("1/1/2024"),]

locs <- rbind(collar44748,
              collar44744,
              collar44746,
              collar44747,
              collar45771)
# remove missed fixes
locs <- locs[!is.na(locs$Latitude),]

#### Calculate home ranges for each cougar ####
#### Create tracks ####

lion_track <- make_track(locs, Easting, Northing,
                    dt, id = AID, crs = 32610)
# make a template raster for making KDEs
trast <- make_trast(lion_track, res = 100)

#### Calculate MCPs for each individual lion ####
# these will be used for kill site RSFs to define available points
IDs <- unique(locs$AID)
polyMCP <- list()
hrs <- list()

# loop through individaul lions and calculate mcp
for(i in 1:length(IDs)){
  temp <- locs[locs$AID == IDs[i],]
  # make track
  lion_track <- make_track(temp, Easting, Northing,
                           dt, id = AID, crs = 32610)
  # make mcp
  lion_mcp <- hr_mcp(lion_track, levels = 1)
  # save to list
  hrs[[i]] <- lion_mcp
  polyMCP[[i]] <- lion_mcp$mcp
  print(i)
}

# Plot mcps over a raster of shrubdist to make sure things look right
shrub <- rast("../ProcessedData/log_shrub_dist.tiff")
plot(shrub)
lines(polyMCP[[1]],col = "red",lwd = 2)
lines(polyMCP[[2]],col = "blue",lwd = 2)
lines(polyMCP[[3]],col = "black",lwd = 2)

legend(x = "topleft", legend = IDs, col = c("red","blue","black"),
       pch = 19)


# save each mcp as a shapefile to use in future analyses
write_sf(polyMCP[[1]],"../ProcessedData/SF1_HomeRange_MCP.shp")
write_sf(polyMCP[[2]],"../ProcessedData/SF9_HomeRange_MCP.shp")
write_sf(polyMCP[[3]],"../ProcessedData/SM6_HomeRange_MCP.shp")
