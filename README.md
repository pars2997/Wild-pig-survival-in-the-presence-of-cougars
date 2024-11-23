## README

Metadata for the manuscript 'Pig Survival and causes of mortality of introduced wild pigs in the presence of cougars' published in Biological Invasions by MA Parsons, KC Vercauteren, JA Dellinger, and JK Young.

The data and R code provided here can be used to recreate results from this manuscript. The data files included several raster layers of topographic and landcover variables, capture and survival monitoring data for collared wild pigs, and a probability of cougar kill occurrence raster.

The R scripts conduct all analyses including survival analysis of wild pigs, a resource selection function of cougar kill sites, and a step selection function for wild pigs.

### DATA FILES

The datafiles include:

1. Raster layers of landcover and topographic data
2. Capture and monitoring data for wild pigs
3. Cougar kill site locations and data

##### log_forest_dist.tiff

This raster layer is a 30x30m resolution raster of the log transformed distance to forest cover. This layer was derived from the 2021 National Landcover Database (https://www.mrlc.gov/data; classes 41, 42, and 43).

##### log_grass_dist.tiff

This raster layer is a 30x30m resolution raster of the log transformed distance to herbaceous cover. This layer was derived from the 2021 National Landcover Database (classes 71, 72, 73, 74, 81, 82).

##### log_riparian_dist.tiff

This raster layer is a 30x30m resolution raster of the log transformed distance to riparian cover. This layer was derived from the 2021 National Landcover Database (classes 11, 12, 90, 95).

##### log_shrub_dist.tiff

This raster layer is a 30x30m resolution raster of the log transformed distance to shrub cover. This layer was derived from the 2021 National Landcover Database (classes 51, 52).

##### propforested_25cell.tif

This raster layer is a 30x30m resolution raster of proportion of forest cover in a 5x5 cell moving window. This layer was derived from the 2021 National Landcover Database (classes 41, 42 and 43).

##### propgrass_25cell.tif

This raster layer is a 30x30m resolution raster of proportion of herbaceous cover in a 5x5 cell moving window. This layer was derived from the 2021 National Landcover Database (classes 71, 72, 73, 74, 81, 82).

##### propriparian_25cell.tif

This raster layer is a 30x30m resolution raster of proportion of riparian cover in a 5x5 cell moving window. This layer was derived from the 2021 National Landcover Database (classes 11, 12, 90, 95).

##### propshrub_25cell.tif

This raster layer is a 30x30m resolution raster of proportion of shrub cover in a 5x5 cell moving window. This layer was derived from the 2021 National Landcover Database (classes 51, 52).

##### is_land.tiff

This raster layer is a 30x30m resolution raster of whether a cell is terrestial land or open water. This layer was derived from the 2021 National Landcover Database (all classes except open water) and was used to mask ocean areas from modeled habitats.

##### resampled_DEM.tiff

This raster layer is a 30x30m resolution raster of elevation. This layer was derived from the USGS 1/3 arc second digital elevation model (https://apps.nationalmap.gov/downloader/)

##### resampled_TRI.tiff

This raster layer is a 30x30m resolution raster of the terrain ruggedness index calculated from resampled_DEM.tiff using the terra package in R.

##### resampled_TPI.tiff

This raster layer is a 30x30m resolution raster of the topographic position index calculated from resampled_DEM.tiff using the terra package in R.

##### LionPredRisk_logdist.tiff

This raster layer is a 30x30m resolution raster of the probability of cougar kill occurrence based on cougar kill site resource selection fucntion modeling. We derived this layer from observed cougar kill sites, but cannot share the raw data because of it's sensitive nature.


##### pig_captures_MJ23_update.csv

This csv file contains the capture information for each indiviudal pig and is used in the survival analysis. Each row represents a capture event and columns provide information on the event and the individual pig. There are 26 columns:

1. pig_id - the ID of the individual pig
2. date - the date of the capture event
3. day - the day of the month of the capture event
4. month - the month of the capture event
5. year - the year of the capture event
6. location - a rough location of the capture event within Fort Hunter Liggett, USA
7. method - type of trap used to capture the animal
8. type - whether the individual pig was a new capture (initial) or a recapture
9. habitat - a rough habitat category of the trapping location
10. drugged - whether the individual pig was chemically immobilized
11. drug_administered - what drug was used if chemical immobilization occurred
12. reversal_administered - what reversal drug was used
13. collar_id - the serial number for the GPS collar deployed on the pig
14. Collar Frequency - the VHF frequency of the deployed collar
15. ear_tag_L_color - the color of the ear tag affixed to the pigs left ear
16. ear_tag_L_number - the number of the ear tag affixed to the pigs left ear
17. ear_tag_R_color - the color of the ear tag affixed to the pigs right ear
18. ear_tag_R_number - the number of the ear tag affixed to the pigs right ear
19. vhf_eartag_ear - whether or not the pig received a VHF eartag
20. vhf_eartag_freq - the VHF frequency of the deployed eartag
21. age_class - the age class of the individual pig
22. nursing/pregnant - whether female pigs were nursing or pregnant
23. current/past - whether nursing/pregnancy was current or previous
24. est_weight - the estimated weight of the pig
25. hair_taken - whether a hair sample was collected
26. blood - whether a blood sample was collected



##### PigEartagLocations3.csv

This csv file includes monitoring data for all tagged and collared pigs used in the survival analysis. Each row represents a detection of a pig and provides information on the pigs location and status. There are 7 columns:

1. Date - the date the pig was detected
2. Pig_ID - the ID of the individual pig
3. Training_Area - the training area (location within Fort Hunter Liggett, CA) that the pig was located 
4. Signal_type - whether the detected collar/tag signal was in normal operation or mortality mode
5. Status - whether the pig was alive (0) or dead (1)
6. Cause - If a mortality occurred, the determined cause of mortality
7. Notes - Any relevant notes related to the detection. Examples include notes on mortalities or comments on tag malfunctions (e.g., tracked a mortality signal to a live pig)

##### Cougar GPS collar data

We are unable to share these raw data because of their sensitive nature.

##### cougar_cluster_investigations

We are unable to share these raw data because of their sensitive nature.

### CODE FILES

1. Estimate pig survival rates
2. Estimate cougar home range areas
3. Estimate RSF of cougar kills for cougar predation risk
4. Estimate SSF of pig GPS collars to understand if they avoid areas of high cougar risk.

##### 01_Pig survival_to share.R 

This file uses pig capture history data to estimate survival rates using Kaplan-Meier esitmation method.
Then uses cox proportional hazard models to estimate if survival varies by age class.
Then creates plots of survival curves, estimated rates +/- 95% CIs
Then runs chi squared test to see if causes of mortality vary by age class.

##### 02_HomeRanges_NewCollars.R

This file uses cougar GPS collar data and creates home ranges.
After inital processing and data cleaning, uses AMT to estiamte 95% autocorrelated kernel density estimates.
And 100% minimum convex polygons for three collared cougars: SF9, SF1, and SM6.
Currently only writes shapefiles for MCPs since that is what is used in cougar RSF.
Note: We cannot provide the raw data to run this code, but provide the code as an illustration of what we did.

##### 03_RSF_multiple-animals.R

This file estimates an resource selection function of cougar kill sites.
It then uses this RSF to create a layers of predicted probability of cougar kill occurrence.
Uses cluster investigation data, cougar MCPs, and GIS data.
Outputs include the parameter estimates for kill site RSF and the cougar risk layers.
Note: We cannot provide the raw data to run this code, but provide the code as an illustration of what we did.
We provide the produced cougar kill occurrene layer which is used in the pig step selection function.

##### 04_SSF_PigCollars.R

This file estimates integrated step selection functions for wild pigs.
It first processes and cleans pig GPS collar data and adds sunrise/sunset times to identify day and night locations.
Then reads in GIS data, including the cougar risk layer.
Finally it estimates an integrated step selection function for pigs both with and without cougar risk included.
The final section compares models with and without cougar risk to understand if it improves model performance.
There is also code for creating plots of pig 50% core area kernel density estimate utilization distributions plotted over cougar risk.