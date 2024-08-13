#####################################################
######## Pig survival modeling with Kaplan-Meier ####
############## Mitchell Parsons #####################
################ February 22, 2024 ##################

#### Load packages ####
library(lubridate)
library(tidyverse)
library(survival)
library(ggplot2)
library(survminer)
library(lmtest)

#### Read in data ####
# this file is pig capture histories with date of detection and status
dat <- read.csv("../RawData/pigeartaglocations3.csv")
dat$Date <- mdy(dat$Date)
head(dat)

# Remove location for recapture of unmonitored pig
dat <- dat[-which(dat$Pig_ID == "MJ23" & dat$Date == mdy("06/01/2023")),]

# Rename one pig to reflect class it was monitored as
dat$Pig_ID[dat$Pig_ID == "MJ23"] <- "MS09"

# this file is actual capture events of pigs, with dates and all covariates
capturedat <- read.csv("../RawData/pig_captures_MJ23_update.csv")
# filter to just those pigs that received VHF or GPS and could be monitored
capturedat <- capturedat[capturedat$pig_id %in% dat$Pig_ID,]

# select only needed columns
capturedat <- capturedat %>% 
  dplyr::select(pig_id, date, age_class,est_weight,collar_id)

# remove the couple recaptures because only need them once
capturedat <- capturedat[-which(duplicated(capturedat$pig_id)),]
# rename columns
colnames(capturedat) <- c("pig_id","capturedate","age","weight","collar_id")


#### Data formating for survival analysis ####

# Create a vector of the last known alive dates for each individual pig
lastdates <- data.frame(cbind(names(tapply(dat[dat$Signal_type == "Normal",]$Date,dat[dat$Signal_type == "Normal",]$Pig_ID,max)),
                 as.numeric(tapply(dat[dat$Signal_type == "Normal",]$Date,dat[dat$Signal_type == "Normal",]$Pig_ID,max))))

# Name columns and correct date format
colnames(lastdates) <- c("pig_id","lastdate")
lastdates$lastdate <- as.Date(as.numeric(lastdates$lastdate),
                              origin = "1970-01-01")

# Combine capture data (i.e. covariates) with last dates
survdata <- capturedat %>% 
  left_join(lastdates, by = "pig_id") %>% 
  mutate(lastdate = ymd(lastdate))

#format dates and rename columns
survdata <- survdata %>% 
  mutate(capturedate = as.Date(as.character(survdata$capturedate),format = "%m/%d/%Y")) %>% 
  mutate(time = difftime(lastdate,capturedate,units = "days"))

# Add in all detections for creating survival histories
colnames(dat) <- c("Date","pig_id","Training_Area","Signal_type","Status","Cause","Notes")
survdata <- survdata %>% 
  left_join(dat)

# Filter to needed columns and create a vector of unique pig IDs for summarizing
survdata <- survdata[,c(1:8,10,11,12)]
survdata$time <- as.numeric(survdata$time)
IDs <- unique(survdata$pig_id)

# Add pig weight estimates from Mayer 2021.
# Body Mass Variation in an Introduced Wild Pig Population
# With Changing Ancestry
pigweights <- data.frame(Ageclass = c("Neonate","Piglet","Juvenile","Yearling","Subadult","Adult"),
                         Age.Months = c(1,5,11,16,27,36),
                         Female.Mass = c(0.754,6.5,26.4,46,62.7,80),
                         Male.Mass = c(0.818,6.7,27.6,49.1,66.9,91.1))

# Calculate average mass and age in days
pigweights <- pigweights %>% 
  mutate(Avg.Mass = (Female.Mass + Male.Mass)/2,
         Age.Days = Age.Months*30)

# Calculate the growth rate during each interval
growth <- c()
for(i in 1:nrow(pigweights)){
  growth[i] <- (pigweights$Avg.Mass[i+1] - pigweights$Avg.Mass[i])/(pigweights$Age.Days[i+1] - pigweights$Age.Days[i])
}
pigweights$growth.rate <- growth

# Plot mass by age and fit linear and quadratic regressions to the data
plot(pigweights$Avg.Mass ~ pigweights$Age.Days, pch = 16)
linmod <- lm(pigweights$Avg.Mass ~ pigweights$Age.Days)
summary(linmod)
abline(linmod)

sqmod <- lm(pigweights$Avg.Mass ~ pigweights$Age.Days + I(pigweights$Age.Days^2))
summary(sqmod)
x <- 0:1500
y <- sqmod$coefficients[1] + sqmod$coefficients[2]*x + sqmod$coefficients[3] * x^2
lines(x,y)

# Check to see which model performs best by AIC
AIC(linmod)
AIC(sqmod)

# Conduct likelihood ratio test of the two models
linLL <- logLik(linmod)
sqLL <- logLik(sqmod)
teststat <- -2 * (as.numeric(linLL) - as.numeric(sqLL))
p.val <- pchisq(teststat,df = 1, lower.tail = F)
p.val
lrtest(linmod,sqmod)

# linear model performs as well as quadratic
# predicts growth of 0.08 kg/day

# Extract coefficient from linear model to use as estimated growth rate
# multiple by 2.2 to convert kg to lbs
growth.rate.lbs <- as.numeric(linmod$coefficients[2]*2.2)

# Add estimated weights to each observation of a pig
# Estimated weight is the capture weight plus the growth rate times the days since capture
# Only do this for juveniles and subadults since adults are less predictable
survdata <- survdata %>% 
  group_by(pig_id) %>% 
  mutate(days_since_capture = as.numeric(Date - capturedate)) %>% 
  ungroup() %>% 
  mutate(est.weight = case_when(
    age %in% c("S","J") ~ weight + days_since_capture*growth.rate.lbs,
    age %in% c("A") ~ weight
  ))

# Create updated data for each pig where each row is an detection
# Detections include the survival period (i.e., previous to current detection)
# And estimated weight at detection
updated <- list()
for(i in 1:length(IDs)){
  starttime <- c()
  starttime[1] <- 0
  temp <- survdata[survdata$pig_id == IDs[i],]
  temp <- temp[order(temp$Date),]
  for(j in 2:nrow(temp)){
    starttime[j] <- temp$days_since_capture[j-1]
  }
  temp$startime <- starttime
  updated[[i]] <- temp
}

# Bind rows into new dataframe
survdata4 <- updated[[1]]
for(i in 2:length(updated)){
  survdata4 <- rbind(survdata4, updated[[i]])
}

# add in end time
survdata4$endtime <- survdata4$days_since_capture

# Add column for age class that updates with changing estimated weight of pigs
survdata4 <- survdata4 %>% 
  mutate(ageclass.updating =case_when(
    est.weight < 60 ~ "J",
    est.weight >= 60 & est.weight < 100 ~ "S",
    est.weight >= 100 ~ "A"))

# Add column for what type of tag each pig had (GPS collar vs VHF eartag)
survdata4 <- survdata4 %>%
  mutate(tag_type = case_when(
    is.na(collar_id) ~ "eartag",
    !is.na(collar_id) ~ "collar"
  ))

head(survdata4)
#summarize tag type data
table(survdata4$tag_type[!duplicated(survdata4$pig_id)])

# number of collars monitored for >6 months
survdata4 %>% 
  group_by(pig_id) %>% 
  summarise(tag_type = first(tag_type),
            time = last(time),
            Status = last(Status)) %>% 
  filter(tag_type == "collar" & (time > 180 | Status == 1)) %>% 
  nrow

# number of eartags monitored for >6 months
survdata4 %>% 
  group_by(pig_id) %>% 
  summarise(tag_type = first(tag_type),
            time = last(time),
            Status = last(Status)) %>% 
  filter(tag_type == "eartag" & (time > 180 | Status == 1)) %>% 
  nrow()

# number of eartags monitored for >6 months
survdata4 %>% 
  group_by(pig_id) %>% 
  summarise(age = first(age),
            time = last(time),
            Status = last(Status)) %>% 
  filter(age == "A" & (time > 180 | Status == 1)) %>% 
  nrow()

# number of eartags monitored for >6 months
survdata4 %>% 
  group_by(pig_id) %>% 
  summarise(age = first(age),
            time = last(time),
            Status = last(Status)) %>% 
  filter(age == "S" & (time > 180 | Status == 1)) %>% 
  nrow()

# number of eartags monitored for >6 months
survdata4 %>% 
  group_by(pig_id) %>% 
  summarise(age = first(age),
            time = last(time),
            Status = last(Status)) %>% 
  filter(age == "J" & (time > 180 | Status == 1)) %>% 
  nrow()


#### Survival analysis ####
# Convert weights from lbs to kg
survdata4$weight <- survdata4$weight/2.2
survdata4$est.weight <- survdata4$est.weight/2.2

survdata4$loc.month <- month(survdata4$Date)
survdata4 <- survdata4 %>% 
  mutate(season = case_when(
    loc.month %in% c(2,3,4) ~ "spring",
    loc.month %in% c(5,6,7) ~ "summer",
    loc.month %in% c(8,9,10) ~ "fall",
    loc.month %in% c(11,12,1) ~ "winter"
  ))

# Create survival object
pigsurv <- Surv(survdata4$startime,survdata4$endtime,survdata4$Status)

mod1 <- survfit(pigsurv ~ ageclass.updating, data = survdata4, id = pig_id)
summary(mod1)

mod2 <- survfit(pigsurv ~ est.weight, data = survdata4, id = pig_id)
summary(mod2)

# look at estimates for annual and 6-month survival rates for each age class
summary(mod1, times = 180, extend = T)


# Run cox proportional hazrd model to determine if survival rate varies by age class
cph1 <- coxph(pigsurv ~ ageclass.updating, data = survdata4, id = pig_id)
summary(cph1)

cph2 <- coxph(pigsurv ~ est.weight + season, data = survdata4, id = pig_id)
summary(cph2)
# Suggests juveniles have a marginally lower survival rate than adults
# but no difference between adults and subadults

# order ages in age order instead of alphabetical order
# note that this makes cph coefficients wonky and I don't know why
# so do this, then rerun kaplan-meier
survdata4$ageclass.updating <- ordered(survdata4$ageclass.updating, levels = c("A","S","J"))
mod1 <- survfit(pigsurv ~ ageclass.updating, data = survdata4, id = pig_id)

# Extract means and confidence intervals for plotting
# do this for both annual and 6 month rates
means <- summary(mod1, times = 180,extend = T)[["surv"]]
lcl <- summary(mod1, times = 180,extend = T)[["lower"]]
ucl <- summary(mod1, times = 180,extend = T)[["upper"]]

#combine estimates into dataframe
cbind(means,lcl,ucl)

#plot estimates +/- 95% CIs for each age class
# par(mfrow = c(1,2))
par(mar = c(5,5,2,2))
plot(means ~ c(1,2,3), xlab = "Age class", ylab = "Six-month survival",
     cex.axis = 1.15, cex.lab = 1.4, pch = 19, cex = 2,xlim = c(0.5,3.5),
     ylim = c(0,1),xaxt = "n",bty = "l")
arrows(x0 = c(1,2,3),y0 = lcl, y1 = ucl, angle = 90, code = 3,length = 0.18)
axis(side = 1, at = c(1,2,3), labels = c("Adult", "Subadult", "Juvenile"))
text(x = 0.5, y = 1, label = "A", cex = 1.5)
# # repeat for annual survival
# means <- summary(mod1, times = 365,extend = T)[["surv"]]
# lcl <- summary(mod1, times = 365,extend = T)[["lower"]]
# ucl <- summary(mod1, times = 365,extend = T)[["upper"]]
# 
# #combine estimates into dataframe
# cbind(means,lcl,ucl)
# par(mar = c(5,4,2,2))
# #plot estimates +/- 95% CIs for each age calss
# plot(means ~ c(1,2,3), xlab = "Age class", ylab = "Annual survival",
#      cex.axis = 1.15, cex.lab = 1.4, pch = 19, cex = 2,xlim = c(0.5,3.5),
#      ylim = c(0,1),xaxt = "n",bty = "l")
# arrows(x0 = c(1,2,3),y0 = lcl, y1 = ucl, angle = 90, code = 3,length = 0.18)
# axis(side = 1, at = c(1,2,3), labels = c("Adult", "Subadult", "Juvenile"))
# text(x = 0.5, y = 1, label = "B", cex = 1.5)

# Make plot of survival curves
# Save as object so you can save to png file for publication
survcurv <- ggsurvplot(mod1, data = survdata4,
                       palette = "grey",
                       legend = "right",
                       legend.title = "Age class",
                       legend.labs = c("Adult", "Subadult", "Juvenile"),
                       xlim = c(0,400),
                       xlab = "Days",
                       font.x = 16,
                       font.y = 16,
                       font.legend = 12,
                       font.tickslab = 14,
                       lwd = 1.5)

png(filename = "../Figures/pigsurvivialcurve.png", width = 5, height = 4,
    units = "in",res = 300)
survcurv
dev.off()

#### Cause specific mortality investigation ####
# create subsetted dataframe only of mortalities 
temp <- survdata4[survdata4$Cause %in% c("Hunter","Lion","Predator",
                                         "Roadkill","Unknown mortality"),]
# code the unknown predator as a lion kill
temp$Cause[temp$Cause == "Predator"] <- "Lion"
# create an age x cause contingency table
mort.tab <- table(temp$Cause,temp$ageclass.updating)
# conduct chisq.test to see if causes of mortality vary by age class
chisq.test(mort.tab)

temp
temp2 <- temp %>% 
  group_by(age,Cause) %>% 
  summarise(TotalMorts = n()) %>% 
  group_by(age) %>% 
  mutate(Total.age = sum(TotalMorts),
         prop.cause = TotalMorts/Total.age) %>% 
  mutate(age2 = case_when(
    age %in% c("A","S") ~ "Adult",
    age %in% c("J") ~ "Juvenile"
  )) %>% 
  mutate(Cause = case_when(
    Cause == "Hunter" ~ "Hunter",
    Cause == "Lion" ~ "Cougar",
    Cause == "Roadkill" ~ "Roadkill",
    Cause == "Unknown mortality" ~ "Unknown"
  ))

ggplot(temp2, aes(x = Cause, y = TotalMorts, fill = age2))+
  geom_bar(stat = "identity")+
  facet_wrap(~ age2)+
  scale_fill_manual(name = "Age class", values = c("burlywood","darkolivegreen4"))+
  theme_classic()+
  xlab("Cause of Mortality")+
  ylab("Frequency")+
  theme(axis.title = element_text(size = 20))+
  theme(axis.text = element_text(size = 16))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 18))+
  theme(strip.background = element_blank())+
  theme(strip.text = element_blank())
