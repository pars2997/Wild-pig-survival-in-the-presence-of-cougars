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
dat <- dat[dat$Date < mdy("02/01/2024"),]

# 
# months <- 1:12
# years <- min(year(dat$Date)):max(year(dat$Date))
# 
# pigcount <- c()
# 
# for(i in 1:length(years)){
#   for(j in 1:length(months)){
#     mindate <- mdy(paste0(months[j],"/01/",years[i]))
#     maxdate <- mdy(paste0(months[j],"/28/",years[i]))
#     temp <- survdata2 %>% 
#       filter(capturedate < mindate & lastdate > maxdate)
#     x <- length(unique(temp$pig_id))
#     pigcount <- c(pigcount,x)
#   }
# }
# 
# pigcount.df <- data.frame(Month = rep(months, 4),
#                           Year = rep(years, each = 12), 
#                           Count = pigcount)
# pigcount.df <- pigcount.df[5:35,]
# 
# pigcount.summer <- pigcount.df %>% 
#   filter(Month %in% c(5,6,7))
# pigcount.winter <- pigcount.df %>% 
#   filter(Month %in% c(11,12,1))
# t.test(pigcount.summer$Count,pigcount.winter$Count)
# 
# table(month(survdata2$lastdate))

dat$Pig_ID[dat$Pig_ID == "MJ23"] <- "MS09"

# this file is actual capture events of pigs, with dates and all covariates
capturedat <- read.csv("../RawData/pig_captures_MJ23_update.csv")
# filter to just those pigs that received VHF or GPS and could be monitored
capturedat <- capturedat[capturedat$pig_id %in% dat$Pig_ID,]
capturedat <- capturedat %>% 
  mutate(location = case_when(
    location %in% c("peanut","stonyres") ~ "core",
    location %in% c("riverrd","TA22") ~ "edge"
  ))

# select only needed columns
capturedat <- capturedat %>% 
  dplyr::select(pig_id, date, location, age_class,est_weight,collar_id)

# remove the couple recaptures because only need them once
capturedat <- capturedat[-which(duplicated(capturedat$pig_id)),]
# rename columns
colnames(capturedat) <- c("pig_id","capturedate", "location","age","weight","collar_id")


#### Data formating for survival analysis ####

# Create a vector of the last known alive dates for each indiviudal pig
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
survdata <- survdata[,c(1:9,11,12,13)]
survdata$time <- as.numeric(survdata$time)
IDs <- unique(survdata$pig_id)


# Create a new object (survdata2) that will have the capture date, last known
# alive date, mortality detection date, and all covariates for each pig
survdata2 <- survdata[2,]

for(i in 2:length(IDs)){
  temp <- survdata[survdata$pig_id == IDs[i],]
  temp2 <- temp[temp$Date == max(temp$Date),]
  survdata2 <- rbind(survdata2,temp2)
  
}

# survdata2 <- survdata2[!duplicated(survdata2$pig_id),]


# Add column for what type of tag each pig had (GPS collar vs VHF eartag)
survdata2 <- survdata2 %>% 
  mutate(tag_type = case_when(
    is.na(collar_id) ~ "eartag",
    !is.na(collar_id) ~ "collar"
  ))

head(survdata2)
#summarize tag type data
table(survdata2$tag_type)

# how many animals were successfully followed for a full year
sum(survdata2$tag_type == "collar" & (survdata2$time > 365 | survdata2$Status == 1))
sum(survdata2$tag_type == "eartag" & (survdata2$time > 365 | survdata2$Status == 1))
sum(survdata2$tag_type == "eartag" & (survdata2$time > 180 | survdata2$Status == 1))
sum(survdata2$tag_type == "collar" & (survdata2$time > 180 | survdata2$Status == 1))

sum(survdata2$age == "A" & (survdata2$time > 365 | survdata2$Status == 1))
sum(survdata2$age == "S" & (survdata2$time > 365 | survdata2$Status == 1))
sum(survdata2$age == "J" & (survdata2$time > 365 | survdata2$Status == 1))
sum(survdata2$age == "A" & (survdata2$time > 180 | survdata2$Status == 1))
sum(survdata2$age == "S" & (survdata2$time > 180 | survdata2$Status == 1))
sum(survdata2$age == "J" & (survdata2$time > 180 | survdata2$Status == 1))

tapply(survdata2$time,survdata2$age,summary)

mean(survdata2$time[survdata2$age == "A" & survdata2$Status == 1])
mean(survdata2$time[survdata2$age == "S" & survdata2$Status == 1])
mean(survdata2$time[survdata2$age == "J" & survdata2$Status == 1])
sd(survdata2$time[survdata2$age == "A" & survdata2$Status == 1])/length(survdata2$time[survdata2$age == "A" & survdata2$Status == 1])
sd(survdata2$time[survdata2$age == "S" & survdata2$Status == 1])/length(survdata2$time[survdata2$age == "S" & survdata2$Status == 1])
sd(survdata2$time[survdata2$age == "J" & survdata2$Status == 1])/length(survdata2$time[survdata2$age == "J" & survdata2$Status == 1])

# survdata2 <- survdata2 %>%
#   mutate(age = case_when(
#     weight < 60 ~ "J",
#     weight >= 60 & weight <= 90 ~ "S",
#     weight > 90 ~ "A"
#   ))

# survdata2 <- survdata2 %>%
#   mutate(ageclass2 = case_when(
#     age == "J" ~ "J",
#     TRUE ~ "A"
#   ))
# 
# 
# 
survdata2 <- survdata2 %>%
  mutate(weightclass = case_when(
    weight <= 33 ~ "S",
    weight > 33 & weight <= 88 ~ "M",
    weight > 88 ~ "L"
  ))


survdata2$weight.kg <- round(survdata2$weight/2.2,2)
table(survdata2$weight.kg)

tapply(survdata2$weight.kg,survdata2$age,mean)

se <- function(x){
  se = sd(x)/sqrt(length(x))
  return(se)
}
tapply(survdata2$weight.kg,survdata2$age,se)

a <- survdata2$weight[survdata2$age == "J" & survdata2$Status == 1]
b <- survdata2$weight[survdata2$age == "J" & survdata2$Status == 0]
t.test(a,b)

# # Add in one pig that was captured but not monitored 
# # and then captured again later, so know survived
# capturedat2 <- read.csv("../RawData/pig_captures_MJ23_update.csv")
# capturedat2 <- capturedat2 %>% 
#   mutate(location = case_when(
#     location %in% c("peanut","stonyres") ~ "core",
#     location %in% c("riverrd","TA22") ~ "edge"
#   ))
# MJ23data <- data.frame(pig_id = "MJ23",
#                        capturedate = mdy(capturedat2$date[capturedat2$pig_id == "MJ23"]),
#                        location = capturedat2$location[capturedat2$pig_id == "MJ23"],
#                        age = capturedat2$age_class[capturedat2$pig_id == "MJ23"],
#                        weight = capturedat2$est_weight[capturedat2$pig_id == "MJ23"],
#                        collar_id = capturedat2$collar_id[capturedat2$pig_id == "MJ23"],
#                        lastdate = mdy(capturedat2$date[capturedat2$pig_id == "MS09"]),
#                        time = NA,
#                        Date = mdy(capturedat2$date[capturedat2$pig_id == "MS09"]),
#                        Signal_type = "Normal",
#                        Status = 0,
#                        Cause = "Collar removed",
#                        ageclass2 = NA,
#                        tag_type = "eartag")
# MJ23data <- MJ23data %>% 
#   mutate(age = "J",
#          time = as.numeric(lastdate - capturedate),
#          ageclass2 = case_when(
#            weight <= 40 ~ "S",
#            weight > 40 & weight <= 80 ~ "M",
#            weight > 80 ~ "L"
#          ))



#### Survival analysis ####

# Create survival object
pigsurv <- Surv(survdata2$time,survdata2$Status)

# Run kaplan-meier survival model with age as a covariate
mod1 <- survfit(pigsurv ~ age, data = survdata2)
summary(mod1)
# look at estimates for annual and 6-month survival rates for each age class
summary(mod1, times = 365, extend = T)
summary(mod1, times = 180, extend = T)



# Run cox proportional hazrd model to determine if survival rate varies by age class
cph1 <- coxph(pigsurv ~ age, data = survdata2)
summary(cph1)
cph2 <- coxph(pigsurv ~ weightclass, data = survdata2)
summary(cph2)

AICcmodavg::AICc(cph1)
AICcmodavg::AICc(cph2)

# Suggests juveniles have a marginally lower survival rate than adults
# but no difference between adults and subadults

# order ages in age order instead of alphabetical order
# note that this makes cph coefficients wonky and I don't know why
# so do this, then rerun kaplan-meier

survdata2$age <- ordered(survdata2$age, levels = c("A","S","J"))
mod1 <- survfit(pigsurv ~ age, data = survdata2)


# Extract means and confidence intervals for plotting
# do this for both annual and 6 month rates
means <- summary(mod1, times = 180,extend = T)[["surv"]]
lcl <- summary(mod1, times = 180,extend = T)[["lower"]]
ucl <- summary(mod1, times = 180,extend = T)[["upper"]]

#combine estimates into dataframe
cbind(means,lcl,ucl)


setEPS()
postscript("../Figures/Parsons-et-al_Figure-1.eps")
#plot estimates +/- 95% CIs for each age calss
# par(mfrow = c(1,2))
par(mar = c(5,5,2,2))
plot(means ~ c(1,2,3), xlab = "Age class", ylab = "Six-month survival",
     cex.axis = 1.15, cex.lab = 1.4, pch = 19, cex = 2,xlim = c(0.5,3.5),
     ylim = c(0,1),xaxt = "n",bty = "l")
arrows(x0 = c(1,2,3),y0 = lcl, y1 = ucl, angle = 90, code = 3,length = 0.18)
axis(side = 1, at = c(1,2,3), labels = c("Adult","Subadult", "Juvenile"))
# text(x = 0.5, y = 1, label = "A", cex = 1.5)

dev.off()

# repeat for annual survival
means <- summary(mod1, times = 365,extend = T)[["surv"]]
lcl <- summary(mod1, times = 365,extend = T)[["lower"]]
ucl <- summary(mod1, times = 365,extend = T)[["upper"]]


#combine estimates into dataframe
cbind(means,lcl,ucl)
par(mar = c(5,4,2,2))
#plot estimates +/- 95% CIs for each age calss
plot(means ~ c(1,2,3), xlab = "Age class", ylab = "Annual survival",
     cex.axis = 1.15, cex.lab = 1.4, pch = 19, cex = 2,xlim = c(0.5,3.5),
     ylim = c(0,1),xaxt = "n",bty = "l")
arrows(x0 = c(1,2,3),y0 = lcl, y1 = ucl, angle = 90, code = 3,length = 0.18)
axis(side = 1, at = c(1,2,3), labels = c("Adult","Subadult", "Juvenile"))
text(x = 0.5, y = 1, label = "B", cex = 1.5)

# Make plot of survival curves
# Save as object so you can save to png file for publication
survcurv <- ggsurvplot(mod1, data = survdata2,
                       palette = "grey",
                       legend = "right",
                       legend.title = "Age class",
                       legend.labs = c("Adult", "Subadult", "Juvenile"),
                       xlim = c(0,200),
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

setEPS()
postscript("../Figures/Parsons-et-al_Figure-2.eps")
survcurv
dev.off()


#### Post hoc survival analysis ####

# Run kaplan-meier survival model with size as a covariate
mod2 <- survfit(pigsurv ~ weightclass, data = survdata2)
summary(mod2)

# Run cox proportional hazrd model to determine if survival rate varies by age class
cph2 <- coxph(pigsurv ~ weightclass, data = survdata2)
summary(cph2)

AICcmodavg::AICc(cph1)
AICcmodavg::AICc(cph2)

# order ages in age order instead of alphabetical order
# note that this makes cph coefficients wonky and I don't know why
# so do this, then rerun kaplan-meier

# survdata2$ageclass2 <- ordered(survdata2$ageclass2, levels = c("L","M","S"))
# mod2 <- survfit(pigsurv ~ ageclass2, data = survdata2)


# Extract means and confidence intervals for plotting
# do this for both annual and 6 month rates
means <- summary(mod2, times = 180,extend = T)[["surv"]]
lcl <- summary(mod2, times = 180,extend = T)[["lower"]]
ucl <- summary(mod2, times = 180,extend = T)[["upper"]]

#combine estimates into dataframe
cbind(means,lcl,ucl)

#plot estimates +/- 95% CIs for each age calss
par(mfrow = c(1,2))
par(mar = c(5,5,2,2))
plot(means ~ c(1,2,3), xlab = "Age class", ylab = "Six-month survival",
     cex.axis = 1.15, cex.lab = 1.4, pch = 19, cex = 2,xlim = c(0.5,3.5),
     ylim = c(0,1),xaxt = "n",bty = "l")
arrows(x0 = c(1,2,3),y0 = lcl, y1 = ucl, angle = 90, code = 3,length = 0.18)
axis(side = 1, at = c(1,2,3), labels = c("Large \n(>40 kg)", "Subadult \n(20-40 kg)", "Juvenile \n(<20kg)"))
text(x = 0.5, y = 1, label = "A", cex = 1.5)
# repeat for annual survival
means <- summary(mod1, times = 365,extend = T)[["surv"]]
lcl <- summary(mod1, times = 365,extend = T)[["lower"]]
ucl <- summary(mod1, times = 365,extend = T)[["upper"]]

#combine estimates into dataframe
cbind(means,lcl,ucl)
par(mar = c(5,4,2,2))
#plot estimates +/- 95% CIs for each age calss
plot(means ~ c(1,2,3), xlab = "Age class", ylab = "Annual survival",
     cex.axis = 1.15, cex.lab = 1.4, pch = 19, cex = 2,xlim = c(0.5,3.5),
     ylim = c(0,1),xaxt = "n",bty = "l")
arrows(x0 = c(1,2,3),y0 = lcl, y1 = ucl, angle = 90, code = 3,length = 0.18)
axis(side = 1, at = c(1,2,3), labels = c("Large \n(>40 kg)", "Subadult \n(20-40 kg)", "Juvenile \n(<20kg)"))
text(x = 0.5, y = 1, label = "B", cex = 1.5)

# Make plot of survival curves
# Save as object so you can save to png file for publication
survcurv <- ggsurvplot(mod2, data = survdata2,
                       palette = "grey",
                       legend = "right",
                       legend.title = "Size",
                       legend.labs = c("Large (>40 kg)", "Subadult (20-40 kg)", "Juvenile (<20kg)"),
                       xlim = c(0,200),
                       xlab = "Days",
                       font.x = 16,
                       font.y = 16,
                       font.legend = 12,
                       font.tickslab = 14,
                       lwd = 1.5)


#### Cause specific mortality investigation ####
# create subsetted dataframe only of mortalities 
temp <- survdata2[survdata2$Cause %in% c("Hunter","Lion","Predator",
                                         "Roadkill","Unknown mortality"),]
# code the unknown predator as a lion kill
temp$Cause[temp$Cause == "Predator"] <- "Unknown mortality"
# create an age x cause contingency table
mort.tab <- table(temp$Cause,temp$age)
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
