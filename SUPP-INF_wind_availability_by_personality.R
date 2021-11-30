
### AIM: EXamine whether wind conditions vary for birds of different personality types

library(lme4); library(MuMIn); library(ggplot2); library(visreg); library(plyr); library(circular);
library(dplyr)


# PREAMBLE ---------------------------------------------------

dat <- data.table::fread("Data_inputs/WAAL_GPS_2010-2021_personality_wind.csv")
dat$DateTime <- fasttime::fastPOSIXct(dat$DateTime)

## Transform ring sex and year to factor
dat[,c(1:3)] <- lapply(dat[,c(1:3)], as.character) 
dat[,c(1:3)] <- lapply(dat[,c(1:3)], as.factor) 


# MODEL WIND SPEED AVAILABILITY -------------------------------------------------

### Look at whether wind availability is predicted by personality and sex

# Remove extremes of personality and rescale
dat.wind <- subset(dat, mean_BLUP_logit < 2 & mean_BLUP_logit > -2)
dat.wind$persScaled <- scales::rescale(dat.wind$mean_BLUP_logit, to = c(-2,2))

## Build global model: interested in interaction of sex and personality (as we known sex has an effect)
wind.lm <- lmerTest::lmer(WindSp ~ Sex*persScaled + Year + (1|Ring/Trip_bout), data = dat.wind, REML = FALSE)
options(na.action = "na.fail") 

## Compare all possible models 
m_set <- dredge(wind.lm)
m_set

#Model selection table 
#(Int) men_BLU_lgt Sex Yer men_BLU_lgt:Sex df   logLik    AICc  delta weight
#16 8.443     0.39980   +   +               + 17 -1284627 2569288   0.00      1
#7  8.277               +   +                 15 -1284653 2569336  48.07      0
#8  8.261    -0.06452   +   +                 16 -1284652 2569336  48.87      0
#12 8.268     0.42800   +                   +  7 -1284667 2569348  60.24      0
#6  8.523    -0.11960       +                 15 -1284678 2569386  98.28      0
#5  8.561                   +                 14 -1284680 2569388 100.11      0
#3  8.310               +                      5 -1284692 2569393 105.57      0
#4  8.315    -0.04718   +                      6 -1284691 2569395 106.96      0
#2  8.756    -0.09559                          5 -1284712 2569434 146.78      0
#1  8.757                                      4 -1284713 2569435 147.13      0


# Best model includes personality, sex, and their interaction

plot(wind.lm)
plot(residuals(wind.lm, type = "pearson") ~ dat$mean_BLUP_logit)
plot(residuals(wind.lm, type = "pearson") ~ dat$Sex)
plot(residuals(wind.lm, type = "pearson") ~ dat$Year)

summary(wind.lm)
emmeans::emmeans(wind.lm, pairwise~Sex)


# MODEL WIND DIRECTION AVAILABILITY ---------------------------------------------------

### Merge in wind data
wind.df <- read.csv("Data_inputs/WAAL_GPS_wind-directions.csv")
wind.df$Ring <- as.factor(as.character((wind.df$Ring)))
wind.df$DateTime <- fasttime::fastPOSIXct(wind.df$DateTime)

# Merge with main dataset
dat <- merge(dat, wind.df, all.x = TRUE, by = c("Ring", "DateTime"))

### Get median wind direction per trip
dat$WindDir2 <- circular(dat$WindDir,
                         type = "angles",
                         units = "degrees", 
                         zero = 0, 
                         rotation = "counter")

dat <- subset(dat, !is.na(WindDir2))
trips.sum <- ddply(dat,.(Sex, Trip_bout, Ring, mean_BLUP_logit), summarize, MedDir = median.circular(WindDir2))

hist(as.numeric(trips.sum$MedDir))
range(as.numeric(trips.sum$MedDir))


## Calculate mean wind direction for all individuals - mean and confidence intervals
mle.vonmises(trips.sum$MedDir)
# mu: 100.3  ( 2.54 )

# kappa: 1.699  ( 0.1008 )


## Bootstrap CIs
mle.vonmises.bootstrap.ci(trips.sum$MedDir, alpha = 0.05)
#Mean Direction:            Low = 95.65   High = 104.72 
#Concentration Parameter:   Low = 1.48   High = 1.98 

# function to convert to where wind is coming FROM on 0-360 scale
dir_conv <- function(x) { if (x > 0) { dir <- 360-(180-x)} else { dir <- 180+x }
  return(dir)
}
dir_conv(100.3) # 280.3
dir_conv(104.72) # 284.72
dir_conv(95.65) # 275.65


### Circular regression for personality

# Circular outcome in radians
rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}

trips.sum$radCirc <- deg2rad(trips.sum$MedDir)

# Bird ID needs to be a numeric index
trips.sum$Ring <- as.numeric(factor(trips.sum$Ring))
trips.mod <- trips.sum[,c(3,6,4,1)]

## Build circular regression model: wind direction ~ personality 

library(bpnreg); library(Rcpp)

system.time(fit.wind <- bpnme(pred.I = radCirc ~ mean_BLUP_logit + Sex + (1|Ring), data = trips.mod,
                              its = 1, burn = 100, n.lag = 3, seed = 101))


print(fit.wind)

#Linear Coefficients 

#Component I: 
#                         mean       mode        sd      LB HPD     UB HPD
#(Intercept)        -0.2608904 -0.2441189 0.1030941 -0.4614688 -0.06009470
#mean_BLUP_logit -0.1528597 -0.1565361 0.1244041 -0.4000871  0.08060275
#SexM               -0.2824122 -0.2756244 0.1491180 -0.5642317  0.01570623

# HPD interval = smallest interval in which 95% posterior distribution located
## HPD intervals for personality include 0, so on average no effect of personality for 
## experienced wind conditions


## Model comparison

fit.wind.null <- bpnme(pred.I = radCirc ~ Sex + (1|birdInd), data = trips.mod,
                       its = 1, burn = 100, n.lag = 5, seed = 181)

## Important coefficients:
# bc = slope of circular regression line at inflection point
# AS = average slope over all x
# SAM = slope at mean of x

#                                mean       mode           sd      LB HPD     UB HPD
#mean_BLUP_logit bc    22.63972436 -0.4531313 2126.0193203 -23.3375034 23.6839226
#mean_BLUP_logit AS    -0.70572370  0.2936290   88.7078513 -10.7990076  9.3712756
#mean_BLUP_logit SAM  -10.77132878  0.3685708  862.4169920 -18.3155895 21.5292140

## All HPD intervals for each of coefficients include 0, so no evidence that preictor
## influences angular error
