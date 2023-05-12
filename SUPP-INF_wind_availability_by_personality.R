
### AIM: Examine whether wind conditions vary for birds of different personality 
### types and whether wind conditions changed over years

# Define the packages
packages <- c("lme4", "MuMIn", "ggplot2", "visreg", "plyr", "circular", "dplyr",
              "data.table", "bpnreg", "Rcpp")

# Install packages not yet installed - change lib to library path
#installed_packages <- packages %in% rownames(installed.packages())

#if (any(installed_packages == FALSE)) {
#  install.packages(packages[!installed_packages], lib = "C:/Users/libraryPath")
#}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))


# PREAMBLE ---------------------------------------------------

dat <- data.table::fread("Data_inputs/WAAL_GPS_2010-2021_personality_wind.csv")
dat$DateTime <- fasttime::fastPOSIXct(dat$DateTime)

## Transform ring sex and year to factor
dat[,c("Ring", "Sex")] <- lapply(dat[,c("Ring", "Sex")], as.character) 
dat[,c("Ring", "Sex")] <- lapply(dat[,c("Ring", "Sex")], as.factor) 

dat$Year <- as.numeric(as.character(dat$Year))

# MODEL WIND SPEED AVAILABILITY BY PERSONALITY -------------------------------------------------

### Look at whether wind availability is predicted by personality and sex

# Remove extremes of personality and rescale
dat.wind <- subset(dat, mean_BLUP_logit < 3 & mean_BLUP_logit > -3)
dat.wind$persScaled <- scales::rescale(dat.wind$mean_BLUP_logit, to = c(-3,3))

## Build global model: interested in interaction of sex and personality (as we known sex has an effect)
wind.lm <- lmerTest::lmer(WindSp ~ Sex*persScaled + Year + (1|Ring/Trip_bout), data = dat.wind, REML = FALSE)
options(na.action = "na.fail") 

## Compare all possible models 
m_set <- dredge(wind.lm)
m_set

#Model selection table 
# (Int) men_BLU_lgt Sex Yer men_BLU_lgt:Sex df   logLik    AICc  delta weight
#3   8.291            +                   5 -1184195 2368400  0.00  0.409
#7  51.640            + -0.02151          6 -1184195 2368402  1.35  0.209
#4   8.295 -0.01428   +                   6 -1184195 2368402  1.93  0.156
#12  8.312 -0.06744   +                +  7 -1184195 2368403  2.88  0.097
#8  52.530 -0.01666   + -0.02195          7 -1184195 2368404  3.25  0.081
#16 50.030 -0.06775   + -0.02070       +  8 -1184194 2368405  4.27  0.048
#1   8.727                                4 -1184214 2368436 35.29  0.000
#2   8.732 -0.06744                       5 -1184213 2368436 35.94  0.000
#5  54.600              -0.02276          5 -1184213 2368437 36.61  0.000
#6  57.990 -0.07000     -0.02445          6 -1184213 2368438 37.17  0.000


# Best model includes sex only
wind.lm.best <- lmerTest::lmer(WindSp ~ Sex + (1|Ring/Trip_bout), data = dat.wind, REML = FALSE)

plot(wind.lm.best)
plot(residuals(wind.lm.best, type = "pearson") ~ dat.wind$Sex)

summary(wind.lm.best)
emmeans::emmeans(wind.lm.best, pairwise~Sex)

# MODEL WIND DIRECTION AVAILABILITY BY PERSONALITY ---------------------------------------------------

### Get median wind direction per trip
dat$WindDir2 <- circular(dat$Dev.wind2,
                         type = "angles",
                         units = "degrees", 
                         zero = 0, 
                         rotation = "counter")

dat <- subset(dat, !is.na(WindDir2))
trips.sum <- ddply(dat,.(Sex, Trip_bout, Ring, mean_BLUP_logit, Year),
                   summarize, MedDir = median.circular(WindDir2))

hist(as.numeric(trips.sum$MedDir))
range(as.numeric(trips.sum$MedDir))


## Calculate mean wind direction for all individuals - mean and confidence intervals
mle.vonmises(trips.sum$MedDir)
# mu: 69.97  ( 0.5137 )

# kappa: 29.37  ( 1.983 )


## Bootstrap CIs
mle.vonmises.bootstrap.ci(trips.sum$MedDir, alpha = 0.05)
# Mean Direction:            Low = 68.97   High = 70.95  
# Concentration Parameter:   Low = 25.67   High = 34.2 

# function to convert to where wind is coming FROM on 0-360 scale
dir_conv <- function(x) { if (x > 0) { dir <- 360-(180-x)} else { dir <- 180+x }
  return(dir)
}

dir_conv(69.97) # 249.97
dir_conv(68.97) # 248.97
dir_conv(70.95) # 250.95


### Circular regression for personality

# Circular outcome in radians
rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}

trips.sum$radCirc <- deg2rad(trips.sum$MedDir)

# Bird ID needs to be a numeric index
trips.sum$Ring <- as.numeric(factor(trips.sum$Ring))
trips.mod <- trips.sum[,c("Ring", "Sex", "Year", "mean_BLUP_logit", "radCirc")]

## Build circular regression model: wind direction ~ personality ---------------------- 

system.time(fit.wind <- bpnme(pred.I = radCirc ~ mean_BLUP_logit + Sex + Year + (1|Ring), data = trips.mod,
                              its = 10000, burn = 1000, n.lag = 3, seed = 101))

beepr::beep(2)

print(fit.wind)

#Linear Coefficients 

#Component I: 
#                       mean         mode          sd        LB HPD      UB HPD
#(Intercept)     -34.71313040 -36.85687194 47.19047145 -127.88854507 56.76676294
#mean_BLUP_logit  -0.05572957  -0.04022718  0.07047042   -0.19213982  0.08661456
#SexM              0.03647552   0.05828197  0.17005781   -0.30191305  0.36985116
#Year              0.01819339   0.01918586  0.02342359   -0.02728164  0.06433881

# HPD interval = smallest interval in which 95% posterior distribution located
## HPD intervals for personality include 0, so on average no effect of personality for 
## experienced wind conditions

coef_circ(fit.wind)

## Important coefficients:
# bc = slope of circular regression line at inflection point
# AS = average slope over all x
# SAM = slope at mean of x

#                                mean       mode           sd      LB HPD     UB HPD
#mean_BLUP_logit bc    3.022932e-04  1.240079e-03 2.475284e-01 -4.676314e-02 5.192100e-02
#Year bc               1.062697e-02 -5.886142e-03 1.834695e+00 -1.946466e-01 1.334120e-01
#mean_BLUP_logit AS    3.296092e-05  2.687187e-07 1.032100e-02  1.417011e-11 1.068164e-03
#Year AS              -3.184135e-04  7.935310e-06 4.858307e-02  4.322466e-08 1.678932e-03
#mean_BLUP_logit SAM   9.527774e-05  2.687055e-07 6.125205e-03  1.417010e-11 1.065363e-03
#Year SAM              2.476959e-04  7.902993e-06 8.322723e-03  4.321891e-08 1.691018e-03
