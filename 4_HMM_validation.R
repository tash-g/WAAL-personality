### CITATION: This code is modified from scripts shared in Clay et al. 2020, J. Anim. Ecol.

### AIM: Validate best HMM using 3 methods: 1) comparison with expert classificaiton; 2) permutation of covariates

# 1. EXPERTLY VALIDATED TRACKS --------------------------------------------

library(momentuHMM); library(data.table)

#### Load best models and label states ####

## Load GPS data and get step lengths/turning angles for viterbi algorithm

all_data <- read.csv("./Data_inputs/WAAL_GPS_2010-2021_personality_wind.csv")

names(all_data)[10] <- "ID"
names(all_data)[8] <- "WindDir"
all_data$DateTime <- as.POSIXct(all_data$DateTime, format = "%Y-%m-%d %H:%M:%S")

all_data <- all_data[order(all_data$ID, all_data$DateTime),]

all_data[,c(1,2,3,9,10)] <- lapply(all_data[,c(1,2,3,9,10)], as.factor) # ring, sex, year, LoD, ID
all_data[,c(4,5,7,8,11)] <- lapply(all_data[,c(4,5,7,8,11)], as.numeric) # lon, lat, windsp, winddir, mean_blup_logit

all_data <- subset(all_data, !is.na(WindDir))

all_data.prepped <- prepData(all_data, type= "LL", 
                   coordNames = c("Longitude", "Latitude")) 

all_data.prepped <- subset(all_data.prepped, !is.na(angle))
all_data.prepped <- subset(all_data.prepped, step < 40)


# Males
file.in <- paste0("./Data_outputs/", paste0("M_mod_", 6, ".RData"))
load(file = file.in)
m.best.mod <- model

male_data <- subset(all_data.prepped, Sex == "M")
male_data$State <- momentuHMM::viterbi(m.best.mod)

male_data$State[male_data$State == 1] <- "Travel"
male_data$State[male_data$State == 2] <- "Search"
male_data$State[male_data$State == 3] <- "Rest"


# Females
file.in <- paste0("./Data_outputs/", paste0("F_mod_", 8, ".RData"))
load(file = file.in)
f.best.mod <- model

female_data <- subset(all_data.prepped, Sex == "F")
female_data$State <- momentuHMM::viterbi(f.best.mod)

female_data$State[female_data$State == 1] <- "Travel"
female_data$State[female_data$State == 2] <- "Search"
female_data$State[female_data$State == 3] <- "Rest"


## Bind together + extract relevant columns
gps_hmm <- rbind(female_data, male_data)
gps_hmm <- gps_hmm[,c(1, 4, 5, 6, 7, 12, 13, 14)] # TripID Ring Sex Year DateTime Longitude Latitude State
colnames(gps_hmm)[8] <- "State_HMM"

## Clear some space
all_data <- NULL
all_data.prepped <- NULL
male_data <- NULL
female_data <- NULL
model <- NULL
gc()


#### Load validation data ####

gps_manual <- read.csv("Data_inputs/WAAL_GPS_2020-2021_manualStates.csv")
gps_manual <- gps_manual %>% 
  mutate(Ring = as.integer(as.factor(Ring))) %>% data.frame()

gps_manual$BirdId <- NULL
gps_manual$Ring <- as.factor(gps_manual$Ring)
gps_manual$DateTime <- as.POSIXct(gps_manual$DateTime, format = "%d/%m/%Y %H:%M")
colnames(gps_manual)[c(6,7)] <- c("Longitude", "Latitude")
colnames(gps_manual)[8] <- "state_manual"

#### Perform the validation ####

## Extract relevant birds
valBirds <- unique(gps_manual$Ring)
gps_hmm.subset <- subset(gps_hmm, Ring %in% valBirds)

## Merge to nearest minute
setkey(setDT(gps_hmm.subset), "Ring", "DateTime")
setkey(setDT(gps_manual), "Ring", "DateTime")
gps_comp <- as.data.frame(gps_hmm.subset[gps_manual, roll = "nearest"])

gps_comp$state_manual[gps_comp$state_manual == 3] <- "Rest"
gps_comp$state_manual[gps_comp$state_manual == 2] <- "Search"
gps_comp$state_manual[gps_comp$state_manual == 1] <- "Travel"

gps_comp$validation <- ifelse(gps_comp$state_manual != gps_comp$State_HMM, 0, 1)
gps_comp <- subset(gps_comp, !is.na(validation))


### Compute accuracy

# Overall
sum(gps_comp$validation)/nrow(gps_comp)  # 78.6 %
sum(gps_comp$validation[gps_comp$State_HMM == "Travel"])/nrow(gps_comp[gps_comp$State_HMM == "Travel",]) # 88.9 %
sum(gps_comp$validation[gps_comp$State_HMM == "Search"])/nrow(gps_comp[gps_comp$State_HMM == "Search",]) # 49.8 %
sum(gps_comp$validation[gps_comp$State_HMM == "Rest"])/nrow(gps_comp[gps_comp$State_HMM == "Rest",]) # 94.6 %





# 2. PERMUTATION -------------------------------------------------------------

### Load in AIC data


#### Plot in comparison to real AIC ####

# AIC of best' model 
AIC.real_F <- AIC(f.best.mod)

# Load in tables for each covariate
windSp <- read.csv("./Data_outputs/AIC_table_F_randomised-WindSp.csv")
winddir <- read.csv("./Data_outputs/AIC_table_F_randomised-WindDir.csv")
lod <- read.csv("./Data_outputs/AIC_table_F_randomised-LoD.csv")
boldness <- read.csv("./Data_outputs/AIC_table_F_randomised-mean_BLUP_logit.csv")

windSp$Covar <- "WindSp"
winddir$Covar <- "WindDir"
lod$Covar <- "LoD"
boldness$Covar <- "mean_BLUP_logit"

merged.table <- rbind(windSp, winddir, lod,  boldness)

# Plot randomised AIC in comparison to 'best' model 
plot.df <- data.frame(AIC = rep(AIC.real_F), Covar <- c("WindSp", "WindDir", "LoD", "mean_BLUP_logit"))

ggplot(plotdf, aes(Covar, AIC)) + geom_boxplot()+geom_hline(yintercept=AIC.real_R, linetype="dashed", 
                                                               color = "red", size = 1.1)+xlab("Dummy covariates")


