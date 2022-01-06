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



# 2. PERMUTATION - RUN MODELS ------------------------------------------------

## Produce datasets where one variable has been randomised and output as .RData file:
## this facilitates parallel processing of the models to speed up computation. This 
## section of the script should be run for each sex separately


# Produce randomised datasets ---------------------------------------------

## Each of the 4 covariates (WindSp, WindDir, LoD, mean_BLUP_logit) is randomised 
## in turn 50 times, giving 4 * 50 = 200 randomised datasets


### Load main dataset
mainDat <- read.csv("Data_inputs/WAAL_GPS_2010-2021_personality_wind.csv", stringsAsFactors = F)
mainDat <- subset(mainDat, Sex == "F")

### Randomise wind speed (WindSp) - outputs 1:50
n.iter <- 50

for (i in 1:n.iter) {
  
  data.R <- mainDat
  data.R$WindSp <- sample(data.R$WindSp)
  
  write.csv(data.R, file = paste0("Data_outputs/data_randomised-", i, ".csv"), row.names = F)
  
}

### Randomise wind direction (WindDir)

n.iter2 <- n.iter + n.iter

for (i in n.iter:n.iter2) {
  
  data.R <- mainDat
  data.R$WindDir <- sample(data.R$WindDir)
  
  write.csv(data.R, file = paste0("Data_outputs/data_randomised-", i, ".csv"), row.names = F)
  
}


### Randomise LoD

n.iter3 <- n.iter + n.iter2

for (i in n.iter2:n.iter3) {
  
  data.R <- mainDat
  data.R$LoD <- sample(data.R$LoD)
  
  write.csv(data.R, file = paste0("Data_outputs/data_randomised-", i, ".csv"), row.names = F)
  
}



### Randomise boldness (mean_BLUP_logit) - needs to be randomised within trip

for (i in n.iter3:n.iter4) {
  
  data.R <- mainDat
  
  # Randomise boldness within each trip
  tripLengths <- rle(as.character(data.R$ID))$lengths # gives runs of trips
  uniqueTrips <- unique(data.R$ID) # gives number of unique trips
  boldSample <- sample(data.R$mean_BLUP_logit, size = length(uniqueTrips)) # sample boldness estimates
  data.R$mean_BLUP_logit <- rep(boldSample, tripLengths) # repeat sampled estimate over run of trip
  
  write.csv(data.R, file = paste0("Data_outputs/data_randomised-", i, ".csv"), row.names = F)
  
}




# Run the randomised models -----------------------------------------------

## This section of the script is set-up to run on parallel machines. To run this in
## R directly, it can be inserted into a for loop to iterate through each dataset. 
## Please note this will require significant computing time. 

### Load model parameters from best model (for males and females separately)

# Load in best models for males and females
file.in <- paste0("./Data_outputs/", paste0("F_mod_", 8, ".RData"))
load(file = file.in)
f.mod <- model

file.in <- paste0("./Data_outputs/", paste0("M_mod_", 6, ".RData"))
load(file = file.in)
m.mod <- model

# Females #
form <- ~WindSp:mean_BLUP_logit + WindSp + mean_BLUP_logit + LoD + WindDir
par <- getPar0(model=f.mod, nbStates=3, formula = form)

# Males #
form <- ~WindSp:mean_BLUP_logit + WindSp + mean_BLUP_logit + LoD + WindDir + WindSp:mean_BLUP_logit:WindDir
par <- getPar0(model=m.mod, nbStates=3, formula = form)


## Assign step lengths, turning angles, and stateNames

shape_0 <- c(12.42, 4.10, 0.33)
scale_0 <- c(3.62, 4.71, 0.17)
stepPar0 <- c(shape_0,scale_0)

location_0 <- c(0.00302, 0.00343, 0.0291)
concentration_0 <- c(50.79, 1.27, 44.02)
anglePar0 <- c(location_0,concentration_0)

stateNames<-c("travel","search", "rest")


### Load randomised data and prepare for HMM

data.R <- read.csv("data_randomised-.csv")

## Prepare data for HMM
permDat <- momentuHMM::prepData(data.R, type= "LL", 
                                coordNames = c("x", "y"))

## Deal with zero lengths by setting 0 to very small numbers
ind_zero <- which(permDat$step == 0)
if (length(ind_zero)>0){
  permDat$step[ind_zero] <- runif(length(ind_zero))/10000
}
ind_zero <- which(permDat$step == "NA")
if (length(ind_zero)>0){
  permDat$step[ind_zero] <- runif(length(ind_zero))/10000
}

permDat <- subset(permDat, !is.na(angle))

### Run the model
mod.p <- momentuHMM::fitHMM(data = permDat, nbStates = 3,
                            dist = list(step = "gamma", angle = "vm"),
                            Par0 = list(step = par$Par$step, angle = par$Par$angle, delat0 = par$delta),
                            estAngleMean = list(angle = T), beta0 = par$beta, 
                            stateNames = stateNames, 
                            formula = form)

### Get AIC and save to file

AIC.p <- data.frame(Iter = i, Formula = as.character(form)[2], AIC = AIC(mod.p))
save(AIC.p, file = "results.RData")



# 2. PERMUTATION - RESULTS -------------------------------------------------------------

library(ggplot2)


# Load AIC data for randomised models -------------------------------------

# Males #
results.M_path <- "Data_outputs/HMM_validation_M-output/"
results.M_files <- list.files(results.M_path)
results.M_list <- vector("list", length = length(results.M_files))

for (i in 1:length(results.M_files)) {
  load(paste0(results.M_path, results.M_files[i]))
  results.M_list[[i]] <- AIC.p
}

results.M <- data.frame(do.call("rbind", results.M_list))

## Order by iteration number
results.M <- results.M[order(results.M$Iter),]

## Match to randomised covariate
covariates <- rep(c("WindSp", "WindDir", "LoD", "mean_BLUP_logit"), each = 50)
results.M$cov <- covariates

## Get AIC of best supported model
file.in <- paste0("./Data_outputs/", paste0("M_mod_", 6, ".RData"))
load(file = file.in)
m.best.mod <- model

AIC.M <- AIC(m.best.mod)


# Females #
results.F_path <- "Data_outputs/HMM_validation_F-output/"
results.F_files <- list.files(results.F_path)
results.F_list <- vector("list", length = length(results.F_files))

for (i in 1:length(results.F_files)) {
  load(paste0(results.F_path, results.F_files[i]))
  results.F_list[[i]] <- AIC.p
}

results.F <- data.frame(do.call("rbind", results.F_list))

## Order by iteration number
results.F <- results.F[order(results.F$Iter),]

## Match to randomised covariate
covariates <- rep(c("WindSp", "WindDir", "LoD", "mean_BLUP_logit"), each = 50)
results.F$cov <- covariates

## Get AIC of best supported model
file.in <- paste0("./Data_outputs/", paste0("F_mod_", 7, ".RData"))
load(file = file.in)
f.best.mod <- model

AIC.F <- AIC(f.best.mod)

# Plot randomised AIC values relative to that of 'best' model -------------

library(ggplot2)

maleplot <- ggplot(aes(x = cov, y = AIC), data = results.M) + geom_boxplot() +
  geom_hline(yintercept = AIC.M, linetype = "dashed", colour = "red", size = 1) +
  labs(x = "Randomised Covariates") +
  theme_bw()


femaleplot <- ggplot(aes(x = cov, y = AIC), data = results.F) + geom_boxplot() +
  geom_hline(yintercept = AIC.F, linetype = "dashed", colour = "red", size = 1) +
  labs(x = "Randomised Covariates") +
  theme_bw()

gridExtra::grid.arrange(femaleplot, maleplot, ncols = 2)
