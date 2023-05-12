### CITATION: This code is modified from scripts shared in Clay et al. 2020, J. Anim. Ecol.

### AIM: Validate best HMM using 3 methods: 1) comparison with expert classificaiton; 2) permutation of covariates

# 1. EXPERTLY VALIDATED TRACKS --------------------------------------------

# Define the packages
packages <- c("data.table", "momentuHMM", "dplyr", "lubridate", "caret", "raster",
              "rgdal", "rnaturalearth", "rnaturalearthdata", "ggpubr")

# Install packages not yet installed - change lib to library path
#installed_packages <- packages %in% rownames(installed.packages())

#if (any(installed_packages == FALSE)) {
#  install.packages(packages[!installed_packages], lib = "C:/Users/libraryPath")
#}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

#### Load validation data #### -------------------------------------------------

gps_manual <- read.csv("Data_inputs/WAAL_gps_2020-2021_manualStates.csv")

gps_manual$DateTime <- as.POSIXct(gps_manual$DateTime, format = "%Y-%m-%d %H:%M:%S")
colnames(gps_manual)[8] <- "state_manual"

# Remove duplicates
gps_manual <- gps_manual[!duplicated(gps_manual[c("Ring", "DateTime")]),]

#### Perform the validation #### -----------------------------------------------

## Covert variables to factors
gps_manual$State_HMM <- tolower(gps_manual$State_HMM)
gps_manual[,c("State_HMM", "behaviour")] <- lapply(gps_manual[,c("State_HMM", "behaviour")], as.factor)

## Sample sizes of each behaviour
nrow(subset(gps_manual, behaviour == "travel")) # 841
nrow(subset(gps_manual, behaviour == "search")) # 695
nrow(subset(gps_manual, behaviour == "rest")) # 791

## Compute overall confusion matrix
confusionMatrix(gps_manual$State_HMM, gps_manual$behaviour)

tab <- table(gps_manual$State_HMM, gps_manual$behaviour)

tab

#         rest search travel
# rest    621     46     33
# search  153    468    155
# travel   17    181    653

## Accuracy
tab[1,1]/colSums(tab)[1] # rest 0.7850822 
tab[2,2]/colSums(tab)[2] # search 0.6733813
tab[3,3]/colSums(tab)[3] # travel 0.7764566 

## Bold birds only
tab_bold <- table(gps_manual$State_HMM[gps_manual$boldCat == "bold"], 
                  gps_manual$behaviour[gps_manual$boldCat == "bold"])

tab_bold

#         rest search travel
# rest    324     14      8
# search   90    240     76
# travel    9    117    370

## Accuracy
tab_bold[1,1]/colSums(tab_bold)[1] # rest 0.7659574
tab_bold[2,2]/colSums(tab_bold)[2] # search 0.6469003 
tab_bold[3,3]/colSums(tab_bold)[3] # travel 0.814978


## Shy birds only

tab_shy <- table(gps_manual$State_HMM[gps_manual$boldCat == "shy"], 
                  gps_manual$behaviour[gps_manual$boldCat == "shy"])

tab_shy

#         rest search travel
# rest    297     32     25
# search   63    228     79
# travel    8     64    283

## Accuracy
tab_shy[1,1]/colSums(tab_shy)[1] # rest 0.8070652 
tab_shy[2,2]/colSums(tab_shy)[2] # search 0.7037037 
tab_shy[3,3]/colSums(tab_shy)[3] # travel 0.7312661 



#### Plot validating vs HMM data #### --------------------------------------------------

### Pre-processing for plotting 

# Load full GPS data for plotting
allgps <- read.csv("Data_inputs/WAAL_GPS_2010-2021_personality_wind.csv")
allgps$DateTime <- as.POSIXct(allgps$DateTime, 
                              format = "%Y-%m-%d %H:%M:%S",
                              tz = "UTC")

# Keep only rings where labelled data
rings.keep <- unique(gps_manual$Ring)
allgps <- subset(allgps, Ring %in% rings.keep)

# Merge in labelled states
gps_labels <- gps_manual[,c("Ring", "DateTime", "State_HMM", "behaviour")]
allgps <- merge(allgps, gps_labels, by = c("Ring", "DateTime"), all.x = T)
allgps$State_HMM <- ifelse(is.na(allgps$State_HMM), "unknown", 
                           as.character(allgps$State_HMM) )
allgps$behaviour <- ifelse(is.na(allgps$behaviour), "unknown", 
                           as.character(allgps$behaviour) )

## Set projections for plotting
proj.dec <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
proj.utm <- "+proj=laea +lat_0=-46.358639 +lon_0=51.706972 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Get world map
world <- ne_countries(scale = "medium", returnclass = "sf")
world2 <- sf::st_transform(world, crs = proj.utm)

# Set behaviour colours 
travel_col <- "#648FFF"
search_col <- "#DC267F"
rest_col <- "#FFB000"
unknown_col <- "#B1B3B8"

# Project data
allgps.proj <- allgps
coordinates(allgps.proj) <- ~ Longitude+Latitude
proj4string(allgps.proj) <- proj.dec
allgps.sp <- spTransform(allgps.proj, CRS(proj.utm))
allgps.t <- as.data.frame(allgps.sp)

# Isolate rings
mykits <- unique(allgps$Ring)

for (i in 1:length(mykits)){
  
  # Select the ring
  mygps <- subset(allgps.t, Ring == mykits[i])
  
  ## Build the manual plot
  manualPlot <- ggplot(data = world2) + 
    geom_sf(fill = "cadetblue", colour = "grey") +
    ggspatial::annotation_scale(location = "bl", width_hint = 0.25, style = "bar") +
    coord_sf(crs = proj.utm, xlim = c(max(mygps$Longitude)+1000, 
                                      min(mygps$Longitude)-1000), 
             ylim = c(min(mygps$Latitude)-1000,
                      max(mygps$Latitude)+1000),
             label_axes = list(top = "E", left = "N", bottom = "E", right = "N")) +
    geom_path(aes(x = Longitude, y = Latitude, colour = behaviour, group = "identity"), 
              size = 1.5, dat = mygps) +
    scale_colour_manual(values = c(rest_col, search_col, travel_col, unknown_col), 
                        labels = c("Rest", "Search", "Travel", "Unlabelled"), 
                        name = "Behaviour") +
    ggtitle(paste("Manual labels:", mygps$Ring[1])) +
    theme(panel.background = element_rect(fill = "white"), 
          panel.grid.major = element_line(colour = "grey80"),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text.x.bottom = element_blank(), axis.title.x.bottom = element_blank(),
          axis.text.y.right = element_blank(), axis.title.y.right = element_blank(),
          axis.title.y.left = element_blank(),
          legend.background = element_rect(fill = "transparent"),
          legend.key = element_blank(),
          legend.position = c(0.92, 0.12),
          plot.tag = element_text(size = 20),
          plot.title = element_text(size = 20, face = "bold"))
  
  assign(paste0("manual_", mygps$Ring[1]),
          manualPlot)
  
  ## Build the HMM plot
  hmmPlot <- ggplot(data = world2) + 
    geom_sf(fill = "cadetblue", colour = "grey") +
    ggspatial::annotation_scale(location = "bl", width_hint = 0.25, style = "bar") +
    coord_sf(crs = proj.utm, xlim = c(max(mygps$Longitude)+1000, 
                                      min(mygps$Longitude)-1000), 
             ylim = c(min(mygps$Latitude)-1000,
                      max(mygps$Latitude)+1000),
             label_axes = list(top = "E", left = "N", bottom = "E", right = "N")) +
    geom_path(aes(x = Longitude, y = Latitude, colour = State_HMM, group = "identity"), 
              size = 1.5, dat = mygps) +
    scale_colour_manual(values = c(rest_col, search_col, travel_col, unknown_col), 
                        labels = c("Rest", "Search", "Travel", "Unlabelled"), 
                        name = "Behaviour") +
    ggtitle(paste("HMM labels:", mygps$Ring[1])) +
    theme(panel.background = element_rect(fill = "white"), 
          panel.grid.major = element_line(colour = "grey80"),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text.x.bottom = element_blank(), axis.title.x.bottom = element_blank(),
          axis.text.y.right = element_blank(), axis.title.y.right = element_blank(),
          axis.title.y.left = element_blank(),
          legend.position = "none",
          plot.tag = element_text(size = 20),
          plot.title = element_text(size = 20, face = "bold")) 
  
  assign(paste0("hmm_", mygps$Ring[1]),
         hmmPlot)
  
}


# #### FIGURE S7 : Validation plots ---------------------------------------

## Choose 2 exemplar plots
manual_A <- manual_WAAL225 + theme(legend.position = c(0.1, 0.82)) + labs(tag = "(a)")
hmm_A <- hmm_WAAL225 + theme(legend.position = "none") + labs(tag = "(b)")
manual_B <- manual_WAAL165 + theme(legend.position = "none") + labs(tag = "(c)")
hmm_B <- hmm_WAAL165 + theme(legend.position = "none") + labs(tag = "(d)")

png(filename = "Figures/FIGS7_validation_plot.png", width = 12, height = 10, units = "in", res = 600)
ggarrange(manual_A, hmm_A,
          manual_B, hmm_B,
          ncol = 2, nrow = 2,
          align = "h")
dev.off()


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

for (i in (n.iter+1):n.iter2) {
  
  data.R <- mainDat
  data.R$WindDir <- sample(data.R$WindDir)
  
  write.csv(data.R, file = paste0("Data_outputs/data_randomised-", i, ".csv"), row.names = F)
  
}


### Randomise LoD

n.iter3 <- n.iter + n.iter2

for (i in (n.iter2+1):n.iter3) {
  
  data.R <- mainDat
  data.R$LoD <- sample(data.R$LoD)
  
  write.csv(data.R, file = paste0("Data_outputs/data_randomised-", i, ".csv"), row.names = F)
  
}



### Randomise boldness (mean_BLUP_logit) - needs to be randomised within trip
n.iter4 <-  n.iter + n.iter3

for (i in (n.iter3+1):n.iter4) {
  
  data.R <- mainDat
  
  # Randomise boldness within each individual
  ringLengths <- rle(as.character(data.R$Ring))$lengths # gives runs of same individual
  uniqueRing <- unique(data.R$Ring) # gives number of unique individuals
  boldSample <- sample(data.R$mean_BLUP_logit, size = length(uniqueRing)) # sample boldness estimates
  data.R$mean_BLUP_logit <- rep(boldSample, ringLengths) # repeat sampled estimate over run of individuals
  
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

library(ggplot2); library(momentuHMM)


# Load AIC data for randomised models -------------------------------------

# Males #
results.M_path <- "Data_outputs/HMM_validation_M-output/"
results.M_files <- list.files(results.M_path)
results.M_list <- vector("list", length = length(results.M_files))

for (i in 1:length(results.M_files)) {
  
  if (file.size(paste0(results.M_path, results.M_files[i])) == 0) next
  load(paste0(results.M_path, results.M_files[i]))
  results.M_list[[i]] <- AIC.p
}

results.M <- data.frame(do.call("rbind", results.M_list))

## Order by iteration number
results.M <- results.M[order(results.M$Iter),]

## Match to randomised covariate
covariates <- data.frame(cov = rep(c("WindSp", "WindDir", "LoD", "mean_BLUP_logit"), each = 50),
                         Iter = seq(1:200)-1)
results.M <- merge(results.M, covariates, by = "Iter", all.x = TRUE)
results.M$cov <- covariates

## Get AIC of best supported model
file.in <- paste0("./Data_outputs/", paste0("M_mod_", 4, ".RData"))
load(file = file.in)
m.best.mod <- mod.p

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
file.in <- paste0("./Data_outputs/", paste0("F_mod_", 9, ".RData"))
load(file = file.in)
f.best.mod <- model

AIC.F <- AIC(f.best.mod)

# Plot randomised AIC values relative to that of 'best' model -------------

maleplot <- ggplot(aes(x = cov, y = AIC), data = results.M) + geom_boxplot() +
  geom_hline(yintercept = AIC.M, linetype = "dashed", colour = "red", size = 1) +
  labs(x = "Randomised Covariates") +
  theme_bw() +
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        legend.position = "none",
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 16)) 


femaleplot <- ggplot(aes(x = cov, y = AIC), data = results.F) + geom_boxplot() +
  geom_hline(yintercept = AIC.F, linetype = "dashed", colour = "red", size = 1) +
  labs(x = "Randomised Covariates") +
  theme_bw() +
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        legend.position = "none",
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 16)) 

png(filename = "Figures/supp_AIC_comparison_FEMALE.png", width = 9, height = 7, units = "in", res = 600)
femaleplot
dev.off()

png(filename = "Figures/supp_AIC_comparison_MALE.png", width = 9, height = 7, units = "in", res = 600)
maleplot
dev.off()




