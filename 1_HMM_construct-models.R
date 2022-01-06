### CITATION: This code is modified from scripts shared in Clay et al. 2020, J. Anim. Ecol.

### AIM: Run HMM combinations based on wind and personality predictors; use AIC to identify the best model


# Preamble ----------------------------------------------------------------

library(momentuHMM); library(ggplot2); library(dplyr)

### Create outputs folder
out.path <- "./Data_outputs/"
if(dir.exists(out.path) == FALSE){
  dir.create(out.path)
}


### LOAD IN TRACKS

gps <- read.csv(file = "./Data_inputs/WAAL_GPS_2010-2021_personality_wind.csv", na.strings = "NA")
names(gps)[10] <- "ID"
names(gps)[8] <- "WindDir"
gps <- gps[order(gps$ID, gps$DateTime),]
gps$DateTime <- NULL

## Isolate columns and specify types
gps[,c(1,2,3,8,9)] <- lapply(gps[,c(1,2,3,8,9)], as.factor) # ring, sex, year, LoD, ID
gps[,c(4:7,10)] <- lapply(gps[,c(4:7,10)], as.numeric) # Longitude, Latitude, WindSp, WindDir

# Some NA wind directions, because we don't have bearings for the last fix of the trip
gps <- subset(gps, !is.na(WindDir))


# INITIALISE HMM DATA ---------------------------------------------------------


dat_OG <- prepData(gps, type= "LL", # longs and lats
               coordNames = c("Longitude", "Latitude")) ## these are our names
head(dat_OG)


# Remove erroneous step lengths that can't be larger than 40 m
hist(dat_OG$step)
nrow(subset(dat_OG, step > 40))/nrow(dat_OG)
dat_OG <- subset(dat_OG, step < 40)

### NOTE: Initial value finding code taken from Clay et al. 2020

## Assign step lengths based on previously identified initial values
shape_0 <- c(12.42, 4.10, 0.33)
scale_0 <- c(3.62, 4.71, 0.17)
stepPar0 <- c(shape_0,scale_0)

## Assigning turning angles based on previously identified initial values
location_0 <- c(0.00302, 0.00343, 0.0291)
concentration_0 <- c(50.79, 1.27, 44.02)
anglePar0 <- c(location_0,concentration_0)


### Set 0s to very small numbers 
dat <- dat_OG

ind_zero <- which(dat$step == 0)
if (length(ind_zero)>0){
  dat$step[ind_zero] <- runif(length(ind_zero))/10000
}
ind_zero <- which(dat$step == "NA")
if (length(ind_zero)>0){
  dat$step[ind_zero] <- runif(length(ind_zero))/10000
}

dat <- na.omit(dat)
dat$LoD <- as.factor(dat$LoD)
dat <- subset(dat, !is.na(angle))
dat <- subset(dat, !is.na(step))

# RUN THE NULL MODEL ------------------------------------------------------

#  first run null models with no covariates on transition probabilities

stateNames <- c("travel","search", "rest")

dat_male <- subset(dat, Sex == "M")
dat_female <- subset(dat, Sex == "F")

m1_M <- fitHMM(data=dat_male, nbStates=3,
             dist=list(step="gamma",angle="vm"),
             Par0=list(step=stepPar0, angle=anglePar0),
             estAngleMean = list(angle=TRUE),
             stateNames = stateNames)

m1_F <- fitHMM(data=dat_female, nbStates=3,
               dist=list(step="gamma",angle="vm"),
               Par0=list(step=stepPar0, angle=anglePar0),
               estAngleMean = list(angle=TRUE),
               stateNames = stateNames)

save(m1_M, file = "Data_outputs/M_mod_1.RData")
save(m1_F, file = "Data_outputs/F_mod_1.RData")

# Plot pseudo-residuals
plotPR(m1_M)
plotPR(m1_F)

# Look at acf of step length over longer time lag
acf(dat$step[!is.na(dat$step)], lag.max = 300)


#load("Data_outputs/M_mod_1.RData")
#load("Data_outputs/F_mod_1.RData")

# SET UP FORMULA FOR ALL MODEL COMBINATIONS -------------------------------

## Run next couple of models for males and females separately - make sure to change model names as appropriate

### CURRENTLY RUNNING FOR FEMALES ###

model_run <- m1_F # will have to manually change for males to m1_M
data_run <- dat_female # will have to manually change for males to dat_male
sex_initial <- "F"

## Set up personality formulae 

formula <- list()	
formula[[2]] <- ~ WindSp:mean_BLUP_logit + WindSp + mean_BLUP_logit + LoD + WindDir + WindSp:mean_BLUP_logit:WindDir + WindSp:WindDir
formula[[3]] <- ~ WindSp:mean_BLUP_logit + WindSp + LoD + WindDir + WindSp:mean_BLUP_logit:WindDir + WindSp:WindDir
formula[[4]] <- ~ WindSp:mean_BLUP_logit + WindSp + mean_BLUP_logit + LoD + WindDir + WindSp:WindDir
formula[[5]] <- ~ WindSp:mean_BLUP_logit + WindSp + LoD + WindDir + WindSp:WindDir
formula[[6]] <- ~ WindSp:mean_BLUP_logit + WindSp + mean_BLUP_logit + LoD + WindDir + WindSp:mean_BLUP_logit:WindDir
formula[[7]] <- ~ WindSp:mean_BLUP_logit + WindSp + LoD + WindDir + WindSp:mean_BLUP_logit:WindDir
formula[[8]] <- ~ WindSp:mean_BLUP_logit + WindSp + mean_BLUP_logit + LoD + WindDir
formula[[9]] <- ~ WindSp:mean_BLUP_logit + WindSp + LoD + WindDir
formula[[10]] <- ~ WindSp + mean_BLUP_logit + LoD + WindDir + WindSp:mean_BLUP_logit:WindDir + WindSp:WindDir
formula[[11]] <- ~ WindSp + LoD + WindDir + WindSp:mean_BLUP_logit:WindDir + WindSp:WindDir
formula[[12]] <- ~ WindSp + mean_BLUP_logit + LoD + WindDir + WindSp:WindDir
formula[[13]] <- ~ WindSp + LoD + WindDir + WindSp:WindDir
formula[[14]] <- ~ WindSp + mean_BLUP_logit + LoD + WindDir
formula[[15]] <- ~ WindSp + LoD + WindDir
formula[[16]] <- ~ WindSp + mean_BLUP_logit + LoD + WindDir + WindSp:mean_BLUP_logit:WindDir
formula[[17]] <- ~ WindSp + LoD + WindDir + WindSp:mean_BLUP_logit:WindDir
formula[[18]] <- ~ WindSp:mean_BLUP_logit + WindSp + mean_BLUP_logit + LoD + WindDir + WindDir:mean_BLUP_logit + WindSp:mean_BLUP_logit:WindDir + WindSp:WindDir
formula[[19]] <- ~ WindSp + mean_BLUP_logit + LoD + WindDir + WindDir:mean_BLUP_logit + WindSp:mean_BLUP_logit:WindDir + WindSp:WindDir
formula[[20]] <- ~ WindSp:mean_BLUP_logit + WindSp + mean_BLUP_logit + LoD + WindDir + WindDir:mean_BLUP_logit + WindSp:WindDir
formula[[21]] <- ~ WindSp + mean_BLUP_logit + LoD + WindDir + WindDir:mean_BLUP_logit + WindSp:WindDir


# this function gets starting values for each model from existing null model fit for a specified covariate formula
Par <- list()
for (i in 2:length(formula)){
  Par[[i]] <- getPar0(model=model_run, nbStates=3, formula = formula[[i]])
}



# RUN ALL THE MODELS ------------------------------------------------------

## the following code iterates through and runs all 40 models, pasting each out in turn

stateNames<-c("travel","search", "rest")

# Specify output lists of length formula
m.list <- vector(mode = "list", length =length(formula))

for (i in 2:length(formula)) {
    print(i)
    m.list[[i]] <- fitHMM(data=data_run, nbStates=3,
                            dist=list(step="gamma",angle="vm"),
                            Par0=list(step=Par[[i]]$Par$step, angle=Par[[i]]$Par$angle,delta0 = Par[[i]]$delta),
                            estAngleMean = list(angle=TRUE), beta0 = Par[[i]]$beta,
                            stateNames = stateNames, 
                            formula = formula[[i]])
    model <- m.list[[i]]
    file.out <- paste0("./Data_outputs/", paste0(sex_initial, "_mod_", i, ".RData"))
    save(model, file = file.out)
   }

# FIND THE BEST MODEL -----------------------------------------------------

# Iterate through each model set, load in models, output autocorrelation plots for step lengths and turning angles, pseudo-residual 
# and qq plots, extract model coefficients for each transition and plot, and paste out AIC table

### Create figures folder
fig.path <- "./Figures/"
if(dir.exists(fig.path) == FALSE){
  dir.create(fig.path)
}

# Specifying output lists of length formula
m.list <- vector(mode = "list", length =length(formula))
out.df <- vector(mode = "list", length =length(formula))

for (i in 1:length(formula)) {
  print(i)
  
  file.in <- paste0("./Data_outputs/", paste0(sex_initial,"_mod_", i, ".RData"))
  load(file = file.in)
  if (i == 1) { m.list[[i]] <- m1_F} else { m.list[[i]] <- model}
  
  ## Extract AIC of model
  if (i == 1) { form_out <- 1} else { form_out <- as.character(formula[[i]])[2]}
  out.df[[i]] <- data.frame(Model = paste0("F_M", i),
                            Formula = form_out, AIC = AIC(m.list[[i]]))
  
  ## Construct autocorrelation plot
  
  pr <- pseudoRes(m.list[[i]])
  
  # Step length
  par(mfrow=c(1,1))
  ylimit <- qnorm((1 + 0.95)/2)/sqrt(length(pr$stepRes[!is.na(pr$stepRes)])) + 1
  acf(pr$stepRes[is.finite(pr$stepRes)],lag.max = 300)
  name.plot <- paste0("./Figures/", paste0(sex_initial, "_mod_", i, "_acf_step.png"))
  dev.copy(png, name.plot,  width = 800, height = 500)
  dev.off()
  
  # Turning angle
  par(mfrow=c(1,1))
  acf(pr$angleRes[!is.na(pr$angleRes)],lag.max = 300)
  name.plot <- paste0("./Figures/", paste0(sex_initial, "_mod_", i, "_acf_angle.png"))
  dev.copy(png, name.plot, width = 800, height = 500)
  dev.off()
  
  # outputting pseudo residual plots
  res.df <- data.frame(stepres = pr$stepRes, angleres = pr$angleRes,
                       sex = m.list[[i]]$data$Sex, ws = m.list[[i]]$data$WindSp, 
                       lod = m.list[[i]]$data$LoD)
  res.df <- subset(res.df, is.finite(stepres))
  
  par(mfrow=c(2,2))
  boxplot(stepres ~ lod, data = res.df)
  plot(stepres ~ ws, data = res.df)
  boxplot(angleres ~ lod, data = res.df)
  plot(angleres ~ ws, data = res.df)
  name.plot <- paste0("./Figures/", paste0(sex_initial, "_mod_", i, "_pseudo-residuals.png"))
  dev.copy(png, name.plot, width = 800, height = 500)
  dev.off()
  
  # outputting qq plots
  plotPR(m.list[[i]])
  name.plot <- paste0("./Figures/", paste0(sex_initial, "_mod_", i, "_acf_qq.png"))
  dev.copy(png, name.plot, width = 800, height = 500)
  dev.off() 
  
  # extracting and plotting model coefficients for each transition
  beta.full <- CIbeta(m.list[[i]])$beta
  beta.full.est <- as.data.frame(beta.full$est)
  beta.full.upr <- as.data.frame(beta.full$upper)
  beta.full.lwr <- as.data.frame(beta.full$lower)
  beta.df <- data.frame(Est = c(beta.full.est$`1 -> 2`, beta.full.est$`1 -> 3`,
                                beta.full.est$`2 -> 1`, beta.full.est$`2 -> 3`,
                                beta.full.est$`3 -> 1`, beta.full.est$`3 -> 2`), 
                        Upr = c(beta.full.upr$`1 -> 2`, beta.full.upr$`1 -> 3`,
                                beta.full.upr$`2 -> 1`, beta.full.upr$`2 -> 3`,
                                beta.full.upr$`3 -> 1`, beta.full.upr$`3 -> 2`), 
                        Lwr = c(beta.full.lwr$`1 -> 2`, beta.full.lwr$`1 -> 3`,
                                beta.full.lwr$`2 -> 1`, beta.full.lwr$`2 -> 3`,
                                beta.full.lwr$`3 -> 1`, beta.full.lwr$`3 -> 2`), 
                        Transitions = rep(colnames(beta.full.est), each = nrow(beta.full.est)),
                        Covariates = rep(rownames(beta.full.est), 3))
  beta.df$Covariates <- as.factor(as.character(beta.df$Covariates))
  beta.df$Transitions <- as.factor(as.character(beta.df$Transitions))
  pd <- position_dodge(width=0.7)
  
  ## Removing intercept to plot
  beta.df2 <- subset(beta.df, Covariates != "(Intercept)",)
  
  pl <- ggplot(beta.df2, aes(Covariates, Est)) + geom_hline(yintercept=0, linetype="dashed", size=1)+
    geom_point(aes(colour = Transitions),position=pd)+
    geom_errorbar(aes(ymin=Lwr, ymax=Upr, colour = Transitions), width=.8, position=pd) +theme_bw()
  print(pl)
  name.plot <- paste0("./Figures/", paste0(sex_initial, "_mod_", i, "_coefficients.png"))
  dev.copy(png, name.plot, width = 800, height = 500)
  dev.off()
}

# Warning messages appear ("removed containing missing values") due to upper and lower CIs which are sometimes "NA"

all.out <- do.call(rbind, out.df)
all.out <- all.out[order(all.out$AIC),]

# Print out AIC table
out.path <- paste0("./Data_outputs/AIC_table_", sex_initial, ".csv")
write.csv(all.out, out.path, row.names=T)



