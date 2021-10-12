library(ggplot2); library(momentuHMM)
source("Functions/gps_functions.R")

ci.stationary <- function(model, cov, alpha) {
  gridLength <- nrow(cov)
  model2 <- momentuHMM:::delta_bc(model) # extended format allowing for backwards compatibility
  nbStates <- length(model2$stateNames)
  Sigma <- model2$mod$Sigma
  
  lci <- matrix(NA,gridLength,nbStates)
  uci <- matrix(NA,gridLength,nbStates)
  
  formula <- model2$conditions$formula # for covariates
  newForm <- momentuHMM:::newFormulas(formula, nbStates) # a formula for each state transition
  newformula <- newForm$newformula
  
  nbCovs <- ncol(model.matrix(newformula, model2$data)) - 1 # model.matrix gives values of covariate terms of the formula for each observation
  # in my opinion, nbCovs gives exactly the same number of terms as formula (-1 is minus intercept)
  
  gamInd <- (length(model2$mod$estimate) - (nbCovs + 1) * nbStates * (nbStates - 1) 
             + 1):(length(model2$mod$estimate)) - ncol(model2$covsDelta) * (nbStates - 1) * (!model2$conditions$stationary) # if the mode is stationary there wouldn't be any covariates
  # here we just enumerate parameters leaving out the first ones which are step and angle related, and the last ones which are delta related
  
  rawCovs <- model2$rawCovs
  # changing order of variables tomatch with raw data
  rawCovs <- rawCovs[,order(names(rawCovs))]
  tempCovs <- cov
  tempCovs <- tempCovs[,sort(names(rawCovs))]
  
  # fixing format problems
  for (i in which(unlist(lapply(rawCovs, is.factor)))) {
    tempCovs[[i]] <- factor(tempCovs[[i]], levels = levels(rawCovs[, 
                                                                   i])) 
  }
  tmpSplineInputs <- momentuHMM:::getSplineFormula(newformula, model2$data, 
                                                   tempCovs) # just a format for spline use later
  desMat <- model.matrix(tmpSplineInputs$formula, data = tmpSplineInputs$covs) # tmpSplineInputs$covs is tempCovs
  # model.matrix gives the design matrix for the formula in the argument
  
  probs <- as.data.frame(stationary(model2, covs=desMat)) # the stationary probability is computed based on desMat, which has only a range of values for one of the covariates
  
  for(state in 1:nbStates) {
    dN <- t(apply(desMat, 1, function(x)
      numDeriv::grad(momentuHMM:::get_stat,model2$mod$estimate[gamInd][unique(c(model2$conditions$betaCons))],covs=matrix(x,nrow=1),nbStates=nbStates,i=state,betaRef=model2$conditions$betaRef,betaCons=model2$conditions$betaCons,workBounds=model2$conditions$workBounds$beta)))
    
    se <- t(apply(dN, 1, function(x)
      suppressWarnings(sqrt(x%*%Sigma[gamInd[unique(c(model2$conditions$betaCons))],gamInd[unique(c(model2$conditions$betaCons))]]%*%x))))
    
    lci[,state] <- 1/(1 + exp(-(log(probs[,state]/(1-probs[,state])) -
                                  qnorm(1-(1-alpha)/2) * (1/(probs[,state]-probs[,state]^2)) * se)))
    uci[,state] <- 1/(1 + exp(-(log(probs[,state]/(1-probs[,state])) +
                                  qnorm(1-(1-alpha)/2) * (1/(probs[,state]-probs[,state]^2)) * se)))
  }
  lci <- as.data.frame(lci)
  names(lci) <- model$stateNames
  uci <- as.data.frame(uci)
  names(uci) <- model$stateNames
  ci.list <- list(lci, uci)
  names(ci.list) <- c("lower", "upper")
  return(ci.list)
}

# AIM OF SCRIPT - to plot predicted stationary probabilities for a range of specified covariates (essentially gives
# time-activity budgets)

# LOAD MODEL OUTPUTS ------------------------------------------------------

# load in best models for males and females
file.in <- paste0("./Data_outputs/", paste0("F_mod_", 8, ".RData"))
load(file = file.in)
f.mod <- model

file.in <- paste0("./Data_outputs/", paste0("M_mod_", 6, ".RData"))
load(file = file.in)
m.mod <- model


##### 1. PLOTTING STATIONARY PROBABILITIES FOR COVARIATE VALUES OF WIND SPEED ###############

min(f.mod$data$WindSp) # 0.05
max(f.mod$data$WindSp) # 20.7

min(m.mod$data$WindSp) # 0.05
max(m.mod$data$WindSp) # 24.2

## Predict wind speed for average wind direction which is 75 degrees relative to bird direction, and for daylight hours when birds are most active. 

# Get covariate data
ws.f <- seq(from=0.05, 20, by = 0.2)
pers.f <- c(round(min(f.mod$data$mean_BLUP_logit), digits = 1), 
            round(max(f.mod$data$mean_BLUP_logit), digits = 1))
comb.f <- expand.grid(ws.f, pers.f)

ws.m <- seq(from=0.05, 24, by = 0.2)
pers.m <- c(round(min(m.mod$data$mean_BLUP_logit), digits = 1), 
            round(max(m.mod$data$mean_BLUP_logit), digits = 1))
comb.m <- expand.grid(ws.m, pers.m)

cov.f <- data.frame(LoD = rep("L"), WindDir = rep(75), WindSp = comb.f$Var1, mean_BLUP_logit = comb.f$Var2)
cov.m <- data.frame(LoD = rep("L"), WindDir = rep(75), WindSp = comb.m$Var1, mean_BLUP_logit = comb.m$Var2)
cov.f$LoD <- as.factor(cov.f$LoD)
cov.m$LoD <- as.factor(cov.m$LoD)

# Get stationary probs (gives for each row of covariates dataset)
s.f <- as.data.frame(stationary(model = f.mod, covs = cov.f))
s.m <- as.data.frame(stationary(model = m.mod, covs = cov.m))

# Create confidence intervals around stationary probabilities
stat.f <- ci.stationary(model = f.mod, cov = cov.f, alpha = 0.95) 
stat.m <- ci.stationary(model = m.mod, cov = cov.m, alpha = 0.95) 


# Combine into dataframe for plotting

prob_cov.f <- cbind(s.f, cov.f)
prob_cov.m <- cbind(s.m, cov.m)

prob.f <- tidyr::gather(prob_cov.f, state, prob, travel:rest, factor_key=TRUE)
prob.m <- tidyr::gather(prob_cov.m, state, prob, travel:rest, factor_key=TRUE)

stat.df <- rbind(prob.f, prob.m)

lwr.f <- tidyr::gather(stat.f[[1]], state, lwr, travel:rest, factor_key = T)$lwr
lwr.m <- tidyr::gather(stat.m[[1]], state, lwr, travel:rest, factor_key = T)$lwr

lwr <- c(lwr.f, lwr.m)

upr.f <- tidyr::gather(stat.f[[2]], state, lwr, travel:rest, factor_key = T)$lwr
upr.m <- tidyr::gather(stat.m[[2]], state, lwr, travel:rest, factor_key = T)$lwr

upr <- c(upr.f, upr.m)

stat.df$lwr <- lwr
stat.df$upr <- upr
stat.df$sex <- c(rep("F", nrow(prob.f)), rep("M", nrow(prob.m)))
stat.df$pers <- ifelse(stat.df$mean_BLUP_logit < -1.5, "shy", "bold")

stat.df[,c(9,10)] <- lapply(stat.df[,c(9,10)], as.factor)

stat.df$pers_state <- as.factor(as.character(paste(stat.df$pers, stat.df$state, sep = "_")))

# create plot

shy_col <- "#00DD2F"
bold_col <- "purple"

stat.df$sex <- factor(stat.df$sex, labels = c("Female", "Male"))

png("Figures/Transition Estimates/Stationary estimates/stationary_transitions_SPEED.png", width = 9, height = 7, units = "in", res = 600)
ggplot(stat.df, aes(WindSp, prob)) + facet_wrap(.~sex, labeller=label_value) + 
  geom_ribbon(size = 1.3, aes(ymin=lwr, ymax=upr, fill = pers_state), alpha=0.15)+
  scale_fill_manual(values = c("grey0", "grey0", "grey0", "grey30", "grey30", "grey30", "grey60", "grey60", "grey60")) +
  geom_line(aes(colour = state, linetype = pers), size = 1.3) + 
  scale_linetype_manual(values=c("solid", "dotted"),labels = c("Bold", "Shy"), name = "Personality") +
  ylim(0, 1) +
  scale_colour_viridis_d(labels = c("Travel", "Search", "Rest"), name = "Behaviour") +
  theme_bw() + #ylab("Stationary probability")+
  scale_x_continuous(limits=c(0, 23)) + xlab(expression(paste("Wind speed (", ms^-1, ")", sep="")))+
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.x=element_text(size=18),
        axis.title.y=element_blank(),
        legend.position = c(0.07,0.77),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 16)) +
  guides(fill = FALSE)
dev.off()




##### 2. PLOTTING STATIONARY PROBABILITIES FOR COVARIATE VALUES OF WIND DIRECTION ###############


min(f.mod$data$WindDir) # 0.001
max(f.mod$data$WindDir) # 180

min(m.mod$data$WindDir) # 0.001
max(m.mod$data$WindDir) # 180

## Predict wind speed for average wind speed which is 9 ms-1, and for daylight hours when birds are most active. 

# Get covariate data
wdir.f <- seq(from=0, 180, by = 1)
pers.f <- c(round(min(f.mod$data$mean_BLUP_logit), digits = 1), 
            round(max(f.mod$data$mean_BLUP_logit), digits = 1))
comb.f <- expand.grid(wdir.f, pers.f)

wdir.m <- seq(from=0, 180, by = 1)
pers.m <- c(round(min(m.mod$data$mean_BLUP_logit), digits = 1), 
            round(max(m.mod$data$mean_BLUP_logit), digits = 1))
comb.m <- expand.grid(wdir.f, pers.m)

cov.f <- data.frame(LoD = rep("L"), WindDir = wdir.f, WindSp = rep(9), mean_BLUP_logit = comb.f$Var2)
cov.m <- data.frame(LoD = rep("L"), WindDir = wdir.f, WindSp = rep(9), mean_BLUP_logit = comb.m$Var2)

# Get stationary probs (gives for each row of covariates dataset)
s.f <- as.data.frame(stationary(model = f.mod, covs = cov.f))
s.m <- as.data.frame(stationary(model = m.mod, covs = cov.m))

# Create confidence intervals around stationary probabilities
stat.f <- ci.stationary(model = f.mod, cov = cov.f, alpha = 0.95) 
stat.m <- ci.stationary(model = m.mod, cov = cov.m, alpha = 0.95) 


# Combine into dataframe for plotting

prob_cov.f <- cbind(s.f, cov.f)
prob_cov.m <- cbind(s.m, cov.m)
prob.f <- tidyr::gather(prob_cov.f, state, prob, travel:rest, factor_key=TRUE)
prob.m <- tidyr::gather(prob_cov.m, state, prob, travel:rest, factor_key=TRUE)

stat.df <- rbind(prob.f, prob.m)

lwr.f <- tidyr::gather(stat.f[[1]], state, lwr, travel:rest, factor_key = T)$lwr
lwr.m <- tidyr::gather(stat.m[[1]], state, lwr, travel:rest, factor_key = T)$lwr
lwr <- c(lwr.f, lwr.m)

upr.f <- tidyr::gather(stat.f[[2]], state, lwr, travel:rest, factor_key = T)$lwr
upr.m <- tidyr::gather(stat.m[[2]], state, lwr, travel:rest, factor_key = T)$lwr
upr <- c(upr.f, upr.m)

stat.df$lwr <- lwr
stat.df$upr <- upr
stat.df$sex <- c(rep("F", nrow(prob.f)), rep("M", nrow(prob.m)))
stat.df$pers <- ifelse(stat.df$mean_BLUP_logit < -1.5, "shy", "bold")

stat.df[,c(9,10)] <- lapply(stat.df[,c(9,10)], as.factor)
stat.df$pers_state <- as.factor(as.character(paste(stat.df$pers, stat.df$state, sep = "_")))

# create plot

shy_col <- "purple"
bold_col <- "#00DD2F"
mid_col <- "#a099a6"

png("Figures/Transition Estimates/Stationary estimates/stationary_transitions_DIR.png", width = 9, height = 7, units = "in", res = 350)
ggplot(stat.df, aes(WindDir, prob)) + facet_wrap(~sex, labeller=labeller(sex=labels)) + 
  geom_ribbon(size = 1.3, aes(ymin=lwr, ymax=upr, fill = pers_state), alpha=0.15)+
  scale_fill_manual(values = c("grey0", "grey0", "grey0", "grey30", "grey30", "grey30", "grey60", "grey60", "grey60")) +
  geom_line(aes(colour = state, linetype = pers), size = 1.3) + 
  scale_linetype_manual(values=c("solid", "dotted"), labels = c("Bold", "Shy"), name = "Personality") +
  scale_colour_viridis_d(labels = c("Travel", "Search", "Rest"), name = "Behaviour") +
  ylim(0, 1)+
  theme_bw() + ylab("Stationary probability")+
  scale_x_continuous(limits=c(0, 180)) + xlab("Wind direction (°)")+
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=18),
        legend.position = c(0.07,0.77),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 16)) +
  guides(fill = FALSE)
dev.off()


# CALCULATE STATIONARY PROBABILITIES FOR EACH WIND DIRECTION -------------------


# Get covariate data
ws.f <- seq(from=0.05, 20, by = 0.2)
wd.f <- c(0, 90, 180)
pers.f <- c(round(min(f.mod$data$mean_BLUP_logit), digits = 1), 
            round(max(f.mod$data$mean_BLUP_logit), digits = 1))
comb.f <- expand.grid(ws.f, pers.f, wd.f)

ws.m <- seq(from=0.05, 24, by = 0.2)
wd.m <- c(0, 90, 180)
pers.m <- c(round(min(m.mod$data$mean_BLUP_logit), digits = 1), 
            round(max(m.mod$data$mean_BLUP_logit), digits = 1))
comb.m <- expand.grid(ws.m, pers.m, wd.m)

cov.f <- data.frame(LoD = rep("L"), WindSp = comb.f$Var1, WindDir = comb.f$Var3, mean_BLUP_logit = comb.f$Var2)
cov.m <- data.frame(LoD = rep("L"), WindSp = comb.m$Var1, WindDir = comb.m$Var3, mean_BLUP_logit = comb.m$Var2)
cov.f$LoD <- as.factor(cov.f$LoD)
cov.m$LoD <- as.factor(cov.m$LoD)

# Get stationary probs (gives for each row of covariates dataset)
s.f <- as.data.frame(stationary(model = f.mod, covs = cov.f))
s.m <- as.data.frame(stationary(model = m.mod, covs = cov.m))

# Create confidence intervals around stationary probabilities
stat.f <- ci.stationary(model = f.mod, cov = cov.f, alpha = 0.95) 
stat.m <- ci.stationary(model = m.mod, cov = cov.m, alpha = 0.95) 


# Combine into dataframe for plotting

prob_cov.f <- cbind(s.f, cov.f)
prob_cov.m <- cbind(s.m, cov.m)

prob.f <- tidyr::gather(prob_cov.f, state, prob, travel:rest, factor_key=TRUE)
prob.m <- tidyr::gather(prob_cov.m, state, prob, travel:rest, factor_key=TRUE)

stat.df <- rbind(prob.f, prob.m)

lwr.f <- tidyr::gather(stat.f[[1]], state, lwr, travel:rest, factor_key = T)$lwr
lwr.m <- tidyr::gather(stat.m[[1]], state, lwr, travel:rest, factor_key = T)$lwr

lwr <- c(lwr.f, lwr.m)

upr.f <- tidyr::gather(stat.f[[2]], state, lwr, travel:rest, factor_key = T)$lwr
upr.m <- tidyr::gather(stat.m[[2]], state, lwr, travel:rest, factor_key = T)$lwr

upr <- c(upr.f, upr.m)

stat.df$lwr <- lwr
stat.df$upr <- upr
stat.df$sex <- c(rep("F", nrow(prob.f)), rep("M", nrow(prob.m)))
stat.df$pers <- ifelse(stat.df$mean_BLUP_logit < -1.5, "shy", "bold")

stat.df[,c(9,10)] <- lapply(stat.df[,c(9,10)], as.factor)

stat.df$pers_state <- as.factor(as.character(paste(stat.df$pers, stat.df$state, sep = "_")))

## Stats for manuscript

# Get switchpoint where travel > search for shy birds
# Crosswind
x <- subset(stat.df, state == "search" & pers == "shy" & sex == "M" & WindDir == 90)$prob
y <- subset(stat.df, state == "travel" & pers == "shy" & sex == "M" & WindDir == 90)$prob

stat.df$WindSp[which(x-y<0)][1] # 13.85

## Females = 9.25

# Headwind
x <- subset(stat.df, state == "search" & pers == "shy" & sex == "M" & WindDir == 180)$prob
y <- subset(stat.df, state == "travel" & pers == "shy" & sex == "M" & WindDir == 180)$prob

stat.df$WindSp[which(x-y<0)][1] # 9.45

## Females = 7.85

# create plot

shy_col <- "#00DD2F"
bold_col <- "purple"

stat.df$sex <- factor(stat.df$sex, labels = c("Female", "Male"))
dirs <- unique(stat.df$WindDir)
dirLabels <- c("tail", "cross", "head")

for (i in 1:length(dirs)){
  
  plot.df <- subset(stat.df, WindDir == dirs[i])
  
  myplot <- ggplot(plot.df, aes(WindSp, prob)) + facet_wrap(.~sex, labeller=label_value) + 
    geom_ribbon(size = 1.3, aes(ymin=lwr, ymax=upr, fill = pers_state), alpha=0.15)+
    scale_fill_manual(values = c("grey0", "grey0", "grey0", "grey30", "grey30", "grey30", "grey60", "grey60", "grey60")) +
    geom_line(aes(colour = state, linetype = pers), size = 1.3) + 
    scale_linetype_manual(values=c("solid", "dotted"),labels = c("Bold", "Shy"), name = "Personality") +
    ylim(0, 1) +
    scale_colour_viridis_d(labels = c("Travel", "Search", "Rest"), name = "Behaviour") +
    theme_bw() + #ylab("Stationary probability")+
    scale_x_continuous(limits=c(0, 23)) + xlab(expression(paste("Wind speed (", ms^-1, ")", sep="")))+
    theme(axis.text.x=element_text(size=16), 
          axis.text.y=element_text(size=16), 
          axis.title.x=element_text(size=18),
          axis.title.y=element_blank(),
          legend.position = c(0.07,0.77),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          strip.placement = "outside",
          strip.background = element_blank(),
          strip.text.x = element_text(size = 16)) +
    guides(fill = FALSE)
  
  png(paste0("Figures/Transition Estimates/Stationary estimates/stationary_transitions_SPEED_", dirLabels[i], "wind.png"), width = 9, height = 7, units = "in", res = 600)
  print(myplot)
  dev.off()
  
}


# PLOT STATIONARY PROBABILITIES FOR EACH WIND DIRECTION -------------------





##### 3. PLOTTING STATIONARY PROBABILITIES BY PERSONALITY ###############

#  for average wind speed and direction

cov.m.bold <- data.frame(lod = "L", dir = 75, ws = 9, pers = "")
cov.m.mid <- data.frame(lod = "L", dir = 75, ws = 9, sex = "F")
cov.m.shy <- data.frame(lod = "L", dir = 75, ws = 9, sex = "M")
cov.m.shy <- data.frame(lod = "L", dir = 75, ws = 9, sex = "F")
cov.m.shy <- data.frame(lod = "L", dir = 75, ws = 9, sex = "F")
cov.m.shy <- data.frame(lod = "L", dir = 75, ws = 9, sex = "F")

s.cro.m <- as.data.frame(stationary(model = m.cro, covs = cov.cro.m))
s.cro.f <- as.data.frame(stationary(model = m.cro, covs = cov.cro.f))
s.sg.m <- as.data.frame(stationary(model = m.sg, covs = cov.sg.m))
s.sg.f <- as.data.frame(stationary(model = m.sg, covs = cov.sg.f))

# creating confidence intervals around stationary probabilities

stat.cro.m <- ci.stationary(model = m.cro, cov = cov.cro.m, alpha = 0.95) 
stat.cro.f <- ci.stationary(model = m.cro, cov = cov.cro.f, alpha = 0.95) 
stat.sg.m <- ci.stationary(model = m.sg, cov = cov.sg.m, alpha = 0.95) 
stat.sg.f <- ci.stationary(model = m.sg, cov = cov.sg.f, alpha = 0.95) 


stat.df <- data.frame(dir = rep(cov.sg.m$dir, 4),
                      prob = c(s.cro.m$travel, s.cro.m$search, s.cro.m$rest,
                               s.cro.f$travel, s.cro.f$search, s.cro.f$rest,
                               s.sg.m$travel, s.sg.m$search, s.sg.m$rest,
                               s.sg.f$travel, s.sg.f$search, s.sg.f$rest),
                      lwr = c(stat.cro.m[[1]]$travel, stat.cro.m[[1]]$search, stat.cro.m[[1]]$rest,
                              stat.cro.f[[1]]$travel, stat.cro.f[[1]]$search, stat.cro.f[[1]]$rest,
                              stat.sg.m[[1]]$travel, stat.sg.m[[1]]$search, stat.sg.m[[1]]$rest,
                              stat.sg.f[[1]]$travel, stat.sg.f[[1]]$search, stat.sg.f[[1]]$rest),
                      upr = c(stat.cro.m[[2]]$travel, stat.cro.m[[2]]$search, stat.cro.m[[2]]$rest,
                              stat.cro.f[[2]]$travel, stat.cro.f[[2]]$search, stat.cro.f[[2]]$rest,
                              stat.sg.m[[2]]$travel, stat.sg.m[[2]]$search, stat.sg.m[[2]]$rest,
                              stat.sg.f[[2]]$travel, stat.sg.f[[2]]$search, stat.sg.f[[2]]$rest),
                      state = rep(rep(c("travel", "search", "rest"), each = length(cov.sg.m$ws)), 4),
                      sex = rep(c("M", "F", "M", "F"), each = length(cov.sg.m$ws)*3),
                      site = rep(c("Crozet", "South Georgia"), each = length(cov.sg.m$ws)*6))

pd <- position_dodge()
ggplot(stat.df, aes(sex, prob)) + facet_wrap(~site)+
  geom_errorbar(aes(ymin=lwr, ymax=upr, group = sex), width=.25, position=pd) +
  geom_point(aes(fill = state), pch = 21, size = 2.7)+
  ylim(0, 1)+theme_bw() + #ylab("Stationary probability")+
  #xlab("Sex")+
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.text.x = element_text(size = 16),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_fill_manual(values = c("yellow", "red", "blue"))





