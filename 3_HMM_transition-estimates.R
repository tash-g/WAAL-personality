### CITATION: This code is modified from scripts shared in Clay et al. 2020, J. Anim. Ecol.

### AIM: Plot predicted transition probabilities for wind speeds, directions and boldness

# PREAMBLE --------------------------------------------------------------

library(ggplot2); library(momentuHMM); library(plyr)


### Create outputs and figures folders if they don't currently exist
out.path <- "./Data_outputs/"
if(dir.exists(out.path) == FALSE){
  dir.create(out.path)
}

figures.path <- "./Figures/"
if(dir.exists(figures.path) == FALSE){
  dir.create(figures.path)
}

## LOAD MODELS

# F best
file.in <- paste0("./Data_outputs/", paste0("F_mod_", 8, ".RData"))
load(file = file.in)
f.best.mod <- model

# M best
file.in <- paste0("./Data_outputs/", paste0("M_mod_", 6, ".RData"))
load(file = file.in)
m.best.mod <- model



# GET TRANSITION ESTIMATES BY SPEED FOR AVERAGE WIND DIRECTION -------------------------------------------

## Average wind direction is 75 degrees relative to bird direction; plots are for daylight only
## Transition predictions are made separately for each sex and for the extremes of boldness score ('shy' and 'bold')

### FEMALES

# Get personality * wind combinations
cov.ws <- seq(from = min(f.best.mod$data$WindSp), max(f.best.mod$data$WindSp), by = 0.2)
cov.pers <- c(min(f.best.mod$data$mean_BLUP_logit), max(f.best.mod$data$mean_BLUP_logit))
cov.comb <- expand.grid(cov.ws, cov.pers)

# Make the dataframe for predictions
cov.f_speed <- data.frame(WindSp=cov.comb$Var1, mean_BLUP_logit = cov.comb$Var2, LoD = "L", 
                        WindDir = 75)
cov.f_speed$LoD <- as.factor(as.character(cov.f_speed$LoD))
cov.f_speed$WindSp <- as.numeric(as.character(cov.f_speed$WindSp))
cov.f_speed$WindDir <- as.numeric(as.character(cov.f_speed$WindDir))
head(cov.f_speed)


# CIreal only allows single row of covariates, so this function iterates through range of covariate values 
# and outputs list of predicted transition probabilities and upper and lower CIs 
# (this takes around 20 minutes to run)

tock <- Sys.time()
ci.list_F <- lapply(1:nrow(cov.f_speed), function(x) {
  print(x)
  cov.sub.df <- cov.f_speed[x,]
  return(CIreal(f.best.mod,covs=cov.sub.df)$gamma)
})
tick <- Sys.time()
tick-tock

#save(ci.list_F, file = "Data_outputs/female_speedtransition_CIs.R")
load("Data_outputs/female_speedtransition_CIs.R")

# Extract means, upper and lower bounds and put into separate lists
means.f_speed <- lapply(ci.list_F, '[[', 1) 
lb.f_speed <- lapply(ci.list_F,'[[',3)
ub.f_speed <- lapply(ci.list_F,'[[',4)



### MALES

# Get all personality * wind combinations
cov.ws <- seq(from = min(m.best.mod$data$WindSp), max(m.best.mod$data$WindSp), by = 0.2)
cov.pers <- c(min(m.best.mod$data$mean_BLUP_logit), max(m.best.mod$data$mean_BLUP_logit))
cov.comb <- expand.grid(cov.ws, cov.pers)

# Construct dataframe for predictions
cov.m_speed <- data.frame(WindSp=cov.comb$Var1, mean_BLUP_logit = cov.comb$Var2, LoD = "L", 
                    WindDir = 75)
cov.m_speed$LoD <- as.factor(as.character(cov.m_speed$LoD))
cov.m_speed$WindSp <- as.numeric(as.character(cov.m_speed$WindSp))
cov.m_speed$WindDir <- as.numeric(as.character(cov.m_speed$WindDir))
head(cov.m_speed)


#  Get list of predicted transition probabilities and upper and lower CIs 
tock <- Sys.time()
ci.list_M <- lapply(1:nrow(cov.m_speed), function(x) {
  print(x)
  cov.sub.df <- cov.m_speed[x,]
  return(CIreal(m.best.mod,covs=cov.sub.df)$gamma)
})
tick <- Sys.time()
tick-tock


#save(ci.list_M, file = "Data_outputs/male_speedtransition_CIs.R")
load("Data_outputs/male_speedtransition_CIs.R")

# Extract means, upper and lower bounds and put into separate lists
means.m_speed <- lapply(ci.list_M, '[[', 1) 
lb.m_speed <- lapply(ci.list_M,'[[',3)
ub.m_speed <- lapply(ci.list_M,'[[',4)


# TABLE 1: EXTRACT TRANSITION ESTIMATE VALUES FOR MIN AND MAX WIND SPEED ----------------------------------------------------

## Rest - search transitions

# Get transition data
state1 <- 3
state2 <- 2
nb_states = 3

m.means <- unlist(lapply(means.m_speed,'[[',nb_states*(state2-1)+state1)) 
m.lb <- unlist(lapply(lb.m_speed,'[[',nb_states*(state2-1)+state1))
m.ub <- unlist(lapply(ub.m_speed,'[[',nb_states*(state2-1)+state1))
df.m <- data.frame(sex = "M", pers = cov.m_speed$mean_BLUP_logit, dir = cov.m_speed$WindDir,
                   wind=cov.m_speed$WindSp, mean = m.means, lower_bound = m.lb, upper_bound = m.ub)
df.m$persCat <- ifelse(df.m$pers == unique(cov.m_speed$mean_BLUP_logit)[1], "shy", "bold")


f.means <- unlist(lapply(means.f_speed,'[[',nb_states*(state2-1)+state1)) 
f.lb <- unlist(lapply(lb.f_speed,'[[',nb_states*(state2-1)+state1))
f.ub <- unlist(lapply(ub.f_speed,'[[',nb_states*(state2-1)+state1))
df.f <- data.frame(sex = "F", pers = cov.f_speed$mean_BLUP_logit, dir = cov.f_speed$WindDir,
                   wind=cov.f_speed$WindSp, mean = f.means, lower_bound = f.lb, upper_bound = f.ub)
df.f$persCat <- ifelse(df.f$pers == unique(cov.f_speed$mean_BLUP_logit)[1], "shy", "bold")

all.df <- rbind(df.m, df.f)


# Isolate probabilities
subset(all.df, sex == "F" & dir == 75 & wind ==  min(subset(all.df, sex == "F")$wind))[1,] # low = 0.1066338  0.08520636   0.1326683     shy
subset(all.df, sex == "F" & dir == 75 & wind ==  max(subset(all.df, sex == "F")$wind))[1,] # high = 0.2305976   0.1689537   0.3064392     shy

subset(all.df, sex == "F" & dir == 75 & wind ==  min(subset(all.df, sex == "F")$wind))[2,] # low = 0.165563    0.127179   0.2127083    bold
subset(all.df, sex == "F" & dir == 75 & wind ==  max(subset(all.df, sex == "F")$wind))[2,] # high = 0.1919329   0.1258836   0.2814778    bold


# Male values
subset(all.df, sex == "M" & dir == 75 & wind ==  min(subset(all.df, sex == "M")$wind))[1,] # low = 0.1473849   0.1182203   0.1822573     shy
subset(all.df, sex == "M" & dir == 75 & wind ==  max(subset(all.df, sex == "M")$wind))[1,] # high = 0.2140871   0.1247433   0.3423888     shy

subset(all.df, sex == "M" & dir == 75 & wind ==  min(subset(all.df, sex == "M")$wind))[2,] # high = 0.1114218  0.08683227   0.1418928   
subset(all.df, sex == "M" & dir == 75 & wind ==  max(subset(all.df, sex == "M")$wind))[2,] # low = 0.1652692  0.09349768   0.2753971    bold




## Search - travel transitions

# Get transition data
state1 <- 2
state2 <- 1
nb_states = 3

m.means <- unlist(lapply(means.m_speed,'[[',nb_states*(state2-1)+state1)) 
m.lb <- unlist(lapply(lb.m_speed,'[[',nb_states*(state2-1)+state1))
m.ub <- unlist(lapply(ub.m_speed,'[[',nb_states*(state2-1)+state1))
df.m <- data.frame(sex = "M", pers = cov.m_speed$mean_BLUP_logit, dir = cov.m_speed$WindDir,
                   wind=cov.m_speed$WindSp, mean = m.means, lower_bound = m.lb, upper_bound = m.ub)
df.m$persCat <- ifelse(df.m$pers == unique(cov.m_speed$mean_BLUP_logit)[1], "shy", "bold")


f.means <- unlist(lapply(means.f_speed,'[[',nb_states*(state2-1)+state1)) 
f.lb <- unlist(lapply(lb.f_speed,'[[',nb_states*(state2-1)+state1))
f.ub <- unlist(lapply(ub.f_speed,'[[',nb_states*(state2-1)+state1))
df.f <- data.frame(sex = "F", pers = cov.f_speed$mean_BLUP_logit, dir = cov.f_speed$WindDir,
                   wind=cov.f_speed$WindSp, mean = f.means, lower_bound = f.lb, upper_bound = f.ub)
df.f$persCat <- ifelse(df.f$pers == unique(cov.f_speed$mean_BLUP_logit)[1], "shy", "bold")

all.df <- rbind(df.m, df.f)

# Isolate probabilities
subset(all.df, sex == "F" & dir == 75 & wind ==  min(subset(all.df, sex == "F")$wind))[1,] # low = 0.09323681   0.0757844   0.1142117     shy
subset(all.df, sex == "F" & dir == 75 & wind ==  max(subset(all.df, sex == "F")$wind))[1,] # high = 0.1995975   0.1535666   0.2552649     shy

subset(all.df, sex == "F" & dir == 75 & wind ==  min(subset(all.df, sex == "F")$wind))[2,] # low =  0.1838433   0.1459631   0.2289191    bold
subset(all.df, sex == "F" & dir == 75 & wind ==  max(subset(all.df, sex == "F")$wind))[2,] # high = 0.1448276    0.101706   0.2021184    bold


# Male values
subset(all.df, sex == "M" & dir == 75 & wind ==  min(subset(all.df, sex == "M")$wind))[1,] # low =  0.05577518  0.04281957  0.07235435     shy
subset(all.df, sex == "M" & dir == 75 & wind ==  max(subset(all.df, sex == "M")$wind))[1,] # high = 0.3070312   0.2006025   0.4389219     shy

subset(all.df, sex == "M" & dir == 75 & wind ==  min(subset(all.df, sex == "M")$wind))[2,] # low = 0.1535877   0.1192771   0.1955761    bold
subset(all.df, sex == "M" & dir == 75 & wind ==  max(subset(all.df, sex == "M")$wind))[2,] # high = 0.09983382  0.05858481   0.1650345    bold


## Travel - search transitions

# Get transition data
state1 <- 1
state2 <- 2
nb_states = 3

m.means <- unlist(lapply(means.m_speed,'[[',nb_states*(state2-1)+state1)) 
m.lb <- unlist(lapply(lb.m_speed,'[[',nb_states*(state2-1)+state1))
m.ub <- unlist(lapply(ub.m_speed,'[[',nb_states*(state2-1)+state1))
df.m <- data.frame(sex = "M", pers = cov.m_speed$mean_BLUP_logit, dir = cov.m_speed$WindDir,
                   wind=cov.m_speed$WindSp, mean = m.means, lower_bound = m.lb, upper_bound = m.ub)
df.m$persCat <- ifelse(df.m$pers == unique(cov.m_speed$mean_BLUP_logit)[1], "shy", "bold")


f.means <- unlist(lapply(means.f_speed,'[[',nb_states*(state2-1)+state1)) 
f.lb <- unlist(lapply(lb.f_speed,'[[',nb_states*(state2-1)+state1))
f.ub <- unlist(lapply(ub.f_speed,'[[',nb_states*(state2-1)+state1))
df.f <- data.frame(sex = "F", pers = cov.f_speed$mean_BLUP_logit, dir = cov.f_speed$WindDir,
                   wind=cov.f_speed$WindSp, mean = f.means, lower_bound = f.lb, upper_bound = f.ub)
df.f$persCat <- ifelse(df.f$pers == unique(cov.f_speed$mean_BLUP_logit)[1], "shy", "bold")

all.df <- rbind(df.m, df.f)

# Isolate probabilities
subset(all.df, sex == "F" & dir == 0 & wind ==  min(subset(all.df, sex == "F")$wind))[1,] # low = 0.1388584   0.1159332   0.1654686     shy
subset(all.df, sex == "F" & dir == 0 & wind ==  max(subset(all.df, sex == "F")$wind))[1,] # high = 0.103329  0.08051051   0.1316886     shy

subset(all.df, sex == "F" & dir == 0 & wind ==  min(subset(all.df, sex == "F")$wind))[2,] # low = 0.1256671  0.09206856   0.1692412    bold
subset(all.df, sex == "F" & dir == 0 & wind ==  max(subset(all.df, sex == "F")$wind))[2,] # high = 0.1243336  0.09902084   0.1550039    bold

# Male values
subset(all.df, sex == "M" & dir == 0 & wind ==  min(subset(all.df, sex == "M")$wind))[1,] # low =  0.1354298   0.1102486    0.165294     shy
subset(all.df, sex == "M" & dir == 0 & wind ==  max(subset(all.df, sex == "M")$wind))[1,] # high = 0.1901425   0.1240988   0.2800944     shy

subset(all.df, sex == "M" & dir == 0 & wind ==  min(subset(all.df, sex == "M")$wind))[2,] # low = 0.04271847  0.02610382  0.06915713    bold
subset(all.df, sex == "M" & dir == 0 & wind ==  max(subset(all.df, sex == "M")$wind))[2,] # high =  0.1433595   0.1153959   0.1767454    bold




# FIGURE 2: PLOT TRANSITIONS ESTIMATES BY SPEED FOR AVERAGE WIND DIRECTION --------------------------------------------

shy_col <- "#00DD2F"
bold_col <- "purple"

behaviour <- data.frame(state = c("Directed travel", "Search", "Rest"), filename = c("travel", "search", "rest"))

nb_states = 3
labels <- c(F = "Female", M = "Male")

for (i in 1:3) {
  
  for (j in 1:3) {
    
    state1 <- i; print(i)
    state2 <- j; print(j) 
    
    ### MALE COEFFECIENTS
    
    m.means <- unlist(lapply(means.m_speed,'[[',nb_states*(state2-1)+state1)) 
    m.lb <- unlist(lapply(lb.m_speed,'[[',nb_states*(state2-1)+state1))
    m.ub <- unlist(lapply(ub.m_speed,'[[',nb_states*(state2-1)+state1))
    df.m <- data.frame(sex = "M", pers = cov.m_speed$mean_BLUP_logit, 
                       wind=cov.m_speed$WindSp, mean = m.means, lower_bound = m.lb, upper_bound = m.ub)
    df.m$persCat <- ifelse(df.m$pers == unique(cov.m_speed$mean_BLUP_logit)[1], "shy", "bold")
    
    ### FEMALE COEFFICIENTS
    
   f.means <- unlist(lapply(means.f_speed,'[[',nb_states*(state2-1)+state1)) 
   f.lb <- unlist(lapply(lb.f_speed,'[[',nb_states*(state2-1)+state1))
   f.ub <- unlist(lapply(ub.f_speed,'[[',nb_states*(state2-1)+state1))
   df.f <- data.frame(sex = "F", pers = cov.f_speed$mean_BLUP_logit, 
                       wind=cov.f_speed$WindSp, mean = f.means, lower_bound = f.lb, upper_bound = f.ub)
   df.f$persCat <- ifelse(df.f$pers == unique(cov.f_speed$mean_BLUP_logit)[1], "shy", "bold")
    
    ### COMBINE DATAFRAMES
    
    all.df <- rbind(df.m, df.f)
    all.df$sex <- as.factor(as.character(all.df$sex))
    all.df$persCat <- as.factor(as.character(all.df$persCat))
    
    
    ### BUILD PLOT
    
    dirPlot <- ggplot(all.df, aes(x = wind, y=mean)) + facet_wrap(~sex, labeller=labeller(sex=labels)) +
      geom_ribbon(size = 1, linetype = "blank", aes(ymin=lower_bound, ymax=upper_bound, col = persCat), alpha=0.15) + 
      geom_line(size = 1, linetype = "dashed", aes(col = persCat)) + 
      theme_bw() + ylab("Transition probability") +
      scale_x_continuous(limits=c(0, 23)) +
      #scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
      xlab("Wind speed (ms-1)") +
      theme(axis.text.x=element_text(size=15), 
            axis.text.y=element_text(size=15), 
            axis.title.x=element_text(size = 16),
            axis.title.y=element_text(size = 16),
            strip.text.x = element_text(size = 16),
            strip.placement = "outside",
            strip.background = element_blank(),
            plot.title = element_text(hjust = 0.5),
            legend.position = "none") +
      scale_color_manual(name = "", labels = c("Bold", "Shy"), values = c(bold_col, shy_col)) +
      ggtitle(paste0(behaviour$state[i], " - > ", behaviour$state[j]))
    
    # Scale y as 0.5 - 1 for same behaviour -> same behaviour; scale 0 - 0.5 for same -> different 
    if (i == j) { dirPlot <- dirPlot + scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0.5, 1)) } else {
      dirPlot <- dirPlot + scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0, 0.5))  }
      
    assign(paste0("speedPlot_", i, "_", j), dirPlot)
    
    png(paste0("Figures/Transition Estimates/SPEED_", behaviour$filename[i], "_", behaviour$filename[j],".png"), 
       width = 9, height = 7, units = "in", res = 350)
    print(dirPlot)
    dev.off()
    
  }
  
}



# GET TRANSITION ESTIMATES BY SPEED FOR CROSS/HEAD/TAIL WIND -------------------------

## Transitions calculated for each of crosswind (90 degrees to trajector), headwind (0 degrees),
## tailwind (180 degrees)
## Transition predictions are made separately for each sex and for the extremes of boldness score ('shy' and 'bold')

### FEMALES 

# Get all personality * wind combinations
cov.ws <- seq(from = min(f.best.mod$data$WindSp), max(f.best.mod$data$WindSp), by = 0.2)
cov.pers <- c(min(f.best.mod$data$mean_BLUP_logit), max(f.best.mod$data$mean_BLUP_logit))
cov.dir <- c(0, 90, 180)
cov.comb <- expand.grid(cov.ws, cov.pers, cov.dir)

# Construct dataset for predictions
cov.f_speed <- data.frame(WindSp=cov.comb$Var1, mean_BLUP_logit = cov.comb$Var2, LoD = "L", 
                          WindDir = cov.comb$Var3)
cov.f_speed$LoD <- as.factor(as.character(cov.f_speed$LoD))
cov.f_speed$WindSp <- as.numeric(as.character(cov.f_speed$WindSp))
cov.f_speed$WindDir <- as.numeric(as.character(cov.f_speed$WindDir))
head(cov.f_speed)


# Output list of predicted transition probabilities and upper and lower CIs 
tock <- Sys.time()
ci.list_F_bydir <- lapply(1:nrow(cov.f_speed), function(x) {
  print(x)
  cov.sub.df <- cov.f_speed[x,]
  return(CIreal(f.best.mod,covs=cov.sub.df)$gamma)
})
tick <- Sys.time()
tick-tock

#save(ci.list_F_bydir, file = "Data_outputs/female_speedtransitionbydir_CIs.R")
load("Data_outputs/female_speedtransitionbydir_CIs.R")

# Extract means, upper and lower bounds and put into separate lists
means.f_speed <- lapply(ci.list_F_bydir, '[[', 1) 
lb.f_speed <- lapply(ci.list_F_bydir,'[[',3)
ub.f_speed <- lapply(ci.list_F_bydir,'[[',4)


# MALES #

# Get all personality * wind possibilities
cov.ws <- seq(from = min(m.best.mod$data$WindSp), max(m.best.mod$data$WindSp), by = 0.2)
cov.pers <- c(min(m.best.mod$data$mean_BLUP_logit), max(m.best.mod$data$mean_BLUP_logit))
cov.dir <- c(0, 90, 180)
cov.comb <- expand.grid(cov.ws, cov.pers, cov.dir)

# Construct dataset for predictions
cov.m_speed <- data.frame(WindSp=cov.comb$Var1, mean_BLUP_logit = cov.comb$Var2, LoD = "L", 
                          WindDir = cov.comb$Var3)
cov.m_speed$LoD <- as.factor(as.character(cov.m_speed$LoD))
cov.m_speed$WindSp <- as.numeric(as.character(cov.m_speed$WindSp))
cov.m_speed$WindDir <- as.numeric(as.character(cov.m_speed$WindDir))
head(cov.m_speed)


# Output list of predicted transition probabilities and upper and lower CIs 
tock <- Sys.time()
ci.list_M_bydir <- lapply(1:nrow(cov.m_speed), function(x) {
  print(x)
  cov.sub.df <- cov.m_speed[x,]
  return(CIreal(m.best.mod,covs=cov.sub.df)$gamma)
})
tick <- Sys.time()
tick-tock


#save(ci.list_M_bydir, file = "Data_outputs/male_speedtransitionbydir_CIs.R")
load("Data_outputs/male_speedtransitionbydir_CIs.R")

# extract means, upper and lower bounds and put into separate lists
means.m_speed <- lapply(ci.list_M_bydir, '[[', 1) 
lb.m_speed <- lapply(ci.list_M_bydir,'[[',3)
ub.m_speed <- lapply(ci.list_M_bydir,'[[',4)






# PLOT TRANSITIONS ESTIMATES BY SPEED FOR CROSS/HEAD/TAILWIND (FIGURES S1-S3)  --------------------------------------------

shy_col <- "#00DD2F"
bold_col <- "purple"
mid_col <- "#a099a6"

behaviour <- data.frame(state = c("Directed travel", "Search", "Rest"), filename = c("travel", "search", "rest"))

nb_states = 3
labels <- c(F = "Female", M = "Male")

for(k in 1:3){
    
    for (i in 1:3) {
      
      for (j in 1:3) {
        
        wind <- k; print(k)
        
        state1 <- i; print(i)
        state2 <- j; print(j) 
        
        ### MALE COEFFECIENTS
        m.means <- unlist(lapply(means.m_speed,'[[',nb_states*(state2-1)+state1)) 
        m.lb <- unlist(lapply(lb.m_speed,'[[',nb_states*(state2-1)+state1))
        m.ub <- unlist(lapply(ub.m_speed,'[[',nb_states*(state2-1)+state1))
        df.m <- data.frame(sex = "M", pers = cov.m_speed$mean_BLUP_logit, dir = cov.m_speed$WindDir,
                           wind=cov.m_speed$WindSp, mean = m.means, lower_bound = m.lb, upper_bound = m.ub)
        df.m$persCat <- ifelse(df.m$pers == unique(cov.m_speed$mean_BLUP_logit)[1], "shy", "bold")
        df.m$windCat <- ifelse(df.m$dir == 0, "tail", "cross")
        df.m$windCat <- ifelse(df.m$dir == 180, "head", as.character(df.m$windCat))
        
        ### FEMALE COEFFICIENTS
        
        f.means <- unlist(lapply(means.f_speed,'[[',nb_states*(state2-1)+state1)) 
        f.lb <- unlist(lapply(lb.f_speed,'[[',nb_states*(state2-1)+state1))
        f.ub <- unlist(lapply(ub.f_speed,'[[',nb_states*(state2-1)+state1))
        df.f <- data.frame(sex = "F", pers = cov.f_speed$mean_BLUP_logit, dir = cov.f_speed$WindDir,
                           wind=cov.f_speed$WindSp, mean = f.means, lower_bound = f.lb, upper_bound = f.ub)
        df.f$persCat <- ifelse(df.f$pers == unique(cov.f_speed$mean_BLUP_logit)[1], "shy", "bold")
        df.f$windCat <- ifelse(df.f$dir == 0, "tail", "cross")
        df.f$windCat <- ifelse(df.f$dir == 180, "head", as.character(df.f$windCat))
        
        ### COMBINE DATAFRAMES
        
        all.df <- rbind(df.m, df.f)
        all.df$sex <- as.factor(as.character(all.df$sex))
        all.df$persCat <- as.factor(as.character(all.df$persCat))
        all.df <- subset(all.df, windCat == unique(all.df$windCat)[k])
        
        ### BUILD PLOT
        
        dirPlot <- ggplot(all.df, aes(x = wind, y=mean)) + facet_wrap(~sex, labeller=labeller(sex=labels)) +
          geom_ribbon(size = 1, linetype = "blank", aes(ymin=lower_bound, ymax=upper_bound, col = persCat), alpha=0.15) + 
          geom_line(size = 1, aes(col = persCat, linetype = persCat)) + 
          theme_bw() + ylab("Transition probability") +
          scale_x_continuous(limits=c(0, 23)) +
          xlab("Wind speed (ms-1)") +
          theme(axis.text.x=element_text(size=15), 
                axis.title.x=element_text(size = 16),
                axis.text.y=element_text(size=15), 
                axis.title.y=element_text(size = 16),
                strip.text.x = element_text(size = 16),
                strip.placement = "outside",
                strip.background = element_blank(),
                plot.title = element_text(hjust = 0.5),
                legend.position = "none") +
          scale_linetype_manual(name = "Personality", values=c("solid", "dashed"), labels = c("Bold", "Shy")) +
          scale_color_manual(name = "Personality", labels = c("Bold", "Shy"), values = c(bold_col, shy_col)) +
          ggtitle(paste0(behaviour$state[i], " - > ", behaviour$state[j]))
        
        # Scale y as 0.5 - 1 for same behaviour -> same behaviour; scale 0 - 0.5 for same -> different 
        if (i == j) { dirPlot <- dirPlot + scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0.5, 1)) } else {
          dirPlot <- dirPlot + scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0, 0.5))  }
        
        assign(paste0("speedPlot_", i, "_", j), dirPlot)
       
        png(paste0("Figures/Transition Estimates/By_wind/SPEED_", behaviour$filename[i], " to ", behaviour$filename[j], " in ", unique(all.df$windCat), "wind.png"), 
            width = 9, height = 7, units = "in", res = 350)
        print(dirPlot)
        dev.off()
        
      }
      
    }

}




