library(ggplot2); library(momentuHMM); library(plyr)


# PLOT PREDICTED TRANSITION PROBABILITIES FOR WIND SPEEDS AND PERSONALITY

# LOADING DATA --------------------------------------------------------------

## LOAD MODELS

# F best
file.in <- paste0("./Data_outputs/", paste0("F_mod_", 8, ".RData"))
load(file = file.in)
f.best.mod <- model

# M best
file.in <- paste0("./Data_outputs/", paste0("M_mod_", 6, ".RData"))
load(file = file.in)
m.best.mod <- model


# GET TRANSITION ESTIMATES - SPEED -------------------------------------------

## Average wind direction 75 degrees relative to bird direction; plot for daylight (TOMMY's METHODS)
## Predict transitions for a range of wind speeds separately for each sex and range of personality tpes

## Get response for 3 levels of personality - very shy, mid, very bold

# Get all personality * wind possibilities
cov.ws <- seq(from = min(f.best.mod$data$WindSp), max(f.best.mod$data$WindSp), by = 0.2)
cov.pers <- c(min(f.best.mod$data$mean_BLUP_logit), max(f.best.mod$data$mean_BLUP_logit))
cov.comb <- expand.grid(cov.ws, cov.pers)

# FEMALES 
cov.f_speed <- data.frame(WindSp=cov.comb$Var1, mean_BLUP_logit = cov.comb$Var2, LoD = "L", 
                          WindDir = 75)
cov.f_speed$LoD <- as.factor(as.character(cov.f_speed$LoD))
cov.f_speed$WindSp <- as.numeric(as.character(cov.f_speed$WindSp))
cov.f_speed$WindDir <- as.numeric(as.character(cov.f_speed$WindDir))
head(cov.f_speed)


# function CIreal only allows single row of covariates, so function iterates through range of covariate values and outputs list of predicted
# transition probabilities and upper and lower CIs 
tock <- Sys.time()
ci.list_F <- lapply(1:nrow(cov.f_speed), function(x) {
  print(x)
  cov.sub.df <- cov.f_speed[x,]
  return(CIreal(f.best.mod,covs=cov.sub.df)$gamma)
})
tick <- Sys.time()
tick-tock

save(ci.list_F, file = "Data_outputs/female_speedtransition_CIs.R")
#load("Data_outputs/female_speedtransition_CIs.R")

# extract means, upper and lower bounds and put into separate lists
means.f_speed <- lapply(ci.list_F, '[[', 1) 
lb.f_speed <- lapply(ci.list_F,'[[',3)
ub.f_speed <- lapply(ci.list_F,'[[',4)



# MALES #

# Get all personality * wind possibilities
cov.ws <- seq(from = min(m.best.mod$data$WindSp), max(m.best.mod$data$WindSp), by = 0.2)
cov.pers <- c(min(m.best.mod$data$mean_BLUP_logit), max(m.best.mod$data$mean_BLUP_logit))
cov.comb <- expand.grid(cov.ws, cov.pers)

cov.m_speed <- data.frame(WindSp=cov.comb$Var1, mean_BLUP_logit = cov.comb$Var2, LoD = "L", 
                          WindDir = 75)
cov.m_speed$LoD <- as.factor(as.character(cov.m_speed$LoD))
cov.m_speed$WindSp <- as.numeric(as.character(cov.m_speed$WindSp))
cov.m_speed$WindDir <- as.numeric(as.character(cov.m_speed$WindDir))
head(cov.m_speed)


# function CIreal only allows single row of covariates, so function iterates through range of covariate values and outputs list of predicted
# transition probabilities and upper and lower CIs 
tock <- Sys.time()
ci.list_M <- lapply(1:nrow(cov.m_speed), function(x) {
  print(x)
  cov.sub.df <- cov.m_speed[x,]
  return(CIreal(m.best.mod,covs=cov.sub.df)$gamma)
})
tick <- Sys.time()
tick-tock


save(ci.list_M, file = "Data_outputs/male_speedtransition_CIs.R")
#load("Data_outputs/male_speedtransition_CIs.R")

# extract means, upper and lower bounds and put into separate lists
means.m_speed <- lapply(ci.list_M, '[[', 1) 
lb.m_speed <- lapply(ci.list_M,'[[',3)
ub.m_speed <- lapply(ci.list_M,'[[',4)



# STATS FOR MANUSCRIPT ----------------------------------------------------

## Probability of landing for bold vs shy females

# Get transition data
state1 <- 2
state2 <- 3
nb_states = 3

f.means <- unlist(lapply(means.f_speed,'[[',nb_states*(state2-1)+state1)) 
f.lb <- unlist(lapply(lb.f_speed,'[[',nb_states*(state2-1)+state1))
f.ub <- unlist(lapply(ub.f_speed,'[[',nb_states*(state2-1)+state1))
df.f <- data.frame(sex = "F", pers = cov.f_speed$mean_BLUP_logit, 
                   wind=cov.f_speed$WindSp, mean = f.means, lower_bound = f.lb, upper_bound = f.ub)
df.f$persCat <- ifelse(df.f$pers == unique(cov.f_speed$mean_BLUP_logit)[1], "shy", "bold")

# Isolate probabilities
subset(df.f, wind == 0.04510045)$mean[2] # bold; 0.1523472
subset(df.f, wind == 0.04510045)$mean[1] # shy; 0.06222135



## Search - search transitions

# Get transition data
state1 <- 2
state2 <- 2
nb_states = 3

m.means <- unlist(lapply(means.m_speed,'[[',nb_states*(state2-1)+state1)) 
m.lb <- unlist(lapply(lb.m_speed,'[[',nb_states*(state2-1)+state1))
m.ub <- unlist(lapply(ub.m_speed,'[[',nb_states*(state2-1)+state1))
df.m <- data.frame(sex = "M", pers = cov.m_speed$mean_BLUP_logit, 
                   wind=cov.m_speed$WindSp, mean = m.means, lower_bound = m.lb, upper_bound = m.ub)
df.m$persCat <- ifelse(df.m$pers == unique(cov.m_speed$mean_BLUP_logit)[1], "shy", "bold")


f.means <- unlist(lapply(means.f_speed,'[[',nb_states*(state2-1)+state1)) 
f.lb <- unlist(lapply(lb.f_speed,'[[',nb_states*(state2-1)+state1))
f.ub <- unlist(lapply(ub.f_speed,'[[',nb_states*(state2-1)+state1))
df.f <- data.frame(sex = "F", pers = cov.f_speed$mean_BLUP_logit, 
                   wind=cov.f_speed$WindSp, mean = f.means, lower_bound = f.lb, upper_bound = f.ub)
df.f$persCat <- ifelse(df.f$pers == unique(cov.f_speed$mean_BLUP_logit)[1], "shy", "bold")

all.df <- rbind(df.m, df.f)

# Isolate probabilities
subset(all.df, sex == "F" & wind == 0.04510045)$mean[2,] # bold = 0.6807008   0.6284974   0.7287347
subset(all.df, sex == "F" & wind == 0.04510045)$mean[1,] # shy = 0.8521937 0.8273743   0.8739880 

subset(all.df, sex == "F" & wind > 20.6)[2,] # bold = 0.8360646   0.7806669   0.8796278
subset(all.df, sex == "F" & wind > 20.6)[1,] # shy = 0.7467889   0.6916119   0.7950195

# Male values
subset(all.df, sex == "M" & wind == 0.05049899)[2,] # bold = 0.777806   0.7346592   0.8156968
subset(all.df, sex == "M" & wind == 0.05049899)[1,] # shy = 0.8376979   0.8046178   0.8661092

subset(all.df, sex == "M" & wind > 24)[2,] # bold = 0.8215681   0.7636542   0.8677485
subset(all.df, sex == "M" & wind > 24)[1,] # shy = 0.7344602    0.654671   0.8014064



## Search - travel transitions

# Get transition data
state1 <- 2
state2 <- 1
nb_states = 3

m.means <- unlist(lapply(means.m_speed,'[[',nb_states*(state2-1)+state1)) 
m.lb <- unlist(lapply(lb.m_speed,'[[',nb_states*(state2-1)+state1))
m.ub <- unlist(lapply(ub.m_speed,'[[',nb_states*(state2-1)+state1))
df.m <- data.frame(sex = "M", pers = cov.m_speed$mean_BLUP_logit, 
                   wind=cov.m_speed$WindSp, mean = m.means, lower_bound = m.lb, upper_bound = m.ub)
df.m$persCat <- ifelse(df.m$pers == unique(cov.m_speed$mean_BLUP_logit)[1], "shy", "bold")


f.means <- unlist(lapply(means.f_speed,'[[',nb_states*(state2-1)+state1)) 
f.lb <- unlist(lapply(lb.f_speed,'[[',nb_states*(state2-1)+state1))
f.ub <- unlist(lapply(ub.f_speed,'[[',nb_states*(state2-1)+state1))
df.f <- data.frame(sex = "F", pers = cov.f_speed$mean_BLUP_logit, 
                   wind=cov.f_speed$WindSp, mean = f.means, lower_bound = f.lb, upper_bound = f.ub)
df.f$persCat <- ifelse(df.f$pers == unique(cov.f_speed$mean_BLUP_logit)[1], "shy", "bold")

all.df <- rbind(df.m, df.f)

# Isolate probabilities
subset(all.df, sex == "F" & wind == 0.04510045)[2,] # bold = 0.166952   0.1322787   0.2085297
subset(all.df, sex == "F" & wind == 0.04510045)[1,] # shy = 0.08558499  0.06965287    0.104751

subset(all.df, sex == "F" & wind > 20.6)[2,] # bold = 0.1343347  0.09417013   0.1880731
subset(all.df, sex == "F" & wind > 20.6)[1,] # shy = 0.1845176   0.1416828   0.2367307

# Male values
subset(all.df, sex == "M" & wind == 0.05049899)[2,] # bold = 0.1426522   0.1107547   0.1818573
subset(all.df, sex == "M" & wind == 0.05049899)[1,] # shy = 0.051035  0.03928624  0.06605572

subset(all.df, sex == "M" & wind > 24)[2,] # bold = 0.1163733   0.0791493   0.1679117
subset(all.df, sex == "M" & wind > 24)[1,] # shy = 0.2389782   0.1731958   0.3200735



# PLOT SPEED TRANSTIONS - LOOP --------------------------------------------

shy_col <- "#00DD2F"
bold_col <- "purple"
mid_col <- "#a099a6"

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
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            strip.placement = "outside",
            strip.background = element_blank(),
            plot.title = element_text(hjust = 0.5),
            legend.position = c(0.06,0.9)) +
      scale_color_manual(name = "Personality", labels = c("Bold", "Shy"), values = c(bold_col, shy_col)) +
      ggtitle(paste0(behaviour$state[i], " - > ", behaviour$state[j]))
    
    if (i == j) { dirPlot <- dirPlot + scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0.5, 1)) } else {
      dirPlot <- dirPlot + scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0, 0.5))  }
    
    assign(paste0("speedPlot_", i, "_", j), dirPlot)
    
    png(paste0("Figures/Transition Estimates/SPEED_", behaviour$filename[i], "_", behaviour$filename[j],".png"), 
        width = 9, height = 7, units = "in", res = 350)
    print(dirPlot)
    dev.off()
    
  }
  
}




# CALCULATE SPEED TRANSITIONS FOR CROSS/HEAD/TAIL WIND -------------------------

# Get all personality * wind possibilities
cov.ws <- seq(from = min(f.best.mod$data$WindSp), max(f.best.mod$data$WindSp), by = 0.2)
cov.pers <- c(min(f.best.mod$data$mean_BLUP_logit), max(f.best.mod$data$mean_BLUP_logit))
cov.dir <- c(0, 90, 180)
cov.comb <- expand.grid(cov.ws, cov.pers, cov.dir)

# FEMALES 
cov.f_speed <- data.frame(WindSp=cov.comb$Var1, mean_BLUP_logit = cov.comb$Var2, LoD = "L", 
                          WindDir = cov.comb$Var3)
cov.f_speed$LoD <- as.factor(as.character(cov.f_speed$LoD))
cov.f_speed$WindSp <- as.numeric(as.character(cov.f_speed$WindSp))
cov.f_speed$WindDir <- as.numeric(as.character(cov.f_speed$WindDir))
head(cov.f_speed)


# function CIreal only allows single row of covariates, so function iterates through range of covariate values and outputs list of predicted
# transition probabilities and upper and lower CIs 
tock <- Sys.time()
ci.list_F_bydir <- lapply(1:nrow(cov.f_speed), function(x) {
  print(x)
  cov.sub.df <- cov.f_speed[x,]
  return(CIreal(f.best.mod,covs=cov.sub.df)$gamma)
})
tick <- Sys.time()
tick-tock

save(ci.list_F_bydir, file = "Data_outputs/female_speedtransitionbydir_CIs.R")
#load("Data_outputs/female_speedtransitionbydir_CIs.R")

# extract means, upper and lower bounds and put into separate lists
means.f_speed <- lapply(ci.list_F_bydir, '[[', 1) 
lb.f_speed <- lapply(ci.list_F_bydir,'[[',3)
ub.f_speed <- lapply(ci.list_F_bydir,'[[',4)


# MALES #

# Get all personality * wind possibilities
cov.ws <- seq(from = min(m.best.mod$data$WindSp), max(m.best.mod$data$WindSp), by = 0.2)
cov.pers <- c(min(m.best.mod$data$mean_BLUP_logit), max(m.best.mod$data$mean_BLUP_logit))
cov.dir <- c(0, 90, 180)
cov.comb <- expand.grid(cov.ws, cov.pers, cov.dir)

cov.m_speed <- data.frame(WindSp=cov.comb$Var1, mean_BLUP_logit = cov.comb$Var2, LoD = "L", 
                          WindDir = cov.comb$Var3)
cov.m_speed$LoD <- as.factor(as.character(cov.m_speed$LoD))
cov.m_speed$WindSp <- as.numeric(as.character(cov.m_speed$WindSp))
cov.m_speed$WindDir <- as.numeric(as.character(cov.m_speed$WindDir))
head(cov.m_speed)


# function CIreal only allows single row of covariates, so function iterates through range of covariate values and outputs list of predicted
# transition probabilities and upper and lower CIs 
tock <- Sys.time()
ci.list_M_bydir <- lapply(1:nrow(cov.m_speed), function(x) {
  print(x)
  cov.sub.df <- cov.m_speed[x,]
  return(CIreal(m.best.mod,covs=cov.sub.df)$gamma)
})
tick <- Sys.time()
tick-tock


save(ci.list_M_bydir, file = "Data_outputs/male_speedtransitionbydir_CIs.R")
#load("Data_outputs/male_speedtransitionbydir_CIs.R")

# extract means, upper and lower bounds and put into separate lists
means.m_speed <- lapply(ci.list_M_bydir, '[[', 1) 
lb.m_speed <- lapply(ci.list_M_bydir,'[[',3)
ub.m_speed <- lapply(ci.list_M_bydir,'[[',4)


# PLOT SPEED TRANSTIONS FOR CROSS/HEAD/TAIL WIND  --------------------------------------------

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
              #axis.text.x=element_blank(),
              axis.title.x=element_text(size = 16),
              #axis.title.x=element_blank(),
              axis.text.y=element_text(size=15), 
              #axis.title.y = element_blank(),
              #axis.text.y = element_blank(),
              axis.title.y=element_text(size = 16),
              strip.text.x = element_text(size = 16),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              strip.placement = "outside",
              strip.background = element_blank(),
              plot.title = element_text(hjust = 0.5),
              #legend.position = c(0.06,0.9)) +
              legend.position = "none") +
        scale_linetype_manual(name = "Personality", values=c("solid", "dashed"), labels = c("Bold", "Shy")) +
        scale_color_manual(name = "Personality", labels = c("Bold", "Shy"), values = c(bold_col, shy_col)) +
        ggtitle(paste0(behaviour$state[i], " - > ", behaviour$state[j]))
      
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






# MULTI PLOT - SPEED ------------------------------------------------------

### Travel/travel + Travel/search

speedsplot_1_1b <- speedsplot_1_1 + theme(axis.text.x=element_blank(), 
                                          axis.title.x=element_blank(),
                                          axis.ticks.x=element_blank(),
                                          axis.title.y=element_text(size = 14, face = "plain"),
                                          title = element_text(size = 16, face = "italic"),
                                          legend.position = c(0.1,0.87))

speedsplot_1_2b <- speedsplot_1_2 + theme(title = element_text(size = 16, face = "italic"), 
                                          axis.title.y=element_text(size = 14, face = "plain"),
                                          axis.title.x=element_text(size = 14, face = "plain"),
                                          legend.position = "none")

png("Figures/transitions_speed_travel.png", width = 9, height = 13, units = "in", res = 350)
gridExtra::grid.arrange(speedsplot_1_1b, speedsplot_1_2b, nrow = 2)
dev.off()


### Rest/rest + Rest/search

speedsplot_3_3b <- speedsplot_3_3 + theme(axis.text.x=element_blank(), 
                                          axis.title.x=element_blank(),
                                          axis.ticks.x=element_blank(),
                                          axis.title.y=element_text(size = 14, face = "plain"),
                                          title = element_text(size = 16, face = "italic"),
                                          legend.position = c(0.1,0.15))

speedsplot_3_2b <- speedsplot_3_2 + theme(title = element_text(size = 16, face = "italic"), 
                                          axis.title.y=element_text(size = 14, face = "plain"),
                                          axis.title.x=element_text(size = 14, face = "plain"),
                                          legend.position = "none")

png("Figures/transitions_speed_rest.png", width = 9, height = 13, units = "in", res = 350)
gridExtra::grid.arrange(speedsplot_3_3b, speedsplot_3_2b, nrow = 2)
dev.off()


### Search/Search + Search/travel + Search/rest

speedsplot_2_2b <- speedsplot_2_2 + theme(axis.text.x=element_blank(), 
                                          axis.title.x=element_blank(),
                                          axis.ticks.x=element_blank(),
                                          axis.title.y=element_text(size = 14, face = "plain"),
                                          title = element_text(size = 16, face = "italic"),
                                          legend.position = "none")

speedsplot_2_1b <- speedsplot_2_1 + theme(axis.text.x=element_blank(), 
                                          axis.title.x=element_blank(),
                                          axis.ticks.x=element_blank(),
                                          title = element_text(size = 16, face = "italic"),
                                          axis.title.y=element_text(size = 14, face = "plain"),
                                          legend.position = c(0.1,0.87))

speedsplot_2_3b <- speedsplot_2_3 + theme(title = element_text(size = 16, face = "italic"), 
                                          axis.title.y=element_text(size = 14, face = "plain"),
                                          axis.title.x=element_text(size = 14, face = "plain"),
                                          legend.position = "none")

png("Figures/transitions_speed_search.png", width = 9, height = 19, units = "in", res = 350)
gridExtra::grid.arrange(speedsplot_2_1b, speedsplot_2_2b, speedsplot_2_3b, nrow = 3)
dev.off()





##### GET TRANSITION ESTIMATES - DIRECTION ###############

# predict wind direction for average wind speed which is 9 ms-1, and for daylight hours when birds are most active. 

# FEMALES
cov.wd <- seq(from = min(f.best.mod$data$WindDir), max(f.best.mod$data$WindDir), by = 1)
cov.pers <- c(min(f.best.mod$data$mean_BLUP_logit), max(f.best.mod$data$mean_BLUP_logit))
cov.comb <- expand.grid(cov.wd, cov.pers)

cov.f_dir <- data.frame(WindDir=cov.comb$Var1, mean_BLUP_logit = cov.comb$Var2, LoD = "L", 
                        WindSp = 9, sex = "F")
cov.f_dir$LoD <- as.factor(as.character(cov.f_dir$LoD))
cov.f_dir$WindSp <- as.numeric(as.character(cov.f_dir$WindSp))
cov.f_dir$WindDir <- as.numeric(as.character(cov.f_dir$WindDir))

# Get CIs

tock <- Sys.time()
ci.list_F_dir <- lapply(1:nrow(cov.f_dir), function(x) {
  print(x)
  cov.sub.df <- cov.f_dir[x,]
  return(CIreal(f.best.mod,covs=cov.sub.df)$gamma)
})
tick <- Sys.time()
tick-tock

save(ci.list_F_dir, file = "Data_outputs/female_dirtransition_CIs.R")
#load("Data_outputs/female_dirtransition_CIs.R")

# extract means, upper and lower bounds and put into separate lists
means.f_dir <- lapply(ci.list_F_dir, '[[', 1) 
lb.f_dir <- lapply(ci.list_F_dir,'[[',3)
ub.f_dir <- lapply(ci.list_F_dir,'[[',4)



# MALES #
cov.wd <- seq(from = min(m.best.mod$data$WindDir), max(m.best.mod$data$WindDir), by = 1)
cov.pers <- c(min(m.best.mod$data$mean_BLUP_logit), max(m.best.mod$data$mean_BLUP_logit))
cov.comb <- expand.grid(cov.wd, cov.pers)

cov.m_dir <- data.frame(WindDir=cov.comb$Var1, mean_BLUP_logit = cov.comb$Var2, LoD = "L", 
                        WindSp = 9)
cov.m_dir$LoD <- as.factor(as.character(cov.m_dir$LoD))
cov.m_dir$WindSp <- as.numeric(as.character(cov.m_dir$WindSp))
cov.m_dir$WindDir <- as.numeric(as.character(cov.m_dir$WindDir))


# Get CIs
tock <- Sys.time()
ci.list_M_dir <- lapply(1:nrow(cov.m_dir), function(x) {
  print(x)
  cov.sub.df <- cov.m_dir[x,]
  return(CIreal(m.best.mod,covs=cov.sub.df)$gamma)
})
tick <- Sys.time()
tick-tock


save(ci.list_M_dir, file = "Data_outputs/male_dirtransition_CIs.R")
#load("Data_outputs/male_dirtransition_CIs.R")

# extract means, upper and lower bounds and put into separate lists
means.m_dir <- lapply(ci.list_M_dir, '[[', 1) 
lb.m_dir <- lapply(ci.list_M_dir,'[[',3)
ub.m_dir <- lapply(ci.list_M_dir,'[[',4)

## LOOP THROUGH EACH TRANSITION

# PLOT DIRECTION TRANSITIONS - LOOP ----------------------------------------------

behaviour <- data.frame(state = c("Directed travel", "Search", "Rest"), filename = c("travel", "search", "rest"))
nb_states = 3
labels <- c(F = "Female", M = "Male")

for (i in 1:3) {
  
  for (j in 1:3) {
    
    state1 <- i; print(i)
    state2 <- j; print(j) 
    
    ### MALE COEFFECIENTS
    
    m.means <- unlist(lapply(means.m_dir,'[[',nb_states*(state2-1)+state1)) 
    m.lb <- unlist(lapply(lb.m_dir,'[[',nb_states*(state2-1)+state1))
    m.ub <- unlist(lapply(ub.m_dir,'[[',nb_states*(state2-1)+state1))
    df.m <- data.frame(sex = "M", pers = cov.m_dir$mean_BLUP_logit, 
                       dir = cov.m_dir$WindDir, mean = m.means, lower_bound = m.lb, upper_bound = m.ub)
    df.m$persCat <- ifelse(df.m$pers == unique(cov.m_dir$mean_BLUP_logit)[1], "shy", "bold")
    
    ### FEMALE COEFFICIENTS
    
    f.means <- unlist(lapply(means.f_dir,'[[',nb_states*(state2-1)+state1)) 
    f.lb <- unlist(lapply(lb.f_dir,'[[',nb_states*(state2-1)+state1))
    f.ub <- unlist(lapply(ub.f_dir,'[[',nb_states*(state2-1)+state1))
    df.f <- data.frame(sex = "F", pers = cov.f_dir$mean_BLUP_logit, 
                       dir = cov.f_dir$WindDir, mean = f.means, lower_bound = f.lb, upper_bound = f.ub)
    df.f$persCat <- ifelse(df.f$pers == unique(cov.f_dir$mean_BLUP_logit)[1], "shy", "bold")
    
    ### COMBINE DATAFRAMES
    
    all.df <- rbind(df.m, df.f)
    all.df$sex <- as.factor(as.character(all.df$sex))
    all.df$persCat <- as.factor(as.character(all.df$persCat))
    
    
    ### BUILD PLOT
    
    dirPlot <- ggplot(all.df, aes(x = dir, y=mean)) + facet_wrap(~sex, labeller=labeller(sex=labels)) +
      geom_ribbon(size = 1, linetype = "blank", aes(ymin=lower_bound, ymax=upper_bound, col = persCat), alpha=0.15) + 
      geom_line(size = 1, linetype = "dashed", aes(col = persCat)) + 
      theme_bw() + ylab("Transition probability") +
      scale_x_continuous(limits=c(0, 180), breaks = seq(0, 180, 60)) +
      xlab("Relative wind direction (°)") +
      #scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0,0.3)) +
      theme(axis.text.x=element_text(size=15), 
            axis.text.y=element_text(size=15), 
            axis.title.x=element_text(size = 16),
            axis.title.y=element_text(size = 16),
            strip.text.x = element_text(size = 16),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            strip.placement = "outside",
            strip.background = element_blank(),
            plot.title = element_text(hjust = 0.5),
            legend.position = c(0.06,0.9)) +
      scale_color_manual(name = "Personality", labels = c("Bold", "Shy"), values = c(bold_col, shy_col)) +
      ggtitle(paste0(behaviour$state[i], " - > ", behaviour$state[j]))
    
    if (i == j) { dirPlot <- dirPlot + scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0.5, 1)) } else {
      dirPlot <- dirPlot + scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0, 0.5))  }
    
    assign(paste0("dirPlot_", i, "_", j), dirPlot)
    
    png(paste0("Figures/Transition Estimates/DIR_", behaviour$filename[i], "_", behaviour$filename[j],".png"), 
        width = 9, height = 7, units = "in", res = 350)
    print(dirPlot)
    dev.off()
    
  }
  
}


# MULTI PLOT - DIRECTION ------------------------------------------------------

### Travel/travel + Travel/search

dirPlot_1_1b <- dirPlot_1_1 + theme(axis.text.x=element_blank(), 
                                    axis.title.x=element_blank(),
                                    axis.ticks.x=element_blank(),
                                    axis.title.y=element_text(size = 14, face = "plain"),
                                    title = element_text(size = 16, face = "italic"),
                                    legend.position = c(0.1,0.12))

dirPlot_1_2b <- dirPlot_1_2 + theme(title = element_text(size = 16, face = "italic"), 
                                    axis.title.y=element_text(size = 14, face = "plain"),
                                    axis.title.x=element_text(size = 14, face = "plain"),
                                    legend.position = "none")

png("Figures/transitions_direction_travel.png", width = 9, height = 13, units = "in", res = 350)
gridExtra::grid.arrange(dirPlot_1_1b, dirPlot_1_2b, nrow = 2)
dev.off()


### Rest/rest + Rest/search

dirPlot_3_3b <- dirPlot_3_3 + theme(axis.text.x=element_blank(), 
                                    axis.title.x=element_blank(),
                                    axis.ticks.x=element_blank(),
                                    axis.title.y=element_text(size = 14, face = "plain"),
                                    title = element_text(size = 16, face = "italic"),
                                    legend.position = c(0.1,0.12))

dirPlot_3_2b <- dirPlot_3_2 + theme(title = element_text(size = 16, face = "italic"), 
                                    axis.title.y=element_text(size = 14, face = "plain"),
                                    axis.title.x=element_text(size = 14, face = "plain"),
                                    legend.position = "none")

png("Figures/transitions_direction_rest.png", width = 9, height = 13, units = "in", res = 350)
gridExtra::grid.arrange(dirPlot_3_3b, dirPlot_3_2b, nrow = 2)
dev.off()


### Search/Search + Search/travel + Search/rest

dirPlot_2_2b <- dirPlot_2_2 + theme(axis.text.x=element_blank(), 
                                    axis.title.x=element_blank(),
                                    axis.ticks.x=element_blank(),
                                    axis.title.y=element_text(size = 14, face = "plain"),
                                    title = element_text(size = 16, face = "italic"),
                                    legend.position = "none")

dirPlot_2_1b <- dirPlot_2_1 + theme(axis.text.x=element_blank(), 
                                    axis.title.x=element_blank(),
                                    axis.ticks.x=element_blank(),
                                    title = element_text(size = 16, face = "italic"),
                                    axis.title.y=element_text(size = 14, face = "plain"),
                                    legend.position = c(0.1,0.12))

dirPlot_2_3b <- dirPlot_2_3 + theme(title = element_text(size = 16, face = "italic"), 
                                    axis.title.y=element_text(size = 14, face = "plain"),
                                    axis.title.x=element_text(size = 14, face = "plain"),
                                    legend.position = "none")

png("Figures/transitions_dir_search.png", width = 9, height = 19, units = "in", res = 350)
gridExtra::grid.arrange(dirPlot_2_1b, dirPlot_2_2b, dirPlot_2_3b, nrow = 3)
dev.off()

##########  3. CALCULATING PREDICTED TRANSITION PROBABILITIES FROM REST TO SEARCH (TAKE OFF) AT WIND SPEED INTERVALS - FOR TABLE 3 ######


#### 3A. CREATE DATAFRAMES OF PREDICTED VALUES  - same as step 1A. #####


# MALE #
cov.ws <- seq(from = min(m.best.mod$data$WindSp), max(m.best.mod$data$WindSp), by = 1)
cov.pers <- c(min(m.best.mod$data$mean_BLUP_logit), median(m.best.mod$data$mean_BLUP_logit), max(m.best.mod$data$mean_BLUP_logit))
cov.comb <- expand.grid(cov.ws, cov.pers)
cov.m <- data.frame(WindSp=cov.comb$Var1, mean_BLUP_logit = cov.comb$Var2, LoD = "L", 
                    sex = "M", dir = 75)

ci.list <- lapply(1:nrow(cov.m), function(x) {
  print(x)
  cov.sub.df <- cov.m[x,]
  return(CIreal(m.best.mod,covs=cov.sub.df)$gamma)
})

means.m <- lapply(ci.list_M, '[[', 1) 
lb.m <- lapply(ci.list_M,'[[',3)
ub.m <- lapply(ci.list_M,'[[',4)

# FEMALE #
cov.ws <- seq(from = min(f.best.mod$data$WindSp), max(f.best.mod$data$WindSp), by = 1)
cov.pers <- c(min(f.best.mod$data$mean_BLUP_logit), median(f.best.mod$data$mean_BLUP_logit), max(f.best.mod$data$mean_BLUP_logit))
cov.comb <- expand.grid(cov.ws, cov.pers)
cov.f <- data.frame(WindSp=cov.comb$Var1, mean_BLUP_logit = cov.comb$Var2, LoD = "L", 
                    sex = "F", dir = 75)

ci.list <- lapply(1:nrow(cov.f), function(x) {
  print(x)
  cov.sub.df <- cov.f[x,]
  return(CIreal(f.best.mod,covs=cov.sub.df)$gamma)
})

means.f <- lapply(ci.list_F, '[[', 1) 
lb.f <- lapply(ci.list_F,'[[',3)
ub.f <- lapply(ci.list_F,'[[',4)



#### 3B. EXTRACTING PROBABILITIES FOR TRANSITION FROM REST TO SEARCH #####

state1 = 3
state2 = 2
nb_states = 3

m.means <- unlist(lapply(means.m,'[[',nb_states*(state2-1)+state1)) 
m.lb <- unlist(lapply(lb.m,'[[',nb_states*(state2-1)+state1))
m.ub <- unlist(lapply(ub.m,'[[',nb_states*(state2-1)+state1))
df.m <- data.frame(sex = "M", pers = cov.m_speed$mean_BLUP_logit, 
                   ws = cov.m_speed$WindSp, dir = cov.m_speed$WindDir, mean = m.means, lower_bound = m.lb, upper_bound = m.ub)
df.m$persCat <- ifelse(df.m$pers == unique(cov.m$mean_BLUP_logit)[1], "shy", "bold")


### FEMALE COEFFICIENTS

f.means <- unlist(lapply(means.f,'[[',nb_states*(state2-1)+state1)) 
f.lb <- unlist(lapply(lb.f,'[[',nb_states*(state2-1)+state1))
f.ub <- unlist(lapply(ub.f,'[[',nb_states*(state2-1)+state1))
df.f <- data.frame(sex = "F", pers = cov.f_speed$mean_BLUP_logit, 
                   ws = cov.f_speed$WindSp, dir = cov.f_speed$WindDir, mean = f.means, lower_bound = f.lb, upper_bound = f.ub)
df.f$persCat <- ifelse(df.f$pers == unique(cov.f$mean_BLUP_logit)[1], "shy", "bold")

all.df <- rbind(df.m, df.f)

## Summarise wind speed categories
all.df$wind_cat <- 1
all.df$wind_cat[all.df$ws > 5 & all.df$ws <= 15] <- 2
all.df$wind_cat[all.df$ws > 15] <- 3

all.df$sex <- as.factor(as.character(all.df$sex))
all.df$persCat <- as.factor(as.character(all.df$persCat))
all.df$wind_cat <- as.factor(as.character(all.df$wind_cat))




