### CITATION: This code is modified from scripts shared in Clay et al. 2020, J. Anim. Ecol.

### AIM: Plot predicted transition probabilities for wind speeds, directions and boldness

# PREAMBLE --------------------------------------------------------------

# Define the packages
packages <- c("ggplot2", "momentuHMM", "plyr", "patchwork", "ggpubr", "gridExtra",
              "grid", "dplyr")

# Install packages not yet installed - change lib to library path
#installed_packages <- packages %in% rownames(installed.packages())

#if (any(installed_packages == FALSE)) {
#  install.packages(packages[!installed_packages], lib = "C:/Users/libraryPath")
#}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# Set colours
shy_col <- "#00DD2F"
bold_col <- "purple"
inter_col <- "#be254e"

## Create outputs and figures folders if they don't currently exist
out.path <- "./Data_outputs/"
if (dir.exists(out.path) == FALSE) {
  dir.create(out.path)
}

figures.path <- "./Figures/"
if (dir.exists(figures.path) == FALSE) {
  dir.create(figures.path)
}


figures.path <- "./Figures/Transition Estimates/"
if (dir.exists(figures.path) == FALSE) {
  dir.create(figures.path)
}


figures.path <- "./Figures/Transition Estimates/By_wind/"
if (dir.exists(figures.path) == FALSE) {
  dir.create(figures.path)
}


# Set boldness quantiles ----------------------------------------------

## Set personality quantiles
upperQ <- 0.75 # 90
lowerQ <- 0.25 # 10

# For labelling outputs - 25% or 10% quantiles
label <- ifelse(upperQ == 0.75, 25, 10)

##  Load data and find quantiles
boldness <- read.csv(file = "./Data_inputs/WAAL_GPS_2010-2021_personality_wind.csv", 
                     na.strings = "NA")

# Isolate individuals
boldness <- boldness %>% 
  dplyr::group_by(Ring) %>%
  dplyr::summarise(boldness = mean_BLUP_logit[1], Sex = Sex[1])

# Extract quantiles
quantiles <- boldness %>% 
  dplyr::group_by(Sex) %>% 
  dplyr::summarise(upperQuan = quantile(boldness, prob = upperQ),
            lowerQuan = quantile(boldness, prob = lowerQ))

# Clear space
rm(boldness); gc()

# Load models ----------------------------------------------------------

# F best
file.in <- paste0("./Data_outputs/", paste0("F_mod_", 4, ".RData"))
load(file = file.in)
f.best.mod <- mod.p

# M best
file.in <- paste0("./Data_outputs/", paste0("M_mod_", 9, ".RData"))
load(file = file.in)
m.best.mod <- mod.p


# FEMALES : GET TRANSITION ESTIMATES BY SPEED FOR CROSSWIND DIRECTION -------------------------------------------

## Get personality * wind combinations
cov.ws <- seq(from = min(f.best.mod$data$WindSp), max(f.best.mod$data$WindSp), by = 0.25)

# Subsample personality to speed up processing
cov.pers <- unique(f.best.mod$data$mean_BLUP_logit)
cov.pers <- cov.pers[order(cov.pers)]
cov.pers <- cov.pers[c(TRUE,FALSE)]

cov.comb <- expand.grid(cov.ws, cov.pers)

# Make the dataframe for predictions
cov.f_speed <- data.frame(WindSp = cov.comb$Var1, mean_BLUP_logit = cov.comb$Var2, LoD = "L", 
                        WindDir = 90)
cov.f_speed$LoD <- as.factor(as.character(cov.f_speed$LoD))
cov.f_speed$WindSp <- as.numeric(as.character(cov.f_speed$WindSp))
cov.f_speed$WindDir <- as.numeric(as.character(cov.f_speed$WindDir))


## CIreal only allows single row of covariates, so this function iterates through range of covariate values 
## and outputs list of predicted transition probabilities and upper and lower CIs 
## (this takes around 2.5 hours to run)

tock <- Sys.time()
ci.list_F <- lapply(1:nrow(cov.f_speed), function(x) {
  print(x)
  cov.sub.df <- cov.f_speed[x,]
  return(CIreal(f.best.mod,covs = cov.sub.df)$gamma)
})
tick <- Sys.time()
tick - tock

# save(ci.list_F, file = paste0("Data_outputs/female_speedtransition_CIs_quantile", label, ".R"))
 load(paste0("Data_outputs/female_speedtransition_CIs_quantile", label, ".R"))

# Extract means, upper and lower bounds and put into separate lists
means.f_speed <- lapply(ci.list_F, '[[', 1) 
lb.f_speed <- lapply(ci.list_F,'[[',3)
ub.f_speed <- lapply(ci.list_F,'[[',4)




# MALES : GET TRANSITION ESTIMATES FOR UPPER AND LOWER BOLDNESS QUANTILES ----------

# Get all personality * wind combinations
cov.pers <- unique(m.best.mod$data$mean_BLUP_logit)
cov.pers <- cov.pers[order(cov.pers)]
cov.pers <- cov.pers[c(TRUE,FALSE)]

cov.m <- data.frame(LoD = "L", WindDir =  90, WindSp = 9, mean_BLUP_logit = cov.pers)

#  Get list of predicted transition probabilities and upper and lower CIs 
tock <- Sys.time()
ci.list_M <- lapply(1:nrow(cov.m), function(x) {
  print(x)
  cov.sub.df <- cov.m[x, ]
  return(CIreal(m.best.mod, covs = cov.sub.df)$gamma)
})
tick <- Sys.time()
tick - tock

# save(ci.list_M, file = paste0("Data_outputs/male_stationaryTransitions_quantile", label, ".R"))
 load(paste0("Data_outputs/male_stationaryTransitions_quantile", label, ".R"))

# Extract means, upper and lower bounds and put into separate lists
means.m <- lapply(ci.list_M, '[[', 1) 
lb.m <- lapply(ci.list_M,'[[',3)
ub.m <- lapply(ci.list_M,'[[',4)


# TABLE 1: EXTRACT TRANSITION ESTIMATE VALUES FOR DIFFERENT WIND SPEEDS ----------------------------------------------------

# Rest (3) -> Search (2) ----------------------------------------------------------

state1 <- 3
state2 <- 2
nb_states = 3

# Get transition data - MALES
m.means <- unlist(lapply(means.m, '[[', nb_states * (state2 - 1) + state1))
m.lb <- unlist(lapply(lb.m, '[[', nb_states * (state2 - 1) + state1))
m.ub <- unlist(lapply(ub.m, '[[', nb_states * (state2 - 1) + state1))

df.m <- data.frame(
  sex = "M",
  pers = cov.m$mean_BLUP_logit,
  mean = m.means,
  lower_bound = m.lb,
  upper_bound = m.ub
)

df.m$persCat <- ifelse(df.m$pers < quantiles$lowerQuan[quantiles$Sex == "M"], "shy", "intermediate")
df.m$persCat <- ifelse(df.m$pers >  quantiles$upperQuan[quantiles$Sex == "M"], "bold", as.character(df.m$persCat))

trans_rest.search.M <- ddply(df.m,.(sex, persCat), summarize, p_mean = mean(mean), p_low = mean(lower_bound), p_upp = mean(upper_bound))
trans_rest.search.M <- subset(trans_rest.search.M, persCat != "intermediate")
trans_rest.search.M

#   sex persCat    p_mean     p_low     p_upp
#1   M    bold 0.1897556 0.1795467 0.2004096
#2   M     shy 0.2037270 0.1936286 0.2142163

# Get transition data - FEMALES
f.means <- unlist(lapply(means.f_speed, '[[', nb_states * (state2 - 1) + state1)) 
f.lb <- unlist(lapply(lb.f_speed, '[[', nb_states * (state2 - 1) + state1))
f.ub <- unlist(lapply(ub.f_speed, '[[', nb_states * (state2 - 1) + state1))

df.f <- data.frame(
  sex = "F",
  pers = cov.f_speed$mean_BLUP_logit,
  dir = cov.f_speed$WindDir,
  wind = cov.f_speed$WindSp,
  mean = f.means,
  lower_bound = f.lb,
  upper_bound = f.ub
)

df.f$persCat <- ifelse(df.f$pers < quantiles$lowerQuan[quantiles$Sex == "F"], "shy", "intermediate")
df.f$persCat <- ifelse(df.f$pers >  quantiles$upperQuan[quantiles$Sex == "F"], "bold", as.character(df.f$persCat))

# Summarise wind into different categories
df.f$windCat <- "low"
df.f$windCat[df.f$wind > 5 & df.f$wind <= 10] <- "mid"
df.f$windCat[df.f$wind > 10] <- "high"


# Summarise transition probabilities
trans_rest.search.F <- ddply(df.f,.(sex, persCat, windCat), summarize, p_mean = mean(mean), p_low = mean(lower_bound), p_upp = mean(upper_bound))
trans_rest.search.F <- subset(trans_rest.search.F, persCat != "intermediate")
trans_rest.search.F

#   sex persCat windCat    p_mean     p_low     p_upp
#1   F    bold    high 0.2457608 0.2216861 0.2716093
#2   F    bold     low 0.1121965 0.1036770 0.1213820
#3   F    bold     mid 0.1538708 0.1464917 0.1615626
#7   F     shy    high 0.2683220 0.2440922 0.2939642
#8   F     shy     low 0.1084081 0.1002293 0.1172137
#9   F     shy     mid 0.1566880 0.1495245 0.1641349

# Search (2) -> Travel (1) ------------------------------------------------

state1 <- 2
state2 <- 1
nb_states <- 3

# Get transition data - MALES
m.means <- unlist(lapply(means.m, '[[', nb_states * (state2 - 1) + state1)) 
m.lb <- unlist(lapply(lb.m, '[[', nb_states * (state2 - 1) + state1))
m.ub <- unlist(lapply(ub.m, '[[', nb_states * (state2 - 1) + state1))

df.m <- data.frame(
  sex = "M",
  pers = cov.m$mean_BLUP_logit,
  mean = m.means,
  lower_bound = m.lb,
  upper_bound = m.ub
)

df.m$persCat <- ifelse(df.m$pers < quantiles$lowerQuan[quantiles$Sex == "M"], "shy", "intermediate")
df.m$persCat <- ifelse(df.m$pers >  quantiles$upperQuan[quantiles$Sex == "M"], "bold", as.character(df.m$persCat))
trans_search.travel.M <- ddply(df.m,.(sex, persCat), summarize, p_mean = mean(mean), p_low = mean(lower_bound), p_upp = mean(upper_bound))
trans_search.travel.M <- subset(trans_search.travel.M, persCat != "intermediate")
trans_search.travel.M

#sex persCat    p_mean     p_low     p_upp
#1   M    bold 0.1422137 0.1352943 0.1494306
#2   M     shy 0.1155982 0.1102372 0.1211883

# Get transition data - FEMALES
f.means <- unlist(lapply(means.f_speed, '[[', nb_states * (state2 - 1) + state1)) 
f.lb <- unlist(lapply(lb.f_speed, '[[', nb_states * (state2 - 1) + state1))
f.ub <- unlist(lapply(ub.f_speed, '[[', nb_states * (state2 - 1) + state1))
df.f <- data.frame(
  sex = "F",
  pers = cov.f_speed$mean_BLUP_logit,
  dir = cov.f_speed$WindDir,
  wind = cov.f_speed$WindSp,
  mean = f.means,
  lower_bound = f.lb,
  upper_bound = f.ub
)

df.f$persCat <- ifelse(df.f$pers < quantiles$lowerQuan[quantiles$Sex == "F"], "shy", "intermediate")
df.f$persCat <- ifelse(df.f$pers >  quantiles$upperQuan[quantiles$Sex == "F"], "bold", as.character(df.f$persCat))

# Summarise wind into different categories
df.f$windCat <- "low"
df.f$windCat[df.f$wind > 5 & df.f$wind <= 10] <- "mid"
df.f$windCat[df.f$wind > 10] <- "high"


# Summarise transition probabilities
trans_search.travel.F <- ddply(df.f,.(sex, persCat, windCat), summarize, p_mean = mean(mean), p_low = mean(lower_bound), p_upp = mean(upper_bound))
trans_search.travel.F <- subset(trans_search.travel.F, persCat != "intermediate")
trans_search.travel.F

#   sex persCat windCat    p_mean      p_low      p_upp
#1   F    bold    high 0.1755042 0.1621072 0.1898878
#2   F    bold     low 0.1444913 0.1350332 0.1545403
#3   F    bold     mid 0.1575283 0.1514317 0.1638322
#7   F     shy    high 0.1510222 0.1396508 0.1632529
#8   F     shy     low 0.1238331 0.1155723 0.1326312
#9   F     shy     mid 0.1346557 0.1294648 0.1400283


# Travel (1) -> Search (2)  -----------------------------------------------

state1 <- 1
state2 <- 2
nb_states <- 3

# Get transition data - MALES
m.means <- unlist(lapply(means.m, '[[', nb_states * (state2 - 1) + state1)) 
m.lb <- unlist(lapply(lb.m, '[[', nb_states * (state2 - 1) + state1))
m.ub <- unlist(lapply(ub.m, '[[', nb_states * (state2 - 1) + state1))

df.m <- data.frame(
  sex = "M",
  pers = cov.m$mean_BLUP_logit,
  mean = m.means,
  lower_bound = m.lb,
  upper_bound = m.ub
)

df.m$persCat <- ifelse(df.m$pers < quantiles$lowerQuan[quantiles$Sex == "M"], "shy", "intermediate")
df.m$persCat <- ifelse(df.m$pers >  quantiles$upperQuan[quantiles$Sex == "M"], "bold", as.character(df.m$persCat))

trans_travel.search.M <- ddply(df.m,.(sex, persCat), summarize, p_mean = mean(mean), p_low = mean(lower_bound), p_upp = mean(upper_bound))
trans_travel.search.M <- subset(trans_travel.search.M, persCat != "intermediate")
trans_travel.search.M

#sex persCat    p_mean     p_low     p_upp
#1   M    bold 0.1095753 0.1044349 0.1149403
#2   M     shy 0.1283524 0.1228351 0.1340826

# Get transition data - FEMALES
f.means <- unlist(lapply(means.f_speed, '[[', nb_states * (state2 - 1) + state1)) 
f.lb <- unlist(lapply(lb.f_speed, '[[', nb_states * (state2 - 1) + state1))
f.ub <- unlist(lapply(ub.f_speed,'[[',nb_states * (state2 - 1) + state1))
df.f <- data.frame(sex = "F", pers = cov.f_speed$mean_BLUP_logit, dir = cov.f_speed$WindDir,
                   wind = cov.f_speed$WindSp, mean = f.means, lower_bound = f.lb, upper_bound = f.ub)

df.f$persCat <- ifelse(df.f$pers < quantiles$lowerQuan[quantiles$Sex == "F"], "shy", "intermediate")
df.f$persCat <- ifelse(df.f$pers >  quantiles$upperQuan[quantiles$Sex == "F"], "bold", as.character(df.f$persCat))

# Summarise wind into different categories
df.f$windCat <- "low"
df.f$windCat[df.f$wind > 5 & df.f$wind <= 10] <- "mid"
df.f$windCat[df.f$wind > 10] <- "high"

# Summarise transition probabilities
trans_travel.search.F <- ddply(df.f,.(sex, persCat, windCat), summarize, p_mean = mean(mean), p_low = mean(lower_bound), p_upp = mean(upper_bound))
trans_travel.search.F <- subset(trans_travel.search.F, persCat != "intermediate")
trans_travel.search.F

#sex persCat windCat    p_mean     p_low     p_upp
#1   F    bold    high 0.1291420 0.1195488 0.1394973
#2   F    bold     low 0.1310245 0.1225873 0.1399945
#3   F    bold     mid 0.1302853 0.1253379 0.1354069
#7   F     shy    high 0.1227754 0.1142623 0.1319170
#8   F     shy     low 0.1380395 0.1293039 0.1472960
#9   F     shy     mid 0.1319140 0.1270367 0.1369551




# Get coefficients for figure 2 ---------------------------------------------------------

# ## Search (2) to travel (1) ## --------------------------------------------------

behaviour <- data.frame(state = c("Travel", "Search", "Rest"), filename = c("travel", "search", "rest"))
nb_states = 3

# Get coefficients
state1 <- 2
state2 <- 1

# MALES #
m.means <- unlist(lapply(means.m,'[[',nb_states*(state2-1)+state1)) 
m.lb <- unlist(lapply(lb.m,'[[',nb_states*(state2-1)+state1))
m.ub <- unlist(lapply(ub.m,'[[',nb_states*(state2-1)+state1))
df.m <- data.frame(sex = "M", pers = cov.m$mean_BLUP_logit, 
                   mean = m.means, lower_bound = m.lb, upper_bound = m.ub)

df.m$persCat <- ifelse(df.m$pers < quantiles$lowerQuan[quantiles$Sex == "M"], "shy", "intermediate")
df.m$persCat <- ifelse(df.m$pers >  quantiles$upperQuan[quantiles$Sex == "M"], "bold", as.character(df.m$persCat))

df.m$state1 <- behaviour$state[state1]
df.m$state2 <- behaviour$state[state2]

# Take mean of responses for each personality/wind combination for plotting
df.m.2_1 <- df.m %>% 
  dplyr::group_by(persCat, state1, state2) %>%
  dplyr::summarise(mean = mean(mean, na.rm = T),
                   lower_bound  = mean(lower_bound , na.rm = T),
                   upper_bound  = mean(upper_bound , na.rm = T)) %>%
  data.frame()

# FEMALES #
f.means <- unlist(lapply(means.f_speed,'[[',nb_states*(state2-1)+state1)) 
f.lb <- unlist(lapply(lb.f_speed,'[[',nb_states*(state2-1)+state1))
f.ub <- unlist(lapply(ub.f_speed,'[[',nb_states*(state2-1)+state1))
df.f <- data.frame(sex = "F", pers = cov.f_speed$mean_BLUP_logit, 
                   wind=cov.f_speed$WindSp, mean = f.means, lower_bound = f.lb, upper_bound = f.ub)

df.f$persCat <- ifelse(df.f$pers < quantiles$lowerQuan[quantiles$Sex == "F"], "shy", "intermediate")
df.f$persCat <- ifelse(df.f$pers >  quantiles$upperQuan[quantiles$Sex == "F"], "bold", as.character(df.f$persCat))

# Take mean of responses for each personality/wind combination for plotting
df.f.2_1 <- df.f %>% 
  dplyr::group_by(wind, persCat) %>%
  dplyr::summarise(mean = mean(mean, na.rm = T),
                   lower_bound  = mean(lower_bound , na.rm = T),
                   upper_bound  = mean(upper_bound , na.rm = T)) %>%
  data.frame()



# ## Travel (1) to search (2) ## --------------------------------------------------

# Get coefficients
state1 <- 1
state2 <- 2

# MALES #
m.means <- unlist(lapply(means.m,'[[',nb_states*(state2-1)+state1)) 
m.lb <- unlist(lapply(lb.m,'[[',nb_states*(state2-1)+state1))
m.ub <- unlist(lapply(ub.m,'[[',nb_states*(state2-1)+state1))
df.m <- data.frame(sex = "M", pers = cov.m$mean_BLUP_logit, 
                   mean = m.means, lower_bound = m.lb, upper_bound = m.ub)

df.m$persCat <- ifelse(df.m$pers < quantiles$lowerQuan[quantiles$Sex == "M"], "shy", "intermediate")
df.m$persCat <- ifelse(df.m$pers >  quantiles$upperQuan[quantiles$Sex == "M"], "bold", as.character(df.m$persCat))

df.m$state1 <- behaviour$state[state1]
df.m$state2 <- behaviour$state[state2]


# Take mean of responses for each personality/wind combination for plotting
df.m.1_2 <- df.m %>% 
  dplyr::group_by(persCat, state1, state2) %>%
  dplyr::summarise(mean = mean(mean, na.rm = T),
                   lower_bound  = mean(lower_bound , na.rm = T),
                   upper_bound  = mean(upper_bound , na.rm = T)) %>%
  data.frame()

# FEMALES #
f.means <- unlist(lapply(means.f_speed,'[[',nb_states*(state2-1)+state1)) 
f.lb <- unlist(lapply(lb.f_speed,'[[',nb_states*(state2-1)+state1))
f.ub <- unlist(lapply(ub.f_speed,'[[',nb_states*(state2-1)+state1))
df.f <- data.frame(sex = "F", pers = cov.f_speed$mean_BLUP_logit, 
                   wind=cov.f_speed$WindSp, mean = f.means, lower_bound = f.lb, upper_bound = f.ub)

df.f$persCat <- ifelse(df.f$pers < quantiles$lowerQuan[quantiles$Sex == "F"], "shy", "intermediate")
df.f$persCat <- ifelse(df.f$pers >  quantiles$upperQuan[quantiles$Sex == "F"], "bold", as.character(df.f$persCat))

# Take mean of responses for each personality/wind combination for plotting
df.f.1_2 <- df.f %>% 
  dplyr::group_by(wind, persCat) %>%
  dplyr::summarise(mean = mean(mean, na.rm = T),
                   lower_bound  = mean(lower_bound , na.rm = T),
                   upper_bound  = mean(upper_bound , na.rm = T)) %>%
  data.frame()


# ## Search (2) to search (2) ## --------------------------------------------------

# Get coefficients
state1 <- 2
state2 <- 2

# MALES #
m.means <- unlist(lapply(means.m,'[[',nb_states*(state2-1)+state1)) 
m.lb <- unlist(lapply(lb.m,'[[',nb_states*(state2-1)+state1))
m.ub <- unlist(lapply(ub.m,'[[',nb_states*(state2-1)+state1))
df.m <- data.frame(sex = "M", pers = cov.m$mean_BLUP_logit, 
                   mean = m.means, lower_bound = m.lb, upper_bound = m.ub)

df.m$persCat <- ifelse(df.m$pers < quantiles$lowerQuan[quantiles$Sex == "M"], "shy", "intermediate")
df.m$persCat <- ifelse(df.m$pers >  quantiles$upperQuan[quantiles$Sex == "M"], "bold", as.character(df.m$persCat))

df.m$state1 <- behaviour$state[state1]
df.m$state2 <- behaviour$state[state2]


# Take mean of responses for each personality/wind combination for plotting
df.m.2_2 <- df.m %>% 
  dplyr::group_by(persCat, state1, state2) %>%
  dplyr::summarise(mean = mean(mean, na.rm = T),
                   lower_bound  = mean(lower_bound , na.rm = T),
                   upper_bound  = mean(upper_bound , na.rm = T)) %>%
  data.frame()

# FEMALES #
f.means <- unlist(lapply(means.f_speed,'[[',nb_states*(state2-1)+state1)) 
f.lb <- unlist(lapply(lb.f_speed,'[[',nb_states*(state2-1)+state1))
f.ub <- unlist(lapply(ub.f_speed,'[[',nb_states*(state2-1)+state1))
df.f <- data.frame(sex = "F", pers = cov.f_speed$mean_BLUP_logit, 
                   wind=cov.f_speed$WindSp, mean = f.means, lower_bound = f.lb, upper_bound = f.ub)

df.f$persCat <- ifelse(df.f$pers < quantiles$lowerQuan[quantiles$Sex == "F"], "shy", "intermediate")
df.f$persCat <- ifelse(df.f$pers >  quantiles$upperQuan[quantiles$Sex == "F"], "bold", as.character(df.f$persCat))

# Take mean of responses for each personality/wind combination for plotting
df.f.2_2 <- df.f %>% 
  dplyr::group_by(wind, persCat) %>%
  dplyr::summarise(mean = mean(mean, na.rm = T),
                   lower_bound  = mean(lower_bound , na.rm = T),
                   upper_bound  = mean(upper_bound , na.rm = T)) %>%
  data.frame()




# Build the male plot -----------------------------------------------------

plotdf.M <- rbind(df.m.2_1, df.m.2_2, df.m.1_2)
plotdf.M$transition <- paste(plotdf.M$state1, "to", plotdf.M$state2, sep = " ")

maleTransitionsPlot <- ggplot(subset(plotdf.M, persCat != "intermediate"), aes(x = transition, y = mean, group = persCat)) +
  geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound, group = persCat), width = 0.25, 
                position = position_dodge(width = 0.6)) +
  geom_point(aes(fill = persCat), pch = 21, size = 2.7, position = position_dodge(width = 0.6)) +
  scale_fill_manual(name = "", values = c(bold_col, shy_col)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(y = "Transition probability", tag = "(d)") +
  ggtitle("Male transition estimates") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.tag = element_text(size = 20),
        plot.title = element_text(size = 20, face = "bold")) 



# Build the female plots --------------------------------------------------

# SEARCH TO SEARCH
plot.f.2_2 <-
  ggplot(subset(df.f.2_2, persCat != "intermediate"), aes(x = wind, y=mean)) +
  geom_ribbon(size = 1, linetype = "blank", aes(ymin=lower_bound, ymax=upper_bound, col = persCat), alpha=0.15) + 
  geom_line(size = 1, aes(col = persCat)) + 
  scale_color_manual(name = "", labels = c("Bolder", "Shyer"), values = c(bold_col, shy_col)) +
  labs(y = "Transition probability", x = expression(paste("Wind speed (", ms^-1, ")", sep = "")),
       tag = "(a)") +
  ggtitle("Search to Search") +
  scale_x_continuous(limits=c(0, 20)) +
  scale_y_continuous( limits = c(0.5,1)) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.x=element_text(size = 18),
        axis.title.y=element_text(size = 18),
        legend.position = c(0.1, 0.95),
        legend.key.size = unit(0.75, "cm"),
        legend.key = element_blank(), 
        legend.text = element_text(size = 12),
        legend.box = "horizontal",
        legend.background = element_rect(fill = "transparent"),
        plot.tag = element_text(size = 20),
        plot.title = element_text(size = 20, face = "bold")) +
  guides(color = guide_legend(override.aes=list(fill=NA)))  

# SEARCH TO TRAVEL

plot.f.2_1 <-
  ggplot(subset(df.f.2_1, persCat != "intermediate"), aes(x = wind, y=mean)) +
  geom_ribbon(size = 1, linetype = "blank", aes(ymin=lower_bound, ymax=upper_bound, col = persCat), alpha=0.15) + 
  geom_line(size = 1, aes(col = persCat)) + 
  scale_color_manual(name = "", labels = c("Bolder", "Shyer"), values = c(bold_col, shy_col)) +
  labs(y = "Transition probability", x = expression(paste("Wind speed (", ms^-1, ")", sep = "")),
       tag = "(b)") +
  ggtitle("Search to Travel") +
  scale_x_continuous(limits=c(0, 20)) +
  scale_y_continuous( limits = c(0,0.3)) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.x=element_text(size = 18),
        axis.title.y=element_blank(),
        legend.position = "none",
        plot.tag = element_text(size = 20),
        plot.title = element_text(size = 20, face = "bold")) +
  guides(color = guide_legend(override.aes=list(fill=NA)))  

# TRAVEL TO SEARCH

plot.f.1_2 <-
  ggplot(subset(df.f.1_2, persCat != "intermediate"), aes(x = wind, y=mean)) +
  geom_ribbon(size = 1, linetype = "blank", aes(ymin=lower_bound, ymax=upper_bound, col = persCat), alpha=0.15) + 
  geom_line(size = 1, aes(col = persCat)) + 
  scale_color_manual(name = "", labels = c("Bolder", "Shyer"), values = c(bold_col, shy_col)) +
  labs(y = "Transition probability", x = expression(paste("Wind speed (", ms^-1, ")", sep = "")),
       tag = "(c)") +
  ggtitle("Travel to Search") +
  scale_x_continuous(limits=c(0, 20)) +
  scale_y_continuous( limits = c(0,0.3)) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.x=element_text(size = 18),
        axis.title.y=element_text(size = 18),
        legend.position = "none",
        plot.tag = element_text(size = 20),
        plot.title = element_text(size = 20, face = "bold")) +
  guides(color = guide_legend(override.aes=list(fill=NA)))  


#### FIGURE 2: TRANSITION ESTIMATES ------------------------------------------

fig2 <- ggarrange(plot.f.2_2, plot.f.2_1,
          plot.f.1_2, maleTransitionsPlot,
          align = "hv",
          nrow = 2,
          ncol = 2)

png(filename = paste0("Figures/FIG2_Transitions_quantile", label, ".png"), width = 12, height = 12, units = "in", res = 700)
fig2
dev.off()



# Plot all female transition estimates by speed, male by boldness for Figure S4 --------------------------------------------

behaviour <- data.frame(state = c("Travel", "Search", "Rest"), filename = c("travel", "search", "rest"))

nb_states = 3
labels <- c(F = "Female", M = "Male")

for (i in 1:3) {
  
  for (j in 1:3) {
    
    state1 <- i; print(i)
    state2 <- j; print(j) 
    
    ### MALE COEFFICIENTS
    
    m.means <- unlist(lapply(means.m,'[[',nb_states*(state2 - 1) + state1)) 
    m.lb <- unlist(lapply(lb.m,'[[',nb_states*(state2 - 1) + state1))
    m.ub <- unlist(lapply(ub.m,'[[',nb_states*(state2 - 1) + state1))
    df.m <- data.frame(sex = "M", pers = cov.m$mean_BLUP_logit, 
                       wind = cov.m$WindSp, mean = m.means, lower_bound = m.lb, upper_bound = m.ub)
    
    df.m$persCat <- ifelse(df.m$pers < quantiles$lowerQuan[quantiles$Sex == "M"], "shy", "intermediate")
    df.m$persCat <- ifelse(df.m$pers >  quantiles$upperQuan[quantiles$Sex == "M"], "bold", as.character(df.m$persCat))
    
    ### FEMALE COEFFICIENTS
    
    f.means <- unlist(lapply(means.f_speed,'[[',nb_states*(state2 - 1) + state1)) 
    f.lb <- unlist(lapply(lb.f_speed,'[[',nb_states*(state2 - 1) + state1))
    f.ub <- unlist(lapply(ub.f_speed,'[[',nb_states*(state2 - 1) + state1))
    df.f <- data.frame(sex = "F", pers = cov.f_speed$mean_BLUP_logit, 
                       wind = cov.f_speed$WindSp, mean = f.means, lower_bound = f.lb, upper_bound = f.ub)
    
    df.f$persCat <- ifelse(df.f$pers < quantiles$lowerQuan[quantiles$Sex == "F"], "shy", "intermediate")
    df.f$persCat <- ifelse(df.f$pers >  quantiles$upperQuan[quantiles$Sex == "F"], "bold", as.character(df.f$persCat))
    
    ### COMBINE DATAFRAMES
    
    all.df <- rbind(df.m, df.f)
    all.df$sex <- as.factor(as.character(all.df$sex))
    all.df$persCat <- as.factor(as.character(all.df$persCat))
    
    # Take mean of responses for each personal/sex/wind combination for plotting
    stat.df <- all.df %>% 
      dplyr::group_by(sex, wind, persCat) %>%
      dplyr::summarise(mean = mean(mean, na.rm = T),
                       lower_bound  = mean(lower_bound , na.rm = T),
                       upper_bound  = mean(upper_bound , na.rm = T)) %>%
      data.frame()
    
    ### BUILD PLOT
    
    female_plot <- ggplot(subset(stat.df, sex == "F"), aes(x = wind, y=mean)) +
      geom_ribbon(size = 1, linetype = "blank", aes(ymin=lower_bound, ymax=upper_bound, col = persCat), alpha=0.15) + 
      geom_line(size = 1, aes(col = persCat)) + 
      scale_color_manual(name = "", labels = c("Bolder", "Intermediate", "Shyer"), 
                         values = c(bold_col, inter_col, shy_col)) +
      labs(y = "Transition probability", x =  expression(paste("Wind speed (", ms^-1, ")", sep = ""))) +
      ggtitle("Females") +
      scale_y_continuous( limits = c(0,0.3)) +
      theme_bw() + 
      theme(axis.text.x=element_text(size=16), 
            axis.text.y=element_text(size=16), 
            axis.title.x=element_text(size = 18),
            axis.title.y=element_text(size = 18),
            legend.position = c(0.25, 0.25),
            legend.key.size = unit(1, "cm"),
            legend.key = element_blank(), 
            legend.text = element_text(size = 16),
            legend.box = "horizontal",
            legend.background = element_rect(fill = "transparent"),
            plot.title = element_text(hjust = 0.5, size = 18)) +
      guides(color = guide_legend(override.aes=list(fill=NA)))
    
    male_plot <- ggplot(subset(stat.df, sex == "M"), aes(x = persCat, y=mean)) +
      geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound, group = persCat), width = 0.25, position = position_dodge(width = 0.6)) +
      geom_point(aes(fill = persCat), pch = 21, size = 2.7, position = position_dodge(width = 0.6)) +
      scale_fill_manual(name = "", labels = c("Bolder", "Intermediate", "Shyer"), 
                        values = c(bold_col, inter_col, shy_col)) +
      scale_x_discrete(labels = c("Bolder", "Shyer")) +
      scale_y_continuous(limits = c(0, 0.3)) +
      xlab("Boldness") +
      ggtitle("Males") +
      theme_bw() +
      theme(axis.text.x = element_text(size=16), 
            axis.text.y = element_blank(), 
            axis.title.x = element_text(size = 18),
            axis.title.y = element_blank(),
            legend.position = "none",
            plot.title = element_text(hjust = 0.5, size = 18)) 
    
    # Scale y as 0.5 - 1 for same behaviour -> same behaviour
    if (i == j) {
      female_plot <-
        female_plot + scale_y_continuous(breaks = scales::pretty_breaks(n = 5),
                                         limits = c(0.5, 1))
      male_plot <-
        male_plot + scale_y_continuous(breaks = scales::pretty_breaks(n = 5),
                                       limits = c(0.5, 1))
    }
    
    # Combine the plots
    ## Run this if printing plots individually
    #combinedFig <- ggarrange(female_plot, male_plot, 
    #                         nrow = 1, 
    #                         align = "h",
    #                         widths = c(1, 0.5))
    
    #finalFig <- annotate_figure(combinedFig, top = textGrob(paste0(behaviour$state[i], " to ", behaviour$state[j]), 
    #                                                       vjust = 0.5, gp = gpar(cex = 1.6)))
    
    # Assign to a filename
    #assign(paste0("speedPlot_", i, "_", j), finalFig)
    
    assign(paste0("F_speedPlot_", i, "_", j), female_plot)
    assign(paste0("M_speedPlot_", i, "_", j), male_plot)
    
    #png(paste0("Figures/Transition Estimates/SPEED_", behaviour$filename[i], "_", behaviour$filename[j], "quantile", label,".png"), 
    #    width = 10, height = 7, units = "in", res = 350)
    #print(finalFig)
    #dev.off()
    
  }
  
}



# Prepare subplots for Figure S4 ---------------------------------------------------------------


# ## Travel row ##  -------------------------------------------------------

# Travel to travel #
M_speedPlot_1_1 <- M_speedPlot_1_1 +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

F_speedPlot_1_1 <- F_speedPlot_1_1 +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

speedPlot_1_1 <- ggarrange(F_speedPlot_1_1, M_speedPlot_1_1, 
                           nrow = 1, 
                           align = "h",
                           widths = c(1, 0.5))

speedPlot_1_1 <- annotate_figure(speedPlot_1_1, top = textGrob(expression(bold("Travel to travel")), 
                                            vjust = 0.5, gp = gpar(cex = 1.6)))

# Travel to search #
M_speedPlot_1_2 <- M_speedPlot_1_2 +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

F_speedPlot_1_2 <- F_speedPlot_1_2 +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

speedPlot_1_2 <- ggarrange(F_speedPlot_1_2, M_speedPlot_1_2, 
                           nrow = 1, 
                           align = "h",
                           widths = c(1, 0.5))

speedPlot_1_2 <- annotate_figure(speedPlot_1_2, top = textGrob(expression(bold("Travel to search")), 
                                              vjust = 0.5, gp = gpar(cex = 1.6)))


# Travel to rest #
M_speedPlot_1_3 <- M_speedPlot_1_3 +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

F_speedPlot_1_3 <- F_speedPlot_1_3 +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

speedPlot_1_3 <- ggarrange(F_speedPlot_1_3, M_speedPlot_1_3, 
                           nrow = 1, 
                           align = "h",
                           widths = c(1, 0.5))

speedPlot_1_3 <- annotate_figure(speedPlot_1_3, top = textGrob(expression(bold("Travel to rest")), 
                                                               vjust = 0.5, gp = gpar(cex = 1.6)))


# ## Search row ##  -------------------------------------------------------

# Search to travel #
M_speedPlot_2_1 <- M_speedPlot_2_1 +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_blank())

F_speedPlot_2_1 <- F_speedPlot_2_1 +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_blank(),
        legend.position = "none")

speedPlot_2_1 <- ggarrange(F_speedPlot_2_1, M_speedPlot_2_1, 
                           nrow = 1, 
                           align = "h",
                           widths = c(1, 0.5))

speedPlot_2_1 <- annotate_figure(speedPlot_2_1, top = textGrob(expression(bold("Search to travel")), 
                                                               vjust = 0.5, gp = gpar(cex = 1.6)))

# Search to search #
M_speedPlot_2_2 <- M_speedPlot_2_2 +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_blank())

F_speedPlot_2_2 <- F_speedPlot_2_2 +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_blank())

speedPlot_2_2 <- ggarrange(F_speedPlot_2_2, M_speedPlot_2_2, 
                           nrow = 1, 
                           align = "h",
                           widths = c(1, 0.5))

speedPlot_2_2 <- annotate_figure(speedPlot_2_2, top = textGrob(expression(bold("Search to search")), 
                                                               vjust = 0.5, gp = gpar(cex = 1.6)))


# Search to rest #
M_speedPlot_2_3 <- M_speedPlot_2_3 +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_blank())

F_speedPlot_2_3 <- F_speedPlot_2_3 +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_blank())

speedPlot_2_3 <- ggarrange(F_speedPlot_2_3, M_speedPlot_2_3, 
                           nrow = 1, 
                           align = "h",
                           widths = c(1, 0.5))

speedPlot_2_3 <- annotate_figure(speedPlot_2_3, top = textGrob(expression(bold("Search to rest")), 
                                                               vjust = 0.5, gp = gpar(cex = 1.6)))

# ## Rest row ##  -------------------------------------------------------

# Rest to travel #
M_speedPlot_3_1 <- M_speedPlot_3_1 +
  theme(plot.title = element_blank())

F_speedPlot_3_1 <- F_speedPlot_3_1 +
  theme(plot.title = element_blank(),
        legend.position = "none")

speedPlot_3_1 <- ggarrange(F_speedPlot_3_1, M_speedPlot_3_1, 
                           nrow = 1, 
                           align = "h",
                           widths = c(1, 0.5))

speedPlot_3_1 <- annotate_figure(speedPlot_3_1, top = textGrob(expression(bold("Rest to travel")), 
                                                               vjust = 0.5, gp = gpar(cex = 1.6)))

# Rest to search #
M_speedPlot_3_2 <- M_speedPlot_3_2 +
  theme(plot.title = element_blank())

F_speedPlot_3_2 <- F_speedPlot_3_2 + ylim(0,0.5) +
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_blank())

speedPlot_3_2 <- ggarrange(F_speedPlot_3_2, M_speedPlot_3_2, 
                           nrow = 1, 
                           align = "h",
                           widths = c(1, 0.5))

speedPlot_3_2 <- annotate_figure(speedPlot_3_2, top = textGrob(expression(bold("Rest to search")), 
                                                               vjust = 0.5, gp = gpar(cex = 1.6)))


# Rest to rest #
M_speedPlot_3_3 <- M_speedPlot_3_3 +
  theme(plot.title = element_blank())

F_speedPlot_3_3 <- F_speedPlot_3_3 +
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_blank())

speedPlot_3_3 <- ggarrange(F_speedPlot_3_3, M_speedPlot_3_3, 
                           nrow = 1, 
                           align = "h",
                           widths = c(1, 0.5))

speedPlot_3_3 <- annotate_figure(speedPlot_3_3, top = textGrob(expression(bold("Rest to rest")), 
                                                               vjust = 0.5, gp = gpar(cex = 1.6)))



#### FIGURE S4 -------------------------------------------------

png(paste0("Figures/FIGS4_quantile", label, ".png"), width = 20, height = 14, units = "in", res = 350)
ggarrange(speedPlot_1_1, speedPlot_1_2, speedPlot_1_3,
          speedPlot_2_1, speedPlot_2_2, speedPlot_2_3,
          speedPlot_3_1, speedPlot_3_2, speedPlot_3_3,
          nrow = 3,
          ncol = 3,
          widths = c(1, 0.925, 0.925))
dev.off()
