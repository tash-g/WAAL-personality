### CITATION: This code is modified from scripts shared in Clay et al. 2020, J. Anim. Ecol.

### AIM: lot predicted stationary probabilities for a range of specified covariates 


# PREAMBLE ----------------------------------------------------------------

# Define the packages
packages <- c("ggplot2", "momentuHMM", "dplyr", "patchwork", "cowplot", "grid",
              "ggpubr")

# Install packages not yet installed - change lib to library path
#installed_packages <- packages %in% rownames(installed.packages())

#if (any(installed_packages == FALSE)) {
#  install.packages(packages[!installed_packages], lib = "C:/Users/libraryPath")
#}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))


### Create outputs and figures folders if they don't currently exist
out.path <- "./Data_outputs/"
if (dir.exists(out.path) == FALSE) {
  dir.create(out.path)
}

figures.path <- "./Figures/Stationary estimates/"
if (dir.exists(figures.path) == FALSE ){
  dir.create(figures.path)
}


source("ci_stationary.R")

# Set colours
shy_col <- "#00DD2F"
bold_col <- "purple"
inter_col <- "#be254e"

# SET BOLDNESS QUANTILES ----------------------------------------------

## Set personality quantiles
upperQ <- 0.9  # 0.75
lowerQ <- 0.1 # 0.25

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


# LOAD MODEL OUTPUTS ------------------------------------------------------

# load in best models for males and females (based on AIC table)
file.in <- paste0("./Data_outputs/", paste0("F_mod_", 4, ".RData"))
load(file = file.in)
f.mod <- mod.p # set to model

file.in <- paste0("./Data_outputs/", paste0("M_mod_", 9, ".RData"))
load(file = file.in)
m.mod <- mod.p # set to model



# FEMALES : CALCULATE STATIONARY PROBABILITIES --------------------------------------

### FEMALES: prob ~ wind * personality
### Calculate for crosswind

# Calculate speed range
min(f.mod$data$WindSp) # 0.05
max(f.mod$data$WindSp) # 20.7

# Get covariate data
ws.f <- seq(from = 0.05, 20, by = 0.25)

# Subsample personality to speed up processing
pers.f <- unique(f.mod$data$mean_BLUP_logit)
pers.f <- pers.f[order(pers.f)]
pers.f <- pers.f[c(TRUE,FALSE)]

comb.f <- expand.grid(ws.f, pers.f)

cov.f <- data.frame(LoD = "L", WindDir = 90, WindSp = comb.f$Var1, mean_BLUP_logit = comb.f$Var2)
cov.f$LoD <- as.factor(cov.f$LoD)

## Get stationary probs (gives for each row of covariates dataset)
s.f <- as.data.frame(stationary(model = f.mod, covs = cov.f))

# Create confidence intervals around stationary probabilities
stat.f <- ci.stationary(model = f.mod, cov = cov.f, alpha = 0.95) 


## Combine into dataframe for plotting
prob_cov.f <- cbind(s.f, cov.f)

prob.f <- tidyr::gather(prob_cov.f, state, prob, travel:rest, factor_key = TRUE)
prob.f$lwr <- tidyr::gather(stat.f[[1]], state, lwr, travel:rest, factor_key = T)$lwr
prob.f$upr <- tidyr::gather(stat.f[[2]], state, lwr, travel:rest, factor_key = T)$lwr
prob.f$sex <- "F"
prob.f$pers <- ifelse(prob.f$mean_BLUP_logit < quantiles$lowerQuan[quantiles$Sex == "F"], "shy", "intermediate")
prob.f$pers <- ifelse(prob.f$mean_BLUP_logit >  quantiles$upperQuan[quantiles$Sex == "F"], "bold", as.character(prob.f$pers))

prob.f[,c("sex","pers")] <- lapply(prob.f[,c("sex","pers")], as.factor)


# Take mean of responses for personality alone
prob.f2 <- prob.f %>% 
  dplyr::group_by(pers, state, sex, LoD) %>%
  dplyr::summarise(prob = mean(prob, na.rm = T),
                   lwr = mean(lwr, na.rm = T),
                   upr = mean(upr, na.rm = T)) %>%
  data.frame()

# Take mean of responses for each personal/sex/wind combination for plotting
prob.f <- prob.f %>% 
  group_by(pers, state, sex, WindSp, WindDir, LoD) %>%
  summarise(prob = mean(prob, na.rm = T),
            lwr = mean(lwr, na.rm = T),
            upr = mean(upr, na.rm = T)) %>%
  data.frame()

prob.f$pers_state <- as.factor(as.character(paste(prob.f$pers, prob.f$state, sep = "_")))


# save(prob.f, file = paste0("Data_outputs/stationary_probs_all_quantile", label, "_F.RData"))
# load(paste0("Data_outputs/stationary_probs_all_quantile", label, "_F.RData"))

# save(prob.f2, file = paste0("Data_outputs/stationary_probs_pers_quantile", label, "_F.RData"))
# load(paste0("Data_outputs/stationary_probs_pers_quantile", label, "_F.RData"))

# Get effect sizes for stationary probability by wind
prob.f$windCat <- "low"
prob.f$windCat[prob.f$WindSp > 5 & prob.f$WindSp <= 10] <- "mid"
prob.f$windCat[prob.f$WindSp > 10] <- "high"

stationary.wind.F <- ddply(prob.f,.(pers, state, windCat), 
                           summarize, p_mean = mean(prob), p_low = mean(lwr), p_upp = mean(upr))
stationary.wind.F

#pers  state windCat     p_mean      p_low      p_upp
#1          bold travel    high 0.54652501 0.51707562 0.57552235
#2          bold travel     low 0.37641369 0.34813811 0.40570327
#3          bold travel     mid 0.46774762 0.45170091 0.48386946
#4          bold search    high 0.39725680 0.37164703 0.42352153
#5          bold search     low 0.32924858 0.31200292 0.34701594
#6          bold search     mid 0.37585715 0.36456586 0.38729010
#7          bold   rest    high 0.05621819 0.04832534 0.06550171
#8          bold   rest     low 0.29433773 0.26555690 0.32476745
#9          bold   rest     mid 0.15639523 0.14516132 0.16833244
#10 intermediate travel    high 0.53016348 0.51260161 0.54760640
#11 intermediate travel     low 0.35463715 0.33794651 0.37174323
#12 intermediate travel     mid 0.44482621 0.43475526 0.45494697
#13 intermediate search    high 0.41500252 0.39947535 0.43072785
#14 intermediate search     low 0.36002012 0.34860243 0.37161434
#15 intermediate search     mid 0.40396802 0.39648513 0.41149688
#16 intermediate   rest    high 0.05483400 0.04971906 0.06050570
#17 intermediate   rest     low 0.28534272 0.26749604 0.30385675
#18 intermediate   rest     mid 0.15120576 0.14384273 0.15887729
#19          shy travel    high 0.51724397 0.49005501 0.54421978
#20          shy travel     low 0.33753176 0.31288270 0.36323257
#21          shy travel     mid 0.42678964 0.41230012 0.44141980
#22          shy search    high 0.42901615 0.40476973 0.45367866
#23          shy search     low 0.38485915 0.36658637 0.40349800
#24          shy search     mid 0.42628329 0.41510023 0.43754341
#25          shy   rest    high 0.05373988 0.04686622 0.06170371
#26          shy   rest     low 0.27760909 0.25136247 0.30541865
#27          shy   rest     mid 0.14692707 0.13706965 0.15736817

stationary.summary.F <- ddply(prob.f,.(pers, state), 
                           summarize, p_mean = mean(prob), p_low = mean(lwr), p_upp = mean(upr))
stationary.summary.F

#           pers  state    p_mean     p_low     p_upp
#1         bold travel 0.4843028 0.4584976 0.5101544
#2         bold search 0.3749048 0.3549657 0.3953373
#3         bold   rest 0.1407923 0.1268422 0.1560258
#4 intermediate travel 0.4649476 0.4494762 0.4804758
#5 intermediate search 0.3984983 0.3860096 0.4111417
#6 intermediate   rest 0.1365541 0.1276942 0.1459364
#7          shy travel 0.4497023 0.4263232 0.4732730
#8          shy search 0.4172937 0.3978065 0.4370997
#9          shy   rest 0.1330040 0.1205411 0.1465486

# MALES : CALCULATE STATIONARY PROBABILITIES --------------------------------------

### MALES: prob ~ personality
### Calculate for crosswind and average speed

# Calculate mean speed
mean(m.mod$data$WindSp) # 9.169809

## Get covariate data
pers.m <- unique(m.mod$data$mean_BLUP_logit)
pers.m <- pers.m[order(pers.m)]
pers.m <- pers.m[c(TRUE,FALSE)]

cov.m <- data.frame(LoD = "L", WindDir =  90, WindSp = 9, mean_BLUP_logit = pers.m)

## Get stationary probabilities
s.m <- as.data.frame(stationary(model = m.mod, covs = cov.m))

# Create confidence intervals around stationary probabilities
stat.m <- ci.stationary(model = m.mod, cov = cov.m, alpha = 0.95) 

## Combine into a dataframe
prob_cov.m <- cbind(s.m, cov.m)
prob.m <- tidyr::gather(prob_cov.m, state, prob, travel:rest, factor_key = TRUE)

prob.m$lwr <- tidyr::gather(stat.m[[1]], state, lwr, travel:rest, factor_key = T)$lwr
prob.m$upr <- tidyr::gather(stat.m[[2]], state, lwr, travel:rest, factor_key = T)$lwr

prob.m$sex <- rep("M", nrow(prob.m))
prob.m$pers <- ifelse(prob.m$mean_BLUP_logit < quantiles$lowerQuan[quantiles$Sex == "M"], "shy", "intermediate")
prob.m$pers <- ifelse(prob.m$mean_BLUP_logit >  quantiles$upperQuan[quantiles$Sex == "M"], "bold", as.character(prob.m$pers))


prob.m[,c(9,10)] <- lapply(prob.m[,c(9,10)], as.factor)

# Take mean of responses for each personality for plotting
prob.m <- prob.m %>% 
  dplyr::group_by(pers, state, sex, WindSp, WindDir, LoD) %>%
  dplyr::summarise(prob = mean(prob, na.rm = T),
            lwr = mean(lwr, na.rm = T),
            upr = mean(upr, na.rm = T)) %>%
  data.frame()

prob.m$pers_state <- as.factor(as.character(paste(prob.m$pers, prob.m$state, sep = "_")))

# save(prob.m, file = paste0("Data_outputs/stationary_probs_quantile", label, "_M.RData"))
# load(paste0("Data_outputs/stationary_probs_quantile", label, "_M.RData"))

## Get male effect sizes

stationary.summary.M <- ddply(prob.m,.(pers, state), 
                              summarize, p_mean = mean(prob), p_low = mean(lwr), p_upp = mean(upr))
stationary.summary.M

#         pers  state     p_mean      p_low     p_upp
#1         bold travel 0.50972669 0.49302726 0.5263968
#2         bold search 0.39171455 0.37860817 0.4049872
#3         bold   rest 0.09855876 0.09052762 0.1072313
#4 intermediate travel 0.46483615 0.45463930 0.4750619
#5 intermediate search 0.43436634 0.42608214 0.4426887
#6 intermediate   rest 0.10079751 0.09519659 0.1066920
#7          shy travel 0.42507827 0.41015641 0.4401450
#8          shy search 0.47283361 0.46026121 0.4854382
#9          shy   rest 0.10208812 0.09459745 0.1101068


# PLOT FEMALE STATIONARY PROBS BY WIND SPEED AND PERSONALITY ---------------------------------

figure3_female <- ggplot(subset(prob.f, pers != "intermediate"), aes(WindSp, prob)) + 
  geom_ribbon(size = 1.3, aes(ymin = lwr, ymax = upr, fill = pers_state), alpha = 0.15) +
  scale_fill_manual(values = c("grey0", "grey0", "grey0", "grey30", "grey30", "grey30", "grey60", "grey60", "grey60")) +
  geom_line(aes(colour = pers, linetype = state), size = 1.3) + 
  scale_linetype_manual(name = "", values = c("solid", "dotted", "dashed"),labels = c("Travel", "Search", "Rest")) +
  scale_color_manual(name = "", labels = c("Bolder", "Shyer"), values = c(bold_col, shy_col)) +
  labs(y = "Stationary probability", x = expression(paste("Wind speed (", ms^-1, ")", sep = "")),
       tag = "(a)") +
  ggtitle("Females") + 
  scale_x_continuous(limits = c(0, 20)) +
  ylim(0,0.7) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 16), 
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.position = c(0.275, 0.87),
        legend.background = element_rect(fill = "transparent"),
        legend.box = "horizontal",
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.key = element_blank(), 
        legend.text = element_text(size = 14),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 16),
        plot.title = element_text(size = 20, face = "bold"),
        plot.tag = element_text(size = 20)) +
  guides(fill = "none",
         colour = guide_legend(order = 2),
         linetype = guide_legend(order = 1))

figure3_female.dens <- ggplot(data = f.mod$data, aes(x = WindSp)) + 
  geom_histogram(color = "black", fill = "white") + 
  theme_void() 

figure3_female.final <- figure3_female + 
  annotation_custom(ggplotGrob(figure3_female.dens), 
                    xmin = -1, xmax = 21, ymin = -0.038, ymax = 0.06)



# PLOT FEMALE STATIONARY PROBS BY PERSONALITY -----------------------------

figure3_female2 <- ggplot(subset(prob.f2, pers != "intermediate"), aes(state, prob)) + 
  geom_errorbar(aes(ymin = lwr, ymax = upr, group = pers), width = 0.25, position = position_dodge(width = 0.6)) +
  geom_point(aes(fill = pers), pch = 21, size = 2.7, position = position_dodge(width = 0.6)) +
  scale_fill_manual(name = "", labels = c("Bolder", "Shyer"), values = c(bold_col, shy_col)) +
  scale_x_discrete(labels = c("Travel", "Search", "Rest")) +
  scale_y_continuous(limits = c(0, 0.7)) +
  labs(y = "Stationary probability", x = "Behaviour", tag = "(b)") +
  ggtitle("Females") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 16), 
        axis.text.y = element_blank(), 
        axis.title.x = element_text(size = 18),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        plot.tag = element_text(size = 20)) 


# PLOT MALE STATIONARY PROBS BY PERSONALITY --------------------------------

figure3_male <- ggplot(subset(prob.m, pers != "intermediate"), aes(state, prob)) + 
  geom_errorbar(aes(ymin = lwr, ymax = upr, group = pers), width = 0.25, position = position_dodge(width = 0.6)) +
  geom_point(aes(fill = pers), pch = 21, size = 2.7, position = position_dodge(width = 0.6)) +
  scale_fill_manual(name = "", labels = c("Bolder", "Shyer"), values = c(bold_col, shy_col)) +
  scale_x_discrete(labels = c("Travel", "Search", "Rest")) +
  scale_y_continuous(limits = c(0, 0.7)) +
  labs(y = "Stationary probability", x = "Behaviour", title = "Males") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 16), 
        axis.text.y = element_blank(), 
        axis.title.x = element_text(size = 18),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 
  


#### FIGURE 4: PLOT STATIONARY PROBABILITiES FOR MALES AND FEMALES ----------------------------------------------------------------

fig4 <- ggarrange(figure3_female.final, figure3_female2, figure3_male, 
                  nrow = 1, 
                  align = "h",
                  widths = c(1, 0.5, 0.44))

png(filename = paste0("Figures/FIG4_stationaryProbs_quantile", label, ".png"),
    width = 12, height = 6, units = "in", res = 700)
fig4
dev.off()




# FEMALES & MALES: STATIONARY PROBS AS A FUNCTION OF WIND SPEED & DIRECTION -----------------------------

### Get covariate data

## Females
wd.f <- c(0, 90, 180)
ws.f <- seq(from = 0.05, 20, by = 0.2)
comb.f <- expand.grid(wd.f, ws.f)

cov.f <- data.frame(LoD = "L", WindDir = comb.f$Var1, WindSp = comb.f$Var2, mean_BLUP_logit = mean(f.mod$data$mean_BLUP_logit))
cov.f$LoD <- as.factor(cov.f$LoD)

# Males
wd.m <- c(0, 90, 180)
ws.m <- seq(from = 0.05, 24, by = 0.2)
comb.m <- expand.grid(wd.m, ws.m)

cov.m <- data.frame(LoD = "L", WindDir = comb.m$Var1, WindSp = comb.m$Var2, mean_BLUP_logit = mean(m.mod$data$mean_BLUP_logit))
cov.m$LoD <- as.factor(cov.m$LoD)

# Get stationary probs (gives for each row of covariates dataset)
s.f <- as.data.frame(stationary(model = f.mod, covs = cov.f))
s.m <- as.data.frame(stationary(model = m.mod, covs = cov.m))

# Create confidence intervals around stationary probabilities
stat.f <- ci.stationary(model = f.mod, cov = cov.f, alpha = 0.95) 
stat.m <- ci.stationary(model = m.mod, cov = cov.m, alpha = 0.95) 


## Combine into dataframe for plotting
prob_cov.f <- cbind(s.f, cov.f)
prob_cov.m <- cbind(s.m, cov.m)

prob.f <- tidyr::gather(prob_cov.f, state, prob, travel:rest, factor_key = TRUE)
prob.m <- tidyr::gather(prob_cov.m, state, prob, travel:rest, factor_key = TRUE)

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
stat.df$sex <- as.factor(stat.df$sex)

# save(stat.df, file = paste0("Data_outputs/stationary_probs_WIND.RData"))
# load(paste0("Data_outputs/stationary_probs_WIND.RData"))


# PLOT STATIONARY PROBABILITIES WITH WIND SPEED AND DIRECTION ---------------------------------

stat.df$sex <- factor(stat.df$sex, labels = c("Females", "Males"))
stat.df$WindDir <- as.factor(as.character(stat.df$WindDir))

travel_col <- "#648FFF"
search_col <- "#DC267F"
rest_col <- "#FFB000"

### HEADWIND

headwindPlot <- ggplot(subset(stat.df, WindDir == 180), aes(x = WindSp, y = prob)) + 
  facet_wrap(.~sex, labeller = label_value) + 
  geom_line(aes(linetype = state, col = state), size = 1.3) + 
  geom_ribbon(size = 1.3, aes(ymin = lwr, ymax = upr, fill = state), alpha = 0.15) +
  scale_fill_manual(values = c("grey30", "grey30", "grey30")) +
  scale_linetype_manual(values = c("solid", "dotted", "dashed"),labels = c("Travel", "Search", "Rest"), name = "Behaviour") +
  scale_colour_manual(values = c(travel_col, search_col, rest_col), labels = c("Travel", "Search", "Rest"), name = "Behaviour") +
  lims(y = c(0,1), x = c(0,23)) +
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = c(0.085, 0.85),
        legend.background = element_rect(fill = "transparent"),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 16)) +
  guides(fill = "none")

headwindPlot <-  ggdraw() +   
  draw_plot(headwindPlot) +
  draw_image(
    "Figures/headwind_schematic.png", x = 0.75, y = 0.9, hjust = 1, vjust = 1, halign = 1, valign = 1,
    width = 0.28, height = 0.28 )

### CROSSWIND

crosswindPlot <- ggplot(subset(stat.df, WindDir == 90), aes(x = WindSp, y = prob)) + 
  facet_wrap(.~sex, labeller = label_value) + 
  geom_line(aes(linetype = state, col = state), size = 1.3) + 
  geom_ribbon(size = 1.3, aes(ymin = lwr, ymax = upr, fill = state), alpha = 0.15) +
  scale_fill_manual(values = c("grey30", "grey30", "grey30")) +
  scale_linetype_manual(values = c("solid", "dotted", "dashed"),labels = c("Travel", "Search", "Rest"), name = "Behaviour") +
  scale_colour_manual(values = c(travel_col, search_col, rest_col), labels = c("Travel", "Search", "Rest"), name = "Behaviour") +
  lims(y = c(0,1), x = c(0,23)) +
  theme_bw() + 
  labs(y = "Stationary probability") +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  guides(fill = "none")

crosswindPlot <-  ggdraw() +   
  draw_plot(crosswindPlot) +
  draw_image(
    "Figures/crosswind_schematic.png", x = 0.76, y = 0.97, hjust = 1, vjust = 1, halign = 1, valign = 1,
    width = 0.3, height = 0.3 )

### TAILWIND

tailwindPlot <- ggplot(subset(stat.df, WindDir == 0), aes(x = WindSp, y = prob)) + 
  facet_wrap(.~sex, labeller = label_value) + 
  geom_line(aes(linetype = state, col = state), size = 1.3) + 
  geom_ribbon(size = 1.3, aes(ymin = lwr, ymax = upr, fill = state), alpha = 0.15) +
  scale_fill_manual(values = c("grey30", "grey30", "grey30")) +
  scale_linetype_manual(values = c("solid", "dotted", "dashed"),labels = c("Travel", "Search", "Rest"), name = "Behaviour") +
  scale_colour_manual(values = c(travel_col, search_col, rest_col), labels = c("Travel", "Search", "Rest"), name = "Behaviour") +
  lims(y = c(0,1), x = c(0,23)) +
  theme_bw() + 
  labs(y = "Stationary probability", x = expression(paste("Wind speed (", ms^-1, ")", sep = ""))) +
  theme(axis.text.x = element_text(size = 16), 
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 18),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  guides(fill = "none")

tailwindPlot <-  ggdraw() +   
  draw_plot(tailwindPlot) +
  draw_image(
    "Figures/tailwind_schematic.png", x = 0.76, y = 0.97, hjust = 1, vjust = 1, halign = 1, valign = 1,
    width = 0.28, height = 0.28 )


#### FIGURE S3 : INTERACTION PLOT OF WIND ON STATIONARY PROBS  ------------------------------------------------------

myfig <- ggarrange(headwindPlot, crosswindPlot, tailwindPlot,
                   heights = c(0.9, 0.85, 1),
                   ncol = 1,
                   nrow = 3)

png("Figures/FIGS3_stationaryProbsByWind.png", width = 7, height = 14, units = "in", res = 350)
annotate_figure(myfig, left = textGrob("Stationary probability", rot = 90, vjust = 0.5, gp = gpar(cex = 1.6)))
dev.off()


#  PREP FOR FIGURE S5 : PLOT FEMALE STATIONARY PROBS BY WIND SPEED AND PERSONALITY ---------------------------------

figureS5_female <- ggplot(prob.f, aes(WindSp, prob)) + 
  geom_ribbon(size = 1.3, aes(ymin = lwr, ymax = upr, fill = pers_state), alpha = 0.15) +
  scale_fill_manual(values = c("grey0", "grey0", "grey0", "grey30", "grey30", "grey30", "grey60", "grey60", "grey60")) +
  geom_line(aes(colour = pers, linetype = state), size = 1.3) + 
  scale_linetype_manual(name = "", values = c("solid", "dotted", "dashed"),labels = c("Travel", "Search", "Rest")) +
  scale_color_manual(name = "", labels = c("Bolder", "Intermediate", "Shyer"), 
                     values = c(bold_col, inter_col, shy_col)) +
  labs(y = "Stationary probability", x = expression(paste("Wind speed (", ms^-1, ")", sep = "")),
       tag = "(a)") +
  scale_x_continuous(limits = c(0, 20)) +
  ylim(0,0.7) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 16), 
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.position = c(0.3, 0.87),
        legend.background = element_rect(fill = "transparent"),
        legend.box = "horizontal",
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.key = element_blank(), 
        legend.text = element_text(size = 14),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 16),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        plot.tag = element_text(size = 20)) +
  guides(fill = "none",
         colour = guide_legend(order = 2),
         linetype = guide_legend(order = 1))

figureS5_female.dens <- ggplot(data = f.mod$data, aes(x = WindSp)) + 
  geom_histogram(color = "black", fill = "white") + 
  theme_void() 

figureS5_female.final <- figureS5_female + 
  annotation_custom(ggplotGrob(figureS5_female.dens), 
                    xmin = -1, xmax = 21, ymin = -0.038, ymax = 0.06)


# PREP FOR FIGURE S5 : PLOT FEMALE STATIONARY PROBS BY PERSONALITY -----------------------------

figureS5_female2 <- ggplot(prob.f2, aes(state, prob)) + 
  geom_errorbar(aes(ymin = lwr, ymax = upr, group = pers), width = 0.25, position = position_dodge(width = 0.6)) +
  geom_point(aes(fill = pers), pch = 21, size = 2.7, position = position_dodge(width = 0.6)) +
  scale_fill_manual(name = "", labels = c("Bolder", "Intermediate", "Shyer"), 
                    values = c(bold_col, inter_col, shy_col)) +
  scale_x_discrete(labels = c("Travel", "Search", "Rest")) +
  scale_y_continuous(limits = c(0, 0.7)) +
  labs(y = "Stationary probability", x = "Behaviour", tag = "(b)") +
  ggtitle("Females") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 16), 
        axis.text.y = element_blank(), 
        axis.title.x = element_text(size = 18),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        plot.tag = element_text(size = 20)) 


# PREP FOR FIGURE S5 : PLOT MALE STATIONARY PROBS BY PERSONALITY --------------------------------

figureS5_male <- ggplot(prob.m, aes(state, prob)) + 
  geom_errorbar(aes(ymin = lwr, ymax = upr, group = pers), width = 0.25, position = position_dodge(width = 0.6)) +
  geom_point(aes(fill = pers), pch = 21, size = 2.7, position = position_dodge(width = 0.6)) +
  scale_fill_manual(name = "", labels = c("Bolder", "Intermediate", "Shyer"), 
                    values = c(bold_col, inter_col, shy_col)) +
  scale_x_discrete(labels = c("Travel", "Search", "Rest")) +
  scale_y_continuous(limits = c(0, 0.7)) +
  labs(y = "Stationary probability", x = "Behaviour", title = "Males") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 16), 
        axis.text.y = element_blank(), 
        axis.title.x = element_text(size = 18),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 



#### FIGURE S5: PLOT STATIONARY PROBABILITiES FOR MALES AND FEMALES WHOLE POPULATION----------------------------------------------------------------

figS5 <- ggarrange(figureS5_female.final, figureS5_female2, figureS5_male, 
                  nrow = 1, 
                  align = "h",
                  widths = c(1, 0.5, 0.44))

png(filename = paste0("Figures/FIGS5_stationaryProbs_quantile", label, ".png"),
    width = 12, height = 6, units = "in", res = 700)
figs5
dev.off()


