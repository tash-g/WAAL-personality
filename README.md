# Boldness predicts individual plasticity in flight responses to winds
Natasha Gillies, Henri Weimerskirch, Jack Thorley, Thomas A. Clay, Lucía Martina Martín López, Rocio Joo, Mathieu Basille, Samantha C. Patrick

## Overview
This repository contains scripts and data to recreate the main results and figures of this paper (currently in prep). 

## Scripts
A short description of each script is given below.

- **1_HMM_construct-models.R** This script is used to fit Hidden Markov Models (HMMs) with every possible combination of boldness and boldness by wind interaction covariates. AIC is used to compare the relative support for each. The code here fits the models, outputs the AIC for each, as well as autocorrelation and goodness-of-fit plots. _Note: depending on processing power, each model may take several hours to compile and so this script will take several days to run to completion._
- **2_HMM_stationary_probabilities.R** This script takes the two best-supported models (for males and females respectively) and uses these to plot stationary probabilities corresponding to Figure 3 in the manuscript and Figure S4 in the Supplementary Information. 
- **3_HMM_transition_estimates.R** This script takes the two best-supported models (for males and females respectively) and uses these to plot transition estimates corresponding to Figure 2 in the manuscript and Figures S1-S3 in the Supplementary Information. 
- **4_HMM_validation.R** This script contains code to run our two validation methods, described in the Methods and Supplementary Information. In part 1)  the behavioural output of our best-supported HMMs is compared with expert classification of behaviours from foraging tracks. In part 2) the value of each covariate in the best supported model is randomly shuffled 50 times, and support for this randomly shuffled model compared to the 'real' model. _Note: depending on processing power, each randomised model in part 2 may take several hours to compile and so this script will take several days to run to completion._
-  **SUPP-INF_plot-GPS.R** This script contains code to analyse the effect of boldness estimates on fitness as a function of bird age. 
- **SUPP-INF_plot-GPS.R** This script contains code to process the GPS data and produce Figures 1 and S1, which show the distribution of birds according to their boldness estimate and across years. 
- **SUPP-INF_wind_availability_by_personality.R** This script contains code to run the analyses described in Supplementary Information, which compare wind speeds and wind directions experienced by birds as a function of their boldness estimate and across years. 

## Data inputs

These data are used in the above scripts. Note that all Rings/BirdIDs have been recoded and so cannot be linked to existing datasets. Please contact the authors if you would like to make use of these datasets, as we may be able to offer additional information, data, or advice. 

- **WAAL_GPS_2010-2021_personality_wind.csv** This is the main dataset for analysis. Each row corresponds to an individual GPS fix. The columns are as follows:
  -  _Ring_: Factor encoding unique ID of bird
  -  _Sex_ : Factor encoding male (M) or female (F)
  -  _Year_ : Factor encoding year of tracking, between 2010-2020
  -  _Longitude_ : Longitudinal position of GPS fix; numeric
  -  _Latitude_ : Latitudunal position of GPS fix; numeric
  -  _DateTime_ : Date and time of fix to nearest second
  -  _WindSp_ : Wind speed in m/s; numeric
  -  _Dev.wind2_ : Wind direction relative to bird trajectory; numeric
  -  _LoD_ : Factor denoting daytime (L) or night (D)
  -  _Trip_bout_ : Factor encoding individual trip identity
  -  _mean_BLUP_logit_ : Numeric boldness score for that individual

- **WAAL_gps_2020-2021_manualStates.csv** Dataset giving expert classification of each GPS fix into a different behavioural state. Used for model validation described in Methods and Supplementary Information. The columns are as follows:
  - _BirdID_ : Factor encoding unique ID of bird
  - _TripID_ : Factor encoding individual trip identity
  - _Sex_ : Factor encoding male (M) or female (F)
  - _Year_ : Factor encoding year of tracking, between 2010-2020
  - _DateTime_ : Date and time of fix to nearest second
  -  _Longitude_ : Longitudinal position of GPS fix; numeric
  - _Latitude_ : Latitudunal position of GPS fix; numeric
  - _STATE_ : Manual classification of behavior, where 1 = travel, 2 = search, 3 = rest
