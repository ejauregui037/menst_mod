# df to eva 10/22/2025

# Variable names and meanings:
# PersonID | cycle_num | cycle_len | bc | 
# cijbasic_lowerCI, cijbasic_mean, cijbasic_upperCI — just provide basic, since c_ijs were essentially the same between the two models
# cijbc_lowerCI, cijbc_mean, cijbc_upperCI
# variance (crude) | tau_i_basic | tau_i_bc | irr (self-reported irregularity) coded as 1 = irregular, 0 = regular)

# note that bc values are missing for individuals without bc data

setwd("~/Desktop/ColumbiaPhD/research_projects/Lauren_menstrualmodelingmethods/menst_mod")
library(dplyr)
library(tidyr)

# source "basicfit"
source("./skiptrack/data_clean/bcfr_basicfit.R")

# self-reported irregularity data
bcfr_irr = read.csv("./data/yws_clue/ywsclue_repro.csv") # BCFR irregularity data (Clue)
irr_id = bcfr_irr |>
  mutate(
    irregular = case_when( # recode irregularity as 0 & 1 -- 1 = irregular, 0 = regular
      IRREGULAR_MENSES == 2 ~ 0, # regular
      IRREGULAR_MENSES == 1 ~ 1 # irregular (Eva said 1 is irregular)
    )
  ) |> 
  select(PersonID, irregular)

# I am not including the BC-adjusted model in this round of the df -- 
# a) when I checked last time, the c_ij estimates were very similar
# b) the really useful variables (e.g. endometriosis, etc.) are not here yet!
# This also means I am leaving out the BC data values

# merge info:
colnames(MCMC_taus_post)[1] = "PersonID" # for merge
bcfr_yws_irr = bcfr_skptr # "bcfr yws irregularity info"
# create variance
bcfr_yws_irr = bcfr_yws_irr |> group_by(PersonID) |> 
  mutate(cycle_var = var(cycle_len)) |> 
  ungroup()
# merge
bcfr_yws_irr = merge(bcfr_yws_irr, MCMC_taus_post, by = "PersonID") # taus
bcfr_yws_irr = cbind(bcfr_yws_irr, cijs_obs) # cijs
bcfr_yws_irr = merge(bcfr_yws_irr, irr_id, by = "PersonID") # self-reported irr

# short version
bcfr_yws_irr_short = bcfr_yws_irr |> 
  select(PersonID, cycle_len, cycle_num, cycle_var, 
         mean_tau, cij_mean, irregular)

# write to csv
write.csv(bcfr_yws_irr, "./skiptrack/data_clean/datasets_tosend/bcfryws_irr.csv")
write.csv(bcfr_yws_irr_short, "./skiptrack/data_clean/datasets_tosend/bcfryws_irr_s.csv")
