
# same format of dataframe as from 10/22/25 "df_to_eva", EXCEPT:
# now, spotting-only reports in Clue have not been counted as periods ("buffer" of one day -- see 'bcfr_skiptrack_nospot.R')

# skipTrack has been refit with these data (no covariates)
# setwd("~/Desktop/ColumbiaPhD/research_projects/Lauren_menstrualmodelingmethods/menst_mod")  
library(dplyr)
library(tidyr)

# fit skipTrack --------------------------------------------------------------------------------
# source "basicfit", but using df with spotting radius applied (i.e. spotting != period, not equivalent signals in data or in calculating cycle length)
source("./skiptrack/data_clean/exploratory/bcfr_basicfit_nospot.R")

# check diagnostics
skipTrack.diagnostics(fitbcfr_basic, param = 'cijs') # trace plot looks... just very very messy??? I think that's just high uncertainty in general...
skipTrack.results(fitbcfr_basic, burnIn = 1500) 

# create csv --------------------------------------------------------------------------------
# self-reported irregularity data
bcfr_irr = read.csv("./data/yws_clue/ywsclue_repro.csv") # BCFR irregularity data (Clue)

# recode & select variables of interest
irr_id = bcfr_irr |>
  mutate(
    # ever-irregular menses
    irr_ever = case_when(
      IRREGULAR_MENSES == 2 ~ 0, # regular
      IRREGULAR_MENSES == 1 ~ 1, # irregular (Eva said 1 is irregular)
      IRREGULAR_MENSES == 3 ~ NA # 3 = "don't know"..... does that mean regular? can't assume
    ), 
    # irregular in the last 12 months
    irr_12mo = case_when(
      MENSES_STARS_AS_EXPECTED == 2 ~ 1, # IRREGULAR -- opposite coding from EVER irregular!!!!
      MENSES_STARS_AS_EXPECTED == 1 ~ 0,
      MENSES_STARS_AS_EXPECTED == 3 ~ NA # 3 = "don't know"..... does that mean regular? can't assume
    )
  ) |> 
  select(PersonID, irr_ever, irr_12mo)

# rename
bcfr_yws_irr = bcfr_skptr # "bcfr yws irregularity info"

# create variance
bcfr_yws_irr = bcfr_yws_irr |> group_by(PersonID) |> 
  mutate(cycle_var = var(cycle_len)) |> 
  ungroup()

# merge
colnames(MCMC_taus_post)[1] = "PersonID" # for merge
bcfr_yws_irr = merge(bcfr_yws_irr, MCMC_taus_post, by = "PersonID") # taus
bcfr_yws_irr = cbind(bcfr_yws_irr, cijs_obs) # cijs
bcfr_yws_irr = merge(bcfr_yws_irr, irr_id, by = "PersonID") # self-reported irr -- ever & last 12 mos

# short version
bcfr_yws_irr_short = bcfr_yws_irr |> 
  select(PersonID, cycle_len, cycle_num, cycle_var, 
         mean_tau, cij_mean, 
         irr_ever, irr_12mo)

# write to csv
write.csv(bcfr_yws_irr, "./skiptrack/data_clean/datasets_tosend/bcfryws_irr_nospot1.csv")
write.csv(bcfr_yws_irr_short, "./skiptrack/data_clean/datasets_tosend/bcfryws_irr_nospot_s1.csv")
 