
# Create a df with the info we want on each individual tau_i, c_ij, & variability measures
# i.e. from skipTrack fit with 10,000 MCMC iterations

# setup
# setwd("~/Desktop/ColumbiaPhD/research_projects/Lauren_menstrualmodelingmethods/menst_mod")  
library(dplyr)
library(tidyr)
library(skipTrack)

# load skipTrack-formatted data -------------------------------------------------------------------------------
source("./skiptrack/clean_clue.R")
# df formatted for skipTrack is named bcfr_skptr, has columns:
print(colnames(bcfr_skptr))

# load fitted model -------------------------------------------------------------------------------
# load fitted skipTrack basic model w/ 10,000 iterations to MCMC:
newdf = readRDS(file = "./skiptrack/fitting/bcfrskip10")
fitbcfr_basic = newdf # rename so that code below does not need to be modified

# calculate taus & posterior stats -------------------------------------------------------------------------------
# get the metrics that are created in bcfr_basicfit10.R, but not saved in the RDS file [which is just the fitted skipTrack object]:

# get precision by Person ID -- output as a df
n_chains = length(fitbcfr_basic$fit) # 6 chains
n_reps = length(fitbcfr_basic$fit[[1]]) # 10,001 reps (+1 starting rep)
burnin = 1502 # +1 for initial guess

store_taulists = vector(mode = "list", length = n_chains*(n_reps-(burnin-1))) # burnin-1 because 1502:3001 starts at 1502, not 1503, & goes to 3001, not 30001

index = 1
for (chain in 1:n_chains){
  for (rep in burnin:n_reps){
    df_chain_rep = fitbcfr_basic$fit[[chain]][[rep]]$iDat
    taus = df_chain_rep[["taus"]]
    store_taulists[[index]] = taus
    index = index + 1
  }
}

taulist_df = as.data.frame(store_taulists, col.names = paste0("run", 1:51000)) # 9000 = 8500 x 6 chains
MCMC_taus = cbind("Individual" = fitbcfr_basic$fit[[1]][[1502]]$iDat[["Individual"]], taulist_df) # accessing IDs from a particular run, actual run accessed not important
# n=56 unique individuals

# now: get some features of the POSTERIOR dist
MCMC_taus_post = MCMC_taus |> rowwise()  |> 
  mutate(
    mean_tau = mean(c_across(2:ncol(MCMC_taus))), # excluding first column, "Individual"
    q025_tau = quantile(c_across(2:ncol(MCMC_taus)), 0.025), # 95% credible interval, ish
    median_tau = quantile(c_across(2:ncol(MCMC_taus)), 0.5),
    q975_tau = quantile(c_across(2:ncol(MCMC_taus)), 0.975)) |> 
  ungroup() |>   
  select(Individual, mean_tau, q025_tau, median_tau, q975_tau)

# ** note: tau_i G-R for tau is calculable with built-in skipTrack functions : ) 
skipTrack.diagnostics(fitbcfr_basic, param = c("tauis"), burnIn = 1500) # super low psrf, super high ESS! excellent

# get c_ijs
cijs_obs = skipTrack.results(fitbcfr_basic)$cijs # one c_ij for each observation, vs. tau_i for each individual
colnames(cijs_obs) = c("cij_lower", "cij_mean", "cij_upper")


# NOW: the main event ---------------------------------------------------------------------------
# create csv --------------------------------------------------------------------------------
# self-reported irregularity data
bcfr_irr = read.csv("./data/yws_clue/ywsclue_repro.csv") # BCFR irregularity data (Clue)
colnames(bcfr_irr)[1] = "personid"

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
  select(personid, irr_ever, irr_12mo)

# rename
bcfr_yws_irr = bcfr_skptr # "bcfr yws irregularity info"

# create variance
bcfr_yws_irr = bcfr_yws_irr |> group_by(personid) |> 
  mutate(cycle_var = var(cycle_len)) |> 
  ungroup()

# merge
colnames(MCMC_taus_post)[1] = "personid" # for merge
bcfr_yws_irr = merge(bcfr_yws_irr, MCMC_taus_post, by = "personid") # taus
bcfr_yws_irr = cbind(bcfr_yws_irr, cijs_obs) # cijs
bcfr_yws_irr = merge(bcfr_yws_irr, irr_id, by = "personid") # self-reported irr -- ever & last 12 mos

# save objects -----------------------------------------------------------------------------------------
# short version
bcfr_yws_irr_short = bcfr_yws_irr |> 
  select(personid, cycle_len, cycle_num, cycle_var, 
         mean_tau, cij_mean, 
         irr_ever, irr_12mo)

# write to csv
write.csv(bcfr_yws_irr, "./skiptrack/datasets_tosend/eva/bcfryws_skiptrack10k.csv")
write.csv(bcfr_yws_irr_short, "./skiptrack/datasets_tosend/eva/bcfryws_skiptrack10k_s.csv")

# # to compare -----------------------------------------------------------------------------------------
# # for visual comparison:
# # comparison of results with 3000 vs 10,000 iterations
# it_3000 = read.csv("./skiptrack/data_clean/datasets_tosend/bcfryws_irr_nospot1.csv")
# it_10000 = read.csv("./skiptrack/data_clean/datasets_tosend/bcfryws_irr_nospot1_10k.csv")
# # scanning through the dfs, results are REMARKABLY similar

