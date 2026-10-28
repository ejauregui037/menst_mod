
# Fit skipTrack to BCFR YWS data with 10,000 (10k) MCMC iterations, no covariates

# setup -----------------------------------------------------------------------------------------------
# setwd("~/Desktop/ColumbiaPhD/research_projects/Lauren_menstrualmodelingmethods/menst_mod")  
library(skipTrack)

# load skipTrack-formatted data -------------------------------------------------------------------------------
source("./skiptrack/clean_clue.R")
# df formatted for skipTrack is named bcfr_skptr, has columns:
print(colnames(bcfr_skptr))

# fit skipTrack -------------------------------------------------------------------------------
# hyper-parameters have been chosen based on previous runs/tests; however, as visible, convergence is not always great.
# burn-in chosen based on plots -- ~3000?
set.seed(9999)
fitbcfr_basic = skipTrack.fit(Y = bcfr_skptr$cycle_len,
                              cluster = bcfr_skptr$personid,
                              chains = 6,
                              reps = 10000)

# save output ----------------------------------------------------------------------
# saveRDS(fitbcfr_basic, file = "./skiptrack/fitting/bcfrskip10.rds")
test = readRDS(file = "./skiptrack/fitting/bcfrskip10.rds")









# post-fit processing ----------------------------------------------------------------------

# diagnostics
skipTrack.diagnostics(fitbcfr_basic, param = 'cijs') # trace plot looks... just very very messy??? I think that's just high uncertainty in general...
skipTrack.results(fitbcfr_basic, burnIn = 1500) 

# methodological note for future: potential issues with G-R statistc & proposed higher threshold than 1.1
# https://arxiv.org/pdf/1812.09384

# get precision by Person ID -- output as a df

n_chains = length(fitbcfr_basic$fit) # 6 chains
n_reps = length(fitbcfr_basic$fit[[1]]) # 3001 reps (starting rep)
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

taulist_df = as.data.frame(store_taulists, col.names = paste0("run", 1:9000)) # 9000 = 1500 x 6 chains
MCMC_taus = cbind("Individual" = fitbcfr_basic$fit[[1]][[1502]]$iDat[["Individual"]], taulist_df)

# now: get some features of the POSTERIOR dist
MCMC_taus_post = MCMC_taus |> rowwise()  |> 
  mutate(
    mean_tau = mean(c_across(2:ncol(MCMC_taus))), # excluding first column, "Individual"
    q025_tau = quantile(c_across(2:ncol(MCMC_taus)), 0.025), # 95% credible interval, ish
    median_tau = quantile(c_across(2:ncol(MCMC_taus)), 0.5),
    q975_tau = quantile(c_across(2:ncol(MCMC_taus)), 0.975)) |> 
  ungroup() |>   
  select(Individual, mean_tau, q025_tau, median_tau, q975_tau)

# ** note: I could calculate my own G-R for tau_i.... ***** could help decide whether these estimates are valid or not!! or look at other potential diagnostics...**

# get c_ijs
cijs_obs = skipTrack.results(fitbcfr_basic)$cijs # one c_ij for each observation, vs. tau_i for each individual
colnames(cijs_obs) = c("cij_lower", "cij_mean", "cij_upper")








