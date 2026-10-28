# confirmation that my code for extraction of tau_i posteriors from MCMC results works: 
# I can reproduce the betas & gammas that come out of their "summary" functions

source("./skiptrack/data_clean/bcfr_basicfit.R")


# results from skipTrack()
skipTrack.results(fitbcfr_basic, burnIn = 1500) 
# same as
summary(fitbcfr_basic, burnIn = 1500)

# vs.

# my code for tau extraction, modified to extract gammas
store_gammalists = vector(mode = "list", length = n_chains*(n_reps-(burnin-1))) # burnin-1 because 1502:3001 starts at 1502, not 1503, & goes to 3001, not 30001
# for-loop
index = 1
for (chain in 1:n_chains){
  for (rep in burnin:n_reps){
    store_gammalists[[index]] = fitbcfr_basic$fit[[chain]][[rep]]$Gamma[1,1]
    index = index + 1
  }
}
# combine
gammalist_df = as.data.frame(store_gammalists, col.names = paste0("run", 1:9000)) # 9000 = 1500 x 6 chains
MCMC_gammas = cbind("Individual" = fitbcfr_basic$fit[[1]][[1502]]$iDat[["Individual"]], gammalist_df)
# get quantiles for credible intervals
MCMC_gammas_post = MCMC_gammas |> rowwise()  |> 
  mutate(
    mean_gamma = mean(c_across(2:ncol(MCMC_gammas))), # excluding first column, "Individual"
    q025_gamma = quantile(c_across(2:ncol(MCMC_gammas)), 0.025), # 95% credible interval, ish
    median_gamma = quantile(c_across(2:ncol(MCMC_gammas)), 0.5),
    q975_gamma = quantile(c_across(2:ncol(MCMC_gammas)), 0.975)) |> 
  ungroup() |>   
  select(Individual, mean_gamma, q025_gamma, median_gamma, q975_gamma)

# returns the same values! as mean, 95% credible interval bounds 



skipTrack.diagnostics(fitbcfr_basic, param = c('tauis'), burnIn = 1500)
skipTrack.diagnostics(fitbcfr_basic, param = c('cijs'), burnIn = 1500)


# tau_i posterior extraction code -- bcfr_basicfit.R

# store_taulists = vector(mode = "list", length = n_chains*(n_reps-(burnin-1))) # burnin-1 because 1502:3001 starts at 1502, not 1503, & goes to 3001, not 30001
# 
# index = 1
# for (chain in 1:n_chains){
#   for (rep in burnin:n_reps){
#     df_chain_rep = fitbcfr_basic$fit[[chain]][[rep]]$iDat
#     taus = df_chain_rep[["taus"]]
#     store_taulists[[index]] = taus
#     index = index + 1
#   }
# }
# 
# taulist_df = as.data.frame(store_taulists, col.names = paste0("run", 1:9000)) # 9000 = 1500 x 6 chains
# MCMC_taus = cbind("Individual" = fitbcfr_basic$fit[[1]][[1502]]$iDat[["Individual"]], taulist_df)
# 
# # now: get some features of the POSTERIOR dist
# MCMC_taus_post = MCMC_taus |> rowwise()  |> 
#   mutate(
#     mean_tau = mean(c_across(2:ncol(MCMC_taus))), # excluding first column, "Individual"
#     q025_tau = quantile(c_across(2:ncol(MCMC_taus)), 0.025), # 95% credible interval, ish
#     median_tau = quantile(c_across(2:ncol(MCMC_taus)), 0.5),
#     q975_tau = quantile(c_across(2:ncol(MCMC_taus)), 0.975)) |> 
#   ungroup() |>   
#   select(Individual, mean_tau, q025_tau, median_tau, q975_tau)


