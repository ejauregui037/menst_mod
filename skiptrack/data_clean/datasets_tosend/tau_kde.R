
# create a spline/fxn to sample from (for MIME sampling to capture uncertainty)

# from tau_i posteriors from MCMC with 10,000 its

# load fitted model -------------------------------------------------------------------------------
# load fitted skipTrack basic model w/ 10,000 iterations to MCMC:
newdf = readRDS(file = "./skiptrack/data_clean/fitbcfr_basic_rep10")
fitbcfr_basic = newdf # rename so that code below does not need to be modified

# calculate taus & posterior stats -------------------------------------------------------------------------------
# get the metrics that are created in bcfr_basicfit10.R, but not saved in the RDS file [which is just the fitted skipTrack object]:

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
# 51001/6 = 8500.167: i.e. 10,000 runs - 1,500 burn-in

# get kernel densities -------------------------------------------------------------------------------
# now: rather than getting FEATURES, use the full draws to estimate smoothed density...
# built-in fxn density() is not easy to sample from. -- test = density(as.numeric(MCMC_taus[1,2:51001])); plot(test)
# However, ks (kde -- kernel density estimation) package *is*

library(ks)
# test = kde(as.numeric(MCMC_taus[1,2:51001])) # normal dist kernels by default; leaving default bandwidth (H,h)
# plot(test)
# hist(as.numeric(MCMC_taus[1,2:51001]), breaks = 100) # visual confirmation: density dist looks correct
# rkde(n=2, test) # & this generates samples!

# note: kde() and density() give the SAME results with default settings! 
# https://bookdown.org/egarpor/NP-UC3M/kde-i-ks.html


# for our purposes:
# save a list of kde objects for sampling, where each is attached to an individual ID (named list): this way, it is easy to iterate through and sample! 
# I will generate a function for this as well

# initialize list for storage of smoothed densities (of tau_i, for each individual): 
MCMC_kdes = vector(mode = "list", length = nrow(MCMC_taus))
# reduce computation in loop
ncol_MCMC = ncol(MCMC_taus)
for (row_num in 1:nrow(MCMC_taus)){
  kde_obj = kde(as.numeric(MCMC_taus[row_num, 2:ncol_MCMC])) # first column is ID #, so start at col #2
  MCMC_kdes[[row_num]] = list("ID" = MCMC_taus[row_num, 1], "kde_obj" = kde_obj)
}
# name list items
names(MCMC_kdes) = paste0("id_", MCMC_taus[,1]) # this is the naming convention! accessed in the function below -- avoid subsetting issues with pure #s as names to list objects

# save this object as RDS [send to Eva]
# saveRDS(MCMC_kdes, "./skiptrack/data_clean/datasets_tosend/tau_kdes")

# set up sampling fxn -----------------------------------------------------------------------------
# Plus: how to sample from this object efficiently... for a given individual 
# i.e. to generate a new "dataset" for MIME s.t. you sample a tau_i value for a given individual i from their kernel density, the kde object
# assume that you have a df with PersonID

# SET A SEED BEFORE RUNNING THIS FUNCTION****
mime_samp_tau_kde = function(n_mime, p_id_vector){
  # mime_samp_tau_kde: Function is designed to promote easy sampling from the kde for MIME estimation. 
    # in particular, this function is tailor-made to support MIME estimation of coefficient uncertainty with # n_mime
    # samples from the tau_i posteriors, since tau_i is the outcome for each individual i in a set of regressions.
  # Arguments:
    # n_mime: Number of samples for MIME estimation (double)
    # p_id_vector: List of IDs corresponding to PersonID in Clue data (atomic vector--integers)
  # Returns:
    # A list of length = n_mime, where each entry is a df with two columns: PersonID & tau_i_samp
      # that is, n_mime samples of tau_i for each individual, from each of their individual smoothed kernel density estimates for their respective tau_i
  
  # first portion: generate n = n_mime samples of tau_i for each individual i, from their own kernel density distribution for tau_i
  p_id_vector1 = paste0("id_", p_id_vector) # avoid subsetting issues with pure #s as names to list objects
  list_individual_samps = vector(mode = "list", length = length(p_id_vector))
  for (ii_p_id in 1:length(p_id_vector1)){
    p_id = p_id_vector1[ii_p_id] # get individual ID (PersonID, with "id_" pasted on at the front)
    kde_obj1 = MCMC_kdes[[p_id]]$kde_obj # get fitted kernel density for individual i
    tau_i_samp = rkde(n=n_mime, kde_obj1) # get n_mime samples
    list_individual_samps[[ii_p_id]] = as.data.frame(list("PersonID" = p_id_vector[ii_p_id], "tau_i_samp" = tau_i_samp))
  }
  
  # current form: a list of length = length(p_id_vector) (i.e. # of individuals i), where each df has nrow() = n_mime; each df in the list corresponds to one individual with n_mime tau_i samples)
  # convert to: a list of length = n_mime (i.e. # of tau_i samples for each individual), where each df has nrow() = length(p_id_vector); each df in the list corresponds to one MIME sample with a tau_i value for each individual i
  list_mime_samples = vector(mode = "list", length = n_mime)
  for (nrow_mime in 1:n_mime){ # this is also the nrow() of each df in the list "list_individual_samps"
    list_mime_samples[[nrow_mime]] = do.call(rbind,
                                             lapply(list_individual_samps, function(x){x[nrow_mime,]}))
  }
  
  return(list_mime_samples)
}

# test = mime_samp_tau_kde(10, c(10001450, 10004453, 10011345))
 
# then: just merge these dfs with the df containing other covars (by = "PersonID") in a loop (over 1:n_mime),
# & you will easily be able to get uncertainty estimates for the coefficients!

# this function ^^ runs very quickly : )

