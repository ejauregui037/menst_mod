# Documentation & function to enable MIME sampling with tau_i posteriors

# ------------------------------------------------------------------------------
# Densities of tau_i posteriors (one for each individual, n=56, in the Clue dataset)

# Note: you can name this object whatever you like
tau_kdes = readRDS(file = "[....]/tau_kdes") # fill in filepath to wherever you downloaded & saved the file
# Format:
  # List of length n = 56
  # Names of list objects correspond to PersonID tag/identifier for individuals in Clue dataset
names(tau_kdes) # this will give you all of the ID #s
  # Each list object contains two elements...
tau_kdes[[1]]$ID # PersonID (integer)
tau_kdes[[1]]$kde_obj # kde object (kde), i.e. the smoothed kernel density of the tau_i posterior distribution for individual i [i.e. for that given PersonID]

# ------------------------------------------------------------------------------
# Custom function to sample from these tau_i densities:

# The goal is to sample from the tau_i density function (smoothed) for each individual.  -- i.e. to generate a new "dataset" for MIME s.t. you sample a tau_i value for a given individual i from their kernel density, the kde object

# 1. Set a seed before running this function. It includes random sampling, so we want it to be reproducible.
set.seed(958049)

# 2. Load in the df you are using to run regressions with the Clue dataset.
XX = read.csv("...")
# That df ("XX") should have PersonID as a column. Then:
# Get your PersonIDs as a vector
p_id_vector = XX$PersonID

# 3. HERE IS THE CUSTOM FUNCTION:

mime_samp_tau_kde = function(n_mime, p_id_vector){
  
  # mime_samp_tau_kde: Function is designed to promote easy sampling from the kde for MIME estimation. 
    # in particular, this function is tailor-made to support MIME estimation of coefficient uncertainty with # n_mime
    # samples from the tau_i posteriors, since tau_i is the outcome for each individual i in a set of regressions.
  # Arguments:
    # n_mime: Number of samples for MIME estimation (type double)
    # p_id_vector: List of IDs corresponding to PersonID in Clue data (type atomic vector--integers)
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
  
  # second portion: reformat
    # current form: a list of length = length(p_id_vector) (i.e. # of individuals i), where each df has nrow() = n_mime; each df in the list corresponds to one individual with n_mime tau_i samples)
    # convert to: a list of length = n_mime (i.e. # of tau_i samples for each individual), where each df has nrow() = length(p_id_vector); each df in the list corresponds to one MIME sample with a tau_i value for each individual i
  list_mime_samples = vector(mode = "list", length = n_mime)
  for (nrow_mime in 1:n_mime){ # this is also the nrow() of each df in the list "list_individual_samps"
    list_mime_samples[[nrow_mime]] = do.call(rbind,
                                             lapply(list_individual_samps, function(x){x[nrow_mime,]}))
  }
  
  return(list_mime_samples)
  
}

# test case which works without having loaded any dataframe:
test = mime_samp_tau_kde(10, c(10001450, 10004453, 10011345)) # these are 3 PersonIDs, & it will generate 10 samples from their respective tau_i smoothed posteriors (kernel density fxns) for each of them

# 4. Using these for MIME:

# Just merge these dfs with the "XX" df containing other covars in a loop (over 1:n_mime)
# This will return n = n_mime dataframes, where all covars other than "tau_i_samp" are the same in EVERY df, and "tau_i_samp" differs for EVERY individual across ALL dfs.
# use (by = "PersonID") 

# Then, if you run your regression over all of these dfs, you can use the resulting list of coefficient estimates to get a 95% pseudo-credible interval, etc. per MIME
