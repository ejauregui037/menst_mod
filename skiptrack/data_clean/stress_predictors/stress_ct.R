
# run skipTrack with stress as a covariate (predictor of cycle irregularity*)
# STRESS_CT = number of stressful events (cont) 

# setup ------------------------------------------------
# setwd("~/Desktop/ColumbiaPhD/research_projects/Lauren_menstrualmodelingmethods/menst_mod")  
library(skipTrack)

# load skipTrack-processed data ------------------------------------------------
source("./skiptrack/data_clean/exploratory/bcfr_skiptrack_nospot.R")
# rename df of interest s.t. code below does not need to be modified:
bcfr_skptr = bcfr_skptr_nospot1
rm(list=setdiff(ls(), "bcfr_skptr")) # remove extraneous objects from global envir

# data prep -------------------------------------------------
# load, merge, & format stressful events data
bcfr_strain = read.csv("./data/yws_clue/ywsclue_strain.csv") # BCFR YWS birth control data(Clue)
bcfr_strain = bcfr_strain[c("PersonID", "STRESS_CT")]
# new df
bcfr_skptr1 = merge(bcfr_skptr, bcfr_strain, by = "PersonID") # 994 obs; everyone has a value for this var! note that each individual has a STRESS_CT value, such that they are repeated for each observed cycle for each individual [like a baseline covar, not time-varying]
# however: there are 5 NAs [5 entries in this df] for STRESS_CT (none for STRESS_TH)... need to drop those
bcfr_skptr1 = bcfr_skptr1 |> drop_na() # 989 obs -- 5 entries were dropped

# fit model -------------------------------------------------
# assume (log)linearity(gamma dist)... will not impose any particular functional form in running this reg

# might take longer to converge now that we have a covar... try with 5000 to start
set.seed(9999)
fitbcfr_stressct = skipTrack.fit(Y = bcfr_skptr1$cycle_len,
                              cluster = bcfr_skptr1$PersonID,
                              # X = as.matrix(bcfr_skptr1[c("other", "pill")]), # covariates for cycle length mean -- left out "none" as "ref" i.e. so 1 df (i.e. intercept/when all covars =0)
                              Z = as.matrix(bcfr_skptr1["STRESS_CT"]), # covariates for cycle length precision i.e. 1/variance
                              chains = 6,
                              reps = 5000)

# diagnostics
skipTrack.diagnostics(fitbcfr_stressct, param = 'Gammas') # trace plot looks... just very very messy??? I think that's just high uncertainty in general...
skipTrack.results(fitbcfr_stressct, burnIn = 1500) 
# these plots look pretty good; chains are mixing... not as much as cij or taui, but...no trends, which is good

# save results
# saveRDS(fitbcfr_stressct, "./skiptrack/data_clean/stress_predictors/fitbcfr_stressct")

# examine beta estimates:
skipTrack.results(fitbcfr_stressct, burnIn = 1500)$Gammas
# shouldn't it be intercept & coeff?
# what is "Coverage tag"
# how to interpret these estimates?








