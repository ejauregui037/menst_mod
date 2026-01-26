# df to Eva 10/15/25

# skipTrack model with birth control as a predictor

# create a df with:
# ID, observed cycle #, observed cycle length, estimated c_ij (lowerCI, mean, higherCI), 
# trimmed to exclude extreme #s, since I think that would impact skipTrack performance (e.g. 1 cycle of 600 days)

# NOTE: tau_i = precision = 1/variance *** but with c_ijs i.e. reporting errors filtered out/inferred through hierarchical model...

# irregular = 1 is irregular!

# setup
setwd("~/Desktop/ColumbiaPhD/research_projects/Lauren_menstrualmodelingmethods/menst_mod")
library(skipTrack)
library(dplyr)
library(tidyr)
# data
bcfr_trk = read.csv("./data/yws_clue/ywsclue_track.csv") # BCFR tracking data (Clue)
bcfr_irr = read.csv("./data/yws_clue/ywsclue_repro.csv") # BCFR irregularity data (Clue)
# add birth control to the predictive model
bcfr_bc = read.csv("./data/yws_clue/ywsclue_bc.csv") # BCFR YWS birth control data(Clue)
bcfr_bc$bc_type = as.factor(bcfr_bc$BC_TYPE)
levels(bcfr_bc$bc_type)
table(bcfr_bc$bc_type) # very few people with implant and other... 1 each

# function to get cycle lengths

cycle_len_only <- function(df){
  
  # exclude bc data
  df_skip = df[df$category == "period" | df$category == "spotting" ,] # none are missing? odd but i guess by definition these are ENTRIES in the clue app -- also note this definition/criteria in writeup for Eva
  
  # get cycle lengths for each individual
  df_skip = df_skip |> group_by(PersonID) |> arrange(track_date_DtDemogr, .by_group = T)  |>    # Cmd Shift M for pipe
    mutate(gap_days = track_date_DtDemogr - lag(track_date_DtDemogr),
           new_cycle = if_else(gap_days > 4, 1, NA)) |>  # if difference between rows >4 -- note this condition in description to Eva
    arrange(track_date_DtDemogr, .by_group = T)  
  
  # ensure that the first entry for each Person ID is a "new cycle"
  df_skip = df_skip |> group_by(PersonID) |> 
    mutate(new_cycle = if_else(is.na(new_cycle) & track_date_DtDemogr == min(track_date_DtDemogr), 1, new_cycle)) 
  
  # calculate cycle length from new_cycle
  df_cycle = df_skip |> group_by(PersonID) |> arrange(track_date_DtDemogr, .by_group = T)  |> 
    filter(new_cycle == 1) |> 
    mutate(cycle_len = lead(track_date_DtDemogr) - track_date_DtDemogr)
  
  # for skiptrack, get cycle_num & eliminate extraneous variables
  df_skptr = df_cycle |>  select(PersonID, cycle_len) |> 
    group_by(PersonID) |> mutate(cycle_num = row_number()) |> 
    filter(!is.na(cycle_len))
  
  # return
  return(df_skptr)
}

# get BCFR clue data formatted for skiptrack
bcfr_skptr = cycle_len_only(bcfr_trk)

# ---------------------------------------------------------------------------------------------
# fit skiptrack with birth control (factor var) as a predictor of cycle length & regularity

# merged = merge(bcfr_skptr, bcfr_bc, by = "PersonID") # ! issue -- some individuals gave multiple values here, but not dates for when those values differed?? -- except "BC_DATE_DtDemogr" ??
# solution:
bcfr_bc1 = bcfr_bc |> group_by(PersonID) |> 
  mutate(bc = bc_type[which.max(BC_DATE_DtDemogr)]) # supposedly choose the most recent value of bc_type?
bcfr_bc1 = bcfr_bc1[c("PersonID", "bc")] # ID & bc type
bcfr_bc1 = unique(bcfr_bc1)
# new df
bcfr_skptr1 = merge(bcfr_skptr, bcfr_bc1, by = "PersonID") # only 945 observations with BC data... some individuals must not have had BC data...

# bc categories are too many--sparse data. Collapse them
table(bcfr_skptr1$bc) # 0 observations with vaginal ring... the rest of them seem fine? but how many INDIVIDUALS
table(unique(bcfr_skptr1[c("PersonID", "bc")])$bc) # this is by individual... I am going to group condoms, implant, & vaginalring in "other"
# recode
bcfr_skptr1 = bcfr_skptr1 |> group_by(PersonID) |> 
  mutate(bc = case_when(
    bc == "condoms" ~ "none", # because, for the purposes of cycle length/regularity, this is nothing hormonal 
    bc == "implant" ~ "other",
    bc == "IUD" ~ "other", # IUD will compose MOST of "other", except 1 "other" + 1 "implant" (there are zero vaginal rings)
    bc == "none" ~ "none",
    bc == "other" ~ "other",
    bc == "pill" ~ "pill",
    bc == "vaginal_ring" ~ "other"
  ))
table(unique(bcfr_skptr1[c("PersonID", "bc")])$bc) # worked

# convert BC category to indicator vars; skipTrack.fit does not take factors
bcfr_skptr1 = bcfr_skptr1 |> mutate(value = 1) |> 
  pivot_wider(names_from = bc, values_from = value, values_fill = 0)

# !!! need to clean out the crazy reporting individuals with nutty cycles/no data!
# exclude...
  # observations with... cycle length > 365
  # individuals with... < 3 cycles observed (not valid for variability/guessing whether cycle was skipped or not; would have to rely heavily on priors from others...TOO much hierarchicalpooling imo)

bcfr_skptr1 = bcfr_skptr1 |> filter(cycle_len < 365)
bcfr_skptr1 = bcfr_skptr1 |> group_by(PersonID) |> filter(n() >= 3) #only eliminates a few ppl! nice : 754 obs

# no other predictors

# -------------------------

# model with BC as a predictor

set.seed(9999)
fitbcfr_bc = skipTrack.fit(Y = bcfr_skptr1$cycle_len,
                           cluster = bcfr_skptr1$PersonID,
                           X = as.matrix(bcfr_skptr1[c("other", "pill")]), # covariates for cycle length mean -- left out "none" as "ref" i.e. so 1 df (i.e. intercept/when all covars =0)
                           Z = as.matrix(bcfr_skptr1[c("other", "pill")]), # covariates for cycle length precision i.e. 1/variance
                           chains = 6,
                           reps = 6000) # well over twice the burn-in... since it seems to be taking a while to converge

# 6 chains 3000 reps ~ 2min to run (longer for the model with predictors X & Z...)
# 6 chains 6000 reps in X & Z model:: 10:02am - 10:12 -- 10 min
# 6 chains 4000 reps in X & Z model ~ 10:19 - 10:28 am -- 9 min...
# back to 6000 reps, 6 chains: 10:30 - XXX -- XX min

# checking fit
# bc
skipTrack.diagnostics(fitbcfr_bc, param = 'cijs') # trace plot looks... just very very messy??? I think that's just high uncertainty in general...
skipTrack.diagnostics(fitbcfr_bc, param = 'Betas') # these chains have not converged to a single answer... still looks bad. I think this is about information in the data, not the actual MCMC fitting process. -- i.e. they are just not that good of predictors for the mean cycle length? which sort of makes sense...
skipTrack.diagnostics(fitbcfr_bc, param = 'Gammas') # I am not seeing gains after 4000, for sure, so we will cut off there... but then this new round with 4000 looks terribly unconverged???
skipTrack.results(fitbcfr_bc, burnIn = 3000) # god freaking dammit LOL  -- we are okay except for the betas... which kind of makes sense to me? I expect these to predict irregularity better than mean...

# methodological note for future: potential issues with G-R statistc & proposed higher threshold than 1.1
# https://arxiv.org/pdf/1812.09384

# let's compare how different the predicted c_ijs are! between the two models
# note that they were run on the same data, even though there might have been a COUPLE observations excluded explicitly because they did not have bc values/data


# basic model ----------------------------------------------------
# data
bcfr_skptr2 = bcfr_skptr
bcfr_skptr2 = bcfr_skptr2 |> filter(cycle_len < 365)
bcfr_skptr2 = bcfr_skptr2 |> group_by(PersonID) |> filter(n() >= 2) #only eliminates a few ppl! nice : 1028 obs

# model
set.seed(9819)
fitbcfr_basic = skipTrack.fit(Y = bcfr_skptr2$cycle_len,
                              cluster = bcfr_skptr2$PersonID,
                              chains = 6, # increasing # chains to help with convergence
                              reps = 3000) # increasing # reps to help with convergence

# checking fit
# basic
skipTrack.diagnostics(fitbcfr_basic, param = 'cijs') # trace plot looks SLIGHTLY better!!!!
skipTrack.results(fitbcfr_basic, burnIn = 750) # "recommended burn-in with real data is at least 750" -- YES & the G-R looks way better!!

# --------------

# combine

# add c_ij results to df
skip_results = bcfr_skptr1
skip_results = cbind(skip_results, skipTrack.results(fitbcfr_basic, burnIn = 750)$cijs)
colnames(skip_results)[colnames(skip_results) == c("Lower", "Mean", "Upper")] <- c("cijbasic_lowerCI", "cijbasic_mean", "cijbasic_upperCI")
skip_results = cbind(skip_results, skipTrack.results(fitbcfr_bc, burnIn = 3000)$cijs)
colnames(skip_results)[colnames(skip_results) == c("Lower", "Mean", "Upper")] <- c("cijbc_lowerCI", "cijbc_mean", "cijbc_upperCI")
# the mean & 95% CI estimates are VERY similar

# --------------

# other variables to the df

# 1. variance
bcfr_sk_cij = bcfr_sk_cij |>  group_by(PersonID) |> mutate(cycle_var = var(cycle_len),
                                                           cycle_sd = sqrt(var(cycle_len)))
# a LOT of this *must* be reporting issues, or else EVERYONE has pretty crazy variation in cycle lengths...


# can see how these align for the 57 ppl who share data:
irr_id = bcfr_irr[c("PersonID","IRREGULAR_MENSES")]
merged = merge(irr_id, bcfr_sk_cij, by = "PersonID")

# recode irregularity as 0 & 1
# 1 = irregular
# 0 = regular
merged = merged |>
  mutate(
    irregular = case_when(
      IRREGULAR_MENSES == 2 ~ 0, # regular
      IRREGULAR_MENSES == 1 ~ 1 # irregular (Eva said 1 is irregular)
    )
  )


# --------------------------------

# tau_i

# !! I would like to add tau_i results to df...
# get precision by Person ID
# get tau_is for each individual.... bc it is the filtered variance minus the skips! I think that is useful...
# iterate through MCMC iterations: exclude burn-in
# lapply() is not internally faster than for-loop

# this is WRONG -- do this later
# empty_list <- vector(mode = "list", length = nrow(fitbcfr_bc$fit[[1]][[540]]$iDat))  # i.e. one entry for each individual
# names(empty_list) = fitbcfr_bc$fit[[1]][[540]]$iDat$Individual
# for (ii in seq(1,length(fitbcfr_bc$fit))){ # 6 -- 6 chains
#   for (jj in seq(3001, length(fitbcfr_bc$fit[[ii]]))){ # 6001 -- 6000 reps, but with 3000 burn-in i.e. they are excluded!
#     MCMC_taus = fitbcfr_bc$fit[[ii]][[jj]]$iDat[c("Individual","taus")]
#     for (kk in seq(1,nrow(MCMC_taus))){
#       individual_id = MCMC_taus[kk, "Individual"]
#       empty_list$individual_id = c(empty_list$individual_id, MCMC_taus[kk, "taus"])
#     }
#   } 
# }
# # combine the list of lists (one list per individual) into a df
# merge(empty_list, by = "Individual") # this will not work...

# now get posterior distribution (mean & 95% credible interval, I would say)

# note: I could calculate my own G-R for tau_i.... ***** could help decide whether these estimates are valid or not!! or look at other potential diagnostics...**



