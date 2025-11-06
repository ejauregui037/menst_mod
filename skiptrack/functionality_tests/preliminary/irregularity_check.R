# examine correlation between self-reported irregularity and tau_i value from basic BCFR algorithm (no covariates, only observed cycles)

# NOTE: tau_i = precision = 1/variance ***

# note: I do not have the codebook for the BCFR data. I do not know whether 2 & 1 = regular or irregular

# setup
# root directory
setwd("~/Desktop/ColumbiaPhD/research_projects/Lauren_menstrualmodelingmethods/menst_mod")
# librarys
library(skipTrack)
library(dplyr)
library(tidyr)
# data
bcfr_trk = read.csv("./data/yws_clue/ywsclue_track.csv") # BCFR tracking data (Clue)
bcfr_irr = read.csv("./data/yws_clue/ywsclue_repro.csv") # BCFR irregularity data (Clue)

# how many individuals in each sample
length(unique(bcfr_trk$PersonID)) # 61 (!) -- & actually cuts down to 57 once i restrict sample to those with period or spotting data
length(unique(bcfr_irr$PersonID)) # 398

bcfr_irr$IRREGULAR_MENSES # coded 2 or 1 -- no idea what either means though... 

# function to get cycle lengths

cycle_len_only <- function(df){
  
  # exclude bc data
  df_skip = df[df$category == "period" | df$category == "spotting" ,] # none are missing? odd but i guess by definition these are ENTRIES in the clue app
  
  # get cycle lengths for each individual
  df_skip = df_skip |> group_by(PersonID) |> arrange(track_date_DtDemogr, .by_group = T)  |>    # Cmd Shift M for pipe
    mutate(gap_days = track_date_DtDemogr - lag(track_date_DtDemogr),
           new_cycle = if_else(gap_days > 4, 1, 0)) |>  # if difference between rows >4
    arrange(track_date_DtDemogr, .by_group = T)  
  df_skip$new_cycle[1] = 1 # doesn't happen with algo, but is implied by meaning
  
  # cycle length
  df_cycle = df_skip |> group_by(PersonID) |> arrange(track_date_DtDemogr, .by_group = T)  |> 
    filter(new_cycle == 1) |> 
    mutate(cycle_len = lead(gap_days))
  
  # for skiptrack
  df_skptr = df_cycle |> select(PersonID, cycle_len) |> 
    group_by(PersonID) |> mutate(cycle_num = row_number()) |> 
    filter(!is.na(cycle_len))
  
  # return
  return(df_skptr)
}

# get BCFR clue data formatted for skiptrack
bcfr_skptr = cycle_len_only(bcfr_trk)
# fit skiptrack, with no covariates -- & no decisions made about the optimal way to actually fit this MCMC model: burn-in, etc.
fitbcfr_basic = skipTrack.fit(Y = bcfr_skptr$cycle_len,
                              cluster = bcfr_skptr$PersonID,
                              chains = 4)
# get precision by Person ID
plot(fitbcfr_basic)
summary(fitbcfr_basic, burnIn = 500)
confint(fitbcfr_basic)
skipTrack.results(fitbcfr_basic)
# mu_is and tau_is are not reported -- unstable??
names(fitbcfr_basic)
skipTrack.diagnostics(fitbcfr_basic) # looks AWFUL lol... could be due to having soooo few data points...

# several issues:
# 1. I am not sure that this model can actually converge with the data provided -- there are very few people, overall, and relatively few sequences...
# 2. tau_i and mu_i are latent variables in the model -- it is likely that the model is not optimizing their estimates! so they might not be good to extract...
# 3. even if I DID extract tau_i, what is it really? it's a measure of variance in an invidiual's cycle lengths, with the 
# extra Bayesian information from other people's cycle lengths. But if the question is whether or not self-reported irregularity matches up with irregularity in 
# the data, then we should just estimate the variability of cycle lengths and regress that on self-reported irregularity directly
   # 4. so then... comparing tau_i & self-reported irregularity would be... for what? to see how they line up in the model? if c_ij is in the model, then wouldn't some 
   # sort of statistic about c_ij be more relevant? but why would we want that...
   # could just be mean c_ij?


# irregular = 1 is irregular!

# new plan: compare...

# get skipTrack c_ijs, to compare
cijs_byrow = skipTrack.results(fitbcfr_basic)$cijs
colnames(cijs_byrow) = c("cij_lower", "cij_mean", "cij_upper")
bcfr_sk_cij = cbind(bcfr_skptr, cijs_byrow)

# 1. variance
bcfr_sk_cij = bcfr_sk_cij |>  group_by(PersonID) |> mutate(cycle_var = var(cycle_len),
                                                         cycle_sd = sqrt(var(cycle_len)))
# a LOT of this *must* be reporting issues, or else EVERYONE has pretty crazy variation in cycle lengths...

# 2. clinical definition ("varies more than 5 days per month" (I'm going to interpret that as "on average" i.e. mean "distance" (mean distance from the mean)))
bcfr_sk_cij = bcfr_sk_cij |>  group_by(PersonID) |> mutate(clin_var = mean(abs(cycle_len - mean(cycle_len))))
# with this criterion, it seems like roughly a third of individuals have clinically "irregular" cycles (if we are using the 5-day (mean distance-ish) cutoff that Lauren mentioned)

# 3. mean c_ijs (i.e. mean "skips", which would imply high variability & c_ij is actually the target estimand for the M-H algorithm so should converge well, etc.)
bcfr_sk_cij = bcfr_sk_cij |>  group_by(PersonID) |> mutate(mean_cijs = mean(cij_mean),
                                                           mean_confint_width = mean(cij_upper - cij_lower)) # i am also interested in the confidence/precision of these c_ij estimates, cause for the bigger one it sucks... need to include covariates!
# !!! recall: these results are with the basic skipTrack, no relevant covariates are in yet... 
# probably super deflates performamnce!! need to include those! 
# This is NOT a minor thing!!!

# regressions
# assume linearity? in relationship between these vars and self-reported irregularity...

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

# 1. variance & sd with self-reported irregularity 
irr_var = lm(merged$irregular ~ merged$cycle_var)
irr_sd = lm(merged$irregular ~ merged$cycle_sd)
# 2. clinical definition i.e. mean distance (kept continuous, did not apply 5-day cutoff)
irr_clin = lm(merged$irregular ~ merged$clin_var)
# 3. mean c_ijs (**issue because of uncertainty... & LACK OF COVARS*** should SIGNIFICANTLY improve performance, I would say! predictive**)
irr_cij = lm(merged$irregular ~ merged$mean_cijs)

# R^2 variables (% variability explained)
summary(irr_var)
summary(irr_sd)
summary(irr_clin)
summary(irr_cij)

# plot
library(ggplot2)
ggplot(data = merged, aes(x = cycle_var, y = irregular)) + geom_point()
ggplot(data = merged, aes(x = cycle_sd, y = irregular)) + geom_point()
ggplot(data = merged, aes(x = clin_var, y = irregular)) + geom_point()
ggplot(data = merged, aes(x = mean_cijs, y = irregular)) + geom_point()

# need to remove ppl with only one cycle reported ! 
# & the people with otherwise weird data (i.e. one cycle reported that was over 500 days !!! she just didn't interact w the app for 1.5 years...)
  # this might improve the fits... outliers can have a BIG impact in linear regression... !esp: see the ones with CRAZY sds!!!
 
 # IMPORTANT DATA-CLEANING STEP!! DO THIS***^^

# I am also personally very interested in the REVERSE, i.e. the difference in the means of the variables according to binary irregular (i.e. the difference in means between those groups)
# 1. variance & sd with self-reported irregularity 
var_irr = lm(merged$cycle_var ~ merged$irregular)
sd_irr = lm(merged$cycle_sd ~ merged$irregular)
# 2. clinical definition i.e. mean distance (kept continuous, did not apply 5-day cutoff)
clin_irr = lm(merged$clin_var ~ merged$irregular)
# 3. mean c_ijs (**issue because of uncertainty... & LACK OF COVARS*** should SIGNIFICANTLY improve performance, I would say! predictive**)
cij_irr = lm(merged$mean_cijs ~ merged$irregular)

# R^2 variables (% variability explained)
summary(var_irr)
summary(sd_irr)
summary(clin_irr)
summary(cij_irr)
# they are indeed all positively associated.....

# "30 observations deleted due to missingness" huh? where? 

# also: these are just not... super related... none of them are GREAT predictors of each other... & the "irregular" as binary is just so -- so much information lost...

# because i'm interested:
cij_sd = lm(merged$mean_cijs ~ merged$cycle_sd)
summary(cij_sd) # R^2 of 0.8359 -- makes sense! but it's not 1:1, which is also nice... 
# that's the bayesian bit, I would presume... could diverge more once I introduce covariates into the skipTrack model to predict both mean cycle & cycle variability!




