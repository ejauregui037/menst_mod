
# examine correlation between self-reported irregularity and tau_i value from basic BCFR algorithm (no covariates, only observed cycles)
# & with skipTrack fittings from models *with* covariates? e.g. BC, endometriosis...






# do I need to update the README for this script?****





# NOTE: tau_i = precision = 1/variance NOT REALLY -- it's the inverse of the covariance matrix... ***







# 1. fit skipTrack with bcfr data, no covariates
# source("./skiptrack/data_clean/bcfr_basicfit.R")

# 2. fit skipTrack with Clue BC data as covariate(s -- indicator vars)
# have I consolidated this into a script? no**
# however, when i did this (not ad hoc, but i didn't set a seed.....) the c_ij estimates were ALMOST identical...

# 3. fit skipTrack with endometriosis as a covariate (leave out other bleeding/potential relevant conditions for now; I think that there are simply note enough unique BCFR YWS # individuals for the model to infer much from such a small dataset...)
# then compare? how is this really abt irregularity tho LOL


# setwd("~/Desktop/ColumbiaPhD/research_projects/Lauren_menstrualmodelingmethods/menst_mod")

# load skipTrack-processed data
source("./skiptrack/data_clean/bcfr_skiptrack.R")
library(skipTrack)

# data
bcfr_irr = read.csv("./data/yws_clue/ywsclue_repro.csv") # BCFR irregularity data (Clue) -- has "ever irregular" variable : IRREGULAR_MENSES

# bcfr_irr$MENSES_STARS_AS_EXPECTED: "1: Yes 2: No 3: Don’t know"
# During the last 12 months, did your menstrual period usually start within 4 days of the day you expected it to start? By “usually” we mean for at least half of the time.

# recode
bcfr_irrdata = bcfr_irr |>
  mutate(
    irr_ever = case_when(
      IRREGULAR_MENSES == 2 ~ 0, # regular
      IRREGULAR_MENSES == 1 ~ 1, # irregular (Eva said 1 is irregular)
      IRREGULAR_MENSES == 3 ~ 0 # 3 = "don't know"..... does that mean regular? too many people to just exclude them (29)
      ), 
    irr_12mo = case_when(
      MENSES_STARS_AS_EXPECTED == 2 ~ 0, # regular
      MENSES_STARS_AS_EXPECTED == 1 ~ 1,
      MENSES_STARS_AS_EXPECTED == 3 ~ 0 # 3 = "don't know"..... does that mean regular? too many people to just exclude them (29)
      )
  )

# (!) there are a lot of 3s ("don't know") -- I let them be "regular".... 
# could just exclude them? ask

# filter for people with these IRR variables
bcfr_irrdata = bcfr_irrdata |> filter(!is.na(IRREGULAR_MENSES) & !is.na(MENSES_STARS_AS_EXPECTED)) # dropped 47 ppl (!) (47/398)

# Check: how correlated are ever & current menses
cor(bcfr_irrdata$irr_ever, bcfr_irrdata$irr_12mo) #  -0.2902531 ?!?!? TOTALLy diff

# how does this align with variance? tau?
bcfr_yws_irrEVA = read.csv("./skiptrack/data_clean/datasets_tosend/bcfryws_irr.csv")
merged = merge(bcfr_irrdata, bcfr_yws_irrEVA, by = 'PersonID')
# check:
merged$irr_ever == merged$irregular # all TRUEs

# now check correlations:

cor(merged$irr_ever, merged$mean_tau)  # negatively associated
# vs.
cor(merged$irr_12mo, merged$mean_tau)  # positively associated!

cor(merged$irr_ever, merged$cycle_var)   # negatively associated
# vs.
cor(merged$irr_12mo, merged$cycle_var) # still negatively associated....?

# this is... strange, but interesting?







rm(list=ls())
# load skipTrack-style data
source("./skiptrack/data_clean/bcfr_skiptrack.R")
# new variables
bcfr_new = read.csv("./data/new_NOTSORTEDYET_needtoconvert_tocsv/Repro_data_new.csv") # BCFR YWS -- medical data -- Clue -- includes:
# same # obs, but new variables available! including endometriosis, ever irregular
# ENDOMETRIOSIS -- "Have you ever been diagnosed with..." -- "1: Yes 2: No 3: Don’t know"	
# ?? how to code don't know
table(bcfr_new$ENDOMETRIOSIS) # only 3 ppl -- exclude them
# recode
bcfr_new1 = bcfr_new |>
  mutate(
    endo = case_when(
      ENDOMETRIOSIS == 2 ~ 0, # regular
      ENDOMETRIOSIS == 1 ~ 1, # irregular (Eva said 1 is irregular)
      ENDOMETRIOSIS == 3 ~ NA # 3 = "don't know"..... does that mean regular? too many people to just exclude them (29)
    )
  ) |> 
  filter(!is.na(endo)) # nice, filtered 3 ppl

# ! note -- only 6 ppl in the sample had endometriosis...

# fit skiptrack -------------------------------------------------
bcfr_skptr1 = merge(bcfr_skptr, bcfr_new1, by = "PersonID") # 1019 obs... EVERYONE had endo data?? how can that be??
# uh oh.... is it that EVERYONE with endo is excluded anyway? i.e. doesn't have clue data?

# check:
setdiff(as.numeric(bcfr_new1$PersonID), as.numeric(bcfr_skptr$PersonID)) # 338 diff.... MOST people do not have clue data... including the ppl with endo :.)
length(unique(bcfr_skptr$PersonID))  # 57
length(unique(bcfr_new1$PersonID)) # 395



# nope: no one in clude has endo


# 
# # fit skiptrack with endo as a predictor of cycle length & regularity
# 
# # hyper-parameters have been chosen based on previous runs/tests; however, as visible, convergence is not always great.
# # burn-in chosen based on plots -- ~3000?
# set.seed(9999)
# fitbcfr_endo = skipTrack.fit(Y = bcfr_new1$cycle_len,
#                            cluster = bcfr_new1$PersonID,
#                            X = as.matrix(bcfr_new1$endo), # covariates for cycle length mean -- left out "none" as "ref" i.e. so 1 df (i.e. intercept/when all covars =0)
#                            Z = as.matrix(bcfr_new1$endo), # covariates for cycle length precision i.e. 1/variance
#                            chains = 6,
#                            reps = 6000) 
# 
# # crude benchmarking:
# # 6 chains 3000 reps ~ 2min to run (longer for the model with predictors X & Z...)
# # 6 chains 6000 reps in X & Z model:: 10:02am - 10:12 -- 10 min
# # 6 chains 4000 reps in X & Z model ~ 10:19 - 10:28 am -- 9 min...
# # back to 6000 reps, 6 chains: 10:30 - XXX -- XX min
# 
# # checking fit
# # bc
# skipTrack.diagnostics(fitbcfr_bc, param = 'cijs') # trace plot looks... just very very messy??? I think that's just high uncertainty in general...
# skipTrack.diagnostics(fitbcfr_bc, param = 'Betas') # these chains have not converged to a single answer... still looks bad. I think this is about information in the data, not the actual MCMC fitting process. -- i.e. they are just not that good of predictors for the mean cycle length? which sort of makes sense...
# skipTrack.diagnostics(fitbcfr_bc, param = 'Gammas') # I am not seeing gains after 4000, for sure, so we will cut off there... but then this new round with 4000 looks terribly unconverged???
# skipTrack.results(fitbcfr_bc, burnIn = 3000) # god freaking dammit LOL  -- we are okay except for the betas... which kind of makes sense to me? I expect these to predict irregularity better than mean...
# 



























# how many individuals in each sample
length(unique(bcfr_trk$PersonID)) # 61 (!) -- & actually cuts down to 57 once i restrict sample to those with period or spotting data

# all-time irregularity
length(unique(bcfr_irr$PersonID)) # 398

bcfr_irr$IRREGULAR_MENSES # coded 2 or 1 -- no idea what either means though... 


# 1. fit skiptrack w no covariates






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




