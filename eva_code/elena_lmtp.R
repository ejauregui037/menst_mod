
# stressmens_updated code.R = where all these objects & variables came from! Eva's code

# use merged_df

# vars: MENSES_STARS_AS_EXPECTED, STRESS_CT/STRESS_TH       &     PREG_EVER + ALCOHOL_EVER + SMOKE_LAST10_EVER + AgeIntv_Demographics

library(lmtp)
# benefits of lmtp for our analysis: flexible exposure & outcome models** i.e. less misspecification
# no time-varying confounding being addressed -- only 2 timepoints
# although this does NOT fix all of our measurement issues!!!



# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# lmtp setup

# to make results more robust:
# add learners?
# add CV folds?

# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------

# ENSURE THAT THESE ARE ML LEARNERS FOR MAXIMAL FLEXIBILITY -- 
# add additional?

learners_trtmodel = c(
  "SL.glm",       # logistic regression
  "SL.ranger",    # rf (ml)
  "SL.glmnet"     # lasso/elastic net (ml)
)

learners_outcomemodel = c(
  "SL.glm",
  "SL.ranger",
  "SL.glmnet"
)

# alllll varnames to lmtp -- what is this for? not sure i needed it lol 
vars_check_count = c("PREG_EVER", "ALCOHOL_EVER", "SMOKE_LAST10_EVER", "AgeIntv_Demographics",
                     "MENSES_STARS_AS_EXPECTED", "STRESS_CT")
vars_check_severity = c("PREG_EVER", "ALCOHOL_EVER", "SMOKE_LAST10_EVER", "AgeIntv_Demographics",
                        "MENSES_STARS_AS_EXPECTED", "STRESS_TH")

# varnames to lmtp
treatment_count = c("STRESS_CT") 
treatment_severity = c("STRESS_TH") 

# basline, fixed confounding vars
baseline_confounders = c("PREG_EVER", "ALCOHOL_EVER", "SMOKE_LAST10_EVER", "AgeIntv_Demographics")

# no time-varying confounding, only two timepoints here
# time_vary_confounders = list(
#   c("r2cesd_BIN"),   # confounders to UE2 (post-baseline, pre-A1)
#   c("r3cesd_BIN")    # confounders to UE3 (post-A1, pre-A2)
# )

# recent irregularity = outcome (binary)
outcome_recent = "MENSES_STARS_AS_EXPECTED"

# DATA
# needs to be complete case... unless we want to handle censroing?? but we'd need to think of how to predict that & time-varying confounding as a result
merged_df_cc_count = merged_df |> drop_na(any_of(vars_check_count)) # can do "all of" or "any of" to drop_NA! COOL
merged_df_cc_severity = merged_df |> drop_na(any_of(vars_check_severity)) # can do "all of" or "any of" to drop_NA! COOL


# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# baseline -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------

# https://www.beyondtheate.com/04_lmtp.html#:~:text=The%20ipw%20and%20gfoRmula%20packages,conducting%20causal%20analyses%20in%20R!

# estimate Y(as observed) [no shift]

shift_null = function(data,trt){
  data[[trt]]  # return treatment values as they are; could also pass NULL to shift fxn in lmtp
}


# run lmtp


# (!) check whether this # of folds is stable **** @ gridsearch ****

set.seed = 90890
exp_Y0 = lmtp_tmle(
  data         = merged_df_cc_count,   
  trt          = treatment_count,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_null,                # static intervention: always set A=1
  outcome_type = "binomial",        # Y is BINARY -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
  
  # poisson might be better in finite samples.... **** TRY IT BOTH WAYS***
  # outcome_type = "binomial",  -- but set to poisson!
  
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5     # default is 10; HOWEVER i got essentially the same results with 5 (& very stable across runs), so i will stick to 5 (10 takes much much longer to run, with no benefit)
  
)

# print results
print(exp_Y0)
summary(exp_Y0) # check what these are telling me

exp_Y00 = lmtp_tmle(
  data         = merged_df_cc_severity,   
  trt          = treatment_severity,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_null,                # static intervention: always set A=1
  outcome_type = "binomial",        # Y is BINARY -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
  
  # poisson might be better in finite samples.... **** TRY IT BOTH WAYS***
  # outcome_type = "binomial",  -- but set to poisson!
  
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5     # default is 10; HOWEVER i got essentially the same results with 5 (& very stable across runs), so i will stick to 5 (10 takes much much longer to run, with no benefit)
  
)

print(exp_Y00)
summary(exp_Y00) # check what these are telling me


# these two should be the SAME^^^^ this is a check!
# they do! nice : )




# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# lmtp -- shift -- stress counts -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------

# 1 less stressful event than observed
shift_count1 = function(data,trt){
  a = data[[trt]]
  ifelse(a >= 1, a - 1, 0) # i.e. zero is a lower bound -- do not shift ppl lower than zero
}

# 5 fewer stressful events than observed
shift_count5 = function(data,trt){
  a = data[[trt]]
  ifelse(a >= 5, a - 5, 0) # i.e. zero is a lower bound -- do not shift ppl lower than zero
}


exp_Y1less = lmtp_tmle(
  data         = merged_df_cc_count,   
  trt          = treatment_count,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_count1,                # static intervention: always set A=1
  outcome_type = "binomial",        # Y is BINARY -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
  # three options: binomial, continuous, survival (which i assume reqs censoring?)
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5     # default is 10; HOWEVER i got essentially the same results with 5 (& very stable across runs), so i will stick to 5 (10 takes much much longer to run, with no benefit)
  
)

exp_Y5less = lmtp_tmle(
  data         = merged_df_cc_count,   
  trt          = treatment_count,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_count5,                # static intervention: always set A=1
  outcome_type = "binomial",        # Y is BINARY -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
  
  # poisson might be better in finite samples.... **** TRY IT BOTH WAYS***
  # outcome_type = "binomial",  -- but set to poisson!
  
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5     # default is 10; HOWEVER i got essentially the same results with 5 (& very stable across runs), so i will stick to 5 (10 takes much much longer to run, with no benefit)
  
)

print(exp_Y1less)
print(exp_Y5less) # oh niceeee that's an effect!!! : D

# calculate RD (or RR)
lmtp_contrast(exp_Y1less, ref = exp_Y0, type = "additive") # use Y00 bc that was on the same sample i.e. same n
lmtp_contrast(exp_Y5less, ref = exp_Y0, type = "additive") # use Y00 bc that was on the same sample i.e. same n
# !!!
# !!!


# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# lmtp -- shift -- stress severity -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------

# 10% less stress severity than observed
shift_severity10 = function(data,trt){
  a = data[[trt]]
  a*0.90 # respects zero bound automatically -- reduce by 10% !
}

shift_severity25 = function(data,trt){
  a = data[[trt]]
  a*0.75 # respects zero bound automatically -- reduce by 10% !
}


exp_Y10perless = lmtp_tmle(
  data         = merged_df_cc_severity,   
  trt          = treatment_severity,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_severity10,                # static intervention: always set A=1
  outcome_type = "continuous",        # Y is CONTINUOUS -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5     # default is 10; HOWEVER i got essentially the same results with 5 (& very stable across runs), so i will stick to 5 (10 takes much much longer to run, with no benefit)
  
  # ! ! !! ! check folds for stability & try other folds!
  
)


exp_Y25perless = lmtp_tmle(
  data         = merged_df_cc_severity,   
  trt          = treatment_severity,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_severity25,                # static intervention: always set A=1
  outcome_type = "continuous",        # Y is CONTINUOUS -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5     # default is 10; HOWEVER i got essentially the same results with 5 (& very stable across runs), so i will stick to 5 (10 takes much much longer to run, with no benefit)
  
  # ! ! !! ! check folds for stability & try other folds!
  
)


print(exp_Y10perless) # uhhh is this a reduction LOL -- not sure about that!
print(exp_Y25perless) # this IS a reduction... wtf HAHA is it nonlinear?! or are we just sooooo underpowered....


# calculate RD (or RR)
lmtp_contrast(exp_Y10perless, ref = exp_Y00, type = "additive") # use Y00 bc that was on the same sample i.e. same n
lmtp_contrast(exp_Y25perless, ref = exp_Y00, type = "additive") # use Y00 bc that was on the same sample i.e. same n
# !!!



# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# check positivity !!! ! ! ! -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------



# # # # #
# see epi 7 ass 2 + slides & syllabi







# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# old eva code @ model-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------


# Stress -- counts
# irregularity -- recent : MENSES_STARS_AS_EXPECTED (why is it not 12mo?)

# mod1_tbl <- glm(MENSES_STARS_AS_EXPECTED ~ STRESS_CT 
#                 + PREG_EVER + ALCOHOL_EVER + SMOKE_LAST10_EVER + AgeIntv_Demographics
#                 , 
#                 family = poisson(link = "log"), 
#                 data = merged_df) %>%
#   tbl_regression(
#     exponentiate = TRUE, 
#     label = list(
#       STRESS_CT ~ "# of Stressful Events"
#       , 
#       PREG_EVER ~ "Ever Pregnant",
#       ALCOHOL_EVER ~ "Ever Alcohol Use",
#       SMOKE_LAST10_EVER ~ "Ever Smoked (Last 10 Years)",
#       AgeIntv_Demographics ~ "Age"
#     ),
#     tidy_fun = partial(tidy_robust, vcov = "HC1")
#   ) %>%
#   add_n(location = "label")

# Stress -- severity
# irregularity -- recent : MENSES_STARS_AS_EXPECTED (why is it not 12mo?)

# mod2_tbl <- glm(MENSES_STARS_AS_EXPECTED ~ STRESS_TH 
#                 + PREG_EVER + ALCOHOL_EVER + SMOKE_LAST10_EVER + AgeIntv_Demographics
#                 ,
#                 family = poisson(link = "log"), 
#                 data = merged_df) %>%
#   tbl_regression(
#     exponentiate = TRUE, 
#     label = list(
#       STRESS_TH ~ "Severity of Stress"
#       , 
#       PREG_EVER ~ "Ever Pregnant",
#       ALCOHOL_EVER ~ "Ever Alcohol Use",
#       SMOKE_LAST10_EVER ~ "Ever Smoked (Last 10 Years)",
#       AgeIntv_Demographics ~ "Age"
#     ),
#     tidy_fun = partial(tidy_robust, vcov = "HC1")
#   ) %>%
#   add_n(location = "label")
# 
# tbl_merged <- tbl_merge(
#   tbls = list(mod1_tbl, mod2_tbl),
#   tab_spanner = c("**# of Stressful Events**", "**Severity of Stress**")
# )