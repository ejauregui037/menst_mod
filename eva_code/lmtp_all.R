# replicate all analyses in paper with lmtp

# first: source stressmens_updated code.R -- code copied below
# I've commented out the models, just so we get the dfs modified as necessary

# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# lmtp setup

# TECHNICAL NOTES

# to make results more robust:

# add learners?
# add CV folds?
  # (!) check whether this # of folds is stable **** @ gridsearch ****
# account for missingness?
  # NOTE: not dealing with missing variables... we don't have CENSORING per se, since it was all collected at once.
  # check whether missingness is random
# Check PS scores etc. for positivity violations



# My conclusion: maybe we shift away from calling these associations causal... we really do NOT have a methodological leg to 
# stand on @ saying we achieved covariate balance i.e. exchangeability


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

# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# table 3 - 1-unit changes in #stressful events, stress severity
# exposures are continuous: code a shift of 1 unit up, causal RR
# -----------------------------------------------------------------------------------------------------------------------------------------------------------

# DATA
# needs to be complete case... unless we want to handle censroing?? but we'd need to think of how to predict that & time-varying confounding as a result
merged_df_cc_count = merged_df |> drop_na(any_of(vars_check_count)) # can do "all of" or "any of" to drop_NA! COOL
merged_df_cc_severity = merged_df |> drop_na(any_of(vars_check_severity)) # can do "all of" or "any of" to drop_NA! COOL

shift_null = function(data,trt){
  data[[trt]]  # return treatment values as they are; could also pass NULL to shift fxn in lmtp
}

set.seed = 90890
ref_count = lmtp_tmle( # no intervention
  data         = merged_df_cc_count,   # count
  trt          = treatment_count,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_null,         # no intervention, i.e. observed values
  outcome_type = "binomial",        # Y is BINARY -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5     
)
ref_sev = lmtp_tmle( # no intervention
  data         = merged_df_cc_severity,   # severity
  trt          = treatment_severity,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_null,               # no intervention, i.e. observed/expected values: diff is the df -- no missingness, with severity var list
  outcome_type = "binomial",        # Y is BINARY -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5    
  
)

# 1 more stressful event than observed -- to correspond to regression results
shift_count1 = function(data,trt){
  a = data[[trt]]
  a + 1
}

up_count1 = lmtp_tmle(
  data         = merged_df_cc_count,   
  trt          = treatment_count,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_count1,     
  outcome_type = "binomial",        # Y is BINARY -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5     # from epi7:  default is 10; HOWEVER i got essentially the same results with 5 (& very stable across runs), so i will stick to 5 (10 takes much much longer to run, with no benefit)
)

# 1 more severity point than observed -- to match up with regression results
shift_severity1 = function(data,trt){
  a = data[[trt]]
  a + 1
}

up_sev1 = lmtp_tmle(
  data         = merged_df_cc_severity,   
  trt          = treatment_severity,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_severity1,          
  outcome_type = "binomial",        
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5     
)

# calculate RRs
table3_count = lmtp_contrast(up_count1, ref = ref_count, type = "rr")
table3_sev = lmtp_contrast(up_sev1, ref = ref_sev, type = "rr") # use Y00 bc that was on the same sample i.e. same n

# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# table 4 -- 5 categories with different exposure timings, one is reference group ()
# categorical: compare setting everyone to ref level vs. to the other 4 levels

# (!) there are definitely positivity violations here... XXX

# also: this ordering is different than in the tables in the paper.... (!)
# -----------------------------------------------------------------------------------------------------------------------------------------------------------

# DATA
vars_check_window = c("PREG_EVER", "ALCOHOL_EVER", "SMOKE_LAST10_EVER", "AgeIntv_Demographics",
                     "MENSES_STARS_AS_EXPECTED", "StressWindow4Cat")
# varnames to lmtp
treatment_window = c("StressWindow4Cat") 
# needs to be complete case... unless we want to handle censroing?? but we'd need to think of how to predict that & time-varying confounding as a result
merged_df_cc_window = merged_df |> drop_na(any_of(vars_check_window)) # can do "all of" or "any of" to drop_NA! COOL
# significantly more ppl lacking these data...

# levels: for my reference
    # mutate(StressWindow4Cat = case_when(
    #   High_Early == 0 & High_Adult == 0 & High_Recent == 0 ~ 0, -- reference level i.e. low in all periods!
    #   High_Early == 1 & High_Adult == 0 & High_Recent == 0 ~ 1,
    #   High_Early == 0 & High_Adult == 1 & High_Recent == 0 ~ 2,
    #   High_Early == 0 & High_Adult == 0 & High_Recent == 1 ~ 3,
    #   High_Early == 1 & High_Adult == 1 & High_Recent == 1 ~ 4,
    #   TRUE ~ NA_real_
    # ))

shift_windowref = function(data,trt){   # low in all periods = reference
  a = data[[trt]]
  factor(rep("0", nrow(data)), levels = levels(a)) # bc it's a factor var
}
shift_window1 = function(data,trt){   # low in all periods = reference
  a = data[[trt]]
  factor(rep("1", nrow(data)), levels = levels(a)) # bc it's a factor var
}
shift_window2 = function(data,trt){   # low in all periods = reference
  a = data[[trt]]
  factor(rep("2", nrow(data)), levels = levels(a)) # bc it's a factor var
}
shift_window3 = function(data,trt){   # low in all periods = reference
  a = data[[trt]]
  factor(rep("3", nrow(data)), levels = levels(a)) # bc it's a factor var
}
shift_window4 = function(data,trt){   # low in all periods = reference
  a = data[[trt]]
  factor(rep("4", nrow(data)), levels = levels(a)) # bc it's a factor var
}


set.seed = 90890
ref_window = lmtp_tmle( 
  data         = merged_df_cc_window,   # severity
  trt          = treatment_window,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_windowref,               # no intervention, i.e. observed/expected values: diff is the df -- no missingness, with severity var list
  outcome_type = "binomial",        # Y is BINARY -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5    
)

ref_window1 = lmtp_tmle( 
  data         = merged_df_cc_window,   # severity
  trt          = treatment_window,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_window1,               # no intervention, i.e. observed/expected values: diff is the df -- no missingness, with severity var list
  outcome_type = "binomial",        # Y is BINARY -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5    
)

ref_window2 = lmtp_tmle( 
  data         = merged_df_cc_window,   # severity
  trt          = treatment_window,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_window2,               # no intervention, i.e. observed/expected values: diff is the df -- no missingness, with severity var list
  outcome_type = "binomial",        # Y is BINARY -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5    
)

ref_window3 = lmtp_tmle( 
  data         = merged_df_cc_window,   # severity
  trt          = treatment_window,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_window3,               # no intervention, i.e. observed/expected values: diff is the df -- no missingness, with severity var list
  outcome_type = "binomial",        # Y is BINARY -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5    
)

ref_window4 = lmtp_tmle( 
  data         = merged_df_cc_window,   # severity
  trt          = treatment_window,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_window4,               # no intervention, i.e. observed/expected values: diff is the df -- no missingness, with severity var list
  outcome_type = "binomial",        # Y is BINARY -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5    
)

# calculate RRs -- ordered as in table
table4_4 = lmtp_contrast(ref_window4, ref = ref_window, type = "rr")
table4_3 = lmtp_contrast(ref_window3, ref = ref_window, type = "rr")
table4_2 = lmtp_contrast(ref_window2, ref = ref_window, type = "rr")
table4_1 = lmtp_contrast(ref_window1, ref = ref_window, type = "rr")
# ref = low in all periods

# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# table 5 -- continuous & dichotomized -- counts + severity -- stratified by menarche timing (early/late)

# is there a reason we don't just check with interaction term?

# some positivity issues... (!)
# lmtp warning message: " one multinomial or binomial class has fewer than 8  observations; dangerous ground"
# -----------------------------------------------------------------------------------------------------------------------------------------------------------

# DATA
merged_df_cc_early_count = df_men_early |> drop_na(any_of(vars_check_count))
merged_df_cc_early_sev = df_men_early |> drop_na(any_of(vars_check_severity)) 
merged_df_cc_late_count = df_men_late |> drop_na(any_of(vars_check_count)) 
merged_df_cc_late_sev = df_men_late |> drop_na(any_of(vars_check_severity)) 
# superrr low n... bc we've split them in half, ofc

# references

# early
set.seed = 90890
ref_count_early = lmtp_tmle( # no intervention
  data         = merged_df_cc_early_count,   # count
  trt          = treatment_count,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_null,         # no intervention, i.e. observed values
  outcome_type = "binomial",        # Y is BINARY -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5     
)
ref_sev_early = lmtp_tmle( # no intervention
  data         = merged_df_cc_early_sev,   # severity
  trt          = treatment_severity,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_null,               # no intervention, i.e. observed/expected values: diff is the df -- no missingness, with severity var list
  outcome_type = "binomial",        # Y is BINARY -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5    
)

# late
ref_count_late = lmtp_tmle( # no intervention
  data         = merged_df_cc_late_count,   # count
  trt          = treatment_count,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_null,         # no intervention, i.e. observed values
  outcome_type = "binomial",        # Y is BINARY -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5     
)
ref_sev_late = lmtp_tmle( # no intervention
  data         = merged_df_cc_late_sev,   # severity
  trt          = treatment_severity,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_null,               # no intervention, i.e. observed/expected values: diff is the df -- no missingness, with severity var list
  outcome_type = "binomial",        # Y is BINARY -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5    
)

# early
up_count1_early = lmtp_tmle(
  data         = merged_df_cc_early_count,   
  trt          = treatment_count,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_count1,     
  outcome_type = "binomial",        # Y is BINARY -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5     # from epi7:  default is 10; HOWEVER i got essentially the same results with 5 (& very stable across runs), so i will stick to 5 (10 takes much much longer to run, with no benefit)
)
up_sev1_early = lmtp_tmle(
  data         = merged_df_cc_early_sev,   
  trt          = treatment_severity,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_severity1,          
  outcome_type = "binomial",        
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5     
)
# late
up_count1_late = lmtp_tmle(
  data         = merged_df_cc_late_count,   
  trt          = treatment_count,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_count1,     
  outcome_type = "binomial",        # Y is BINARY -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5     # from epi7:  default is 10; HOWEVER i got essentially the same results with 5 (& very stable across runs), so i will stick to 5 (10 takes much much longer to run, with no benefit)
)
up_sev1_late = lmtp_tmle(
  data         = merged_df_cc_late_sev,   
  trt          = treatment_severity,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_severity1,          
  outcome_type = "binomial",        
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5     
)

# calculate RRs
# early
table6_count_early = lmtp_contrast(up_count1_early, ref = ref_count_early, type = "rr")
table6_sev_early = lmtp_contrast(up_sev1_early, ref = ref_sev_early, type = "rr") 
# late
table6_count_late = lmtp_contrast(up_count1_late, ref = ref_count_late, type = "rr")
table6_sev_late = lmtp_contrast(up_sev1_late, ref = ref_sev_late, type = "rr") 

# (!) these looks CRAZY... but there are positivity violations so maybe that makes sense?!

# binary
treatment_count_bin = c("STRESS_CT_CAT") 
treatment_severity_bin = c("STRESS_TH_CAT") 
# new lists
vars_check_count1 = c("PREG_EVER", "ALCOHOL_EVER", "SMOKE_LAST10_EVER", "AgeIntv_Demographics",
                      "MENSES_STARS_AS_EXPECTED", "STRESS_CT_CAT")
vars_check_severity1 = c("PREG_EVER", "ALCOHOL_EVER", "SMOKE_LAST10_EVER", "AgeIntv_Demographics",
                         "MENSES_STARS_AS_EXPECTED", "STRESS_TH_CAT")

# DATA
merged_df_cc_early_count1 = df_men_early |> drop_na(any_of(vars_check_count1))
merged_df_cc_early_sev1 = df_men_early |> drop_na(any_of(vars_check_severity1)) 
merged_df_cc_late_count1 = df_men_late |> drop_na(any_of(vars_check_count1)) 
merged_df_cc_late_sev1 = df_men_late |> drop_na(any_of(vars_check_severity1)) 

# new shifts
shift_null_bin_count = function(data,trt){   # low in all periods = reference
  a = data[[trt]]
  factor(rep("Low Stress", nrow(data)), levels = levels(a)) # bc it's a factor var
}
shift_1_bin_count = function(data,trt){   # low in all periods = reference
  a = data[[trt]]
  factor(rep("High Stress", nrow(data)), levels = levels(a)) # bc it's a factor var
}
shift_null_bin_sev = function(data,trt){   # low in all periods = reference
  a = data[[trt]]
  factor(rep("Low Severity", nrow(data)), levels = levels(a)) # bc it's a factor var
}
shift_1_bin_sev = function(data,trt){   # low in all periods = reference
  a = data[[trt]]
  factor(rep("High Severity", nrow(data)), levels = levels(a)) # bc it's a factor var
}

# early
set.seed = 90890
ref_count_early_bin = lmtp_tmle( # no intervention
  data         = merged_df_cc_early_count1,   # count
  trt          = treatment_count_bin,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_null_bin_count,         # no intervention, i.e. observed values
  outcome_type = "binomial",        # Y is BINARY -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5     
)
ref_sev_early_bin = lmtp_tmle( # no intervention
  data         = merged_df_cc_early_sev1,   # severity
  trt          = treatment_severity_bin,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_null_bin_sev,               # no intervention, i.e. observed/expected values: diff is the df -- no missingness, with severity var list
  outcome_type = "binomial",        # Y is BINARY -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5    
)

# late
ref_count_late_bin = lmtp_tmle( # no intervention
  data         = merged_df_cc_late_count1,   # count
  trt          = treatment_count_bin,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_null_bin_count,         # no intervention, i.e. observed values
  outcome_type = "binomial",        # Y is BINARY -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5     
)
ref_sev_late_bin = lmtp_tmle( # no intervention
  data         = merged_df_cc_late_sev1,   # severity
  trt          = treatment_severity_bin,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_null_bin_sev,               # no intervention, i.e. observed/expected values: diff is the df -- no missingness, with severity var list
  outcome_type = "binomial",        # Y is BINARY -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5    
)

# early
up_count1_early_bin = lmtp_tmle(
  data         = merged_df_cc_early_count1,   
  trt          = treatment_count_bin,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_1_bin_count,     
  outcome_type = "binomial",        # Y is BINARY -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5     # from epi7:  default is 10; HOWEVER i got essentially the same results with 5 (& very stable across runs), so i will stick to 5 (10 takes much much longer to run, with no benefit)
)
up_sev1_early_bin = lmtp_tmle(
  data         = merged_df_cc_early_sev1,   
  trt          = treatment_severity_bin,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_1_bin_sev,          
  outcome_type = "binomial",        
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5     
)
# late
up_count1_late_bin = lmtp_tmle(
  data         = merged_df_cc_late_count1,   
  trt          = treatment_count_bin,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_1_bin_count,     
  outcome_type = "binomial",        # Y is BINARY -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5     # from epi7:  default is 10; HOWEVER i got essentially the same results with 5 (& very stable across runs), so i will stick to 5 (10 takes much much longer to run, with no benefit)
)
up_sev1_late_bin = lmtp_tmle(
  data         = merged_df_cc_late_sev1,   
  trt          = treatment_severity_bin,   
  outcome      = outcome_recent,  
  baseline     = baseline_confounders,    
  # time_vary    = time_vary_confounders,      # none 
  shift        = shift_1_bin_sev,          
  outcome_type = "binomial",        
  learners_trt = learners_trtmodel,      # for PS model
  learners_outcome = learners_outcomemodel, # for outcome model 
  folds        = 5     
)

# calculate RRs
# effects are stronger in ppl with early menarche...but we are SO underpowered it's crazy
# early
table6_count_early_bin = lmtp_contrast(up_count1_early_bin, ref = ref_count_early_bin, type = "rr")
table6_sev_early_bin = lmtp_contrast(up_sev1_early_bin, ref = ref_sev_early_bin, type = "rr") 
# late
table6_count_late_bin = lmtp_contrast(up_count1_late_bin, ref = ref_count_late_bin, type = "rr")
table6_sev_late_bin = lmtp_contrast(up_sev1_late_bin, ref = ref_sev_late_bin, type = "rr") 

# WE ARE SOOOO CRAZY UNDERPOWERED ?!?!?

# table(merged_df_cc_late_sev1$STRESS_TH_CAT, merged_df_cc_late_sev1$MENSES_STARS_AS_EXPECTED)
# 
# 0  1
# Low Severity  66  4
# High Severity 64  7

# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# table 6 -- Clue subsample is too small (n=56), so I'll just do it in the YWS sample

# positivity violations are not raising ERRORs in these analyses -- although I haven't actually looked at the data to see how bad it may be...
# -----------------------------------------------------------------------------------------------------------------------------------------------------------

# DATA

# merged_df <- merged_df %>%
#   mutate(MENSES_DAYS = factor(MENSES_DAYS) %>% 
#            relevel(ref = "26-32 days")) 

# collapse to binary to use with lmtp: category X vs NOT category X, since multinomial not supported natively
for (lev in levels(merged_df$MENSES_DAYS)) {
  merged_df[[paste0("Y_", lev)]] <- as.integer(merged_df$MENSES_DAYS == lev)
}

# DATA
# alllll varnames to lmtp -- what is this for? not sure i needed it lol 
vars_check_length_count = c("PREG_EVER", "ALCOHOL_EVER", "SMOKE_LAST10_EVER", "AgeIntv_Demographics",
                     "Y_26-32 days", "Y_<26 days", "Y_>32 days", "Y_Too variable to say", "STRESS_CT")
vars_check_length_sev = c("PREG_EVER", "ALCOHOL_EVER", "SMOKE_LAST10_EVER", "AgeIntv_Demographics",
                      "Y_26-32 days", "Y_<26 days", "Y_>32 days", "Y_Too variable to say", "STRESS_TH")
# needs to be complete case... unless we want to handle censroing?? but we'd need to think of how to predict that & time-varying confounding as a result
merged_df_cc_count_length = merged_df |> drop_na(any_of(vars_check_length_count)) # can do "all of" or "any of" to drop_NA! COOL
merged_df_cc_severity_length = merged_df |> drop_na(any_of(vars_check_length_sev)) # can do "all of" or "any of" to drop_NA! COOL
# n essentially preserved : )

# refs
length_reflmtp_count = function(length){
  ref_count = lmtp_tmle( # no intervention
    data         = merged_df_cc_count_length,   # count
    trt          = treatment_count,   
    outcome      = length,  
    baseline     = baseline_confounders,    
    # time_vary    = time_vary_confounders,      # none 
    shift        = shift_null,         # no intervention, i.e. observed values
    outcome_type = "binomial",        # Y is BINARY -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
    learners_trt = learners_trtmodel,      # for PS model
    learners_outcome = learners_outcomemodel, # for outcome model 
    folds        = 5     
  )
}
length1ref_count = length_reflmtp_count("Y_<26 days")
length2ref_count = length_reflmtp_count("Y_26-32 days")
length3ref_count = length_reflmtp_count("Y_>32 days")
length4ref_count = length_reflmtp_count("Y_Too variable to say")

length_reflmtp_sev = function(length){
  ref_count = lmtp_tmle( # no intervention
    data         = merged_df_cc_severity_length,   # count
    trt          = treatment_severity,   
    outcome      = length,  
    baseline     = baseline_confounders,    
    # time_vary    = time_vary_confounders,      # none 
    shift        = shift_null,         # no intervention, i.e. observed values
    outcome_type = "binomial",        # Y is BINARY -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
    learners_trt = learners_trtmodel,      # for PS model
    learners_outcome = learners_outcomemodel, # for outcome model 
    folds        = 5     
  )
}
length1ref_sev = length_reflmtp_sev("Y_<26 days")
length2ref_sev = length_reflmtp_sev("Y_26-32 days")
length3ref_sev = length_reflmtp_sev("Y_>32 days")
length4ref_sev = length_reflmtp_sev("Y_Too variable to say")

# up1
length_up1lmtp_count = function(length){
  up1_count = lmtp_tmle( # no intervention
    data         = merged_df_cc_count_length,   # count
    trt          = treatment_count,   
    outcome      = length,  
    baseline     = baseline_confounders,    
    # time_vary    = time_vary_confounders,      # none 
    shift        = shift_count1,         # no intervention, i.e. observed values
    outcome_type = "binomial",        # Y is BINARY -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
    learners_trt = learners_trtmodel,      # for PS model
    learners_outcome = learners_outcomemodel, # for outcome model 
    folds        = 5     
  )
}
length1up1_count = length_up1lmtp_count("Y_<26 days")
length2up1_count = length_up1lmtp_count("Y_26-32 days")
length3up1_count = length_up1lmtp_count("Y_>32 days")
length4up1_count = length_up1lmtp_count("Y_Too variable to say")

length_up1lmtp_sev = function(length){
  up1_count = lmtp_tmle( # no intervention
    data         = merged_df_cc_severity_length,   # count
    trt          = treatment_severity,   
    outcome      = length,  
    baseline     = baseline_confounders,    
    # time_vary    = time_vary_confounders,      # none 
    shift        = shift_severity1,         # no intervention, i.e. observed values
    outcome_type = "binomial",        # Y is BINARY -- but also, doesn't matter for lmtp, since these are NUISANCE models, not the target E[Y] -- i.e. should be the same if i say "binomial" vs some sort of Poission log link
    learners_trt = learners_trtmodel,      # for PS model
    learners_outcome = learners_outcomemodel, # for outcome model 
    folds        = 5     
  )
}
length1up1_sev = length_up1lmtp_sev("Y_<26 days")
length2up1_sev = length_up1lmtp_sev("Y_26-32 days")
length3up1_sev = length_up1lmtp_sev("Y_>32 days")
length4up1_sev = length_up1lmtp_sev("Y_Too variable to say")

# calculate RRs
# count
table6_count1 = lmtp_contrast(length1up1_count, ref = length1ref_count, type = "rr")
table6_count2 = lmtp_contrast(length2up1_count, ref = length2ref_count, type = "rr")
table6_count3 = lmtp_contrast(length3up1_count, ref = length3ref_count, type = "rr")
table6_count4 = lmtp_contrast(length4up1_count, ref = length4ref_count, type = "rr")
# severity
table6_sev1 = lmtp_contrast(length1up1_sev, ref = length1ref_sev, type = "rr")
table6_sev2 = lmtp_contrast(length2up1_sev, ref = length2ref_sev, type = "rr")
table6_sev3 = lmtp_contrast(length3up1_sev, ref = length3ref_sev, type = "rr")
table6_sev4 = lmtp_contrast(length4up1_sev, ref = length4ref_sev, type = "rr")

# these are not directly comparable to what are currently in our tables... lmtp does not natively support a reference level like lmtp

# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# table 7 is deeply undepowered/small sample... X -- using skipTrack measures of variance as outcomes...
# -----------------------------------------------------------------------------------------------------------------------------------------------------------


