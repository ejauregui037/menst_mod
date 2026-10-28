# fit skipTrack with bcfr data, covariates -- birth control data (by group)

setwd("~/Desktop/ColumbiaPhD/research_projects/Lauren_menstrualmodelingmethods/menst_mod")
library(skipTrack)

# load skipTrack-style data
source("./skiptrack/data_clean/bcfr_skiptrack.R")

# data prep -------------------------------------------------

# load, merge, & format BC data
bcfr_bc = read.csv("./data/yws_clue/ywsclue_bc.csv") # BCFR YWS birth control data(Clue)
bcfr_bc$bc_type = as.factor(bcfr_bc$BC_TYPE)

# merged = merge(bcfr_skptr, bcfr_bc, by = "PersonID") # ! issue -- some individuals gave multiple values here, but not dates for when those values differed?? -- except "BC_DATE_DtDemogr" ??
# solution: supposedly choose the most recent value of bc_type?
bcfr_bc1 = bcfr_bc |> group_by(PersonID) |> 
  mutate(bc = bc_type[which.max(BC_DATE_DtDemogr)]) # supposedly choose the most recent value of bc_type?
bcfr_bc1 = bcfr_bc1[c("PersonID", "bc")] # ID & bc type
bcfr_bc1 = unique(bcfr_bc1)

# new df
bcfr_skptr1 = merge(bcfr_skptr, bcfr_bc1, by = "PersonID") # only 945 observations with BC data... some individuals must not have had BC data...

# bc categories are too many--sparse data. Collapse them -- n=0 for vaginal ring, condoms == no hormonal influence/"none", then let implant/IUD be "other" (no details on copper vs. hormonal)
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
# check
table(unique(bcfr_skptr1[c("PersonID", "bc")])$bc) 

# convert BC category to indicator vars; skipTrack.fit does not take factors
bcfr_skptr1 = bcfr_skptr1 |> mutate(value = 1) |> 
  pivot_wider(names_from = bc, values_from = value, values_fill = 0)

# fit skiptrack -------------------------------------------------

# fit skiptrack with birth control (factor var) as a predictor of cycle length & regularity

# hyper-parameters have been chosen based on previous runs/tests; however, as visible, convergence is not always great.
# burn-in chosen based on plots -- ~3000?
set.seed(9999)
fitbcfr_bc = skipTrack.fit(Y = bcfr_skptr1$cycle_len,
                           cluster = bcfr_skptr1$PersonID,
                           X = as.matrix(bcfr_skptr1[c("other", "pill")]), # covariates for cycle length mean -- left out "none" as "ref" i.e. so 1 df (i.e. intercept/when all covars =0)
                           Z = as.matrix(bcfr_skptr1[c("other", "pill")]), # covariates for cycle length precision i.e. 1/variance
                           chains = 6,
                           reps = 6000) 

# crude benchmarking:
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

# get precision by Person ID





# modifying step by step to isolate bug
# ALSO::: once I get this working, log & push to github (after sending df to Eva, with updated documentation)
# then do MDRTB edits



