setwd("~/Desktop/ColumbiaPhD/research_projects/Lauren_menstrualmodelingmethods/menst_mod")

# Clue tracking data (? who are in this subset--"Columbia data")
cbc = read.csv("./data/ctrack.csv")
# BCFR tracking data (Clue)
bcfr = read.csv("./data/yws_clue/ywsclue_track.csv")


# test run: get bcfr into skiptrack format & fit
library(skipTrack)
library(dplyr)
library(tidyr)

# id column is fine
# I have to make a decision: here I will say that spotting == part of period
unique(bcfr$category, ) # period, spotting, bc pill, iud, bc ring -- we will ignore everything except period, for now
# have to make decsisions about what constitutes a "cycle" to turn these reporting data into [length of cycle] data
# I will say that any gap of more than 4 days constitutes a new cycle
typeof(bcfr_skip$track_date_DtDemogr) # nice
# check that with Lauren

# unique(bcfr$category)
# [1] "period"             "spotting"           "birth_control_pill" "iud"                "birth_control_ring"
# exclude bc data
bcfr_skip = bcfr[bcfr$category == "period" | bcfr$category == "spotting" ,] # none are missing? odd but i guess by definition these are ENTRIES in the clue app

# get cycle lengths for each individual
bcfr_skip = bcfr_skip |> group_by(PersonID) |> arrange(track_date_DtDemogr, .by_group = T)  |>    # Cmd Shift M for pipe
  mutate(gap_days = track_date_DtDemogr - lag(track_date_DtDemogr),
         new_cycle = if_else(gap_days > 4, 1, 0)) |>  # if difference between rows >4
  arrange(track_date_DtDemogr, .by_group = T)  
bcfr_skip$new_cycle[1] = 1 # doesn't happen with algo, but is implied by meaning

# cycle length
bcfr_cycle = bcfr_skip |> group_by(PersonID) |> arrange(track_date_DtDemogr, .by_group = T)  |> 
  filter(new_cycle == 1) |> 
  mutate(cycle_len = lead(gap_days))

# for skiptrack
bcfr_skptr = bcfr_cycle |> select(PersonID, cycle_len) |> 
  group_by(PersonID) |> mutate(cycle_num = row_number()) |> 
  filter(!is.na(cycle_len))

# checking number of people in these data
length(unique(bcfr_skptr$PersonID))  # 57... is that what Lauren said to expect? that feels quite LOW

# plot data
hist(bcfr_skptr$cycle_len, breaks = seq(0,600,10))
# I think I should just exclude entries where the # days is over 300; those people clearly just never opened the app..
hist(bcfr_skptr$cycle_len, breaks = seq(0,600,2), xlim=c(0,150))
hist(bcfr_skptr$cycle_len, breaks = seq(0,600,2), xlim=c(0,70))

# ! didn't exclude the people with absurdly outlying values

# fit skiptrack
fitbcfr_basic = skipTrack.fit(Y = bcfr_skptr$cycle_len,
                              cluster = bcfr_skptr$PersonID,
                              chains = 4)
                              # no X (cycle length mean) covars -- YET -- Beta
                              # no Z (cycle length precision) covars -- YET -- Gamma
# runtime -- less than 1 min, about 40-50sec with this basic model
# results
fitbcfr_basic
fitbcfr_basic$fit # overwhelming amt of stuff -- would need to look at what all of this is...
# what is [1]][[207]]$ijDat this? is it literally running me through ALL of the MCMC chains used/guessed in fitting???
skipTrack.diagnostics(fitbcfr_basic) # note: doesn't look super converged...
skipTrack.results(fitbcfr_basic) # but G-R is alright? maybe it's overconfident... why woudl that be

# I need to think about what Lauren actually WANTs out of this -- 
# prediction? mean cycle length for each individual? prob that a given cycle was "skipped" i.e. its c_ij?? ** for cleaning the data... 
  # I thought max # skips was 10, how is there a c_ij which was estimated to be >10? -- c_ij is the number of TRUE cycles, so cij = 10 is *9* skips (jesus)
# let's go with that last one

# bcfr -- c_ijs to get actual cycle lengths
cijs_byrow = skipTrack.results(fitbcfr_basic)$cijs
colnames(cijs_byrow) = c("cij_lower", "cij_mean", "cij_upper")
bcfr_sk_inf = cbind(bcfr_skptr, cijs_byrow)
# "true" cycle length = cycle length / cij
bcfr_sk_inf$cyclelen_cij = bcfr_sk_inf$cycle_len / bcfr_sk_inf$cij_mean
# !!! note that for the highest cycle_len (& therefore the highest cij_mean),
# cij_mean are LESS CERTAIN... but I think that is more a product of sparse data than fitting...


# could rerun all of this with JUST periods, not spotting;
# also: I have not addressed whether some bleeding could be part of the "same" cycle when it's 7 days apart, etc. -- but i don't think that makes sense regardless


# ! have not dealt with the precision of these cij estimates! at some level, this might just be a restriction due to the relatively few # individuals 
# and if people have very little reported cycle data, we can't learn much about their individual cycle patterns...
# ex: Person 20075437 only reported one 2-yr long cycle -- her cycle length estimate is *entirely* derived from the information of the general sample i.e. mu rather than mu_i, since there is no data to really estimate a separate mu_i

    
# !! main thing that's missing -- skiptrack can take predictors for variation, length of cycle 
# I would need to decide what those are... discuss with lauren & look at dataset excel sheet/key for what makes sense

# -------------------------------------------------------------------------------------------------------
# Now: think about relevant covariates for...

# mean cycle length
# cycle variability!!

# Covariates available for the bcfr data:
# -- in this dataset: 
# -- in other datasets:
# note: will need to format to match bcfr_skptr formatting (not a big deal)




