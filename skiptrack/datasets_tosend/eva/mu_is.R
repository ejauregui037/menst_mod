
# Extract inferred mu_is (median cycle length for each individual i) from n=56 skipTrack individuals
# & plot as a histogram

# ***** THIS SCRIPT IS SOURCED BY mu_cyclelen_dists.Rmd for mockups/all histograms next to each other for comparison!!







# load fitted model -------------------------------------------------------------------------------

# Goal:
# extract mu_is i.e. cycle lengths -- estimated within skipTrack! but they are part of the posterior...

  # some notes on units for muijs:
  # 1. this is a MEDIAN cycle length!! to the lognormal dist [bc it is so skewed, parameterized with median]
    # https://www.theoj.org/joss-papers/joss.06928/10.21105.joss.06928.pdf
    # https://en.wikipedia.org/wiki/Log-normal_distribution
  # 2. it is in fact the NAT LOG of the median cycle length! ^^ so we need to antilog it : ) exp()

# (1) crude way: divide all observed cycle lengths by c_ij

# (2) cleaner way: extract mu_i for each individual...
  # GOT IT -- NOTE: the names are nasty, but the documentation reveals what is going on:
  # https://github.com/LukeDuttweiler/skipTrack/blob/228d39b08c83fb125a68b6c59085f96a81c19dcb/R/MCMC.R#L106
  # line 106 :    'muijs' = sapply(cluster, function(ind){initialParams$muis[iDat$Individual == ind]}))
    # mui is misleadingly named as muij, BUT! we can see here that it is actually the extraction of the muis!!!!
# AND: while this fitbcfr_basic$fit[[1]][[1]] makes it seem like muijs = mu overall, THAT IS J BC all individuals were initialized with the same mu_i for MCMC!
# if we instead look at fitbcfr_basic$fit[[1]][[100]], we can see that a) muij changes for each INDIVIDUAL, but is constant WITHIN each individual! & b) mijs change for each observation! which, from the same github link, you can see is just the cij-correct individual observation-level inferred cycle length



# setup -------------------------------------------------------------------------------

# setwd("~/Desktop/ColumbiaPhD/research_projects/Lauren_menstrualmodelingmethods/menst_mod")  
library(dplyr)
library(tidyr)




# calculate mu_is & posterior stats -------------------------------------------------------------------------------
# get the metrics that are created in bcfr_basicfit10.R, but not saved in the RDS file [which is just the fitted skipTrack object]:

# ran & saved -- takes a long time to run: 

        # # load fitted skipTrack basic model w/ 10,000 iterations to MCMC:
        # newdf = readRDS(file = "./skiptrack/fitting/bcfrskip10")
        # fitbcfr_basic = newdf # rename for easy tracking

        # # get mean cycle length(mu_i) by Person ID -- output as a df
        # n_chains = length(fitbcfr_basic$fit) # 6 chains
        # n_reps = length(fitbcfr_basic$fit[[1]]) # 10,001 reps (+1 starting rep)
        # burnin = 1502 # +1 for initial guess [i.e. skip initial guess + burn-in reps]
        # 
        # store_muilists = vector(mode = "list", length = n_chains*(n_reps-(burnin-1))) # burnin-1 because 1502:3001 starts at 1502, not 1503, & goes to 3001, not 30001
        # 
        # index = 1
        # for (chain in 1:n_chains){
        #   for (rep in burnin:n_reps){ # skip burnin: posterior dist comes once chains have settled
        #     df_chain_rep = fitbcfr_basic$fit[[chain]][[rep]]$ijDat # for some reason these are stord in ijDat rather than iDat -- iDat would have been cleaner! like taus
        #     # separate by individual
        #     df_chain_rep = unique(df_chain_rep[c("Individual", "muijs")]) # remove "duplicate" rows i.e. get one row per individual! since muijs are constant for each obs j for individual i, within a given MCMC iteration
        #     # order of individuals is still preserved^^ so next step is fine
        #     # now extract
        #     muis = df_chain_rep[["muijs"]] # misleadingly named!
        #     store_muilists[[index]] = muis
        #     index = index + 1
        #   }
        # }
        # 
        # muilist_df = as.data.frame(store_muilists, col.names = paste0("run", 1:51000)) # 9000 = 8500 x 6 chains
        # MCMC_muis = cbind("Individual" = fitbcfr_basic$fit[[1]][[1502]]$iDat[["Individual"]], muilist_df) # accessing IDs from a particular run, actual run accessed not important
        # # n=56 unique individuals
        # 
        # # now: get some features of the POSTERIOR dist
        # MCMC_taus_post = MCMC_muis |> rowwise()  |> 
        #   mutate(
        #     mean_mui_ln = mean(c_across(2:ncol(MCMC_muis))), # excluding first column, "Individual"
        #     q025_mui_ln = quantile(c_across(2:ncol(MCMC_muis)), 0.025), # 95% credible interval, ish
        #     median_mui_ln = quantile(c_across(2:ncol(MCMC_muis)), 0.5),
        #     q975_mui_ln = quantile(c_across(2:ncol(MCMC_muis)), 0.975)) |> 
        #   ungroup() |>   
        #   select(Individual, mean_mui_ln, q025_mui_ln, median_mui_ln, q975_mui_ln)
        # 
        # # add the antilogged versions!
        # MCMC_taus_post1 = MCMC_taus_post |> mutate(
        #   mean_mui = exp(mean_mui_ln),
        #   q025_mui = exp(q025_mui_ln), 
        #   median_mui = exp(median_mui_ln), # median & mean are extremely similar! posterior dists OF EACH MU_i mu_i must not be too skewed. (!= the mean/median of cycle lengths for each individual, these are mean/median of the posteriors for the MEDIAN CYCLE LENGTH of each individual i!!)
        #   q975_mui = exp(q975_mui_ln)
        # )  

# the above takes SO long to run -- too annoying to run repeatedly for plots. Save it out.
# saveRDS(MCMC_taus_post1, "./skiptrack/datasets_tosend/eva/muis_post.rds")

MCMC_taus_post1 = readRDS(file = "./skiptrack/datasets_tosend/eva/muis_post.rds")



# --------------------------------------------------------------------------------------------------------------------------------
# bar plots
# catgories, from eva plots:
levels_1to7 = c("<21 days", "21-25 days", "25-33 days", "33-36 days", "36-61 days", "61-90 days", ">90 days") # no gaps!
levels_1to7_gaps = c("<21 days", "21-25 days", "26-32 days", "33-35 days", "36-60 days", "61-90 days", ">90 days")
# --------------------------------------------------------------------------------------------------------------------------------






# plot histogram of skipTrack-inferred cycle lengths -----------------------------------------------------------------------------------------------------------------
# for plotting: all we need is the mean antilogged mu_i (i.e. inferred median cycle length for each of n=56 individuals)

# hist(MCMC_taus_post1$mean_mui, breaks = c(0, 21, 25, 36, 61), freq = T)
# hist(MCMC_taus_post1$mean_mui, main = "median cycle len dist") # this distribution looks VERY DIFFERENT from self-report or uncorrected skipTrack ! ! 
# # but it also makes no assumptions & is essentially smoothing based on other patterns in reported cycles...

# median cycle lengths
mu_skp = MCMC_taus_post1 # copy to keep
mu_skp = mu_skp[c("Individual", "mean_mui")] # use mean of posterior for each inferred MEDIAN CYCLE LENGTH mu_i

# generate bins/factors -- LOWER BOUND INCLUSIVE
mu_skp = mu_skp |> 
  mutate(median_cycle_bin = case_when(mean_mui<21 ~ 1,
                                      mean_mui>=21 & mean_mui<25 ~ 2,
                                      mean_mui>=25 & mean_mui<33 ~ 3,
                                      mean_mui>=33 & mean_mui<36 ~ 4,
                                      mean_mui>=36 & mean_mui<61 ~ 5,
                                      mean_mui>=61 & mean_mui<=90 ~ 6, # this one is upper inclusive just to match the questionnaire!
                                      mean_mui>90 ~ 7))

# note bound inclusivity in subtitle/captions
mu_skp$median_cycle_bin = factor(mu_skp$median_cycle_bin, levels = c(1,2,3,4,5,6,7), labels = levels_1to7)


# mean cycle lengths  -----------------------------------------

# (!) I could extract the inferred mean from the lognormal dist... there is a closed form expression for mean of lognormal, using sd & median !
# (!!!!)
# ln(mean) = ln(median) + (sd^2)/2   

# get tau_is, already calculated elsewhere : )
bcfryws_skiptrack10k_s = read.csv(file = "./skiptrack/datasets_tosend/eva/bcfryws_skiptrack10k_s.csv")
bcfryws_tau = unique(bcfryws_skiptrack10k_s[c("personid", "mean_tau")])
colnames(bcfryws_tau)[1] = "Individual"
mu_skp = merge(mu_skp, bcfryws_tau)
# calculate mean
mu_skp$mean_cyclelen = mu_skp$mean_mui * exp(1/(2*mu_skp$mean_tau))

mu_skp = mu_skp |> 
  mutate(mean_cycle_bin = case_when(mean_cyclelen<21 ~ 1,
                                    mean_cyclelen>=21 & mean_cyclelen<25 ~ 2,
                                    mean_cyclelen>=25 & mean_cyclelen<33 ~ 3,
                                    mean_cyclelen>=33 & mean_cyclelen<36 ~ 4,
                                    mean_cyclelen>=36 & mean_cyclelen<61 ~ 5,
                                    mean_cyclelen>=61 & mean_cyclelen<=90 ~ 6, # this one is upper inclusive just to match the questionnaire!
                                    mean_cyclelen>90 ~ 7))
# MODIFIED BINS -- to match labels in survey DESPITE continuous measure in skipTrack
mu_skp = mu_skp |> 
  mutate(mean_cycle_bin1 = case_when(mean_cyclelen<21 ~ 1, # "<21"
                                    mean_cyclelen>=21 & mean_cyclelen<25.5 ~ 2, # "21-25"
                                    mean_cyclelen>=25.5 & mean_cyclelen<32.5 ~ 3, # "26-32"
                                    mean_cyclelen>=32.5 & mean_cyclelen<35.5 ~ 4, # "33-35"
                                    mean_cyclelen>=35.5 & mean_cyclelen<60.5 ~ 5, # "36-60"
                                    mean_cyclelen>=60.5 & mean_cyclelen<=90 ~ 6, # "61-90"
                                    mean_cyclelen>90 ~ 7)) # ">90"

mu_skp$mean_cycle_bin = factor(mu_skp$mean_cycle_bin, levels = c(1,2,3,4,5,6,7), labels = levels_1to7)
mu_skp$mean_cycle_bin1 = factor(mu_skp$mean_cycle_bin1, levels = c(1,2,3,4,5,6,7), labels = levels_1to7_gaps)

# clean environment:
rm(list = setdiff(ls(), c("mu_skp", "levels_1to7", "levels_1to7_gaps")))












# --------------------------------------------------------------------------------------------------------------------------------

# comparison; crude correction (i.e. c_ij --> mean & median cycle lengths for each individual i)

bcfryws_skiptrack10k_s = read.csv(file = "./skiptrack/datasets_tosend/eva/bcfryws_skiptrack10k_s.csv")

cij_corrected = bcfryws_skiptrack10k_s # copy

# median AND mean -- cause maybe that's driving this as well... mean skew s& median doesn't!!
# in fact this is a CRITICAL diff jaja
cij_corrected$cyclelen_corrected = cij_corrected$cycle_len/cij_corrected$cij_mean
cij_corrected = cij_corrected[c("personid", "cyclelen_corrected")]
cij_corrected = cij_corrected |> group_by(personid) |> mutate(median_corrected = median(cyclelen_corrected),
                                                                                mean_corrected = mean(cyclelen_corrected))

cij_corrected = unique(cij_corrected[c("personid", "median_corrected", "mean_corrected")]) # 56 obs -- nice

# now look at these:
# hist(cij_corrected$median_corrected)
# hist(cij_corrected$mean_corrected) # definitely has more of those extreme values, although not a TON!

# 36-60 are still MUCh rarer ^^ in the skipTrack corrected data! it is definitely clocking those as accidental skips jaja
# since they are so close to 2x cycle length fr!

# generate bins/factors -- LOWER BOUND INCLUSIVE
cij_corrected = cij_corrected |> 
  mutate(median_corrected_bin = case_when(median_corrected<21 ~ 1,
                                      median_corrected>=21 & median_corrected<25 ~ 2,
                                      median_corrected>=25 & median_corrected<33 ~ 3,
                                      median_corrected>=33 & median_corrected<36 ~ 4,
                                      median_corrected>=36 & median_corrected<61 ~ 5,
                                      median_corrected>=61 & median_corrected<=90 ~ 6, # this one is upper inclusive just to match the questionnaire!
                                      median_corrected>90 ~ 7),
         mean_corrected_bin = case_when(mean_corrected<21 ~ 1,
                                        mean_corrected>=21 & mean_corrected<25 ~ 2,
                                        mean_corrected>=25 & mean_corrected<33 ~ 3,
                                        mean_corrected>=33 & mean_corrected<36 ~ 4,
                                        mean_corrected>=36 & mean_corrected<61 ~ 5,
                                        mean_corrected>=61 & mean_corrected<=90 ~ 6, # this one is upper inclusive just to match the questionnaire!
                                        mean_corrected>90 ~ 7),)
# same thing: 'match' bins in survey with continuous results...
cij_corrected = cij_corrected |> 
  mutate(mean_corrected_bin1 = case_when(mean_corrected<21 ~ 1,
                                        mean_corrected>=21 & mean_corrected<25.5 ~ 2,
                                        mean_corrected>=25.5 & mean_corrected<32.5 ~ 3,
                                        mean_corrected>=32.5 & mean_corrected<35.5 ~ 4,
                                        mean_corrected>=35.5 & mean_corrected<60.5 ~ 5,
                                        mean_corrected>=60.5 & mean_corrected<=90 ~ 6, 
                                        mean_corrected>90 ~ 7),)

cij_corrected$median_corrected_bin = factor(cij_corrected$median_corrected_bin, levels = c(1,2,3,4,5,6,7), labels = levels_1to7)
cij_corrected$mean_corrected_bin = factor(cij_corrected$mean_corrected_bin, levels = c(1,2,3,4,5,6,7), labels = levels_1to7)
cij_corrected$mean_corrected_bin1 = factor(cij_corrected$mean_corrected_bin1, levels = c(1,2,3,4,5,6,7), labels = levels_1to7_gaps)

# clean environment:
rm(list = setdiff(ls(), c("mu_skp", "levels_1to7", "levels_1to7_gaps",
                          "cij_corrected", "bcfryws_skiptrack10k_s")))











# --------------------------------------------------------------------------------------------------------------------------------

# skipTrack-reported cycle lengths, not corrected

bcfryws_skiptrack_meanmed = bcfryws_skiptrack10k_s |> group_by(personid) |> summarise(mean_cycle = mean(cycle_len),
                                                                                      median_cycle = median(cycle_len))
hist(bcfryws_skiptrack_meanmed$median_cycle, breaks = c(21, 26, 33, 36, 61, 90), freq = T)
# what happened to the in-between means in Eva's figure? that doesn't make sense... i.e. 35-36 days?
hist(bcfryws_skiptrack_meanmed$mean_cycle, breaks = c(21, 26, 33, 36, 61, 90), freq = T)

# generate bins/factors -- LOWER BOUND INCLUSIVE
bcfryws_skiptrack_meanmed = bcfryws_skiptrack_meanmed |> 
                                  mutate(median_cycle_bin = case_when(median_cycle<21 ~ 1,
                                                                 median_cycle>=21 & median_cycle<25 ~ 2,
                                                                 median_cycle>=25 & median_cycle<33 ~ 3,
                                                                 median_cycle>=33 & median_cycle<36 ~ 4,
                                                                 median_cycle>=36 & median_cycle<61 ~ 5,
                                                                 median_cycle>=61 & median_cycle<=90 ~ 6, # this one is upper inclusive just to match the questionnaire!
                                                                 median_cycle>90 ~ 7),
                                    mean_cycle_bin = case_when(mean_cycle<21 ~ 1,
                                                                 mean_cycle>=21 & mean_cycle<25 ~ 2,
                                                                 mean_cycle>=25 & mean_cycle<33 ~ 3,
                                                                 mean_cycle>=33 & mean_cycle<36 ~ 4,
                                                                 mean_cycle>=36 & mean_cycle<61 ~ 5,
                                                                 mean_cycle>=61 & mean_cycle<=90 ~ 6, # this one is upper inclusive just to match the questionnaire!
                                                                 mean_cycle>90 ~ 7),)
# 'match' bins to BCFR YWS survey
bcfryws_skiptrack_meanmed = bcfryws_skiptrack_meanmed |> 
  mutate(mean_cycle_bin1 = case_when(mean_cycle<21 ~ 1,
                                    mean_cycle>=21 & mean_cycle<25.5 ~ 2,
                                    mean_cycle>=25.5 & mean_cycle<32.5 ~ 3,
                                    mean_cycle>=32.5 & mean_cycle<35.5 ~ 4,
                                    mean_cycle>=35.5 & mean_cycle<60.5 ~ 5,
                                    mean_cycle>=60.5 & mean_cycle<=90 ~ 6, #
                                    mean_cycle>90 ~ 7),)

# note bound inclusivity in subtitle/captions
bcfryws_skiptrack_meanmed$median_cycle_bin = factor(bcfryws_skiptrack_meanmed$median_cycle_bin, levels = c(1,2,3,4,5,6,7), labels = levels_1to7)
bcfryws_skiptrack_meanmed$mean_cycle_bin = factor(bcfryws_skiptrack_meanmed$mean_cycle_bin, levels = c(1,2,3,4,5,6,7), labels = levels_1to7)
bcfryws_skiptrack_meanmed$mean_cycle_bin1 = factor(bcfryws_skiptrack_meanmed$mean_cycle_bin1, levels = c(1,2,3,4,5,6,7), labels = levels_1to7_gaps)

# do these match up with Eva's plots? -- yes -- but fix it to make it a a barplot like what she has : )

# clean environment:
rm(list = setdiff(ls(), c("mu_skp", "levels_1to7", "levels_1to7_gaps",
                          "cij_corrected", "bcfryws_skiptrack_meanmed")))










# --------------------------------------------------------------------------------------------------------------------------------

# self-reported cycle lengths 

bcfr_cyclen = read.csv(file = "~/Desktop/ColumbiaPhD/research_projects/Lauren_menstrualmodelingmethods/menst_mod/data/yws_clue/ywsclue_repro.csv")
# unique(bcfr_cyclen$MENSES_DAYS)
# bcfr_cyclen = bcfr_cyclen |> mutate(cycle_len = case_when(MENSES_DAYS == X ~ XX))

# Eva already made a figure for this... 
# ! ! ! -- more practically, I can't do this because I don't have a key for this variable... so all th levels are meaningless lol

# I can infer...

# bcfr_cyclen_tb = table(bcfr_cyclen$MENSES_DAYS) # excludes NAs
# barplot(bcfr_cyclen_tb)
# for ggplot
bcfr_cyclen$MENSES_DAYS_fact = factor(bcfr_cyclen$MENSES_DAYS, levels = c(1,2,3,4,5,6,7,8,9), 
                                      labels = c(levels_1to7_gaps, "8__code", "9__code")) # there are gaps, in the questionnaire

# bcfr_cyclen_tb = as.data.frame(table(bcfr_cyclen$MENSES_DAYS)) # excludes NAs
# colnames(bcfr_cyclen_tb) = c("Cycle_length", "Frequency")
# I guess that 8 & 9 are don't know/missing codes?? but these relative frequencies correspond to Lauren's 

# & in Clue subset! -----------
bcfryws_skiptrack10k_s = read.csv(file = "./skiptrack/datasets_tosend/eva/bcfryws_skiptrack10k_s.csv")
ids_clue = unique(bcfryws_skiptrack10k_s$personid)
bcfr_cyclen_clue = bcfr_cyclen[bcfr_cyclen$PersonID %in% ids_clue,] # n=56 -- good

# # clean environment:
# rm(list = setdiff(ls(), c("mu_skp", "levels_1to7", "levels_1to7_gaps",
#                           "cij_corrected", "bcfryws_skiptrack_meanmed",
#                           "bcfr_cyclen", "bcfr_cyclen_clue")))


# --------------------------------------------------------------------------------------------------------------------------------


# combine to send to eva -- REMOVE AFTER

# combine1 = bcfr_cyclen_clue[c("PersonID", "MENSES_DAYS", "MENSES_DAYS_fact")]
# combine2 = bcfryws_skiptrack_meanmed
# colnames(combine2) = c("PersonID", "meancycle_Clue", "mediancycle_Clue", "mediancycle_Clue_bin", "meancycle_Clue_bin")
# combine3 = cij_corrected
# colnames(combine3) = c("PersonID", "mediancycle_cijcorrected", "meancycle_cijcorrected", "mediancycle_cijcorrected_bin", "meancycle_cijcorrected_bin")
# combine4 = mu_skp
# colnames(combine4) = c("PersonID", "mediancycle_inferred", "mediancycle_inferred_bin", "tau", "meancycle_inferred", "meancycle_inferred_bin")
# combine4$tau = NULL
# 
# combined = merge(combine1, combine2, by = "PersonID")
# combined = merge(combined, combine3, by = "PersonID")
# combined = merge(combined, combine4, by = "PersonID")
# 
# colnames(combined)
# write.csv(combined, "~/Downloads/cyclelen_vars.csv")
# 
# rm(combine1, combine2, combine3, combine4, combined)

# --------------------------------------------------------------------------------------------------------------------------------

# check: how many cycles for each individual, on average:

# count_enties = bcfryws_skiptrack10k_s |> group_by(personid) |> summarise(num_cycles = max(cycle_num))
# mean(count_enties$num_cycles)

