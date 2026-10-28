# fxn to clean/transform Clue data

# reformatting fxn-----------------------------------------------------------------------------------
# get cycle lengths
clue_cyclelens = function(df_clue_in, spotting_radius){
  # -----------------------------------------------------------
  # arguments: df_clue_in, spotting_radius
    # df_clue_in: Clue data tracking sheet df with cols=c("personid", "category", "track_date)
    # spotting_radius (days): integer boundary for how far from a "period day" a spotting day entry must be to be considered *not* part of that period.
     # e.g. spotting_radius = 1 means that "spotting" days are only counted as contiguous "period" days if they immediately precede or succeed a "period" day (e.g. "period" 12, "spotting" 13: excluded. "period" 12, "spotting" 14: excluded.)
  # returns: a df with columns: cycle_len, cycle_num
    # personid
    # cycle_num: the # j of cycle observed for individual i [index]
    # cycle_len: the length of cycle j for individual i (days)
    # cycle_start: the 1st day of each cycle, as logged by user. Good to check that algo is working correctly (compare visually against original dataset.)
  # -----------------------------------------------------------
  
  # exclude birth control data in YWS -- only period or spotting counted as entries
  df_clue = df_clue_in[df_clue_in$category == "period" | df_clue_in$category == "spotting" ,] # none are missing? odd but i guess by definition these are ENTRIES in the clue app
  
  # Exclude entries as "counted" if spotting occurs in isolation, i.e. non-contiguous with period data
  # for loop: calculate distance (for each person) in time between entry & nearest-in-time "period" entry
  store_ppl = list()
  counter = 1
  for (ii in unique(df_clue$personid)){
    df_person = df_clue[df_clue$personid == ii,]
    # distance (in days) from a "period day" -- useful for spotting days
    df_person$min_diff_to_period = unlist(lapply(df_person$track_date, function(x){min(abs(x - df_person[df_person$category == "period",]$track_date))}))
    store_ppl[[counter]] = df_person
    counter = counter + 1
  }
  df_clue = as.data.frame(do.call(rbind, store_ppl))
  
  # keep periods, and spotting within X [spotting_radius] days of a period
  df_clue = df_clue |> filter(!(category == "spotting" & min_diff_to_period > spotting_radius))
  
  # generate cycle lengths with dplyr
  df_clue_out = df_clue |> 
    arrange(personid, track_date) |> # arrange in increasing order
    group_by(personid) |> # group by person
    # create "chunks", where each "chunk" is entries corresponding to a separate period
    # we consider those to be entries with a gap of >7 days; if they are within 7 days, they are considered to be part of the same period
    mutate(new_cycle = if_else(row_number() == 1L | track_date - lag(track_date) > 7, 1L, 0L), # flag the first row of each chunk with 1 or 0 (L = integer)
           cycle_num1 = cumsum(new_cycle) ) |>  #; then produce a continuous flag (1 1 1 2 2 2 3 3 3 etc...) with cumsum() indicating which # period this is for the individual
    group_by(personid, cycle_num1) |> # group by ID & cycle_num (i.e. each cycle as individual dfs, until we ungroup)
    summarise(cycle_start = first(track_date),  # for each given cycle chunk (sub-df), the start of the cycle is the first day
              .groups = "drop_last") |> # remove grouping by cycle, just keeps by personid
    mutate(cycle_len = lead(cycle_start) - cycle_start, # cycle len = day1 of current cycle - day1 of previous cycle
           cycle_num = row_number()) |> # more robust way to re-calculate cycle_num, even though it will be the same as cycle_num1 (above) in this code
    filter(!is.na(cycle_len)) |> # if lead() does not exist (i.e. final observed "chunk"/cycle), returns NA; we just don't know how long that cycle is, since we only know the first day of it & not the first day of the next cycle.
    ungroup() # remove dplyr grouping
  
  # # cleaning -- cycle variability measures can only be calculated if you observe >=2 cycles
  # exclude a) overly long cycles, b) individuals with only 1 observed cycle
  df_clue_out = df_clue_out |> filter(cycle_len < 365) 
  df_clue_out = df_clue_out |> group_by(personid) |> filter(n() >= 2)
  # redo cycle_num in light of potential exclusions
  df_clue_out = df_clue_out |> arrange(personid, cycle_num) |>
    group_by(personid) |> mutate(cycle_num = row_number()) |>
    filter(!is.na(cycle_len))
  
  # return df
  df_clue_out$cycle_num1 = NULL # redundant
  return(df_clue_out)
}

# ---------------------------------------------------
# apply to main skipTrack dataset -- with "spotting" buffer of 1 day--i.e. spotting must be consecutive with "period" day to be counted as part of a period:
library(dplyr)
library(tidyr) 
# setwd
# setwd("~/Desktop/ColumbiaPhD/research_projects/Lauren_menstrualmodelingmethods/menst_mod")
bcfr_trk = read.csv("./data/yws_clue/ywsclue_track.csv") 
colnames(bcfr_trk) = c("personid", "category", "type", "track_date")
bcfr_skptr = clue_cyclelens(bcfr_trk, 1)  # buffer of 1 day


