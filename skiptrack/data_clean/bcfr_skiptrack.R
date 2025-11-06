# generate a dataset with 3 variables: PersonID, cycle number, cycle length

# From: full YWS BCFR dataset, no covariates, no other irregularity measures, cleaned & prepped for analysis

# To run with: skipTrack 

setwd("~/Desktop/ColumbiaPhD/research_projects/Lauren_menstrualmodelingmethods/menst_mod")
library(dplyr)
library(tidyr) 

# load data
bcfr_trk = read.csv("./data/yws_clue/ywsclue_track.csv") # BCFR tracking data (Clue)

# function to calculate cycle lengths with no bugs
cycle_len_only <- function(df){
  
  # exclude bc data in YWS -- either period or spotting counted
  df_skip = df[df$category == "period" | df$category == "spotting" ,] # none are missing? odd but i guess by definition these are ENTRIES in the clue app
  
  # get cycle lengths for each individual
  df_skip = df_skip |> group_by(PersonID) |> arrange(track_date_DtDemogr, .by_group = T)  |>    # Cmd Shift M for pipe
    mutate(gap_days = track_date_DtDemogr - lag(track_date_DtDemogr),
           new_cycle = if_else(gap_days > 7, 1, NA)) |>  # if difference between rows >7 -- note this condition in description to Eva
    arrange(track_date_DtDemogr, .by_group = T)
  
  # ensure that the first entry for each Person ID is a "new cycle"
  df_skip = df_skip |> group_by(PersonID) |>
    mutate(new_cycle = if_else(is.na(new_cycle) & track_date_DtDemogr == min(track_date_DtDemogr), 1, new_cycle))
  
  # cycle length
  df_cycle = df_skip |> group_by(PersonID) |> arrange(track_date_DtDemogr, .by_group = T)  |> 
    filter(new_cycle == 1) |> 
    mutate(cycle_len = lead(gap_days))
  
  # filters for cycle_num & cycle_len (future steps will mess up cycle_num unless they drop full PersonID)
  df_cycle = df_cycle |> filter(cycle_len < 365) # we are just eliminating these from observations--could have implications?
  df_cycle = df_cycle |> filter(cycle_len > 0) # No zero-length "cycles" -- could happen if they reported spotting & period on the same day...
  df_cycle = df_cycle |> group_by(PersonID) |> filter(n() >= 2) # only eliminates a few ppl! nice 
  
  # ! note: step above must happen BEFORE cycle_num is generated, or you can get "skips" where the # cycles observed are not continuous (e.g. cycle_num = 1,3,4,5...)
  
  # for skiptrack
  df_skptr = df_cycle |> select(PersonID, cycle_len) |> 
    group_by(PersonID) |> mutate(cycle_num = row_number()) |> 
    filter(!is.na(cycle_len))
  
  # return
  return(df_skptr)
}

# get BCFR clue data formatted for skiptrack
bcfr_skptr = cycle_len_only(bcfr_trk)


