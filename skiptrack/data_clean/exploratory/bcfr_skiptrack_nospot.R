
# bcfr_skiptrack.R -- BUT: here, we treat spotting as a separate event
# SPECIFICALLY: if there are *isolated* spotting events, those are NOT considered to be new cycles (although if "spotting" is reported in advance of a period, that is considered the start of that period)
# & Check how much that choice changes the outputted df



# need to create readmes for the data/ folder, including the subfolder where I have put this script!!


# copied from 'bcfr_skiptrack.R' -- just with changes at lines ****** 

# -----------------------------------------------------------------------------------------------

# setwd("~/Desktop/ColumbiaPhD/research_projects/Lauren_menstrualmodelingmethods/menst_mod") 
library(dplyr)
library(tidyr) 

bcfr_trk = read.csv("./data/yws_clue/ywsclue_track.csv") 

cycle_len_only1 <- function(df, spotting_radius){
  
  # exclude bc data in YWS -- either period or spotting counted
  df_skip = df[df$category == "period" | df$category == "spotting" ,] # none are missing? odd but i guess by definition these are ENTRIES in the clue app
 
   # ******
  
  # Change: exclude entries as "counted" if spotting occurs in isolation, i.e. non-contiguous with period data: currently using a bound of 1 day (i.e. if spotting is buffered by nothing on either side, then I will exclude it)
  # for loop: calculate distance (for each person) in time between entry & nearest-in-time "period" entry
  store_ppl = list()
  counter = 1
  for (ii in unique(df_skip$PersonID)){
    df_person = df_skip[df_skip$PersonID == ii,]
    # distance (in days) from a "period day" -- useful for spotting days
    df_person$min_diff_to_period = unlist(lapply(df_person$track_date_DtDemogr, function(x){min(abs(x - df_person[df_person$category == "period",]$track_date_DtDemogr))}))
    store_ppl[[counter]] = df_person
    counter = counter + 1
  }
  df_skip = as.data.frame(do.call(rbind, store_ppl))
  
  # & if spotting is >1 day away (i.e. not followed or preceded by a period by more than a 1 day gap (e.g. spotting - nothing - period is *excluded*)): exclude
  # UNLESS ! it is ALSO not part of a string of spotting leading up to a period...
  
  # simple version: let's just filter by X days; if you are X days away i.e. reporting spotting >X days away from (in advance or after) a period (even if there are other spotting reports in-between), we are gonna filter that day of spotting data
  # keep periods, and spotting within X [spotting_radius] days of a period
  df_skip = df_skip |> filter(!(category == "spotting" & min_diff_to_period > spotting_radius))
  
  # ******
  
  df_skip = df_skip |> group_by(PersonID) |> arrange(track_date_DtDemogr, .by_group = T)  |> 
    mutate(gap_days = track_date_DtDemogr - lag(track_date_DtDemogr),
           new_cycle = if_else(gap_days > 7, 1, NA)) |>  
    arrange(track_date_DtDemogr, .by_group = T)
  
  df_skip = df_skip |> group_by(PersonID) |>
    mutate(new_cycle = if_else(is.na(new_cycle) & track_date_DtDemogr == min(track_date_DtDemogr), 1, new_cycle))
  
  df_cycle = df_skip |> group_by(PersonID) |> arrange(track_date_DtDemogr, .by_group = T)  |> 
    filter(new_cycle == 1) |> 
    mutate(cycle_len = lead(gap_days))
  
  df_cycle = df_cycle |> filter(cycle_len < 365) 
  df_cycle = df_cycle |> filter(cycle_len > 0) 
  df_cycle = df_cycle |> group_by(PersonID) |> filter(n() >= 2) 
  
  df_skptr = df_cycle |> select(PersonID, cycle_len) |> 
    group_by(PersonID) |> mutate(cycle_num = row_number()) |> 
    filter(!is.na(cycle_len))
  
  return(df_skptr)
}

# -----------------------------------------------------------------------------------------------

# & Check how much spotting choice choice changes the outputted df -- depending on the "radius" applied to spotting

# vary spotting_radius

# load skipTrack-processed data
source("./skiptrack/data_clean/bcfr_skiptrack.R") 
library(daff)

# Compare --  bcfr_skptr_nospot  VS.  bcfr_skptr
bcfr_skptr_nospot3 = cycle_len_only1(bcfr_trk, 3)   # 995 obs, not 1019 -- there are a few missing entries!
# compare buffers (i.e. "spotting_radius")
bcfr_skptr_nospot1 = cycle_len_only1(bcfr_trk, 1)  # 994 (essentially no difference -- one row is modified, one is deleted)
diff_result <- diff_data(bcfr_skptr_nospot3, bcfr_skptr_nospot1)
render_diff(diff_result)

# compare dfs with different spotting handling
diff_result <- diff_data(bcfr_skptr, bcfr_skptr_nospot1) # 32 "added", 57 "deleted", 201 cycle lengths modified 
render_diff(diff_result)

# worth re-running skipTrack with this ^^ i think
 
