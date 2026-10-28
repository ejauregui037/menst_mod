# ---
  #  title: "Stress and menstruation analyses edited"
  #  author: "Eva Siegel"
  #  date: "2026-04-14"
  
library(haven)
library(tidyverse)
library(dplyr)
library(readxl)
library(knitr)
library(labelled)
library(gmodels)
library(gtsummary)
library(kableExtra)
library(janitor)
library(stringr)
library(rstatix)
library(broom)
library(broom.helpers)
library(parameters)
library(mgcv)
library(nnet)
library(tableone) 
library(labelled)
library(purrr)
library(forcats)
library(sandwich)

# setwd('/users/evasiegel/desktop/lauren/data')
setwd('~/Desktop/ColumbiaPhD/research_projects/Lauren_menstrualmodelingmethods/menst_mod/eva_code/')

file_path <- "250407 yws & clue data for lauren (1).xlsx"

sheets_to_clean <- c(
  "ALCOHOL_TOBACCO", "MEDICAL_HISTORY", "PHYSICAL_ACTIVITY", 
  "REPRODUCTIVE", "Occupational_History", "Sleep", "Derived"
)

cleaned_data <- map(set_names(sheets_to_clean), ~ {
  read_excel(file_path, sheet = .x) %>% select(-any_of("CENTER_NO"))
})

list2env(cleaned_data, .GlobalEnv)

yws       <- read_excel(file_path, sheet = "STRAIN")

yws_clue  <- read_excel(file_path, sheet = "CLUE_Participants")
yws_clues <- read_excel(file_path, sheet = "CLUE_Sessions")
yws_clueb <- read_excel(file_path, sheet = "CLUE_Birth_Control")
ywscluet  <- read_excel(file_path, sheet = "CLUE_Tracking")

ywsnew1   <- read_excel("Copy of 251011 Additional YWS Data For Lauren.xlsx", sheet = "MEDICAL_HISTORY")
ywsnew2   <- read_excel("Copy of 251011 Additional YWS Data For Lauren.xlsx", sheet = "REPRODUCTIVE")

skip<-read_csv("bcfryws_skiptrack10k_s.csv")
skip <- skip %>%
  group_by(PersonID) %>%
  mutate(mean_cycle_len = mean(cycle_len, na.rm = TRUE)) 
skip1 <- skip %>%
  group_by(PersonID) %>%
  slice(1)

skip1$cycle_len <- NULL
skip1$cycle_num <- NULL
skip1$cij_mean <- NULL

skip2<-read_csv("cyclelen_vars.csv")

list_of_dfs <- list(yws, ALCOHOL_TOBACCO, MEDICAL_HISTORY, PHYSICAL_ACTIVITY, 
                    REPRODUCTIVE, Occupational_History, Sleep, Derived, ywsnew1, ywsnew2, skip1, skip2)

merged_df <- list_of_dfs[[1]]

for (i in 2:length(list_of_dfs)) {
  merged_df <- merge(merged_df, list_of_dfs[[i]], by = "PersonID", all = TRUE)
}

#######data cleaning#######

merged_df <- merged_df %>%
  mutate(
    FIRST_MENSES_AGE   = na_if(FIRST_MENSES_AGE, 99),
    MENSES_NUMBER_LAST12 = na_if(MENSES_NUMBER_LAST12, 6),
    MENSES_DAYS        = na_if(MENSES_DAYS, 9),
    HEIGHT_CM          = na_if(HEIGHT_CM, 6),
    HEIGHT_M           = na_if(HEIGHT_M, 2),
    WEIGHT_POUND       = na_if(WEIGHT_POUND, 996)
  ) %>%
  mutate(
    menarche = coalesce(FIRST_MENSES_AGE, FIRST_MENSES_GRADE),
    menarched = recode(menarche, `28`=9.5, `29`=10.5, `30`=11.5, `31`=12.5, 
                       `32`=13.5, `33`=14.5, `34`=15.5, `35`=16.5, `36`=17.5),
    men_group = if_else(menarched > median(menarched, na.rm = TRUE), "Later", "Earlier"),
    men_group_v2 = case_when(
      menarched < 12.5  ~ "Earlier (< 12.5)",
      menarched >= 12.5 ~ "Later (>= 12.5)",
      TRUE              ~ NA_character_
    )
  ) %>%
  mutate(
    Menses_Recent = factor(
      recode(
        as.character(MENSES_STARS_AS_EXPECTED), 
        "1" = "Regular",
        "2" = "Irregular",
        "3" = "Don't Know",
        .default = NA_character_
      ),
      levels = c("Irregular", "Regular", "Don't Know")
    )
  )   %>%
      mutate(
        Menses_Lifetime = factor(
          recode(
            as.character(IRREGULAR_MENSES),
            "1" = "Irregular",
            "2" = "Regular",
            "3" = "Don't Know",
            .default = NA_character_
          ),
          levels = c("Irregular", "Regular", "Don't Know")
      )
    ) %>%
  mutate(
   Race_D = factor(Race_D, levels = 1:5, labels = c("White", "Black", "Latina", "Asian", "Other")),
    across(c(PREG_EVER, SMOKE_LAST10_EVER, ALCOHOL_EVER), 
           ~ factor(.x, levels = c(2, 1), labels = c("No", "Yes")))
  ) %>%
  mutate(
    tot_height = coalesce(HEIGHT_M * 100 + HEIGHT_CM, (HEIGHT_FEET * 12 + HEIGHT_INCH) * 2.54) / 100,
    tot_weight = coalesce(WEIGHT_KG, WEIGHT_POUND / 2.2),
    bmi        = tot_weight / (tot_height^2),
    bmi_cat    = cut(bmi, breaks = c(0, 18.5, 25, 30, Inf), 
                     labels = c("Underweight", "Normal weight", "Overweight", "Obesity"), 
                     right = FALSE)
  ) %>%
  mutate(
    STRESS_CT_CAT = if_else(STRESS_CT > median(STRESS_CT, na.rm = TRUE), "High Stress", "Low Stress"),
    STRESS_TH_CAT = if_else(STRESS_TH > median(STRESS_TH, na.rm = TRUE), "High Severity", "Low Severity")
  ) %>%
  mutate(
    MENSES_DAYS = fct_collapse(as.character(MENSES_DAYS),
                               "<26 days"            = c("1", "2"),
                               "26-32 days"          = "3",
                               ">32 days"            = c("4", "5", "6", "7"),
                               "Too variable to say" = "8"
    ) %>% fct_relevel("<26 days", "26-32 days", ">32 days", "Too variable to say"),
    MISSED_MENSES_60DAYS = factor(
      recode(as.character(MISSED_MENSES_60DAYS),
             "1" = "Yes",
             "2" = "No",
             "3" = "Don't Know"
      ),
      levels = c("Yes", "No", "Don't Know")
    )
  ) %>%
  mutate(
    MENSES_NUMBER_LAST12 = fct_collapse(as.character(MENSES_NUMBER_LAST12),
                                 "1-3 periods"   = "1",
                                 "4-6 periods"   = "2",
                                 "7-10 periods"  = "3",
                                 "11-14 periods" = "4",
                                 ">14 periods"   = "5"
    )
  )  

var_label(merged_df) <- list(
  Race_D = "Race/Ethnicity",
  bmi_cat = "BMI",
  PREG_EVER = "Ever pregnant",
  menarched = "Menarche",
  AgeIntv_Demographics = "Age",
  IRREGULAR_MENSES = "Lifetime irregularity",
  MISSED_MENSES_60DAYS = "Missed menses > 60 days", 
  MENSES_STARS_AS_EXPECTED = "Recent irregularity", 
  MENSES_NUMBER_LAST12 = "Menses in past year",
  MENSES_DAYS = "Average cycle length"
  
)

myVars <- c("AgeIntv_Demographics", "Race_D", "bmi_cat", "PREG_EVER", "SMOKE_LAST10_EVER", "ALCOHOL_EVER", "menarched", "men_group", "STRESS_CT", "EATotCT", "AdultCT", "RecTotCT", "STRESS_TH", "STRESS_CT_CAT")

catVars <- c("Race_D","bmi_cat","PREG_EVER", "SMOKE_LAST10_EVER", "ALCOHOL_EVER", "men_group", "STRESS_CT_CAT")

########################################Table 1##################################################################
# table1 <- merged_df %>%
#   select(all_of(myVars)) %>% 
#   tbl_summary(
#     by = STRESS_CT_CAT,           
#     missing = "ifany",       
#     statistic = list(
#       all_continuous() ~ "{mean} ({sd})",    
#       all_categorical() ~ "{n} ({p}%)"        
#     )
#   ) %>%
#   add_p() %>%                
#   add_overall() %>%          
#   bold_labels()              
# 
# table1

########################################Table 2##################################################################
myVars <- c("Menses_Lifetime", "MISSED_MENSES_60DAYS", "Menses_Recent", "MENSES_NUMBER_LAST12", "MENSES_DAYS", "STRESS_CT_CAT")

catVars <- c("Menses_Lifetime", "MISSED_MENSES_60DAYS", "Menses_Recent", "MENSES_NUMBER_LAST12", "MENSES_DAYS", "STRESS_CT_CAT")

# table2 <- merged_df %>%
#   select(all_of(myVars)) %>% 
#   tbl_summary(
#     by = STRESS_CT_CAT,           
#     missing = "ifany",       
#     statistic = list(
#       all_continuous() ~ "{mean} ({sd})",    
#       all_categorical() ~ "{n} ({p}%)"        
#     )
#   ) %>%
#   add_p() %>%                
#   add_overall() %>%          
#   bold_labels()       
# 
# table2

#########################################Table 3###################################################################
##recent menses
merged_df$MENSES_STARS_AS_EXPECTED<- as.numeric(merged_df$MENSES_STARS_AS_EXPECTED)

merged_df$MENSES_STARS_AS_EXPECTED[merged_df$MENSES_STARS_AS_EXPECTED == 3] <- NA

merged_df <- merged_df %>%
  mutate(
    MENSES_STARS_AS_EXPECTED = case_when(
      MENSES_STARS_AS_EXPECTED == 1 ~ 0, # Regular
      MENSES_STARS_AS_EXPECTED == 2 ~ 1, # Irregular
      TRUE ~ NA_real_
    )
  )

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
# 
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
# 
# tbl_merged %>%
#   modify_caption("**Adjusted Relative Risk Ratios, Lifetime Stress and Recent Irregular Menses**")

#binary exposures
merged_df$STRESS_CT_CAT <- relevel(factor(merged_df$STRESS_CT_CAT), ref = "Low Stress")
merged_df$STRESS_TH_CAT <- relevel(factor(merged_df$STRESS_TH_CAT), ref = "Low Severity")

# mod1_tbl <- glm(MENSES_STARS_AS_EXPECTED ~ STRESS_CT_CAT 
#                 + PREG_EVER + ALCOHOL_EVER + SMOKE_LAST10_EVER + AgeIntv_Demographics
#                 , 
#                 family = poisson(link = "log"), 
#                 data = merged_df) %>%
#   tbl_regression(
#     exponentiate = TRUE, 
#     label = list(
#       STRESS_CT_CAT ~ "High Stress"
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
# mod2_tbl <- glm(MENSES_STARS_AS_EXPECTED ~ STRESS_TH_CAT 
#                 + PREG_EVER + ALCOHOL_EVER + SMOKE_LAST10_EVER + AgeIntv_Demographics
#                 ,
#                 family = poisson(link = "log"), 
#                 data = merged_df) %>%
#   tbl_regression(
#     exponentiate = TRUE, 
#     label = list(
#       STRESS_TH_CAT ~ "High Stress Severity"
#       , 
#       PREG_EVER ~ "Ever Pregnant",
#       ALCOHOL_EVER ~ "Ever Alcohol Use",
#       SMOKE_LAST10_EVER ~ "Ever Smoked (Last 10 Years)",
#       AgeIntv_Demographics ~ "Age"
#     ),
#     tidy_fun = partial(tidy_robust, vcov = "HC1")
#   ) %>%
#   add_n(location = "label")

# tbl_merged <- tbl_merge(
#   tbls = list(mod1_tbl, mod2_tbl),
#   tab_spanner = c("**# of Stressful Events**", "**Severity of Stress**")
# )
# 
# tbl_merged %>%
#   modify_caption("**Adjusted Relative Risk Ratios, Lifetime Stress and Recent Irregular Menses**")

##same for lifetime menses
merged_df$IRREGULAR_MENSES<- as.numeric(merged_df$IRREGULAR_MENSES)
merged_df$IRREGULAR_MENSES[merged_df$IRREGULAR_MENSES == 3] <- NA

merged_df <- merged_df %>%
  mutate(
    IRREGULAR_MENSES = case_when(
      IRREGULAR_MENSES == 1 ~ 1, # Irregular
      IRREGULAR_MENSES == 2 ~ 0, # Regular
      TRUE ~ NA_real_
    )
  )

# mod1_tbl <- glm(IRREGULAR_MENSES ~ STRESS_CT 
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
# 
# mod2_tbl <- glm(IRREGULAR_MENSES ~ STRESS_TH 
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
# 
# tbl_merged %>%
#   modify_caption("**Adjusted Relative Risk Ratios, Lifetime Stress and Lifetime Irregular Menses**")

#binary
# mod1_tbl <- glm(IRREGULAR_MENSES ~ STRESS_CT_CAT 
#                 + PREG_EVER + ALCOHOL_EVER + SMOKE_LAST10_EVER + AgeIntv_Demographics
#                 , 
#                 family = poisson(link = "log"), 
#                 data = merged_df) %>%
#   tbl_regression(
#     exponentiate = TRUE, 
#     label = list(
#       STRESS_CT_CAT ~ "High Stress"
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
# mod2_tbl <- glm(IRREGULAR_MENSES ~ STRESS_TH_CAT 
#                 + PREG_EVER + ALCOHOL_EVER + SMOKE_LAST10_EVER + AgeIntv_Demographics
#                 ,
#                 family = poisson(link = "log"), 
#                 data = merged_df) %>%
#   tbl_regression(
#     exponentiate = TRUE, 
#     label = list(
#       STRESS_TH_CAT ~ "High Stress Severity"
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
# 
# tbl_merged %>%
#   modify_caption("**Adjusted Relative Risk Ratios, Lifetime Stress and Lifetime Irregular Menses**")

#############################################Table 4#################################################
merged_df <- merged_df %>%
  mutate(
    High_Early = if_else(EATotCT >= median(EATotCT, na.rm = TRUE), 1, 0),
    High_Adult = if_else(AdultCT >= median(AdultCT, na.rm = TRUE), 1, 0),
    High_Recent = if_else(RecTotCT >= median(RecTotCT, na.rm = TRUE), 1, 0)
  ) %>%
  mutate(StressWindow4Cat = case_when(
    High_Early == 0 & High_Adult == 0 & High_Recent == 0 ~ 0,
    High_Early == 1 & High_Adult == 0 & High_Recent == 0 ~ 1,
    High_Early == 0 & High_Adult == 1 & High_Recent == 0 ~ 2,
    High_Early == 0 & High_Adult == 0 & High_Recent == 1 ~ 3,
    High_Early == 1 & High_Adult == 1 & High_Recent == 1 ~ 4,
    TRUE ~ NA_real_
  ))

merged_df$StressWindow4Cat <- as.factor(merged_df$StressWindow4Cat)

merged_df %>%
  count(StressWindow4Cat) %>%
  mutate(percentage = n / sum(n) * 100)


# smodel1_tbl <- glm(MENSES_STARS_AS_EXPECTED ~ StressWindow4Cat + PREG_EVER + 
#                      ALCOHOL_EVER + SMOKE_LAST10_EVER + AgeIntv_Demographics, 
#                    family = poisson(link = "log"), 
#                    data = merged_df) %>% 
#   tbl_regression(
#     exponentiate = TRUE, 
#     label = list(
#       StressWindow4Cat ~ "Stress Timing", 
#       PREG_EVER ~ "Ever Pregnant",
#       ALCOHOL_EVER ~ "Ever Alcohol Use",
#       SMOKE_LAST10_EVER ~ "Ever Smoked (Last 10 Years)",
#       AgeIntv_Demographics ~ "Age"
#     ),
#     tidy_fun = partial(tidy_robust, vcov = "HC1")
#   ) %>%
#   add_n(location = "label")
# 
# smodel1_tbl %>%
#   modify_caption("**Adjusted Relative Risk Ratios, Stress Timing and Recent Regular Menses**")

# smodel1_tbl <- glm(IRREGULAR_MENSES ~ StressWindow4Cat + PREG_EVER + 
#                      ALCOHOL_EVER + SMOKE_LAST10_EVER + AgeIntv_Demographics, 
#                    family = poisson(link = "log"), 
#                    data = merged_df) %>% 
#   tbl_regression(
#     exponentiate = TRUE, 
#     label = list(
#       StressWindow4Cat ~ "Stress Timing", 
#       PREG_EVER ~ "Ever Pregnant",
#       ALCOHOL_EVER ~ "Ever Alcohol Use",
#       SMOKE_LAST10_EVER ~ "Ever Smoked (Last 10 Years)",
#       AgeIntv_Demographics ~ "Age"
#     ),
#     tidy_fun = partial(tidy_robust, vcov = "HC1")
#   ) %>%
#   add_n(location = "label")
# 
# smodel1_tbl %>%
#   modify_caption("**Adjusted Relative Risk Ratios, Stress Timing and Lifetime Regular Menses**")


########################################################Table 5####################################################
# Create new dataframes for early and late menarche
df_men_early <- merged_df %>%
  filter(men_group == "Earlier")

df_men_late <- merged_df %>%
  filter(men_group == "Later")

df_men_early$STRESS_CT_CAT <- relevel(factor(df_men_early$STRESS_CT_CAT), ref = "Low Stress")
df_men_late$STRESS_CT_CAT <- relevel(factor(df_men_late$STRESS_CT_CAT), ref = "Low Stress")

df_men_early$STRESS_TH_CAT <- relevel(factor(df_men_early$STRESS_TH_CAT), ref = "Low Severity")
df_men_late$STRESS_TH_CAT <- relevel(factor(df_men_late$STRESS_TH_CAT), ref = "Low Severity")

#early
# mod1_tbl <- glm(MENSES_STARS_AS_EXPECTED ~ STRESS_CT 
#                 + PREG_EVER + ALCOHOL_EVER + SMOKE_LAST10_EVER + AgeIntv_Demographics
#                 , 
#                 family = poisson(link = "log"), 
#                 data = df_men_early) %>%
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
# 
# mod1_tbl %>%
#   modify_caption("**Adjusted Relative Risk Ratios, Early Menarche, # of Stressful Events and Recent Regular Menses**")

# mod1_tbl <- glm(MENSES_STARS_AS_EXPECTED ~ STRESS_TH 
#                 + PREG_EVER + ALCOHOL_EVER + SMOKE_LAST10_EVER + AgeIntv_Demographics
#                 , 
#                 family = poisson(link = "log"), 
#                 data = df_men_early) %>%
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
# mod1_tbl %>%
#   modify_caption("**Adjusted Relative Risk Ratios, Early Menarche, Severity of Stress and Recent Regular Menses**")


# mod1_tbl <- glm(MENSES_STARS_AS_EXPECTED ~ STRESS_CT_CAT 
#                 + PREG_EVER + ALCOHOL_EVER + SMOKE_LAST10_EVER + AgeIntv_Demographics
#                 , 
#                 family = poisson(link = "log"), 
#                 data = df_men_early) %>%
#   tbl_regression(
#     exponentiate = TRUE, 
#     label = list(
#       STRESS_CT_CAT ~ "High Stress"
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
# mod1_tbl %>%
#   modify_caption("**Adjusted Relative Risk Ratios, Early Menarche, High Stress and Recent Regular Menses**")

# mod1_tbl <- glm(MENSES_STARS_AS_EXPECTED ~ STRESS_TH_CAT 
#                 + PREG_EVER + ALCOHOL_EVER + SMOKE_LAST10_EVER + AgeIntv_Demographics
#                 , 
#                 family = poisson(link = "log"), 
#                 data = df_men_early) %>%
#   tbl_regression(
#     exponentiate = TRUE, 
#     label = list(
#       STRESS_TH_CAT ~ "High Severity of Stress"
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
# mod1_tbl %>%
#   modify_caption("**Adjusted Relative Risk Ratios, Early Menarche, High Stress Severity and Recent Regular Menses**")
#use same code for late menarche and for lifetime irregularity

######################################################Table 6############################################
#YWS
merged_df <- merged_df %>%
  mutate(MENSES_DAYS = factor(MENSES_DAYS) %>% 
           relevel(ref = "26-32 days")) 

# smodel1_tbl <- multinom(
#   MENSES_DAYS ~ STRESS_CT + PREG_EVER + 
#     ALCOHOL_EVER + SMOKE_LAST10_EVER + AgeIntv_Demographics, 
#   data = merged_df
# ) %>% 
#   tbl_regression(
#     exponentiate = TRUE, 
#     label = list(
#       STRESS_CT ~ "# of Stressful Events", 
#       PREG_EVER ~ "Ever Pregnant",
#       ALCOHOL_EVER ~ "Ever Alcohol Use",
#       SMOKE_LAST10_EVER ~ "Ever Smoked (Last 10 Years)",
#       AgeIntv_Demographics ~ "Age"
#     )) %>% 
#   bold_labels() %>%
#   add_n(location = "label")
# 
# smodel1_tbl %>%
#   modify_caption("**Adjusted Relative Risk Ratios, Stressful Events and YWS Cycle Length**")

merged_df <- merged_df %>%
  mutate(MENSES_DAYS = factor(MENSES_DAYS) %>% 
           relevel(ref = "26-32 days")) 

# smodel1_tbl <- multinom(
#   MENSES_DAYS ~ STRESS_TH + PREG_EVER + 
#     ALCOHOL_EVER + SMOKE_LAST10_EVER + AgeIntv_Demographics, 
#   data = merged_df
# ) %>% 
#   tbl_regression(
#     exponentiate = TRUE, 
#     label = list(
#       STRESS_TH ~ "Stress severity", 
#       PREG_EVER ~ "Ever Pregnant",
#       ALCOHOL_EVER ~ "Ever Alcohol Use",
#       SMOKE_LAST10_EVER ~ "Ever Smoked (Last 10 Years)",
#       AgeIntv_Demographics ~ "Age"
#     )) %>% 
#   bold_labels() %>%
#   add_n(location = "label")
# 
# smodel1_tbl %>%
#   modify_caption("**Adjusted Relative Risk Ratios, Stress Severity and YWS Cycle Length**")

#skiptrack
merged_df <- merged_df %>% 
  mutate(
    meancycle_Clue_cat = case_when(
      meancycle_Clue < 26 ~ "<26 days",
      meancycle_Clue >= 26 & meancycle_Clue <= 32 ~ "26-32 days",
      meancycle_Clue > 32 ~ ">32 days",
      TRUE ~ NA_character_
    ), 
    meancycle_Clue_cat = factor(meancycle_Clue_cat, 
         
                                                       levels = c("26-32 days", "<26 days", ">32 days"))
  )
merged_df <- merged_df %>%
  mutate(meancycle_Clue_cat = factor(meancycle_Clue_cat) %>% 
           relevel(ref = "26-32 days")) 
  
# smodel1_tbl <- multinom(
#   meancycle_Clue_cat ~ STRESS_CT + PREG_EVER + 
#     ALCOHOL_EVER + SMOKE_LAST10_EVER + AgeIntv_Demographics, 
#   data = merged_df
# ) %>% 
#   tbl_regression(
#     exponentiate = TRUE, 
#     label = list(
#       STRESS_CT ~ "# of Stressful Events", 
#       PREG_EVER ~ "Ever Pregnant",
#       ALCOHOL_EVER ~ "Ever Alcohol Use",
#       SMOKE_LAST10_EVER ~ "Ever Smoked (Last 10 Years)",
#       AgeIntv_Demographics ~ "Age"
#     )) %>% 
#   bold_labels() %>%
#   add_n(location = "label")
# 
# smodel1_tbl %>%
#   modify_caption("**Adjusted Relative Risk Ratios, Stressful Events and Clue Cycle Length**")

# smodel1_tbl <- multinom(
#   meancycle_Clue_cat ~ STRESS_TH + PREG_EVER + 
#     ALCOHOL_EVER + SMOKE_LAST10_EVER + AgeIntv_Demographics, 
#   data = merged_df
# ) %>% 
#   tbl_regression(
#     exponentiate = TRUE, 
#     label = list(
#       STRESS_TH ~ "Stress Severity", 
#       PREG_EVER ~ "Ever Pregnant",
#       ALCOHOL_EVER ~ "Ever Alcohol Use",
#       SMOKE_LAST10_EVER ~ "Ever Smoked (Last 10 Years)",
#       AgeIntv_Demographics ~ "Age"
#     )) %>% 
#   bold_labels() %>%
#   add_n(location = "label")
# 
# smodel1_tbl %>%
#   modify_caption("**Adjusted Relative Risk Ratios, Stress Severity and Clue Cycle Length**")

#again with mean cycle_inferred
merged_df <- merged_df %>% 
  mutate(
    meancycle_inferred_cat = case_when(
      meancycle_inferred < 26 ~ "<26 days",
      meancycle_inferred >= 26 & meancycle_inferred <= 32 ~ "26-32 days",
      meancycle_inferred > 32 ~ ">32 days",
      TRUE ~ NA_character_
    ), 
    meancycle_inferred_cat = factor(meancycle_inferred_cat, 
                                
                                levels = c("26-32 days", "<26 days", ">32 days"))
  )
merged_df <- merged_df %>%
  mutate(meancycle_inferred_cat = factor(meancycle_inferred_cat) %>% 
           relevel(ref = "26-32 days")) 

# smodel1_tbl <- multinom(
#   meancycle_inferred_cat ~ STRESS_CT + PREG_EVER + 
#     ALCOHOL_EVER + SMOKE_LAST10_EVER + AgeIntv_Demographics, 
#   data = merged_df
# ) %>% 
#   tbl_regression(
#     exponentiate = TRUE, 
#     label = list(
#       STRESS_CT ~ "# of Stressful Events", 
#       PREG_EVER ~ "Ever Pregnant",
#       ALCOHOL_EVER ~ "Ever Alcohol Use",
#       SMOKE_LAST10_EVER ~ "Ever Smoked (Last 10 Years)",
#       AgeIntv_Demographics ~ "Age"
#     )) %>% 
#   bold_labels() %>%
#   add_n(location = "label")
# 
# smodel1_tbl %>%
#   modify_caption("**Adjusted Relative Risk Ratios, Stressful Events and Inferred Clue Cycle Length**")

# smodel1_tbl <- multinom(
#   meancycle_inferred_cat ~ STRESS_TH + PREG_EVER + 
#     ALCOHOL_EVER + SMOKE_LAST10_EVER + AgeIntv_Demographics, 
#   data = merged_df
# ) %>% 
#   tbl_regression(
#     exponentiate = TRUE, 
#     label = list(
#       STRESS_TH ~ "Stress Severity", 
#       PREG_EVER ~ "Ever Pregnant",
#       ALCOHOL_EVER ~ "Ever Alcohol Use",
#       SMOKE_LAST10_EVER ~ "Ever Smoked (Last 10 Years)",
#       AgeIntv_Demographics ~ "Age"
#     )) %>% 
#   bold_labels() %>%
#   add_n(location = "label")
# 
# smodel1_tbl %>%
#   modify_caption("**Adjusted Relative Risk Ratios, Stress Severity and Inferred Cycle Length**")

########################################Table 7############################################
# mod1_tbl <- glm(mean_tau ~ STRESS_CT, family = gaussian(link = "identity"), data = merged_df, na.action = na.omit) %>%
#   tbl_regression( 
#     label = list(STRESS_CT ~ "# of Stressful Events", 
#                  PREG_EVER = "Ever Pregnant",
#                  ALCOHOL_EVER = "Ever Alcohol Use",
#                  SMOKE_LAST10_EVER = "Ever Smoked (Last 10 Years)",
#                  AgeIntv_Demographics = "Age"
#     ))%>% 
#   add_n() 
# 
# print(mod1_tbl)

# mod1_tbl <- glm(mean_tau ~ STRESS_TH, family = gaussian(link = "identity"), data = merged_df, na.action = na.omit) %>%
#   tbl_regression( 
#     label = list(STRESS_TH ~ "Severity of Stress", 
#                  PREG_EVER = "Ever Pregnant",
#                  ALCOHOL_EVER = "Ever Alcohol Use",
#                  SMOKE_LAST10_EVER = "Ever Smoked (Last 10 Years)",
#                  AgeIntv_Demographics = "Age"
#     ))%>% 
#   add_n() 
# 
# print(mod1_tbl)
  