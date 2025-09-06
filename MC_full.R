setwd("//cbsp.nl/Productie/primair/bs_samen_sec1/Werk/DATA/Stage_Furkan/0-Bronbestanden")

library(data.table)
library(haven)
library(plyr)
library(MASS)
library(tidyverse)
library(readxl)
library(nFactors)
library(rlist)
library(magrittr)
library(sampling)
library(scales)
library(gridExtra)

# Select only data from certain year
year   <- 2022

# All study variables by name
var_selection <- c("OMZETPS210000", "PERSONS100000",
                   "LOONSOM110002", "RESULTS120000",
                   "OUT_PRODUCT200000")

# All auxiliary variables by name
auxiliary <- c("Netto_Omzet_minus_accijnzen",
               "Arbeidskosten",
               "Lonen",
               "Sociale_lasten",
               "Pensioenlasten",
               "Overige_personeelskosten",
               "Inkoopwaarde_omzet",
               "Uitbesteed_werk",
               "Inkoopwaarde_omzet_excl_uitbesteed",
               "Marge_omzet")

# All kerncels by code
kc_selection = c("10500",
                 "33120",
                 "47110",
                 "62000",
                 "80000",
                 "85510")

#------------------------------------------------------------------------------#

################ Obtain all files ################


# Obtain all files from all folders within working directory
file.names <- list.files(recursive=T)

# Extract all PS files from 2022
get_pattern <- paste0("(?=.*main)(?=.*", year, ")")
files_year <- grep(get_pattern, file.names, perl=T, value=T)

# Load in all the PS files
files_list <- lapply(files_year, read_sav)


#------------------------------------------------------------------------------#


################ Load in data ################


# Combine PS files and load in remaining files
PS_complete <- do.call(rbind, files_list)
BEW_full <- read_xlsx(list.files(pattern = "bewerking"))
WIA_full <- read.csv(grep("(?=.*WIA)(?=.*2.0.csv)", file.names, perl=T, value=T),
                     sep=";")
kerncel_data <- list.files(pattern = "kerncel")
KC <- read.csv(kerncel_data[grep("sbi", kerncel_data)], sep=";")


#----------------#
### CBS DATA INFO ###
# PS = Productiestatistiek, annual survey sample of companies about economic activities
# WIA = Profit declaration register, declaration of annual profits of companies to tax authorities
# KC = Information about which SBI codes belongs to which kerncel code
# BEW = For some kerncels, there is underrepresentation/nonresponse in the PS, and
#       supplementary estimate is produced. BEW tells which kerncels it concerns. 
#----------------#


#------------------------------------------------------------------------------#


################ Data wrangling ################


PS_fil <- PS_complete %>%
  
  filter(!is.na(EindGewicht)) %>%   # Filter out NA sampling weights 
  mutate_if(is.labelled, as.character) %>%  # Convert lbl class to chr
  replace(is.na(.), 0) %>%  # Missing value = nonresponse (0)
  filter(kerncelcode %in% kc_selection) %>%  # Keep selection of kerncels
  # Keep selection of relevant variables
  mutate(BE_ID = BE_Id,
         SBI = sbi,
         KERNCEL = toupper(kerncelcode),
         year = STATJAAR,
         RECHTSVORM = RechtsvormID,
         GK = gk_SBS,
         Eindgewicht = EindGewicht,
         IMPUTATIE = ImputatieGebruiken,
         DUMMY = IsDummy,
         # Combine study variables to create auxiliary variables that match
         # variables from WIA
         Netto_Omzet_minus_accijnzen = OMZETPS210000 - BEDRLST349610,
         Arbeidskosten = LOONSOM110002 + LOONSOM121200 + LOONSOM121300 +
           BEDRLST345500 + BEDRLST345600 + BEDRLST345700 + BEDRLST345900 +
           SUBSIDI140000 + INKWRDE120310 + INKWRDE120320 + INKWRDE120330,
         Lonen = LOONSOM110002 + ONTVANG120000 + SUBSIDI140000 + INKWRDE120310,
         Sociale_lasten = LOONSOM121300 + INKWRDE120320,
         Pensioenlasten = LOONSOM121200 + INKWRDE120330,
         Overige_personeelskosten = BEDRLST345500 + BEDRLST345600 +
           BEDRLST345700 + BEDRLST345900 + SUBSIDI140000,
         Inkoopwaarde_omzet = INKWRDE100000,
         Uitbesteed_werk = INKWRDE132000,
         Inkoopwaarde_omzet_excl_uitbesteed = INKWRDE100000 - INKWRDE132000,
         Marge_omzet = OMZETPS210000 - BEDRLST349610 - INKWRDE100000,
         .keep = "none") %>%
  # Keep selection of study variables
  left_join(PS_complete %>% select(BE_Id, all_of(var_selection)),
            by = join_by(BE_ID == BE_Id))



WIA_fil <- WIA_full %>%
  # WIA is one full df with all years, filter only relevant rows and get the 
  # kerncel code corresponding to the SBI code (kerncel unavailable in WIA).
  mutate(TIJDVAK = as.numeric(substr(TIJDVAK, 6, 9)),
         KERNCEL = KC$kerncelPS_1[match(WIA_full$SBI_COOR, KC$sbi2008_2)]) %>%
  filter(TIJDVAK == year,
         X2030   == "A" | X2030 == "C",
         X2031   == 100,
         X2032   == 100) %>%
  replace(is.na(.), 0) %>%   # Missing value = nonresponse (0)
  filter(KERNCEL %in% kc_selection) %>%  # Keep selection of kerncels
  # Keep selection of relevant variables
  mutate(BE_ID = BE_ID,
         SBI = SBI_COOR,
         KERNCEL = toupper(as.character(KERNCEL)),
         GK = GK_SBS,
         year = TIJDVAK,
         RECHTSVORM = RECHTSVORM,
         # Combine variables to create auxiliary variables that match those in PS
         Netto_Omzet_minus_accijnzen = X1886 / 1000,
         Arbeidskosten = (X1926 + X1930 + X1932 + X1934) / 1000,
         Lonen = X1926 / 1000,
         Sociale_lasten = X1930 / 1000,
         Pensioenlasten = X1932 / 1000,
         Overige_personeelskosten = X1934 / 1000,
         Inkoopwaarde_omzet = (X1920 + X1922) / 1000,
         Uitbesteed_werk = X1922 / 1000,
         Inkoopwaarde_omzet_excl_uitbesteed = X1920 / 1000,
         Marge_omzet = (X1886 - X1920 - X1922) / 1000,
         .keep = "none") %>%
  distinct(BE_ID, .keep_all=T) # Remove duplicate companies


# Only keep a selection of kerncels in BEW 
BEW <- BEW_full %>%
  mutate(KERNCEL = toupper(kerncel)) %>%
  filter(Verslagjaar == year,
         KERNCEL %in% kc_selection) %>%
  select(-kerncel)


# Size class (GK) is determinant of kerncels for which supplementary estimations 
# are produced. The first integer can be matched in BEW
PS_fil$GK1 <- as.numeric(substr(PS_fil$GK, 1, 1))
WIA_fil$GK1 <- as.numeric(substr(WIA_fil$GK, 1, 1))

filter_obsv <- function(df)
  # Creates a new variable ("check") that tells for which units/kerncels CBS 
  # has (not) produced supplementary estimations.
{
  df %>%
    left_join(BEW, by="KERNCEL", relationship="many-to-many") %>%
    mutate(check = if_else(GK1 >= GKMin & GK1 <= GKMax & GKBereikTypeID %in% 4,
                           1, 0, 0)) %>%
    group_by(BE_ID, KERNCEL) %>%
    summarise(check = max(check), .groups="drop") %>%
    right_join(df, by=c("BE_ID", "KERNCEL")) %>%
    mutate(check = coalesce(check, 0))
}
PS_obs <- filter_obsv(PS_fil)
WIA_obs <- filter_obsv(WIA_fil)

# Save those observations that should be filtered out of PS separately. 
PS_plus <- PS_obs %>%
  filter(IMPUTATIE == "True" | DUMMY == "True" | check == 0) %>%
  select(-GK1)

# Keep check = 1 and match units with PS based on BE_ID (delta = 1)
WIA_keep <- WIA_obs %>%
  filter(check == 1) %>%  # No supplementary estimates from CBS 
  mutate(delta = 1,       # Units that are also included in PS
         Z = 1) %>%       # Units in WIA are Z=1
  relocate(delta, .before = Netto_Omzet_minus_accijnzen) %>%
  select(-c(check))
PS_keep  <- PS_obs  %>%
  anti_join(PS_plus) %>%  # No supplementary estimates from CBS (check=1)
  # Indicator variable: if included in WIA, delta=1, if not, delta=0
  mutate(delta = ifelse(BE_ID %in% WIA_keep$BE_ID, 1, 0), 
         w_i = 0,         # Predefine column for calibration weights
         Z = 0) %>%       # Units in PS are Z=0
  relocate(c(w_i, delta), .before = Netto_Omzet_minus_accijnzen) %>%
  select(-c(check, IMPUTATIE, DUMMY))


inclusion <- function(df, PS)
  # Redefine the sampling weights such that it is equal for all units inside a
  # kerncel
{
  INCL <- PS %>%
    group_by(KERNCEL) %>%
    summarise(uni = n(), sums = sum(Eindgewicht), d_i = sums / uni, .groups = "drop")
  
  new_df <- df %>%
    left_join(INCL, by=join_by(KERNCEL)) %>%
    select(-c(uni, sums, GK1)) %>%
    filter(d_i > 0)
  
  return(new_df)
}
PS  <- inclusion(PS_keep, PS_obs)
WIA <- inclusion(WIA_keep, PS_obs)


## Checkpoint
rm(list = setdiff(ls(), c("PS",
                          "WIA",
                          "auxiliary",
                          "kc_selection",
                          "var_selection"
)))

get_all <- ls()
save(list = get_all, file = "MC2022.RData")
load("MC2022.RData")

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

############# Recreate population data with simulation scenarios ###############



setwd("F:/Documents/Rdata/MCsimulation")
source("MC_functions.R")
set.seed(66)


PS_WIA <- PS %>% 
  # Row-bind WIA-auxiliary data of PS units with delta=1 to PS
  left_join(WIA %>% select(BE_ID, all_of(auxiliary)), by = join_by(BE_ID),
            suffix = c(".PS", ".WIA")) %>% 
  # Determine if sampling weight should be rounded up/down for duplication 
  # based on a coin flip
  mutate(highlow = rbinom(n(),1,0.5),
         rounded_di = if_else(highlow==1, ceiling(d_i), floor(d_i))) 


# Duplicate observations based on the rounded sampling weight
dupe_pop <- PS_WIA[rep(row.names(PS_WIA), times = PS_WIA$rounded_di), ]
dupe_pop %<>% select(-c(highlow, Eindgewicht)) %>% 
  group_by(BE_ID) %>%
  # Set aside the original observations, because these will not be perturbed
  slice(-1) %>%
  ungroup() 

## From here:
# B-sample = WIA
# A-sample = PS

# Perturbation of individual duplicated X&Y values (this can take hours)
all_vars <- c(paste0(auxiliary, ".PS"), paste0(auxiliary, ".WIA"), var_selection)
WIA_vars <- paste0(auxiliary, ".WIA")

dupe_pop_dt1  <- as.data.table(copy(dupe_pop))
dupe_pop_dt10 <- as.data.table(copy(dupe_pop))
dupe_pop_dt25 <- as.data.table(copy(dupe_pop))

# 2 percent quantile range (base)
dupe_pop_dt1[, (all_vars) := lapply(all_vars, function(var){
  sapply(.SD[[var]], function(x) perturb_x(.SD, var, x))
}), by = KERNCEL]

# 20 percent quantile range (B10)
dupe_pop_dt10[, (WIA_vars) := lapply(WIA_vars, function(var){
  sapply(.SD[[var]], function(x) perturb_x(.SD, var, x, delta = 0.1))
}), by = KERNCEL]

# 50 percent quantile range (B25)
dupe_pop_dt25[, (WIA_vars) := lapply(WIA_vars, function(var){
  sapply(.SD[[var]], function(x) perturb_x(.SD, var, x, delta = 0.25))
}), by = KERNCEL]


# For scenario without measurement errors in B
noisy_pop0  <- as.data.frame(dupe_pop_dt1)
noisy_pop0$Netto_Omzet_minus_accijnzen.WIA <- noisy_pop0$Netto_Omzet_minus_accijnzen.PS 
noisy_pop0$Arbeidskosten.WIA <- noisy_pop0$Arbeidskosten.PS 
noisy_pop0$Lonen.WIA <- noisy_pop0$Lonen.PS 
noisy_pop0$Sociale_lasten.WIA <- noisy_pop0$Sociale_lasten.PS 
noisy_pop0$Pensioenlasten.WIA <- noisy_pop0$Pensioenlasten.PS 
noisy_pop0$Overige_personeelskosten.WIA <- noisy_pop0$Overige_personeelskosten.PS 
noisy_pop0$Inkoopwaarde_omzet.WIA <- noisy_pop0$Inkoopwaarde_omzet.PS 
noisy_pop0$Uitbesteed_werk.WIA <- noisy_pop0$Uitbesteed_werk.PS 
noisy_pop0$Inkoopwaarde_omzet_excl_uitbesteed.WIA <- noisy_pop0$Inkoopwaarde_omzet_excl_uitbesteed.PS 
noisy_pop0$Marge_omzet.WIA <- noisy_pop0$Marge_omzet.PS 

PS_WIA0 <- PS_WIA 
PS_WIA0$Netto_Omzet_minus_accijnzen.WIA <- PS_WIA0$Netto_Omzet_minus_accijnzen.PS 
PS_WIA0$Arbeidskosten.WIA <- PS_WIA0$Arbeidskosten.PS 
PS_WIA0$Lonen.WIA <- PS_WIA0$Lonen.PS 
PS_WIA0$Sociale_lasten.WIA <- PS_WIA0$Sociale_lasten.PS 
PS_WIA0$Pensioenlasten.WIA <- PS_WIA0$Pensioenlasten.PS 
PS_WIA0$Overige_personeelskosten.WIA <- PS_WIA0$Overige_personeelskosten.PS 
PS_WIA0$Inkoopwaarde_omzet.WIA <- PS_WIA0$Inkoopwaarde_omzet.PS 
PS_WIA0$Uitbesteed_werk.WIA <- PS_WIA0$Uitbesteed_werk.PS 
PS_WIA0$Inkoopwaarde_omzet_excl_uitbesteed.WIA <- PS_WIA0$Inkoopwaarde_omzet_excl_uitbesteed.PS 
PS_WIA0$Marge_omzet.WIA <- PS_WIA0$Marge_omzet.PS 

noisy_pop1  <- as.data.frame(dupe_pop_dt1)  # base
noisy_pop10 <- as.data.frame(dupe_pop_dt10) # B10
noisy_pop25 <- as.data.frame(dupe_pop_dt25) # B25


create_sampleB <- function(noisy_pop, PS_WIA)
  # Separates sample B from the population and merges original units back to the
  # data.
{
  noisy_pop %<>%  
    filter(!is.na(BE_ID)) %>%                       
    mutate(BE_ID = row_number() + 100) %>%   # Create unique ID's for duplicated units
    full_join(PS_WIA) %>%  # Join back original units
    select(-c(Eindgewicht, highlow, d_i)) %>% 
    rename(d_i = rounded_di) # Use rounded sampling weights as new sampling weights
  
  sampleB <- noisy_pop %>% 
    filter(delta == 1) %>% 
    select(BE_ID, KERNCEL, RECHTSVORM, GK, d_i, delta, ends_with(".WIA"), 
           -ends_with(".PS"), -all_of(var_selection)) %>% 
    rename_with(~ str_remove(., ".WIA"), ends_with(".WIA")) %>% 
    mutate(Z = 1)
  
  return(list(pop = noisy_pop,
              B   = sampleB))
}

current_pop <- create_sampleB(noisy_pop1, PS_WIA)$pop  # Full population with 0.01 perturbation
B_0  <- create_sampleB(noisy_pop0, PS_WIA0)$B          # B without measurement errors
B_1  <- create_sampleB(noisy_pop1, PS_WIA)$B           # B with 0.01 perturbation
B_10 <- create_sampleB(noisy_pop10, PS_WIA)$B          # B with 0.1 perturbation
B_25 <- create_sampleB(noisy_pop25, PS_WIA)$B          # B with 0.25 perturbation

# Remove additional B-auxiliary variables from population 
full_pop <- current_pop %>% 
  select(-ends_with(".WIA")) %>% 
  rename_with(~ str_remove(., ".PS"), ends_with(".PS")) 


# Obtain true population totals of study variables in list and matrix form
T_y_list <- full_pop %>% 
  group_by(KERNCEL) %>%
  group_map(~ colSums(.x %>% select(all_of(var_selection)))) %>% 
  set_names(sort(as.numeric(kc_selection))) %>% 
  map( ~ list(T_y = rbind(.x)))

T_y_mat <- full_pop %>% 
  group_by(KERNCEL) %>%
  group_map(~ colSums(.x %>% select(all_of(var_selection)))) %>% 
  bind_rows() %>% 
  as.matrix() %>% 
  set_rownames(sort(as.numeric(kc_selection)))


## Checkpoint
save.image(file = "MC_prep.RData")
load("MC_prep.RData")

#------------------------------------------------------------------------------#

############# Monte Carlo Simulation ###############

setwd("F:/Documents/Rdata/MCsimulation")
source("MC_functions.R")


# Total Monte Carlo iterations
N <- 1000

start_time <- Sys.time()  # Keep track of time

# Scenario 1: Baseline/realistic scenario, perturbation tau=0.01
hat_T_y_base  <- monte_carlo(full_pop, B_1,  N, T_y_mat, T_y_list, auxiliary, var_selection)

# Scenario 2: Baseline scenario, but values from A=B (no measurement error in B)
hat_T_y_B0  <- monte_carlo(full_pop, B_0,  N, T_y_mat, T_y_list, auxiliary, var_selection)

# Scenario 3: Baseline scenario, but measurement errors in B are increased (tau=0.1, tau=0.25)
hat_T_y_B10 <- monte_carlo(full_pop, B_10, N, T_y_mat, T_y_list, auxiliary, var_selection)
hat_T_y_B25 <- monte_carlo(full_pop, B_25, N, T_y_mat, T_y_list, auxiliary, var_selection)

# Scenario 4: Shrinking sample A by:
## 80%
hat_T_y_Aprop0.8 <- monte_carlo(full_pop, B_1,  N, T_y_mat, T_y_list, auxiliary, var_selection, A_prop = 0.8, A_slicing = T) 
## 50%
hat_T_y_Aprop0.5 <- monte_carlo(full_pop, B_1,  N, T_y_mat, T_y_list, auxiliary, var_selection, A_prop = 0.5, A_slicing = T) 

# Scenario 5: Shrinking sample B by:
## 80%
hat_T_y_Bprop0.8 <- monte_carlo(full_pop, B_1,  N, T_y_mat, T_y_list, auxiliary, var_selection, B_prop = 0.8, B_slicing = T) 
## 50%
hat_T_y_Bprop0.5 <- monte_carlo(full_pop, B_1,  N, T_y_mat, T_y_list, auxiliary, var_selection, B_prop = 0.5, B_slicing = T) 

end_time <- Sys.time()
end_time - start_time  # What is total runtime? (~ 4 hours)



## Checkpoint
save(list = c("hat_T_y_base",
              "hat_T_y_B0",
              "hat_T_y_B10",
              "hat_T_y_B25",
              "hat_T_y_Aprop0.8",
              "hat_T_y_Aprop0.5",
              "hat_T_y_Bprop0.8", 
              "hat_T_y_Bprop0.5"), file = "RESULTS2022.RData")
load("RESULTS2022.RData")



#------------------------------------------------------------------------------#

############# Post-Hoc Analysis ###############

kc_selection_small <- c("10500", "80000")
full_pop_small <- full_pop %>% filter(KERNCEL %in% kc_selection_small)
B_1_small <- B_1 %>% filter(KERNCEL %in% kc_selection_small)

  
# Obtain true population totals of study variables in list and matrix for delta=1
T_y_list1 <- full_pop_small %>%
  filter(delta == 1) %>% 
  group_by(KERNCEL) %>%
  group_map(~ colSums(.x %>% select(all_of(var_selection)))) %>% 
  set_names(sort(as.numeric(kc_selection_small))) %>% 
  map( ~ list(T_y = rbind(.x)))

T_y_mat1 <- full_pop_small %>% 
  filter(delta == 1) %>% 
  group_by(KERNCEL) %>%
  group_map(~ colSums(.x %>% select(all_of(var_selection)))) %>% 
  bind_rows() %>% 
  as.matrix() %>% 
  set_rownames(sort(as.numeric(kc_selection_small)))
  
# Obtain true population totals of study variables in list and matrix for delta=0
T_y_list0 <- full_pop_small %>% 
  filter(delta == 0) %>% 
  group_by(KERNCEL) %>%
  group_map(~ colSums(.x %>% select(all_of(var_selection)))) %>% 
  set_names(sort(as.numeric(kc_selection_small))) %>% 
  map( ~ list(T_y = rbind(.x)))

T_y_mat0 <- full_pop_small %>% 
  filter(delta == 0) %>% 
  group_by(KERNCEL) %>%
  group_map(~ colSums(.x %>% select(all_of(var_selection)))) %>% 
  bind_rows() %>% 
  as.matrix() %>% 
  set_rownames(sort(as.numeric(kc_selection_small)))


setwd("~/Rdata/MCsimulation/posthoc")
source("delta0_1functions.R")

set.seed(66)

# Total Monte Carlo iterations
N <- 1000

# Get estimates
estimates_delta1  <- monte_carlo(full_pop_small, B_1_small,  N, T_y_mat1, T_y_list1, auxiliary, var_selection, Delta=1)
estimates_delta0  <- monte_carlo(full_pop_small, B_1_small,  N, T_y_mat0, T_y_list0, auxiliary, var_selection, Delta=0)


# Get calibration weights
source("delta0_1weights.R")

full_pop10500 <- full_pop_small %>% filter(KERNCEL == "10500") 
full_pop80000 <- full_pop_small %>% filter(KERNCEL == "80000") 
B_1_10500 <- B_1_small %>% filter(KERNCEL == "10500")
B_1_80000 <- B_1_small %>% filter(KERNCEL == "80000")

weights_delta1_10500 <- monte_carlo(full_pop10500, B_1_10500,  N, auxiliary, Delta=1)
weights_delta1_80000 <- monte_carlo(full_pop80000, B_1_80000,  N, auxiliary, Delta=1)
weights_delta0_10500 <- monte_carlo(full_pop10500, B_1_10500,  N, auxiliary, Delta=0)
weights_delta0_80000 <- monte_carlo(full_pop80000, B_1_80000,  N, auxiliary, Delta=0)
