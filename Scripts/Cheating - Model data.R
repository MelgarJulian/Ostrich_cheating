##**************************##
# Frequency-dependent cheating in a cooperative breeding bird - Model-ready data----
##**************************##
#
##*****************##
## Load packages----
##*****************##
#
pacman::p_load("magrittr", "dplyr", "doBy", "data.table", "reshape", "tidyr", "tidyverse", "date")
#
#
##*************##
## Load data----
##*************##
#
load("R_Output_files/Cheating_data.RData")
#
#
##****************************##
## Guide to data frame names----
##****************************##
#
#Data frames are named df_<topic>, using the same topic names as the models in
#Cheating_MCMCglmmModels.R, so that a frame and the model fitted to it can be
#read off each other. Neither carries a number, so models and sections can be
#reordered without renaming anything.
#
#  Data frame                  Model fitted to it        Contents
#  df_strategy_counts          m_strategy_counts         Number of cheats and cooperators per group
#  df_cheat_prop               m_cheat_prop              Proportion of cheats per group
#  df_cheat_prob               m_cheat_prob              Probability of being a cheat, females seen in >1 year
#  df_eggs                     m_eggs                    Reproductive investment of individual females
#  df_eggs_excl_1coop          m_eggs_excl_1coop         As df_eggs, without groups with one cooperator
#  df_mass                     m_mass                    Body mass at the start and end of the season
#  df_eggs_within_female       m_eggs_within_female      Females recorded both cheating and cooperating
#  df_incubation_individual    m_incubation_individual   Individual incubation effort
#  df_incubation_group         m_incubation_group        Group-level incubation effort
#  df_aggression               m_aggression              Aggression received, by breeding strategy
#  df_laying_window            m_laying_window           Probability of laying within a two-day window
#  df_laying_group             m_laying_group            Group-level egg laying across time
#  df_incubation_time          m_incubation_time         Group-level incubation effort across time
#  df_activity_allocation      m_activity_nocheat,       Daily split between laying and incubation
#                              m_activity_cheat
#  df_nest_failure             m_nest_failure            Whether a nest hatched any eggs (built in the models script)
#  df_hatch_success            m_hatch_success           Proportion hatching, nests that hatched at least one
#
#The remaining frames are stepping stones rather than model inputs:
#
#  df_females             base female-level frame, filtered to cheats and cooperators
#  df_females_n           counts of each breeding strategy per group
#  df_females_breeding    base frame for the egg models
#  df_cooperators         base frame for individual incubation
#  df_groups              base group-level frame
#  df_hatching            group hatching outcomes; parent of df_nest_failure
#  df_activity            whether a nest was active on a given day; parent of df_activity_allocation
#  df_eggs_all            df_eggs before excluding groups with 0 or 6 cooperators
#  df_groups_all          df_groups before excluding groups with 0 or 6 cooperators
#
#df_eggs_all and df_groups_all carry an InModel column marking which rows the
#models keep, for figures that need to show the excluded groups.
#
#
#******#
### Frequency of cooperators and cheats----
#******#
#
#Subset data to only include females
df_females <- individuals2012.2018[individuals2012.2018$sex == "F",]
#
#Subset data to only include cheats and cooperators
table(df_females$BS, exclude = NULL)
df_females <- df_females[df_females$BS != "Other",]
#
#Subset data to exclude groups with only one female
df_females <- df_females[df_females$nofemales != 1,]
#
#Subset data to only include camps where eggs were laid during the period with natural incubation
#and eggs were genotyped during the period with natural incubation and 4 weeks before
#Match data
df_females$totGenoeggsCampInc4 <- camps2012.2018$totGenoeggsCampInc4[match(df_females$yearcamp, camps2012.2018$yearcamp)]
df_females$no.eggsCampInc <- camps2012.2018$no.eggsCampInc[match(df_females$yearcamp, camps2012.2018$yearcamp)]
df_females <- df_females[which(df_females$totGenoeggsCampInc4 != 0 & df_females$no.eggsCampInc != 0),]
#
#Exclude individuals in camps with no nests
#Match data
df_females$nonests <- camps2012.2018$nonests[match(df_females$yearcamp, camps2012.2018$yearcamp)]
table(df_females$nonests, exclude = NULL)#No camps with no nests
#
#Exclude individuals in camps with no female incubation
df_females <- df_females[df_females$femSitdurCamp != 0,]
#
#Subset data to only include relevant columns
names(df_females)
df_females <- df_females[c("year", "camp", "yearcamp", "nomales", "nofemales", "noFActive", "nonests", "BS", "age", "id")]
#
#Match data on number of nests
df_females$nonests = camps2012.2018$nonests[match(df_females$yearcamp, camps2012.2018$nonests)]
#
#
####df_strategy_counts: Number of cheats and cooperators per group----
#****************************************************************************
#
#Compute the number of of individuals of each breeding strategy present in each group
df_females_n <- df_females %>% dplyr::select(yearcamp, nomales, nofemales, BS)%>%
  dplyr::group_by(yearcamp, BS, nomales, nofemales) %>% dplyr::summarize(n=dplyr::n())
df_females_n$yrcamp_cc <- paste(df_females_n$yearcamp,df_females_n$BS,sep="_")
#
#Create a list of camps with a row for each breeding strategy
df_strategy_counts <- data.frame(yearcamp = rep(unique(df_females$yearcamp),2),
                           BS = c(rep("Cheat", length(unique(df_females$yearcamp))),
                                  rep("Cooperator",length(unique(df_females$yearcamp)))))
#Create a unique id for every camp and breeding strategy
df_strategy_counts$yrcamp_cc <- paste(df_strategy_counts$yearcamp,df_strategy_counts$BS,sep="_")
#Match relevant information
df_strategy_counts$year <- camps2012.2018$year[match(df_strategy_counts$yearcamp,camps2012.2018$yearcamp)]
df_strategy_counts$camp <- camps2012.2018$camp[match(df_strategy_counts$yearcamp,camps2012.2018$yearcamp)]
df_strategy_counts$nomales <- df_females_n$nomales[match(df_strategy_counts$yearcamp,df_females_n$yearcamp)]
df_strategy_counts$nofemales <- df_females_n$nofemales[match(df_strategy_counts$yearcamp,df_females_n$yearcamp)]
df_strategy_counts$noFCoops <- individuals2012.2018$noFCoopsCamp[match(df_strategy_counts$yearcamp,individuals2012.2018$yearcamp)]
df_strategy_counts$noFCheats <- individuals2012.2018$noFCheatsCamp[match(df_strategy_counts$yearcamp,individuals2012.2018$yearcamp)]
df_strategy_counts$noFActive <- individuals2012.2018$noFActive[match(df_strategy_counts$yearcamp,individuals2012.2018$yearcamp)]
#Number of nests: a few groups had a second nest, so this is fitted as a covariate
df_strategy_counts$nonests <- camps2012.2018$nonests[match(df_strategy_counts$yearcamp,camps2012.2018$yearcamp)]
#Z-transform number of males, females and nests
df_strategy_counts$Znomales <- scale(df_strategy_counts$nomales)
df_strategy_counts$Znofemales <- scale(df_strategy_counts$nofemales)
df_strategy_counts$Znonests <- scale(df_strategy_counts$nonests)
#Compute quadratics
df_strategy_counts$Znomales2 <- df_strategy_counts$Znomales * df_strategy_counts$Znomales
df_strategy_counts$Znofemales2 <- df_strategy_counts$Znofemales * df_strategy_counts$Znofemales
#Match number of each breeding strategy per group
df_strategy_counts$n <- df_females_n$n[match(df_strategy_counts$yrcamp_cc, df_females_n$yrcamp_cc)]
df_strategy_counts$n[is.na(df_strategy_counts$n)] <- 0
#
#
####df_cheat_prop: Proportion of cheats per group----
#************************************************************************
#
#Subset data to only include camps where eggs were laid during the period with natural incubation
#and eggs were genotyped during the period with natural incubation and 4 weeks before
df_cheat_prop <- camps2012.2018[which(camps2012.2018$totGenoeggsCampInc4 != 0 & camps2012.2018$no.eggsCampInc != 0),]
#
#Subset data to exclude groups with only one female
df_cheat_prop <- df_cheat_prop[df_cheat_prop$nofemales != 1,]
#
#Exclude individuals in camps with no female incubation
df_cheat_prop <- df_cheat_prop[df_cheat_prop$Ftotsit != 0,]
#
#Subset the data to include only columns of interest
names(df_cheat_prop)
df_cheat_prop <- df_cheat_prop[c("yearcamp", "camp",  "year", "noFCoopsCamp", "noFCheatsCamp",
                           "noFOtherCamp", "nomales", "nofemales", "FCheatinCamp", "noFActive", "prop.Cheats", "nonests")]
#
#Z-transform number of cooperators, number of other females and number of nests
df_cheat_prop$ZnoFCoopsCamp <- scale(df_cheat_prop$noFCoopsCamp)
df_cheat_prop$ZnoFOtherCamp <- scale(df_cheat_prop$noFOtherCamp)
df_cheat_prop$Znonests <- scale(df_cheat_prop$nonests)
#Compute quadratics
df_cheat_prop$ZnoFCoopsCamp2 <- df_cheat_prop$ZnoFCoopsCamp * df_cheat_prop$ZnoFCoopsCamp
df_cheat_prop$ZnoFOtherCamp2 <- df_cheat_prop$ZnoFOtherCamp * df_cheat_prop$ZnoFOtherCamp
#
#
####df_cheat_prob: Probability of being a cheat, for females seen in more than one year----
#************************************************************************************
#
#Subset data to only include females
df_cheat_prob <- individuals2012.2018[individuals2012.2018$sex == "F",]
#
#Subset data to exclude groups with only one female
df_cheat_prob <- df_cheat_prob[df_cheat_prob$nofemales != 1,]
#
#Subset data to only include camps where eggs were laid during the period with natural incubation
#and eggs were genotyped during the period with natural incubation and 4 weeks before
#Match data
df_cheat_prob$totGenoeggsCampInc4 <- camps2012.2018$totGenoeggsCampInc4[match(df_cheat_prob$yearcamp, camps2012.2018$yearcamp)]
df_cheat_prob$no.eggsCampInc <- camps2012.2018$no.eggsCampInc[match(df_cheat_prob$yearcamp, camps2012.2018$yearcamp)]
df_cheat_prob <- df_cheat_prob[which(df_cheat_prob$totGenoeggsCampInc4 != 0 & df_cheat_prob$no.eggsCampInc != 0),]
#
#Exclude individuals in camps with no female incubation
df_cheat_prob <- df_cheat_prob[df_cheat_prob$femSitdurCamp != 0,]
#
#Subset data to only include individuals that appear in more than one year
length(unique(df_cheat_prob$id))#There are 159 females
length(unique(df_cheat_prob$id[duplicated(df_cheat_prob$id)]))#108 females appear in more than one year
df_cheat_prob <- df_cheat_prob[df_cheat_prob$id %in% df_cheat_prob$id[duplicated(df_cheat_prob$id)],]
#
#Subset data to only include relevant columns
names(df_cheat_prob)
df_cheat_prob <- df_cheat_prob[c("year", "camp", "yearcamp", "nomales", "nofemales", "BS", "age", "id")]
#
#Compute probability of cheating
df_cheat_prob$prob.Cheat <- ifelse(df_cheat_prob$BS == "Cheat", 1, 0)
#
#Add info on number of cheats, cooperators and other
df_cheat_prob = df_cheat_prob |> mutate(cooperators=camps2012.2018$noFCoopsCamp[match(yearcamp,camps2012.2018$yearcamp)],
                                    cheats=camps2012.2018$noFCheatsCamp[match(yearcamp,camps2012.2018$yearcamp)],
                                    other=camps2012.2018$noFOtherCamp[match(yearcamp,camps2012.2018$yearcamp)]) 
#
#Add info on mass
df_cheat_prob = df_cheat_prob |> mutate(year_id=paste(year,id,sep="_"),
                                    massin=mass$massin[match(year_id,mass$year_id)],
                                    massout=mass$massout[match(year_id,mass$year_id)],
                                    massloss=massout-massin)
#
#Filter out individuals with missing mass data, and compute focal number of cooperators
df_cheat_prob = df_cheat_prob |> dplyr::filter(!is.na(massin),!is.na(massout)) |> group_by(id,BS) |> 
  mutate(year=factor(year),
         massin=massin,
         massout=massout,
         massloss=massloss,
         cooperators=cooperators,
         cheats=cheats,
         cooperators_focal=ifelse(BS == "Cooperator",cooperators-1,cooperators),
         other_focal=ifelse(BS == "Other",other-1,other))
#
#Z-transform numerical fixed effects
df_cheat_prob$Zage <- scale(df_cheat_prob$age)
df_cheat_prob$Zmassin <- scale(df_cheat_prob$massin)
df_cheat_prob$Zmassout <- scale(df_cheat_prob$massout)
df_cheat_prob$Zcooperators_focal <- scale(df_cheat_prob$cooperators_focal)
df_cheat_prob$Zother_focal <- scale(df_cheat_prob$other_focal)

#Compute quadratics
df_cheat_prob$Zage2 <- df_cheat_prob$Zage * df_cheat_prob$Zage
df_cheat_prob$Zmassin2 <- df_cheat_prob$Zmassin * df_cheat_prob$Zmassin
df_cheat_prob$Zmassout2 <- df_cheat_prob$Zmassout * df_cheat_prob$Zmassout
df_cheat_prob$Zcooperators_focal2 <- df_cheat_prob$Zcooperators_focal * df_cheat_prob$Zcooperators_focal
df_cheat_prob$Zother_focal2 <- df_cheat_prob$Zother_focal * df_cheat_prob$Zother_focal
#
#Save data as dataframe
df_cheat_prob <- as.data.frame(df_cheat_prob)
#
#
#******#
### Consequences of cheating and cooperation for the individual----
#******#
#
#Subset data to only include females
df_females_breeding <- individuals2012.2018[individuals2012.2018$sex == "F",]
#
#Subset data to only include camps where eggs were laid during the period with natural incubation
#and eggs were genotyped during the period with natural incubation and 4 weeks before
#Match data
df_females_breeding$totGenoeggsCampInc4 <- camps2012.2018$totGenoeggsCampInc4[match(df_females_breeding$yearcamp, camps2012.2018$yearcamp)]
df_females_breeding$no.eggsCampInc <- camps2012.2018$no.eggsCampInc[match(df_females_breeding$yearcamp, camps2012.2018$yearcamp)]
df_females_breeding <- df_females_breeding[which(df_females_breeding$totGenoeggsCampInc4 != 0 & df_females_breeding$no.eggsCampInc != 0),]
#
#Subset data to only include Cooperators and cheats
df_females_breeding <- df_females_breeding[which(df_females_breeding$BS == "Cooperator" |
                                                 df_females_breeding$BS == "Cheat"),]
#
#Subset data to exclude groups with only one female
df_females_breeding <- df_females_breeding[df_females_breeding$nofemales != 1,]
#
#Exclude individuals in camps with no female incubation
df_females_breeding <- df_females_breeding[df_females_breeding$femSitdurCamp != 0,]
#
#Subset data to only include relevant columns
names(df_females_breeding)
df_females_breeding <- df_females_breeding[c("year", "camp", "yearcamp", "nomales", "nofemales", "noFCoopsCamp",
                           "noFCheatsCamp", "noFOtherCamp", "id", "age", "RSInc", "RSInc4",
                           "HatchIndInc", "relativeRSInc4", "BS", "noFActive", "prop.Cheats")]
#
#Create a column with presence/absence of cheats
df_females_breeding$FCheatinCamp <- ifelse(df_females_breeding$noFCheatsCamp == 0, "No", "Yes")
#
#Create a column with information about whether Cooperators have a cheat in their group or not
df_females_breeding$BS_2 <- ifelse(df_females_breeding$BS == "Cooperator" &
                                    df_females_breeding$FCheatinCamp == "Yes", "Coop_w_Cheat",
                                  ifelse(df_females_breeding$BS == "Cooperator" &
                                           df_females_breeding$FCheatinCamp == "No", "Coop_w/o_Cheat",
                                         df_females_breeding$BS))
#Set the level order explicitly. MCMCglmm builds the factor by sorting, and the sort order of
#"Coop_w/o_Cheat" and "Coop_w_Cheat" depends on the collation locale, which would change the order of
#the fixed effects without changing the estimates. The Extended Data tables name fixed effects by
#position, so an unpinned order would silently swap the labels of the two cooperator categories.
df_females_breeding$BS_2 <- factor(df_females_breeding$BS_2, levels = c("Cheat", "Coop_w/o_Cheat", "Coop_w_Cheat"))
#
#
####df_mass: Body mass at the start and end of the season----
#************************************************************
#
#Add info on mass
df_mass = df_females_breeding |> mutate(year_id=paste(year,id,sep="_"),
                                        massin=mass$massin[match(year_id,mass$year_id)],
                                        massout=mass$massout[match(year_id,mass$year_id)],
                                        massloss=massout-massin)
#Check new columns
table(df_mass$massin, exclude = NULL)#There are three NAs. Remove.
df_mass <- df_mass[which(!is.na(df_mass$massin)),]
table(df_mass$massin, exclude = NULL)#Good!
table(df_mass$massout, exclude = NULL)#There is one NA. Remove.
df_mass <- df_mass[which(!is.na(df_mass$massout)),]
table(df_mass$massout, exclude = NULL)#Good!
table(df_mass$massloss, exclude = NULL)#No NAs. Good!
#
#Compute z-transformed initial mass
df_mass$Zmassin <- scale(df_mass$massin)
#Compute quadratic of z-transformed reproductive investment
df_mass$Zmassin2 <- df_mass$Zmassin * df_mass$Zmassin
#
#Compute z-transformed number of cooperators
df_mass$ZnoFCoops <- scale(df_mass$noFCoopsCamp)
#Compute quadratic of z-transformed reproductive investment
df_mass$ZnoFCoops2 <- df_mass$ZnoFCoops * df_mass$ZnoFCoops
#
#Compute z-transformed age
df_mass$Zage <- scale(df_mass$age)
#Compute quadratic of z-transformed reproductive investment
df_mass$Zage2 <- df_mass$Zage * df_mass$Zage
#
#Compute z-transformed number of other females
df_mass$ZnoFOther <- scale(df_mass$noFOtherCamp)
#Compute quadratics of Zno.FOther
df_mass$ZnoFOther2 <- df_mass$ZnoFOther * df_mass$ZnoFOther
#
#Convert to data frame
df_mass <- as.data.frame(df_mass)
#
#
####df_eggs and df_eggs_all: Reproductive investment of individual females----
#****************************************************************************
#
#Create a data frame to use when drawing the figure
df_eggs_all <- df_females_breeding
#Create a column with info of whether a given row is included in the model data
df_eggs_all$InModel <- ifelse(df_eggs_all$noFCoopsCamp == 6 | df_eggs_all$noFCoopsCamp == 0, "No", "Yes")
#
#Subset data to exclude individuals in groups with no, or six female cooperators
#(no comparisons between different breeding strategies are possible in these groups since all cheat or all cooperate)
df_eggs <- df_females_breeding[df_females_breeding$noFCoopsCamp != 6 & df_females_breeding$noFCoopsCamp != 0,]
#
#Match proportion of broken eggs
df_eggs$prop.BrokenInc <- camps2012.2018$prop.BrokenInc[match(df_eggs$yearcamp, camps2012.2018$yearcamp)]
#
#Format variables
str(df_eggs)
#Format year as factor
df_eggs$year <- factor(df_eggs$year)
#Format id as factor
df_eggs$id <- factor(df_eggs$id)
#Compute Z-transformed numerical variables
df_eggs$Znomales <- scale(df_eggs$nomales)
df_eggs$Znofemales <- scale(df_eggs$nofemales)
df_eggs$ZnoFOther <- scale(df_eggs$noFOtherCamp)
df_eggs$ZnoFActive <- scale(df_eggs$noFActive)
df_eggs$ZnoFCoops <- scale(df_eggs$noFCoopsCamp)
#Compute logit values for proportions
df_eggs$logit_prop.BrokenInc <- car::logit(df_eggs$prop.BrokenInc)
#
#Compute quadratics
df_eggs$Znomales2 <- df_eggs$Znomales * df_eggs$Znomales
df_eggs$Znofemales2 <- df_eggs$Znofemales * df_eggs$Znofemales
df_eggs$ZnoFOther2 <- df_eggs$ZnoFOther * df_eggs$ZnoFOther
df_eggs$ZnoFActive2 <- df_eggs$ZnoFActive * df_eggs$ZnoFActive
df_eggs$ZnoFCoops2 <- df_eggs$ZnoFCoops * df_eggs$ZnoFCoops
#
#Save data as dataframe
df_eggs <- as.data.frame(df_eggs)
#
#m_eggs_excl_1coop asks whether the result holds once groups with a single
#cooperator are dropped
df_eggs_excl_1coop <- df_eggs[df_eggs$noFCoopsCamp != 1,]
#
#
####df_eggs_within_female: Females recorded both cheating and cooperating----
#******************************************************************************************************
#
#Subset data to only include females
df_eggs_within_female <- df_females_breeding
#
#Subset data to only include relevant columns
names(df_eggs_within_female)
df_eggs_within_female <- df_eggs_within_female[c("year", "camp", "yearcamp", "nomales", "nofemales", "BS", "age", "id","RSInc", "RSInc4", "relativeRSInc4")]
#
#Compute probability of cheating each year
df_eggs_within_female$prob.Cheat <- ifelse(df_eggs_within_female$BS == "Cheat", 1, 0)
#
#Find out whether an individual is a cheat, a cooperator or an opportunist
df_eggs_within_female <- df_eggs_within_female %>% group_by(id) %>% mutate(BS_id = case_when(all(prob.Cheat == 1) ~ "Cheat",
                                                                         all(prob.Cheat == 0) ~ "Cooperator",
                                                                         any(prob.Cheat == 0) & any(prob.Cheat == 1) ~ "Opportunist"))
df_eggs_within_female <- as.data.frame(df_eggs_within_female)
#
#Subset data to only include Opportunists
df_eggs_within_female <- df_eggs_within_female[df_eggs_within_female$BS_id == "Opportunist",]
#
#Re-format variables
str(df_eggs_within_female)
#Format year as factor
df_eggs_within_female$year <- factor(df_eggs_within_female$year)
#Format id as factor
df_eggs_within_female$id <- factor(df_eggs_within_female$id)
#
#Compute Z-transformed numerical variables
df_eggs_within_female$Znomales <- scale(df_eggs_within_female$nomales)
df_eggs_within_female$Znofemales <- scale(df_eggs_within_female$nofemales)
df_eggs_within_female$Zage <- scale(df_eggs_within_female$age)
#Compute quadratics
df_eggs_within_female$Znomales2 <- df_eggs_within_female$Znomales * df_eggs_within_female$Znomales
df_eggs_within_female$Znofemales2 <- df_eggs_within_female$Znofemales * df_eggs_within_female$Znofemales
df_eggs_within_female$Zage2 <- df_eggs_within_female$Zage * df_eggs_within_female$Zage
#
#
#******#
### Investment in cooperation----
#******#
#
#Subset data to only include females
df_cooperators <- individuals2012.2018[individuals2012.2018$sex == "F",]
#
#Subset data to only include cooperators
df_cooperators <- df_cooperators[df_cooperators$BS == "Cooperator",]
#
#Subset data to only include camps where eggs were laid during the period with natural incubation
#and eggs were genotyped during the period with natural incubation and 4 weeks before
#Match data
df_cooperators$totGenoeggsCampInc4 <- camps2012.2018$totGenoeggsCampInc4[match(df_cooperators$yearcamp, camps2012.2018$yearcamp)]
df_cooperators$no.eggsCampInc <- camps2012.2018$no.eggsCampInc[match(df_cooperators$yearcamp, camps2012.2018$yearcamp)]
#Subset data
df_cooperators <- df_cooperators[which(df_cooperators$totGenoeggsCampInc4 != 0 & df_cooperators$no.eggsCampInc != 0),]
#
#Subset data to exclude groups with only one female
df_cooperators <- df_cooperators[df_cooperators$nofemales != 1,]
#
#Subset data to only include relevant columns
names(df_cooperators)
df_cooperators <- df_cooperators[c("year", "camp", "yearcamp", "BS", "nomales", "nofemales", "noFActive", "id",
                           "yearID", "noFCheatsCamp", "noFCoopsCamp", "noFOtherCamp", "prop.Cheats", 
                           "RSInc", "RSInc4", "relativeRSInc4", "sitdurbird", "totobsdurYear")]
#
#Check the frequency of individuals in camps with x number of cheats, cooperators and other
table(df_cooperators$noFCheatsCamp, exclude = NULL)#Have an extra column with presence/absence of cheats
df_cooperators$FCheatinCamp <- ifelse(df_cooperators$noFCheatsCamp == 0, "No", "Yes")
table(df_cooperators$FCheatinCamp)
table(df_cooperators$noFCoopsCamp)#Keep as it is
table(df_cooperators$noFOtherCamp)#Have an extra column with presence/absence
df_cooperators$FOtherinCamp <- ifelse(df_cooperators$noFOtherCamp == 0, "No", "Yes")
table(df_cooperators$FOtherinCamp)
#
#Compute incubation effort as a proportion of total observation time
df_cooperators$Incubation_effort <- df_cooperators$sitdurbird / df_cooperators$totobsdurYear
#
#Create a column where number of female cooperators is merged into two categories
df_cooperators$noFCoops_cat2 <- ifelse(df_cooperators$noFCoopsCamp < 3, "Small",
                                         ifelse(df_cooperators$noFCoopsCamp == 6, NA, "Big"))
table(df_cooperators$noFCoops_cat2)
#
#Subset data to exclude individuals in groups with six cooperators (no cheats are present in these groups,
#so no comparison is possible between groups with and without cheats)
df_cooperators <- df_cooperators[df_cooperators$noFCoopsCamp != 6,]
#
#
####df_incubation_individual: Individual incubation effort----
#********************************************************************************************************************
#
df_incubation_individual <- df_cooperators
#
#Match proportion of broken eggs
df_incubation_individual$prop.BrokenInc <- camps2012.2018$prop.BrokenInc[match(df_incubation_individual$yearcamp, camps2012.2018$yearcamp)]
#
#Format variables
str(df_incubation_individual)
#Format year as factor
df_incubation_individual$year <- factor(df_incubation_individual$year)
#Format id as factor
df_incubation_individual$id <- factor(df_incubation_individual$id)
#Compute Z-transformed numerical variables
df_incubation_individual$Znomales <- scale(df_incubation_individual$nomales)
df_incubation_individual$Znofemales <- scale(df_incubation_individual$nofemales)
df_incubation_individual$ZnoFCoops <- scale(df_incubation_individual$noFCoopsCamp)
df_incubation_individual$ZnoFCheats <- scale(df_incubation_individual$noFCheatsCamp)
df_incubation_individual$ZnoFOther <- scale(df_incubation_individual$noFOtherCamp)
df_incubation_individual$ZrelativeRSInc4 <- scale(df_incubation_individual$relativeRSInc4)
#Compute logit values for proportions
df_incubation_individual$logit_prop.BrokenInc <- car::logit(df_incubation_individual$prop.BrokenInc)
#
#Compute quadratics
df_incubation_individual$Znomales2 <- df_incubation_individual$Znomales * df_incubation_individual$Znomales
df_incubation_individual$Znofemales2 <- df_incubation_individual$Znofemales * df_incubation_individual$Znofemales
df_incubation_individual$ZnoFCoops2 <- df_incubation_individual$ZnoFCoops * df_incubation_individual$ZnoFCoops
df_incubation_individual$ZnoFCheats2 <- df_incubation_individual$ZnoFCheats *  df_incubation_individual$ZnoFCheats
df_incubation_individual$ZnoFOther2 <- df_incubation_individual$ZnoFOther *  df_incubation_individual$ZnoFOther
df_incubation_individual$ZrelativeRSInc42 <- df_incubation_individual$ZrelativeRSInc4 *  df_incubation_individual$ZrelativeRSInc4
#
#Format data as dataframe
df_incubation_individual <- as.data.frame(df_incubation_individual)
#
#
#******#
### Hatching success----
#******#
#
#Hatching success in relation to group composition and group size
#****************************************************************
#
#Subset data to only include camps where eggs were laid during the period with natural incubation
#and eggs were genotyped during the period with natural incubation and 4 weeks before
df_groups <- camps2012.2018[which(camps2012.2018$totGenoeggsCampInc4 != 0 & camps2012.2018$no.eggsCampInc != 0),]
#
#Subset data to exclude groups with only one female
df_groups <- df_groups[df_groups$nofemales != 1,]
#
#Exclude individuals in camps with no female incubation
df_groups <- df_groups[df_groups$Ftotsit != 0,]
#
#Subset the data to include only columns of interest
names(df_groups)
df_groups <- df_groups[c("yearcamp", "camp",  "year", "nonests", "noFCoopsCamp", "noFCheatsCamp",
                           "noFOtherCamp", "nomales", "nofemales", "no.eggsCampInc", "no.BrokenInc",
                           "prop.BrokenInc", "no.HatchedInc", "FCheatinCamp", "noFActive", "prop.Cheats")]
#
#Check the frequency of number of helpers in camps
table(df_groups$noFOtherCamp)#Have an extra column with presence/absence
df_groups$FOtherinCamp <- ifelse(df_groups$noFOtherCamp == 0, "No", "Yes")
table(df_groups$FOtherinCamp)
#
#Compute hatching success
df_groups$Hatching_success <- df_groups$no.HatchedInc / df_groups$no.eggsCampInc
#
#
####df_hatching and df_groups_all: Group-level hatching outcomes----
#********************************************************************
#Create a data frame to use when drawing the figure
df_groups_all <- df_groups
#Create a column with info of whether a given row is included in the model data
df_groups_all$InModel <- ifelse(df_groups_all$noFCoopsCamp != 0 & df_groups_all$noFCoopsCamp != 6,
                                      "Yes", "No")
#
#Subset data to exclude individuals in groups with no, and with six cooperators
#(no cheats are present in these groups, so no comparison is possible between groups with and without cheats)
df_hatching <- df_groups_all[df_groups_all$noFCoops != 0 & df_groups_all$noFCoops != 6,]
#
#Format variables
str(df_hatching)
df_hatching$year <- factor(df_hatching$year)
#
#Compute Z-transformed numerical variables
df_hatching$Znomales <- scale(df_hatching$nomales)
df_hatching$Znofemales <- scale(df_hatching$nofemales)
df_hatching$ZnoFCheats <- scale(df_hatching$noFCheatsCamp)
df_hatching$ZnoFOther <- scale(df_hatching$noFOtherCamp)
df_hatching$ZnoFCoops <- scale(df_hatching$noFCoopsCamp)
#
#Compute logit values for proportions
df_hatching$logit_prop.BrokenInc <- car::logit(df_hatching$prop.BrokenInc)
#
#Get quadratics
df_hatching$Znomales2 <- df_hatching$Znomales * df_hatching$Znomales
df_hatching$Znofemales2 <- df_hatching$Znofemales * df_hatching$Znofemales
df_hatching$ZnoFCheats2 <- df_hatching$ZnoFCheats * df_hatching$ZnoFCheats
df_hatching$ZnoFOther2 <- df_hatching$ZnoFOther * df_hatching$ZnoFOther
df_hatching$ZnoFCoops2 <- df_hatching$ZnoFCoops * df_hatching$ZnoFCoops
df_hatching$logit_prop.BrokenInc2 <- df_hatching$logit_prop.BrokenInc * df_hatching$logit_prop.BrokenInc
#
#
####df_incubation_group: Group-level incubation effort----
#***********************************************************************************
df_incubation_group <- df_hatching
#
#Match female incubation at primary nest (accounting for overlapping incubatioon records)
df_incubation_group$Finc_N1 <- incNest_F$Finc_N1[match(df_incubation_group$yearcamp, incNest_F$yearcamp)]
#
#Match total observation time
df_incubation_group$totobstime <- camps2012.2018$totobstime[match(df_incubation_group$yearcamp, camps2012.2018$yearcamp)]
#
#
#******#
### Female aggression----
#******#
#
####df_aggression: Aggression received by females of each breeding strategy----
#************************************************************************************
#
df_aggression <- Femfem_coop_agg_count_care4_grp_long
#
#Subset data to only include camps where eggs were laid during the period with natural incubation
#and eggs were genotyped during the period with natural incubation and 4 weeks before
df_aggression <- df_aggression[which(df_aggression$totGenoeggsCampInc4 != 0 & df_aggression$no.eggsCampInc != 0),]
#
#Subset data to exclude groups with only one female
df_aggression <- df_aggression[df_aggression$nofemales != 1,]
#
#Subset data to exclude groups without cheats
df_aggression <- df_aggression[df_aggression$no.FCheats != 0,]
#
#Subset data to exclude groups with only one cooperator (no aggression possible towards self)
df_aggression <- df_aggression[df_aggression$no.FCoops != 1,]
#
#Subset data to only include Cheats and cooperators
df_aggression <- df_aggression[df_aggression$BSvic != "Other",]
#
#Match camp info
df_aggression$camp <- camps2012.2018$camp[match(df_aggression$yearcamp, camps2012.2018$yearcamp)]
names(df_aggression)
#
#Compute number of females with the same breeding strategy as the victim
df_aggression = df_aggression %>%
  mutate(noBSvic = ifelse(BSvic == "Cooperator", no.FCoops - 1, no.FCheats))
#
#Compute Z-transformed numerical variables
df_aggression$Znomales <- scale(df_aggression$nomales)
df_aggression$Znofemales <- scale(df_aggression$nofemales)
df_aggression$Zno.FCoops <- scale(df_aggression$no.FCoops)
df_aggression$Zno.FOther <- scale(df_aggression$no.FOther)
df_aggression$ZnoBSvic <- scale(df_aggression$noBSvic)
df_aggression$ZtotObstime  <- scale(df_aggression$totObstime)
#
#Compute quadratics
df_aggression$Znomales2 <- df_aggression$Znomales * df_aggression$Znomales
df_aggression$Znofemales2 <- df_aggression$Znofemales * df_aggression$Znofemales
df_aggression$Zno.FCoops2 <- df_aggression$Zno.FCoops * df_aggression$Zno.FCoops
df_aggression$Zno.FOther2 <- df_aggression$Zno.FOther * df_aggression$Zno.FOther
df_aggression$ZnoBSvic2 <- df_aggression$ZnoBSvic * df_aggression$ZnoBSvic
df_aggression$ZtotObstime2 <- df_aggression$ZtotObstime * df_aggression$ZtotObstime
#
#
#******#
###df_laying_window: Probability of an individual laying within a two-day window----
#******#
#
#Compute data
#Include only cheats and cooperators
df_laying_window <- indLayingwZeros[indLayingwZeros$BS != "Other",]
#
#Create two-day windows
df_laying_window$Window_2days <- df_laying_window$daysfromEggsleft - (df_laying_window$daysfromEggsleft %% 2)
#Compute probability of laying within a two-day window
df_laying_window <- df_laying_window %>%
  group_by(Window_2days, yearID) %>%
  mutate(probEgg_2days = as.integer(any(probEgg == 1)))
#
#Compute a unique id for every individual and 2-day window.
df_laying_window$yearID_2days <- paste(df_laying_window$yearID, df_laying_window$Window_2days, sep = "_")  
#Keep only unique values for yearID_2days
df_laying_window <- df_laying_window[!duplicated(df_laying_window$yearID_2days), ]
#Get rid of unnecessary columns
names(df_laying_window)
df_laying_window <- df_laying_window %>% dplyr::select(yearID, year, yearcamp, BS, Window_2days, probEgg_2days)
#
#Include number of females in group
df_laying_window$nofemales <- camps2012.2018$nofemales[match(df_laying_window$yearcamp, camps2012.2018$yearcamp)]
#
#Subset data to exclude groups with only one female, as in every other analysis.
#A female breeding alone has no one to cooperate with or cheat on, so her
#breeding strategy is not comparable with that of a female in a group.
df_laying_window <- df_laying_window[df_laying_window$nofemales != 1,]
#
#Apply the same group-level exclusions as the other analyses, so that all of
#them describe the same set of groups: camps where eggs were laid and genotyped
#during the period with natural incubation, that had a nest, and where females
#incubated.
df_laying_window$totGenoeggsCampInc4 <- camps2012.2018$totGenoeggsCampInc4[match(df_laying_window$yearcamp, camps2012.2018$yearcamp)]
df_laying_window$no.eggsCampInc <- camps2012.2018$no.eggsCampInc[match(df_laying_window$yearcamp, camps2012.2018$yearcamp)]
df_laying_window <- df_laying_window[which(df_laying_window$totGenoeggsCampInc4 != 0 & df_laying_window$no.eggsCampInc != 0),]
#
df_laying_window$nonests <- camps2012.2018$nonests[match(df_laying_window$yearcamp, camps2012.2018$yearcamp)]
df_laying_window <- df_laying_window[df_laying_window$nonests != 0,]
#
df_laying_window$Ftotsit <- camps2012.2018$Ftotsit[match(df_laying_window$yearcamp, camps2012.2018$yearcamp)]
df_laying_window <- df_laying_window[df_laying_window$Ftotsit != 0,]
#
#Add info on camp
df_laying_window$camp <- camps2012.2018$camp[match(df_laying_window$yearcamp, camps2012.2018$yearcamp)]
#Add info on id
df_laying_window$id <- treatments2012.2018$id[match(df_laying_window$yearID, treatments2012.2018$yearID)]
#
#Include only females
df_laying_window$sex <- treatments2012.2018$sex[match(df_laying_window$id, treatments2012.2018$id)]
df_laying_window <- df_laying_window[df_laying_window$sex != "M",]
#
#Compute z-values, after the exclusions so that the scaling reflects the
#rows that are modelled
df_laying_window$ZWindow_2days <- scale(df_laying_window$Window_2days)
df_laying_window$Znofemales <- scale(df_laying_window$nofemales)
#Compute quadratics
df_laying_window$ZWindow_2days2 <- df_laying_window$ZWindow_2days * df_laying_window$ZWindow_2days
#
#Convert to df
df_laying_window <- as.data.frame(df_laying_window)
#
#
#******#
###df_laying_group: Group-level egg laying across time----
#******#
#
#Create a data frame with the data
df_laying_group <- eggs_date_camp
#
#Subset data to only include camps where eggs were laid during the period with natural incubation
#and eggs were genotyped during the period with natural incubation and 4 weeks before
#Match data
df_laying_group$totGenoeggsCampInc4 <- camps2012.2018$totGenoeggsCampInc4[match(df_laying_group$yearcamp, camps2012.2018$yearcamp)]
df_laying_group$no.eggsCampInc <- camps2012.2018$no.eggsCampInc[match(df_laying_group$yearcamp, camps2012.2018$yearcamp)]
df_laying_group <- df_laying_group[which(df_laying_group$totGenoeggsCampInc4 != 0 & df_laying_group$no.eggsCampInc != 0),]
#
#Exclude individuals in camps with no nests
#Match data
df_laying_group$nonests <- camps2012.2018$nonests[match(df_laying_group$yearcamp, camps2012.2018$yearcamp)]
table(df_laying_group$nonests, exclude = NULL)
df_laying_group <- df_laying_group[df_laying_group$nonests != 0,]
#
#Exclude individuals in camps with no female incubation
df_laying_group$Ftotsit <- camps2012.2018$Ftotsit[match(df_laying_group$yearcamp, camps2012.2018$yearcamp)]
df_laying_group <- df_laying_group[df_laying_group$Ftotsit != 0,]
#
#Match info on whether a camp has a cheat or not
df_laying_group$CheatinCamp <- camps2012.2018$FCheatinCamp[match(df_laying_group$yearcamp, camps2012.2018$yearcamp)]
#Remove NAs
df_laying_group <- na.omit(df_laying_group)
#
#Compute remaining unlaid eggs
df_laying_group <- df_laying_group %>%
  arrange(yearcamp, dayorderEggs) %>%
  group_by(yearcamp) %>%   
  mutate(remaining = tot_no.eggs - dplyr::lag(cumulative_eggs, default = 0)) %>%
  ungroup() %>%
  as.data.frame()
#
#Match info on group size, year and camp
df_laying_group <- df_laying_group %>%
  left_join(camps2012.2018 %>%
              dplyr::select(yearcamp, nofemales, noCoops = noFCoopsCamp, noOther = noFOtherCamp, year, camp),
            by = "yearcamp")
#
#Subset data to exclude groups with only one female, as in every other analysis.
#A female nesting alone cannot have a cheat in her group, so these nests would
#otherwise contribute only to the cheat-absent category while following a very
#different laying trajectory from a shared nest.
df_laying_group <- df_laying_group[df_laying_group$nofemales != 1,]
#
#Scale time- and group size variables, after the exclusions so that the scaling
#reflects the rows that are modelled
df_laying_group$ZDay <- scale(df_laying_group$dayorderEggs)
df_laying_group$Znofemales <- scale(df_laying_group$nofemales)
df_laying_group$ZnoCoops <- scale(df_laying_group$noCoops)
df_laying_group$ZnoOther <- scale(df_laying_group$noOther)
#
#Compute quadratic and cubic time variables
df_laying_group$ZDay2 <- df_laying_group$ZDay^2
df_laying_group$ZDay3 <- df_laying_group$ZDay^3
#
#
#******#
###df_incubation_time: Group-level incubation effort across time----
#******#
#
df_incubation_time <- sitdur_date_filtered
#Exclude individuals in camps with no nests
#Match data
df_incubation_time$nonests <- camps2012.2018$nonests[match(df_incubation_time$yearcamp, camps2012.2018$yearcamp)]
table(df_incubation_time$nonests, exclude = NULL)#No camps with no nests
#
#Exclude individuals in camps with no female incubation
df_incubation_time$Ftotsit <- camps2012.2018$Ftotsit[match(df_incubation_time$yearcamp, camps2012.2018$yearcamp)]
df_incubation_time <- df_incubation_time[df_incubation_time$Ftotsit != 0,]
#
#Match info on whether a camp has a cheat or not
df_incubation_time$CheatinCamp <- camps2012.2018$FCheatinCamp[match(df_incubation_time$yearcamp, camps2012.2018$yearcamp)]
#
#Match info on group size and camp
df_incubation_time <- df_incubation_time %>%
  left_join(camps2012.2018 %>%
              dplyr::select(yearcamp, noCoops = noFCoopsCamp, noOther = noFOtherCamp, camp),
            by = "yearcamp")
#
#Scale time- and group size variables
df_incubation_time$ZDay <- scale(df_incubation_time$dayordereggsleft)
df_incubation_time$Znofemales <- scale(df_incubation_time$nofemales)
df_incubation_time$ZnoCoops <- scale(df_incubation_time$noCoops)
df_incubation_time$ZnoOther <- scale(df_incubation_time$noOther)
#Compute quadratic and cubic time variables
df_incubation_time$ZDay2 <- df_incubation_time$ZDay^2
df_incubation_time$ZDay3 <- df_incubation_time$ZDay^3
df_incubation_time <- as.data.frame(df_incubation_time)
#
#
#******#
###df_activity: Whether a nest was active on a given day----
#******#
#
str(df_laying_group)
str(df_incubation_time)
df_activity = df_incubation_time %>% #This is the data with incubation per day
  dplyr::select(yearcamp, date, sitdur_date, obsdur_date, dayordereggsleft, propsit, year, totGenoeggsCampInc4, no.eggsCampInc,
         nofemales, nonests, Ftotsit, CheatinCamp, noCoops, noOther, camp) %>%
  #Match eggs laid in a given day
  left_join(df_laying_group %>% dplyr::select(no.eggs, cumulative_eggs, laydate_date, yearcamp), by = c("date" = "laydate_date", "yearcamp")) %>%
  #Create a column with activity, scale time- and group size variables, and compute quadratic of time variable
  mutate(activity = ifelse(sitdur_date == 0 & no.eggs == 0, 0, 1), ZDay = scale(dayordereggsleft),
         Znofemales = scale(nofemales), ZnoCoops = scale(noCoops), ZnoOther = scale(noOther),
         ZDay2 = ZDay^2)
#
#
#******#
###df_activity_allocation: Split of daily activity between laying and incubation----
#******#
#
#For days with activity, check how activity is allocated between egg-laying and incubating
df_activity_allocation = df_activity %>%
  dplyr::filter(activity == 1) %>%
  #Compute propotion of eggs laid (a proportion of the population maximum per day = 6), and activity allocation: 1 = only incubation, 0 = only egg-laying.
  ##qlogis transformation, avoids exact 0/1
  mutate(prop_eggs = no.eggs / 6, act_alloc = propsit / (propsit + prop_eggs), logit_alloc = qlogis(pmin(pmax(act_alloc, 0.01), 0.99)))
#
##**********##
#Save data----
##**********##
#
#Only the frames that something downstream reads are saved. The rest are
#stepping stones that exist while this script runs.
save(df_females, df_strategy_counts, df_cheat_prop, df_cheat_prob,
     df_eggs, df_eggs_excl_1coop, df_eggs_all, df_mass, df_eggs_within_female,
     df_incubation_individual, df_hatching, df_incubation_group, df_aggression,
     df_laying_window, df_laying_group, df_incubation_time,
     df_activity, df_activity_allocation,
     treatments2012.2018, sitdur_date_ind,
     file = "R_Output_files/Cheating_Model_Data.RData")
#
#
#