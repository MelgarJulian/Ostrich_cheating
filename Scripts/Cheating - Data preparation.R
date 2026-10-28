###*************************###
# Frequency-dependent cheating in a cooperative breeding bird - Data Preparation ----
###*************************###
#
#
##*****************##
## Load packages----
##*****************##
#
pacman::p_load(magrittr, dplyr, doBy, data.table, reshape, tidyr, tidyverse, date, segmented, purrr, ggplot2, conflicted)
#
#
##*************##
## Load data----
##*************##
#
treatments2012.2018 <-
  (read.csv("R_Input_files/treatments2012.2018 - 20.05.2020.csv"))
sitduration2012.2018 <-
  read.csv("R_Input_files/sitduration2012.2018 - 13.07.2024.csv", sep = "")
eggs2012.2018 <-
  read.csv("R_Input_files/eggs2012.2018 - 23.04.2024.csv")
parentage2012.2018 <- 
  read.csv("R_Input_files/parentage2012.2018.csv")
copulation2012.2018 <-
  read.csv("R_Input_files/copulation2012.2018 - 31.08.2020.csv")
femfem2012.2018 <- 
  read.csv("R_Input_files/femfem2012.2018 - 29.07.2022.csv")
malefem2012.2018 <- 
  read.csv("R_Input_files/malefem2012.2018 - 29.07.2022.csv")
mass <- read.csv("R_Input_files/mass.csv")
#
#
##*********************************************##
## General data formatting and data filtering----
##*********************************************##
#
#******#
### Treatments----
#******#
#
#Format date columns as date
treatments2012.2018$eggsleftDate <- as.Date(treatments2012.2018$eggsleftDate, "%Y-%m-%d")
treatments2012.2018$seasonstart_date <- as.Date(treatments2012.2018$seasonstart_date, "%Y-%m-%d")
treatments2012.2018$seasonend_date <- as.Date(treatments2012.2018$seasonend_date, "%Y-%m-%d")
#
#******#
### sitduration2012.2018----
#******#
#
#Format date columns as date
sitduration2012.2018$date <- as.Date(sitduration2012.2018$date, "%Y-%m-%d")
#
#Filter data from after day 32 of natural incubation
#sitduration2012.2018 <- sitduration2012.2018 %>% dplyr::filter(dayorderEggs < 33)
#
#
#******#
### eggs2012.2018----
#******#
#
#Filter data from after day 32 of natural incubation
#eggs2012.2018 <- eggs2012.2018 %>% dplyr::filter(daysfromEggsleft < 33)
#
#******#
### parentage2012.2018----
#******#
#
#Format date columns as date
str(parentage2012.2018)
parentage2012.2018$laydate <- as.Date(parentage2012.2018$laydate, "%Y-%m-%d")
parentage2012.2018$eggsleftDate <- as.Date(parentage2012.2018$eggsleftDate, "%Y-%m-%d")
#
#Filter data from after day 32 of natural incubation
#parentage2012.2018 <- parentage2012.2018 %>% dplyr::filter(daysfromEggsleft < 33)
#
#******#
### copulation2012.2018----
#******#
#
#Filter data from after day 32 of natural incubation
#copulation2012.2018 <- copulation2012.2018 %>% dplyr::filter(dayorderEggs < 33)
#
#******#
### femfem2012.2018----
#******#
#
#Filter data from after day 32 of natural incubation
#femfem2012.2018 <- femfem2012.2018 %>% dplyr::filter(dayorderEggs < 33)
#
#******#
### malefem2012.2018----
#******#
#
#Filter data from after day 32 of natural incubation
#malefem2012.2018 <- malefem2012.2018 %>% dplyr::filter(dayorderEggs < 33)
#
#******#
### Mass----
#******#
#
#Format year and create a year_id column
mass = mass |> dplyr::mutate(year=yr+1900,
                      year_id=paste(year,id,sep="_"))
#
#
#***************************#
## Create new data frames----
#***************************#
#
#******#
### Summary data per individual and year----
#******#
#
#create a df with one row per each individual present in the period with care
individuals2012.2018 <- treatments2012.2018[which(treatments2012.2018$incampeggsleft == "Yes"),]
#
#Check whether there are individuals in camps with removal score of 4
#This means that an individual was replaced less than 14 days before the onset of natural incubation,
#and that the sex-specific group size was reduced to less that 5/6 of its original size.
table(individuals2012.2018$Removal.replacement.score, exclude = NULL)#18 individuals affected
#check
individuals2012.2018[which(individuals2012.2018$Removal.replacement.score == 4),]#3 camps affected:
#2015_25: 1 male (out of three males in camp)
#removed late in season (12/01)
#2017_22: 1 female (out of 3)
#removed before period with care started.
#2017_29: The only female removed
#early in period with care (11/13)
#Remove individuals in camps with a removal score of 4
#Replace NAs with 0s
individuals2012.2018$Removal.replacement.score[which(is.na(individuals2012.2018$Removal.replacement.score))] <- 0
#Remove individuals
individuals2012.2018 <- individuals2012.2018[which(individuals2012.2018$Removal.replacement.score != 4),]
#
#Remove irrelevant columns
names(individuals2012.2018)
individuals2012.2018 <- individuals2012.2018[c("year", "camp", "yearcamp", "grpsize", "nomales", "nofemales", "id", "age", "sex", "idcode",
                                               "yearID")]
#
#Compute total observation time per year during the period with care
#Compute the data
totobsdurYear <- 
  summaryBy(obsduration ~ year, data = sitduration2012.2018 %>%
              distinct(session, .keep_all = TRUE), FUN = sum)#OBS! This only gives observations bouts when incubation was observed!
#Add info on total observation duration per year
individuals2012.2018$totobsdurYear <- totobsdurYear$obsduration.sum[match(individuals2012.2018$year, totobsdurYear$year)]
#
#
#### Incubation----
#
#Compute total incubation time for every group
sitdurcamp <- summaryBy(sitduration ~ yearcamp, data = sitduration2012.2018, FUN = sum)
#Match data to the individuals df
individuals2012.2018$totSitdurCamp <- sitdurcamp$sitduration.sum[match(individuals2012.2018$yearcamp, sitdurcamp$yearcamp)]
#Change NAs to zeros
table(individuals2012.2018$totSitdurCamp, exclude = NULL)
individuals2012.2018$totSitdurCamp[which(is.na(individuals2012.2018$totSitdurCamp))] <- 0
table(individuals2012.2018$totSitdurCamp, exclude = NULL)
#
#Compute total female incubation time for every group
Femsitdurcamp <- summaryBy(sitduration ~ yearcamp, data = sitduration2012.2018[sitduration2012.2018$sex == "F",], FUN = sum)
#Match data to the individuals df
individuals2012.2018$femSitdurCamp <- Femsitdurcamp$sitduration.sum[match(individuals2012.2018$yearcamp, Femsitdurcamp$yearcamp)]
#Change NAs to zeros
table(individuals2012.2018$femSitdurCamp, exclude = NULL)
individuals2012.2018$femSitdurCamp[which(is.na(individuals2012.2018$femSitdurCamp))] <- 0
table(individuals2012.2018$femSitdurCamp, exclude = NULL)
#
#Visualise female incubation duration
plot(table(sitduration2012.2018$sitduration[sitduration2012.2018$sex == "F"]))
#
#There are some incubations that are only one minute long. These could be females at the nest laying eggs, or soliciting. Get rid!
sitduration2012.2018 <- sitduration2012.2018[sitduration2012.2018$sitduration != 1,]
#
#Compute time that every individual spent sitting
sitdurbird <- summaryBy(sitduration ~ yearID + yearcamp + sex, data = sitduration2012.2018, FUN = sum)
#
#Add info on time spent sitting per individual in the individuals2012.2018 df
individuals2012.2018$sitdurbird <- sitdurbird$sitduration.sum[match(individuals2012.2018$yearID, sitdurbird$yearID)]
#Change NAs to zeros
individuals2012.2018$sitdurbird[is.na(individuals2012.2018$sitdurbird)] <- 0
#
#Compute time that individuals spent incubating relative to total observation time
individuals2012.2018$tot.prop.IncInd <- individuals2012.2018$sitdurbird / individuals2012.2018$totobsdurYear
#
#Compute the number of incubations per individuals
incsind <- as.data.frame(table(sitduration2012.2018$yearID)) %>%
  dplyr::rename(yearID = Var1, times_inc = Freq)
#Match data
individuals2012.2018$noIncs_year <- incsind$times_inc[match(individuals2012.2018$yearID, incsind$yearID)]
#Change NAs to zeros
individuals2012.2018$noIncs_year[is.na(individuals2012.2018$noIncs_year)] <- 0
#Visualize
plot(table(individuals2012.2018$noIncs_year))
#
#
####Parentage----
#Create a column with parentage score for individual during the four weeks previous to onset of incubation 
#Compute data
#Females
parentageW4F <- summaryBy(sampleID ~ yearID_Mother, data = parentage2012.2018[which(parentage2012.2018$daysfromEggsleft < 0 &
                                                                                      parentage2012.2018$daysfromEggsleft >= -28),],
                          FUN = length)
#Males
parentageW4M <- summaryBy(sampleID ~ yearID_Father, data = parentage2012.2018[which(parentage2012.2018$daysfromEggsleft < 0 &
                                                                                      parentage2012.2018$daysfromEggsleft >= -28),],
                          FUN = length)
#Bind dfs
#Re-name columns
parentageW4F <- dplyr::rename(parentageW4F, yearID = yearID_Mother)
parentageW4M <- dplyr::rename(parentageW4M, yearID = yearID_Father)
#Bind
parentageW4 <- rbind(parentageW4F, parentageW4M)
#Match data
individuals2012.2018$parentageW4 <- parentageW4$sampleID[match(individuals2012.2018$yearID, parentageW4$yearID)]
#Change NAs to zeros
individuals2012.2018$parentageW4[is.na(individuals2012.2018$parentageW4)] <- 0
#
#Create a column with parentage score for camps during the four weeks previous to onset of incubation 
#Compute data
CampParentageW4 <- summaryBy(sampleID ~ yearcamp, data = parentage2012.2018[which(parentage2012.2018$daysfromEggsleft < 0 &
                                                                                    parentage2012.2018$daysfromEggsleft >= -28),],
                             FUN = length)
#Match data
individuals2012.2018$CampParentageW4 <- CampParentageW4$sampleID[match(individuals2012.2018$yearcamp, CampParentageW4$yearcamp)]
#Change NAs to zeros
individuals2012.2018$CampParentageW4[is.na(individuals2012.2018$CampParentageW4)] <- 0
#
#Compute number of eggs produced per individual during the period with care
#Females
parentageIncF <- summaryBy(sampleID ~ yearID_Mother, data = parentage2012.2018[which(parentage2012.2018$daysfromEggsleft < 0 &
                                                                                       parentage2012.2018$daysfromEggsleft >= -28),],
                           FUN = length)
#Males
parentageIncM <- summaryBy(sampleID ~ yearID_Father, data = parentage2012.2018[which(parentage2012.2018$daysfromEggsleft < 0 &
                                                                                       parentage2012.2018$daysfromEggsleft >= -28),],
                           FUN = length)
#Bind dfs
#Re-name columns
parentageIncF <- dplyr::rename(parentageIncF, yearID = yearID_Mother)
parentageIncM <- dplyr::rename(parentageIncM, yearID = yearID_Father)
#Bind
parentageInc <- rbind(parentageIncF, parentageIncM)
#Match data
individuals2012.2018$parentageInc <- parentageInc$sampleID[match(individuals2012.2018$yearID, parentageInc$yearID)]
#Change NAs to zeros
individuals2012.2018$parentageInc[is.na(individuals2012.2018$parentageInc)] <- 0
#
####Copulation----
#Compute the total number of copulations for every male
#Compute year_ID for the males in the copulation df
copulation2012.2018$year_IDM <- paste(copulation2012.2018$year, copulation2012.2018$idMale, sep = "_")
copsMale <- as.data.frame(table(copulation2012.2018$year_IDM))
#Compute the total number of copulations for every female
#Compute year_ID for the females in the copulation df
copulation2012.2018$year_IDF <- paste(copulation2012.2018$year, copulation2012.2018$idFemale, sep = "_")
copsFemale <-  as.data.frame(table(copulation2012.2018$year_IDF))
#Merge the male and the female copulation dfs
copsInd <- rbind(copsMale, copsFemale)
#Match copulations per individual in the individuals df
individuals2012.2018$no.CopsInd <- copsInd$Freq[match(individuals2012.2018$yearID, copsInd$Var1)]
#Change NAs to 0s
individuals2012.2018$no.CopsInd[is.na(individuals2012.2018$no.CopsInd)] <- 0
#
#Compute number of copulations during the period with care and 28 days earlier (see Malecki et al 2004)
#Compute the total number of copulations during the period with care for every male
copsMaleInc <- as.data.frame(table(copulation2012.2018$year_IDM[copulation2012.2018$cop.period == "P2"]))
#Compute the total number of copulations during the period with care for every female
copsFemaleInc <-  as.data.frame(table(copulation2012.2018$year_IDF[copulation2012.2018$cop.period == "P2"]))
#Merge the male and the female copulation dfs
copsIndInc <- rbind(copsMaleInc, copsFemaleInc)
#Match copulations per individual in the individuals df
individuals2012.2018$no.CopsIndInc <- copsIndInc$Freq[match(individuals2012.2018$yearID, copsIndInc$Var1)]
#Change NAs to 0s
individuals2012.2018$no.CopsIndInc[is.na(individuals2012.2018$no.CopsIndInc)] <- 0
#
#Compute successful copulations
#Compute the total number of successful copulations for every male
succ.copsMale <- as.data.frame(table(copulation2012.2018$year_IDM[copulation2012.2018$succ == "Y"]))
#Compute the total number of successful copulations for every female
succ.copsFemale <-  as.data.frame(table(copulation2012.2018$year_IDF[copulation2012.2018$succ == "Y"]))
#Merge the male and the female successful copulation dfs
succ.copsInd <- rbind(succ.copsMale, succ.copsFemale)
#Match copulations per individual in the individuals df
individuals2012.2018$succ.CopsInd <- succ.copsInd$Freq[match(individuals2012.2018$yearID, succ.copsInd$Var1)]
#Change NAs to 0s
individuals2012.2018$succ.CopsInd[is.na(individuals2012.2018$succ.CopsInd)] <- 0
#
#Compute successful copulations during the period with care care and 28 days earlier
#Compute the total number of successful copulations during the period with care for every male
succCopsMaleInc <- as.data.frame(table(copulation2012.2018$year_IDM[copulation2012.2018$succ == "Y" &
                                                                      copulation2012.2018$cop.period == "P2"]))
#Compute the total number of successful copulations during the period with care for every female
succCopsFemaleInc <-  as.data.frame(table(copulation2012.2018$year_IDF[copulation2012.2018$succ == "Y" &
                                                                         copulation2012.2018$cop.period == "P2"]))
#Merge the male and the female successful copulation dfs
succCopsIndInc <- rbind(succCopsMaleInc, succCopsFemaleInc)
#Match successful copulations per individual in the individuals df
individuals2012.2018$no.succCopsIndInc <- succCopsIndInc$Freq[match(individuals2012.2018$yearID, succCopsIndInc$Var1)]
#Change NAs to 0s
individuals2012.2018$no.succCopsIndInc[is.na(individuals2012.2018$no.succCopsIndInc)] <- 0
#
#Compute successful copulations for every individual during the period without care (up to 28 days before onset of care)
#Compute the total number of successful copulations during the period without care for every male
succCopsMaleNC <- as.data.frame(table(copulation2012.2018$year_IDM[copulation2012.2018$succ == "Y" &
                                                                     copulation2012.2018$cop.period == "P1"]))
#Compute the total number of successful copulations during the period without care for every female
succCopsFemaleNC <-  as.data.frame(table(copulation2012.2018$year_IDF[copulation2012.2018$succ == "Y" &
                                                                        copulation2012.2018$cop.period == "P1"]))
#Merge the male and the female successful copulation dfs
succCopsIndNC <- rbind(succCopsMaleNC, succCopsFemaleNC)
#Match successful copulations per individual in the individuals df
individuals2012.2018$no.succCopsIndNC <- succCopsIndNC$Freq[match(individuals2012.2018$yearID, succCopsIndNC$Var1)]
#Change NAs to 0s
individuals2012.2018$no.succCopsIndNC[is.na(individuals2012.2018$no.succCopsIndNC)] <- 0
#
#Create a column with copulation probability for every individual
#(1 for individuals that were observed copulating, and 0 for individuals that never were)
individuals2012.2018$Cop.prob <- ifelse(individuals2012.2018$no.CopsInd != 0, 1, 0)
#
#Create a column with copulation probability for every individual during the period with care
#(1 for individuals that were observed copulating, and 0 for individuals that never were)
individuals2012.2018$Cop.probInc <- ifelse(individuals2012.2018$no.CopsIndInc != 0, 1, 0)
#
#Create a column with successful copulation probability for every individual
#(1 for individuals that were observed copulating successfully, and 0 for individuals that never were)
individuals2012.2018$SuccCop.prob <- ifelse(individuals2012.2018$succ.CopsInd != 0, 1, 0)
#
#Create a column with successful copulation probability for every individual during the period with care
#(1 for individuals that were observed copulating successfully,
#and 0 for individuals that never were)
individuals2012.2018$SuccCop.probInc <- ifelse(individuals2012.2018$no.succCopsIndInc != 0, 1, 0)
#
####Reproductive success----
#####Eggs----
######Total----
#Compute total reproductive success for every individual
#Compute the total number of eggs with genotype for every male
GenoMale <- as.data.frame(table(parentage2012.2018$yearID_Father))
#Compute the total number of eggs with genotype for every female
GenoFemale <-  as.data.frame(table(parentage2012.2018$yearID_Mother))
#Merge the male and the female Geno dfs
GenoInd <- rbind(GenoMale, GenoFemale)
#Match number of eggs with genotype per individual in the individuals df
individuals2012.2018$RS <- GenoInd$Freq[match(individuals2012.2018$yearID, GenoInd$Var1)]
#Change NAs to 0s
individuals2012.2018$RS[is.na(individuals2012.2018$RS)] <- 0
#
######Period with care----
#Compute reproductive success for every individual during the period with care
#Compute the total number of eggs with genotype during the period with care for every male
GenoMaleInc <- as.data.frame(table(parentage2012.2018$yearID_Father[parentage2012.2018$eggsleft == "Yes"]))
#Compute the total number of eggs with genotype during the period with care for every female
GenoFemaleInc <- as.data.frame(table(parentage2012.2018$yearID_Mother[parentage2012.2018$eggsleft == "Yes"]))
#Merge the male and the female Geno dfs
GenoIndInc <- rbind(GenoMaleInc, GenoFemaleInc)
#Match number of eggs with genotype per individual in the individuals df
individuals2012.2018$RSInc <- GenoIndInc$Freq[match(individuals2012.2018$yearID, GenoIndInc$Var1)]
#Change NAs to 0s
individuals2012.2018$RSInc[is.na(individuals2012.2018$RSInc)] <- 0
#
######Period with care and four weeks before----
#Compute reproduction success for every individual during the period with care and four weeks before.
individuals2012.2018$RSInc4 <- individuals2012.2018$parentageW4 + individuals2012.2018$RSInc
#
#Compute sex-specific mean RSInc4 per year
individuals2012.2018 <- individuals2012.2018 %>%
  group_by(year, sex) %>%
  dplyr::mutate(meanRSInc4 = mean(RSInc4, na.rm = TRUE)) %>% ungroup()
#Compute RSInc4 relative to year- and sex-specific mean
individuals2012.2018$relativeRSInc4 <- individuals2012.2018$RSInc4 / individuals2012.2018$meanRSInc4
#
#### Hatching success----
#For every male, compute the number of hatched eggs out of the eggs with genotype laid during the period with care
HatchMale_Inc <- as.data.frame(table(parentage2012.2018$yearID_Father[which(parentage2012.2018$eggsleft == "Yes")]))
#For every female, compute the number of hatched eggs out of the eggs with genotype laid during the period with care and four weeks before
HatchFemale_Inc <- as.data.frame(table(parentage2012.2018$yearID_Mother[which(parentage2012.2018$eggsleft == "Yes")]))
#Merge the male and the female hatching dfs
HatchIndInc <- rbind(HatchMale_Inc, HatchFemale_Inc)
#Match number of hatched eggs with genotype per individual in the individuals df
individuals2012.2018$HatchIndInc <- HatchIndInc$Freq[match(individuals2012.2018$yearID, HatchIndInc$Var1)]
#Change NAs to 0s
individuals2012.2018$HatchIndInc[is.na(individuals2012.2018$HatchIndInc)] <- 0
#
#For every male, compute the number of hatched eggs out of the eggs with genotype laid during the period with care and four weeks before
HatchMale_Inc4 <- as.data.frame(table(parentage2012.2018$yearID_Father[which(parentage2012.2018$daysfromEggsleft >= -28 &
                                                                               parentage2012.2018$c_hatched == 1)]))
#For every female, compute the number of hatched eggs out of the eggs with genotype laid during the period with care and four weeks before
HatchFemale_Inc4 <- as.data.frame(table(parentage2012.2018$yearID_Mother[which(parentage2012.2018$daysfromEggsleft >= -28 &
                                                                                 parentage2012.2018$c_hatched == 1)]))
#Merge the male and the female hatching dfs
HatchIndInc4 <- rbind(HatchMale_Inc4, HatchFemale_Inc4)
#Match number of hatched eggs with genotype per individual in the individuals df
individuals2012.2018$HatchIndInc4 <- HatchIndInc4$Freq[match(individuals2012.2018$yearID, HatchIndInc4$Var1)]
#Change NAs to 0s
individuals2012.2018$HatchIndInc4[is.na(individuals2012.2018$HatchIndInc4)] <- 0
#
#### Breeding scores----
#
#Compute the number of reproductive females in each group (using RS during period with care and four weeks before for cheats)
#Cooperator = incubates for a total of more than 3 minutes, or has been observed incubating more than once
#Cheat = produces offspring but does not classify to be a cooperator
#Other = does not produce offspring nor classifies to be a cooperator
#
##OBS! Breeding score/strategy was changed from Breed_2RSInc4 to BS 03/11/2023
individuals2012.2018 <- individuals2012.2018 %>%
  dplyr::mutate(BS = ifelse(sitdurbird > 3 | noIncs_year > 1, "Cooperator",
                            ifelse(RSInc4 != 0 & sitdurbird <= 3 & noIncs_year <= 1, "Cheat", "Other")))
table(individuals2012.2018$BS, exclude = NULL)#There are no "extra" birds. Good!
range(individuals2012.2018$sitdurbird[individuals2012.2018$BS == "Cooperator"])#All cooperators incubated. Good!
range(individuals2012.2018$sitdurbird[individuals2012.2018$BS == "Cheat"])#No Cheats incubated for more than 3 minutes. Good!
range(individuals2012.2018$RSInc4[individuals2012.2018$BS == "Cheat"])#All cheats reproduced. Good!
range(individuals2012.2018$sitdurbird[individuals2012.2018$BS == "Other"])#No "other" incubated. Good!
range(individuals2012.2018$RSInc4[individuals2012.2018$BS == "Other"])#No "other" reproduced. Good!
#
######Number of cooperators----
#Compute the number of cooperators
noFCoops <-
  summaryBy(idcode ~ yearcamp, FUN = length, data = individuals2012.2018[which(individuals2012.2018$BS == "Cooperator" &
                                                                                 individuals2012.2018$sex == "F"),])
#Match number of cooperators to the individuals df
individuals2012.2018$noFCoopsCamp <- noFCoops$idcode.length[match(individuals2012.2018$yearcamp, noFCoops$yearcamp)]
#Change NAs to zeros
individuals2012.2018$noFCoopsCamp[is.na(individuals2012.2018$noFCoopsCamp)] <- 0
#
######Number of cheats----
#Compute the number of Cheats
noFCheats <-
  summaryBy(idcode ~ yearcamp, FUN = length, data = individuals2012.2018[which(individuals2012.2018$BS == "Cheat" &
                                                                                 individuals2012.2018$sex == "F"),])
#Match number of Cheats to the individuals df
individuals2012.2018$noFCheatsCamp <- noFCheats$idcode.length[match(individuals2012.2018$yearcamp, noFCheats$yearcamp)]
#Change NAs to zeros
individuals2012.2018$noFCheatsCamp[is.na(individuals2012.2018$noFCheatsCamp)] <- 0
#
######Number of other----
#Compute the number of other
noFOther <-
  summaryBy(idcode ~ yearcamp, FUN = length, data = individuals2012.2018[which(individuals2012.2018$BS == "Other" &
                                                                                 individuals2012.2018$sex == "F"),])
#Match number of other to the individuals df
individuals2012.2018$noFOtherCamp <- noFOther$idcode.length[match(individuals2012.2018$yearcamp, noFOther$yearcamp)]
#Change NAs to zeros
individuals2012.2018$noFOtherCamp[is.na(individuals2012.2018$noFOtherCamp)] <- 0
#
######Number of active females----
#Compute the number of females in breeding group (cooperators + cheats)
individuals2012.2018$noFActive <- individuals2012.2018$noFCoopsCamp + individuals2012.2018$noFCheatsCamp
#
######Proportion of cheats----
#Compute proportion of cheats
individuals2012.2018$prop.Cheats <- individuals2012.2018$noFCheatsCamp / individuals2012.2018$noFActive
#
#
#******#
### Summary data by camp----
#******#
#
#Create a df with one row per camp
camps2012.2018 <- treatments2012.2018 %>% distinct(yearcamp, .keep_all = T)
#
#Make yearcamp the first column
camps2012.2018 <- camps2012.2018 %>% dplyr::select(yearcamp, everything())
#
#Keep only relevant columns
colnamescamps <- data.frame(names = names(camps2012.2018), number = 1:length(names(camps2012.2018)))
colnamescamps
camps2012.2018 <- camps2012.2018[c("yearcamp", "year", "camp", "grpsize", "nomales", "nofemales", "eggsleftDate", "Removal.replacement.score")]
#
#Remove camps with removal/replacement score of 4
camps2012.2018 <- camps2012.2018[which(camps2012.2018$Removal.replacement.score != 4 | is.na(camps2012.2018$Removal.replacement.score)),]
#
####Nest number----
#Create a column with info on whether a camp had a nest in which sitting was recorded (1 = yes, 0 = no)
camps2012.2018$nestincamp <- "0"
camps2012.2018$nestincamp[camps2012.2018$yearcamp %in% sitduration2012.2018$yearcamp] <- "1"
table(camps2012.2018$nestincamp)#There are 5 camps with no nest in which incubation was recorded. Check these camps!
camps2012.2018[camps2012.2018$nestincamp == 0,]
#According to field notes and egg collection notes:
#-Camp 2012_29 had 1 nest
#-Camp 2012_30 had 1 nest
#-Camp 2012_68 had no nests
#-Camp 2013_55 had no nests
#-Camp 2016_30 had 1 nest
#Did these camps produce any eggs during the period with natural incubation?
sum(eggs2012.2018[eggs2012.2018$yearcamp == "2012_29" &
                    eggs2012.2018$eggs == "eggsleft",]$egglaid)#15 eggs were laid in camp 2012_29.
sum(eggs2012.2018[eggs2012.2018$yearcamp == "2012_30" &
                    eggs2012.2018$eggs == "eggsleft",]$egglaid)#14 eggs were laid in camp 2012_30.
sum(eggs2012.2018[eggs2012.2018$yearcamp == "2012_68" &
                    eggs2012.2018$eggs == "eggsleft",]$egglaid)#No eggs were laid in camp 2012_68.
sum(eggs2012.2018[eggs2012.2018$yearcamp == "2013_55" &
                    eggs2012.2018$eggs == "eggsleft",]$egglaid)#8 eggs were laid in camp 2013_55.
sum(eggs2012.2018[eggs2012.2018$yearcamp == "2016_30" &
                    eggs2012.2018$eggs == "eggsleft",]$egglaid)#8 eggs were laid in camp 2016_30.
#Record 1 nest for the camps which had a recorded nest that was active during the period with natural incubation.
#The eggs in camps 2013_55 and 2017_22, which had no recorded nest, are assumed to be loose eggs.
camps2012.2018[camps2012.2018$yearcamp == "2012_29" | camps2012.2018$yearcamp == "2012_30" |
                 camps2012.2018$yearcamp == "2016_30",]$nestincamp <- 1
#
#Create a column with info on whether a camp had a secondary nest
camps2012.2018$nest2incamp <- "0"
camps2012.2018$nest2incamp[camps2012.2018$yearcamp %in% sitduration2012.2018$yearcamp[which(sitduration2012.2018$nestrank == 2)]] <- "1"
#
#Create a column with info on whether a camp had a tertiary nest
camps2012.2018$nest3incamp <- "0"
camps2012.2018$nest3incamp[camps2012.2018$yearcamp %in% sitduration2012.2018$yearcamp[which(sitduration2012.2018$nestrank == 3)]] <- "1"
#
#Make nestincamp columns numeric
nestincamp <- c("nestincamp", "nest2incamp", "nest3incamp")
camps2012.2018 <- camps2012.2018 %>% dplyr::mutate_at(nestincamp, list(~ as.numeric(.)))
str(camps2012.2018)
#
#Compute number of nests per camp
camps2012.2018$nonests <- camps2012.2018$nestincamp + camps2012.2018$nest2incamp + camps2012.2018$nest3incamp
#
#### Egg production----
#####Period without care----
#Create a column with number of eggs laid in camp during period with no care
#Compute data
eggsCampNoCare <- summaryBy(egglaid ~ yearcamp, data = eggs2012.2018[which(eggs2012.2018$eggs == "eggsout"),], FUN = sum)
#Match data
camps2012.2018$no.eggsCampNoCare <- eggsCampNoCare$egglaid.sum[match(camps2012.2018$yearcamp, eggsCampNoCare$yearcamp)]
#
#####Period with care----
#Create a column with number of eggs laid in camp during period with care
#Compute data
eggsCampInc <- summaryBy(egglaid ~ yearcamp, data = eggs2012.2018[which(eggs2012.2018$eggs == "eggsleft"),], FUN = sum)
#Match data
camps2012.2018$no.eggsCampInc <- eggsCampInc$egglaid.sum[match(camps2012.2018$yearcamp, eggsCampInc$yearcamp)]
#Change NAs to zeros
camps2012.2018$no.eggsCampInc[is.na(camps2012.2018$no.eggsCampInc)] <- 0
#
#####Period with care and four weeks before----
#Compute data
eggsCampInc4 <- summaryBy(egglaid ~ yearcamp, data = eggs2012.2018[which(eggs2012.2018$daysfromEggsleft > -27),], FUN = sum)
#Match data
camps2012.2018$no.eggsCampInc4 <- eggsCampInc4$egglaid.sum[match(camps2012.2018$yearcamp, eggsCampInc4$yearcamp)]
#
#### Broken eggs----
#####Period with care----
#Create a column with number of eggs broken in camps during period with care
#Compute data
brokenCampInc <- summaryBy(egglaid ~ yearcamp, data = eggs2012.2018[which(eggs2012.2018$e_hatchoutcome == 12 &
                                                                            eggs2012.2018$eggs == "eggsleft"),],
                           FUN = sum)
#Match data
camps2012.2018$no.BrokenInc <- brokenCampInc$egglaid.sum[match(camps2012.2018$yearcamp, brokenCampInc$yearcamp)]
#Change NAs to zeros
camps2012.2018$no.BrokenInc[which(is.na(camps2012.2018$no.BrokenInc))] <- 0
#
#Compute proportion of broken eggs
camps2012.2018$prop.BrokenInc <- camps2012.2018$no.BrokenInc / camps2012.2018$no.eggsCampInc
#
#### Hatching----
#####Period with care----
#Create a column with number of eggs hatched in camp during period with care
#Compute data
hatchedCampInc <- summaryBy(c_hatched ~ yearcamp, data = eggs2012.2018[which(!is.na(eggs2012.2018$c_hatched) &
                                                                               eggs2012.2018$eggs == "eggsleft"),],
                            FUN = sum)
#Match data
camps2012.2018$no.HatchedInc <- hatchedCampInc$c_hatched.sum[match(camps2012.2018$yearcamp, hatchedCampInc$yearcamp)]
#
#####Period with care and four weeks before----
#Compute data
hatchedCampInc4 <- summaryBy(c_hatched ~ yearcamp, data = eggs2012.2018[which(eggs2012.2018$daysfromEggsleft > -27),], FUN = sum)
#Match data
camps2012.2018$no.HatchedInc4 <- hatchedCampInc4$c_hatched.sum[match(camps2012.2018$yearcamp, hatchedCampInc4$yearcamp)]
#Change NAs to zeros
camps2012.2018$no.HatchedInc4[is.na(camps2012.2018$no.HatchedInc4)] <- 0
#
#### Genotyped eggs----
#####Period with care----
GenoCampInc <- as.data.frame(table(parentage2012.2018$yearcamp[parentage2012.2018$eggsleft == "Yes"]))
#Match total number of eggs with genotype per camp during the period with care to the camps df
camps2012.2018$totGenoeggsCampInc <- GenoCampInc$Freq[match(camps2012.2018$yearcamp, GenoCampInc$Var1)]
#
#####Period with care and four weeks before----
GenoCampInc4 <- as.data.frame(table(parentage2012.2018$yearcamp[parentage2012.2018$daysfromEggsleft > -27]))
#Match total number of eggs with genotype per camp to the camps df
camps2012.2018$totGenoeggsCampInc4 <- GenoCampInc4$Freq[match(camps2012.2018$yearcamp, GenoCampInc4$Var1)]
#Change NAs to zeros
camps2012.2018$totGenoeggsCampInc4[which(is.na(camps2012.2018$totGenoeggsCampInc4))] <- 0
#
####Incubation----
#####Total----
#Create a column with total incubation duration in each camp
#Compute the data
sitpercamp <- summaryBy(sitduration ~ yearcamp, data = sitduration2012.2018, FUN = sum)
#Create a column where to insert the data, giving NA as default value
camps2012.2018$totsit <- NA
#Insert data to camp df
camps2012.2018$totsit <- sitpercamp$sitduration.sum[match(camps2012.2018$yearcamp, sitpercamp$yearcamp)]
#Convert NAs to 0s
camps2012.2018[is.na(camps2012.2018$totsit),]$totsit <- 0
#
#####Females----
#Create a column with total incubation duration by females in each camp
#Compute the data
Fsitpercamp <- summaryBy(sitduration ~ yearcamp, data = sitduration2012.2018[sitduration2012.2018$sex == "F",], FUN = sum)
#Create a column where to insert the data, giving NA as default value
camps2012.2018$Ftotsit <- NA
#Insert data to camp df
camps2012.2018$Ftotsit <- Fsitpercamp$sitduration.sum[match(camps2012.2018$yearcamp, Fsitpercamp$yearcamp)]
#Convert NAs to 0s
camps2012.2018[is.na(camps2012.2018$Ftotsit),]$Ftotsit <- 0
#
#####Proportion----
######Total----
#Create a column with proportion of time incubating
#Insert data on total observation time per year
camps2012.2018$totobstime <- sitduration2012.2018$obsdurationMinYear[match(camps2012.2018$year, sitduration2012.2018$year)]
#Compute proportion of time incubating
camps2012.2018$propsit <- camps2012.2018$totsit / camps2012.2018$totobstime
#
######Females----
#Create a column with proportion of time incubating by females
#Compute proportion of time incubating by females
camps2012.2018$Fpropsit <- camps2012.2018$Ftotsit / camps2012.2018$totobstime
#
####Copulation----
#Compute total number of successful copulations in each group during period with care and 28 days before
succCopsCampInc <- as.data.frame(table(copulation2012.2018$yearcamp[copulation2012.2018$succ == "Y" & copulation2012.2018$cop.period == "P2"]))
#Match data to camps df
camps2012.2018$succCops_P2 <- succCopsCampInc$Freq[match(camps2012.2018$yearcamp, succCopsCampInc$Var1)]
#Change NAs to zeros
camps2012.2018$succCops_P2[which(is.na(camps2012.2018$succCops_P2))] <- NA
#
#
####Breeding scores----
#Match number of inds. of a given breeding strategy per group (using RS during period with care and four weeks before)
camps2012.2018$noFCheatsCamp <- individuals2012.2018$noFCheatsCamp[match(camps2012.2018$yearcamp, individuals2012.2018$yearcamp)]
camps2012.2018$noFCoopsCamp <- individuals2012.2018$noFCoopsCamp[match(camps2012.2018$yearcamp, individuals2012.2018$yearcamp)]
camps2012.2018$noFOtherCamp <- individuals2012.2018$noFOtherCamp[match(camps2012.2018$yearcamp, individuals2012.2018$yearcamp)]
#Convert NAs to zeros
camps2012.2018$noFCheatsCamp[is.na(camps2012.2018$noFCheatsCamp)] <- 0
camps2012.2018$noFCoopsCamp[is.na(camps2012.2018$noFCoopsCamp)] <- 0
camps2012.2018$noFOtherCamp[is.na(camps2012.2018$noFOtherCamp)] <- 0
#
#Create column with information on whether there is a female cheat in camp or not
camps2012.2018$FCheatinCamp <- ifelse(camps2012.2018$noFCheatsCamp == 0, "No", "Yes")
#
######Number of active females----
#Compute the number of females in breeding group (cooperators + cheats)
camps2012.2018$noFActive <- camps2012.2018$noFCoopsCamp + camps2012.2018$noFCheatsCamp
#
######Proportion of cheats----
#Compute proportion of cheats
camps2012.2018$prop.Cheats <- camps2012.2018$noFCheatsCamp / camps2012.2018$noFActive
#
#
#******#
### Female-female aggression----
#******#
#
#Select relevant columns
names(femfem2012.2018)
Femfem <- femfem2012.2018[c("year", "date", "obstime", "camp", "tagAgg", "codeAgg", "idAgg", "tagVic", "codeVic", "idVic", "aggscore",
                            "vicscore", "obstime.minutes", "yearcamp", "nomales", "nofemales", "dayorderStart", "dayorderEggs", "eggs")]
#
#Compute yearID
Femfem$yearIDagg <- paste(Femfem$year, Femfem$idAgg, sep = "_")
Femfem$yearIDvic <- paste(Femfem$year, Femfem$idVic, sep = "_")
#
#Match breeding strategy
Femfem$BSagg <- individuals2012.2018$BS[match(Femfem$yearIDagg, individuals2012.2018$yearID)]
Femfem$BSvic <- individuals2012.2018$BS[match(Femfem$yearIDvic, individuals2012.2018$yearID)]
#Check
table(Femfem$BSagg, exclude = NULL)#There are 10 NAs in the BSagg-column
Femfem[which(is.na(Femfem$BSagg)),]#All individuals are from camps 2015_25, 2017_22 or 2018_25.
#Camps 2015_25 & 2017_22 were excluded from the data because of late removals of individuals. Exclude camps from the femfem data as well.
Femfem <- Femfem[Femfem$yearcamp != "2015_25" & Femfem$yearcamp != "2017_22",]
#Check camp 2018_25
treatments2012.2018[treatments2012.2018$yearcamp == "2018_25",]#The individual with the NA died before the end of the season and was never
#assigned a breeding strategy.
#Get rid of the interactions
Femfem <- Femfem[which(!is.na(Femfem$BSagg)),]
table(Femfem$BSagg, exclude = NULL)#Good!
#Check BSvic
table(Femfem$BSvic, exclude = NULL)#There is 1 NA
Femfem[which(is.na(Femfem$BSvic)),]#Its is the same individual that appeared as aggressor above.
#Get rid of interaction
Femfem <- Femfem[which(!is.na(Femfem$BSvic)),]
table(Femfem$BSvic, exclude = NULL)#Good!
#
#Match data on whether there was a cheat in the camp or not
Femfem$FCheatincamp <- camps2012.2018$FCheatinCamp[match(Femfem$yearcamp, camps2012.2018$yearcamp)]
#
#Create a df with per individual data on female-female aggression
#The only aggressor behaviour that is not defined as aggressive is clucking (0). Get rid of such interactions
table(Femfem$aggscore, exclude = NULL)
Femfem_aggr <- Femfem[which(Femfem$aggscore != 0),]
#
#
####Period with care and four weeks before----
#Create a df with per individual data on aggression during period with care and four weeks before
#subset data to only include period with care and four weeks before
Femfem_aggr_care4 <- Femfem_aggr[which(Femfem_aggr$dayorderEggs > -29),]
#Create a df with number of interactions per individual, where the individual is the aggressor
Femfem_agg_count_care4 <-
  summaryBy(aggscore ~ yearIDagg + BSagg + yearcamp + year, data = Femfem_aggr_care4, FUN = length)
Femfem_agg_count_care4 <- Femfem_agg_count_care4 %>% dplyr::rename(Tot.no.ints = "aggscore.length")
#Compute data set with number of aggressive interactions by breeding strategy of the victim
Femfem_agg_count_care4_BSvic <-
  summaryBy(aggscore ~ yearIDagg + BSagg + BSvic + yearcamp + year, data = Femfem_aggr_care4, FUN = length)
#
#Create a df with a row per individual, including individuals that were not seen interacting
#Create a list of all females present in camps during period with care
Femfem_agg_count_care4_BSvic_all <- 
  as.data.frame(treatments2012.2018$yearID[which(treatments2012.2018$sex == "F" & treatments2012.2018$incampeggsleft == "Yes")])
#Rename column
names(Femfem_agg_count_care4_BSvic_all)
Femfem_agg_count_care4_BSvic_all <- 
  Femfem_agg_count_care4_BSvic_all %>%
  dplyr::rename(yearIDagg = "treatments2012.2018$yearID[which(treatments2012.2018$sex == \"F\" & treatments2012.2018$incampeggsleft == \"Yes\")]")
#Add relevant columns
#BSagg
Femfem_agg_count_care4_BSvic_all$BSagg <-
  individuals2012.2018$BS[match(Femfem_agg_count_care4_BSvic_all$yearID, individuals2012.2018$yearID)]
table(Femfem_agg_count_care4_BSvic_all$BSagg, exclude = NULL)#There are 9 NAs. Check!
#yearcamp
Femfem_agg_count_care4_BSvic_all$yearcamp <-
  treatments2012.2018$yearcamp[match(Femfem_agg_count_care4_BSvic_all$yearID, treatments2012.2018$yearID)]
#Check the NAs for BSagg
Femfem_agg_count_care4_BSvic_all[which(is.na(Femfem_agg_count_care4_BSvic_all$BSagg)),]#All from excluded camps. Remove!
Femfem_agg_count_care4_BSvic_all <- Femfem_agg_count_care4_BSvic_all[which(!is.na(Femfem_agg_count_care4_BSvic_all$BSagg)),]
#year
Femfem_agg_count_care4_BSvic_all$year <-
  treatments2012.2018$year[match(Femfem_agg_count_care4_BSvic_all$yearID, treatments2012.2018$yearID)]
#Tot.no.ints
Femfem_agg_count_care4_BSvic_all$Tot.no.ints <- 
  Femfem_agg_count_care4$Tot.no.ints[match(Femfem_agg_count_care4_BSvic_all$yearID, Femfem_agg_count_care4$yearID)]
#Change Nas to 0 in the Tot.no.ints-column
Femfem_agg_count_care4_BSvic_all$Tot.no.ints <- ifelse(is.na(Femfem_agg_count_care4_BSvic_all$Tot.no.ints), 0,
                                                       Femfem_agg_count_care4_BSvic_all$Tot.no.ints)
#no.intsCooperators
Femfem_agg_count_care4_BSvic_all$no.intsCooperators <- 
  Femfem_agg_count_care4_BSvic$aggscore.length[which(Femfem_agg_count_care4_BSvic$BSvic == "Cooperator")][match(Femfem_agg_count_care4_BSvic_all$yearIDagg,
                                                                                                                Femfem_agg_count_care4_BSvic$yearIDagg[which(Femfem_agg_count_care4_BSvic$BSvic == "Cooperator")])]
#Change NAs to 0s
Femfem_agg_count_care4_BSvic_all$no.intsCooperators <- ifelse(is.na(Femfem_agg_count_care4_BSvic_all$no.intsCooperators), 0,
                                                              Femfem_agg_count_care4_BSvic_all$no.intsCooperators)
#no.intsCheats
Femfem_agg_count_care4_BSvic_all$no.intsCheats <-
  Femfem_agg_count_care4_BSvic$aggscore.length[which(Femfem_agg_count_care4_BSvic$BSvic == "Cheat")][match(Femfem_agg_count_care4_BSvic_all$yearIDagg,
                                                                                                           Femfem_agg_count_care4_BSvic$yearIDagg[which(Femfem_agg_count_care4_BSvic$BSvic == "Cheat")])]
#Change NAs to 0s
Femfem_agg_count_care4_BSvic_all$no.intsCheats <- ifelse(is.na(Femfem_agg_count_care4_BSvic_all$no.intsCheats), 0,
                                                         Femfem_agg_count_care4_BSvic_all$no.intsCheats)
#no.intsOther
Femfem_agg_count_care4_BSvic_all$no.intsOther <- 
  Femfem_agg_count_care4_BSvic$aggscore.length[which(Femfem_agg_count_care4_BSvic$BSvic == "Other")][match(Femfem_agg_count_care4_BSvic_all$yearIDagg,
                                                                                                           Femfem_agg_count_care4_BSvic$yearIDagg[which(Femfem_agg_count_care4_BSvic$BSvic == "Other")])]
#Change NAs to 0s
Femfem_agg_count_care4_BSvic_all$no.intsOther <- ifelse(is.na(Femfem_agg_count_care4_BSvic_all$no.intsOther), 0,
                                                        Femfem_agg_count_care4_BSvic_all$no.intsOther)

#
#add information on number of females with a given BS
Femfem_agg_count_care4_BSvic_all$no.FCheats <-
  individuals2012.2018$noFCheatsCamp[match(Femfem_agg_count_care4_BSvic_all$yearcamp, individuals2012.2018$yearcamp)]
Femfem_agg_count_care4_BSvic_all$no.FCoops <-
  individuals2012.2018$noFCoopsCamp[match(Femfem_agg_count_care4_BSvic_all$yearcamp, individuals2012.2018$yearcamp)]
Femfem_agg_count_care4_BSvic_all$no.FOther <-
  individuals2012.2018$noFOtherCamp[match(Femfem_agg_count_care4_BSvic_all$yearcamp, individuals2012.2018$yearcamp)]
#
#Compute number of interactions corrected for number of individuals of a given BS
Femfem_agg_count_care4_BSvic_all$no.intsCooperators.corr <- 
  Femfem_agg_count_care4_BSvic_all$no.intsCooperators / Femfem_agg_count_care4_BSvic_all$no.FCoops
Femfem_agg_count_care4_BSvic_all$no.intsCheats.corr <- 
  Femfem_agg_count_care4_BSvic_all$no.intsCheats / Femfem_agg_count_care4_BSvic_all$no.FCheats
Femfem_agg_count_care4_BSvic_all$no.intsOther.corr <- 
  Femfem_agg_count_care4_BSvic_all$no.intsOther / Femfem_agg_count_care4_BSvic_all$no.FOther
#
#Add info on group size
Femfem_agg_count_care4_BSvic_all$nomales <-
  treatments2012.2018$nomales[match(Femfem_agg_count_care4_BSvic_all$yearcamp, treatments2012.2018$yearcamp)]
Femfem_agg_count_care4_BSvic_all$nofemales <-
  treatments2012.2018$nofemales[match(Femfem_agg_count_care4_BSvic_all$yearcamp, treatments2012.2018$yearcamp)]
#Add info on total observation time
#Compute this data from the male-female interaction data.
#Create a unique session id
malefem2012.2018$session <- paste(malefem2012.2018$date, malefem2012.2018$timestart, sep = "_")
#Create a list of sessions
obsSessions2012_2018 <- as.data.frame(unique(malefem2012.2018$session))
obsSessions2012_2018 <- obsSessions2012_2018 %>% dplyr::rename(session = "unique(malefem2012.2018$session)")
#Add data on session duration
obsSessions2012_2018$obsduration <-
  malefem2012.2018$obsduration[match(obsSessions2012_2018$session, malefem2012.2018$session)]
#Add data on year
obsSessions2012_2018$year <-
  malefem2012.2018$year[match(obsSessions2012_2018$session, malefem2012.2018$session)]
#Add data on the day order from the onset of the period with care
obsSessions2012_2018$dayorderEggs <-
  malefem2012.2018$dayorderEggs[match(obsSessions2012_2018$session, malefem2012.2018$session)]
#Include only the period with care and four weeks before
obsSessions2012_2018_care4 <- obsSessions2012_2018[obsSessions2012_2018$dayorderEggs > -29,]
#Summarize data per year
Obsdur_year_care4 <- summaryBy(obsduration ~ year, data = obsSessions2012_2018_care4, FUN = sum)
#Match data on total observation time during period with care
Femfem_agg_count_care4_BSvic_all$totObstime <-
  Obsdur_year_care4$obsduration.sum[match(Femfem_agg_count_care4_BSvic_all$year, Obsdur_year_care4$year)]
#
#
#Create a long df with aggression during period with care by BS of victim
#Stack victim BS
names(Femfem_agg_count_care4_BSvic_all)
Femfem_agg_count_care4_long <- cbind(Femfem_agg_count_care4_BSvic_all[1:5], stack(Femfem_agg_count_care4_BSvic_all[6:8]), row.names = NULL)
#Change names
names(Femfem_agg_count_care4_long)
Femfem_agg_count_care4_long <- Femfem_agg_count_care4_long %>% dplyr::rename(no.ints = 6, BSvic = 7)
#Change BSvic column
Femfem_agg_count_care4_long$BSvic <- gsub("no.intsCooperators", "Cooperator", Femfem_agg_count_care4_long$BSvic)
Femfem_agg_count_care4_long$BSvic <- gsub("no.intsCheats", "Cheat", Femfem_agg_count_care4_long$BSvic)
Femfem_agg_count_care4_long$BSvic <- gsub("no.intsOther", "Other", Femfem_agg_count_care4_long$BSvic)
#Compute number of interactions corrected for number of inds with a given BS
#add information on number of females with a given BS
Femfem_agg_count_care4_long$no.FBS <-
  ifelse(Femfem_agg_count_care4_long$BSvic == "Cooperator",
         individuals2012.2018$noFCoopsCamp[match(Femfem_agg_count_care4_long$yearcamp,
                                                 individuals2012.2018$yearcamp)],
         ifelse(Femfem_agg_count_care4_long$BSvic == "Cheat",
                individuals2012.2018$noFCheatsCamp[match(Femfem_agg_count_care4_long$yearcamp,
                                                         individuals2012.2018$yearcamp)],
                individuals2012.2018$noFOtherCamp[match(Femfem_agg_count_care4_long$yearcamp,
                                                        individuals2012.2018$yearcamp)]))
#Compute corrected no. interactions
Femfem_agg_count_care4_long$no.ints.corr <- Femfem_agg_count_care4_long$no.ints / Femfem_agg_count_care4_long$no.FBS
#
#Add info on group size
Femfem_agg_count_care4_long$nomales <-
  treatments2012.2018$nomales[match(Femfem_agg_count_care4_long$yearcamp, treatments2012.2018$yearcamp)]
Femfem_agg_count_care4_long$nofemales <-
  treatments2012.2018$nofemales[match(Femfem_agg_count_care4_long$yearcamp, treatments2012.2018$yearcamp)]
#
#Add info on whether there were cheats in the group or not
Femfem_agg_count_care4_long$cheatinCamp <- camps2012.2018$FCheatinCamp[match(Femfem_agg_count_care4_long$yearcamp,
                                                                             camps2012.2018$yearcamp)]
#
#
#####Cooperator aggression----
#Create a df with cooperator aggression data summarized per group
#Subset data to only include interactions where the aggressor is a cooperator
Femfem_coop_agg_count_care4 <- Femfem_agg_count_care4_BSvic_all[which(Femfem_agg_count_care4_BSvic_all$BSagg == "Cooperator"),]
#Compute per camp data
Femfem_coop_agg_count_care4_grp <-
  summaryBy(Tot.no.ints + no.intsCooperators + no.intsCheats + no.intsOther ~ yearcamp + year + no.FCheats +
              no.FCoops + no.FOther + nomales + nofemales + totObstime, data = Femfem_coop_agg_count_care4, FUN = sum)
#Compute number of interactions corrected for number of inds with a given BS
Femfem_coop_agg_count_care4_grp$no.intsCooperators.corr <- Femfem_coop_agg_count_care4_grp$no.intsCooperators.sum /
  Femfem_coop_agg_count_care4_grp$no.FCoops
Femfem_coop_agg_count_care4_grp$no.intsCheats.corr <- Femfem_coop_agg_count_care4_grp$no.intsCheats.sum /
  Femfem_coop_agg_count_care4_grp$no.FCheats
Femfem_coop_agg_count_care4_grp$no.intsother.corr <- Femfem_coop_agg_count_care4_grp$no.intsOther.sum /
  Femfem_coop_agg_count_care4_grp$no.FOther
#
#Create a long group level df with aggression during period with care by BS of victim
#Stack victim BS and number of interactions per BS
names(Femfem_coop_agg_count_care4_grp)
Femfem_coop_agg_count_care4_grp_long <- cbind(Femfem_coop_agg_count_care4_grp[1:9],
                                              stack(Femfem_coop_agg_count_care4_grp[10:12]),
                                              stack(Femfem_coop_agg_count_care4_grp[13:15]), row.names = NULL)
head(Femfem_coop_agg_count_care4_grp_long)
#Change names
names(Femfem_coop_agg_count_care4_grp_long)
Femfem_coop_agg_count_care4_grp_long <- Femfem_coop_agg_count_care4_grp_long %>% dplyr::rename(no.ints = 10, BSvic = 11,
                                                                                               no.ints.corr = 12)
head(Femfem_coop_agg_count_care4_grp_long)
#Delete redundant column
Femfem_coop_agg_count_care4_grp_long$ind <- NULL
#
#Change BSvic column
Femfem_coop_agg_count_care4_grp_long$BSvic <- gsub("no.intsCooperators.sum", "Cooperator", Femfem_coop_agg_count_care4_grp_long$BSvic)
Femfem_coop_agg_count_care4_grp_long$BSvic <- gsub("no.intsCheats.sum", "Cheat", Femfem_coop_agg_count_care4_grp_long$BSvic)
Femfem_coop_agg_count_care4_grp_long$BSvic <- gsub("no.intsOther.sum", "Other", Femfem_coop_agg_count_care4_grp_long$BSvic)
#
#Match egg data
Femfem_coop_agg_count_care4_grp_long$totGenoeggsCampInc4 <- 
  camps2012.2018$totGenoeggsCampInc4[match(Femfem_coop_agg_count_care4_grp_long$yearcamp, camps2012.2018$yearcamp)]
Femfem_coop_agg_count_care4_grp_long$no.eggsCampInc <- 
  camps2012.2018$no.eggsCampInc[match(Femfem_coop_agg_count_care4_grp_long$yearcamp, camps2012.2018$yearcamp)]
#
#Match incubation data
Femfem_coop_agg_count_care4_grp_long$Ftotsit <- 
  camps2012.2018$Ftotsit[match(Femfem_coop_agg_count_care4_grp_long$yearcamp, camps2012.2018$yearcamp)]
#
#
#******#
###Incubation and Egg laying across time#####----
#******#
#
####Nest-level incubation----
#****************************
#Data of incubation time per nest and session. Accounting for overlapping incubation records.
#
#Merge overlapping time periods at nest level (when more than one bird has been recorded sitting on the same nest at the same time)
#There are some incubations that are only one minute long. These could be females at the nest laying eggs, or soliciting. Get rid!
sittimesNest <- sitduration2012.2018[sitduration2012.2018$sitduration != 1,] %>%
  #Order the observations chronologically (with incubation start as reference) for every nest in every observation session
  dplyr::arrange(sessionNest, timesitCorr.min) %>%
  #Group together all observations from the same nest and observation session
  dplyr::group_by(sessionNest) %>%
  #Add an index variable. If two incubation bouts overlap temporally, they are given the same index value.
  #The first observation in the ordered and grouped vector is given a 0.
  dplyr::mutate(indx = c(0, cumsum(as.numeric(lead(timesitCorr.min)) > #This line computes the cumulative sum for the start of the next observation
                                     cummax(as.numeric(timefinishCorr.min)))[-n()])) %>% #This line computes the cumulative maximum value for the end of an
  #incubation. [-n()] removes the last item in the cumsum to accommodate
  #the added 0 at the beginning of the vector.
  #Group by nest, observation session and index
  group_by(sessionNest, indx) %>%
  #Merge overlapping incubations. 
  dplyr::summarise(timesitCorr.min = dplyr::first(timesitCorr.min),
                   timefinishCorr.min = dplyr::last(timefinishCorr.min), .groups = 'drop') %>%
  #Duration per collapsed (merged) nest-sitting event
  dplyr::mutate(sitdurationNest = timefinishCorr.min - timesitCorr.min)
#
#Compute incubations at nest level by sex
sittimesNest_sex <- sitduration2012.2018[sitduration2012.2018$sitduration != 1,] %>%
  #Order the observations chronologically (with incubation start as reference) for every nest in every observation session
  dplyr::arrange(sessionNest, timesitCorr.min) %>%
  #Group together all observations from the same nest and observation session
  dplyr::group_by(sessionNest, sex) %>%
  #Add an index variable. If two incubation bouts overlap temporally, they are given the same index value.
  #The first observation in the ordered and grouped vector is given a 0.
  dplyr::mutate(indx = c(0, cumsum(as.numeric(lead(timesitCorr.min)) > #This line computes the cumulative sum for the start of the next observation
                                     cummax(as.numeric(timefinishCorr.min)))[-n()])) %>% #This line computes the cumulative maximum value for the end of an
  #incubation. [-n()] removes the last item in the cumsum to accommodate
  #the added 0 at the beginning of the vector.
  #Group by nest, observation session and index
  group_by(sessionNest, sex, indx) %>%
  #Merge overlapping incubations. 
  dplyr::summarise(timesitCorr.min = dplyr::first(timesitCorr.min),
                   timefinishCorr.min = dplyr::last(timefinishCorr.min), .groups = 'drop') %>%
  #Duration per collapsed (merged) nest-sitting event
  dplyr::mutate(sitdurationNest = timefinishCorr.min - timesitCorr.min)
#
#Compute sitduration per nest and session and add relevant columns
incNest <- sittimesNest %>%
  group_by(sessionNest) %>%
  dplyr::summarise(sitdurationNestSess = sum(sitdurationNest), .groups = "drop") %>%
  left_join(sitduration2012.2018 %>%
              distinct(sessionNest, .keep_all = TRUE) %>%
              dplyr::select(sessionNest, session, nest, yearcampNest, yearcamp, camp, year,
                            date, timestart, timeend, obsduration,
                            nestrank, nomales, nofemales),
            by = "sessionNest")
#
#Compute sitduration per nest, sex and session and add relevant columns
incNest_sex <- sittimesNest_sex %>%
  group_by(sessionNest, sex) %>%
  dplyr::summarise(sitdurationNestSess_sex = sum(sitdurationNest), .groups = "drop") %>%
  left_join(sitduration2012.2018 %>%
              distinct(sessionNest, .keep_all = TRUE) %>%
              dplyr::select(sessionNest, session, nest, yearcampNest, yearcamp, camp, year,
                            date, timestart, timeend, obsduration,
                            nestrank, nomales, nofemales),
            by = "sessionNest")
#
#Compute female incubation in primary nests
incNest_F <- incNest_sex[incNest_sex$sex == "F" & incNest_sex$nestrank == 1,] %>%
  group_by(yearcamp) %>%
  dplyr::summarise(Finc_N1 = sum(sitdurationNestSess_sex), totObstime = sum(obsduration), .groups = "drop")
#
#Compute number of days since eggs were left
incNest <- incNest %>%
  #Add eggs-left date from treatments data
  left_join(treatments2012.2018 %>% dplyr::select(year, eggsleftDate) %>% distinct(),
            by = "year") %>%
  #Compute number of days since eggs left
  dplyr::mutate(dayordereggsleft = as.numeric(date - eggsleftDate)) %>%
  #Clean up unnecessary column
  dplyr::select(-eggsleftDate)
#

#Compute total observation time in a given date
#Create a data frame with unique sessions
UniqueSessions <- incNest %>% distinct(session, .keep_all = T)
#Compute the data
totobsdurDate <- summaryBy(obsduration ~ date, data = UniqueSessions, FUN = sum)
#Add data to sittimesNest df
incNest$totobsdurDate <- totobsdurDate$obsduration.sum[match(incNest$date, totobsdurDate$date)]
#
#
####Summary data at nest level----
#*********************************
#
#Create a new df with one row per nest
nests <- sitduration2012.2018 %>%
  distinct(yearcampNest, .keep_all = TRUE) %>%#Keep one row per nest
  dplyr::select(yearcampNest, year, camp, nest, yearcamp, nestrank,
                grpsize, nofemales, nomales)#Keep only relevant columns
#
#Add nest that are not detected by the incubation data.
#The nests where there was no incubation were:
#-Camp 2012_29, N1
#-Camp 2012_30, N1
#-Camp 2016_30, N1
missing_nests <- tibble::tibble(yearcampNest = c("2012_29N1", "2012_30N1", "2016_30N1"),
                                yearcamp = c("2012_29", "2012_30", "2016_30"),
                                nest = "N1",
                                nestrank = 1,
                                year = c(2012, 2012, 2016),
                                camp = c("29", "30", "30"))
#Add missing nests
nests <- bind_rows(nests, missing_nests)
#Add group size data by joining treatments info
nests <- nests %>%
  left_join(treatments2012.2018 %>%
              distinct(yearcamp, grpsize, nofemales, nomales), 
            by = "yearcamp", suffix = c("", ".treat")) %>%
  dplyr::mutate(grpsize = coalesce(grpsize, grpsize.treat),
                nofemales = coalesce(nofemales, nofemales.treat),
                nomales = coalesce(nomales, nomales.treat)) %>%
  dplyr::select(-ends_with(".treat"))#Delete redundant columns
#
#Compute total sitting duration per nest-year
totsitNestYear <- incNest %>%
  group_by(yearcampNest) %>%
  dplyr::summarise(totsitNestYear = sum(sitdurationNestSess), .groups = "drop")
#Add data to nest dataset, set NA to 0 for nests without incubation
nests <- nests %>%
  left_join(totsitNestYear, by = "yearcampNest") %>%
  dplyr::mutate(totsitNestYear = replace_na(totsitNestYear, 0))
#
#Compute total observation time per year
totobsdurYear <- sitduration2012.2018 %>%
  distinct(session, .keep_all = TRUE) %>%
  group_by(year) %>%
  dplyr::summarise(totobsdurYear = sum(obsduration), .groups = "drop")
#Add to nest data
nests <- nests %>%
  left_join(totobsdurYear, by = "year")
#
#Compute proportion of time nests where protected
nests <- nests %>%
  dplyr::mutate(propsitNestYear = totsitNestYear / totobsdurYear)
#
#
####Summary data at nest level - during time with natural incubation----
#***********************************************************************
#
#Create a new df with one row per nest
nestsInc <- nests %>%
  dplyr::select(yearcampNest, year, camp, nest, yearcamp, nestrank,
                grpsize, nofemales, nomales)#Keep only relevant columns
#
#Compute sitting duration per nest-year during time with natural incubation
totsitNestInc <- incNest[incNest$dayordereggsleft >= 0,] %>%
  group_by(yearcampNest) %>%
  dplyr::summarise(totsitNestInc = sum(sitdurationNestSess), .groups = "drop")
#Add data to nest dataset, set NA to 0 for nests without incubation
nestsInc <- nestsInc %>%
  left_join(totsitNestInc, by = "yearcampNest") %>%
  dplyr::mutate(totsitNestInc = replace_na(totsitNestInc, 0))
#
#Compute observation time per year during time with natural incubation
totobsdurInc <- sitduration2012.2018[sitduration2012.2018$eggs == "eggsleft",] %>%
  distinct(session, .keep_all = TRUE) %>%
  group_by(year) %>%
  dplyr::summarise(totobsdurInc = sum(obsduration), .groups = "drop")
#Add to nest data
nestsInc <- nestsInc %>%
  left_join(totobsdurInc, by = "year")
#
#Compute proportion of time nests where protected
nestsInc <- nestsInc %>%
  dplyr::mutate(propsitNestInc = totsitNestInc / totobsdurInc)
#
#
####Nest-level incubation, including no-sits----
#***********************************************
#OBS! This data frame takes into account sessions when some incubation (in any nest) was observed.
#Sessions with no observations at all are not included here.
#
#For each nest, add the sessions when no observations were recorded in the nest as zero sittings
#Create somewhere to put the data
incNestZeros <- as.numeric()
#Create a df with unique sessions
UniqueSessions <- incNest %>% distinct(session, .keep_all = T)
#Select unique nests
YrCampNests <- unique(nests[c("yearcampNest", "year")])
#Run a loop to compute the data
for(i in YrCampNests$yearcampNest){
  #i = "2012_65N1" #Run this to check in case of errors.
  #Find the year in which the nest was observed
  Year <- YrCampNests[YrCampNests$yearcampNest == i,]$year
  #Find sessions in which nest i was guarded
  SessionPresentForNest <- unique(incNest$session[incNest$yearcampNest %in% i])
  #Find sessions in which nest i was not guarded
  SessionsMissing <- UniqueSessions[which(!UniqueSessions$session %in% SessionPresentForNest & UniqueSessions$year == Year),]
  #Convert session to character
  SessionsMissing$session <- as.character(SessionsMissing$session)
  #Add data to entries with no incubation
  if(length(SessionsMissing$session) > 0){
    OutNest <- data.frame(cbind("sitdurationNestSess" = 0, "session" = SessionsMissing$session, "yearcampNest" = i))
    #Add nest level data
    ToAdd <- merge(OutNest, unique(nests[,c("nest", "nestrank", "year", "camp", "yearcamp", "nomales", "nofemales",
                                            "yearcampNest")]),
                   by = "yearcampNest", all.x = TRUE)
    #Add session level data
    ToAdd <- merge(ToAdd, UniqueSessions[,c("timestart","timeend", "obsduration", "date", "session")], by = "session", all.x = TRUE)
    incNestZeros <- rbind(incNestZeros,ToAdd)
  }
}
#
#Create a df with incubation duration and records for when no incubation was observed in a given nest
#Add columns present in incNest but not present in incNestZeros and give the NA values
#Identify missing columns in incNestZeros
missing_cols <- setdiff(names(incNest), names(incNestZeros))
#Add missing columns as NA to `incNestZeros`, make sure columns are of the same type
str(incNest[, intersect(names(incNest), missing_cols), drop = FALSE])
for (col in missing_cols) {
  template <- incNest[[col]]
  incNestZeros[[col]] <- NA
  # Coerce to correct type (based on `incNest`)
  if (is.numeric(template)) {
    incNestZeros[[col]] <- as.numeric(NA)
  } else if (is.character(template)) {
    incNestZeros[[col]] <- as.character(NA)
  } 
}
#Make sure the rest of the columns are also of the same type in both dataframes
str(incNest)
for (col in names(incNest)) {
  if (is.Date(incNest[[col]])) {
    incNestZeros[[col]] <- as.Date(incNestZeros[[col]])
  } else if (is.numeric(incNest[[col]]) || is.integer(incNest[[col]])) {
    incNestZeros[[col]] <- as.numeric(incNestZeros[[col]])
  } else if (is.character(incNest[[col]])) {
    incNestZeros[[col]] <- as.character(incNestZeros[[col]])
  }
}
#
#Bind dataframes and compute/match relevant columns
incNestAll <- bind_rows(incNest, incNestZeros %>% dplyr::select(all_of(names(incNest)))) %>%
  mutate(sessionNest = paste(session, yearcampNest, sep = "_"),
         dayordereggsleft = incNest$dayordereggsleft[match(date, incNest$date)],
         totobsdurDate = incNest$totobsdurDate[match(date, incNest$date)])
#
#
####Incubation per day on primary nest----
#*****************************************
#Nests that the groups incubated most when eggs were left in the camps will be classified as primary
#
#Identify primary nests
#Current nestrank does not account for overlapping incubation, and includes data from when eggs were collected
#Function to select primary nest per group and warn in case of ties
select_primary_nests <- function(data) {
  data_ties <- data %>%
    group_by(yearcamp) %>%
    dplyr::mutate(max_duration = max(propsitNestInc, na.rm = TRUE),
           n_max = sum(propsitNestInc == max_duration)) %>%
    ungroup()
  
  tied_groups <- unique(data_ties$yearcamp[data_ties$n_max > 1])
  
  if (length(tied_groups) > 0) {
    warning(sprintf(
      "Tie detected in: %s",
      paste(tied_groups, collapse = ", ")
    ))
  }
  
  primary_nests <- data_ties %>%
    group_by(yearcamp) %>%
    dplyr::filter(propsitNestInc == max(propsitNestInc, na.rm = TRUE)) %>%
    slice(1) %>%# In case of tie, takes first
    ungroup() %>%
    dplyr::select(-max_duration, -n_max)
  
  return(primary_nests)
}
#Run function
primary_nests <- select_primary_nests(nestsInc)
table(primary_nests$nestrank)#The new ranking does not change anything. Use nestrank.
#
#Compute incubation duration per day for every group (in primary nest)
sitdur_date <- incNestAll[incNestAll$dayordereggsleft >= 0 & incNestAll$nestrank == 1,] %>%
  group_by(yearcamp, date) %>%                                   
  dplyr::summarise(sitdur_date = sum(sitdurationNestSess, na.rm = TRUE), .groups = "drop") %>% 
  left_join(incNestAll[incNestAll$dayordereggsleft >= 0,] %>% 
              distinct(session, .keep_all = TRUE) %>%                   
              group_by(date) %>% 
              dplyr::summarise(obsdur_date = sum(obsduration, na.rm = TRUE), .groups = "drop"), by = "date") %>% 
  left_join(incNestAll %>% 
              distinct(date, dayordereggsleft), by = "date") %>% 
  dplyr::mutate(propsit = sitdur_date / obsdur_date)
#
#
####Trimming data at end of season for consistency across years----
#******************************************************************
#
#Add information on year in the sitdur_date data frame
sitdur_date$year <- substr(sitdur_date$yearcamp, start = 1, stop = 4)
#
#Check season length per year (when was the latest incubation recorded?)
sitdur_date %>%
  group_by(year) %>%
  summarise(latest_day = max(dayordereggsleft, na.rm = TRUE))#The shortest year is 2018, with 21 days, the longest is 2013, with 60.
#
#Plot mean proportion of time incubating per year
#Summarize data:
propsit_date_sum_year <- sitdur_date %>%
  group_by(dayordereggsleft, year) %>% 
  dplyr::summarise(mean = mean(propsit, na.rm = T),
                   lcl = mean(propsit, na.rm = T) - (sd(propsit)/sqrt(n())),
                   ucl = mean(propsit, na.rm = T) + (sd(propsit)/sqrt(n())), .groups = "drop")
#
ggplot(sitdur_date,
       aes(x = dayordereggsleft, y = propsit)) +
  geom_point(alpha = 0.2) +
  geom_line(data = propsit_date_sum_year,
            aes(x = dayordereggsleft, y = mean), linewidth = 1, alpha = 2) +
  geom_pointrange(data = propsit_date_sum_year,
                  aes(x = dayordereggsleft, y = mean, ymin = lcl, ymax = ucl),
                  fatten = 5, size = 0.6, shape = 21, alpha = 1) +
  geom_vline(xintercept = c(21, 26, 30, 32, 33), linetype = "dotted") +
  labs(x = "Days from onset of experimental incubation",
       y = "Proportion of eggs laid") + facet_wrap(~ year) + theme_minimal()
#
#Plot mean proportion of time incubating across all years
#Summarize data:
propsit_date_sum <- sitdur_date %>%
  group_by(dayordereggsleft) %>% 
  dplyr::summarise(mean = mean(propsit, na.rm = T),
                   lcl = mean(propsit, na.rm = T) - (sd(propsit)/sqrt(n())),
                   ucl = mean(propsit, na.rm = T) + (sd(propsit)/sqrt(n())), .groups = "drop")
#
ggplot(sitdur_date,
       aes(x = dayordereggsleft, y = propsit)) +
  geom_point(alpha = 0.2) +
  geom_line(data = propsit_date_sum,
            aes(x = dayordereggsleft, y = mean), linewidth = 1, alpha = 2) +
  geom_pointrange(data = propsit_date_sum,
                  aes(x = dayordereggsleft, y = mean, ymin = lcl, ymax = ucl),
                  fatten = 5, size = 0.6, shape = 21, alpha = 1) +
  geom_vline(xintercept = c(21, 26, 30, 32, 33), linetype = "dotted") +
  labs(x = "Days from onset of experimental incubation",
       y = "Proportion of eggs laid") + theme_minimal()
#Error around the mean seems to increase at 33 days, and after 32 days, most years are not represented in the data.
#Filter the data to only include up to 32 days after onset of natural incubation.
#
sitdur_date_filtered <- sitdur_date[sitdur_date$dayordereggsleft < 33,]
#
####Incubation thresholds----
#****************************
#Here, I define when a group started incubating consistently.
#For this I use breakpoint analysis
#
#Write function to perform segmented regression on each group
seg_regression <- function(group_data) {
  #Check if we have enough data points
  if (nrow(group_data) < 6) {
    return(list(
      breakpoint = NA,
      mean_before = NA,
      mean_after = NA,
      change = NA,
      converged = FALSE,
      error = "Insufficient data"
    ))
  }
  #
  #Sort by day order
  group_data <- group_data %>% arrange(dayordereggsleft)
  #
  #Try segmented regression
  tryCatch({
    # Fit initial linear model
    lm_fit <- lm(propsit ~ dayordereggsleft, data = group_data)
    #
    #Try segmented regression with automatic breakpoint detection
    seg_fit <- segmented(lm_fit, seg.Z = ~dayordereggsleft, npsi = 1)
    #
    #Extract breakpoint
    breakpoint <- seg_fit$psi[1, "Est."]
    #
    #Calculate means before and after breakpoint
    mean_before <- mean(group_data$propsit[group_data$dayordereggsleft <= breakpoint], na.rm = TRUE)
    mean_after <- mean(group_data$propsit[group_data$dayordereggsleft > breakpoint], na.rm = TRUE)
    #
    return(list(
      breakpoint = round(breakpoint),
      mean_before = mean_before,
      mean_after = mean_after,
      change = mean_after - mean_before,
      converged = TRUE,
      error = NA
    ))
  }, error = function(e) {
    return(list(
      breakpoint = NA,
      mean_before = NA,
      mean_after = NA,
      change = NA,
      converged = FALSE,
      error = as.character(e)
    ))
  })
}
#
#Not all groups are relevant. Filter data to include only relevant groups
sitdur_date_filtered <- sitdur_date_filtered %>%
  left_join(camps2012.2018 %>%#Add relevant data
              dplyr::select(yearcamp, totGenoeggsCampInc4, no.eggsCampInc, nofemales), by = "yearcamp") %>%
  dplyr::filter(totGenoeggsCampInc4 != 0, no.eggsCampInc != 0, nofemales != 1) %>%#Filter data to match criteria.
  group_by(yearcamp) %>% #Filter also out groups with no incubation during the relevant period
  dplyr::filter(sum(propsit != 0) > 0) %>% ungroup()
#
#Apply segmented regression to each group
breakpoints <- sitdur_date_filtered[sitdur_date_filtered$dayordereggsleft >= 0,] %>%
  group_by(yearcamp) %>%
  nest() %>%
  dplyr::mutate(seg_breakpoints = map(data, seg_regression)) %>%
  unnest_wider(seg_breakpoints)#Warnings suggest that breakpoints are not supported in group 2013_2018.
#
#Filter successful breakpoints
successful_breakpoints <- breakpoints %>%
  dplyr::filter(converged == TRUE, !is.na(breakpoint))#Two groups were filtered out
#
#Summary statistics
yearly_stats <- successful_breakpoints %>%
  dplyr::mutate(year = as.integer(str_sub(yearcamp, 1, 4))) %>%
  group_by(year) %>%
  dplyr::summarise(n_groups = n(),
                   mean_breakpoint = mean(breakpoint, na.rm = TRUE),
                   sd_breakpoint = sd(breakpoint, na.rm = TRUE),
                   mean_change = mean(change, na.rm = TRUE),
                   .groups = "drop")
yearly_stats#There is strong variation in breakpoints across years
#
#Detailed results table
breakpoints_group <- successful_breakpoints %>%
  dplyr::select(yearcamp, breakpoint, mean_before, mean_after, change) %>%
  dplyr::mutate(
    mean_before = round(mean_before, 3),
    mean_after = round(mean_after, 3),
    change = round(change, 3)) %>%
  arrange(yearcamp)
#
#Visualize the data
#
#1. Histogram of breakpoints
breakpoint_hist <- successful_breakpoints %>%
  ggplot(aes(x = breakpoint)) +
  geom_histogram(bins = 15, fill = "steelblue", alpha = 0.7, color = "black") +
  geom_vline(xintercept = median(successful_breakpoints$breakpoint), 
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = "Distribution of Incubation Onset Days",
    subtitle = paste("Median =", median(successful_breakpoints$breakpoint), "days"),
    x = "Day when consistent incubation begins",
    y = "Number of groups"
  ) +
  theme_minimal()
breakpoint_hist
#
#2. Year-by-year boxplot
year_boxplot <- successful_breakpoints %>%
  dplyr::mutate(year = as.integer(str_sub(yearcamp, 1, 4))) %>%
  ggplot(aes(x = factor(year), y = breakpoint)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  labs(
    title = "Incubation Onset Timing by Year",
    x = "Year",
    y = "Day when consistent incubation begins") +
  theme_minimal()
year_boxplot
#
#3. Time series plots
#Create one multi-panel plot for each year. One panel per group. 
#create a list with years
years <- 2012:2018
#Create an empty list to store plot objects
multi_panel_plots <- list()
#Run loop to create plots
for (yr in years) {
  year_str <- as.character(yr)
  #
  # Filter data for the current year
  df_year <- sitdur_date_filtered %>%
    dplyr::filter(dayordereggsleft >= 0, grepl(year_str, yearcamp))
  breakpoints_year <- breakpoints_group %>% 
    dplyr::filter(grepl(year_str, yearcamp))
  #
  # Build the ggplot
  p <- ggplot(df_year, aes(x = dayordereggsleft, y = propsit)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_line(alpha = 0.7) +
    facet_wrap(~ yearcamp) +
    labs(title = paste("Proportion of Time Incubating by Group -", year_str),
         x = "Day in incubation period",
         y = "Proportion of time incubating") +
    theme_minimal()
  
  # Add breakpoints if any exist for this year
  if (nrow(breakpoints_year) > 0) {
    p <- p +
      geom_vline(data = breakpoints_year,
                 aes(xintercept = breakpoint),
                 color = "red", linetype = "dashed", linewidth = 1) +
      geom_text(data = breakpoints_year,
                aes(x = breakpoint + 2, y = 0.8, label = paste("Breakpoint:", breakpoint)),
                color = "red",
                inherit.aes = FALSE)
  }
  #
  # Store the plot in the list
  multi_panel_plots[[year_str]] <- p
}
#Visually control plots
multi_panel_plots[["2012"]]
multi_panel_plots[["2013"]]
multi_panel_plots[["2014"]]
multi_panel_plots[["2015"]]
multi_panel_plots[["2016"]]
multi_panel_plots[["2017"]]
multi_panel_plots[["2018"]]
#The breakpoints for many of the groups do not look reasonable
#Breakpoint analysis does not seem to be reliable at a group level. Try to do it at a population level.
#
#Population-level analysis
#Compute day-by-day data
population_data <- sitdur_date_filtered[sitdur_date_filtered$dayordereggsleft >= 0,] %>%
  group_by(dayordereggsleft) %>%
  dplyr::summarise(mean_propsit = mean(propsit, na.rm = TRUE),
                   se_propsit = sd(propsit, na.rm = TRUE) / sqrt(n()),
                   n_obs = n(), .groups = "drop")
#
#Fit population-level segmented regression
pop_lm <- lm(mean_propsit ~ dayordereggsleft, data = population_data)
pop_seg <- segmented(pop_lm, seg.Z = ~dayordereggsleft, npsi = 1)
summary(pop_seg)#Breakpoint at 15 days
slope(pop_seg)#The slope goes from positive, to flat
mean(population_data$mean_propsit[population_data$dayordereggsleft == 15])#Men inc. at breakpoint=76%
mean(population_data$mean_propsit[population_data$dayordereggsleft <= pop_seg$psi[1, "Est."]])#Mean inc. before=34%
mean(population_data$mean_propsit[population_data$dayordereggsleft >= pop_seg$psi[1, "Est."]])#Mean inc. efter=69%
min(population_data$mean_propsit[population_data$dayordereggsleft >= pop_seg$psi[1, "Est."]])#Min inc. after=59%
#
#Visualize the data
population_plot <- population_data %>%
  ggplot(aes(x = dayordereggsleft, y = mean_propsit)) +
  geom_point(alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_propsit - se_propsit, 
                    ymax = mean_propsit + se_propsit), 
                width = 0.5, alpha = 0.5) +
  geom_vline(xintercept = pop_seg$psi[1, "Est."], 
             color = "red", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 0.59, 
             color = "black", linetype = "dashed", linewidth = 0.7) +
  labs(title = "Population-Level Incubation Pattern",
       subtitle = paste("Population breakpoint:", round(pop_seg$psi[1, "Est."], 1), "days"),
       x = "Day in incubation period",
       y = "Mean proportion incubating") +
  theme_minimal()
population_plot#The population-level breakpoints looks reasonable.
#
#Using breakpoint analysis at population level, we found that across years and groups, a breakpoint in
#incubation pattern is reached at day 15. After this day, nests are incubated in average 69% of the time.
#Using this average incubation effort, we define the onset of consistent incubation as the first day of
#the earliest three consecutive days in which a group invests at least 69% of their time in incubation.
#
#Function to find the earliest occurrence where the incubation effort meets the threshold (>=0.69)
#and is immediately followed by two more consecutive y-values also meeting the threshold.:
find_inc.Threshold <- function(df, y_col, x_col, threshold = 0.69) {
  y <- df[[y_col]]
  x <- df[[x_col]]
  n <- length(y)
  for (i in 1:(n - 2)) {
    if (y[i] >= threshold && y[i + 1] >= threshold && y[i + 2] >= threshold) {
      return(x[i])#Return x (date/time) value for the first qualifying run
    }
  }
  return(NA)
}
#
#Find incubation thresholds
inc_thresholds <- sitdur_date_filtered[sitdur_date_filtered$dayordereggsleft >= 0,] %>%
  arrange(yearcamp, dayordereggsleft) %>%
  group_by(yearcamp) %>%
  group_modify(~ tibble(inc_threshold = find_inc.Threshold(.x, y_col = "propsit", x_col = "dayordereggsleft"))) %>%
  ungroup()
#
#Create one multi-panel plot for each year with incubation thresholds for each group. 
#create a list with years
years <- 2012:2018
#Create an empty list to store plot objects
threshold_plots <- list()
#Run loop to create plots
for (yr in years) {
  year_str <- as.character(yr)
  #
  # Filter data for the current year
  df_year <- sitdur_date_filtered %>%
    dplyr::filter(dayordereggsleft >= 0, grepl(year_str, yearcamp))
  thresholds_year <- inc_thresholds %>% 
    dplyr::filter(grepl(year_str, yearcamp))
  #
  #Remove NA inc_thresholds
  thresholds_year_nonNA <- thresholds_year %>%
    dplyr::filter(!is.na(inc_threshold))
  #
  # Build the ggplot
  p <- ggplot(df_year, aes(x = dayordereggsleft, y = propsit)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_line(alpha = 0.7) +
    geom_hline(yintercept = 0.69, linetype = "dashed") +
    facet_wrap(~ yearcamp) +
    labs(title = paste("Proportion of Time Incubating by Group -", year_str),
         x = "Day in incubation period",
         y = "Proportion of time incubating") +
    theme_minimal()
  #
  #Add thresholds if any exist for this year
  if (nrow(thresholds_year_nonNA) > 0) {
    p <- p +
      geom_vline(data = thresholds_year_nonNA,
                 aes(xintercept = inc_threshold),
                 color = "red", linetype = "dashed", linewidth = 1) +
      geom_text(data = thresholds_year_nonNA,
                aes(x = inc_threshold, y = 0.4, label = paste("Day", inc_threshold)),
                color = "red", size = 3,
                inherit.aes = FALSE)
  }
  #
  # Store the plot in the list
  threshold_plots[[year_str]] <- p
}
#Visually control plots
threshold_plots[["2012"]]
threshold_plots[["2013"]]
threshold_plots[["2014"]]
threshold_plots[["2015"]]
threshold_plots[["2016"]]
threshold_plots[["2017"]]
threshold_plots[["2018"]]
#Looks good!
mean(inc_thresholds$inc_threshold[which(!is.na(inc_thresholds$inc_threshold))])
pop_seg$psi[1, "Est."]
#Groups reach the incubation threshold in average at day 15.
#This is in accordance with the population-level data.
#
#
####Individual-level incubation----
#**********************************
#NOTE: here, incubation is at any of the camp's nests as we are interested in individual effort, not nest protection.
#
#Create a vector of unique sessions
UniqueSessions <- unique(sitduration2012.2018$session)
#
#Are all the yearcamps in the incubation df in the treatments df?
unique(sitduration2012.2018$yearcamp[!sitduration2012.2018$yearcamp %in% treatments2012.2018$yearcamp])
#factor(0) -> Yes! All of them are in there.
#
#Create somewhere to collect output data
sitdurationAllInds <- as.numeric()
#
#Run a for loop to give an incubation time for every session to all individuals in the camps
for(i in UniqueSessions){
  #
  #i = "2012-10-27_06:50" #(A normal observation). Run this line in case of errors.
  #
  #Finding incubation observations of time and year (i)
  sitduration_temp <- sitduration2012.2018[sitduration2012.2018$session %in% i,]
  #
  #find potential individuals in the camp based on year
  PotHit <- treatments2012.2018[treatments2012.2018$year %in% sitduration_temp$year,]
  #
  #find individuals which were in the camps at the present date
  date_temp <- unique(sitduration_temp$date)
  Hit <- PotHit[PotHit$seasonstart_date <= date_temp & PotHit$seasonend_date >= date_temp,]
  #
  #Add sitduration and time of observation to the relevant camp
  ToAdd <- merge(Hit, sitduration_temp[,c("date", "observer", "idcode", "session", "obsduration", "timesitCorr", "timefinishCorr", "timestart", "timeend",
                                          "sitduration", "nest", "nestrank", "eggs")],
                 by = "idcode", all = TRUE)
  #
  #Tranfer session specific data to individuals that did not sit.
  ToAdd$date <- date_temp
  ToAdd$observer <- unique(na.omit(ToAdd$observer))
  ToAdd$session <- unique(na.omit(ToAdd$session))
  ToAdd$obsduration <- unique(na.omit(ToAdd$obsduration))
  ToAdd$timestart <- unique(na.omit(ToAdd$timestart))
  ToAdd$timeend <- unique(na.omit(ToAdd$timeend))
  ToAdd$eggs <- unique(na.omit(ToAdd$eggs))
  #
  #Check if last row has year = Na, which indicates a mismatch/missing info
  if(is.na(ToAdd[nrow(ToAdd),"year"]))
  {
    print("ERROR MISSING DATA")
    print(i)
  }
  #
  #Add to output dataframe
  sitdurationAllInds <- rbind(sitdurationAllInds,ToAdd)
  #
}
#
#Any idcodes from incubation data not matched in the treatments data?
sitdurationAllInds$idcode[is.na(sitdurationAllInds$year)] #
#factor(0) = Good!
#
#Change NAs to 0s in the sitduration column
sitdurationAllInds[is.na(sitdurationAllInds$sitduration),]$sitduration <- 0
#
#Compute sitduration per date for every individual
sitdur_date_ind <- sitdurationAllInds %>% group_by(yearID, date) %>%
  dplyr::summarize(sitdur_date = sum(sitduration, na.rm = TRUE), .groups = "drop")
#
#Add observation duration per day
#Compute observation duration per day
obsdur_date <- sitduration2012.2018 %>%  distinct(session, .keep_all = TRUE) %>%
  group_by(date) %>% dplyr::mutate(obsdur_date = sum(obsduration)) %>%
  distinct(date, .keep_all = TRUE) %>% dplyr::select(date, obsdur_date) %>% ungroup()
#Match data
sitdur_date_ind$obsdur_date <- obsdur_date$obsdur_date[match(sitdur_date_ind$date, obsdur_date$date)]
#
#Compute proportion of time incubating in a given camp and date
sitdur_date_ind$propsit <- sitdur_date_ind$sitdur_date / sitdur_date_ind$obsdur_date
#
#Match day from experimental onset of natural incubation
sitdur_date_ind$dayorderEggs <- sitduration2012.2018$dayorderEggs[match(sitdur_date_ind$date, sitduration2012.2018$date)]
#
#
####Camp-level egg laying----
#****************************
#
#Create a list and a data frame with all dates within the breeding season of every year
#Make sure date columns are Date type
treatments2012.2018 <- treatments2012.2018 %>%
  dplyr::mutate(seasonstart_date = as.Date(seasonstart_date),
         seasonend_date = as.Date(seasonend_date))
#Define years
years <- 2012:2018
#Create a data frame with all days in the breeding season for each year
all_days_df <- map_df(years, function(y) {
  #Filter for the year
  df_year <- dplyr::filter(treatments2012.2018, year == y)
  #Determine min start date and max end date for that year
  start_date <- min(df_year$seasonstart_date)
  end_date <- max(df_year$seasonend_date)
  #Create a sequence of dates for the breeding season in this year
  tibble(date = seq(start_date, end_date, by = "day"), year = y)
})
#Add a column where date is a factor
all_days_df <- all_days_df %>%
  dplyr::mutate(date_factor = factor(date))
#
#Create a columns with laydate_Date in factor format in the egg df
eggs2012.2018 <- eggs2012.2018 %>%
  dplyr::mutate(laydate_date = as.Date(laydate_date),
         laydate_Date_Factor = factor(laydate_date))
#
#Create a data frame with all the camps for each year
allcamps2012.2018 <- treatments2012.2018 %>% distinct(yearcamp, .keep_all = T) %>% dplyr::select(yearcamp, year)
#
#Run a for loop to add zeros to camps that did not produce any eggs in a given date
#Create new data frame to collect output data
eggs2012.2018wZeros <- as.numeric()
#Run the loop
for(i in all_days_df$date_factor){
  #
  #i = "2012-05-15" #An observation without any eggs laid. Run this to see whether days without eggs are being added properly.
  #i = "2012-07-23" #A normal day with four eggs laid. Run this line in case of errors.
  #
  #Finding eggs laid in date (i)
  eggs_temp <- eggs2012.2018[eggs2012.2018$laydate_Date_Factor %in% i,]
  #
  #Find camps in "i" based on year
  Hit <- allcamps2012.2018[allcamps2012.2018$year %in% all_days_df$year[which(all_days_df$date_factor == i)],]
  #
  #Add relevant information from the egg df
  ToAdd <- merge(Hit, eggs_temp[,c("egglaid", "e_hatchoutcome", "c_hatched", "grpsize","nomales", "nofemales", 
                                   "laydate_date", "yearcamp", "eggid")], by = "yearcamp", all = TRUE)
  #
  #Transfer data on group and date to camps without any eggs in  a given date
  ToAdd$grpsize <- treatments2012.2018$grpsize[match(ToAdd$yearcamp, treatments2012.2018$yearcamp)]
  ToAdd$nomales <- treatments2012.2018$nomales[match(ToAdd$yearcamp, treatments2012.2018$yearcamp)]
  ToAdd$nofemales <- treatments2012.2018$nofemales[match(ToAdd$yearcamp, treatments2012.2018$yearcamp)]
  ToAdd$laydate_date <- as.Date(i)
  #Give camps without any eggs in a given date a value of 0 for egglaid
  ToAdd$egglaid[which(is.na(ToAdd$egglaid))] <- 0
  #
  #Add to output dataframe
  eggs2012.2018wZeros <- rbind(eggs2012.2018wZeros,ToAdd)
}
#Check the new data
table(eggs2012.2018wZeros$egglaid, exclude = NULL)#14875 new camp-dates with no eggs added.
eggs2012.2018[!eggs2012.2018$eggid %in% eggs2012.2018wZeros$eggid,]#There are some eggs missing,
#but all are eggs with no lay date
#
#Order data by laydate
eggs2012.2018wZeros <- eggs2012.2018wZeros[order(eggs2012.2018wZeros$laydate_date),]
#
#Create a df with dates for start of incubation period
StartInc <- treatments2012.2018 %>% distinct(year, .keep_all = T) %>% dplyr::select(year, eggsleftDate)
#
#Match data to egg df
eggs2012.2018wZeros$eggsleftDate <- StartInc$eggsleftDate[match(eggs2012.2018wZeros$year, StartInc$year)]
#
#Create a column with information about whether the eggs were collected from camps or not
eggs2012.2018wZeros <- eggs2012.2018wZeros %>% dplyr::mutate(eggsleft = ifelse(laydate_date < eggsleftDate, "No", "Yes"))
#
#Create a column with days from start of incubation period
eggs2012.2018wZeros$dayorderEggs <- as.numeric(eggs2012.2018wZeros$laydate_date - eggs2012.2018wZeros$eggsleftDate)
#
#Sub-set data to only include period with natural incubation, up to 32 days.
eggs <- eggs2012.2018wZeros[eggs2012.2018wZeros$eggsleft == "Yes" & eggs2012.2018wZeros$dayorderEggs < 33,]
#
#Compute total number of eggs laid per camp during period with natural incubation, up to 32 days
eggs_camp <- eggs %>% group_by(yearcamp) %>%
  dplyr::summarize(tot_no.eggs = sum(egglaid, na.rm = TRUE)) %>%
  as.data.frame()
#
#Compute eggs laid per date and camp
eggs_date_camp <- eggs %>% group_by(yearcamp, laydate_date) %>%
  dplyr::summarize(no.eggs = sum(egglaid, na.rm = TRUE), .groups = "drop")
#
#Compute cumulative number of eggs for each nest by date
eggs_date_camp <- eggs_date_camp %>%
  dplyr::group_by(yearcamp) %>%
  dplyr::arrange(laydate_date) %>%
  dplyr::mutate(cumulative_eggs = cumsum(no.eggs)) %>%
  ungroup() %>% as.data.frame()
#
#Match total number of eggs in camp
eggs_date_camp$tot_no.eggs <- eggs_camp$tot_no.eggs[match(eggs_date_camp$yearcamp, eggs_camp$yearcamp)]
#
#Are there any camps with no eggs?
unique(eggs_date_camp$yearcamp[eggs_date_camp$tot_no.eggs == 0])#Yes, "2012_68" & "2016_66"
#Get rid!
eggs_date_camp <- eggs_date_camp[eggs_date_camp$tot_no.eggs != 0,]
#Check if these camps are in the incubation dataset
sitdur_date_filtered[sitdur_date_filtered$yearcamp == "2012_68" | sitdur_date_filtered$yearcamp == "2016_66",]#No, good!
#
#Compute proportion of eggs in camp out of the total number at the end of incubation
eggs_date_camp$propeggs <- eggs_date_camp$cumulative_eggs / eggs_date_camp$tot_no.eggs
#
#Match day from experimental onset of natural incubation
eggs_date_camp$dayorderEggs <- eggs$dayorderEggs[match(eggs_date_camp$laydate_date, eggs$laydate_date)]
#
#
####Individual-level egg laying----
#**********************************
#
#Create a data frame with all dates within the period with care and 28 days weeks before
#Define years
years <- 2012:2018
# Make sure dates are Date class
treatments2012.2018 <- treatments2012.2018 %>%
  dplyr::mutate(eggsleftDate = as.Date(eggsleftDate, "%Y-%m-%d"),
         seasonend_date = as.Date(seasonend_date, "%Y-%m-%d"))
#Generate sequences of dates with all dates within the period with care and 28 days weeks before
incDays_list <- lapply(years, function(yr) {
  start_date <- min(treatments2012.2018$eggsleftDate[treatments2012.2018$year == yr]) - 28
  end_date <- max(treatments2012.2018$seasonend_date[treatments2012.2018$year == yr])
  data.frame(date = seq(start_date, end_date, by = "day"), year = yr)
})
# Combine all into one data frame
incDays2012.2018 <- bind_rows(incDays_list)
#Format columns
incDays2012.2018 <- incDays2012.2018 %>%
  dplyr::mutate(date_factor = factor(date), year = as.factor(year))
#
#Create a data frame with all the individuals present during incubation
allinds <- treatments2012.2018 %>%
  dplyr::filter(incampeggsleft == "Yes") %>%
  dplyr::select(yearID, year)
#
#Run a for loop to add zeros to individuals that did not lay any eggs in a given date
#Create new data frame to collect output data
indLayingwZeros <- as.numeric()
#Run the loop
for(i in incDays2012.2018$date_factor){
  #
  #i = "2012-10-13"#Example day
  #
  #Finding eggs laid in date (i) in the parentage
  parentage_temp <- parentage2012.2018[parentage2012.2018$laydate_date %in% i,]
  #
  #Find individuals in "i" based on year
  Hit <- allinds[allinds$year %in% incDays2012.2018$year[which(incDays2012.2018$date_factor == i)],]
  #
  #Add relevant information from the parentage df
  ToAdd <- merge(Hit, parentage_temp[,c("eggid", "c_hatched", "laydate_date", "eggsleftDate", "daysfromEggsleft",
                                        "yearID_Mother", "yearcamp")], by.x = "yearID", by.y = "yearID_Mother",
                 all = TRUE)
  #
  #Transfer data on date and camp to individuals without any eggs in  a given date
  ToAdd$eggsleftDate <- treatments2012.2018$eggsleftDate[match(ToAdd$year, treatments2012.2018$year)]
  ToAdd$daysfromEggsleft <- as.numeric(as.Date(i) - ToAdd$eggsleftDate)
  ToAdd$yearcamp <- treatments2012.2018$yearcamp[match(ToAdd$yearID, treatments2012.2018$yearID)]
  ToAdd$laydate_date <- i
  #Create a column for probability of laying an egg 
  ToAdd$probEgg <- ifelse(is.na(ToAdd$eggid), 0, 1)
  #
  #Add to output dataframe
  indLayingwZeros <- rbind(indLayingwZeros,ToAdd)
}
#
#Add breeding strategy and order by laydate
indLayingwZeros <- indLayingwZeros %>%
  #Order by laydate_date
  arrange(laydate_date) %>%
  #Get BS fron individuals data
  left_join(dplyr::select(individuals2012.2018, yearID, BS), by = "yearID") %>%
  dplyr::filter(!is.na(BS))#Remove individuals with no BS
#
#Filter out days after day 32 of natural incubation
indLayingwZeros <- indLayingwZeros[indLayingwZeros$daysfromEggsleft < 33,]
#
#
#******#
### Export data----
#******#
#
save(individuals2012.2018, camps2012.2018, Femfem_coop_agg_count_care4_grp_long, treatments2012.2018, mass,
     indLayingwZeros, eggs_date_camp, sitdur_date_filtered, sitdur_date_ind, incNest_F,
     file = "R_Output_files/Cheating_data.RData")
#
#
#