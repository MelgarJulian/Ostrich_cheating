##*********************************
#Frequency-dependent cheating in a cooperative breeding bird - MCMCglmm models----
##*********************************
#Script by Julian Melgar, Charlie Cornwallis & Mads Schou
#
##*****************
##Preliminaries----
##*****************
#
###Load packages----
#*******************
pacman::p_load(MCMCglmm, ggplot2, dplyr, doParallel, foreach, segmented, nlme)
#
###Load data----
#***************
load("R_Output_files/CheatingModels.RData")
#
##*******************************
##Guide to model object names----
##*******************************
#
#Model objects are named for what they estimate, not for their position in the
#manuscript, so that models can be reordered, added or dropped without renaming
#anything. Each name is built as:
#
#  m_<topic>_<structure>_<component>
#
#  <topic>      what the model is about, e.g. m_eggs, m_nest_failure
#  <structure>  which set of terms was fitted:
#                 main  main effects only
#                 int   main effects plus two-way interactions
#                 int3  main effects plus two- and three-way interactions
#                       (only m_incubation_individual goes this far)
#                 quad  as int, plus interactions with quadratic terms
#  <component>  c1, c2, c3  the three MCMC chains of that fit
#               sol, vcv    the mcmc.list combining the three chains, used for
#                           the convergence checks
#
#So m_eggs_quad_c2 is the second chain of the egg-number model fitted with
#quadratic interactions, and m_eggs_quad_sol holds all three chains of its Sol.
#Models reported at only one level of complexity have no _int or _quad objects.
#
#The manuscript and Extended Data tables number models in reading order
#(manuscript Model N -> Extended Data Table N). That numbering changes whenever the paper is
#restructured, so it is recorded here rather than in the object names:
#
#  Topic                      Manuscript   Table   Subject
#  m_strategy_counts          Model 1      S1      Group size and the occurrence of cheats and cooperators
#  m_cheat_prop               Model 2      S2      Group size and the proportion of cheats
#  m_cheat_prob               Model 3      S3      Predictors of the probability of being a cheat
#  m_mass                     Model 3      S3b     Body mass at the end of the season by breeding strategy
#  m_eggs                     Model 4      S4      Cheating, group composition and the number of eggs laid
#  m_eggs_excl_1coop          Model 4b     S4b     As Model 4, excluding groups with one cooperator
#  m_eggs_within_female       Model 5      S5      Reproductive success when the same female cheats vs cooperates
#  m_laying_group             Model 6      S6      Group-level egg laying across time
#  m_laying_window            Model 6      S6b     Probability of laying within a two-day window
#  m_aggression               Model 7      S7      Cooperator aggression at group level
#  m_incubation_group         Model 8      S8      Cheating, group composition and group-level incubation effort
#  m_incubation_time          Model 9      S9      Cheating and group-level incubation effort across time
#  m_incubation_individual    Model 10     S10     Group composition, reproductive success and individual incubation
#  m_nest_failure             Model 11     S11     Cheating and complete nest failure
#  m_hatch_success            Model 11     S11b    Hatching success of nests that hatched at least one egg
#  m_activity_nocheat,        Model 12     S12     Cheating and the activity allocation of groups
#    m_activity_cheat
#
#m_mass and m_laying_window are described at the end of the Model 3 and Model 6
#sections of the manuscript methods rather than under headings of their own.
#Two objects fall outside the scheme because they are one-off checks rather than
#part of a series: m_nest_failure_broken (nest failure refitted with the
#proportion of eggs broken) and fit_activity_nocheat / fit_activity_cheat (the
#lme fits that the segmented activity models are built from).
#
#The data frames follow the same topic names with a df_ prefix, so m_eggs is
#fitted to df_eggs. They are built in Cheating - Model data.R, which carries the
#full list at the top.
###Model settings----
#********************
#
#Every model is fitted at these settings, giving 3600 stored samples per chain
#and 10800 across the three chains. They are the settings reported in the
#manuscript methods, and they are not a default to be adjusted casually:
#changing them means the committed chains in CheatingModels.RData no longer
#correspond to the script, and the reported settings become wrong.
nitts<-2000000
burn<-200000
thin<-500
#
#To try a change quickly, set these by hand in the console rather than editing
#the lines above, so that what is committed always states what was actually run:
#  nitts<-1100; burn<-100; thin<-1     #test run, results are meaningless
#
##**********************
##Random number seeds----
##**********************
#
#Each fit is preceded by its own set.seed(), so any model can be re-run on its
#own and reproduce exactly. The seed is written at the call site, so it travels
#with the model if the script is reordered.
#
#The value encodes the model, so duplicates are easy to spot:
#  seed = topic * 100 + structure * 10 + chain
#with topics numbered in Extended Data table order (m_strategy_counts 1 through
#m_hatch_success 15), structure main 1, int 2, quad or int3 3, and chain 1 to 3.
#So m_eggs_quad_c2 is 532. The one-off m_nest_failure_broken uses 1409.
#
#The model image is saved with compress = "xz" throughout. Posterior samples are
#dense doubles, which the default gzip barely compresses: the same image is
#122 MB under gzip and 26 MB under xz. That keeps it under the 100 MB per-file
#limit GitHub enforces, at the cost of a few seconds per save.
#
#Set plot margins to 1, to avoid plotting errors
par(mar=c(1,1,1,1))
#
#
start <- Sys.time()
#
##*********##
##Priors----
##*********##
#
#The priors are defined together because most are shared between models:
#prior1 by m_strategy_counts, m_laying_group, m_aggression and m_incubation_time,
#prior3 by m_cheat_prob and m_laying_window, prior4 by the five individual-level
#models, and prior5 by m_incubation_group and m_hatch_success. Defining them in
#one place means any model block can be run on its own.
#
prior1 = list(R = list(V=1, nu=0.002), G = list(G1=list(V=1, nu=0.002), G2=list(V=1, nu=0.002), G3=list(V=1, nu=0.002)))
#
prior2 = list(R = list(V=1, nu=0.002), G = list(G1=list(V=1, nu=0.002), G2=list(V=1, nu=0.002)))
#
prior3 = list(R = list(V=1, nu=0.002,fix=1), G = list(G1=list(V=1, nu=0.002), G2=list(V=1, nu=0.002), G3=list(V=1, nu=0.002), G4=list(V=1, nu=0.002)))
#
prior4 = list(R = list(V=1, nu=0.002), G = list(G1=list(V=1, nu=0.002), G2=list(V=1, nu=0.002), G3=list(V=1, nu=0.002), G4=list(V=1, nu=0.002)))
#
prior5 = list(R = list(V=1, nu=0.002), G = list(G1=list(V=1, nu=0.002), G2=list(V=1, nu=0.002)))
#
#Residual variance must be fixed for a binary response as it is not identifiable
prior8 = list(R = list(V=1, nu=0.002, fix=1), G = list(G1=list(V=1, nu=0.002), G2=list(V=1, nu=0.002)))
#

#***************************************
##Frequency of cooperators and cheats----
#***************************************
#
###m_strategy_counts: Effect of group size on number of cheats and cooperators (Model 1, Extended Data Table 1)----
##***********************************************************************
#
#
#
#Run main effects
set.seed(111)
m_strategy_counts_main_c1 <- MCMCglmm(n  ~  BS-1 + Znonests,
  random= ~ year + camp + yearcamp,
  family ="poisson", data = df_strategy_counts,
  prior = prior1, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
summary(m_strategy_counts_main_c1)
set.seed(112)
m_strategy_counts_main_c2 <- MCMCglmm(n  ~  BS-1 + Znonests,
  random= ~ year + camp + yearcamp,
  family ="poisson", data = df_strategy_counts,
  prior = prior1, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
set.seed(113)
m_strategy_counts_main_c3 <- MCMCglmm(n  ~  BS-1 + Znonests,
  random= ~ year + camp + yearcamp,
  family ="poisson", data = df_strategy_counts,
  prior = prior1, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
#
#Check convergence of models#
#
par(mar=c(2,2,2,2))#To avoid errors with figure margins being too large
m_strategy_counts_main_sol<-mcmc.list(m_strategy_counts_main_c1$Sol,m_strategy_counts_main_c2$Sol,m_strategy_counts_main_c3$Sol);plot(m_strategy_counts_main_sol)
m_strategy_counts_main_vcv<-mcmc.list(m_strategy_counts_main_c1$VCV,m_strategy_counts_main_c2$VCV,m_strategy_counts_main_c3$VCV);plot(m_strategy_counts_main_vcv)
gelman.diag(m_strategy_counts_main_sol)
gelman.diag(m_strategy_counts_main_vcv, multivariate = FALSE)
#
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#
#Then run two-way interactions
set.seed(121)
m_strategy_counts_int_c1 <- MCMCglmm(n  ~  BS-1 + Znonests +
    BS:Znofemales,
  random= ~ year + camp + yearcamp,
  family ="poisson", data = df_strategy_counts,
  prior = prior1, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
summary(m_strategy_counts_int_c1)
set.seed(122)
m_strategy_counts_int_c2 <- MCMCglmm(n  ~  BS-1 + Znonests +
    BS:Znofemales,
  random= ~ year + camp + yearcamp,
  family ="poisson", data = df_strategy_counts,
  prior = prior1, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
set.seed(123)
m_strategy_counts_int_c3 <- MCMCglmm(n  ~  BS-1 + Znonests +
    BS:Znofemales,
  random= ~ year + camp + yearcamp,
  family ="poisson", data = df_strategy_counts,
  prior = prior1, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
#
#Check convergence of models#
#
par(mar=c(2,2,2,2))#To avoid errors with figure margins being too large
m_strategy_counts_int_sol<-mcmc.list(m_strategy_counts_int_c1$Sol,m_strategy_counts_int_c2$Sol,m_strategy_counts_int_c3$Sol);plot(m_strategy_counts_int_sol)
m_strategy_counts_int_vcv<-mcmc.list(m_strategy_counts_int_c1$VCV,m_strategy_counts_int_c2$VCV,m_strategy_counts_int_c3$VCV);plot(m_strategy_counts_int_vcv)
gelman.diag(m_strategy_counts_int_sol)
gelman.diag(m_strategy_counts_int_vcv, multivariate = FALSE)
#
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
###m_cheat_prop: Effect of number of cooperators on proportion of cheats (Model 2, Extended Data Table 2)----
##**********************************************************************
#
#
#Run model
set.seed(211)
m_cheat_prop_main_c1 <- MCMCglmm(cbind(noFCheatsCamp, noFActive - noFCheatsCamp)  ~  ZnoFCoopsCamp + ZnoFOtherCamp + Znonests,
  random= ~ year + camp,
  family ="multinomial2", data = df_cheat_prop,
  prior = prior2, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
summary(m_cheat_prop_main_c1)
set.seed(212)
m_cheat_prop_main_c2 <- MCMCglmm(cbind(noFCheatsCamp, noFActive - noFCheatsCamp)  ~  ZnoFCoopsCamp + ZnoFOtherCamp + Znonests,
  random= ~ year + camp,
  family ="multinomial2", data = df_cheat_prop,
  prior = prior2, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
set.seed(213)
m_cheat_prop_main_c3 <- MCMCglmm(cbind(noFCheatsCamp, noFActive - noFCheatsCamp)  ~  ZnoFCoopsCamp + ZnoFOtherCamp + Znonests,
  random= ~ year + camp,
  family ="multinomial2", data = df_cheat_prop,
  prior = prior2, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
#
#Check convergence of models#
#
par(mar=c(2,2,2,2))#To avoid errors with figure margins being too large
m_cheat_prop_main_sol<-mcmc.list(m_cheat_prop_main_c1$Sol,m_cheat_prop_main_c2$Sol,m_cheat_prop_main_c3$Sol);plot(m_cheat_prop_main_sol)#**Fit does not look very good (nitts=1100000)*
m_cheat_prop_main_vcv<-mcmc.list(m_cheat_prop_main_c1$VCV,m_cheat_prop_main_c2$VCV,m_cheat_prop_main_c3$VCV);plot(m_cheat_prop_main_vcv)
gelman.diag(m_cheat_prop_main_sol)
gelman.diag(m_cheat_prop_main_vcv, multivariate = FALSE)
#
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
###m_cheat_prob: How does probability of being a cheat change with group structure and age? (Model 3, Extended Data Table 3)----
#******************************************************************************************
#
#Run model
set.seed(311)
m_cheat_prob_main_c1<- MCMCglmm(prob.Cheat ~ Zcooperators_focal + Zother_focal + Zage + Zmassin,
  random=~year + camp + yearcamp + id,
  family ="threshold", data = df_cheat_prob,
  prior = prior3, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE,trunc=T)
summary(m_cheat_prob_main_c1)
set.seed(312)
m_cheat_prob_main_c2<- MCMCglmm(prob.Cheat ~ Zcooperators_focal + Zother_focal + Zage + Zmassin,
  random=~year + camp + yearcamp + id,
  family ="threshold", data = df_cheat_prob,
  prior = prior3, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE,trunc=T)
set.seed(313)
m_cheat_prob_main_c3<- MCMCglmm(prob.Cheat ~ Zcooperators_focal + Zother_focal + Zage + Zmassin,
  random=~year + camp + yearcamp + id,
  family ="threshold", data = df_cheat_prob,
  prior = prior3, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE,trunc=T)
#
#Check convergence of models
par(mar=c(2,2,2,2))#To avoid errors with figure margins being too large
m_cheat_prob_main_sol<-mcmc.list(m_cheat_prob_main_c1$Sol,m_cheat_prob_main_c2$Sol,m_cheat_prob_main_c3$Sol);plot(m_cheat_prob_main_sol)
m_cheat_prob_main_vcv<-mcmc.list(m_cheat_prob_main_c1$VCV,m_cheat_prob_main_c2$VCV,m_cheat_prob_main_c3$VCV);plot(m_cheat_prob_main_vcv)
gelman.diag(m_cheat_prob_main_sol)
gelman.diag(m_cheat_prob_main_vcv, multivariate = FALSE)
#
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#***********************************************
##Consequences of cheating for the individual----
#***********************************************
#
###m_mass: How is mass affected by breeding strategy, age and social context (Model 4, Extended Data Table 4)----
#*****************************************************************************
#
#Run model with main effects
set.seed(411)
m_mass_main_c1 <- MCMCglmm(massout ~ BS-1,
                           random= ~ camp + year + yearcamp + id,
                           family ="gaussian", data = df_mass,
                           prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
summary(m_mass_main_c1)
set.seed(412)
m_mass_main_c2 <- MCMCglmm(massout ~ BS-1,
                           random= ~ camp + year + yearcamp + id,
                           family ="gaussian", data = df_mass,
                           prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
set.seed(413)
m_mass_main_c3 <- MCMCglmm(massout ~ BS-1,
                           random= ~ camp + year + yearcamp + id,
                           family ="gaussian", data = df_mass,
                           prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
#
#Check convergence of models#
#
par(mar=c(2,2,2,2))#To avoid errors with figure margins being too large
m_mass_main_sol<-mcmc.list(m_mass_main_c1$Sol,m_mass_main_c2$Sol,m_mass_main_c3$Sol);plot(m_mass_main_sol)
m_mass_main_vcv<-mcmc.list(m_mass_main_c1$VCV,m_mass_main_c2$VCV,m_mass_main_c3$VCV);plot(m_mass_main_vcv)
gelman.diag(m_mass_main_sol)
gelman.diag(m_mass_main_vcv, multivariate = FALSE)
#
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#
#Then run two-way interactions
set.seed(421)
m_mass_int_c1 <- MCMCglmm(massout ~ BS-1 + BS:ZnoFCoops + BS:ZnoFOther + BS:Zmassin,
                          random= ~ camp + year + yearcamp + id,
                          family ="gaussian", data = df_mass,
                          prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
summary(m_mass_int_c1)
set.seed(422)
m_mass_int_c2 <- MCMCglmm(massout ~ BS-1 + BS:ZnoFCoops + BS:ZnoFOther + BS:Zmassin,
                          random= ~ camp + year + yearcamp + id,
                          family ="gaussian", data = df_mass,
                          prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
set.seed(423)
m_mass_int_c3 <- MCMCglmm(massout ~ BS-1 + BS:ZnoFCoops + BS:ZnoFOther + BS:Zmassin,
                          random= ~ camp + year + yearcamp + id,
                          family ="gaussian", data = df_mass,
                          prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
#
#Check convergence of models#
#
par(mar=c(2,2,2,2))#To avoid errors with figure margins being too large
m_mass_int_sol<-mcmc.list(m_mass_int_c1$Sol,m_mass_int_c2$Sol,m_mass_int_c3$Sol);plot(m_mass_int_sol)
m_mass_int_vcv<-mcmc.list(m_mass_int_c1$VCV,m_mass_int_c2$VCV,m_mass_int_c3$VCV);plot(m_mass_int_vcv)
gelman.diag(m_mass_int_sol)
gelman.diag(m_mass_int_vcv, multivariate = FALSE)
#
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
###m_eggs: The effect of the presence of cheats on reproductive investment (Model 5, Extended Data Table 5)----
##******************************************************************************
#
#Start with the main effects
set.seed(511)
m_eggs_main_c1 <- MCMCglmm(RSInc4 ~ BS_2-1 + ZnoFOther + logit_prop.BrokenInc,
  random= ~ camp + year + yearcamp + id,
  family ="poisson", data = df_eggs,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
summary(m_eggs_main_c1)
set.seed(512)
m_eggs_main_c2 <- MCMCglmm(RSInc4 ~ BS_2-1 + ZnoFOther + logit_prop.BrokenInc,
  random= ~ camp + year + yearcamp + id,
  family ="poisson", data = df_eggs,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
set.seed(513)
m_eggs_main_c3 <- MCMCglmm(RSInc4 ~ BS_2-1 + ZnoFOther + logit_prop.BrokenInc,
  random= ~ camp + year + yearcamp + id,
  family ="poisson", data = df_eggs,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
#
#Check convergence of models
par(mar=c(2,2,2,2))#To avoid errors with figure margins being too large
m_eggs_main_sol<-mcmc.list(m_eggs_main_c1$Sol,m_eggs_main_c2$Sol,m_eggs_main_c3$Sol);plot(m_eggs_main_sol)
m_eggs_main_vcv<-mcmc.list(m_eggs_main_c1$VCV,m_eggs_main_c2$VCV,m_eggs_main_c3$VCV);plot(m_eggs_main_vcv)
gelman.diag(m_eggs_main_sol)
gelman.diag(m_eggs_main_vcv, multivariate = FALSE)
#
#Then include two-way interactions
set.seed(521)
m_eggs_int_c1 <- MCMCglmm(RSInc4 ~ BS_2-1 + ZnoFOther + logit_prop.BrokenInc +
    BS_2:ZnoFCoops,
  random= ~ camp + year + yearcamp + id,
  family ="poisson", data = df_eggs,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
summary(m_eggs_int_c1)
set.seed(522)
m_eggs_int_c2 <- MCMCglmm(RSInc4 ~ BS_2-1 + ZnoFOther + logit_prop.BrokenInc +
    BS_2:ZnoFCoops,
  random= ~ camp + year + yearcamp + id,
  family ="poisson", data = df_eggs,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
set.seed(523)
m_eggs_int_c3 <- MCMCglmm(RSInc4 ~ BS_2-1 + ZnoFOther + logit_prop.BrokenInc +
    BS_2:ZnoFCoops,
  random= ~ camp + year + yearcamp + id,
  family ="poisson", data = df_eggs,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
#
#Check convergence of models
par(mar=c(2,2,2,2))#To avoid errors with figure margins being too large
m_eggs_int_sol<-mcmc.list(m_eggs_int_c1$Sol,m_eggs_int_c2$Sol,m_eggs_int_c3$Sol);plot(m_eggs_int_sol)
m_eggs_int_vcv<-mcmc.list(m_eggs_int_c1$VCV,m_eggs_int_c2$VCV,m_eggs_int_c3$VCV);plot(m_eggs_int_vcv)
gelman.diag(m_eggs_int_sol)
gelman.diag(m_eggs_int_vcv, multivariate = FALSE)
#
#Finally include interactions with quadratics
set.seed(531)
m_eggs_quad_c1 <- MCMCglmm(RSInc4 ~ BS_2-1 + ZnoFOther + logit_prop.BrokenInc +
    BS_2:ZnoFCoops +
    BS_2:ZnoFCoops2,
  random= ~ camp + year + yearcamp + id,
  family ="poisson", data = df_eggs,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
summary(m_eggs_quad_c1)
set.seed(532)
m_eggs_quad_c2 <- MCMCglmm(RSInc4 ~ BS_2-1 + ZnoFOther + logit_prop.BrokenInc +
    BS_2:ZnoFCoops +
    BS_2:ZnoFCoops2,
  random= ~ camp + year + yearcamp + id,
  family ="poisson", data = df_eggs,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
set.seed(533)
m_eggs_quad_c3 <- MCMCglmm(RSInc4 ~ BS_2-1 + ZnoFOther + logit_prop.BrokenInc +
    BS_2:ZnoFCoops +
    BS_2:ZnoFCoops2,
  random= ~ camp + year + yearcamp + id,
  family ="poisson", data = df_eggs,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)

#
#Check convergence of models
m_eggs_quad_sol<-mcmc.list(m_eggs_quad_c1$Sol,m_eggs_quad_c2$Sol,m_eggs_quad_c3$Sol);plot(m_eggs_main_sol)
m_eggs_quad_vcv<-mcmc.list(m_eggs_quad_c1$VCV,m_eggs_quad_c2$VCV,m_eggs_quad_c3$VCV);plot(m_eggs_quad_vcv)
gelman.diag(m_eggs_quad_sol)
gelman.diag(m_eggs_quad_vcv, multivariate = FALSE)
#
#
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
###m_eggs_excl_1coop: as m_eggs, but excluding groups with one cooperator (Model 5b, Extended Data Table 5)----
##********************************************************************
#
#Start with the main effects
set.seed(611)
m_eggs_excl_1coop_main_c1 <- MCMCglmm(RSInc4 ~ BS_2-1 + ZnoFOther + logit_prop.BrokenInc,
  random= ~ camp + year + yearcamp + id,
  family ="poisson", data = df_eggs_excl_1coop,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
summary(m_eggs_excl_1coop_main_c1)
set.seed(612)
m_eggs_excl_1coop_main_c2 <- MCMCglmm(RSInc4 ~ BS_2-1 + ZnoFOther + logit_prop.BrokenInc,
  random= ~ camp + year + yearcamp + id,
  family ="poisson", data = df_eggs_excl_1coop,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
set.seed(613)
m_eggs_excl_1coop_main_c3 <- MCMCglmm(RSInc4 ~ BS_2-1 + ZnoFOther + logit_prop.BrokenInc,
  random= ~ camp + year + yearcamp + id,
  family ="poisson", data = df_eggs_excl_1coop,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
#
#Check convergence of models
par(mar=c(2,2,2,2))#To avoid errors with figure margins being too large
m_eggs_excl_1coop_main_sol<-mcmc.list(m_eggs_excl_1coop_main_c1$Sol,m_eggs_excl_1coop_main_c2$Sol,m_eggs_excl_1coop_main_c3$Sol);plot(m_eggs_excl_1coop_main_sol)
m_eggs_excl_1coop_main_vcv<-mcmc.list(m_eggs_excl_1coop_main_c1$VCV,m_eggs_excl_1coop_main_c2$VCV,m_eggs_excl_1coop_main_c3$VCV);plot(m_eggs_excl_1coop_main_vcv)
gelman.diag(m_eggs_excl_1coop_main_sol)
gelman.diag(m_eggs_excl_1coop_main_vcv, multivariate = FALSE)
#
#Then include two-way interactions
set.seed(621)
m_eggs_excl_1coop_int_c1 <- MCMCglmm(RSInc4 ~ BS_2-1 + ZnoFOther + logit_prop.BrokenInc +
    BS_2:ZnoFCoops,
  random= ~ camp + year + yearcamp + id,
  family ="poisson", data = df_eggs_excl_1coop,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
summary(m_eggs_excl_1coop_int_c1)
set.seed(622)
m_eggs_excl_1coop_int_c2 <- MCMCglmm(RSInc4 ~ BS_2-1 + ZnoFOther + logit_prop.BrokenInc +
    BS_2:ZnoFCoops,
  random= ~ camp + year + yearcamp + id,
  family ="poisson", data = df_eggs_excl_1coop,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
set.seed(623)
m_eggs_excl_1coop_int_c3 <- MCMCglmm(RSInc4 ~ BS_2-1 + ZnoFOther + logit_prop.BrokenInc +
    BS_2:ZnoFCoops,
  random= ~ camp + year + yearcamp + id,
  family ="poisson", data = df_eggs_excl_1coop,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
#
#Check convergence of models
par(mar=c(2,2,2,2))#To avoid errors with figure margins being too large
m_eggs_excl_1coop_int_sol<-mcmc.list(m_eggs_excl_1coop_int_c1$Sol,m_eggs_excl_1coop_int_c2$Sol,m_eggs_excl_1coop_int_c3$Sol);plot(m_eggs_excl_1coop_int_sol)
m_eggs_excl_1coop_int_vcv<-mcmc.list(m_eggs_excl_1coop_int_c1$VCV,m_eggs_excl_1coop_int_c2$VCV,m_eggs_excl_1coop_int_c3$VCV);plot(m_eggs_excl_1coop_int_vcv)
gelman.diag(m_eggs_excl_1coop_int_sol)
gelman.diag(m_eggs_excl_1coop_int_vcv, multivariate = FALSE)
#
#Finally include interactions with quadratics
set.seed(631)
m_eggs_excl_1coop_quad_c1 <- MCMCglmm(RSInc4 ~ BS_2-1 + ZnoFOther + logit_prop.BrokenInc +
    BS_2:ZnoFCoops +
    BS_2:ZnoFCoops2,
  random= ~ camp + year + yearcamp + id,
  family ="poisson", data = df_eggs_excl_1coop,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
summary(m_eggs_excl_1coop_quad_c1)
set.seed(632)
m_eggs_excl_1coop_quad_c2 <- MCMCglmm(RSInc4 ~ BS_2-1 + ZnoFOther + logit_prop.BrokenInc +
    BS_2:ZnoFCoops +
    BS_2:ZnoFCoops2,
  random= ~ camp + year + yearcamp + id,
  family ="poisson", data = df_eggs_excl_1coop,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
set.seed(633)
m_eggs_excl_1coop_quad_c3 <- MCMCglmm(RSInc4 ~ BS_2-1 + ZnoFOther + logit_prop.BrokenInc +
    BS_2:ZnoFCoops +
    BS_2:ZnoFCoops2,
  random= ~ camp + year + yearcamp + id,
  family ="poisson", data = df_eggs_excl_1coop,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
#
#Check convergence of models
m_eggs_excl_1coop_quad_sol<-mcmc.list(m_eggs_excl_1coop_quad_c1$Sol,m_eggs_excl_1coop_quad_c2$Sol,m_eggs_excl_1coop_quad_c3$Sol);plot(m_eggs_excl_1coop_quad_sol)
m_eggs_excl_1coop_quad_vcv<-mcmc.list(m_eggs_excl_1coop_quad_c1$VCV,m_eggs_excl_1coop_quad_c2$VCV,m_eggs_excl_1coop_quad_c3$VCV);plot(m_eggs_excl_1coop_quad_vcv)
gelman.diag(m_eggs_excl_1coop_quad_sol)
gelman.diag(m_eggs_excl_1coop_quad_vcv, multivariate = FALSE)
#
#
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
###m_eggs_within_female: Do opportunistic individuals have higher reproductive success when they cheat? (Model 6, Extended Data Table 6)----
##*****************************************************************************************
#
#Run model
set.seed(711)
m_eggs_within_female_main_c1 <- MCMCglmm(RSInc4  ~ BS-1 + Znofemales + Zage,
  random= ~ camp + year + yearcamp + id,
  family ="poisson", data = df_eggs_within_female,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
summary(m_eggs_within_female_main_c1)
set.seed(712)
m_eggs_within_female_main_c2 <- MCMCglmm(RSInc4  ~ BS-1 + Znofemales + Zage,
  random= ~ camp + year + yearcamp + id,
  family ="poisson", data = df_eggs_within_female,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
set.seed(713)
m_eggs_within_female_main_c3 <- MCMCglmm(RSInc4  ~ BS-1 + Znofemales + Zage,
  random= ~ camp + year + yearcamp + id,
  family ="poisson", data = df_eggs_within_female,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
summary(m_eggs_within_female_main_c1)
#
#This comparison tells us about the intercepts (mean RS of cheats vs cooperators)
HPDinterval(m_eggs_within_female_main_c1$Sol[,"BSCheat"]-m_eggs_within_female_main_c1$Sol[,"BSCooperator"])
#
#
#Check convergence of models#
#
par(mar=c(2,2,2,2))#To avoid errors with figure margins being too large
m_eggs_within_female_main_sol<-mcmc.list(m_eggs_within_female_main_c1$Sol,m_eggs_within_female_main_c2$Sol,m_eggs_within_female_main_c3$Sol);plot(m_eggs_within_female_main_sol)
m_eggs_within_female_main_vcv<-mcmc.list(m_eggs_within_female_main_c1$VCV,m_eggs_within_female_main_c2$VCV,m_eggs_within_female_main_c3$VCV);plot(m_eggs_within_female_main_vcv)
gelman.diag(m_eggs_within_female_main_sol)
gelman.diag(m_eggs_within_female_main_vcv, multivariate = FALSE)
#
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#************************
##Timing of egg laying----
#************************
#
###m_laying_group: How does nest-buildup differ between groups with and without cheats? (Model 7, Extended Data Table 7)----
#************************************************************************************
#
set.seed(811)
m_laying_group_main_c1 <- MCMCglmm(cbind(no.eggs, remaining - no.eggs) ~ CheatinCamp-1 + ZnoCoops + ZnoOther + ZDay,
  random = ~  year + camp + yearcamp,
  family = "multinomial2",
  data = df_laying_group,
  prior = prior1,
  nitt = nitts,
  burnin = burn,
  thin = thin,
  verbose = F)
summary(m_laying_group_main_c1)
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#
set.seed(812)
m_laying_group_main_c2 <- MCMCglmm(cbind(no.eggs, remaining - no.eggs) ~ CheatinCamp-1 + ZnoCoops + ZnoOther + ZDay,
  random = ~  year + camp + yearcamp,
  family = "multinomial2",
  data = df_laying_group,
  prior = prior1,
  nitt = nitts,
  burnin = burn,
  thin = thin,
  verbose = F)
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#
set.seed(813)
m_laying_group_main_c3 <- MCMCglmm(cbind(no.eggs, remaining - no.eggs) ~ CheatinCamp-1 + ZnoCoops + ZnoOther + ZDay,
  random = ~  year + camp + yearcamp,
  family = "multinomial2",
  data = df_laying_group,
  prior = prior1,
  nitt = nitts,
  burnin = burn,
  thin = thin,
  verbose = F)
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#
#Check convergence of models#
par(mar=c(2,2,2,2))#To avoid errors with figure margins being too large
m_laying_group_main_sol<-mcmc.list(m_laying_group_main_c1$Sol,m_laying_group_main_c2$Sol,m_laying_group_main_c3$Sol);plot(m_laying_group_main_sol)
m_laying_group_main_vcv<-mcmc.list(m_laying_group_main_c1$VCV,m_laying_group_main_c2$VCV,m_laying_group_main_c3$VCV);plot(m_laying_group_main_vcv)
gelman.diag(m_laying_group_main_sol)
gelman.diag(m_laying_group_main_vcv, multivariate = FALSE)
#
#Two-way interactions
set.seed(821)
m_laying_group_int_c1 <- MCMCglmm(cbind(no.eggs, remaining - no.eggs) ~ CheatinCamp-1 + ZnoCoops + ZnoOther + ZDay +
    CheatinCamp:ZDay-1 + ZnoCoops:ZDay,
  random = ~  year + camp + yearcamp,
  family = "multinomial2",
  data = df_laying_group,
  prior = prior1,
  nitt = nitts,
  burnin = burn,
  thin = thin,
  verbose = F)
summary(m_laying_group_int_c1)
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
set.seed(822)
m_laying_group_int_c2 <- MCMCglmm(cbind(no.eggs, remaining - no.eggs) ~ CheatinCamp-1 + ZnoCoops + ZnoOther + ZDay +
    CheatinCamp:ZDay-1 + ZnoCoops:ZDay,
  random = ~  year + camp + yearcamp,
  family = "multinomial2",
  data = df_laying_group,
  prior = prior1,
  nitt = nitts,
  burnin = burn,
  thin = thin,
  verbose = F)
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#
set.seed(823)
m_laying_group_int_c3 <- MCMCglmm(cbind(no.eggs, remaining - no.eggs) ~ CheatinCamp-1 + ZnoCoops + ZnoOther + ZDay +
    CheatinCamp:ZDay-1 + ZnoCoops:ZDay,
  random = ~  year + camp + yearcamp,
  family = "multinomial2",
  data = df_laying_group,
  prior = prior1,
  nitt = nitts,
  burnin = burn,
  thin = thin,
  verbose = F)
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#
#Check convergence of models#
par(mar=c(2,2,2,2))#To avoid errors with figure margins being too large
m_laying_group_int_sol<-mcmc.list(m_laying_group_int_c1$Sol,m_laying_group_int_c2$Sol,m_laying_group_int_c3$Sol);plot(m_laying_group_int_sol)
m_laying_group_int_vcv<-mcmc.list(m_laying_group_int_c1$VCV,m_laying_group_int_c2$VCV,m_laying_group_int_c3$VCV);plot(m_laying_group_int_vcv)
gelman.diag(m_laying_group_int_sol)
gelman.diag(m_laying_group_int_vcv, multivariate = FALSE)
#
#Quadratics
set.seed(831)
m_laying_group_quad_c1 <- MCMCglmm(cbind(no.eggs, remaining - no.eggs) ~ CheatinCamp-1 + ZnoCoops + ZnoOther + ZDay +
    CheatinCamp:ZDay-1 + ZnoCoops:ZDay +
    CheatinCamp:ZDay2-1 + ZnoCoops:ZDay2,
  random = ~  year + camp + yearcamp,
  family = "multinomial2",
  data = df_laying_group,
  prior = prior1,
  nitt = nitts,
  burnin = burn,
  thin = thin,
  verbose = F)
summary(m_laying_group_quad_c1)
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#
set.seed(832)
m_laying_group_quad_c2 <- MCMCglmm(cbind(no.eggs, remaining - no.eggs) ~ CheatinCamp-1 + ZnoCoops + ZnoOther + ZDay +
    CheatinCamp:ZDay-1 + ZnoCoops:ZDay +
    CheatinCamp:ZDay2-1 + ZnoCoops:ZDay2,
  random = ~  year + camp + yearcamp,
  family = "multinomial2",
  data = df_laying_group,
  prior = prior1,
  nitt = nitts,
  burnin = burn,
  thin = thin,
  verbose = F)
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#
set.seed(833)
m_laying_group_quad_c3 <- MCMCglmm(cbind(no.eggs, remaining - no.eggs) ~ CheatinCamp-1 + ZnoCoops + ZnoOther + ZDay +
    CheatinCamp:ZDay-1 + ZnoCoops:ZDay +
    CheatinCamp:ZDay2-1 + ZnoCoops:ZDay2,
  random = ~  year + camp + yearcamp,
  family = "multinomial2",
  data = df_laying_group,
  prior = prior1,
  nitt = nitts,
  burnin = burn,
  thin = thin,
  verbose = F)
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#
#
#Check convergence of models#
par(mar=c(2,2,2,2))#To avoid errors with figure margins being too large
m_laying_group_quad_sol<-mcmc.list(m_laying_group_quad_c1$Sol,m_laying_group_quad_c2$Sol,m_laying_group_quad_c3$Sol);plot(m_laying_group_quad_sol)
m_laying_group_quad_vcv<-mcmc.list(m_laying_group_quad_c1$VCV,m_laying_group_quad_c2$VCV,m_laying_group_quad_c3$VCV);plot(m_laying_group_quad_vcv)
gelman.diag(m_laying_group_quad_sol)
gelman.diag(m_laying_group_quad_vcv, multivariate = FALSE)
#
###m_laying_window: Egg laying probability through time (Model 8, Extended Data Table 8)----
#***************************************************
#
#Start with main effects
set.seed(911)
m_laying_window_main_c1 <- MCMCglmm(probEgg_2days ~ BS-1 + Znofemales,
  random=~year + camp + yearcamp + id,
  family ="categorical", data = df_laying_window,
  prior = prior3, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE,trunc=T)
summary(m_laying_window_main_c1)
#
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#
set.seed(912)
m_laying_window_main_c2 <- MCMCglmm(probEgg_2days ~ BS-1 + Znofemales,
  random=~year + camp + yearcamp + id,
  family ="categorical", data = df_laying_window,
  prior = prior3, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE,trunc=T)
#
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#
set.seed(913)
m_laying_window_main_c3 <- MCMCglmm(probEgg_2days ~ BS-1 + Znofemales,
  random=~year + camp + yearcamp + id,
  family ="categorical", data = df_laying_window,
  prior = prior3, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE,trunc=T)
#
#Check convergence of models#
par(mar=c(2,2,2,2))#To avoid errors with figure margins being too large
m_laying_window_main_sol<-mcmc.list(m_laying_window_main_c1$Sol,m_laying_window_main_c2$Sol,m_laying_window_main_c3$Sol);plot(m_laying_window_main_sol)
m_laying_window_main_vcv<-mcmc.list(m_laying_window_main_c1$VCV,m_laying_window_main_c2$VCV,m_laying_window_main_c3$VCV);plot(m_laying_window_main_vcv)
gelman.diag(m_laying_window_main_sol)
gelman.diag(m_laying_window_main_vcv, multivariate = FALSE)
#
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#
#Then run two-way interactions
set.seed(921)
m_laying_window_int_c1 <- MCMCglmm(probEgg_2days ~ BS-1 + Znofemales +
    BS:ZWindow_2days,
  random=~year + camp + yearcamp + id,
  family ="categorical", data = df_laying_window,
  prior = prior3, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE,trunc=T)
summary(m_laying_window_int_c1)
#
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#
set.seed(922)
m_laying_window_int_c2 <- MCMCglmm(probEgg_2days ~ BS-1 + Znofemales +
    BS:ZWindow_2days,
  random=~year + camp + yearcamp + id,
  family ="categorical", data = df_laying_window,
  prior = prior3, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE,trunc=T)
#
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#
set.seed(923)
m_laying_window_int_c3 <- MCMCglmm(probEgg_2days ~ BS-1 + Znofemales +
    BS:ZWindow_2days,
  random=~year + camp + yearcamp + id,
  family ="categorical", data = df_laying_window,
  prior = prior3, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE,trunc=T)
#
#Check convergence of models#
par(mar=c(2,2,2,2))#To avoid errors with figure margins being too large
m_laying_window_int_sol<-mcmc.list(m_laying_window_int_c1$Sol,m_laying_window_int_c2$Sol,m_laying_window_int_c3$Sol);plot(m_laying_window_int_sol)
m_laying_window_int_vcv<-mcmc.list(m_laying_window_int_c1$VCV,m_laying_window_int_c2$VCV,m_laying_window_int_c3$VCV);plot(m_laying_window_int_vcv)
gelman.diag(m_laying_window_int_sol)
gelman.diag(m_laying_window_int_vcv, multivariate = FALSE)
#
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#
#And finally run interactions with quadratics
set.seed(931)
m_laying_window_quad_c1 <- MCMCglmm(probEgg_2days ~ BS-1 + Znofemales +
    BS:ZWindow_2days +
    BS:ZWindow_2days2,
  random=~year + camp + yearcamp + id,
  family ="categorical", data = df_laying_window,
  prior = prior3, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE,trunc=T)
summary(m_laying_window_quad_c1)
#
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#
set.seed(932)
m_laying_window_quad_c2 <- MCMCglmm(probEgg_2days ~ BS-1 + Znofemales +
    BS:ZWindow_2days +
    BS:ZWindow_2days2,
  random=~year + camp + yearcamp + id,
  family ="categorical", data = df_laying_window,
  prior = prior3, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE,trunc=T)
#
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#
set.seed(933)
m_laying_window_quad_c3 <- MCMCglmm(probEgg_2days ~ BS-1 + Znofemales +
    BS:ZWindow_2days +
    BS:ZWindow_2days2,
  random=~year + camp + yearcamp + id,
  family ="categorical", data = df_laying_window,
  prior = prior3, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE,trunc=T)
#
#Check convergence of models#
par(mar=c(2,2,2,2))#To avoid errors with figure margins being too large
m_laying_window_quad_sol<-mcmc.list(m_laying_window_quad_c1$Sol,m_laying_window_quad_c2$Sol,m_laying_window_quad_c3$Sol);plot(m_laying_window_quad_sol)
m_laying_window_quad_vcv<-mcmc.list(m_laying_window_quad_c1$VCV,m_laying_window_quad_c2$VCV,m_laying_window_quad_c3$VCV);plot(m_laying_window_quad_vcv)
gelman.diag(m_laying_window_quad_sol)
gelman.diag(m_laying_window_quad_vcv, multivariate = FALSE)
#
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#******************************************
##Aggression between breeding strategies----
#******************************************
#
###m_aggression: Cooperator aggression (Model 9, Extended Data Table 9)----
#*************************************
#
#Start with main effects
set.seed(1011)
m_aggression_main_c1 <- MCMCglmm(no.ints  ~ BSvic-1 + ZtotObstime + ZnoBSvic + Znofemales,
                                 random= ~ year + camp + yearcamp,
                                 family ="poisson", data = df_aggression,
                                 prior = prior1, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
summary(m_aggression_main_c1)
set.seed(1012)
m_aggression_main_c2 <- MCMCglmm(no.ints  ~ BSvic-1 + ZtotObstime + ZnoBSvic + Znofemales,
                                 random= ~ year + camp + yearcamp,
                                 family ="poisson", data = df_aggression,
                                 prior = prior1, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
set.seed(1013)
m_aggression_main_c3 <- MCMCglmm(no.ints  ~ BSvic-1 + ZtotObstime + ZnoBSvic + Znofemales,
                                 random= ~ year + camp + yearcamp,
                                 family ="poisson", data = df_aggression,
                                 prior = prior1, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
#
#Check convergence of models#
par(mar=c(2,2,2,2))#To avoid errors with figure margins being too large
m_aggression_main_sol<-mcmc.list(m_aggression_main_c1$Sol,m_aggression_main_c2$Sol,m_aggression_main_c3$Sol);plot(m_aggression_main_sol)
m_aggression_main_vcv<-mcmc.list(m_aggression_main_c1$VCV,m_aggression_main_c2$VCV,m_aggression_main_c3$VCV);plot(m_aggression_main_vcv)
gelman.diag(m_aggression_main_sol)
gelman.diag(m_aggression_main_vcv, multivariate = FALSE)
#
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#
#Then run two-way interactions
set.seed(1021)
m_aggression_int_c1 <- MCMCglmm(no.ints  ~ BSvic-1 + ZtotObstime + ZnoBSvic + Znofemales +
                                  Znofemales:BSvic,
                                random= ~ year + camp + yearcamp,
                                family ="poisson", data = df_aggression,
                                prior = prior1, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
summary(m_aggression_int_c1)
set.seed(1022)
m_aggression_int_c2 <- MCMCglmm(no.ints  ~ BSvic-1 + ZtotObstime + ZnoBSvic + Znofemales +
                                  Znofemales:BSvic,
                                random= ~ year + camp + yearcamp,
                                family ="poisson", data = df_aggression,
                                prior = prior1, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
set.seed(1023)
m_aggression_int_c3 <- MCMCglmm(no.ints  ~ BSvic-1 + ZtotObstime + ZnoBSvic + Znofemales +
                                  Znofemales:BSvic,
                                random= ~ year + camp + yearcamp,
                                family ="poisson", data = df_aggression,
                                prior = prior1, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
#
#Check convergence of models#
par(mar=c(2,2,2,2))#To avoid errors with figure margins being too large
m_aggression_int_sol<-mcmc.list(m_aggression_int_c1$Sol,m_aggression_int_c2$Sol,m_aggression_int_c3$Sol);plot(m_aggression_int_sol)
m_aggression_int_vcv<-mcmc.list(m_aggression_int_c1$VCV,m_aggression_int_c2$VCV,m_aggression_int_c3$VCV);plot(m_aggression_int_vcv)
gelman.diag(m_aggression_int_sol)
gelman.diag(m_aggression_int_vcv, multivariate = FALSE)
#
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#****************************
##Investment in incubation----
#****************************
#
###m_incubation_group: Group-level incubation effort in relation to presence of cheats (Model 10, Extended Data Table 10)----
#***************************************************************************
#Main effects
set.seed(1111)
m_incubation_group_main_c1 <- MCMCglmm(cbind(Finc_N1, totobstime - Finc_N1) ~ FCheatinCamp-1,
  random= ~ year + camp,
  family ="multinomial2", data = df_incubation_group,
  prior = prior5, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
summary(m_incubation_group_main_c1)
set.seed(1112)
m_incubation_group_main_c2 <- MCMCglmm(cbind(Finc_N1, totobstime - Finc_N1) ~ FCheatinCamp-1,
  random= ~ year + camp,
  family ="multinomial2", data = df_incubation_group,
  prior = prior5, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
set.seed(1113)
m_incubation_group_main_c3 <- MCMCglmm(cbind(Finc_N1, totobstime - Finc_N1) ~ FCheatinCamp-1,
  random= ~ year + camp,
  family ="multinomial2", data = df_incubation_group,
  prior = prior5, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
#Check convergence of models
par(mar=c(2,2,2,2))#To avoid errors with figure margins being too large
m_incubation_group_main_sol<-mcmc.list(m_incubation_group_main_c1$Sol,m_incubation_group_main_c2$Sol,m_incubation_group_main_c3$Sol)
plot(m_incubation_group_main_sol)
m_incubation_group_main_vcv<-mcmc.list(m_incubation_group_main_c1$VCV,m_incubation_group_main_c2$VCV,m_incubation_group_main_c3$VCV)
plot(m_incubation_group_main_vcv)
gelman.diag(m_incubation_group_main_sol)
gelman.diag(m_incubation_group_main_vcv, multivariate = FALSE)
#
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#
#Two-way interactions
set.seed(1121)
m_incubation_group_int_c1 <- MCMCglmm(cbind(Finc_N1, totobstime - Finc_N1) ~ FCheatinCamp-1 +
    FCheatinCamp:ZnoFCoops + FCheatinCamp:ZnoFOther,
  random= ~ year + camp,
  family ="multinomial2", data = df_incubation_group,
  prior = prior5, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
summary(m_incubation_group_int_c1)
set.seed(1122)
m_incubation_group_int_c2 <- MCMCglmm(cbind(Finc_N1, totobstime - Finc_N1) ~ FCheatinCamp-1 +
    FCheatinCamp:ZnoFCoops + FCheatinCamp:ZnoFOther,
  random= ~ year + camp,
  family ="multinomial2", data = df_incubation_group,
  prior = prior5, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
set.seed(1123)
m_incubation_group_int_c3 <- MCMCglmm(cbind(Finc_N1, totobstime - Finc_N1) ~ FCheatinCamp-1 +
    FCheatinCamp:ZnoFCoops + FCheatinCamp:ZnoFOther,
  random= ~ year + camp,
  family ="multinomial2", data = df_incubation_group,
  prior = prior5, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
#Check convergence of models
par(mar=c(2,2,2,2))#To avoid errors with figure margins being too large
m_incubation_group_int_sol<-mcmc.list(m_incubation_group_int_c1$Sol,m_incubation_group_int_c2$Sol,m_incubation_group_int_c3$Sol)
plot(m_incubation_group_int_sol)
m_incubation_group_int_vcv<-mcmc.list(m_incubation_group_int_c1$VCV,m_incubation_group_int_c2$VCV,m_incubation_group_int_c3$VCV)
plot(m_incubation_group_int_vcv)
gelman.diag(m_incubation_group_int_sol)
gelman.diag(m_incubation_group_int_vcv, multivariate = FALSE)
#
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
###m_incubation_time: Group-level incubation effort over time (Model 11, Extended Data Table 11)----
####**************************************************
start_time <- Sys.time()
#Main effects
set.seed(1211)
m_incubation_time_main_c1 <- MCMCglmm(cbind(sitdur_date, obsdur_date - sitdur_date) ~ CheatinCamp-1 + ZnoCoops + ZnoOther + ZDay,
  random = ~  year + camp + yearcamp,
  family = "multinomial2",
  data = df_incubation_time,
  prior = prior1,
  nitt = nitts,
  burnin = burn,
  thin = thin,
  verbose = F)
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
summary(m_incubation_time_main_c1)
#
set.seed(1212)
m_incubation_time_main_c2 <- MCMCglmm(cbind(sitdur_date, obsdur_date - sitdur_date) ~ CheatinCamp-1 + ZnoCoops + ZnoOther + ZDay,
  random = ~  year + camp + yearcamp,
  family = "multinomial2",
  data = df_incubation_time,
  prior = prior1,
  nitt = nitts,
  burnin = burn,
  thin = thin,
  verbose = F)
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#
set.seed(1213)
m_incubation_time_main_c3 <- MCMCglmm(cbind(sitdur_date, obsdur_date - sitdur_date) ~ CheatinCamp-1 + ZnoCoops + ZnoOther + ZDay,
  random = ~  year + camp + yearcamp,
  family = "multinomial2",
  data = df_incubation_time,
  prior = prior1,
  nitt = nitts,
  burnin = burn,
  thin = thin,
  verbose = F)
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#
#Check convergence of models#
par(mar=c(2,2,2,2))#To avoid errors with figure margins being too large
m_incubation_time_main_sol<-mcmc.list(m_incubation_time_main_c1$Sol,m_incubation_time_main_c2$Sol,m_incubation_time_main_c3$Sol);plot(m_incubation_time_main_sol)
m_incubation_time_main_vcv<-mcmc.list(m_incubation_time_main_c1$VCV,m_incubation_time_main_c2$VCV,m_incubation_time_main_c3$VCV);plot(m_incubation_time_main_vcv)
gelman.diag(m_incubation_time_main_sol)
gelman.diag(m_incubation_time_main_vcv, multivariate = FALSE)
#
#Two-way interactions
set.seed(1221)
m_incubation_time_int_c1 <- MCMCglmm(cbind(sitdur_date, obsdur_date - sitdur_date) ~ CheatinCamp-1 + ZnoCoops + ZnoOther + ZDay +
    CheatinCamp:ZDay-1 + ZnoCoops:ZDay,
  random = ~  year + camp + yearcamp,
  family = "multinomial2",
  data = df_incubation_time,
  prior = prior1,
  nitt = nitts,
  burnin = burn,
  thin = thin,
  verbose = F)
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
summary(m_incubation_time_int_c1)
#
set.seed(1222)
m_incubation_time_int_c2 <- MCMCglmm(cbind(sitdur_date, obsdur_date - sitdur_date) ~ CheatinCamp-1 + ZnoCoops + ZnoOther + ZDay +
    CheatinCamp:ZDay-1 + ZnoCoops:ZDay,
  random = ~  year + camp + yearcamp,
  family = "multinomial2",
  data = df_incubation_time,
  prior = prior1,
  nitt = nitts,
  burnin = burn,
  thin = thin,
  verbose = F)
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#
set.seed(1223)
m_incubation_time_int_c3 <- MCMCglmm(cbind(sitdur_date, obsdur_date - sitdur_date) ~ CheatinCamp-1 + ZnoCoops + ZnoOther + ZDay +
    CheatinCamp:ZDay-1 + ZnoCoops:ZDay,
  random = ~  year + camp + yearcamp,
  family = "multinomial2",
  data = df_incubation_time,
  prior = prior1,
  nitt = nitts,
  burnin = burn,
  thin = thin,
  verbose = F)
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#Quadratics
set.seed(1231)
m_incubation_time_quad_c1 <- MCMCglmm(cbind(sitdur_date, obsdur_date - sitdur_date) ~ CheatinCamp-1 + ZnoCoops + ZnoOther + ZDay +
    CheatinCamp:ZDay-1 + ZnoCoops:ZDay + 
    CheatinCamp:ZDay2-1 + ZnoCoops:ZDay2,
  random = ~  year + camp + yearcamp,
  family = "multinomial2",
  data = df_incubation_time,
  prior = prior1,
  nitt = nitts,
  burnin = burn,
  thin = thin,
  verbose = F)
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
summary(m_incubation_time_quad_c1)
#
set.seed(1232)
m_incubation_time_quad_c2 <- MCMCglmm(cbind(sitdur_date, obsdur_date - sitdur_date) ~ CheatinCamp-1 + ZnoCoops + ZnoOther + ZDay +
    CheatinCamp:ZDay-1 + ZnoCoops:ZDay + 
    CheatinCamp:ZDay2-1 + ZnoCoops:ZDay2,
  random = ~  year + camp + yearcamp,
  family = "multinomial2",
  data = df_incubation_time,
  prior = prior1,
  nitt = nitts,
  burnin = burn,
  thin = thin,
  verbose = F)
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#
set.seed(1233)
m_incubation_time_quad_c3 <- MCMCglmm(cbind(sitdur_date, obsdur_date - sitdur_date) ~ CheatinCamp-1 + ZnoCoops + ZnoOther + ZDay +
    CheatinCamp:ZDay-1 + ZnoCoops:ZDay + 
    CheatinCamp:ZDay2-1 + ZnoCoops:ZDay2,
  random = ~  year + camp + yearcamp,
  family = "multinomial2",
  data = df_incubation_time,
  prior = prior1,
  nitt = nitts,
  burnin = burn,
  thin = thin,
  verbose = F)
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
end_time <- Sys.time()
#How long did it take to run all the m_incubation_time models?
elapsed <- as.numeric(difftime(end_time, start_time, units="secs"))
#Write function to get time as hh:mm
format_hours <- function(decimal_hours) {
  hours <- floor(decimal_hours)
  minutes <- round((decimal_hours - hours) * 60)
  sprintf("%d:%02d", hours, minutes)
}
format_hours(elapsed/3600)#
###m_incubation_individual: Incubation effort in relation to cheating, group composition and reproductive success, square root-converting RS (Model 12, Extended Data Table 12)----
#***********************************************************************************************************************
#
#Main effects
set.seed(1311)
m_incubation_individual_main_c1 <- MCMCglmm(cbind(sitdurbird, totobsdurYear - sitdurbird) ~ FCheatinCamp-1,
  random= ~ year + camp + yearcamp + id,
  family ="multinomial2", data = df_incubation_individual,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
summary(m_incubation_individual_main_c1)
set.seed(1312)
m_incubation_individual_main_c2 <- MCMCglmm(cbind(sitdurbird, totobsdurYear - sitdurbird) ~ FCheatinCamp-1,
  random= ~ year + camp + yearcamp + id,
  family ="multinomial2", data = df_incubation_individual,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
set.seed(1313)
m_incubation_individual_main_c3 <- MCMCglmm(cbind(sitdurbird, totobsdurYear - sitdurbird) ~ FCheatinCamp-1,
  random= ~ year + camp + yearcamp + id,
  family ="multinomial2", data = df_incubation_individual,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
#
#Check convergence of models#
par(mar=c(2,2,2,2))#To avoid errors with figure margins being too large
m_incubation_individual_main_sol<-mcmc.list(m_incubation_individual_main_c1$Sol,m_incubation_individual_main_c2$Sol,m_incubation_individual_main_c3$Sol);plot(m_incubation_individual_main_sol)
m_incubation_individual_main_vcv<-mcmc.list(m_incubation_individual_main_c1$VCV,m_incubation_individual_main_c2$VCV,m_incubation_individual_main_c3$VCV);plot(m_incubation_individual_main_vcv)
gelman.diag(m_incubation_individual_main_sol)
gelman.diag(m_incubation_individual_main_vcv, multivariate = FALSE)
#
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#
#Two-way interactions
set.seed(1321)
m_incubation_individual_int_c1 <- MCMCglmm(cbind(sitdurbird, totobsdurYear - sitdurbird) ~ FCheatinCamp-1 +
    FCheatinCamp:sqrt(relativeRSInc4) + FCheatinCamp:ZnoFCoops + FCheatinCamp:ZnoFOther,
  random= ~ year + camp + yearcamp + id,
  family ="multinomial2", data = df_incubation_individual,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
summary(m_incubation_individual_int_c1)
set.seed(1322)
m_incubation_individual_int_c2 <- MCMCglmm(cbind(sitdurbird, totobsdurYear - sitdurbird) ~ FCheatinCamp-1 +
    FCheatinCamp:sqrt(relativeRSInc4) + FCheatinCamp:ZnoFCoops + FCheatinCamp:ZnoFOther,
  random= ~ year + camp + yearcamp + id,
  family ="multinomial2", data = df_incubation_individual,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
set.seed(1323)
m_incubation_individual_int_c3 <- MCMCglmm(cbind(sitdurbird, totobsdurYear - sitdurbird) ~ FCheatinCamp-1 +
    FCheatinCamp:sqrt(relativeRSInc4) + FCheatinCamp:ZnoFCoops + FCheatinCamp:ZnoFOther,
  random= ~ year + camp + yearcamp + id,
  family ="multinomial2", data = df_incubation_individual,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
#
#Check convergence of models#
par(mar=c(2,2,2,2))#To avoid errors with figure margins being too large
m_incubation_individual_int_sol<-mcmc.list(m_incubation_individual_int_c1$Sol,m_incubation_individual_int_c2$Sol,m_incubation_individual_int_c3$Sol);plot(m_incubation_individual_int_sol)
m_incubation_individual_int_vcv<-mcmc.list(m_incubation_individual_int_c1$VCV,m_incubation_individual_int_c2$VCV,m_incubation_individual_int_c3$VCV);plot(m_incubation_individual_int_vcv)
gelman.diag(m_incubation_individual_int_sol)
gelman.diag(m_incubation_individual_int_vcv, multivariate = FALSE)
#
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#
#Three-way interactions
set.seed(1331)
m_incubation_individual_int3_c1 <- MCMCglmm(cbind(sitdurbird, totobsdurYear - sitdurbird) ~ FCheatinCamp-1 +
    FCheatinCamp:sqrt(relativeRSInc4) + FCheatinCamp:ZnoFCoops + FCheatinCamp:ZnoFOther +
    FCheatinCamp:sqrt(relativeRSInc4):ZnoFCoops,
  random= ~ year + camp + yearcamp + id,
  family ="multinomial2", data = df_incubation_individual,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
summary(m_incubation_individual_int3_c1)
set.seed(1332)
m_incubation_individual_int3_c2 <- MCMCglmm(cbind(sitdurbird, totobsdurYear - sitdurbird) ~ FCheatinCamp-1 +
    FCheatinCamp:sqrt(relativeRSInc4) + FCheatinCamp:ZnoFCoops + FCheatinCamp:ZnoFOther +
    FCheatinCamp:sqrt(relativeRSInc4):ZnoFCoops,
  random= ~ year + camp + yearcamp + id,
  family ="multinomial2", data = df_incubation_individual,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
set.seed(1333)
m_incubation_individual_int3_c3 <- MCMCglmm(cbind(sitdurbird, totobsdurYear - sitdurbird) ~ FCheatinCamp-1 +
    FCheatinCamp:sqrt(relativeRSInc4) + FCheatinCamp:ZnoFCoops + FCheatinCamp:ZnoFOther +
    FCheatinCamp:sqrt(relativeRSInc4):ZnoFCoops,
  random= ~ year + camp + yearcamp + id,
  family ="multinomial2", data = df_incubation_individual,
  prior = prior4, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
#
#Check convergence of models
par(mar=c(2,2,2,2))#To avoid errors with figure margins being too large
m_incubation_individual_int3_sol<-mcmc.list(m_incubation_individual_int3_c1$Sol,m_incubation_individual_int3_c2$Sol,m_incubation_individual_int3_c3$Sol);plot(m_incubation_individual_int3_sol)
m_incubation_individual_int3_vcv<-mcmc.list(m_incubation_individual_int3_c1$VCV,m_incubation_individual_int3_c2$VCV,m_incubation_individual_int3_c3$VCV);plot(m_incubation_individual_int3_vcv)
gelman.diag(m_incubation_individual_int3_sol)
gelman.diag(m_incubation_individual_int3_vcv, multivariate = FALSE)#Year not converging very well
#
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#****************************************
##Group-level consequences of cheating----
#****************************************
#
###m_nest_failure and m_hatch_success: Hatching success in relation to presence of cheats (Models 13 and 13b, Extended Data Table 13)----
##*****************************************************************
#Hatching success is modelled in two parts, because the effect of cheating acts on whether a nest
#hatched any eggs at all rather than on the proportion of eggs that hatched. m_nest_failure models
#complete failure (Extended Data Table 11) and m_hatch_success models the proportion hatching among nests that
#hatched at least one egg (Extended Data Table 11b).
#The number of cooperators is not fitted as a quadratic: every group containing a cheat and a single
#cooperator failed completely, and because those groups lie at the end of the range of group sizes they
#dominate the estimated curvature.
#The proportion of eggs broken is deliberately not fitted: egg breakage is itself a route to nest
#failure, so conditioning on it would remove part of the effect being estimated.
#
#Data for the two parts: all groups, and those that hatched at least one egg
df_nest_failure <- df_hatching
df_nest_failure$failed <- as.numeric(df_nest_failure$no.HatchedInc == 0)
df_hatch_success <- df_nest_failure[df_nest_failure$failed == 0,]
#
##Part 1: did the nest hatch any eggs at all? (Model 13)
#Start with main effects
set.seed(1411)
m_nest_failure_main_c1 <- MCMCglmm(failed ~ FCheatinCamp-1 + ZnoFOther,
  random=~year + camp, family ="categorical", data = df_nest_failure,
  prior = prior8, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE, trunc = TRUE)
summary(m_nest_failure_main_c1)
set.seed(1412)
m_nest_failure_main_c2 <- MCMCglmm(failed ~ FCheatinCamp-1 + ZnoFOther,
  random=~year + camp, family ="categorical", data = df_nest_failure,
  prior = prior8, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE, trunc = TRUE)
set.seed(1413)
m_nest_failure_main_c3 <- MCMCglmm(failed ~ FCheatinCamp-1 + ZnoFOther,
  random=~year + camp, family ="categorical", data = df_nest_failure,
  prior = prior8, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE, trunc = TRUE)
#
#Check convergence of models
m_nest_failure_main_sol<-mcmc.list(m_nest_failure_main_c1$Sol,m_nest_failure_main_c2$Sol,m_nest_failure_main_c3$Sol)
m_nest_failure_main_vcv<-mcmc.list(m_nest_failure_main_c1$VCV,m_nest_failure_main_c2$VCV,m_nest_failure_main_c3$VCV)
gelman.diag(m_nest_failure_main_sol)
gelman.diag(m_nest_failure_main_vcv, multivariate = FALSE)
#
#Then add the interaction between presence of cheats and number of cooperators
set.seed(1421)
m_nest_failure_int_c1 <- MCMCglmm(failed ~ FCheatinCamp-1 + ZnoFOther + FCheatinCamp:ZnoFCoops,
  random=~year + camp, family ="categorical", data = df_nest_failure,
  prior = prior8, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE, trunc = TRUE)
summary(m_nest_failure_int_c1)
set.seed(1422)
m_nest_failure_int_c2 <- MCMCglmm(failed ~ FCheatinCamp-1 + ZnoFOther + FCheatinCamp:ZnoFCoops,
  random=~year + camp, family ="categorical", data = df_nest_failure,
  prior = prior8, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE, trunc = TRUE)
set.seed(1423)
m_nest_failure_int_c3 <- MCMCglmm(failed ~ FCheatinCamp-1 + ZnoFOther + FCheatinCamp:ZnoFCoops,
  random=~year + camp, family ="categorical", data = df_nest_failure,
  prior = prior8, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE, trunc = TRUE)
#
#Check convergence of models
m_nest_failure_int_sol<-mcmc.list(m_nest_failure_int_c1$Sol,m_nest_failure_int_c2$Sol,m_nest_failure_int_c3$Sol)
m_nest_failure_int_vcv<-mcmc.list(m_nest_failure_int_c1$VCV,m_nest_failure_int_c2$VCV,m_nest_failure_int_c3$VCV)
gelman.diag(m_nest_failure_int_sol)
gelman.diag(m_nest_failure_int_vcv, multivariate = FALSE)
#
#Table 13 reports a negative posterior mode for the year variance.
#Check why this happens.
#1)Collect the three independently fitted chains of the nest-failure model.
fits = list( c1 = m_nest_failure_int_c1,
              c2 = m_nest_failure_int_c2,
              c3 = m_nest_failure_int_c3)
#2)Apply the same checks to each chain and display one row per chain.
t(vapply(fits, function(m) {
  #Extract posterior draws of the year variance.
  #drop = FALSE preserves the matrix structure for posterior.mode().
  x = m$VCV[, "year", drop = FALSE]
  #Smallest sampled variance: should be positive.
  c(minimum = min(x),
    #Number of negative variance draws: should be zero.
    negative_draws = sum(x < 0),
    #Number of NA, NaN or infinite draws: should be zero.
    nonfinite_draws = sum(!is.finite(x)),
    #Recalculate the posterior mode using the same function as in reporting.
    #A negative mode despite entirely positive draws identifies the issue
    #as arising during mode estimation, rather than posterior sampling.
    mode = as.numeric(MCMCglmm::posterior.mode(x)))
}, numeric(4)))
#Positive minima and zero negative/non-finite draws confirm that the sampled
#variances are valid. The negative mode is therefore a summarisation
#artefact, not a negative fitted variance.
#
#
#Does the effect of cheats run through egg breakage? Refitting with the proportion of eggs broken
#included does not weaken it, so breakage and cheating are separate routes to nest failure
set.seed(1409)
m_nest_failure_broken <- MCMCglmm(failed ~ FCheatinCamp-1 + ZnoFOther + logit_prop.BrokenInc + FCheatinCamp:ZnoFCoops,
  random=~year + camp, family ="categorical", data = df_nest_failure,
  prior = prior8, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE, trunc = TRUE)
summary(m_nest_failure_broken)
#
##Part 2: of the nests that hatched eggs, what proportion hatched? (Model 13b)
#Start with main effects
set.seed(1511)
m_hatch_success_main_c1 <- MCMCglmm(cbind(no.HatchedInc, no.eggsCampInc - no.HatchedInc) ~ FCheatinCamp-1 + ZnoFOther,
  random=~year + camp, family ="multinomial2", data = df_hatch_success,
  prior = prior5, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
summary(m_hatch_success_main_c1)
set.seed(1512)
m_hatch_success_main_c2 <- MCMCglmm(cbind(no.HatchedInc, no.eggsCampInc - no.HatchedInc) ~ FCheatinCamp-1 + ZnoFOther,
  random=~year + camp, family ="multinomial2", data = df_hatch_success,
  prior = prior5, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
set.seed(1513)
m_hatch_success_main_c3 <- MCMCglmm(cbind(no.HatchedInc, no.eggsCampInc - no.HatchedInc) ~ FCheatinCamp-1 + ZnoFOther,
  random=~year + camp, family ="multinomial2", data = df_hatch_success,
  prior = prior5, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
#
#Check convergence of models
m_hatch_success_main_sol<-mcmc.list(m_hatch_success_main_c1$Sol,m_hatch_success_main_c2$Sol,m_hatch_success_main_c3$Sol)
m_hatch_success_main_vcv<-mcmc.list(m_hatch_success_main_c1$VCV,m_hatch_success_main_c2$VCV,m_hatch_success_main_c3$VCV)
gelman.diag(m_hatch_success_main_sol)
gelman.diag(m_hatch_success_main_vcv, multivariate = FALSE)
#
#Then add the interaction between presence of cheats and number of cooperators
set.seed(1521)
m_hatch_success_int_c1 <- MCMCglmm(cbind(no.HatchedInc, no.eggsCampInc - no.HatchedInc) ~ FCheatinCamp-1 + ZnoFOther + FCheatinCamp:ZnoFCoops,
  random=~year + camp, family ="multinomial2", data = df_hatch_success,
  prior = prior5, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
summary(m_hatch_success_int_c1)
set.seed(1522)
m_hatch_success_int_c2 <- MCMCglmm(cbind(no.HatchedInc, no.eggsCampInc - no.HatchedInc) ~ FCheatinCamp-1 + ZnoFOther + FCheatinCamp:ZnoFCoops,
  random=~year + camp, family ="multinomial2", data = df_hatch_success,
  prior = prior5, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
set.seed(1523)
m_hatch_success_int_c3 <- MCMCglmm(cbind(no.HatchedInc, no.eggsCampInc - no.HatchedInc) ~ FCheatinCamp-1 + ZnoFOther + FCheatinCamp:ZnoFCoops,
  random=~year + camp, family ="multinomial2", data = df_hatch_success,
  prior = prior5, nitt=nitts, burnin=burn, thin=thin, verbose = FALSE)
#
#Check convergence of models
m_hatch_success_int_sol<-mcmc.list(m_hatch_success_int_c1$Sol,m_hatch_success_int_c2$Sol,m_hatch_success_int_c3$Sol)
m_hatch_success_int_vcv<-mcmc.list(m_hatch_success_int_c1$VCV,m_hatch_success_int_c2$VCV,m_hatch_success_int_c3$VCV)
gelman.diag(m_hatch_success_int_sol)
gelman.diag(m_hatch_success_int_vcv, multivariate = FALSE)
#
#
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#

###m_activity_nocheat and m_activity_cheat: Fragmented regression of group-level activity allocation (egg-laying vs incubation) across time (Model 14, Extended Data Table 14)----
####*********************************************************************************************************
#
# Fit separately by group, with yearcamp as a random effect
fit_activity_nocheat  = lme(logit_alloc ~ ZDay + ZnoCoops + ZnoOther, random = ~ 1 | yearcamp,
  data = subset(df_activity_allocation, CheatinCamp == "No"))
m_activity_nocheat  = segmented.lme(fit_activity_nocheat, seg.Z = ~ ZDay, random = list(yearcamp = pdDiag(~1)),
  psi = 0)
summary(m_activity_nocheat)
#
fit_activity_cheat = lme(logit_alloc ~ ZDay + ZnoCoops + ZnoOther, random = ~ 1 | yearcamp,
  data = subset(df_activity_allocation, CheatinCamp == "Yes"))
m_activity_cheat = segmented.lme(fit_activity_cheat, seg.Z = ~ ZDay, random = list(yearcamp = pdDiag(~1)),
  psi = 0)
summary(m_activity_cheat)
#
confint(m_activity_nocheat)#Breakpoint CI for groups without cheats
confint(m_activity_cheat)#Breakpoint for groups with cheats
save.image("R_Output_files/CheatingModels.RData", compress = "xz")
#
#
