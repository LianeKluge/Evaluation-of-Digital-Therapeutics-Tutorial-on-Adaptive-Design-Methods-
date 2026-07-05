# Example 1 memodio DiGA based on https://diga.bfarm.de/de/verzeichnis/02866/fachkreise
#(5 July 2026) 
# authors: Liane Kluge, Laura Maaß, Felix Holl, Werner Brannath

# Please  note that for the code below for the two-stage adaptive design, the code 
# from the file Code_Eval DTx_Kluge_Maaß_Holl_Brannath must be run which, in turn,
# requires that the code from Part 1-4 (NOT part 5!) of 
# https://github.com/LianeKluge/Adaptive-Designs-in-Fast-Track-Registration-Processes
# is run

# Please note: all sample sizes are per group!

# analyse results from directory (pilot study) ---------------------------------

n_1 <- 70 # per group sample size of pilot study

# estimate se from results in DiGA directory:
L <- 2.941- 0.858
c <- qt(0.975,df=2*n_1-2)
se <- L*sqrt(n_1/2)/(2*c)

I_1 <- n_1/(2*se^2) # information of pilot study

delrel <- 1 # assumption on MCID

alpha_f <- min(1-pnorm(delrel*sqrt(I_1)), 0.05) # one-sided significance level
# that needs to be undercut for a successful conditional registration, resulting
# from the demand of having a relevant effect estimation and statistical
# tendency at the one-sided level of 0.05

alpha <- 0.025 # one-sided level for permanent registration

delta_apriori <- (qnorm(1-0.15)+qnorm(1-alpha_f))/sqrt(I_1) # calculate the 
# a priori (before pilot study) effect assumption under which the study could
# have been powered for 85%


est_1 <- 1.899 # observed effect after conducting pilot study

# separate Two-stage design-------------------------------------------------------

# Calculate for separate 2-stage design based on the estimate of the primary 
# endpoint the confirmatory study sample size, for 80%, 85% and 90% power:

n_2_80PercentPower_separat <- 
  ceiling((2*se^2*(qnorm(1-0.2)+qnorm(1-alpha))^2)/est_1^2) # sample size for 80% power at est_1

n_2_85PercentPower_separat <- 
  ceiling((2*se^2*(qnorm(1-0.15)+qnorm(1-alpha))^2)/est_1^2) # sample size for 85% power at est_1

n_2_90PercentPower_separat <- 
  ceiling((2*se^2*(qnorm(1-0.1)+qnorm(1-alpha))^2)/est_1^2) # sample size for 90% power at est_1

# One-stage Confirmatory design -------------------------------------------------

# calculate sample sizes for One-Stage Confirmatory Design:

# sample size for 85% power at delta_apriori
n_85PercentPower_at_delta_ariori_singleStage <- ceiling((2*se^2*(qnorm(1-0.15)+qnorm(1-alpha))^2)/delta_apriori^2)
# sample size for 80% power at delta_apriori
n_80PercentPower_at_delta_ariori_singleStage <- ceiling((2*se^2*(qnorm(1-0.2)+qnorm(1-alpha))^2)/delta_apriori^2)
# sample size for 85% power at delrel
n_85PercentPower_at_delrel_singleStage <- ceiling((2*se^2*(qnorm(1-0.15)+qnorm(1-alpha))^2)/delrel^2)
# sample size for 80% power at delrel
n_80PercentPower_at_delrel_singleStage <- ceiling((2*se^2*(qnorm(1-0.2)+qnorm(1-alpha))^2)/delrel^2)

# Adaptive Two-stage design ------------------------------------------------------
# run code from Code_Eval DTx_Kluge_Maaß_Holl_Brannath.R:

# caclulate minimal trial study sample size for adaptive design for overall power
# of 80% at delta_apriori and constant c_AF1 for conditional error function using code 
# from Code_Eval DTx_Kluge_Maaß_Holl_Brannath.R. This yields:
# n_2_min_Fisher_cp80, n_2_min_Fisher_cp85, n_2_min_Fisher_cp90 and c_AF1


A_F <- function(z1, c) { # Fishers product test
  A_def <- c/(1-pnorm(z1))
  return(min(A_def, 0.5))
}

A_f <- A_F(z1=sqrt(I_1)*est_1, c_AF1) # second stage significance level with
# c_AF1 from Code_Eval DTx_Kluge_Maaß_Holl_Brannath.R

# calculate for adaptive 2-stage design confirmatory stage sample size for
# 80%, 85% and 90% conditional power under overall power of 80%

n_2_80Percent_without_OverallPower <- ceiling(2*se^2*(qnorm(1-0.2)+qnorm(1-A_f))^2/est_1^2)
n_2_ConditionalPower80Percent_adaptiv_Fisher <- max(n_2_80Percent_without_OverallPower, n_2_min_Fisher_cp80)

n_2_85Percent_without_OverallPower <- ceiling(2*se^2*(qnorm(1-0.15)+qnorm(1-A_f))^2/est_1^2)
n_2_ConditionalPower85Percent_adaptiv_Fisher <- max(n_2_85Percent_without_OverallPower, n_2_min_Fisher_cp85)

n_2_90Percent_without_OverallPower <- ceiling(2*se^2*(qnorm(1-0.10)+qnorm(1-A_f))^2/est_1^2)
n_2_ConditionalPower90Percent_adaptiv_Fisher <- max(n_2_90Percent_without_OverallPower, n_2_min_Fisher_cp90)

# extra: separate 2-stage design with overall power control at delta_apriori-----


# caclulate minimal trial study sample size for separate design for overall power
# of 80%, 85% and 90% at delta_apriori using Code from Code_Eval DTx_Kluge_Maaß_Holl_Brannath. 
# This yields:
# n_2_min_separate_cp80. n_2_min_separate_cp85, n_2_min_separate_cp90 

# now calculate for separate 2-stage design confirmatory sample size
n_2_80PercentPower_separate_withOverallPower <- ceiling(max(n_2_80PercentPower_separat, n_2_min_separate_cp80))

n_2_85PercentPower_separate_withOverallPower <- ceiling(max(n_2_85PercentPower_separat, n_2_min_separate_cp85))

n_2_90PercentPower_separate_withOverallPower <- ceiling(max(n_2_90PercentPower_separat, n_2_min_separate_cp90))
