# Example 2: Adaptive designs with multiple endpoints
# based on DiGA Orthopy (see https://diga.bfarm.de/de/verzeichnis/00902/fachkreise)
#(5 July 2026)
# authors: Liane Kluge, Laura Maaß, Felix Holl, Werner Brannath

# Please  note that for the code below for the two-stage adaptive design, the code 
# from the file Code_Eval DTx_Kluge_Maaß_Holl_Brannath must be run which, in turn,
# requires that the code from Part 1-4 (NOT part 5!) of 
# https://github.com/LianeKluge/Adaptive-Designs-in-Fast-Track-Registration-Processes
# is run

# Please note: all sample sizes are per group!
# analyse results from directory (pilot study) ---------------------------------

n_1 <- 40 # per group sample size of pilot study

# estimate se for the endpoints pain and symptoms from results in DiGA directory:
L_pain <- 12.40 - (-1.7)
c_pain <- qt(0.975,df=2*n_1-2)
se_pain <- L_pain*sqrt(n_1/2)/(2*c_pain)
I_1_pain <- n_1/(2*se_pain^2) # information of pilot study

L_symptoms <- 14.29 - (0.42)
c_symptoms <- qt(0.975,df=2*n_1-2)
se_symptoms <- L_symptoms*sqrt(n_1/2)/(2*c_symptoms)
I_1_symptoms <- n_1/(2*se_symptoms^2) # information of pilot study


# calculate constant for conditional error function of Fisher:------------------
# PLEASE NOTE: run code from Code_INNO_Kluge_Maaß_Holl_Brannath.R to obtain c_AF2

c_AF2 <- 0.004351979

A_F_p <- function(p, c) { # Conditional error function of Fishers product test
  A_def <- c/p
  return(min(A_def, 0.5))
}

# calculate sample sizes for confirmatory study:--------------------------------
est_1_pain <- 5.67
est_1_symptoms <- 7.35
p_1_pain <- 0.0986/2 # one-sided p-value
p_1_symptoms <- 0.0377/2 # one-sided p-value

# Decision: KOOS sub scale symptoms is primary endpoint for confirmatory study 

# for separate study design: calculate sample size for confirmatory study under 
# the aim of 80%, 85% and 90% power for the primary endpoint symptoms:
n_2_80PercentPower_symptoms_separate <- 
  ceiling(2*se_symptoms^2*(qnorm(1-0.2)+qnorm(1-0.025))^2/est_1_symptoms^2)
n_2_80PercentPower_symptoms_separate

n_2_85PercentPower_symptoms_separate <- 
  ceiling(2*se_symptoms^2*(qnorm(1-0.15)+qnorm(1-0.025))^2/est_1_symptoms^2)
n_2_85PercentPower_symptoms_separate

n_2_90PercentPower_symptoms_separate <- 
  ceiling(2*se_symptoms^2*(qnorm(1-0.10)+qnorm(1-0.025))^2/est_1_symptoms^2)
n_2_90PercentPower_symptoms_separate

# for adaptive design: calculate sample size for second stage under 
# the aim of 80%, 85% and 90% conditional power for the primary endpoint 
# symptoms taking into account the combined hypothesis


A_f_symptoms <- A_F_p(p=p_1_symptoms, c_AF2)
A_f_combined <- A_F_p(min(2*p_1_pain, 2*p_1_symptoms), c_AF2)

n_2_80PercentPower_symptoms_adaptive <- 
  ceiling(2*se_symptoms^2*(qnorm(1-0.2)+qnorm(1-min(A_f_symptoms,A_f_combined)))^2/est_1_symptoms^2)
n_2_80PercentPower_symptoms_adaptive

n_2_85PercentPower_symptoms_adaptive <- 
  ceiling(2*se_symptoms^2*(qnorm(1-0.15)+qnorm(1-min(A_f_symptoms,A_f_combined)))^2/est_1_symptoms^2)
n_2_85PercentPower_symptoms_adaptive

n_2_90PercentPower_symptoms_adaptive <- 
  ceiling(2*se_symptoms^2*(qnorm(1-0.10)+qnorm(1-min(A_f_symptoms,A_f_combined)))^2/est_1_symptoms^2)
n_2_90PercentPower_symptoms_adaptive
