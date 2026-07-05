# Calculate values needed for example 1 and 2
# authors: Liane Kluge, Laura Maaß, Felix Holl, Werner Brannath
# PLEASE NOTE:
# FIRST run Part 1-4 (NOT part 5!) of the code from:
# https://github.com/LianeKluge/Adaptive-Designs-in-Fast-Track-Registration-Processes
# and then run code below!
###############################################################################

# common input parameters: -----------------------------------------------------
alpha <- 0.025
beta <- 0.2
tol <- 1 / 10 ^ 3
eta_f <- qnorm(1 - beta) + qnorm(1 - alpha)
stabConstIntegral <- 38.56802# for greater values, dnorm becomes zero
uppBoundSearch_I_2Min <- 10
tolBisecDeterm_I_2Min <- 1/10^3
tolBisecDetermLevelConst <- 1/10^6


A_F <- function(z1, c) { # Fishers product test
  A_def <- c/(1-pnorm(z1))
  return(min(A_def, 0.5))
}

# conditional error function for design with separate studies:
A_C <- function(z1, c){
  return(alpha)
}


# Example 1: memodio------------------------------------------------------------
xi_1 = 1.545918/1
alpha_c <- 0.05
alpha <- 0.025
rel_I1_min_xi_1<- determine_I1min(xi_1, beta=0.2, alpha_c)
rel_I1_max_xi_1 <- determine_I1max(xi_1, alpha)
relativeI1vec_1 <- seq(rel_I1_min_xi_1+1/10^6, rel_I1_max_xi_1 -1/10^6,
                       by = 0.001) 
numI1 <- length(relativeI1vec_1)

# set important boundaries: 
z_er1 <- Inf # no early rejection
z_f1 <- pmax(sqrt(relativeI1vec_1)*eta_f/xi_1, qnorm(1-alpha_c))
z_s1 <- z_f1  # assumed binding futility stop

# determine level constants:
c_AF1 <- determineLevelConst(A_F, relativeI1vec_1, z_s1, z_er1,
                            searchIntervalDeterLevelConst= c(0.0038,1),
                            tolBisecDetermLevelConst,
                            stabConstIntegral)
c_AC1 <- rep(alpha,numI1)

# set further input parameter:
p_1_1 <- 0.8 # aimed overall sucess probability

#################### determination of I_2,min/I_delta: ########################

relativeI_2MinFisher1_cp80 <- determinationI_2Min(relativeI1vec_1, p_1_1, xi_1, A_F,
                                            c_AF1, z_er1, z_f1, z_s1, cp=0.2,
                                            uppBoundSearch_I_2Min,
                                            tolBisecDeterm_I_2Min,
                                            stabConstIntegral)
relativeI_2MinFisher1_cp85 <- determinationI_2Min(relativeI1vec_1, p_1_1, xi_1, A_F,
                                                  c_AF1, z_er1, z_f1, z_s1, cp=0.15,
                                                  uppBoundSearch_I_2Min,
                                                  tolBisecDeterm_I_2Min,
                                                  stabConstIntegral)
relativeI_2MinFisher1_cp90 <- determinationI_2Min(relativeI1vec_1, p_1_1, xi_1, A_F,
                                                  c_AF1, z_er1, z_f1, z_s1, cp=0.1,
                                                  uppBoundSearch_I_2Min,
                                                  tolBisecDeterm_I_2Min,
                                                  stabConstIntegral)
relativeI_2MinClas1_cp80 <- determinationI_2Min(relativeI1vec_1, p_1_1, xi_1, A_C,
                                           c_AC1, z_er1, z_f1, z_s1, cp=0.2,
                                           uppBoundSearch_I_2Min, 
                                           tolBisecDeterm_I_2Min,
                                           stabConstIntegral)
relativeI_2MinClas1_cp85 <- determinationI_2Min(relativeI1vec_1, p_1_1, xi_1, A_C,
                                                c_AC1, z_er1, z_f1, z_s1, cp=0.15,
                                                uppBoundSearch_I_2Min, 
                                                tolBisecDeterm_I_2Min,
                                                stabConstIntegral)
relativeI_2MinClas1_cp90 <- determinationI_2Min(relativeI1vec_1, p_1_1, xi_1, A_C,
                                                c_AC1, z_er1, z_f1, z_s1, cp=0.1,
                                                uppBoundSearch_I_2Min, 
                                                tolBisecDeterm_I_2Min,
                                                stabConstIntegral)

# caclulate minimal trial study sample size for adaptive design 
# (and separate design) for overall power of 80% at delta_apriori: 
n_1 <- 70
n_80PercentPower_at_delta_ariori_singleStage <- 74 # calculated with Example 1 file

index <- max(which(relativeI1vec_1<=n_1/n_80PercentPower_at_delta_ariori_singleStage))
n_2_min_Fisher_cp80 <- ceiling(relativeI_2MinFisher1_cp80[index]*n_80PercentPower_at_delta_ariori_singleStage)
n_2_min_Fisher_cp85 <- ceiling(relativeI_2MinFisher1_cp85[index]*n_80PercentPower_at_delta_ariori_singleStage)
n_2_min_Fisher_cp90 <- ceiling(relativeI_2MinFisher1_cp90[index]*n_80PercentPower_at_delta_ariori_singleStage)

c_AF <- c_AF1[index] # needed for conditional error function


n_2_min_separate_cp80 <- ceiling(relativeI_2MinClas1_cp80[index]*n_80PercentPower_at_delta_ariori_singleStage)
n_2_min_separate_cp85 <- ceiling(relativeI_2MinClas1_cp85[index]*n_80PercentPower_at_delta_ariori_singleStage)
n_2_min_separate_cp90 <- ceiling(relativeI_2MinClas1_cp90[index]*n_80PercentPower_at_delta_ariori_singleStage)

# Example 2: Orthopy -----------------------------------------------------------
# set important boundaries: 
z_er2 <- Inf # no early rejection
z_s2 <- -Inf # no early binding futility bound
relativeI1vec_2 <- 1 # can be set arbitrary because of no binding futility stop

# determine level constants:
c_AF2 <- determineLevelConst(A_F, relativeI1vec_2, z_s2, z_er2,
                             searchIntervalDeterLevelConst= c(0.0038,1),
                             tolBisecDetermLevelConst,
                             stabConstIntegral)
