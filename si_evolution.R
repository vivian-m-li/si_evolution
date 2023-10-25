# Social information evolution model
rm(list = ls())
#ptm <- proc.time()
require(ggplot2)
require(gridExtra)
source("analyze.R")

# # assumptions (parameters in function below):
# Ni = 500, # number of individuals
# tf = 30, #300 # time steps for a given generation
# e_gain = 1, # energy units gained each time step in which the ith individual does not flee
# coef_false = 0.20, # coefficient that determines the prob of a false alarm (a smaller value than f_pred)
# # we start with X individuals, whose evolvable traits of interest ('jumpiness', social faith, density dependence in social faith) are determined by random draws from Gaussain distributions
# # the jumpiness trait determines an individual's probability of: randomly fleeing (NOT due to a real threat, i.e., a false alarm) AND an individual's probability of detecting an attacking predator. As a first pass, let's simply make the false alarm probability much lower than the response to predators probability. Perhaps 10-20% of the response to predators probability:
# ######################### CAN CUT THIS (repeated in code below)
# coef_false = 0.20,
# f_pred = runif(Ni, min=0, max=1), # using uniform distribution to constrain sampling between 0 and 1
# # the 'social faith' parm is an intercept term that represents the probability that a given individual will flee in a time step, given that another individual in the group fled. We can make the ordering of individuals random, so that we don't make it consistently less or more likely that certain individuals will flee
# s_faith = runif(Ni, min=0, max=0.5),
# s_dd = runif(Ni, min=-2, max=2),
# #########################
# fit = matrix(NA, ncol=tf,nrow=Ni),
# fit[ , 1] = 1, # setting fitness of all individuals at the start at 1
# maxf = 100 # number of generations to run the model through


evo.fun1 <- function(Ni = 100, tf = 30, e_gain = 1, coef_false = 0.20, maxf = 100, prob_pred = 0.2, max_group_size = 25) {
  
  f_pred = runif(Ni, min = 0, max = 1)
  s_faith = runif(Ni, min = 0, max = 0.5)
  s_dd = runif(Ni, min = -2, max = 2)
  fit = matrix(NA, ncol = tf, nrow = Ni)
  fit[ , 1] = 1
  
  # f_false = coef_false*f_pred
  flights <- list() # collecting all flight info for a generation
  attacks_all <- list()
  eaten_detect_all <- list()
  eaten_nodetect_all <- list()
  fit_traits_gen0 <- list() # list to store final fitness and traits of all individuals (alive and dead) in each generation
  fit_traits_gen <- list() # list to store final fitness and traits of SURVIVING individuals in each generation
  trait_mean <- matrix(NA, nrow = maxf, ncol = 3) # to store mean values of survivor traits for each generation
  trait_sd <- matrix(NA, nrow = maxf, ncol = 3) # to store sd of survivor traits for each generation
  f_all <- list() # this list stores the new set of traits for each generation (draw from parents in prev gen). This will end up with maxf-1 items (since there was not one for the first gen).
  flights_master <- list()
  eaten_detect_master <- list()
  eaten_nodetect_master <- list()
  
  for (f in 1:maxf){
    if (f == 1){ #first gen traits are randomly determined
      f_pred <- runif(Ni, min = 0, max = 1) # using uniform distribution to constrain sampling bewteen 0 and 1
      f_false <- coef_false*f_pred
      s_faith <- runif(Ni, min = 0, max = 0.5)
      s_dd <- runif(Ni, min = -2, max = 2)
    } else{ # from the 2nd gen onward, traits are inherited based on parent fitness
      f_pred <- f_all[[f]][ , 1]
      f_false <- coef_false*f_pred
      s_faith <- f_all[[f]][ , 2]
      s_dd <- f_all[[f]][ , 3]
    }
    
    # For each time step, the population gets reassembled into groups based on a uniform distribution of group sizes. Then, each group is potentially subjected to a predator attack (based on a background predation level, probability set to 0.2 by default below). 
    
    for (t in 1:(tf - 1)){ # running through all the time steps of a generation
      group_sizes = seq(1:max_group_size)
      groups <- list()
      j <- 1
      while (sum(unlist(groups)) < Ni){
        j <- j + 1
        groups[j] = sample(group_sizes, size = 1, replace = FALSE)
      }
      groups = unlist(groups)
      remainder = sum(groups) - Ni
      groups[length(groups)] <- groups[length(groups)] - remainder
      group_vec <- list()
      for (k in 1:length(groups)){
        vec <- rep(k, groups[k])
        group_vec[[k]] <- vec
      }	
      group_vec <- unlist(group_vec) # expanding group to list of individuals with group assignments
      seq0 <- sample(seq(1:Ni), size = Ni, replace = FALSE) # setting a random sequence of individuals that can each respond to the predator or to non threats (false alarm))
      groups_df <- data.frame(individual = seq0, groupID = group_vec)
      flights0 <- matrix(NA, nrow=length(unique(group_vec)), ncol = 4) # collecting flight info for each group in a time step
      eaten_detect0 <- matrix(NA, nrow=length(unique(group_vec)), ncol = 2)
      eaten_nodetect0 <- matrix(NA, nrow=length(unique(group_vec)), ncol = 2)
      attacks_vec <- rep(NA, length(unique(group_vec)))
      for (g in 1:length(unique(group_vec))){
        pred <- rbinom(1, 1, prob_pred) #0.1 # determining if a predator attacks the group for the given time step
        attacks_vec[g] <- pred
        prev_flee <- 0 # setting the initial number of observable neighbor flights for a given group in a given time step
        subgroup <- subset(groups_df, groups_df$groupID == g)
        ddensity <- sum(as.numeric((fit[subgroup$individual, t] > 0))) # determining the density of LIVING individuals in the group
        prev_detect <- 0 # setting the initial number of observable neighbor flights from a predator for a given group in a given time step in which a predator attacks
        eaten_detect_vec <- rep(0, nrow(subgroup))
        eaten_nodetect_vec <- rep(0, nrow(subgroup))
        for (i in 1:length(subgroup$individual)){
          ii <- subgroup$individual[i] # determining the order of prey decisions randomly for each time step
          #density <- length(subgroup$individual)
          eaten_detect <- 0; eaten_nodetect <- 0 # setting defaults/placeholders
          if (fit[ii, t] == 0){ # if prey were EATEN in a previous time step, their fitness remains at zero, otherwise their fitness is >0 (they start out with a fitness value of 1)
            fit[ii, t + 1] <- 0
          } else { # if focal was NOT eaten in a previous time step
            if (prev_flee > 0) { # if one or more neighbors has already fled in this time step
              p_flee_s <- prev_flee*s_faith[ii] + (ddensity - prev_flee)*s_dd[ii] # I think this density should be the density of REMAINING individuals
              if (p_flee_s > 1){
                p_flee_s <- 1
              }
              if (p_flee_s < 0){
                p_flee_s <- 0
              }
              flee0 <- rbinom(1, 1, f_false[ii]) + rbinom(1, 1, p_flee_s) # flight of focal depends on both random flight and flight of neighbors
              flee <- as.numeric(flee0 >= 1)
            }else { # if no neighbors have fled yet in this time step
              flee <- rbinom(1, 1, f_false[ii])
            }
            prev_flee <- prev_flee + flee # updating flights from false alarms
            if (flee == 1 | pred == 1){ # if focal flees OR the predator attacks, the focal does not gain energy (eat) in this time step
              fit[ii, t + 1] <- fit[ii, t] 
            }
            if (pred==1 & flee==0){ # if a predator attacks AND the focal did NOT flee on this time step...
              p_detect_s <- sum(prev_detect)*s_faith[ii] + (ddensity-prev_flee-prev_detect)*s_dd[ii]
              if (p_detect_s > 1) {
                p_detect_s <- 1
              }
              if (p_detect_s < 0) {
                p_detect_s <- 0
              }
              detect0 <- rbinom(1, 1, f_pred[ii]) + rbinom(1, 1, p_detect_s) # move to prob!?!!!!!!!!!!!!!
              detect <- as.numeric(detect0 >= 1)
              prev_detect <- prev_detect + detect # updating flights from predator
              if (detect == 1){ # if the focal detects the predator
                peaten_detect <- 1/((length(subgroup$individual) - prev_flee) + 10) # probability of life or death (the +10 is a bonus for early detection)
                eaten_detect <- rbinom(1, 1, peaten_detect)
                eaten_detect_vec[i] <- eaten_detect
              } else { # if the focal does not detect the predator
                peaten_nodetect <- 1/((length(subgroup$individual) - prev_flee)) # LOWER prob life/death (since pred was not detected)
                eaten_nodetect <- rbinom(1, 1, peaten_nodetect)
                eaten_nodetect_vec[i] <- eaten_nodetect
                # add up each eaten for each group and store in a list at each time point
              }
              if (eaten_detect == 1 | eaten_nodetect == 1){
                fit[ii, t + 1] <- 0 # if eaten, the focal's fitness goes to 0				
              }
            }
            if (pred == 0 & flee == 0){ # if predator does NOT attack AND focal does NOT flee (from a false alarm)
              fit[ii, t + 1] <- fit[ii, t] + e_gain
            }
          }
          #print(c(t,g,i))
        } # end of i loop
        # make a matrix that is col=var, row=# groups, and all this matrix to a list each time step
        flights0[g, ] <- c(t, ddensity, prev_flee, prev_detect) # recording all flights, from false alarms AND predators, for a given group and time step (groups change each time step)
        eaten_detect0[g, ] <- c(t, sum(eaten_detect_vec))
        eaten_nodetect0[g, ] <- c(t, sum(eaten_nodetect_vec))
      } # end of g loop
      attacks_all[[t]] <- attacks_vec
      eaten_detect_all[[t]] <- eaten_detect0
      eaten_nodetect_all[[t]] <- eaten_nodetect0
      flights[[t]] <- flights0
    } # end of t loop
    flights_all <- do.call(rbind, flights)
    flights_master[[f]] <- flights_all
    eaten_detect_mat <- do.call(rbind, eaten_detect_all)
    eaten_detect_master[[f]] <- eaten_detect_mat
    eaten_nodetect_mat <- do.call(rbind, eaten_nodetect_all)
    eaten_nodetect_master[[f]] <- eaten_nodetect_mat	
    
    survive0 <- cbind(fit[,tf], f_pred, s_faith, s_dd) # combining final fitness and traits for this generation
    fit_traits_gen0[[f]] <- survive0
    survive <- survive0[survive0[,1]>0,] # subsetting only the survivors (those not eaten by predators)
    fit_traits_gen[[f]] <- survive
    trait_mean[f,] <- c(mean(survive[,2]), mean(survive[,3]), mean(survive[,4]))
    trait_sd[f,] <- c(sd(survive[,2]), sd(survive[,3]), sd(survive[,4]))
    survive[,1] <- survive[,1]/sum(survive[,1]) # calculating proportional fitness (fitness of population sums to zero)
    surv_df <- data.frame(survive) 
    surv_df$index <- seq(1:length(surv_df$V1)) # creating a column of row numbers to sample from for the traits of next gen
    f_index <- sample(surv_df$index, Ni, prob=surv_df$V1, replace=TRUE) # sampling next gen traits weighted by parent fitness
    f_all[[f + 1]] <- survive[f_index,2:4] # traits of all Ni individuals for next gen	
  } # end of f loop
  
  # to select only survivors from the complete 'fit_traits_gen' list (which contains fitnesses for all individuals in every generation), I'd modify this line:
  #dev.new()
  #par(mfrow=c(1,3))
  traits <- c("jumpiness", "sociality", "density dependence in sociality")
  gg_trait_ls <- list()
  out_ls <- list()
  mean_ff <- round(mean(fit_traits_gen[[maxf]][,1]), 2) # mean final fitness for last generation
  sd_ff <- round(sd(fit_traits_gen[[maxf]][,1]), 2) # sd of final fitness for last generation
  for (p in 1:length(traits)){
    dff <- data.frame(generation = seq(1:nrow(trait_mean)), trait_mean = trait_mean[ , p], lb = trait_mean[ , p] - trait_sd[ , p], ub = trait_mean[ , p] + trait_sd[ , p])
    out_ls[[p]] <- dff
    gg_trait_ls[[p]] <- ggplot(dff, aes(generation)) + 
      geom_line(aes(y=trait_mean), colour="black") + 
      geom_ribbon(aes(ymin=lb, ymax=ub), alpha=0.2) +
      labs(y = traits[p])
  }
  #proc.time() - ptm
  dev.new()
  grid.arrange(gg_trait_ls[[1]] + ggtitle(paste("Fitness=",mean_ff,"+/-",sd_ff,"SD")), gg_trait_ls[[2]], gg_trait_ls[[3]], ncol = 1)
  
  # multiplot(gg_trait_ls[[1]] + ggtitle(paste("Fitness=",mean_ff,"+/-",sd_ff,"SD")),gg_trait_ls[[2]],gg_trait_ls[[3]])
  
  return(out_ls)
}
