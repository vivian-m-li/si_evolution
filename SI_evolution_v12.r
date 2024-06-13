# after Hein meeting, replace second term with alpha*e^(beta*N)  

# NOTE: it took 6 hours to run 4 pred levels and 10 runs for 2 models (there was a 3rd attempted model run, which ended in error, but I'm not sure how long that took. So, it's possible that the time to run the two models only was significantly less than 6 hours...

# Social information evolution model
rm(list=ls())

require(ggplot2)
require(grid)

# multiplot for multi-panel ggplots
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# creating cluster 
require(doSNOW)
require(foreach)

cl<-makeCluster(4) # number of CPU cores
registerDoSNOW(cl)

ptm <- proc.time()

s1_ON=1; s2_ON=1; s_dd_ON=1
flee_fun = function(s1_ON=1, s2_ON=1, s_dd_ON=1) {
	p_pred <- seq(from=0.02, to=0.2, length.out = 4) #10
	master_mat_N1 <- list()
	master_mat_N2 <- list()
	false_flights_dlist <- list()
	pred_flights_dlist <- list()

	dList = foreach (d=1:length(p_pred)) %dopar% { # d is for 'danger'!
		runs = 10
		false_flights_rlist <- list()
		pred_flights_rlist <- list()	
		for (r in 1:runs){
			# assumptions:
			N1i <- 100 # number of individuals from sp 1
			N2i <- 100
			tf <- 250 #250 time steps for a given generation
			maxg <- 25 # max group size
			e_gain <- 1 # energy units gained each time step in which the ith individual does not flee
			coef_false <- 0.20 # coefficient that determines the prob of a false alarm (a smaller value than f_pred)

			# we start with X individuals, whose evolvable traits of interest ('jumpiness', social faith, density dependence in social faith) are determined by random draws from Gaussain distributions
			# the jumpiness trait determines an individual's probability of: randomly fleeing (NOT due to a real threat, i.e., a false alarm) AND an individual's probability of detecting an attacking predator. As a first pass, let's simply make the false alarm probability much lower than the response to predators probability. Perhaps 10-20% of the response to predators probability:
			fit_N1 <- matrix(NA, ncol=tf,nrow=N1i)
			fit_N2 <- fit_N1
			fit_N1[,1] <- 1 # setting fitness of all individuals at the start at 1
			fit_N2[,1] <- 1
			flights <- list() # collecting all flight info for a generation
			attacks_all <- list()
			eaten_detect_all <- list()
			eaten_nodetect_all <- list()
			maxf <- 200 # number of generations to run the model through
			fit_traits_gen0_N1 <- list() # list to store final fitness and traits of all individuals (alive and dead) in each generation, for sp. 1
			fit_traits_gen0_N2 <- list() # list to store final fitness and traits of all individuals (alive and dead) in each generation, for sp. 1
			fit_traits_gen_N1 <- list() # list to store final fitness and traits of SURVIVING individuals in each generation
			fit_traits_gen_N2 <- list()
			trait_mean_N1 <- matrix(NA, nrow=maxf, ncol=4) # to store mean values of survivor traits for each generation
			trait_mean_N2 <- matrix(NA, nrow=maxf, ncol=4) 
			trait_sd_N1 <- matrix(NA, nrow=maxf, ncol=4) # to store sd of survivor traits for each generation
			trait_sd_N2 <- matrix(NA, nrow=maxf, ncol=4) # to store sd of survivor traits for each generation
			f_all_N1 <- list() # this list stores the new set of traits for each generation (draw from parents in prev gen). This will end up with maxf-1 items (since there was not one for the first gen).
			f_all_N2 <- list() 
			flights_master_f <- list()
			eaten_detect_master <- list()
			eaten_nodetect_master <- list()
			false_flights_flist <- list()
			pred_flights_flist <- list()

			for (f in 1:maxf){
				if (f==1){ #first gen traits are randomly determined
					f_pred_N1 <- runif(N1i, min=0, max=0.5)
					f_false_N1 <- coef_false*f_pred_N1
					f_pred_N2 <- runif(N2i, min=0, max=0.5)
					f_false_N2 <- coef_false*f_pred_N2
					if (s1_ON != 0) {
						s1_N1 <- runif(N1i, min=0, max=0.2)
						s1_N2 <- runif(N2i, min=0, max=0.2)
					} else {
						s1_N1 <- rep(0,100)
						s1_N2 <- rep(0,100)
					}
					if (s2_ON != 0) {
						s2_N1 <- runif(N1i, min=0, max=0.2)
						s2_N2 <- runif(N2i, min=0, max=0.2)
					} else {
						s2_N1 <- rep(0,100)
						s2_N2 <- rep(0,100)
					}
					if (s_dd_ON != 0) {
						s_dd_N1 <- runif(N1i, min=0, max=0.2)
						s_dd_N2 <- runif(N2i, min=0, max=0.2)
					} else {
						s_dd_N1 <- rep(0,100)
						s_dd_N2 <- rep(0,100)
					}
				} else{ # from the 2nd gen onward, traits are inherited based on parent fitness
					f_pred_N1 <- f_all_N1[[f]][,1]
					f_false_N1 <- coef_false*f_pred_N1
					f_pred_N2 <- f_all_N2[[f]][,1]
					f_false_N2 <- coef_false*f_pred_N2
					s1_N1 <- f_all_N1[[f]][,2]
					s1_N2 <- f_all_N2[[f]][,2]
					s2_N1 <- rep(0,100) #f_all_N1[[f]][,3]
					s2_N2 <- rep(0,100) #f_all_N2[[f]][,3]
					s_dd_N1 <- f_all_N1[[f]][,4]
					s_dd_N2 <- f_all_N1[[f]][,4]
				}
				false_flights_tlist <- list()
				pred_flights_tlist <- list()
				for (t in 1:(tf-1)){ # running through all the time steps of a generation
					group_sizes = seq(1:maxg)
					groups <- list()
					j = 1
					while (sum(unlist(groups))<(N1i+N2i)){ # sampling groups until both species run out
						j <- j+1
						groups[j] = sample(group_sizes, size=1, replace=TRUE)
					}
					groups = unlist(groups)
					remainder = sum(groups)-(N1i+N2i)
					groups[length(groups)] <- groups[length(groups)]-remainder # subtracting remainder from last group
					group_vec <- list()
					for (k in 1:length(groups)){ # creating a vector that will assign each individual (arranged in a random order, below) to a group number
						vec <- rep(k,groups[k])
						group_vec[[k]] <- vec
					}	
					group_vec <- unlist(group_vec)
					seq0 <- sample(seq(1:(N1i+N2i)), size=N1i+N2i, replace=FALSE) # setting a random sequence of individuals that can each respond to the predator or to non threats (false alarm))
					groups_df <- data.frame(individual=seq0, groupID=group_vec)
					flights0 <- matrix(NA, nrow=length(unique(group_vec)), ncol=4) # collecting flight info for each group in a time step
					eaten_detect0 <- matrix(NA, nrow=length(unique(group_vec)), ncol=2)
					eaten_nodetect0 <- matrix(NA, nrow=length(unique(group_vec)), ncol=2)
					attacks_vec <- rep(NA, length(unique(group_vec)))
					for (g in 1:length(unique(group_vec))){ # running through each group
						pred <- rbinom(1, 1, p_pred[[d]]) #0.1 # determining if a predator attacks the group for the given time step
						attacks_vec[g] <- pred
						prev_flee <- 0 # setting the initial number of observable neighbor flights for a given group in a given time step
						subgroup0 <- subset(groups_df, groups_df$groupID==g)
						subgroup1 <- subset(subgroup0, subgroup0$individual<101)
						subgroup2 <- subset(subgroup0, subgroup0$individual>100)
						subgroup2$individual <- subgroup2$individual - 100 
						density_N1 <- sum(as.numeric((fit_N1[subgroup1$individual,t]>0)))
						density_N2 <- sum(as.numeric((fit_N2[subgroup2$individual,t]>0)))
						density <- density_N1 + density_N2 # determining density of LIVING individuals in group
						prev_detect <- 0 # setting the initial number of observable neighbor flights from a predator for a given group in a given time step in which a predator attacks
						eaten_detect_vec <- rep(0,length(subgroup0$individual))
						eaten_nodetect_vec <- rep(0,length(subgroup0$individual))
						for (i in 1:length(subgroup0$individual)){ # running through each individual in gth group
							if (subgroup0$individual[i]<101) { # if the individual is from sp 1
								ii <- subgroup0$individual[i]
								sp <- 1
							} else {
								ii0 <- subgroup0$individual[i]
								ii <- ii0 - 100
								sp <- 2
							}
							eaten_detect <- 0; eaten_nodetect <- 0 # setting defaults/placeholders
							if (sp == 1){
								if (fit_N1[ii,t]==0){ # if prey were EATEN in a previous time step, their fitness remains at zero, otherwise their fitness is >0 (they start out with a fitness value of 1)
									fit_N1[ii,t+1] <- 0
								} else { # if focal was NOT eaten in a previous time step
									p_flee <- f_false_N1[ii]*exp(s_dd_N1[ii]*(density-prev_flee)) + s1_N1[ii]*prev_flee*exp(s2_N1[ii]*(density-prev_flee))  # density of REMAINING individuals
									if (p_flee>1){
										p_flee=1
									}
									if (p_flee<0){
										p_flee=0
									}
									flee <- rbinom(1,1,p_flee) # flight of focal depends on both random flight and flight of neighbors
									prev_flee <- prev_flee + flee # updating flights from false alarms
									if (flee==1 | pred==1){ # if focal flees OR the predator attacks, the focal does not gain energy (eat) in this time step
										fit_N1[ii,t+1] <- fit_N1[ii,t] 
									}
									if (pred==1 & flee==0){ # if a predator attacks AND the focal did NOT already flee (via false alarm) on this time step...
										p_detect <- f_pred_N1[ii]*exp(s_dd_N1[ii]*(density-prev_flee-prev_detect))+ s1_N1[ii]*prev_detect*exp(s2_N1[ii]*(density-prev_flee-prev_detect))
										if (p_detect>1) {
											p_detect=1
										}
										if (p_detect<0) {
											p_detect=0
										}
										detect <- rbinom(1,1,p_detect)
										prev_detect <- prev_detect + detect # updating flights from predator
										if (detect==1){ # if the focal detects the predator
											peaten_detect <- 1/((length(subgroup0$individual)-prev_flee) + 10) # probability of life or death (the +10 is a bonus for early detection)
											eaten_detect <- rbinom(1,1,peaten_detect)
											eaten_detect_vec[i] <- eaten_detect
										} else { # if the focal does not detect the predator
											peaten_nodetect <- 1/((length(subgroup0$individual)-prev_flee-prev_detect)) # LOWER prob life/death (since pred was not detected) # NOTE**: could also adjust density further here: subtracting individuals who fled from the predator
											eaten_nodetect <- rbinom(1,1,peaten_nodetect)
											eaten_nodetect_vec[i] <- eaten_nodetect
											# add up each eaten for each group and store in a list at each time point
										}
										if (eaten_detect==1 | eaten_nodetect==1){
											fit_N1[ii,t+1] <- 0 # if eaten, the focal's fitness goes to 0				
										}
									}
									if (pred==0 & flee==0){ # if predator does NOT attack AND focal does NOT flee (from a false alarm)
										fit_N1[ii,t+1] <- fit_N1[ii,t] + e_gain
									}
								}
							} else {
								if (fit_N2[ii,t]==0){ # if prey were EATEN in a previous time step, their fitness remains at zero, otherwise their fitness is >0 (they start out with a fitness value of 1)
									fit_N2[ii,t+1] <- 0
								} else { # if focal was NOT eaten in a previous time step
									p_flee <- f_false_N1[ii]*exp(s_dd_N1[ii]*(density-prev_flee)) + s1_N1[ii]*prev_flee*exp(s2_N1[ii]*(density-prev_flee))  # density of REMAINING individuals
									if (p_flee>1){
										p_flee=1
									}
									if (p_flee<0){
										p_flee=0
									}
									flee <- rbinom(1,1,p_flee) # flight of focal depends on both random flight and flight of neighbors
									prev_flee <- prev_flee + flee # updating flights from false alarms
									if (flee==1 | pred==1){ # if focal flees OR the predator attacks, the focal does not gain energy (eat) in this time step
										fit_N2[ii,t+1] <- fit_N2[ii,t] 
									}
									if (pred==1 & flee==0){ # if a predator attacks AND the focal did NOT flee on this time step...
										p_detect <- f_pred_N2[ii]*exp(s_dd_N2[ii]*(density-prev_flee-prev_detect))+ s1_N2[ii]*prev_detect*exp(s2_N2[ii]*(density-prev_flee-prev_detect))									
										if (p_detect>1) {
											p_detect=1
										}
										if (p_detect<0) {
											p_detect=0
										}
										detect <- rbinom(1,1,p_detect)
										prev_detect <- prev_detect + detect # updating flights from predator
										if (detect==1){ # if the focal detects the predator
											peaten_detect <- 0 #1/((length(subgroup0$individual)-prev_flee) + 10) # probability of life or death (the +10 is a bonus for early detection)
											eaten_detect <- rbinom(1,1,peaten_detect)
											eaten_detect_vec[i] <- eaten_detect
										} else { # if the focal does not detect the predator
											peaten_nodetect <- 1/((length(subgroup0$individual)-prev_flee-prev_detect)) # LOWER prob life/death (since pred was not detected) # NOTE**: could also adjust density further here: subtracting individuals who fled from the predator
											eaten_nodetect <- rbinom(1,1,peaten_nodetect)
											eaten_nodetect_vec[i] <- eaten_nodetect
											# add up each eaten for each group and store in a list at each time point
										}
										if (eaten_detect==1 | eaten_nodetect==1){
											fit_N2[ii,t+1] <- 0 # if eaten, the focal's fitness goes to 0				
										}
									}
									if (pred==0 & flee==0){ # if predator does NOT attack AND focal does NOT flee (from a false alarm)
										fit_N2[ii,t+1] <- fit_N2[ii,t] + e_gain
									}
								}
							}
						#print(c(t,g,i))
						} # end of i loop
						# make a matrix that is col=var, row=# groups, and all this matrix to a list each time step
						flights0[g,] <- c(t, density, prev_flee, prev_detect) # recording all flights, from false alarms AND predators, for a given group and time step (groups change each time step)
						eaten_detect0[g,] <- c(t, sum(eaten_detect_vec))
						eaten_nodetect0[g,] <- c(t, sum(eaten_nodetect_vec))
					} # end of g loop
					
					flights_eaten_df <- data.frame(p_pred = rep(p_pred[d], length(flights0[,1])), run = r, gen = rep(f, length(flights0[,1])), time= flights0[,1], dens_false= flights0[,2], false_flee= flights0[,3], attack= attacks_vec, dens_pred=flights0[,2]-flights0[,3], pred_flee= flights0[,4], eaten_detect= eaten_detect0[,2], eaten_nodetect= eaten_nodetect0[,2])
					
					false_flights <- subset(flights_eaten_df, flights_eaten_df$false_flee>0)
					false_flights_tlist[[t]] <- false_flights
					pred_flights <- subset(flights_eaten_df, flights_eaten_df$pred_flee>0)
					pred_flights_tlist[[t]] <- pred_flights
					####

				} # end of t loop (completed a time step for the series of random groups of individuals
						
				false_flights_flist[[f]] <- do.call(rbind, false_flights_tlist)
				pred_flights_flist[[f]] <- do.call(rbind, pred_flights_tlist)
				
				# sp 1:
				survive_N1_0 <- cbind(fit_N1[,tf], f_pred_N1, s1_N1, s2_N1, s_dd_N1) # combining final fitness and traits for this generation
				fit_traits_gen0_N1[[f]] <- survive_N1_0
				survive_N1 <- survive_N1_0[survive_N1_0[,1]>0,] # subsetting only the survivors (those not eaten by predators)
				fit_traits_gen_N1[[f]] <- survive_N1
				trait_mean_N1[f,] <- c(mean(survive_N1[,2]), mean(survive_N1[,3]), mean(survive_N1[,4]), mean(survive_N1[,5]))
				trait_sd_N1[f,] <- c(sd(survive_N1[,2]), sd(survive_N1[,3]), sd(survive_N1[,4]), sd(survive_N1[,5]))
				survive_N1[,1] <- survive_N1[,1]/sum(survive_N1[,1]) # calculating proportional fitness (fitness of population sums to zero)
				surv_N1_df <- data.frame(survive_N1) 
				surv_N1_df$index <- seq(1:length(surv_N1_df$V1)) # creating a column of row numbers to sample from for the traits of next gen
				f_index_N1 <- sample(surv_N1_df$index, N1i, prob=surv_N1_df$V1, replace=TRUE) # sampling next gen traits weighted by parent fitness
				f_all_N1_0 <- survive_N1[f_index_N1,2:5] # traits of all Ni individuals for next gen	
				f_mut_N1_0 <- cbind(rnorm(100, mean=0, sd=0.04*abs(mean(f_all_N1_0[,1]))), rnorm(100, mean=0, sd=0.04*abs(mean(f_all_N1_0[,2]))), rnorm(100, mean=0, sd=0.04*abs(mean(f_all_N1_0[,3]))), rnorm(100, mean=0, sd=0.04*abs(mean(f_all_N1_0[,4]))))	# introducing mutation
				if (s1_ON == 0) { 
					f_mut_N1_0[,2] = rep(0, length(f_mut_N1_0[,2])) 
				}
				if (s2_ON == 0) { 
					f_mut_N1_0[,3] = rep(0, length(f_mut_N1_0[,3])) 
				}
				if (s_dd_ON == 0) { 
					f_mut_N1_0[,4] = rep(0, length(f_mut_N1_0[,4])) 
				}
				f_mut_N1 <- f_all_N1_0 + f_mut_N1_0
				f_mut_N1[,1] <- replace(f_mut_N1[,1], f_mut_N1[,1] < 0, 0)
				f_mut_N1[,2] <- replace(f_mut_N1[,2], f_mut_N1[,2] < 0, 0)
				f_all_N1[[f+1]] <- f_mut_N1
				
				# sp 2 (repeating above)
				survive_N2_0 <- cbind(fit_N2[,tf], f_pred_N2, s1_N2, s2_N2, s_dd_N2) # combining final fitness and traits for this generation
				fit_traits_gen0_N2[[f]] <- survive_N2_0
				survive_N2 <- survive_N2_0[survive_N2_0[,1]>0,] # subsetting only the survivors (those not eaten by predators)
				fit_traits_gen_N2[[f]] <- survive_N2
				trait_mean_N2[f,] <- c(mean(survive_N2[,2]), mean(survive_N2[,3]), mean(survive_N2[,4]), mean(survive_N2[,5]))
				trait_sd_N2[f,] <- c(sd(survive_N2[,2]), sd(survive_N2[,3]), sd(survive_N2[,4]), sd(survive_N2[,5]))
				survive_N2[,1] <- survive_N2[,1]/sum(survive_N2[,1]) # calculating proportional fitness (fitness of population sums to zero)
				surv_N2_df <- data.frame(survive_N2) 
				surv_N2_df$index <- seq(1:length(surv_N2_df$V1)) # creating a column of row numbers to sample from for the traits of next gen
				f_index_N2 <- sample(surv_N2_df$index, N2i, prob=surv_N2_df$V1, replace=TRUE) # sampling next gen traits weighted by parent fitness
				f_all_N2_0 <- survive_N2[f_index_N2,2:5] # traits of all Ni individuals for next gen	
				f_mut_N2_0 <- cbind(rnorm(100, mean=0, sd=0.04*abs(mean(f_all_N2_0[,1]))), rnorm(100, mean=0, sd=0.04*abs(mean(f_all_N2_0[,2]))), rnorm(100, mean=0, sd=0.04*abs(mean(f_all_N2_0[,3]))), rnorm(100, mean=0, sd=0.04*abs(mean(f_all_N2_0[,4])))) # introducing mutation
				if (s1_ON == 0) { 
					f_mut_N1_0[,2] = rep(0, length(f_mut_N1_0[,2])) 
				}
				if (s2_ON == 0) { 
					f_mut_N1_0[,3] = rep(0, length(f_mut_N1_0[,3])) 
				}
				if (s_dd_ON == 0) { 
					f_mut_N1_0[,4] = rep(0, length(f_mut_N1_0[,4])) 
				}
				f_mut_N2 <- f_all_N2_0 + f_mut_N2_0
				f_mut_N2[,1] <- replace(f_mut_N2[,1], f_mut_N2[,1] < 0, 0)
				f_mut_N2[,2] <- replace(f_mut_N2[,2], f_mut_N2[,2] < 0, 0)
				f_all_N2[[f+1]] <- f_mut_N2
			
			} # end of f loop
			
			#flights_master_r[[r]] <- flights_master_f
			false_flights_rlist[[r]] <- do.call(rbind, false_flights_flist)
			pred_flights_rlist[[r]] <- do.call(rbind, pred_flights_flist)
			
			# sp 1
			mean_fit_maxf_N1 <- mean(fit_traits_gen_N1[[maxf]][,1]) # mean final fitness for last generation
			sd_fit_maxf_N1 <- sd(fit_traits_gen_N1[[maxf]][,1]) # sd of final fitness for last generation
			mean_f_pred_maxf_N1 <- mean(fit_traits_gen_N1[[maxf]][,2]) # mean and sd of each trait for last generation...
			sd_f_pred_maxf_N1 <- sd(fit_traits_gen_N1[[maxf]][,2]) 
			mean_s1_maxf_N1 <- mean(fit_traits_gen_N1[[maxf]][,3]) 
			sd_s1_maxf_N1 <- sd(fit_traits_gen_N1[[maxf]][,3])
			mean_s2_maxf_N1 <- mean(fit_traits_gen_N1[[maxf]][,4]) 
			sd_s2_maxf_N1 <- sd(fit_traits_gen_N1[[maxf]][,4]) 
			mean_s_dd_maxf_N1 <- mean(fit_traits_gen_N1[[maxf]][,5]) 
			sd_s_dd_maxf_N1 <- sd(fit_traits_gen_N1[[maxf]][,5]) 
			
			# sp 2
			mean_fit_maxf_N2 <- mean(fit_traits_gen_N2[[maxf]][,1]) # mean final fitness for last generation
			sd_fit_maxf_N2 <- sd(fit_traits_gen_N2[[maxf]][,1]) # sd of final fitness for last generation
			mean_f_pred_maxf_N2 <- mean(fit_traits_gen_N2[[maxf]][,2]) # mean and sd of each trait for last generation...
			sd_f_pred_maxf_N2 <- sd(fit_traits_gen_N2[[maxf]][,2]) 
			mean_s1_maxf_N2 <- mean(fit_traits_gen_N2[[maxf]][,3]) 
			sd_s1_maxf_N2 <- sd(fit_traits_gen_N2[[maxf]][,3])
			mean_s2_maxf_N2 <- mean(fit_traits_gen_N2[[maxf]][,4]) 
			sd_s2_maxf_N2 <- sd(fit_traits_gen_N2[[maxf]][,4]) 
			mean_s_dd_maxf_N2 <- mean(fit_traits_gen_N2[[maxf]][,5]) 
			sd_s_dd_maxf_N2 <- sd(fit_traits_gen_N2[[maxf]][,5]) 	
			
			if (r==1){
				master_mat_N1 <- cbind(pred=p_pred[d], run = r, fit_mean_N1 = mean_fit_maxf_N1, fit_sd_N1 = sd_fit_maxf_N1, f_pred_mean_N1 = mean_f_pred_maxf_N1, f_pred_lb_N1 = mean_f_pred_maxf_N1 - sd_f_pred_maxf_N1, f_pred_ub_N1 = mean_f_pred_maxf_N1 + sd_f_pred_maxf_N1, s1_mean_N1 = mean_s1_maxf_N1, s1_lb = mean_s1_maxf_N1 - sd_s1_maxf_N1, s1_ub_N1 = mean_s1_maxf_N1 + sd_s1_maxf_N1, s2_mean_N1 = mean_s2_maxf_N1, s2_lb_N1 = mean_s2_maxf_N1 - sd_s2_maxf_N1, s2_ub_N1 = mean_s2_maxf_N1 + sd_s2_maxf_N1, s_dd_mean_N1 = mean_s_dd_maxf_N1, s_dd_lb_N1 = mean_s_dd_maxf_N1 - sd_s_dd_maxf_N1, s_dd_ub_N1 = mean_s_dd_maxf_N1 + sd_s_dd_maxf_N1) #[[d]]
				
				master_mat_N2 <- cbind(pred=p_pred[d], run = r, fit_mean_N2 = mean_fit_maxf_N2, fit_sd_N2 = sd_fit_maxf_N2, f_pred_mean_N2 = mean_f_pred_maxf_N2, f_pred_lb_N2 = mean_f_pred_maxf_N2 - sd_f_pred_maxf_N2, f_pred_ub_N2 = mean_f_pred_maxf_N2 + sd_f_pred_maxf_N2, s1_mean_N2 = mean_s1_maxf_N2, s1_lb = mean_s1_maxf_N2 - sd_s1_maxf_N2, s1_ub_N2 = mean_s1_maxf_N2 + sd_s1_maxf_N2, s2_mean_N2 = mean_s2_maxf_N2, s2_lb_N2 = mean_s2_maxf_N2 - sd_s2_maxf_N2, s2_ub_N2 = mean_s2_maxf_N2 + sd_s2_maxf_N2, s_dd_mean_N2 = mean_s_dd_maxf_N2, s_dd_lb_N2 = mean_s_dd_maxf_N2 - sd_s_dd_maxf_N2, s_dd_ub_N2 = mean_s_dd_maxf_N2 + sd_s_dd_maxf_N2) #[[d]]
			} else{
				master_mat_N1_0 <- cbind(pred=p_pred[d], run = r, fit_mean_N1 = mean_fit_maxf_N1, fit_sd_N1 = sd_fit_maxf_N1, f_pred_mean_N1 = mean_f_pred_maxf_N1, f_pred_lb_N1 = mean_f_pred_maxf_N1 - sd_f_pred_maxf_N1, f_pred_ub_N1 = mean_f_pred_maxf_N1 + sd_f_pred_maxf_N1, s1_mean_N1 = mean_s1_maxf_N1, s1_lb = mean_s1_maxf_N1 - sd_s1_maxf_N1, s1_ub_N1 = mean_s1_maxf_N1 + sd_s1_maxf_N1, s2_mean_N1 = mean_s2_maxf_N1, s2_lb_N1 = mean_s2_maxf_N1 - sd_s2_maxf_N1, s2_ub_N1 = mean_s2_maxf_N1 + sd_s2_maxf_N1, s_dd_mean_N1 = mean_s_dd_maxf_N1, s_dd_lb_N1 = mean_s_dd_maxf_N1 - sd_s_dd_maxf_N1, s_dd_ub_N1 = mean_s_dd_maxf_N1 + sd_s_dd_maxf_N1) 
				
				master_mat_N2_0 <- cbind(pred=p_pred[d], run = r, fit_mean_N2 = mean_fit_maxf_N2, fit_sd_N2 = sd_fit_maxf_N2, f_pred_mean_N2 = mean_f_pred_maxf_N2, f_pred_lb_N2 = mean_f_pred_maxf_N2 - sd_f_pred_maxf_N2, f_pred_ub_N2 = mean_f_pred_maxf_N2 + sd_f_pred_maxf_N2, s1_mean_N2 = mean_s1_maxf_N2, s1_lb = mean_s1_maxf_N2 - sd_s1_maxf_N2, s1_ub_N2 = mean_s1_maxf_N2 + sd_s1_maxf_N2, s2_mean_N2 = mean_s2_maxf_N2, s2_lb_N2 = mean_s2_maxf_N2 - sd_s2_maxf_N2, s2_ub_N2 = mean_s2_maxf_N2 + sd_s2_maxf_N2, s_dd_mean_N2 = mean_s_dd_maxf_N2, s_dd_lb_N2 = mean_s_dd_maxf_N2 - sd_s_dd_maxf_N2, s_dd_ub_N2 = mean_s_dd_maxf_N2 + sd_s_dd_maxf_N2)
				
				master_mat_N1 <- rbind(master_mat_N1, master_mat_N1_0) #[[d]] [[d]]
				master_mat_N2 <- rbind(master_mat_N2, master_mat_N2_0) #[[d]] [[d]]
			}
			#plot(master_mat[,3]~jitter(rep(1,length(master_mat[,3])),1),xlab=c("f_pred","s_faith","s_dd"))
			
			#master_mat[,6],master_mat[,9]
			
			#points(
			print(r)
		} # end of r loop
		false_flights_dlist <- do.call(rbind, false_flights_rlist) #[[d]]
		pred_flights_dlist <- do.call(rbind, pred_flights_rlist) #[[d]]
		
		return(list(master_mat_N1, master_mat_N2, false_flights_dlist, pred_flights_dlist))
		print(d)
	} # end of d loop
	
	stopCluster(cl) # stopping the cluster, since the foreach loop is complete
	
	false_flights_all = do.call(rbind, sapply(dList,`[`,3)) #do.call(rbind, false_flights_dlist)
	pred_flights_all = do.call(rbind, sapply(dList,`[`,4)) #do.call(rbind, pred_flights_dlist)

	master_mat_N1_ALL <- do.call(rbind, sapply(dList,`[`,1)) #do.call(rbind, master_mat_N1)
	master_mat_N2_ALL <- do.call(rbind, sapply(dList,`[`,2)) #do.call(rbind, master_mat_N2)

	return(list(false_flights_all, pred_flights_all, master_mat_N1_ALL, master_mat_N2_ALL))
}

#trait1_2 = flee_fun(s1_ON=0, s2_ON=0, s_dd_ON=0)
#save(trait1_2, file="trait1_v2.rdata")
#proc.time() - ptm


cl<-makeCluster(4) # number of CPU cores
registerDoSNOW(cl)
trait2_2 = flee_fun(s1_ON=0, s2_ON=0, s_dd_ON=1)
save(trait2_2, file="trait2_v2.rdata")

cl<-makeCluster(4) # number of CPU cores
registerDoSNOW(cl)
trait3_2 = flee_fun(s1_ON=1, s2_ON=0, s_dd_ON=1)
save(trait3_2, file="trait3_v2.rdata")

cl<-makeCluster(4) # number of CPU cores
registerDoSNOW(cl)
trait4_2 = flee_fun(s1_ON=1, s2_ON=1, s_dd_ON=1)
save(trait4_2, file="trait4_v2.rdata")
proc.time() - ptm

##################################################
# Dec. 17, 2018: RUN EVERYTHING ABOVE OVERNIGHT!




s1_ON=0; s2_ON=0; s_dd_ON=0
# OK, trying to run trait 1 without foreach:
flee_fun_nopar = function(s1_ON=1, s2_ON=1, s_dd_ON=1) {
	p_pred <- seq(from=0.02, to=0.2, length.out = 6) #10
	master_mat_N1 <- list()
	master_mat_N2 <- list()
	false_flights_dlist <- list()
	pred_flights_dlist <- list()

	for (d in 1:length(p_pred)) { # d is for 'danger'!
		runs = 6
		false_flights_rlist <- list()
		pred_flights_rlist <- list()	
		for (r in 1:runs){
			# assumptions:
			N1i <- 100 # number of individuals from sp 1
			N2i <- 100
			tf <- 250 #250 time steps for a given generation
			maxg <- 25 # max group size
			e_gain <- 1 # energy units gained each time step in which the ith individual does not flee
			coef_false <- 0.20 # coefficient that determines the prob of a false alarm (a smaller value than f_pred)

			# we start with X individuals, whose evolvable traits of interest ('jumpiness', social faith, density dependence in social faith) are determined by random draws from Gaussain distributions
			# the jumpiness trait determines an individual's probability of: randomly fleeing (NOT due to a real threat, i.e., a false alarm) AND an individual's probability of detecting an attacking predator. As a first pass, let's simply make the false alarm probability much lower than the response to predators probability. Perhaps 10-20% of the response to predators probability:
			fit_N1 <- matrix(NA, ncol=tf,nrow=N1i)
			fit_N2 <- fit_N1
			fit_N1[,1] <- 1 # setting fitness of all individuals at the start at 1
			fit_N2[,1] <- 1
			flights <- list() # collecting all flight info for a generation
			attacks_all <- list()
			eaten_detect_all <- list()
			eaten_nodetect_all <- list()
			maxf <- 200 # number of generations to run the model through
			fit_traits_gen0_N1 <- list() # list to store final fitness and traits of all individuals (alive and dead) in each generation, for sp. 1
			fit_traits_gen0_N2 <- list() # list to store final fitness and traits of all individuals (alive and dead) in each generation, for sp. 1
			fit_traits_gen_N1 <- list() # list to store final fitness and traits of SURVIVING individuals in each generation
			fit_traits_gen_N2 <- list()
			trait_mean_N1 <- matrix(NA, nrow=maxf, ncol=4) # to store mean values of survivor traits for each generation
			trait_mean_N2 <- matrix(NA, nrow=maxf, ncol=4) 
			trait_sd_N1 <- matrix(NA, nrow=maxf, ncol=4) # to store sd of survivor traits for each generation
			trait_sd_N2 <- matrix(NA, nrow=maxf, ncol=4) # to store sd of survivor traits for each generation
			f_all_N1 <- list() # this list stores the new set of traits for each generation (draw from parents in prev gen). This will end up with maxf-1 items (since there was not one for the first gen).
			f_all_N2 <- list() 
			flights_master_f <- list()
			eaten_detect_master <- list()
			eaten_nodetect_master <- list()
			false_flights_flist <- list()
			pred_flights_flist <- list()

			for (f in 1:maxf){
				if (f==1){ #first gen traits are randomly determined
					f_pred_N1 <- runif(N1i, min=0, max=0.5)
					f_false_N1 <- coef_false*f_pred_N1
					f_pred_N2 <- runif(N2i, min=0, max=0.5)
					f_false_N2 <- coef_false*f_pred_N2
					if (s1_ON != 0) {
						s1_N1 <- runif(N1i, min=0, max=0.2)
						s1_N2 <- runif(N2i, min=0, max=0.2)
					} else {
						s1_N1 <- rep(0,100)
						s1_N2 <- rep(0,100)
					}
					if (s2_ON != 0) {
						s2_N1 <- runif(N1i, min=0, max=0.2)
						s2_N2 <- runif(N2i, min=0, max=0.2)
					} else {
						s2_N1 <- rep(0,100)
						s2_N2 <- rep(0,100)
					}
					if (s_dd_ON != 0) {
						s_dd_N1 <- runif(N1i, min=0, max=0.2)
						s_dd_N2 <- runif(N2i, min=0, max=0.2)
					} else {
						s_dd_N1 <- rep(0,100)
						s_dd_N2 <- rep(0,100)
					}
				} else{ # from the 2nd gen onward, traits are inherited based on parent fitness
					f_pred_N1 <- f_all_N1[[f]][,1]
					f_false_N1 <- coef_false*f_pred_N1
					f_pred_N2 <- f_all_N2[[f]][,1]
					f_false_N2 <- coef_false*f_pred_N2
					s1_N1 <- f_all_N1[[f]][,2]
					s1_N2 <- f_all_N2[[f]][,2]
					s2_N1 <- rep(0,100) #f_all_N1[[f]][,3]
					s2_N2 <- rep(0,100) #f_all_N2[[f]][,3]
					s_dd_N1 <- f_all_N1[[f]][,4]
					s_dd_N2 <- f_all_N1[[f]][,4]
				}
				false_flights_tlist <- list()
				pred_flights_tlist <- list()
				for (t in 1:(tf-1)){ # running through all the time steps of a generation
					group_sizes = seq(1:maxg)
					groups <- list()
					j = 1
					while (sum(unlist(groups))<(N1i+N2i)){ # sampling groups until both species run out
						j <- j+1
						groups[j] = sample(group_sizes, size=1, replace=TRUE)
					}
					groups = unlist(groups)
					remainder = sum(groups)-(N1i+N2i)
					groups[length(groups)] <- groups[length(groups)]-remainder # sibstracting remainder from last group
					group_vec <- list()
					for (k in 1:length(groups)){ # creating a vector that will assign each individual (arranged in a random order, below) to a group number
						vec <- rep(k,groups[k])
						group_vec[[k]] <- vec
					}	
					group_vec <- unlist(group_vec)
					seq0 <- sample(seq(1:(N1i+N2i)), size=N1i+N2i, replace=FALSE) # setting a random sequence of individuals that can each respond to the predator or to non threats (false alarm))
					groups_df <- data.frame(individual=seq0, groupID=group_vec)
					flights0 <- matrix(NA, nrow=length(unique(group_vec)), ncol=4) # collecting flight info for each group in a time step
					eaten_detect0 <- matrix(NA, nrow=length(unique(group_vec)), ncol=2)
					eaten_nodetect0 <- matrix(NA, nrow=length(unique(group_vec)), ncol=2)
					attacks_vec <- rep(NA, length(unique(group_vec)))
					for (g in 1:length(unique(group_vec))){ # running through each group
						pred <- rbinom(1, 1, p_pred[[d]]) #0.1 # determining if a predator attacks the group for the given time step
						attacks_vec[g] <- pred
						prev_flee <- 0 # setting the initial number of observable neighbor flights for a given group in a given time step
						subgroup0 <- subset(groups_df, groups_df$groupID==g)
						subgroup1 <- subset(subgroup0, subgroup0$individual<101)
						subgroup2 <- subset(subgroup0, subgroup0$individual>100)
						subgroup2$individual <- subgroup2$individual - 100 
						density_N1 <- sum(as.numeric((fit_N1[subgroup1$individual,t]>0)))
						density_N2 <- sum(as.numeric((fit_N2[subgroup2$individual,t]>0)))
						density <- density_N1 + density_N2 # determining density of LIVING individuals in group
						prev_detect <- 0 # setting the initial number of observable neighbor flights from a predator for a given group in a given time step in which a predator attacks
						eaten_detect_vec <- rep(0,length(subgroup0$individual))
						eaten_nodetect_vec <- rep(0,length(subgroup0$individual))
						for (i in 1:length(subgroup0$individual)){ # running through each individual in gth group
							if (subgroup0$individual[i]<101) { # if the individual is from sp 1
								ii <- subgroup0$individual[i]
								sp <- 1
							} else {
								ii0 <- subgroup0$individual[i]
								ii <- ii0 - 100
								sp <- 2
							}
							eaten_detect <- 0; eaten_nodetect <- 0 # setting defaults/placeholders
							if (sp == 1){
								if (fit_N1[ii,t]==0){ # if prey were EATEN in a previous time step, their fitness remains at zero, otherwise their fitness is >0 (they start out with a fitness value of 1)
									fit_N1[ii,t+1] <- 0
								} else { # if focal was NOT eaten in a previous time step
									p_flee <- f_false_N1[ii]*exp(s_dd_N1[ii]*(density-prev_flee)) + s1_N1[ii]*prev_flee*exp(s2_N1[ii]*(density-prev_flee))  # density of REMAINING individuals
									if (p_flee>1){
										p_flee=1
									}
									if (p_flee<0){
										p_flee=0
									}
									flee <- rbinom(1,1,p_flee) # flight of focal depends on both random flight and flight of neighbors
									prev_flee <- prev_flee + flee # updating flights from false alarms
									if (flee==1 | pred==1){ # if focal flees OR the predator attacks, the focal does not gain energy (eat) in this time step
										fit_N1[ii,t+1] <- fit_N1[ii,t] 
									}
									if (pred==1 & flee==0){ # if a predator attacks AND the focal did NOT already flee (via false alarm) on this time step...
										p_detect <- f_pred_N1[ii]*exp(s_dd_N1[ii]*(density-prev_flee-prev_detect))+ s1_N1[ii]*prev_detect*exp(s2_N1[ii]*(density-prev_flee-prev_detect))
										if (p_detect>1) {
											p_detect=1
										}
										if (p_detect<0) {
											p_detect=0
										}
										detect <- rbinom(1,1,p_detect)
										prev_detect <- prev_detect + detect # updating flights from predator
										if (detect==1){ # if the focal detects the predator
											peaten_detect <- 1/((length(subgroup0$individual)-prev_flee) + 10) # probability of life or death (the +10 is a bonus for early detection)
											eaten_detect <- rbinom(1,1,peaten_detect)
											eaten_detect_vec[i] <- eaten_detect
										} else { # if the focal does not detect the predator
											peaten_nodetect <- 1/((length(subgroup0$individual)-prev_flee-prev_detect)) # LOWER prob life/death (since pred was not detected) # NOTE**: could also adjust density further here: subtracting individuals who fled from the predator
											eaten_nodetect <- rbinom(1,1,peaten_nodetect)
											eaten_nodetect_vec[i] <- eaten_nodetect
											# add up each eaten for each group and store in a list at each time point
										}
										if (eaten_detect==1 | eaten_nodetect==1){
											fit_N1[ii,t+1] <- 0 # if eaten, the focal's fitness goes to 0				
										}
									}
									if (pred==0 & flee==0){ # if predator does NOT attack AND focal does NOT flee (from a false alarm)
										fit_N1[ii,t+1] <- fit_N1[ii,t] + e_gain
									}
								}
							} else {
								if (fit_N2[ii,t]==0){ # if prey were EATEN in a previous time step, their fitness remains at zero, otherwise their fitness is >0 (they start out with a fitness value of 1)
									fit_N2[ii,t+1] <- 0
								} else { # if focal was NOT eaten in a previous time step
									p_flee <- f_false_N1[ii]*exp(s_dd_N1[ii]*(density-prev_flee)) + s1_N1[ii]*prev_flee*exp(s2_N1[ii]*(density-prev_flee))  # density of REMAINING individuals
									if (p_flee>1){
										p_flee=1
									}
									if (p_flee<0){
										p_flee=0
									}
									flee <- rbinom(1,1,p_flee) # flight of focal depends on both random flight and flight of neighbors
									prev_flee <- prev_flee + flee # updating flights from false alarms
									if (flee==1 | pred==1){ # if focal flees OR the predator attacks, the focal does not gain energy (eat) in this time step
										fit_N2[ii,t+1] <- fit_N2[ii,t] 
									}
									if (pred==1 & flee==0){ # if a predator attacks AND the focal did NOT flee on this time step...
										p_detect <- f_pred_N2[ii]*exp(s_dd_N2[ii]*(density-prev_flee-prev_detect))+ s1_N2[ii]*prev_detect*exp(s2_N2[ii]*(density-prev_flee-prev_detect))									
										if (p_detect>1) {
											p_detect=1
										}
										if (p_detect<0) {
											p_detect=0
										}
										detect <- rbinom(1,1,p_detect)
										prev_detect <- prev_detect + detect # updating flights from predator
										if (detect==1){ # if the focal detects the predator
											peaten_detect <- 0 #1/((length(subgroup0$individual)-prev_flee) + 10) # probability of life or death (the +10 is a bonus for early detection)
											eaten_detect <- rbinom(1,1,peaten_detect)
											eaten_detect_vec[i] <- eaten_detect
										} else { # if the focal does not detect the predator
											peaten_nodetect <- 1/((length(subgroup0$individual)-prev_flee-prev_detect)) # LOWER prob life/death (since pred was not detected) # NOTE**: could also adjust density further here: subtracting individuals who fled from the predator
											eaten_nodetect <- rbinom(1,1,peaten_nodetect)
											eaten_nodetect_vec[i] <- eaten_nodetect
											# add up each eaten for each group and store in a list at each time point
										}
										if (eaten_detect==1 | eaten_nodetect==1){
											fit_N2[ii,t+1] <- 0 # if eaten, the focal's fitness goes to 0				
										}
									}
									if (pred==0 & flee==0){ # if predator does NOT attack AND focal does NOT flee (from a false alarm)
										fit_N2[ii,t+1] <- fit_N2[ii,t] + e_gain
									}
								}
							}
						#print(c(t,g,i))
						} # end of i loop
						# make a matrix that is col=var, row=# groups, and all this matrix to a list each time step
						flights0[g,] <- c(t, density, prev_flee, prev_detect) # recording all flights, from false alarms AND predators, for a given group and time step (groups change each time step)
						eaten_detect0[g,] <- c(t, sum(eaten_detect_vec))
						eaten_nodetect0[g,] <- c(t, sum(eaten_nodetect_vec))
					} # end of g loop
					
					flights_eaten_df <- data.frame(p_pred = rep(p_pred[d], length(flights0[,1])), run = r, gen = rep(f, length(flights0[,1])), time= flights0[,1], dens_false= flights0[,2], false_flee= flights0[,3], attack= attacks_vec, dens_pred=flights0[,2]-flights0[,3], pred_flee= flights0[,4], eaten_detect= eaten_detect0[,2], eaten_nodetect= eaten_nodetect0[,2])
					
					false_flights <- subset(flights_eaten_df, flights_eaten_df$false_flee>0)
					false_flights_tlist[[t]] <- false_flights
					pred_flights <- subset(flights_eaten_df, flights_eaten_df$pred_flee>0)
					pred_flights_tlist[[t]] <- pred_flights
					####

				} # end of t loop (completed a time step for the series of random groups of individuals
						
				false_flights_flist[[f]] <- do.call(rbind, false_flights_tlist)
				pred_flights_flist[[f]] <- do.call(rbind, pred_flights_tlist)
				
				# sp 1:
				survive_N1_0 <- cbind(fit_N1[,tf], f_pred_N1, s1_N1, s2_N1, s_dd_N1) # combining final fitness and traits for this generation
				fit_traits_gen0_N1[[f]] <- survive_N1_0
				survive_N1 <- survive_N1_0[survive_N1_0[,1]>0,] # subsetting only the survivors (those not eaten by predators)
				fit_traits_gen_N1[[f]] <- survive_N1
				trait_mean_N1[f,] <- c(mean(survive_N1[,2]), mean(survive_N1[,3]), mean(survive_N1[,4]), mean(survive_N1[,5]))
				trait_sd_N1[f,] <- c(sd(survive_N1[,2]), sd(survive_N1[,3]), sd(survive_N1[,4]), sd(survive_N1[,5]))
				survive_N1[,1] <- survive_N1[,1]/sum(survive_N1[,1]) # calculating proportional fitness (fitness of population sums to zero)
				surv_N1_df <- data.frame(survive_N1) 
				surv_N1_df$index <- seq(1:length(surv_N1_df$V1)) # creating a column of row numbers to sample from for the traits of next gen
				f_index_N1 <- sample(surv_N1_df$index, N1i, prob=surv_N1_df$V1, replace=TRUE) # sampling next gen traits weighted by parent fitness
				f_all_N1_0 <- survive_N1[f_index_N1,2:5] # traits of all Ni individuals for next gen	
				f_mut_N1_0 <- cbind(rnorm(100, mean=0, sd=0.04*abs(mean(f_all_N1_0[,1]))), rnorm(100, mean=0, sd=0.04*abs(mean(f_all_N1_0[,2]))), rnorm(100, mean=0, sd=0.04*abs(mean(f_all_N1_0[,3]))), rnorm(100, mean=0, sd=0.04*abs(mean(f_all_N1_0[,4]))))	# introducing mutation
				if (s1_ON == 0) { 
					f_mut_N1_0[,2] = rep(0, length(f_mut_N1_0[,2])) 
				}
				if (s2_ON == 0) { 
					f_mut_N1_0[,3] = rep(0, length(f_mut_N1_0[,3])) 
				}
				if (s_dd_ON == 0) { 
					f_mut_N1_0[,4] = rep(0, length(f_mut_N1_0[,4])) 
				}
				f_mut_N1 <- f_all_N1_0 + f_mut_N1_0
				f_mut_N1[,1] <- replace(f_mut_N1[,1], f_mut_N1[,1] < 0, 0)
				f_mut_N1[,2] <- replace(f_mut_N1[,2], f_mut_N1[,2] < 0, 0)
				f_all_N1[[f+1]] <- f_mut_N1
				
				# sp 2 (repeating above)
				survive_N2_0 <- cbind(fit_N2[,tf], f_pred_N2, s1_N2, s2_N2, s_dd_N2) # combining final fitness and traits for this generation
				fit_traits_gen0_N2[[f]] <- survive_N2_0
				survive_N2 <- survive_N2_0[survive_N2_0[,1]>0,] # subsetting only the survivors (those not eaten by predators)
				fit_traits_gen_N2[[f]] <- survive_N2
				trait_mean_N2[f,] <- c(mean(survive_N2[,2]), mean(survive_N2[,3]), mean(survive_N2[,4]), mean(survive_N2[,5]))
				trait_sd_N2[f,] <- c(sd(survive_N2[,2]), sd(survive_N2[,3]), sd(survive_N2[,4]), sd(survive_N2[,5]))
				survive_N2[,1] <- survive_N2[,1]/sum(survive_N2[,1]) # calculating proportional fitness (fitness of population sums to zero)
				surv_N2_df <- data.frame(survive_N2) 
				surv_N2_df$index <- seq(1:length(surv_N2_df$V1)) # creating a column of row numbers to sample from for the traits of next gen
				f_index_N2 <- sample(surv_N2_df$index, N2i, prob=surv_N2_df$V1, replace=TRUE) # sampling next gen traits weighted by parent fitness
				f_all_N2_0 <- survive_N2[f_index_N2,2:5] # traits of all Ni individuals for next gen	
				f_mut_N2_0 <- cbind(rnorm(100, mean=0, sd=0.04*abs(mean(f_all_N2_0[,1]))), rnorm(100, mean=0, sd=0.04*abs(mean(f_all_N2_0[,2]))), rnorm(100, mean=0, sd=0.04*abs(mean(f_all_N2_0[,3]))), rnorm(100, mean=0, sd=0.04*abs(mean(f_all_N2_0[,4])))) # introducing mutation
				if (s1_ON == 0) { 
					f_mut_N1_0[,2] = rep(0, length(f_mut_N1_0[,2])) 
				}
				if (s2_ON == 0) { 
					f_mut_N1_0[,3] = rep(0, length(f_mut_N1_0[,3])) 
				}
				if (s_dd_ON == 0) { 
					f_mut_N1_0[,4] = rep(0, length(f_mut_N1_0[,4])) 
				}
				f_mut_N2 <- f_all_N2_0 + f_mut_N2_0
				f_mut_N2[,1] <- replace(f_mut_N2[,1], f_mut_N2[,1] < 0, 0)
				f_mut_N2[,2] <- replace(f_mut_N2[,2], f_mut_N2[,2] < 0, 0)
				f_all_N2[[f+1]] <- f_mut_N2
			
			} # end of f loop
			
			#flights_master_r[[r]] <- flights_master_f
			false_flights_rlist[[r]] <- do.call(rbind, false_flights_flist)
			pred_flights_rlist[[r]] <- do.call(rbind, pred_flights_flist)
			
			# sp 1
			mean_fit_maxf_N1 <- mean(fit_traits_gen_N1[[maxf]][,1]) # mean final fitness for last generation
			sd_fit_maxf_N1 <- sd(fit_traits_gen_N1[[maxf]][,1]) # sd of final fitness for last generation
			mean_f_pred_maxf_N1 <- mean(fit_traits_gen_N1[[maxf]][,2]) # mean and sd of each trait for last generation...
			sd_f_pred_maxf_N1 <- sd(fit_traits_gen_N1[[maxf]][,2]) 
			mean_s1_maxf_N1 <- mean(fit_traits_gen_N1[[maxf]][,3]) 
			sd_s1_maxf_N1 <- sd(fit_traits_gen_N1[[maxf]][,3])
			mean_s2_maxf_N1 <- mean(fit_traits_gen_N1[[maxf]][,4]) 
			sd_s2_maxf_N1 <- sd(fit_traits_gen_N1[[maxf]][,4]) 
			mean_s_dd_maxf_N1 <- mean(fit_traits_gen_N1[[maxf]][,5]) 
			sd_s_dd_maxf_N1 <- sd(fit_traits_gen_N1[[maxf]][,5]) 
			
			# sp 2
			mean_fit_maxf_N2 <- mean(fit_traits_gen_N2[[maxf]][,1]) # mean final fitness for last generation
			sd_fit_maxf_N2 <- sd(fit_traits_gen_N2[[maxf]][,1]) # sd of final fitness for last generation
			mean_f_pred_maxf_N2 <- mean(fit_traits_gen_N2[[maxf]][,2]) # mean and sd of each trait for last generation...
			sd_f_pred_maxf_N2 <- sd(fit_traits_gen_N2[[maxf]][,2]) 
			mean_s1_maxf_N2 <- mean(fit_traits_gen_N2[[maxf]][,3]) 
			sd_s1_maxf_N2 <- sd(fit_traits_gen_N2[[maxf]][,3])
			mean_s2_maxf_N2 <- mean(fit_traits_gen_N2[[maxf]][,4]) 
			sd_s2_maxf_N2 <- sd(fit_traits_gen_N2[[maxf]][,4]) 
			mean_s_dd_maxf_N2 <- mean(fit_traits_gen_N2[[maxf]][,5]) 
			sd_s_dd_maxf_N2 <- sd(fit_traits_gen_N2[[maxf]][,5]) 	
			
			if (r==1){
				master_mat_N1[[d]] <- cbind(pred=p_pred[d], run = r, fit_mean_N1 = mean_fit_maxf_N1, fit_sd_N1 = sd_fit_maxf_N1, f_pred_mean_N1 = mean_f_pred_maxf_N1, f_pred_lb_N1 = mean_f_pred_maxf_N1 - sd_f_pred_maxf_N1, f_pred_ub_N1 = mean_f_pred_maxf_N1 + sd_f_pred_maxf_N1, s1_mean_N1 = mean_s1_maxf_N1, s1_lb = mean_s1_maxf_N1 - sd_s1_maxf_N1, s1_ub_N1 = mean_s1_maxf_N1 + sd_s1_maxf_N1, s2_mean_N1 = mean_s2_maxf_N1, s2_lb_N1 = mean_s2_maxf_N1 - sd_s2_maxf_N1, s2_ub_N1 = mean_s2_maxf_N1 + sd_s2_maxf_N1, s_dd_mean_N1 = mean_s_dd_maxf_N1, s_dd_lb_N1 = mean_s_dd_maxf_N1 - sd_s_dd_maxf_N1, s_dd_ub_N1 = mean_s_dd_maxf_N1 + sd_s_dd_maxf_N1) #[[d]]
				
				master_mat_N2[[d]] <- cbind(pred=p_pred[d], run = r, fit_mean_N2 = mean_fit_maxf_N2, fit_sd_N2 = sd_fit_maxf_N2, f_pred_mean_N2 = mean_f_pred_maxf_N2, f_pred_lb_N2 = mean_f_pred_maxf_N2 - sd_f_pred_maxf_N2, f_pred_ub_N2 = mean_f_pred_maxf_N2 + sd_f_pred_maxf_N2, s1_mean_N2 = mean_s1_maxf_N2, s1_lb = mean_s1_maxf_N2 - sd_s1_maxf_N2, s1_ub_N2 = mean_s1_maxf_N2 + sd_s1_maxf_N2, s2_mean_N2 = mean_s2_maxf_N2, s2_lb_N2 = mean_s2_maxf_N2 - sd_s2_maxf_N2, s2_ub_N2 = mean_s2_maxf_N2 + sd_s2_maxf_N2, s_dd_mean_N2 = mean_s_dd_maxf_N2, s_dd_lb_N2 = mean_s_dd_maxf_N2 - sd_s_dd_maxf_N2, s_dd_ub_N2 = mean_s_dd_maxf_N2 + sd_s_dd_maxf_N2) #[[d]]
			} else{
				master_mat_N1_0 <- cbind(pred=p_pred[d], run = r, fit_mean_N1 = mean_fit_maxf_N1, fit_sd_N1 = sd_fit_maxf_N1, f_pred_mean_N1 = mean_f_pred_maxf_N1, f_pred_lb_N1 = mean_f_pred_maxf_N1 - sd_f_pred_maxf_N1, f_pred_ub_N1 = mean_f_pred_maxf_N1 + sd_f_pred_maxf_N1, s1_mean_N1 = mean_s1_maxf_N1, s1_lb = mean_s1_maxf_N1 - sd_s1_maxf_N1, s1_ub_N1 = mean_s1_maxf_N1 + sd_s1_maxf_N1, s2_mean_N1 = mean_s2_maxf_N1, s2_lb_N1 = mean_s2_maxf_N1 - sd_s2_maxf_N1, s2_ub_N1 = mean_s2_maxf_N1 + sd_s2_maxf_N1, s_dd_mean_N1 = mean_s_dd_maxf_N1, s_dd_lb_N1 = mean_s_dd_maxf_N1 - sd_s_dd_maxf_N1, s_dd_ub_N1 = mean_s_dd_maxf_N1 + sd_s_dd_maxf_N1) 
				
				master_mat_N2_0 <- cbind(pred=p_pred[d], run = r, fit_mean_N2 = mean_fit_maxf_N2, fit_sd_N2 = sd_fit_maxf_N2, f_pred_mean_N2 = mean_f_pred_maxf_N2, f_pred_lb_N2 = mean_f_pred_maxf_N2 - sd_f_pred_maxf_N2, f_pred_ub_N2 = mean_f_pred_maxf_N2 + sd_f_pred_maxf_N2, s1_mean_N2 = mean_s1_maxf_N2, s1_lb = mean_s1_maxf_N2 - sd_s1_maxf_N2, s1_ub_N2 = mean_s1_maxf_N2 + sd_s1_maxf_N2, s2_mean_N2 = mean_s2_maxf_N2, s2_lb_N2 = mean_s2_maxf_N2 - sd_s2_maxf_N2, s2_ub_N2 = mean_s2_maxf_N2 + sd_s2_maxf_N2, s_dd_mean_N2 = mean_s_dd_maxf_N2, s_dd_lb_N2 = mean_s_dd_maxf_N2 - sd_s_dd_maxf_N2, s_dd_ub_N2 = mean_s_dd_maxf_N2 + sd_s_dd_maxf_N2)
				
				master_mat_N1[[d]] <- rbind(master_mat_N1[[d]], master_mat_N1_0) #[[d]] [[d]]
				master_mat_N2[[d]] <- rbind(master_mat_N2[[d]], master_mat_N2_0) #[[d]] [[d]]
			}
			#plot(master_mat[,3]~jitter(rep(1,length(master_mat[,3])),1),xlab=c("f_pred","s_faith","s_dd"))
			
			#master_mat[,6],master_mat[,9]
			
			#points(

		} # end of r loop
		false_flights_dlist[[d]] <- do.call(rbind, false_flights_rlist) #[[d]]
		pred_flights_dlist[[d]] <- do.call(rbind, pred_flights_rlist) #[[d]]
		print(d)
		#return(list(master_mat_N1, master_mat_N2, false_flights_dlist, pred_flights_dlist))
	} # end of d loop
	
	#stopCluster(cl) # stopping the cluster, since the foreach loop is complete
	
	false_flights_all = do.call(rbind, false_flights_dlist)
	pred_flights_all = do.call(rbind, pred_flights_dlist)

	master_mat_N1_ALL <- do.call(rbind, master_mat_N1)
	master_mat_N2_ALL <- do.call(rbind, master_mat_N2)

	return(list(false_flights_all, pred_flights_all, master_mat_N1_ALL, master_mat_N2_ALL))
}

ptm <- proc.time()
trait1 = flee_fun_nopar(s1_ON=0, s2_ON=0, s_dd_ON=0)
save(trait1, file="trait1_v1.rdata")
proc.time() - ptm








# Oct 8: 1 idea would be to regress each trait from one species against that trait in the other. What we may see is that traits are positively correlated when predation is low, and negatively correlated when predation is high (meaning strategies are the same

par(mfrow=c(1,3))
rbPal <- colorRampPalette(c('blue','red'))
cols = rbPal(runs)[as.numeric(cut(master_mat_ALL$pred, breaks = runs))]
# f_pred:
plot(master_mat_ALL$f_pred_mean_N1~master_mat_ALL$f_pred_mean_N2, col=cols, pch= c(rep(1,runs),rep(2,runs),rep(3,runs),rep(4,runs)), xlab="sp 2 responsiveness to personal info (f_pred)", ylab="sp 1 responsiveness to personal info (f_pred)", ylim=c(min(c(master_mat_ALL$f_pred_mean_N1,master_mat_ALL$f_pred_mean_N2)),max(c(master_mat_ALL$f_pred_mean_N1,master_mat_ALL$f_pred_mean_N2))), xlim=c(min(c(master_mat_ALL$f_pred_mean_N1,master_mat_ALL$f_pred_mean_N2)), max(c(master_mat_ALL$f_pred_mean_N1,master_mat_ALL$f_pred_mean_N2))))
abline(h=0,lty=3)
abline(v=0,lty=3)
# s1:
plot(master_mat_ALL$s1_mean_N1~master_mat_ALL$s1_mean_N2, col=cols, pch= c(rep(1,runs),rep(2,runs),rep(3,runs),rep(4,runs)), xlab="sp 2 base responsiveness to neighbor flight (s1)", ylab="sp 1 base responsiveness to neighbor flight (s1)", ylim=c(min(c(master_mat_ALL$s1_mean_N1,master_mat_ALL$s1_mean_N2)), max(c(master_mat_ALL$s1_mean_N1,master_mat_ALL$s1_mean_N2))), xlim=c(min(c(master_mat_ALL$s1_mean_N1,master_mat_ALL$s1_mean_N2)), max(c(master_mat_ALL$s1_mean_N1,master_mat_ALL$s1_mean_N2))))
abline(h=0,lty=3)
abline(v=0,lty=3)
# s2:
#plot(master_mat_ALL$s2_mean_N1~master_mat_ALL$s2_mean_N2, col=cols, pch= c(rep(1,runs),rep(2,runs),rep(3,runs),rep(4,runs)), xlab="sp 2 dd in response to neighbor flights (s2)", ylab="sp 1 dd in response to neighbor flights (s2)", ylim=c(min(c(master_mat_ALL$s2_mean_N1,master_mat_ALL$s2_mean_N2)), max(c(master_mat_ALL$s2_mean_N1,master_mat_ALL$s2_mean_N2))), xlim=c(min(c(master_mat_ALL$s2_mean_N1,master_mat_ALL$s2_mean_N2)), max(c(master_mat_ALL$s2_mean_N1,master_mat_ALL$s2_mean_N2))))
#abline(h=0,lty=3)
#abline(v=0,lty=3)
# s_dd:
plot(master_mat_ALL$s_dd_mean_N1~master_mat_ALL$s_dd_mean_N2, col=cols, pch= c(rep(1,runs),rep(2,runs),rep(3,runs),rep(4,runs)), xlab="sp 2 dd in flight response (s_dd)", ylab="sp 1 dd in flight response (s_dd)", ylim=c(min(c(master_mat_ALL$s_dd_mean_N1,master_mat_ALL$s_dd_mean_N2)), max(c(master_mat_ALL$s_dd_mean_N1,master_mat_ALL$s_dd_mean_N2))), xlim=c(min(c(master_mat_ALL$s_dd_mean_N1,master_mat_ALL$s_dd_mean_N2)), max(c(master_mat_ALL$s_dd_mean_N1,master_mat_ALL$s_dd_mean_N2))))
abline(h=0,lty=3)
abline(v=0,lty=3)



# saving output
save(master_mat_N1, file="evo_p02_2_r2_N1_3TRAITS.rdata")
save(master_mat_N2, file="evo_p02_2_r2_N2_3TRAITS.rdata")

load("evo_p02_2_r2_N1_3TRAITS.rdata")
head(master_mat_N1[[1]])



# collecting flight info:
false_flights_ALL <- do.call(rbind, false_flights_dlist)
save(false_flights_ALL,file="falsef_evo_p02_2_r2_N1N2_3TRAITS.rdata")

pred_flights_ALL <- do.call(rbind, pred_flights_dlist)
save(pred_flights_ALL,file="predf_evo_p02_2_r2_N1N2_3TRAITS.rdata")

#master_mat
proc.time() - ptm







dev.new()
par(mfrow=c(1,4))
traits <- c("jumpiness (from personal info; 'f_flight')","response to neighbor flights (s1)","dd in response to neighbor flights (s2)", "dd in focal flight (s_dd)")
rbPal <- colorRampPalette(c('red','blue'))
ylo <- c(0,0,-5,-1.3)#c(0,0,-2,-0.2)
yhi <- c(2.5,0.5,4,1) #c(0.5,0.2,2,0.2)
for (p in 1:4){ # number of evolvable traits
	plot(0, type='n', xlim=c(min(p_pred),max(p_pred)), ylim=c(ylo[p],yhi[p]), ylab="",xlab="")
	for  (d in 1:length(p_pred)){
		pred_prob <- p_pred[d]
		#cols = rbPal(runs)[as.numeric(cut(master_mat[[d]][,3],breaks = runs))]
		cols = rbPal(runs)[as.numeric(cut(seq(1:runs),breaks = runs))]
		points(master_mat[[d]][,p+4+2*(p-1)]~rep(pred_prob,length(master_mat[[d]][,p+4+2*(p-1)])), col=cols, pch=seq(1:runs), cex=2)
		abline(h=0,lty=3)
		#points(master_mat[[d]][,p+5+2*(p-1)]~rep(pred_prob,length(master_mat[[d]][,p+3+2*(p-1)])), col=cols, pch="-")
		#points(master_mat[[d]][,p+6+2*(p-1)]~rep(pred_prob,length(master_mat[[d]][,p+3+2*(p-1)])), col=cols, pch="-")
	}
	title(xlab="predation probability", ylab=traits[p])
}

dev.new()
par(mfrow=c(1,4))
for (p in 1:4){
	plot(0, type='n', xlim=c(min(p_pred),max(p_pred[length(p_pred)])), ylim=c(ylo[p],yhi[p]), ylab="",xlab="")
	for  (d in 1:length(p_pred)){
		pred_prob <- p_pred[d]
		#cols = rbPal(10)[as.numeric(cut(master_mat[[d]][,2],breaks = 10))]
		#cols = seq(1:runs)
		points(mean(master_mat[[d]][,p+4+2*(p-1)])~pred_prob, pch=16, cex=2)
		points(mean(master_mat[[d]][,p+4+2*(p-1)])-sd(master_mat[[d]][,p+5+2*(p-1)])~pred_prob, col=cols, pch="-")
		points(mean(master_mat[[d]][,p+4+2*(p-1)])+sd(master_mat[[d]][,p+5+2*(p-1)])~pred_prob, col=cols, pch="-")
		abline(h=0,lty=3)
	}
	title(xlab="predation probability", ylab=traits[p])
}

# Oct. 7: Now thinking of a way to uncover distinct strategies that may exist in the data. One way to do this: look at correlation between s2 and s_dd:

dev.new()
par(mfrow=c(2,5))
trait_compare = c(1,4)
third_trait = 3
for (d in 1:length(p_pred)){
	plot(0, type='n',  ylab="",xlab="", ylim=c(min(master_mat[[d]][,trait_compare[1]+4+2*(trait_compare[1]-1)]),max(master_mat[[d]][,trait_compare[1]+4+2*(trait_compare[1]-1)])), xlim=c(min(master_mat[[d]][,trait_compare[2]+4+2*(trait_compare[2]-1)]),max(master_mat[[d]][,trait_compare[2]+4+2*(trait_compare[2]-1)]))) #ylim=c(0,0.5), xlim=c(-4,4)
	pred_prob <- p_pred[d]
	cols = rbPal(runs)[as.numeric(cut(master_mat[[d]][,third_trait+4+2*(third_trait-1)], breaks = runs))] # a 3rd trait
	# fitness: cols = rbPal(runs)[as.numeric(cut(master_mat[[d]][,3],breaks = runs))]
	#cols = rbPal(runs)[as.numeric(cut(seq(1:runs),breaks = runs))]
	points(master_mat[[d]][,trait_compare[1]+4+2*(trait_compare[1]-1)]~master_mat[[d]][,trait_compare[2]+4+2*(trait_compare[2]-1)], col=cols, cex=2)
	abline(h=0,lty=3)
	abline(v=0,lty=3)
	title(main=paste("p_pred=",round(p_pred[d],3)),xlab=traits[trait_compare[2]], ylab=traits[trait_compare[1]])

	#points(master_mat[[d]][,p+5+2*(p-1)]~rep(pred_prob,length(master_mat[[d]][,p+3+2*(p-1)])), col=cols, pch="-")
	#points(master_mat[[d]][,p+6+2*(p-1)]~rep(pred_prob,length(master_mat[[d]][,p+3+2*(p-1)])), col=cols, pch="-")
}

# These trait comparison plots didn't seem to reveal clear distinct strategies. I should think of how else I could identify these...

# Next: 2 spp.




# for this plot, I should combine all matrices of the list, calculate the fitness color range from the full set, and provide a legend of the color values - ACTUALLY, i think this would just show that fitness declines with higher predation, so i may leave colors as is
# Sept. 21, 2018: Still need to do above, also need to:
# 1) run to more generations (to reduce error bars) and maybe play with initial parm value ranges
# 2) once I feel good about above, work on the accompanying animations

# Sept. 24: plot box plots or means with error bars. Also want to plot some examples of positive versus negative dd (4th parm) to see if the two different strategies I'd revealed before re-emerge.
# I should create a plot that instead of coloring by fitness, colors by order (so we can see which values of the different traits were associated with one another).
# Next: stochastic death process model?
# Using above model to allow social individuals to group more


### Oct. 7, 2018: Making 'response to false alarm' X 'g1, gf_PLo, gf_Phi' plot:
falsef_gen_i <- subset(false_flights_ALL, false_flights_ALL$gen==1 & false_flights_ALL$time==1) #false_flights_ALL$time==1
falsef_gen_f <- subset(false_flights_ALL, false_flights_ALL$gen>maxf-10 & false_flights_ALL$time==1)
subset(pred_flights_ALL, pred_flights_ALL$gen<11 & pred_flights_ALL$p_pred==0.1)
predf_gen_if0 <- subset(pred_flights_ALL, pred_flights_ALL$gen<11 | pred_flights_ALL$gen>maxf-10) # sampling 1st and last 10 generations
predf_gen_if_20 <- subset(predf_gen_if0, predf_gen_if0$dens_pred>18)
predf_gen_if_05 <- subset(predf_gen_if0, predf_gen_if0$dens_pred== 5)

# first, responses to a false alarm will just be the number of false alarm flights in a time point -1 (minus the initial false alarmer):
falsef_geni_20_PLo <- subset(false_flights_ALL, false_flights_ALL$gen==1 & false_flights_ALL$dens_false==20 & false_flights_ALL$p_pred==0.05)
falsef_geni_20_PHi <- subset(false_flights_ALL, false_flights_ALL$gen==1 & false_flights_ALL$dens_false==20 & false_flights_ALL$p_pred==0.1)

falsef_genf_20_PLo <- subset(false_flights_ALL, false_flights_ALL$gen==maxf & false_flights_ALL$dens_false==20 & false_flights_ALL$p_pred==0.05)
falsef_genf_20_PHi <- subset(false_flights_ALL, false_flights_ALL$gen==maxf & false_flights_ALL$dens_false==20 & false_flights_ALL$p_pred==0.1)

c(mean(falsef_geni_20_PLo[,6]-1), mean(falsef_geni_20_PHi[,6]-1), mean(falsef_genf_20_PLo[,6]-1), mean(falsef_genf_20_PHi[,6]-1))
c(length(falsef_geni_20_PLo[,6]), length(falsef_geni_20_PHi[,6]), length(falsef_genf_20_PLo[,6]), length(falsef_genf_20_PHi[,6]))

# now repeating above for predator flights, to see if this can help explain false alarm patterns:
predf_geni_20_PLo <- subset(pred_flights_ALL, pred_flights_ALL$gen==1 & pred_flights_ALL$dens_pred>9 & pred_flights_ALL$p_pred==0.05)
predf_geni_20_PHi <- subset(pred_flights_ALL, pred_flights_ALL$gen==1 & pred_flights_ALL$dens_pred>9 & pred_flights_ALL$p_pred==0.1)

predf_genf_20_PLo <- subset(pred_flights_ALL, pred_flights_ALL$gen==maxf & pred_flights_ALL$dens_pred==20 & pred_flights_ALL$p_pred==0.05)
predf_genf_20_PHi <- subset(pred_flights_ALL, pred_flights_ALL$gen==maxf & pred_flights_ALL$dens_pred==20 & pred_flights_ALL$p_pred==0.1)

c(mean(predf_geni_20_PLo[,9]-1), mean(predf_geni_20_PHi[,9]-1), mean(predf_genf_20_PLo[,9]-1), mean(predf_genf_20_PHi[,9]-1))
c(length(predf_geni_20_PLo[,9]), length(predf_geni_20_PHi[,9]), length(predf_genf_20_PLo[,9]), length(predf_genf_20_PHi[,9]))



#-------------------------------------------------------------------------------------------------------------

# Next, 'response to a false/REAL flight' X 'prey density' plots:
par(mfrow=c(2,6))

# false flights first
falsef_geni_PLo <- subset(false_flights_ALL, false_flights_ALL$gen==1 & false_flights_ALL$p_pred==0.02)
with(falsef_geni_PLo, boxplot((false_flee-1)~dens_false, ylim=c(0,25), ylab="# of flights following false alarm flight"))
points(aggregate(falsef_geni_PLo[, 6]-1, list(falsef_geni_PLo$dens_false), mean), col="red")
title(paste("gen=1 p_pred=", p_pred[1]))
text(10, 20, paste("# of groups w/ >=1 \nfalse flights=",length(falsef_geni_PLo[,1])))
abline(a=-1,b=1, col="grey")

falsef_geni_PHi <- subset(false_flights_ALL, false_flights_ALL$gen==1 & false_flights_ALL$p_pred==0.2)
with(falsef_geni_PHi, boxplot((false_flee-1)~dens_false, ylim=c(0,25)))
points(aggregate(falsef_geni_PHi[, 6]-1, list(falsef_geni_PHi$dens_false), mean), col="red")
title(paste("gen=",1," p_pred=", p_pred[length(p_pred)]))
text(10, 20, paste("# of groups w/ >=1 \nfalse flights=",length(falsef_geni_PHi[,1])))
abline(a=-1,b=1, col="grey")

falsef_genf_PLo <- subset(false_flights_ALL, false_flights_ALL$gen==maxf & false_flights_ALL$p_pred==0.02)
with(falsef_genf_PLo, boxplot((false_flee-1)~dens_false, ylim=c(0,25)))
points(aggregate(falsef_genf_PLo[, 6]-1, list(falsef_genf_PLo$dens_false), mean), col="red")
title(paste("gen=",maxf," p_pred=", p_pred[1]))
text(10, 20, paste("# of groups w/ >=1 \nfalse flights=",length(falsef_genf_PLo[,1])))
abline(a=-1,b=1, col="grey")

falsef_genf_Pmid1 <- subset(false_flights_ALL, false_flights_ALL$gen==maxf & false_flights_ALL$p_pred==0.08)
with(falsef_genf_Pmid1, boxplot((false_flee-1)~dens_false, ylim=c(0,25)))
points(aggregate(falsef_genf_Pmid1[, 6]-1, list(falsef_genf_Pmid1$dens_false), mean), col="red")
title(paste("gen=",maxf," p_pred=", p_pred[2]))
text(10, 20, paste("# of groups w/ >=1 \nfalse flights=",length(falsef_genf_Pmid1[,1])))
abline(a=-1,b=1, col="grey")

falsef_genf_Pmid2 <- subset(false_flights_ALL, false_flights_ALL$gen==maxf & false_flights_ALL$p_pred==0.14)
with(falsef_genf_Pmid2, boxplot((false_flee-1)~dens_false, ylim=c(0,25)))
points(aggregate(falsef_genf_Pmid2[, 6]-1, list(falsef_genf_Pmid2$dens_false), mean), col="red")
title(paste("gen=",maxf," p_pred=", p_pred[3]))
text(10, 20, paste("# of groups w/ >=1 \nfalse flights=",length(falsef_genf_Pmid2[,1])))
abline(a=-1,b=1, col="grey")

falsef_genf_PHi <- subset(false_flights_ALL, false_flights_ALL$gen==maxf & false_flights_ALL$p_pred==0.2)
with(falsef_genf_PHi, boxplot((false_flee-1)~dens_false, ylim=c(0,25)))
points(aggregate(falsef_genf_PHi[, 6]-1, list(falsef_genf_PHi$dens_false), mean), col="red")
title(paste("gen=",maxf," ,p_pred=", p_pred[length(p_pred)]))
text(10, 20, paste("# of groups w/ >=1 \nfalse flights=",length(falsef_genf_PHi[,1])))
abline(a=-1,b=1, col="grey")

# pred flights second
predf_geni_PLo <- subset(pred_flights_ALL, pred_flights_ALL$gen<11 & pred_flights_ALL$p_pred==0.02)
with(predf_geni_PLo, boxplot((pred_flee-1)~dens_pred, ylim=c(0,25), ylab="# of flights following predator avoidance flight (true alarm)", xlab="pre-flight prey density"))
points(aggregate(predf_geni_PLo[, 9]-1, list(predf_geni_PLo$dens_pred), mean), col="red")
title("gen=1-10")
abline(a=-1,b=1, col="grey")
#text(8, 20, paste("# of groups attacked=",length(predf_geni_PLo[,1])))

predf_geni_PHi <- subset(pred_flights_ALL, pred_flights_ALL$gen<11 & pred_flights_ALL$p_pred==0.2)
with(predf_geni_PHi, boxplot((pred_flee-1)~dens_pred, ylim=c(0,25), xlab="pre-flight prey density"))
points(aggregate(predf_geni_PHi[, 9]-1, list(predf_geni_PHi$dens_pred), mean), col="red")
title("gen=1-10")
abline(a=-1,b=1, col="grey")
#text(8, 20, paste("# of groups attacked=",length(predf_geni_PHi[,1])))

predf_genf_PLo <- subset(pred_flights_ALL, pred_flights_ALL$gen==maxf & pred_flights_ALL$p_pred==0.02)
with(predf_genf_PLo, boxplot((pred_flee-1)~dens_pred, ylim=c(0,25), xlab="pre-flight prey density"))
points(aggregate(predf_genf_PLo[, 9]-1, list(predf_genf_PLo$dens_pred), mean), col="red")
text(8, 20, paste("# of groups \nattacked=",length(predf_genf_PLo[,1])))
abline(a=-1,b=1, col="grey")

predf_genf_Pmid1 <- subset(pred_flights_ALL, pred_flights_ALL$gen==maxf & pred_flights_ALL$p_pred==0.08)
with(predf_genf_Pmid1, boxplot((pred_flee-1)~dens_pred, ylim=c(0,25), xlab="pre-flight prey density"))
points(aggregate(predf_genf_Pmid1[, 9]-1, list(predf_genf_Pmid1$dens_pred), mean), col="red")
text(8, 20, paste("# of groups \nattacked=",length(predf_genf_Pmid1[,1])))
abline(a=-1,b=1, col="grey")

predf_genf_Pmid2 <- subset(pred_flights_ALL, pred_flights_ALL$gen==maxf & pred_flights_ALL$p_pred==0.14)
with(predf_genf_Pmid2, boxplot((pred_flee-1)~dens_pred, ylim=c(0,25), xlab="pre-flight prey density"))
points(aggregate(predf_genf_Pmid2[, 9]-1, list(predf_genf_Pmid2$dens_pred), mean), col="red")
text(8, 20, paste("# of groups \nattacked=",length(predf_genf_Pmid2[,1])))
abline(a=-1,b=1, col="grey")

predf_genf_PHi <- subset(pred_flights_ALL, pred_flights_ALL$gen==maxf & pred_flights_ALL$p_pred==0.2)
with(predf_genf_PHi, boxplot((pred_flee-1)~dens_pred, ylim=c(0,25), xlab="pre-flight prey density"))
points(aggregate(predf_genf_PHi[, 9]-1, list(predf_genf_PHi$dens_pred), mean), col="red")
text(8, 20, paste("# of groups \nattacked=",length(predf_genf_PHi[,1])))
abline(a=-1,b=1, col="grey")


pushViewport(viewport())
grid.lines(x = 0.34, y = c(0,1), gp = gpar(col = "black"))
popViewport()





# to select only survivors from the complete 'fit_traits_gen' list (which contains fitnesses for all individuals in every generation), I'd modify this line:
dev.new()
par(mfrow=c(1,4))
traits <- c("jumpiness","sociality","dd1","dd2")
gg_trait_ls <- list()
mean_ff <- round(mean(fit_traits_gen[[maxf]][,1]), 2) # mean final fitness for last generation
sd_ff <- round(sd(fit_traits_gen[[maxf]][,1]), 2) # sd of final fitness for last generation
for (p in 1:4){
	df <- data.frame(generation=seq(1:length(trait_mean[,p])), trait_mean=trait_mean[,p], lb=trait_mean[,p]-trait_sd[,p], ub=trait_mean[,p]+trait_sd[,p])
	gg_trait_ls[[p]] <- ggplot(df, aes(generation)) + 
	geom_line(aes(y=trait_mean), colour="black") + 
	geom_ribbon(aes(ymin=lb, ymax=ub), alpha=0.2) +
	labs(y = traits[p])
}
#proc.time() - ptm
multiplot(gg_trait_ls[[1]] +	ggtitle(paste("Fitness=",mean_ff,"+/-",sd_ff,"SD")),gg_trait_ls[[2]],gg_trait_ls[[3]])
# SWEET! THIS WORKED
# Looks like alternative strategies emerge. But I should test the fitness (and variance thereof) of each.
# 153, 166 for pos dd trait, 121 (with high pos dd!), 
# 146 for neg dd trait
# NOTE: when predation is high, there appears to be more variability in the resulting strategies and final fitness (perhaps because there is a much greater chance that even individuals with ideal traits are randomly executed and have less an ability to guide the population).

# Sept. 20, 2018: Goals:
# 1) Run above a bunch and log the stable traits. Plot mean stable traits (color them by fitness) for each run over predation probability. Look for bifurcation.









# NEXT: randomly sampling a large group: 20 individuals and looking at response to false alarms and predators between gen1 and genf:


# for plotting flights, I could plot all individuals in a group in a row and have them move to a new row if they fled, and I could do this sequentially, so that you can see the cascade. I could also color individuals cold to warm based on their: 1) jumpiness, 2) social sensitivity (this would require keeping track of individuals, which I've not done -- I don't think this will be worth doing, since we're primarily interested in the evolution, which I can show in more convincing ways).
# to code this, I want to first break the 'flights' dataframe up into 'false alarm' and 'predator' flights:


# Plots for animations of flights (false and true)
flights_all <- do.call(rbind, flights)
eaten_detect_mat <- do.call(rbind, eaten_detect_all)
eaten_nodetect_mat <- do.call(rbind, eaten_nodetect_all)

flights_eaten_df <- data.frame(time= flights_all[,1], dens_false= flights_all[,2], false_flee= flights_all[,3], attack= unlist(attacks_all), dens_pred=flights_all[,2]-flights_all[,3], pred_flee= flights_all[,4], eaten_detect= eaten_detect_mat[,2], eaten_nodetect= eaten_nodetect_mat[,2])

false_flights <- subset(flights_eaten_df, flights_eaten_df$false_flee>0)
pred_flights <- subset(flights_eaten_df, flights_eaten_df$pred_flee>0)

# plotting: first, false alarms, THEN flight from predators, including death!
# false alarms
init_stack <- 3
fin_stack <- 5
sub <- head(false_flights)
for (i in 1:length(sub$time)){ #length(false_flights$time)
	dens <- sub$dens_false[i]
	false_flee <- sub$false_flee[i]	
	for (p in 11:(init_stack+10)){ # adding 10 to keep image file numbering in order
		png(paste("fflee_",i,"_",p,".png"))
		plot(rep(1,dens),ylim=c(1,3), pch=16, xaxt='n', yaxt='n', ann=FALSE)
		dev.off()
	}
	for (j in 1:(false_flee)){
		png(paste("fflee_",i,"_",j+init_stack+10,".png"))
		plot(c(rep(1,dens-j),rep(2,1),rep(3,(j-1))),ylim=c(1,3), pch=16, xaxt='n', yaxt='n', ann=FALSE)
		dev.off()
	}
	for (m in 1:(fin_stack)){
		png(paste("fflee_",i,"_",false_flee+init_stack+m+10,".png"))
		plot(c(rep(1,dens-false_flee), rep(3,false_flee)),ylim=c(1,3), pch=16, xaxt='n', yaxt='n', ann=FALSE)
		dev.off()
	}	
}

# one idea: randomly sample different densities and plot them, perhaps with densities as rows and generations as columns (perhaps f=1,5,10,20,etc)
# Could also animate more by having more than 3 positions (maybe 5?)
# Oct. 2: how I'll plot false_flights_dlist
# Could first subset time 1 from the first and last generation. It would be good to have some reps within a predation level to look for different strategies
# 

# show PREDATION: I could turn the color of the individual that dies (from black to red), I could even make it so that predation that happens for prey that detect the predator happens graphically during the flight (this would require that I have three or more positions to move individuals -- i.e., they flee from an initial position to an intermediate position and then a final position (an they turn red and stop in the intermediate position).
# I think the idea will be to show that before traits evolve, false alarms and predator attacks look similarly. BUT, perhaps with evolution, they will look very different (particularly at high densities). 
# To do this, I think I should randomly sample a large group size (maybe 25) from an early vs. late generation. Then I can compare what false alarms vs. predator attacks look like between these two. I would do this in conjunction with plotting the evolved traits: jumpiness, sociality, and density dependence in sociality
# When I come back to create the pred animation, I'll borrow from the code above

# update: I may want to create three types of animated plots: strategy 1, strategy 2, and random (from traits determined by random draws).

# To add mutation, I'd add a random number, drawn from a Gaussian distribution with mean 0 and variance gamma

# I could expand the model in a few obvious ways:
# Make sociality trait also raise prob of being in a larger group (this could put a scaling weight on group sizes sampled for assignment)
# make energy to fitness, and density dependence functions nonlinear. Energy to fitness should be asymptotic, could do the same for dd sociality.
# can use exponential decay and mirror of that for density dependence:
# (1 + or - exp(-s_dd*density))/2 - intercept 

