source("si_evolution.R")


run_sims <- function() {
  maxf0 <- 500
  evo.fun1(maxf = maxf0)

  # max_groups <- c(5, 10, 20)

  # g1 <- evo.fun1(max_group_size = max_groups[1], maxf = maxf0)
  # g2 <- evo.fun1(max_group_size = max_groups[2], maxf = maxf0)
  # g3 <- evo.fun1(max_group_size = max_groups[3], maxf = maxf0)
  #  
  # sub1 <- data.frame(max_group_size = max_groups, soc_mean = c(mean(g1[[2]][which(g1[[2]][ , "generation"] == maxf), "trait_mean"]), mean(g2[[2]][which(g2[[2]][ , "generation"] == maxf), "trait_mean"]), mean(g3[[2]][which(g3[[2]][ , "generation"] == maxf), "trait_mean"])),
  # 	soc_lb = c(mean(g1[[2]][which(g1[[2]][ , "generation"] == maxf), "lb"]), mean(g2[[2]][which(g2[[2]][ , "generation"] == maxf), "lb"]), mean(g3[[2]][which(g3[[2]][ , "generation"] == maxf), "lb"])),
  # 	soc_ub = c(mean(g1[[2]][which(g1[[2]][ , "generation"] == maxf), "ub"]), mean(g2[[2]][which(g2[[2]][ , "generation"] == maxf), "ub"]), mean(g3[[2]][which(g3[[2]][ , "generation"] == maxf), "ub"])))
  # 	
  # dev.new();
  # with(sub1, plot(soc_mean ~ max_group_size, ylim = c(min(sub1$soc_lb), max(sub1$soc_ub))))
  # with(sub1, points(soc_lb ~ max_group_size, pch = "-"))
  # with(sub1, points(soc_ub ~ max_group_size, pch = "-"))



  # ptm <- proc.time()
  # 
  # lopred1 <- evo.fun1(prob_pred = 0.2, maxf = 1000)
  # lopred2 <- evo.fun1(prob_pred = 0.2, maxf = 1000)
  # lopred3 <- evo.fun1(prob_pred = 0.2, maxf = 1000)
  # lopred4 <- evo.fun1(prob_pred = 0.2, maxf = 1000)
  # lopred5 <- evo.fun1(prob_pred = 0.2, maxf = 1000)
  # 
  # mdpred1 <- evo.fun1(prob_pred = 0.4, maxf = 1000)
  # mdpred2 <- evo.fun1(prob_pred = 0.4, maxf = 1000)
  # mdpred3 <- evo.fun1(prob_pred = 0.4, maxf = 1000)
  # mdpred4 <- evo.fun1(prob_pred = 0.4, maxf = 1000)
  # mdpred5 <- evo.fun1(prob_pred = 0.4, maxf = 1000)
  # 
  # proc.time() - ptm

}

run_sims()
