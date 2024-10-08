# si_evolution
This project analyzes the evolution of social behavior in foraging populations using an agent-based model.

# Running the simulations
1. Fill in `sim_params` in run_simulations.py
2. Run analyze.ipynb to plot

# Changes made to si_evolution.R
- Parallelized the simulations to run on a computing cluster
- Changed prob_pred (previously: probability each group is attacked) to lambda (now: avg # of groups attacked according to a poisson distribution)
  - `num_groups_attacked = min(num_groups, np.random.poisson(prob_pred))`
  - Using prob_pred biased the number of groups attacked based on group size
- Scaled the effects of s_faith and s_dd (previously, the effects of s_dd were dominating the effects of s_faith) on `p_flee_s` and `p_detect_s`:
   - ex: `p_flee_s = prev_flee * s_faith[indiv_idx] + (ddensity - prev_flee) / 25 * s_dd[indiv_idx]`
   - Added the `/25` but this scaling factor is still affected by group size
- Changed some of the code on keeping track of `eaten_detect/eaten_nodetect`
- Added noise to simulate mutation of trait values
- Kept track of metrics at each timestep/generation and saved in an output file for each simulation
- New analysis code (looks through thousands of relevant simulation results) + plots in `analyze.py` and `plot.py`

# Future Work
1. Mutation rate - If we run our simulation with predation = 0, then we expect to see jumpiness decrease to 0, which it never does so fully. Alternatively, if we run our simulation with predation = 0 and initial trait values = 0, then we also expect the trait values to remain at 0. Without mutation, we don't know whether traits are reaching their optimal values or conveniently reaching a steady state value based on the trait values in the population. The mutation rate would have to be relatively high, or we would have to increase the number of generations that the simulation is run for.
2. Rescale the effects of s_dd and s_faith - Currently, if you run the simulations and plot mean trait values by average group size, we see that sociality increases with decreasing average/max group size. This is likely an artifact of the code when we calculate p_flee and p_detect (lines 174 and 191) by dividing the effect of s_dd by max_group_size. We should divide by a constant (like 10 or 25) instead, aiming for s_faith and s_dd to have equal effects on p_flee and p_detect.
3. Social learning - Allowing individuals to reach more optimal trait values during a generation would be an interesting avenue to explore.
4. Competition - We're currently assuming that there are infinite resources, which can be thought of as true for certain ecosystems, such as herbivorous reef fish that feed on algae. However, implementing competition would force individuals to form smaller groups that are still large enough so as to not diminish their chances of survival. 
5. Leadership - In many species, even in fission-fusion societies, social hierarchy exists where some individuals may have more influence over the decisions of others. We could add a new trait for social influence and explore how the trait impacts group dynamics and decision making.
