import warnings
import numpy as np
import pandas as pd
from pandas.errors import SettingWithCopyWarning
import random
from analyze import write_output
from si_types import *
from constants import *
from helper import *
from plot import plot_traits
from collections import defaultdict
from typing import List, DefaultDict, Set

warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
pd.set_option("display.max_rows", None)


def evo_fun(
    directory: str,
    params: Parameters,
    *,
    save_output: bool = True,
    plot: bool = False,
):
    Ni = params.Ni
    tf = params.tf
    e_gain = params.e_gain
    coef_false = params.coef_false
    maxf = params.maxf
    prob_pred = params.prob_pred
    max_group_size = params.max_group_size

    # Initialize trait values
    f_pred = np.random.uniform(0, 1, Ni)
    s_faith = np.random.uniform(0, 0.5, Ni)
    s_dd = np.random.uniform(-2, 2, Ni)

    # Initialize data structures for storing outputs per timestep/generation
    attacks_all = []
    fit_traits_gen0 = []
    fit_traits_gen = []
    trait_mean = np.full((maxf, 3), np.nan)
    trait_sd = np.full((maxf, 3), np.nan)
    trait_var = np.full((maxf, 3), np.nan)
    f_all = []
    flights_master = []
    eaten_detect_master = []
    eaten_nodetect_master = []

    output = init_outputs(params)

    for f in range(maxf):
        init_outputs_per_generation(output)
        flights = []
        eaten_detect_all = []
        eaten_nodetect_all = []

        if f == 0:
            f_pred = np.random.uniform(0, 1, Ni)
            s_faith = np.random.uniform(0, 0.5, Ni)
            s_dd = np.random.uniform(-2, 2, Ni)
        else:
            f_pred = f_all[-1]["f_pred"].values
            s_faith = f_all[-1]["s_faith"].values
            s_dd = f_all[-1]["s_dd"].values
        f_false = coef_false * f_pred
        fit = np.full((Ni, tf), np.nan)
        fit[:, 0] = 1

        # For each time step, the population gets reassembled into groups based on a uniform distribution of group sizes. Then, each group is potentially subjected to a predator attack (based on a background predation level, probability set to 0.2 by default below).
        group_sizes: List[int] = []
        false_flights_by_group_size: DefaultDict[int, List[float]] = defaultdict(list)
        true_flights_by_group_size: DefaultDict[int, List[float]] = defaultdict(list)
        indivs_dead: Set[int] = set()
        num_pred_attacks = 0
        num_pred_by_group_size: DefaultDict[int, int] = defaultdict(lambda: 0)
        num_caught_by_group_size: DefaultDict[int, int] = defaultdict(lambda: 0)
        for t in range(1, tf):
            indivs_alive = [i for i in list(range(1, Ni + 1)) if i not in indivs_dead]
            num_groups, groups_df = assign_groups(indivs_alive, max_group_size)

            flights0 = np.full((num_groups, 4), np.nan)
            eaten_detect0 = np.full((num_groups, 2), np.nan)
            eaten_nodetect0 = np.full((num_groups, 2), np.nan)
            attacks_vec = np.full(num_groups, np.nan)

            for indiv_id in indivs_dead:
                fit[indiv_id - 1, t] = 0

            num_groups_attacked = min(num_groups, np.random.poisson(prob_pred))
            output.prop_groups_attacked[-1].append(num_groups_attacked / num_groups)
            subgroups = list(range(num_groups))
            groups_attacked = random.choices(subgroups, k=num_groups_attacked)
            for group_idx in subgroups:
                pred = group_idx in groups_attacked
                num_pred_attacks += int(pred)
                subgroup = groups_df[groups_df["group_id"] == group_idx + 1]
                ddensity = len(subgroup)
                group_sizes.append(ddensity)
                attacks_vec[group_idx] = pred
                prev_flee = 0  # false flights
                prev_detect = 0  # true flights
                eaten_detect_vec = []
                eaten_nodetect_vec = []
                poss_deaths_weighted: Dict[int, float] = {}  # indiv_id: prob_death
                detected: Dict[int, bool] = {}  # indiv_id: detected

                for i in range(len(subgroup)):
                    indiv_id = subgroup["individual"].values[i]
                    indiv_idx = indiv_id - 1
                    eaten_detect = 0
                    eaten_nodetect = 0
                    if prev_flee > 0:
                        p_flee_s = (
                            prev_flee * s_faith[indiv_idx]
                            + (ddensity - prev_flee) / 25 * s_dd[indiv_idx]
                        )
                        p_flee_s = min(1, p_flee_s)
                        p_flee_s = max(0, p_flee_s)
                        flee = random_binomial(f_false[indiv_idx]) or random_binomial(
                            p_flee_s
                        )
                    else:
                        flee = random_binomial(f_false[indiv_idx])
                    prev_flee += flee

                    if flee or pred:
                        fit[indiv_idx, t] = fit[indiv_idx, t - 1]

                    if pred and not flee:
                        p_detect_s = (
                            prev_detect * s_faith[indiv_idx]
                            + (ddensity - prev_flee - prev_detect)
                            / 25
                            * s_dd[indiv_idx]
                        )
                        p_detect_s = min(1, p_detect_s)
                        p_detect_s = max(0, p_detect_s)
                        detect = random_binomial(f_pred[indiv_idx]) or random_binomial(
                            p_detect_s
                        )
                        prev_detect += detect
                        detected[indiv_id] = detect

                        if detect:
                            peaten_detect = 1 / ((len(subgroup) - prev_flee) + 10)
                            poss_deaths_weighted[indiv_id] = peaten_detect
                        else:
                            peaten_nodetect = 1 / (len(subgroup) - prev_flee)
                            poss_deaths_weighted[indiv_id] = peaten_nodetect

                    if not pred and not flee:
                        fit[indiv_idx, t] = fit[indiv_idx, t - 1] + e_gain

                indiv_eaten = False
                sorted_poss_deaths = dict(
                    sorted(
                        poss_deaths_weighted.items(),
                        key=lambda item: item[1],
                        reverse=True,
                    )
                )
                for indiv_id, p_eaten in sorted_poss_deaths.items():
                    indiv_idx = indiv_id - 1
                    if not indiv_eaten:
                        eaten = random_binomial(p_eaten)
                        if eaten:
                            indiv_eaten = True
                            fit[indiv_idx, t] = 0
                            indivs_dead.add(indiv_id)
                            if detected[indiv_id]:
                                output.detected_pred_deaths[-1] += 1
                            else:
                                output.nondetected_pred_deaths[-1] += 1
                    fit[indiv_idx, t] = fit[indiv_idx, t - 1]

                num_pred_by_group_size[ddensity] += pred
                num_caught_by_group_size[ddensity] += int(indiv_eaten)

                flights0[group_idx, :] = [t, ddensity, prev_flee, prev_detect]
                eaten_detect0[group_idx, :] = [
                    t,
                    sum(eaten_detect_vec),
                ]
                eaten_nodetect0[group_idx, :] = [t, sum(eaten_nodetect_vec)]
                if prev_flee > 0:
                    false_flights_by_group_size[ddensity].append(prev_flee / ddensity)
                if prev_detect > 0:
                    true_flights_by_group_size[ddensity].append(
                        prev_detect / (ddensity - prev_flee)
                    )

            attacks_all.append(attacks_vec)
            eaten_detect_all.append(eaten_detect0)
            eaten_nodetect_all.append(eaten_nodetect0)
            flights.append(flights0)

        output.total_deaths.append(len(indivs_dead))
        output.group_size.append(calc_stat(group_sizes))
        output.all_group_sizes.append(group_sizes)
        output.false_flights.append(
            {
                group_size: calc_mean(freq_false_flights)
                for group_size, freq_false_flights in false_flights_by_group_size.items()
            }
        )
        output.true_flights.append(
            {
                group_size: calc_mean(freq_true_flights)
                for group_size, freq_true_flights in true_flights_by_group_size.items()
            }
        )
        output.pred_catch_rate.append(
            len(indivs_dead) / max(num_pred_attacks, 1)
        )  # TODO: save None if num_pred_attacks is none

        output.pred_catch_by_group_size.append(
            {
                group_size: num_caught_by_group_size.get(group_size, 0) / num_preds
                for group_size, num_preds in num_pred_by_group_size.items()
                if num_preds > 0
            }
        )
        flights_master.append(np.concatenate(flights))
        eaten_detect_master.append(np.concatenate(eaten_detect_all))
        eaten_nodetect_master.append(np.concatenate(eaten_nodetect_all))

        survive0 = pd.DataFrame(
            {"fit": fit[:, tf - 1], "f_pred": f_pred, "s_faith": s_faith, "s_dd": s_dd}
        )
        fit_traits_gen0.append(survive0)
        survive = survive0[survive0["fit"] > 0]
        fit_traits_gen.append(survive.copy())

        energetic_states = calc_stat(survive0.iloc[:, 0].to_numpy())
        output.energetic_states.append(energetic_states)
        fitness_stat = calc_stat(survive.iloc[:, 0].to_numpy())
        output.fitness.append(fitness_stat)

        survive["fit"] = survive["fit"] / np.sum(survive["fit"])

        surv_df = pd.DataFrame(survive)
        surv_df["index"] = range(len(surv_df["fit"]))
        f_index = np.random.choice(surv_df["index"], Ni, p=surv_df["fit"])
        gen_final_traits = surv_df.iloc[f_index, 1:4].copy()

        # add a little bit of noise as mutation
        for trait in gen_final_traits.columns:
            sd = 0.02 * abs(np.mean(gen_final_traits[trait]))
            noise = np.random.normal(0, sd, gen_final_traits[trait].shape)
            gen_final_traits[trait] += noise
        f_all.append(gen_final_traits)

        for lst, op in [(trait_mean, np.mean), (trait_sd, np.std), (trait_var, np.var)]:
            lst[f, :] = [
                op(surv_df["f_pred"]),
                op(surv_df["s_faith"]),
                op(surv_df["s_dd"]),
            ]
        output.trait_values.append(
            [
                Stat(mean=trait_mean[f][trait_idx], variance=trait_var[f, trait_idx])
                for trait_idx in range(3)
            ]
        )

    if save_output:
        write_output(directory, output)

    if plot:
        plot_traits(trait_mean, trait_sd, fit_traits_gen)
