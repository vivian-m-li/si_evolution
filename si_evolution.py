import warnings
import numpy as np
import pandas as pd
from pandas.errors import SettingWithCopyWarning
import random
from analyze import write_output
from si_types import *
from constants import *
from plot import plot_traits
from collections import defaultdict
from typing import List, Tuple, DefaultDict, Set

warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)


def random_binomial(prob: float):
    return 1 if random.random() < prob else 0


def calc_stat(data) -> Stat:
    return Stat(mean=sum(data) / len(data), variance=np.var(data))


def assign_groups(
    indivs_alive: List[int], max_group_size: int
) -> Tuple[int, pd.DataFrame]:
    groups = []
    indivs = indivs_alive.copy()
    random.shuffle(indivs)
    while indivs:
        group_size = random.randint(1, min(max_group_size, len(indivs)))
        current_group = indivs[:group_size]
        groups.append(current_group)
        indivs = indivs[group_size:]
    group_lookup = {}
    for i, members in enumerate(groups):
        group_id = i + 1
        for member in members:
            group_lookup[member] = group_id

    group_vec = []
    for i in indivs_alive:
        group_vec.append(group_lookup[i])
    groups_df = pd.DataFrame({"individual": indivs_alive, "group_id": group_vec})
    num_groups = len(set(group_vec))
    return num_groups, groups_df


def init_outputs(params: Parameters) -> SimOutput:
    return SimOutput(
        parameters=params,
        total_deaths=[],
        false_flights=[],
        true_flights=[],
        detected_pred_deaths=[],
        nondetected_pred_deaths=[],
        trait_values=[],
        energetic_states=[],
        fitness=[],
        group_size=[],
    )


# Mutates the original output object
def init_outputs_per_generation(output: SimOutput) -> None:
    output.detected_pred_deaths.append(0)
    output.nondetected_pred_deaths.append(0)


def evo_fun(
    directory: str,
    params: Parameters,
    *,
    sim_id: int,
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
    fit = np.full((Ni, tf), np.nan)
    fit[:, 0] = 1

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
            f_false = coef_false * f_pred
            s_faith = np.random.uniform(0, 0.5, Ni)
            s_dd = np.random.uniform(-2, 2, Ni)
        else:
            f_pred = f_all[f - 1]["f_pred"].values
            f_false = coef_false * f_pred
            s_faith = f_all[f - 1]["s_faith"].values
            s_dd = f_all[f - 1]["s_dd"].values

        # For each time step, the population gets reassembled into groups based on a uniform distribution of group sizes. Then, each group is potentially subjected to a predator attack (based on a background predation level, probability set to 0.2 by default below).

        group_sizes: List[int] = []
        false_flights_by_group_size: DefaultDict[int, List[float]] = defaultdict(list)
        true_flights_by_group_size: DefaultDict[int, List[float]] = defaultdict(list)
        indivs_dead: Set[int] = set()
        for t in range(1, tf):
            indivs_alive = [i for i in list(range(1, Ni + 1)) if i not in indivs_dead]
            num_groups, groups_df = assign_groups(indivs_alive, max_group_size)
            group_sizes.append(Ni / num_groups)

            flights0 = np.full((num_groups, 4), np.nan)
            eaten_detect0 = np.full((num_groups, 2), np.nan)
            eaten_nodetect0 = np.full((num_groups, 2), np.nan)
            attacks_vec = np.full(num_groups, np.nan)

            for indiv_id in indivs_dead:
                fit[indiv_id - 1, t] = 0

            for group_idx in range(num_groups):
                pred = random_binomial(prob_pred)
                attacks_vec[group_idx] = pred
                prev_flee = 0
                subgroup = groups_df[groups_df["group_id"] == group_idx + 1]
                ddensity = len(subgroup)
                prev_detect = 0
                eaten_detect_vec = []
                eaten_nodetect_vec = []

                num_false_flee = 0
                num_true_flee = 0
                for i in range(len(subgroup)):
                    indiv_id = subgroup["individual"].values[i]
                    indiv_idx = indiv_id - 1
                    eaten_detect = 0
                    eaten_nodetect = 0
                    if prev_flee > 0:
                        p_flee_s = (
                            prev_flee * s_faith[indiv_idx]
                            + (ddensity - prev_flee) * s_dd[indiv_idx]
                        )
                        p_flee_s = min(1, p_flee_s)
                        p_flee_s = max(0, p_flee_s)
                        flee = random_binomial(f_false[indiv_idx]) or random_binomial(
                            p_flee_s
                        )
                    else:
                        flee = random_binomial(f_false[indiv_idx])
                    prev_flee += flee

                    if flee == 1:
                        fit[indiv_idx, t] = fit[indiv_idx, t - 1]
                        if pred == 1:  # true flight
                            num_true_flee += 1
                        else:  # false flight
                            num_false_flee += 1
                    elif pred == 1:
                        p_detect_s = (
                            prev_detect * s_faith[indiv_idx]
                            + (ddensity - prev_flee - prev_detect) * s_dd[indiv_idx]
                        )
                        p_detect_s = min(1, p_detect_s)
                        p_detect_s = max(0, p_detect_s)
                        detect = random_binomial(f_pred[indiv_idx]) or random_binomial(
                            p_detect_s
                        )
                        prev_detect += detect

                        if detect == 1:
                            peaten_detect = 1 / ((len(subgroup) - prev_flee) + 10)
                            eaten_detect = random_binomial(peaten_detect)
                            eaten_detect_vec.append(eaten_detect)
                        else:
                            peaten_nodetect = 1 / (len(subgroup) - prev_flee)
                            eaten_nodetect = random_binomial(peaten_nodetect)
                            eaten_nodetect_vec.append(eaten_nodetect)

                        if eaten_detect == 1 or eaten_nodetect == 1:
                            fit[indiv_idx, t] = 0
                            indivs_dead.add(indiv_id)

                        if eaten_detect == 1:
                            output.detected_pred_deaths[-1] += 1

                        elif eaten_nodetect == 1:
                            output.nondetected_pred_deaths[-1] += 1

                        else:
                            fit[indiv_idx, t] = fit[indiv_idx, t - 1]
                    else:
                        fit[indiv_idx, t] = fit[indiv_idx, t - 1] + e_gain

                flights0[group_idx, :] = [t, ddensity, prev_flee, prev_detect]
                eaten_detect0[group_idx, :] = [
                    t,
                    sum(eaten_detect_vec),
                ]
                eaten_nodetect0[group_idx, :] = [t, sum(eaten_nodetect_vec)]
                false_flights_by_group_size[ddensity].append(num_false_flee / ddensity)
                if pred:
                    true_flights_by_group_size[ddensity].append(
                        num_true_flee / ddensity
                    )

            attacks_all.append(attacks_vec)
            eaten_detect_all.append(eaten_detect0)
            eaten_nodetect_all.append(eaten_nodetect0)
            flights.append(flights0)

        output.total_deaths.append(len(indivs_dead))
        output.group_size.append(calc_stat(group_sizes))
        output.false_flights.append(
            {
                group_size: sum(freq_false_flights) / len(freq_false_flights)
                for group_size, freq_false_flights in false_flights_by_group_size.items()
            }
        )
        output.true_flights.append(
            {
                group_size: sum(freq_true_flights) / len(freq_true_flights)
                for group_size, freq_true_flights in true_flights_by_group_size.items()
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
        f_all.append(survive.iloc[f_index, 1:4])

        for lst, op in [(trait_mean, np.mean), (trait_sd, np.std), (trait_var, np.var)]:
            lst[f, :] = [
                op(survive["f_pred"]),
                op(survive["s_faith"]),
                op(survive["s_dd"]),
            ]
        output.trait_values.append(
            [
                Stat(mean=trait_mean[f][trait_idx], variance=trait_var[f, trait_idx])
                for trait_idx in range(3)
            ]
        )

    if save_output:
        write_output(directory, sim_id, output)

    if plot:
        plot_traits(trait_mean, trait_sd, fit_traits_gen)
