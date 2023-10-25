import warnings
import numpy as np
import pandas as pd
from pandas.errors import SettingWithCopyWarning
import random
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from typing import List, Dict, Tuple

warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)


def assign_groups(Ni: int, max_group_size: int) -> Tuple[int, pd.DataFrame]:
    groups = []
    indivs = list(range(1, Ni + 1))
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

    indivs = list(range(1, Ni + 1))
    group_vec = []
    for i in indivs:
        group_vec.append(group_lookup[i])

    groups_df = pd.DataFrame({"individual": indivs, "group_id": group_vec})
    num_groups = len(set(group_vec))
    return num_groups, groups_df


def evo_fun(
    *,
    Ni: int = 100,  # number of individuals
    tf: int = 30,  # time steps for a given generation
    e_gain: int = 1,  # energy units gained each time step in which the ith individual does not flee
    coef_false: float = 0.2,  # coefficient that determines the prob of a false alarm (a smaller value than f_pred)
    maxf: int = 100,  # number of generations to run the model through
    prob_pred: float = 0.2,
    max_group_size: int = 25,
):
    f_pred = np.random.uniform(0, 1, Ni)
    s_faith = np.random.uniform(0, 0.5, Ni)
    s_dd = np.random.uniform(-2, 2, Ni)
    fit = np.full((Ni, tf), np.nan)
    fit[:, 0] = 1

    attacks_all = []
    fit_traits_gen0 = []
    fit_traits_gen = []
    trait_mean = np.full((maxf, 3), np.nan)
    trait_sd = np.full((maxf, 3), np.nan)
    f_all = []
    flights_master = []
    eaten_detect_master = []
    eaten_nodetect_master = []

    for f in range(maxf):
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

        for t in range(1, tf):
            num_groups, groups_df = assign_groups(Ni, max_group_size)

            flights0 = np.full((num_groups, 4), np.nan)
            eaten_detect0 = np.full((num_groups, 2), np.nan)
            eaten_nodetect0 = np.full((num_groups, 2), np.nan)
            attacks_vec = np.full(num_groups, np.nan)

            for group_id in range(num_groups):
                pred = np.random.binomial(1, prob_pred)
                attacks_vec[group_id] = pred
                prev_flee = 0
                subgroup = groups_df[groups_df["group_id"] == group_id + 1]
                ddensity = len(subgroup)
                prev_detect = 0
                eaten_detect_vec = []
                eaten_nodetect_vec = []

                for i in range(len(subgroup)):
                    ii = subgroup["individual"].values[i] - 1
                    eaten_detect = 0
                    eaten_nodetect = 0
                    if fit[ii, t - 1] == 0:
                        fit[ii, t] = 0
                    else:
                        if prev_flee > 0:
                            p_flee_s = (
                                prev_flee * s_faith[ii]
                                + (ddensity - prev_flee) * s_dd[ii]
                            )
                            p_flee_s = min(1, p_flee_s)
                            p_flee_s = max(0, p_flee_s)
                            flee0 = np.random.binomial(
                                1, f_false[ii]
                            ) + np.random.binomial(1, p_flee_s)
                            flee = int(flee0 >= 1)
                        else:
                            flee = np.random.binomial(1, f_false[ii])

                        prev_flee += flee

                        if flee == 1:
                            fit[ii, t] = fit[ii, t - 1]
                        elif pred == 1:
                            p_detect_s = (
                                prev_detect * s_faith[ii]
                                + (ddensity - prev_flee - prev_detect) * s_dd[ii]
                            )
                            p_detect_s = min(1, p_detect_s)
                            p_detect_s = max(0, p_detect_s)
                            detect0 = np.random.binomial(
                                1, f_pred[ii]
                            ) + np.random.binomial(1, p_detect_s)
                            detect = int(detect0 >= 1)
                            prev_detect += detect

                            if detect == 1:
                                peaten_detect = 1 / ((len(subgroup) - prev_flee) + 10)
                                eaten_detect = np.random.binomial(1, peaten_detect)
                                eaten_detect_vec.append(eaten_detect)
                            else:
                                peaten_nodetect = 1 / (len(subgroup) - prev_flee)
                                eaten_nodetect = np.random.binomial(1, peaten_nodetect)
                                eaten_nodetect_vec.append(eaten_nodetect)

                            if eaten_detect == 1 or eaten_nodetect == 1:
                                fit[ii, t] = 0
                            else:
                                fit[ii, t] = fit[ii, t - 1]
                        else:
                            fit[ii, t] = fit[ii, t - 1] + e_gain

                flights0[group_id, :] = [t, ddensity, prev_flee, prev_detect]
                eaten_detect0[group_id, :] = [t, len(eaten_detect_vec)]
                eaten_nodetect0[group_id, :] = [t, len(eaten_nodetect_vec)]

            attacks_all.append(attacks_vec)
            eaten_detect_all.append(eaten_detect0)
            eaten_nodetect_all.append(eaten_nodetect0)
            flights.append(flights0)

        flights_master.append(np.concatenate(flights))
        eaten_detect_master.append(np.concatenate(eaten_detect_all))
        eaten_nodetect_master.append(np.concatenate(eaten_nodetect_all))

        survive0 = pd.DataFrame(
            {"fit": fit[:, tf - 1], "f_pred": f_pred, "s_faith": s_faith, "s_dd": s_dd}
        )
        fit_traits_gen0.append(survive0)
        survive = survive0[survive0["fit"] > 0]
        fit_traits_gen.append(survive.copy())
        trait_mean[f, :] = [
            np.mean(survive["f_pred"]),
            np.mean(survive["s_faith"]),
            np.mean(survive["s_dd"]),
        ]
        trait_sd[f, :] = [
            np.std(survive["f_pred"]),
            np.std(survive["s_faith"]),
            np.std(survive["s_dd"]),
        ]

        survive["fit"] = survive["fit"] / np.sum(survive["fit"])

        surv_df = pd.DataFrame(survive)
        surv_df["index"] = range(len(surv_df["fit"]))
        f_index = np.random.choice(surv_df["index"], Ni, p=surv_df["fit"])
        f_all.append(survive.iloc[f_index, 1:4])

    traits = ["jumpiness", "sociality", "density dependence in sociality"]
    out_ls = []

    fig = plt.figure(figsize=(6, 10))
    gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 1])
    for i in range(len(traits)):
        dff = pd.DataFrame(
            {
                "generation": range(1, len(trait_mean) + 1),
                "trait_mean": trait_mean[:, i],
                "lb": trait_mean[:, i] - trait_sd[:, i],
                "ub": trait_mean[:, i] + trait_sd[:, i],
            }
        )
        out_ls.append(dff)

        ax = plt.subplot(gs[i, 0])
        ax.plot(dff["generation"], dff["trait_mean"], color="black", label=traits[i])
        ax.fill_between(dff["generation"], dff["lb"], dff["ub"], alpha=0.2)
        ax.set_ylabel(traits[i])

    mean_ff = round(np.mean(fit_traits_gen[-1]["fit"]), 2)
    sd_ff = round(np.std(fit_traits_gen[-1]["fit"]), 2)

    plt.suptitle(f"Fitness = {mean_ff} +/- {sd_ff} SD")
    plt.xlabel("Generation")
    plt.tight_layout()
    plt.show()

    return out_ls
