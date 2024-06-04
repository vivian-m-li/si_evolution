import numpy as np
import pandas as pd
import random
from typing import List, Tuple
from si_types import *


def calc_mean(data: List[float]) -> List[float]:
    return np.mean(np.array(data))


def random_binomial(prob: float):
    return 1 if random.random() < prob else 0


def calc_stat(data) -> Stat:
    return Stat(mean=calc_mean(data), variance=np.var(data))


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
        all_group_sizes=[],
        pred_catch_rate=[],
        pred_catch_by_group_size=[],
        prop_groups_attacked=[],
    )


# Mutates the original output object
def init_outputs_per_generation(output: SimOutput) -> None:
    output.detected_pred_deaths.append(0)
    output.nondetected_pred_deaths.append(0)
    output.prop_groups_attacked.append([])


def build_output_path(params: Parameters):
    return f"{params.Ni}_{params.tf}_{params.e_gain}_{params.coef_false}_{params.maxf}_{params.prob_pred}_{params.max_group_size}"
