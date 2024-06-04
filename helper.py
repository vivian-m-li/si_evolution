import os
import csv
import ast
import uuid
from scipy import stats
import numpy as np
import pandas as pd
import random
from typing import List, Tuple, Any
from dataclasses import fields, asdict
from si_types import *

CONFIDENCE_LEVEL = 0.95


def calc_mean(data: List[float]) -> List[float]:
    return np.mean(np.array(data))


def random_binomial(prob: float):
    return 1 if random.random() < prob else 0


def calc_stat(data) -> Stat:
    return Stat(mean=calc_mean(data), variance=np.var(data))


def eval_int_float(x):
    if type(x) != str:
        return x
    if "." not in x or x.endswith(".0"):
        return int(float(x))
    if "." in x:
        return float(x)
    return x


def cast_data_types(row: List[str]) -> List[Any]:
    data = []
    for x in row:
        try:
            value = ast.literal_eval(x)
            if isinstance(value, dict):
                data.append(value)
            elif isinstance(value, list):
                data.append(value)
            else:
                data.append(eval_int_float(x))
        except (ValueError, SyntaxError):
            data.append(eval_int_float(x))
    return data


def cast_dataclass_types(instance: Any) -> Any:
    for field in fields(instance):
        value = getattr(instance, field.name)
        casted_value = field.type(value)
        setattr(instance, field.name, casted_value)
    return instance


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


def cast_param_dataclass(output_params: OutputParameters) -> Parameters:
    b_fields = asdict(output_params)
    a_fields = {field.name: b_fields[field.name] for field in fields(Parameters)}
    return Parameters(**a_fields)


def build_output_path(params: Parameters):
    params = cast_dataclass_types(params)
    param_lst = [
        params.Ni,
        params.tf,
        params.e_gain,
        params.coef_false,
        params.maxf,
        params.prob_pred,
        params.max_group_size,
    ]
    param_lst = [str(x) for x in param_lst]
    return "_".join(param_lst)


def calc_confidence_interval(means: List[float]) -> Tuple[float, float]:
    mean_of_means = calc_mean(means)
    std_dev_of_sample_means = np.std(means, ddof=1) / np.sqrt(len(means))
    z_critical = stats.norm.ppf((1 + CONFIDENCE_LEVEL) / 2)
    margin_of_error = z_critical * std_dev_of_sample_means
    confidence_interval = (
        mean_of_means - margin_of_error,
        mean_of_means + margin_of_error,
    )
    return confidence_interval


def get_sim_id(file_name: str) -> int:
    if not os.path.isfile(file_name):
        return 0

    df = pd.read_csv(file_name)
    df.dropna(axis=0, how="all", inplace=True)
    df.to_csv(file_name, index=False)

    sim_id = 0
    if not df.empty:
        sim_id = int(df.iloc[:, 0].max()) + 1

    return sim_id


def write_new_all_file(file_name):
    results_file = open(file_name, "w")
    writer_object = csv.writer(results_file, lineterminator="\n")
    writer_object.writerow(
        [
            "id",
            "Ni",
            "tf",
            "e_gain",
            "coef_false",
            "maxf",
            "prob_pred",
            "max_group_size",
        ]
    )
    results_file.close()


def write_sim_file(directory: str, sim_id: str, output: SimOutput):
    sim_file = open(f"{directory}/{sim_id}.csv", "a")
    writer_object = csv.writer(sim_file, lineterminator="\n")
    headers = [
        "generation",
        "total_deaths",
        "freq_false_flights",
        "freq_true_flights",
        "freq_detected_pred_deaths",
        "freq_nondetected_pred_deaths",
        "group_size_mean",
        "group_size_var",
        "all_group_sizes",
        "energetic_states_mean",
        "energetic_states_var",
        "fitness_mean",
        "fitness_var",
        "f_pred_mean",
        "f_pred_var",
        "s_faith_mean",
        "s_faith_var",
        "s_dd_mean",
        "s_dd_var",
        "pred_catch_rate",
        "pred_catch_by_group_size",
        "prop_groups_attacked",
    ]
    writer_object.writerow(headers)

    sim_results = []
    for g in range(output.parameters.maxf):
        sim_results.append(
            [
                g + 1,
                output.total_deaths[g],
                str(output.false_flights[g]),
                str(output.true_flights[g]),
                output.detected_pred_deaths[g] / max(1, output.total_deaths[g]),
                output.nondetected_pred_deaths[g] / max(1, output.total_deaths[g]),
                output.group_size[g].mean,
                output.group_size[g].variance,
                str(output.all_group_sizes[g]),
                output.energetic_states[g].mean,
                output.energetic_states[g].variance,
                output.fitness[g].mean,
                output.fitness[g].variance,
                output.trait_values[g][0].mean,
                output.trait_values[g][0].variance,
                output.trait_values[g][1].mean,
                output.trait_values[g][1].variance,
                output.trait_values[g][2].mean,
                output.trait_values[g][2].variance,
                output.pred_catch_rate[g],
                str(output.pred_catch_by_group_size[g]),
                str(output.prop_groups_attacked[g]),
            ]
        )
    writer_object.writerows(sim_results)

    sim_file.close()


def write_output(directory: str, output: SimOutput):
    if not os.path.exists(directory):
        os.makedirs(directory)
    output_path = build_output_path(output.parameters)
    output_dir = f"{directory}/{output_path}"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    unique_sim_id = uuid.uuid4()
    write_sim_file(output_dir, unique_sim_id, output)


def write_output_old(directory: str, sim_id: int, output: SimOutput):
    results_file_name = f"{directory}/all.csv"
    if not os.path.exists(directory):
        os.makedirs(directory)
        write_new_all_file(results_file_name)

    results_file = open(results_file_name, "a")
    writer_object = csv.writer(results_file, lineterminator="\n")
    writer_object.writerow(
        [
            sim_id,
            output.parameters.Ni,
            output.parameters.tf,
            output.parameters.e_gain,
            output.parameters.coef_false,
            output.parameters.maxf,
            output.parameters.prob_pred,
            output.parameters.max_group_size,
        ]
    )
    results_file.close()

    write_sim_file(directory, sim_id, output)


def get_all_outputs(
    out_file_path: str, all_params: List[OutputParameters]
) -> List[List[pd.DataFrame]]:
    sims: List[List[pd.DataFrame]] = [[] for _ in range(len(all_params))]
    for i, params in enumerate(all_params):
        output_dir = (
            f"{out_file_path}/{build_output_path(cast_param_dataclass(params))}"
        )
        try:
            files = os.listdir(output_dir)
        except Exception:
            continue
        for file_name in files:
            df = pd.read_csv(f"{output_dir}/{file_name}").to_numpy()
            sims[i].append(df)
    return sims


def get_all_outputs_old(
    out_file_path: str, all_params: List[OutputParameters]
) -> List[List[pd.DataFrame]]:
    sims: List[List[pd.DataFrame]] = [[] for _ in range(len(all_params))]
    all_sims_file = open(f"{out_file_path}/all.csv", "r")
    reader_object = csv.reader(all_sims_file, delimiter=",")
    next(reader_object)
    for row in reader_object:
        (
            sim_id,
            Ni,
            tf,
            e_gain,
            coef_false,
            maxf,
            prob_pred,
            max_group_size,
        ) = cast_data_types(row)
        for i, params in enumerate(all_params):
            if (
                params.Ni == Ni
                and params.tf == tf
                and params.e_gain == e_gain
                and params.coef_false == coef_false
                and params.maxf == maxf
                and params.prob_pred == prob_pred
                and params.max_group_size == max_group_size
            ):
                df = pd.read_csv(f"{out_file_path}/{sim_id}.csv").to_numpy()
                sims[i].append(df)

    all_sims_file.close()
    return sims
