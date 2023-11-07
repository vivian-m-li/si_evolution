import os
import csv
import pandas as pd
import numpy as np
import ast
from scipy import stats
from plot import *
from si_types import *
from constants import *
from collections import defaultdict
from typing import Optional, List, DefaultDict, Any

CONFIDENCE_LEVEL = 0.95


def calc_confidence_interval(means: List[float]) -> Tuple[float, float]:
    means = np.array(means)
    mean_of_means = np.mean(means)
    std_dev_of_sample_means = np.std(means, ddof=1) / np.sqrt(len(means))
    z_critical = stats.norm.ppf((1 + CONFIDENCE_LEVEL) / 2)
    margin_of_error = z_critical * std_dev_of_sample_means
    confidence_interval = (
        mean_of_means - margin_of_error,
        mean_of_means + margin_of_error,
    )
    return confidence_interval


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
            else:
                data.append(eval_int_float(x))
        except (ValueError, SyntaxError):
            try:
                data.append(eval_int_float(x))
            except Exception:
                import pdb

                pdb.set_trace()
    return data


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


def write_output(directory: str, sim_id: int, output: SimOutput):
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
                output.detected_pred_deaths[g] / output.total_deaths[g],
                output.nondetected_pred_deaths[g] / output.total_deaths[g],
                output.group_size[g].mean,
                output.group_size[g].variance,
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
            ]
        )
    writer_object.writerows(sim_results)

    sim_file.close()


def get_outputs(
    out_file_path: str, params: Optional[OutputParameters]
) -> List[pd.DataFrame]:
    sims = []

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
        if params is not None and (
            params.Ni != Ni
            or params.tf != tf
            or params.e_gain != e_gain
            or params.coef_false != coef_false
            or params.maxf != maxf
            or params.prob_pred != prob_pred
            or params.max_group_size != max_group_size
        ):
            continue
        df = pd.read_csv(f"{out_file_path}/{sim_id}.csv").to_numpy()
        sims.append(df)

    all_sims_file.close()

    return sims


def process_results(
    sim_outputs: List[pd.DataFrame], params: OutputParameters
) -> Results:
    num_generations = params.maxf

    freq_false_flight_by_group_size: List[DefaultDict[int, List[float]]] = []
    freq_true_flight_by_group_size: List[DefaultDict[int, List[float]]] = []
    fitness_stat: List[Stat] = []
    trait_values: List[List[Stat]] = [[], [], []]
    for i in range(num_generations):
        freq_false_flight_by_group_size.append(defaultdict(list))
        freq_true_flight_by_group_size.append(defaultdict(list))
        group_stats: Dict[str, GroupStats] = {
            "fitness": GroupStats(means=[], vars=[]),
            "f_pred": GroupStats(means=[], vars=[]),
            "s_faith": GroupStats(means=[], vars=[]),
            "s_dd": GroupStats(means=[], vars=[]),
        }
        for sim in sim_outputs:
            (
                gen,
                total_deaths,
                freq_false_flights,
                freq_true_flights,
                freq_detected_pred_deaths,
                freq_nondetected_pred_deaths,
                group_size_mean,
                group_size_var,
                energetic_states_mean,
                energetic_states_var,
                fitness_mean,
                fitness_var,
                f_pred_mean,
                f_pred_var,
                s_faith_mean,
                s_faith_var,
                s_dd_mean,
                s_dd_var,
            ) = cast_data_types(sim[i])
            for key, mean, var in [
                ["fitness", fitness_mean, fitness_var],
                ["f_pred", f_pred_mean, f_pred_var],
                ["s_faith", s_faith_mean, s_faith_var],
                ["s_dd", s_dd_mean, s_dd_var],
            ]:
                group_stats[key].means.append(mean)
                group_stats[key].vars.append(var)
            for group_size, freq_false_flight in freq_false_flights.items():
                freq_false_flight_by_group_size[-1][group_size].append(
                    freq_false_flight
                )
            for group_size, freq_true_flight in freq_true_flights.items():
                freq_true_flight_by_group_size[-1][group_size].append(freq_true_flight)
        fitness_stat.append(
            Stat(
                mean=sum(group_stats["fitness"].means)
                / len(group_stats["fitness"].means),
                variance=np.var(group_stats["fitness"].vars),
                confidence_interval=calc_confidence_interval(
                    group_stats["fitness"].means
                ),
            )
        )
        for j, trait in enumerate(["f_pred", "s_faith", "s_dd"]):
            trait_values[j].append(
                Stat(
                    mean=sum(group_stats[trait].means) / len(group_stats[trait].means),
                    variance=np.var(group_stats[trait].vars),
                    confidence_interval=calc_confidence_interval(
                        group_stats[trait].means
                    ),
                )
            )

    freq_false_flights_unbinned: List[float] = []
    freq_true_flights_unbinned: List[float] = []
    freq_false_flights_binned: List[List[Optional[float]]] = []
    freq_true_flights_binned: List[List[Optional[float]]] = []
    for i in range(num_generations):
        freq_false_flights_binned.append([])
        freq_true_flights_binned.append([])
        all_false_flights_by_group = freq_false_flight_by_group_size[i]
        all_true_flights_by_group = freq_true_flight_by_group_size[i]

        all_false_flights = [x for y in all_false_flights_by_group.values() for x in y]
        all_true_flights = [x for y in all_true_flights_by_group.values() for x in y]

        freq_false_flights_unbinned.append(
            sum(all_false_flights) / len(all_false_flights)
        )
        freq_true_flights_unbinned.append(sum(all_true_flights) / len(all_true_flights))

        for j in range(1, params.max_group_size + 1, params.group_bin_size):
            freq_false_flights = []
            freq_true_flights = []
            for group_size in range(j, j + params.group_bin_size):
                if group_size in all_false_flights_by_group:
                    freq_false_flights.extend(all_false_flights_by_group[group_size])
                if group_size in all_true_flights_by_group:
                    freq_true_flights.extend(all_true_flights_by_group[group_size])

            freq_false_flights_binned[-1].append(
                None
                if len(freq_false_flights) == 0
                else sum(freq_false_flights) / len(freq_false_flights)
            )
            freq_true_flights_binned[-1].append(
                None
                if len(freq_true_flights) == 0
                else sum(freq_true_flights) / len(freq_true_flights)
            )

    return Results(
        freq_false_flights_binned,
        freq_true_flights_binned,
        freq_false_flights_unbinned,
        freq_true_flights_unbinned,
        fitness_stat,
        trait_values,
    )


def mult_sim_analysis(
    *, out_file_path: str, params: Optional[OutputParameters], plots: List[str]
) -> None:
    sim_outputs = get_outputs(out_file_path, params)
    results = process_results(sim_outputs, params)

    if "flight_freq_by_group_size" in plots:
        plot_false_flight_freq(results, params)
        plot_true_flight_freq(results, params)

    if "fitness" in plots:
        plot_fitness(results)

    if "all_mean_trait_values" in plots:
        plot_all_mean_trait_values(results)
