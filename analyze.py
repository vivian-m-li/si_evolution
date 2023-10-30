import os
import csv
import pandas as pd
from si_types import *


def get_sim_id(file_name: str) -> int:
    df = pd.read_csv(file_name)
    df.dropna(axis=0, how="all", inplace=True)
    df.to_csv(file_name, index=False)

    sim_id = 0
    if not df.empty:
        sim_id = int(df.iloc[:, 0].max()) + 1

    return sim_id


def write_output(directory: str, sim_id: int, output: SimOutput):
    results_file_name = f"{directory}/all.csv"

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
