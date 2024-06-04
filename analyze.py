import pandas as pd
import numpy as np
from plot import *
from si_types import *
from constants import *
from helper import *
from collections import defaultdict
from typing import Optional, List, DefaultDict


def get_params_to_analyze(analysis):
    params = []
    param = ""

    if analysis == "pop_size":
        for max_group_size in [15, 25, 50]:
            for Ni in [100, 500]:
                params.append(
                    OutputParameters(
                        max_group_size=max_group_size,
                        Ni=Ni,
                        group_bin_size=int(max_group_size / 5),
                    )
                )
        return [params], "pop_size"

    if analysis == "group_size":
        for prob_pred in [0.2, 2]:
            params.append([])
            for max_group_size in [5, 10, 15, 20, 25, 30, 40, 50, 100]:
                params[-1].append(
                    OutputParameters(
                        prob_pred=prob_pred,
                        max_group_size=max_group_size,
                        group_bin_size=int(max_group_size / 5),
                    )
                )
        return params, "avg_group_size"

    if analysis == "pred":
        for prob_pred in [
            0,
            0.002,
            0.005,
            0.01,
            0.02,
            0.04,
            0.06,
            0.08,
            0.1,
            0.12,
            0.14,
            0.16,
            0.18,
            0.2,
        ]:
            params.append(OutputParameters(prob_pred=prob_pred))
        return [params], "prob_pred"

    if analysis == "lambda":
        for p_lambda in [
            0,
            0.02,
            0.05,
            0.1,
            0.2,
            0.4,
            0.6,
            0.8,
            1,
            1.2,
            1.4,
            1.6,
            1.8,
            2.0,
        ]:
            params.append(OutputParameters(prob_pred=p_lambda))
        return [params], "lambda"

    if analysis == "e_gain":
        for prob_pred in [0.2, 2.0]:
            params.append([])
            for e_gain in [0.5, 1, 1.5, 2]:
                params[-1].append(OutputParameters(prob_pred=prob_pred, e_gain=e_gain))
        return params, "e_gain"

    return params, param


def process_results(
    sim_outputs: List[pd.DataFrame], params: OutputParameters
) -> Results:
    num_generations = params.maxf

    freq_false_flight_by_group_size: List[DefaultDict[int, List[float]]] = []
    freq_true_flight_by_group_size: List[DefaultDict[int, List[float]]] = []
    freq_detected_pred_deaths_all: List[float] = []
    freq_nondetected_pred_deaths_all: List[float] = []
    fitness_stat: List[Stat] = []
    trait_values: List[List[Stat]] = [[], [], []]
    deaths_stat: List[Stat] = []
    num_groups_per_gen: List[Stat] = []
    pred_catch_stat: List[Stat] = []
    all_group_sizes_stat: List[Stat] = []
    all_group_sizes: List[int] = []
    all_prop_pred_visits: List[List[float]] = []
    for i in range(num_generations):
        freq_false_flight_by_group_size.append(defaultdict(list))
        freq_true_flight_by_group_size.append(defaultdict(list))
        freq_detected_pred_deaths_gen: List[float] = []
        freq_nondetected_pred_deaths_gen: List[float] = []
        group_stats: Dict[str, GroupStats] = {
            "fitness": GroupStats(means=[], vars=[]),
            "f_pred": GroupStats(means=[], vars=[]),
            "s_faith": GroupStats(means=[], vars=[]),
            "s_dd": GroupStats(means=[], vars=[]),
        }
        num_groups: List[float] = []
        all_deaths: List[int] = []
        all_catches: List[float] = []
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
                group_sizes,
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
                pred_catch_rate,
                pred_catch_by_group_size,
                prop_groups_attacked,
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
            freq_detected_pred_deaths_gen.append(freq_detected_pred_deaths)
            freq_nondetected_pred_deaths_gen.append(freq_nondetected_pred_deaths)
            all_group_sizes_stat.append(
                Stat(mean=group_size_mean, variance=group_size_var)
            )
            all_group_sizes.extend(group_sizes)
            num_groups.append(params.Ni / group_size_mean)
            all_deaths.append(total_deaths)
            all_catches.append(pred_catch_rate)
            all_prop_pred_visits.append(prop_groups_attacked)

        freq_detected_pred_deaths_all.append(calc_mean(freq_detected_pred_deaths_gen))
        freq_nondetected_pred_deaths_all.append(
            calc_mean(freq_nondetected_pred_deaths_gen)
        )
        fitness_stat.append(
            Stat(
                mean=calc_mean(group_stats["fitness"].means),
                variance=np.var(group_stats["fitness"].vars),
                confidence_interval=calc_confidence_interval(
                    group_stats["fitness"].means
                ),
            )
        )
        for j, trait in enumerate(["f_pred", "s_faith", "s_dd"]):
            trait_values[j].append(
                Stat(
                    mean=calc_mean(group_stats[trait].means),
                    variance=np.var(group_stats[trait].vars),
                    confidence_interval=calc_confidence_interval(
                        group_stats[trait].means
                    ),
                )
            )
        num_groups_per_gen.append(
            Stat(mean=calc_mean(num_groups), variance=np.var(num_groups))
        )
        deaths_stat.append(
            Stat(mean=calc_mean(all_deaths), variance=np.var(all_deaths))
        )
        pred_catch_stat.append(
            Stat(mean=calc_mean(all_catches), variance=np.var(all_catches))
        )

    all_group_size_means = np.array([x.mean for x in all_group_sizes_stat])
    all_group_size_vars = np.array([x.variance for x in all_group_sizes_stat])
    avg_group_size = Stat(
        mean=round(np.mean(all_group_size_means), 2),
        variance=np.mean(all_group_size_vars),
        confidence_interval=(
            np.percentile(all_group_sizes, 2.5),
            np.percentile(all_group_sizes, 97.5),
        ),
    )
    all_prop_pred_visits = np.array(all_prop_pred_visits)
    avg_prop_pred_visits = Stat(
        np.mean(all_prop_pred_visits), np.var(all_prop_pred_visits)
    )
    prop_pred_visit_means = np.mean(all_prop_pred_visits, axis=0)
    prop_pred_visit_vars = np.var(all_prop_pred_visits, axis=0)
    prop_pred_visits_by_timestep: List[Stat] = [
        Stat(mean=m, variance=v)
        for m, v in zip(prop_pred_visit_means, prop_pred_visit_vars)
    ]

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

        freq_false_flights_unbinned.append(calc_mean(all_false_flights))
        freq_true_flights_unbinned.append(calc_mean(all_true_flights))

        if params.group_bin_size is None:
            continue

        for j in range(1, params.max_group_size + 1, params.group_bin_size):
            freq_false_flights = []
            freq_true_flights = []
            for group_size in range(j, j + params.group_bin_size):
                if group_size in all_false_flights_by_group:
                    freq_false_flights.extend(all_false_flights_by_group[group_size])
                if group_size in all_true_flights_by_group:
                    freq_true_flights.extend(all_true_flights_by_group[group_size])

            freq_false_flights_binned[-1].append(
                None if len(freq_false_flights) == 0 else calc_mean(freq_false_flights)
            )
            freq_true_flights_binned[-1].append(
                None if len(freq_true_flights) == 0 else calc_mean(freq_true_flights)
            )

    return Results(
        freq_false_flights_binned,
        freq_true_flights_binned,
        freq_false_flights_unbinned,
        freq_true_flights_unbinned,
        freq_detected_pred_deaths_all,
        freq_nondetected_pred_deaths_all,
        fitness_stat,
        trait_values,
        avg_group_size,
        num_groups_per_gen,
        deaths_stat,
        pred_catch_stat,
        avg_prop_pred_visits,
        prop_pred_visits_by_timestep,
    )


def mult_sim_analysis(
    *,
    out_file_path: str,
    all_params: List[OutputParameters],
    plots: List[str],
    param: Optional[str] = None,
) -> None:
    all_outputs = get_all_outputs(out_file_path, all_params)
    all_results: List[MultResults] = []
    for i, params in enumerate(all_params):
        sim_outputs = all_outputs[i]
        if len(sim_outputs) == 0:
            continue
        results = process_results(sim_outputs, params)
        all_results.append(MultResults(params, results))

        if "flight_freq_by_group_size" in plots:
            plot_false_flight_freq(results, params)
            plot_true_flight_freq(results, params)

        if "fitness" in plots:
            plot_fitness(results)

        if "all_mean_trait_values" in plots:
            plot_all_mean_trait_values(results)

    if param is not None:
        analysis_param = PARAM_FUNCS[param]

        if "avg_flight" in plots:
            plot_avg_false_flight(all_results, analysis_param)
            plot_avg_true_flight(all_results, analysis_param)

        if "detected_nondetected_pred_deaths" in plots:
            plot_detected_nondetected_pred_deaths(all_results, analysis_param)

        if "total_deaths_per_gen" in plots:
            plot_total_deaths_per_gen(all_results, analysis_param)

        if "final_fitness" in plots:
            plot_final_fitness(all_results, analysis_param)

        if "final_trait_values" in plots:
            plot_final_trait_values(all_results, analysis_param)

        if "final_flight_freq" in plots:
            plot_final_flight_freq(all_results, analysis_param)

        if "kills_per_visits_per_gen" in plots:
            plot_kills_per_visits_per_gen(all_results, analysis_param)

        if "final_kills_per_visits" in plots:
            plot_final_kills_per_visits(all_results, analysis_param)

        if "traits_by_gen" in plots:
            plot_traits_by_gen(all_results, analysis_param)

        if "prob_pred_by_lambda" in plots:
            plot_prob_pred_by_lambda(all_results, analysis_param)

        if "prob_pred_by_lambda_per_timestep" in plots:
            plot_prob_pred_by_lambda_per_timestep(all_results, analysis_param)


def run_mult_sim_analysis(params, param):
    for param_set in params:
        mult_sim_analysis(
            out_file_path=OUT_FILE_DIR,
            all_params=param_set,
            plots=[
                # "flight_freq_by_group_size",
                # "fitness",
                # "all_mean_trait_values",
                # "avg_flight",
                # "detected_nondetected_pred_deaths",
                # "total_deaths_per_gen",
                # "final_fitness",
                # "final_trait_values",
                "final_flight_freq",
                # "kills_per_visits_per_gen",
                # "final_kills_per_visits",
                "traits_by_gen",
                # "prob_pred_by_lambda",
                # "prob_pred_by_lambda_per_timestep",
            ],
            param=param,
        )


if __name__ == "__main__":
    params, param = get_params_to_analyze("lambda")
    run_mult_sim_analysis(params, param)
