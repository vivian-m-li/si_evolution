import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.stats import linregress
from collections import defaultdict
from si_types import *
from constants import *
from helper import *
from typing import Union


def get_color(index: int, num_colors: int) -> Tuple[float, float, float]:
    cmap_blues = plt.get_cmap("Blues")
    return cmap_blues(0.3 + index / num_colors)


def plot_traits(trait_mean, trait_sd, fit_traits_gen):
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
    plt.savefig(f"{PLOT_FILE_DIR}/traits.png")

    return out_ls


def plot_flight_freq_by_group_size(
    flights_binned: List[List[Optional[float]]],
    flights_unbinned: List[float],
    params: OutputParameters,
    flight_type: str,
) -> None:
    plt.figure(figsize=(10, 6))

    x = list(range(1, params.maxf + 1))
    num_bins = len(flights_binned[0])
    for bin_idx in range(num_bins):
        y = [gen[bin_idx] for gen in flights_binned]
        bin_lower_bound = bin_idx * params.group_bin_size
        plt.plot(
            x,
            y,
            label=f"{bin_lower_bound + 1}-{bin_lower_bound + params.group_bin_size}",
            color=get_color(bin_idx, num_bins),
        )
    plt.plot(x, flights_unbinned, label="All group sizes", color="#555555")

    plt.legend(title="Group Size", loc="upper right")
    plt.title(
        f"Frequency of {flight_type} Flights Across Generations, Max Group Size = {params.max_group_size}"
    )
    plt.xlabel("Generation")
    plt.ylabel(f"Proportion of Group Fleeing ({flight_type})")
    plt.savefig(f"{PLOT_FILE_DIR}/flight_freq_by_group_size.png")


def plot_false_flight_freq(r: Results, params: OutputParameters) -> None:
    plot_flight_freq_by_group_size(
        r.freq_false_flights_binned, r.freq_false_flights_unbinned, params, "False"
    )


def plot_true_flight_freq(r: Results, params: OutputParameters) -> None:
    plot_flight_freq_by_group_size(
        r.freq_true_flights_binned, r.freq_true_flights_unbinned, params, "True"
    )


def plot_mean_trait_values(trait: str, values: List[Stat]) -> None:
    x = list(range(1, len(values) + 1))
    y = [val.mean for val in values]

    plt.figure(figsize=(10, 6))
    plt.plot(x, y)
    plt.fill_between(
        x,
        [val.confidence_interval[0] for val in values],
        [val.confidence_interval[1] for val in values],
        alpha=0.2,
    )

    plt.xlabel("Generation")
    plt.ylabel(trait)
    plt.title(f"Mean {trait} Across Generations")
    plt.savefig(f"{PLOT_FILE_DIR}/mean_trait_values.png")


def plot_all_mean_trait_values(r: Results) -> None:
    traits = ["jumpiness", "sociality", "density dependence in sociality"]
    for i, trait in enumerate(traits):
        plot_mean_trait_values(trait, r.trait_values[i])


def plot_fitness(r: Results) -> None:
    x = list(range(1, len(r.fitness_stat) + 1))
    y = [x.mean for x in r.fitness_stat]

    plt.figure(figsize=(10, 6))
    plt.plot(x, y)
    plt.fill_between(
        x,
        [val.confidence_interval[0] for val in r.fitness_stat],
        [val.confidence_interval[1] for val in r.fitness_stat],
        alpha=0.2,
    )

    plt.xlabel("Generation")
    plt.ylabel("Fitness")
    plt.title("Mean Fitness Across Generations")
    plt.savefig(f"{PLOT_FILE_DIR}/fitness.png")


def plot_fitness_by_gen(results: List[MultResults], param: AnalysisParam) -> None:
    legend_elements = []
    x = list(range(1, results[0].params.maxf + 1))

    plot_data = []
    for j, r in enumerate(results):
        color = get_color(j, len(results))
        legend_elements.append(
            Line2D([0], [0], color=color, label=param.func(r.params, r.results))
        )
        means = [y.mean for y in r.results.fitness_stat]
        plot_data.append(means)
        ci_lower = [y.confidence_interval[0] for y in r.results.fitness_stat]
        ci_upper = [y.confidence_interval[1] for y in r.results.fitness_stat]
        plt.fill_between(x, ci_lower, ci_upper, color=color, alpha=0.3)

    # ensures lines show up in front of confidence intervals
    for j, r in enumerate(results):
        color = get_color(j, len(results))
        plt.plot(
            x,
            plot_data[j],
            color=color,
        )

    plt.ylabel("Fitness")
    plt.legend(
        title=param.label,
        handles=legend_elements,
        loc="upper right",
    )
    plt.title(f"Mean Fitness Across Generations at Varying {param.label}")
    plt.xlabel("Generation")
    plt.tight_layout()
    plt.savefig(f"{PLOT_FILE_DIR}/fitness_by_gen_{param.label}.png")


def plot_avg_flight(
    results: List[Tuple[Parameters, List[float]]],
    flight_type: str,
    param: AnalysisParam,
) -> None:
    plt.figure(figsize=(10, 6))

    x = list(range(1, results[0][0].maxf + 1))
    for i, r in enumerate(results):
        params, flights = r
        plt.plot(
            x,
            flights,
            label=param.func(params),
            color=get_color(
                i,
                len(results),
            ),
        )

    plt.legend(title=param.label, loc="upper right")
    plt.title(
        f"Average Frequency of {flight_type} Flights Across Generations At Varying {param.label}"
    )
    plt.xlabel("Generation")
    plt.ylabel(f"Proportion of Group Fleeing ({flight_type})")
    plt.savefig(f"{PLOT_FILE_DIR}/avg_flight_{param.label}.png")


def plot_avg_false_flight(results: List[MultResults], param: AnalysisParam) -> None:
    plot_avg_flight(
        [(r.params, r.results.freq_false_flights_unbinned) for r in results],
        "False",
        param,
    )


def plot_avg_true_flight(results: List[MultResults], param: AnalysisParam) -> None:
    plot_avg_flight(
        [(r.params, r.results.freq_true_flights_unbinned) for r in results],
        "True",
        param,
    )


def plot_detected_nondetected_pred_deaths(
    results: List[MultResults], param: AnalysisParam
) -> None:
    plt.figure(figsize=(10, 6))

    legend_elements = []
    x = list(range(1, results[0].params.maxf + 1))
    for i, r in enumerate(results):
        color = get_color(i, len(results))
        plt.plot(
            x,
            r.results.freq_detected_pred_deaths_all,
            color=color,
        )
        plt.plot(
            x,
            r.results.freq_nondetected_pred_deaths_all,
            color=color,
            linestyle=(0, (1, 5)),
        )
        legend_elements.append(
            Line2D([0], [0], color=color, label=param.func(r.params, r.results))
        )

    legend_elements.extend(
        [
            Line2D([0], [0], label=" ", color="white", markersize=0),
            Line2D([0], [0], color="black", label="Detected Predator Deaths"),
            Line2D(
                [0],
                [0],
                color="black",
                label="Nondetected Predator Deaths",
                linestyle=(0, (1, 5)),
            ),
        ]
    )

    plt.legend(title=param.label, handles=legend_elements, loc="upper right")
    plt.title(
        f"Detected and Nondetected Deaths from Predators At Varying {param.label}"
    )
    plt.xlabel("Generation")
    plt.ylim(-0.025, 1.025)
    plt.ylabel(f"Freq Death")
    plt.savefig(f"{PLOT_FILE_DIR}/detected_nondetected_pred_deaths_{param.label}.png")


def get_steady_state_value(values: Union[List[Stat], List[float]]) -> float:
    last_fifty_vals = values[-50:]
    if type(values[0]) == Stat:
        last_fifty_vals = [x.mean for x in last_fifty_vals]
    return np.mean(np.array(last_fifty_vals))


def plot_final_fitness(results: List[MultResults], param: AnalysisParam) -> None:
    plt.figure(figsize=(10, 6))
    data = {}
    labels = {}
    for r in results:
        x_val = param.func(r.params, r.results)
        data[x_val] = get_steady_state_value(r.results.fitness_stat)
        labels[x_val] = (
            param.label_func(r.params, r.results)
            if param.label_func is not None
            else ""
        )

    x_vals = []
    y_vals = []
    for x_val in sorted(data):
        x_vals.append(x_val)
        y_vals.append(data[x_val])

    plt.scatter(x_vals, y_vals)
    plt.plot(x_vals, y_vals, linestyle="dashed")
    for i, x_val in enumerate(x_vals):
        plt.annotate(
            labels[x_val],
            (x_val, y_vals[i]),
            textcoords="offset points",
            xytext=(5, 5),
            ha="center",
        )

    plt.title(f"Fitness at Varying {param.label}")
    plt.xlabel(param.label)
    plt.ylabel("Fitness")
    plt.savefig(f"{PLOT_FILE_DIR}/final_fitness_{param.label}.png")


def plot_final_trait_values(results: List[MultResults], param: AnalysisParam) -> None:
    plt.figure(figsize=(10, 6))
    data = {}
    for r in results:
        x_val = param.func(r.params, r.results)
        data[x_val] = []
        for i, trait_vals in enumerate(r.results.trait_values):
            data[x_val].append(get_steady_state_value(trait_vals))

    x_vals = []
    y_vals = [[], [], []]
    for x_val in sorted(data):
        x_vals.append(x_val)
        for i, y_val in enumerate(data[x_val]):
            y_vals[i].append(y_val)

    for i, label in enumerate(
        [
            "jumpiness",
            "sociality",
            "density dependence in sociality",
        ]
    ):
        plt.scatter(x_vals, y_vals[i], label=label, color=COLOR_MAP[i])
        plt.plot(x_vals, y_vals[i], linestyle="dashed", color=COLOR_MAP[i])

    plt.legend()
    plt.title(f"Trait Values at Varying {param.label}")
    plt.xlabel(param.label)
    plt.ylabel("Trait Value")
    plt.savefig(f"{PLOT_FILE_DIR}/final_trait_values_{param.label}.png")


def plot_total_deaths_per_gen(results: List[MultResults], param: AnalysisParam) -> None:
    plt.figure(figsize=(10, 6))

    legend_elements = []
    x = list(range(1, results[0].params.maxf + 1))
    for i, r in enumerate(results):
        color = get_color(i, len(results))
        plt.plot(
            x,
            [y.mean for y in r.results.deaths_stat],
            color=color,
        )
        legend_elements.append(
            Line2D([0], [0], color=color, label=param.func(r.params, r.results))
        )

    plt.legend(title=param.label, handles=legend_elements, loc="upper right")
    plt.title(f"Total Deaths from Predators At Varying {param.label}")
    plt.xlabel("Generation")
    plt.ylabel(f"# Deaths")
    plt.savefig(f"{PLOT_FILE_DIR}/total_deaths_per_gen_{param.label}.png")


def plot_final_flight_freq(results: List[MultResults], param: AnalysisParam) -> None:
    plt.figure(figsize=(10, 6))
    data = defaultdict(dict)
    labels = {}
    x_bounds = {}
    for r in results:
        x_val = param.func(r.params, r.results)
        data[x_val]["false_flights"] = get_steady_state_value(
            r.results.freq_false_flights_unbinned
        )
        data[x_val]["true_flights"] = get_steady_state_value(
            r.results.freq_true_flights_unbinned
        )
        labels[x_val] = (
            param.label_func(r.params, r.results)
            if param.label_func is not None
            else ""
        )
        x_bounds[x_val] = (
            param.error_func(r.params, r.results)
            if param.error_func is not None
            else None
        )

    x_vals = []
    freq_false_flights = []
    freq_true_flights = []
    for x_val in sorted(data):
        x_vals.append(x_val)
        freq_false_flights.append(data[x_val]["false_flights"])
        freq_true_flights.append(data[x_val]["true_flights"])

    false_best_fit = linregress(filter_nan(x_vals, freq_false_flights))
    true_best_fit = linregress(filter_nan(x_vals, freq_true_flights))
    print(
        f"Freq false flights={round(false_best_fit.slope, 2)}x + {round(false_best_fit.intercept, 2)}"
    )
    print(
        f"Freq true flights={round(true_best_fit.slope, 2)}x + {round(true_best_fit.intercept, 2)}"
    )

    for i, x_val in enumerate(x_vals):
        if param.label_func:
            plt.annotate(
                labels[x_val],
                (x_val, freq_false_flights[i]),
                textcoords="offset points",
                xytext=(5, 5),
                ha="center",
            )
        # if param.error_func:
        #     plt.plot(
        #         [x_bounds[x_val][0], x_bounds[x_val][1]],
        #         [freq_false_flights[i], freq_false_flights[i]],
        #         color="gray",
        #         alpha=0.5,
        #     )

    plt.scatter(
        x_vals, freq_false_flights, label="freq false flights", color=COLOR_MAP[0]
    )
    plt.plot(x_vals, freq_false_flights, linestyle="dashed", color=COLOR_MAP[0])

    plt.scatter(
        x_vals, freq_true_flights, label="freq true flights", color=COLOR_MAP[1]
    )
    plt.plot(x_vals, freq_true_flights, linestyle="dashed", color=COLOR_MAP[1])

    plt.legend()
    plt.title(f"Frequency of False & True Flights", fontsize=18)
    plt.xlabel(param.label, fontsize=16)
    plt.ylim(-0.025, 1.025)
    plt.ylabel("Freq. Flight", fontsize=16)
    plt.savefig(f"{PLOT_FILE_DIR}/final_flight_freq_{param.label}.png")


def plot_kills_per_visits_per_gen(results: List[MultResults], param: AnalysisParam):
    plt.figure(figsize=(10, 6))

    legend_elements = []
    x = list(range(1, results[0].params.maxf + 1))
    for i, r in enumerate(results):
        color = get_color(i, len(results))
        plt.plot(
            x,
            [
                y.mean
                # * r.params.prob_pred * (r.params.tf - 1)
                for y in r.results.pred_catch_stat
            ],
            color=color,
        )
        legend_elements.append(
            Line2D([0], [0], color=color, label=param.func(r.params, r.results))
        )

    plt.legend(title=param.label, handles=legend_elements, loc="upper right")
    plt.title(f"Predator Catch Rate At Varying {param.label}")
    plt.xlabel("Generation")
    plt.ylim(-0.025, 1.025)
    plt.ylabel(f"Catch Rate (# Kills/# Visits)")
    plt.savefig(f"{PLOT_FILE_DIR}/kills_per_visits_per_gen_{param.label}.png")


def plot_final_kills_per_visits(
    results: List[MultResults], param: AnalysisParam
) -> None:
    plt.figure(figsize=(10, 6))
    data = {}
    labels = {}
    for r in results:
        x_val = param.func(r.params, r.results)
        data[x_val] = (
            get_steady_state_value(r.results.pred_catch_stat)
            # * r.params.prob_pred
            # * (r.params.tf - 1)
        )
        labels[x_val] = (
            param.label_func(r.params, r.results)
            if param.label_func is not None
            else ""
        )

    x_vals = []
    y_vals = []
    for x_val in sorted(data):
        x_vals.append(x_val)
        y_vals.append(data[x_val])

    plt.scatter(x_vals, y_vals)
    plt.plot(x_vals, y_vals, linestyle="dashed")
    for i, x_val in enumerate(x_vals):
        plt.annotate(
            labels[x_val],
            (x_val, y_vals[i]),
            textcoords="offset points",
            xytext=(5, 5),
            ha="center",
        )

    plt.title(f"Predator Catch Rate at Varying {param.label}")
    plt.xlabel(param.label)
    plt.ylim(-0.025, 1.025)
    plt.ylabel("Catch Rate (# Kills/# Visits)")
    plt.savefig(f"{PLOT_FILE_DIR}/final_kills_per_visits_{param.label}.png")


def plot_traits_by_gen(results: List[MultResults], param: AnalysisParam) -> None:
    plt.figure(figsize=(7, 12))
    gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 1])
    traits = ["jumpiness", "sociality", "density dependence in sociality"]
    for i, trait in enumerate(traits):
        ax = plt.subplot(gs[i, 0])

        legend_elements = []
        x = list(range(1, results[0].params.maxf + 1))

        plot_data = []
        for j, r in enumerate(results):
            color = get_color(j, len(results))
            legend_elements.append(
                Line2D([0], [0], color=color, label=param.func(r.params, r.results))
            )
            means = [y.mean for y in r.results.trait_values[i]]
            plot_data.append(means)
            ci_lower = [y.confidence_interval[0] for y in r.results.trait_values[i]]
            ci_upper = [y.confidence_interval[1] for y in r.results.trait_values[i]]
            plt.fill_between(x, ci_lower, ci_upper, color=color, alpha=0.3)

            is_stable, diff = calc_trait_stability(means)
            if not is_stable:
                print("not stable", r.params, diff)

        # ensures lines show up in front of confidence intervals
        for j, r in enumerate(results):
            color = get_color(j, len(results))
            plt.plot(
                x,
                plot_data[j],
                color=color,
            )

        ax.set_ylabel(trait)
        if i == 0:
            ax.legend(
                title=param.label,
                handles=legend_elements,
                loc="upper right",
            )
    plt.suptitle(f"Mean Trait Values Across Generations at Varying {param.label}")
    plt.xlabel("Generation")
    plt.tight_layout()
    plt.savefig(f"{PLOT_FILE_DIR}/traits_by_gen_{param.label}.png")


def plot_traits_by_gen_separated(
    results: List[MultResults], param: AnalysisParam
) -> None:
    plt.figure(figsize=(10, 6))
    legend_elements = []
    x = list(range(1, results[0].params.maxf + 1))
    for j, r in enumerate(results):
        color = get_color(j, len(results))
        plt.plot(
            x,
            [y.mean for y in r.results.trait_values[0]],
            color=color,
        )
        legend_elements.append(
            Line2D([0], [0], color=color, label=param.func(r.params, r.results))
        )

    plt.ylabel("jumpiness", fontsize=16)
    plt.legend(
        title=param.label,
        handles=legend_elements,
        loc="upper right",
    )
    plt.suptitle(f"Mean Trait Values Across Generations", fontsize=18)
    plt.xlabel("Generation", fontsize=16)
    plt.savefig(f"{PLOT_FILE_DIR}/jumpiness_{param.label}.png")

    plt.figure(figsize=(8, 10))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
    traits = ["sociality", "density dependence in sociality"]
    for i, trait in enumerate(traits):
        ax = plt.subplot(gs[i, 0])

        legend_elements = []
        x = list(range(1, results[0].params.maxf + 1))
        for j, r in enumerate(results):
            color = get_color(j, len(results))
            plt.plot(
                x,
                [y.mean for y in r.results.trait_values[i + 1]],
                color=color,
            )
            legend_elements.append(
                Line2D([0], [0], color=color, label=param.func(r.params, r.results))
            )

        ax.set_ylabel(trait, fontsize=16)
        if i == 0:
            ax.legend(
                title=param.label,
                handles=legend_elements,
                loc="upper right",
            )
    plt.suptitle(f"Mean Trait Values Across Generations", fontsize=18)
    plt.xlabel("Generation", fontsize=16)
    plt.tight_layout()
    plt.savefig(f"{PLOT_FILE_DIR}/sociality_sdd_{param.label}.png")


def plot_prob_pred_by_lambda(results: List[MultResults], param: AnalysisParam) -> None:
    plt.figure(figsize=(10, 6))
    data = {}
    labels = {}
    for r in results:
        x_val = param.func(r.params, r.results)
        data[x_val] = r.results.avg_prop_pred_visits.mean
        labels[x_val] = (
            param.label_func(r.params, r.results)
            if param.label_func is not None
            else ""
        )

    x_vals = []
    y_vals = []
    for x_val in sorted(data):
        x_vals.append(x_val)
        y_vals.append(data[x_val])

    plt.scatter(x_vals, y_vals)
    plt.plot(x_vals, y_vals, linestyle="dashed")
    for i, x_val in enumerate(x_vals):
        plt.annotate(
            labels[x_val],
            (x_val, y_vals[i]),
            textcoords="offset points",
            xytext=(5, 5),
            ha="center",
        )

    plt.title(f"% of Groups Attacked at Varying {param.label}")
    plt.xlabel(param.label)
    plt.ylabel(f"% of Groups Attacked")
    plt.savefig(f"{PLOT_FILE_DIR}/prob_pred_by_lambda_{param.label}.png")


def plot_prob_pred_by_lambda_per_timestep(
    results: List[MultResults], param: AnalysisParam
) -> None:
    plt.figure(figsize=(10, 6))

    legend_elements = []
    x = list(range(1, results[0].params.tf))
    for i, r in enumerate(results):
        color = get_color(i, len(results))
        plt.plot(
            x,
            [y.mean for y in r.results.prop_pred_visits_by_timestep],
            color=color,
        )
        legend_elements.append(
            Line2D([0], [0], color=color, label=param.func(r.params, r.results))
        )

    plt.legend(title=param.label, handles=legend_elements, loc="upper right")
    plt.title(f"% of Groups Attacked At Varying {param.label}")
    plt.xlabel("Timestep")
    plt.ylim(-0.025, 1.025)
    plt.ylabel(f"% of Groups Attacked")
    plt.savefig(f"{PLOT_FILE_DIR}/prob_pred_by_lambda_per_timestep_{param.label}.png")


def plot_trait_hist(trait: str, values: List[float]):
    plt.figure(figsize=(10, 6))
    plt.hist(values, bins=15)
    plt.title(f"Steady State Values for {trait}")
    plt.xlabel("Trait Value")
    plt.ylabel("Count")
    plt.show()


def plot_traits_scatter(trait_values: List[List[float]]):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(trait_values[:, 0], trait_values[:, 1], trait_values[:, 2])

    ax.set_xlabel("Jumpiness")
    ax.set_ylabel("Sociality")
    ax.set_zlabel("Density Dependence in Sociality")
    plt.show()
