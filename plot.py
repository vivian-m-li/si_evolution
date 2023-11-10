import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from si_types import *
from constants import *


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
    plt.show()

    return out_ls


def plot_flight_freq_by_group_size(
    flights_binned: List[List[Optional[float]]],
    flights_unbinned: List[float],
    params: OutputParameters,
    flight_type: str,
) -> None:
    plt.figure(figsize=(10, 6))

    x = list(range(1, params.maxf + 1))
    for bin_idx in range(len(flights_binned[0])):
        y = [gen[bin_idx] for gen in flights_binned]
        bin_lower_bound = bin_idx * params.group_bin_size
        plt.plot(
            x,
            y,
            label=f"{bin_lower_bound + 1}-{bin_lower_bound + params.group_bin_size}",
            color=COLOR_MAP[bin_idx],
        )
    plt.plot(x, flights_unbinned, label="All group sizes", color="#555555")

    plt.legend(title="Group Size", loc="upper right")
    plt.title(
        f"Frequency of {flight_type} Flights Across Generations, Max Group Size = {params.max_group_size}"
    )
    plt.xlabel("Generation")
    plt.ylabel(f"Freq {flight_type} Flight")
    plt.show()


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
    plt.show()


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
    plt.show()


def plot_avg_flight(
    results: List[Tuple[Parameters, List[float]]],
    flight_type: str,
    param: AnalysisParam,
) -> None:
    plt.figure(figsize=(10, 6))

    x = list(range(1, results[0][0].maxf + 1))
    for i, r in enumerate(results):
        params, flights = r
        plt.plot(x, flights, label=param.func(params), color=COLOR_MAP[i])

    plt.legend(title=param.label, loc="upper right")
    plt.title(
        f"Average Frequency of {flight_type} Flights Across Generations At Varying {param.label}"
    )
    plt.xlabel("Generation")
    plt.ylabel(f"Freq {flight_type} Flight")
    plt.show()


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
        color = COLOR_MAP[i]
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
            Line2D([0], [0], color=color, label=param.func(r.params))
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
    plt.ylabel(f"Freq Death")
    plt.show()
