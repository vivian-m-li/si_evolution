import matplotlib.pyplot as plt
from si_types import *
from constants import *


def plot_flight_freq_by_group_size(
    flights_binned: List[List[Optional[float]]],
    flights_unbinned: List[float],
    params: Parameters,
    flight_type: str,
) -> None:
    plt.figure(figsize=(10, 6))

    x = list(range(1, params.maxf + 1))
    for bin_idx in range(len(flights_binned[0])):
        y = [gen[bin_idx] for gen in flights_binned]
        bin_lower_bound = bin_idx * GROUP_BIN_SIZE
        plt.plot(
            x,
            y,
            label=f"{bin_lower_bound + 1}-{bin_lower_bound + GROUP_BIN_SIZE}",
            color=COLOR_MAP[bin_idx],
        )
    plt.plot(x, flights_unbinned, label="All group sizes", color="#555555")

    plt.legend(title="Group Size", loc="upper right")
    plt.title(
        f"Frequency of {flight_type} Flights Across Generations, Max Group Size = {params.max_group_size}"
    )
    plt.xlabel("Generation")
    plt.ylabel(f"Freq {flight_type} Flight")


def plot_false_flight_freq(r: Results, params: Parameters) -> None:
    plot_flight_freq_by_group_size(
        r.freq_false_flights_binned, r.freq_false_flights_unbinned, params, "False"
    )


def plot_true_flight_freq(r: Results, params: Parameters) -> None:
    plot_flight_freq_by_group_size(
        r.freq_true_flights_binned, r.freq_true_flights_unbinned, params, "True"
    )


def plot_mean_trait_values(trait: str) -> None:
    return


def plot_all_mean_trait_values(r: Results) -> None:
    traits = ["jumpiness", "ociality", "density dependence in sociality"]
    for i, trait in enumerate(traits):
        plot_mean_trait_values(trait)


def plot_fitness(r: Results) -> None:
    x = list(range(1, len(r.fitness_stat) + 1))
    y = [x.mean for x in r.fitness_stat]
    y_err = [x.variance for x in r.fitness_stat]

    plt.figure(figsize=(10, 6))
    plt.plot(x, y)
    # plt.errorbar(
    #     x,
    #     y,
    #     yerr=y_err,
    #     ecolor="red",
    #     fmt="o-",
    # )

    plt.xlabel("Generation")
    plt.ylabel("Fitness")
    plt.title("Mean Fitness Across Generations")

    plt.show()
