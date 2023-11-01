import matplotlib.pyplot as plt
from si_types import *
from constants import *


def plot_flight_freq_by_group_size(
    flights_binned: List[List[Optional[float]]], params: Parameters, flight_type: str
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

    plt.legend(title="Group Size", loc="upper right")
    plt.title(f"Frequency of {flight_type} Flights Across Generations")
    plt.xlabel("Generation")
    plt.ylabel(f"Freq {flight_type} Flight")


def plot_false_flight_freq(r: Results, params: Parameters) -> None:
    plot_flight_freq_by_group_size(r.freq_false_flights_binned, params, "False")


def plot_true_flight_freq(r: Results, params: Parameters) -> None:
    plot_flight_freq_by_group_size(r.freq_true_flights_binned, params, "True")
