from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple, Callable


@dataclass
class Parameters:
    Ni: int = 100  # number of individuals
    tf: int = 30  # time steps for a given generation
    e_gain: int = 1  # energy units gained each time step in which the ith individual does not flee
    coef_false: float = 0.2  # coefficient that determines the prob of a false alarm (a smaller value than f_pred)
    maxf: int = 500  # number of generations to run the model through
    prob_pred: float = 0.02
    max_group_size: int = 25


@dataclass
class OutputParameters(Parameters):
    group_bin_size: Optional[int] = None


@dataclass
class Stat:
    mean: float
    variance: float
    confidence_interval: Optional[Tuple[float, float]] = None


@dataclass
class SimOutput:
    parameters: Parameters
    total_deaths: List[int]  # per generation
    false_flights: List[Dict[int, float]]  # per generation, per group size
    true_flights: List[Dict[int, float]]  # per generation, per group size
    detected_pred_deaths: List[int]  # per generation
    nondetected_pred_deaths: List[int]  # per generation
    trait_values: List[
        List[Stat]
    ]  # per generation - list of stats for f_pred, s_faith, s_dd
    energetic_states: List[Stat]  # per generation
    fitness: List[Stat]  # per generation
    group_size: List[Stat]  # per generation
    all_group_sizes: List[int]
    pred_catch_rate: List[float]  # per generation
    pred_catch_by_group_size: List[
        Dict[int, float]
    ]  # num_caught/num_attacks mean per group size per per generation


@dataclass
class Flights:
    false_true_flights: int
    total_flights: int


@dataclass
class Results:
    freq_false_flights_binned: List[List[Optional[float]]]
    freq_true_flights_binned: List[List[Optional[float]]]
    freq_false_flights_unbinned: List[float]
    freq_true_flights_unbinned: List[float]
    freq_detected_pred_deaths_all: List[float]
    freq_nondetected_pred_deaths_all: List[float]
    fitness_stat: List[
        Stat
    ]  # mean, var of the means/vars of each simulation, per generation
    trait_values: List[List[Stat]]
    avg_group_size: Stat
    num_groups_per_gen: List[Stat]
    deaths_stat: List[Stat]
    pred_catch_stat: List[Stat]


@dataclass
class GroupStats:
    means: List[float]
    vars: List[float]


@dataclass
class MultResults:
    params: Parameters
    results: Results


@dataclass
class AnalysisParam:
    label: str
    func: Callable[[Parameters, Optional[Results]], float]
    label_func: Optional[Callable[[Parameters, Optional[Results]], str]] = None
    error_func: Optional[
        Callable[[Parameters, Optional[Results]], Tuple[float, float]]
    ] = None
