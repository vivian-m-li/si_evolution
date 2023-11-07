from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple


@dataclass
class Parameters:
    Ni: int = 100  # number of individuals
    tf: int = 30  # time steps for a given generation
    e_gain: int = 1  # energy units gained each time step in which the ith individual does not flee
    coef_false: float = 0.2  # coefficient that determines the prob of a false alarm (a smaller value than f_pred)
    maxf: int = 500  # number of generations to run the model through
    prob_pred: float = 0.2
    max_group_size: int = 25


@dataclass
class OutputParameters(Parameters):
    group_bin_size: int = 5


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
    fitness_stat: List[Stat]
    trait_values: List[List[Stat]]


@dataclass
class GroupStats:
    means: List[float]
    vars: List[float]
