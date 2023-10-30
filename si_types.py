from dataclasses import dataclass
from typing import List, Dict


@dataclass
class Parameters:
    Ni: int
    tf: int
    e_gain: int
    coef_false: float
    maxf: int
    prob_pred: float
    max_group_size: int


@dataclass
class Stat:
    mean: float
    variance: float


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
    group_size: List[Stat]  # per generation


@dataclass
class Flights:
    false_true_flights: int
    total_flights: int
