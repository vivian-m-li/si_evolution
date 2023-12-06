import numpy as np
import matplotlib.cm as cm
from si_types import *

out_file_path = "results"

DEFAULT_PARAMS = Parameters()

COLOR_MAP = [x for y in [cm.Set2.colors, cm.Set1.colors] for x in y]

PARAM_FUNCS = {
    "prob_pred": AnalysisParam(
        label="Prob of Visitation", func=lambda p, r=None: p.prob_pred
    ),
    "max_group_size": AnalysisParam(
        label="Max Group Size", func=lambda p, r=None: p.max_group_size
    ),
    "e_gain": AnalysisParam(label="Energetic Gain", func=lambda p, r=None: p.e_gain),
    "e_gain/prob_pred": AnalysisParam(
        label="Energetic Gain/Probability of Predation",
        func=lambda p, r=None: p.e_gain / p.prob_pred,
        label_func=lambda p, r=None: f"{p.e_gain}/{p.prob_pred}",
    ),
    "avg_group_size": AnalysisParam(
        label="Average Group Size",
        func=lambda p, r: r.avg_group_size.mean,
    ),
    "pop_size": AnalysisParam(
        label="Average Group Size",
        func=lambda p, r: r.avg_group_size.mean,
        error_func=lambda p, r: (
            max(r.avg_group_size.mean - (1.96 * np.sqrt(r.avg_group_size.variance)), 1),
            min(
                r.avg_group_size.mean + (1.96 * np.sqrt(r.avg_group_size.variance)),
                p.max_group_size,
            ),
        ),
        label_func=lambda p, r=None: f"Ni={p.Ni}",
    ),
}

prev_sim_params = [
    # varying prob pred
    Parameters(prob_pred=0),
    Parameters(prob_pred=0.002),
    Parameters(prob_pred=0.005),
    Parameters(prob_pred=0.01),
    Parameters(prob_pred=0.02),
    Parameters(prob_pred=0.04),
    Parameters(prob_pred=0.06),
    Parameters(prob_pred=0.08),
    Parameters(prob_pred=0.1),
    Parameters(prob_pred=0.12),
    Parameters(prob_pred=0.14),
    Parameters(prob_pred=0.16),
    Parameters(prob_pred=0.18),
    Parameters(prob_pred=0.2),
    # varying max group size, prob pred = 0.02
    Parameters(prob_pred=0.02, max_group_size=5),
    Parameters(prob_pred=0.02, max_group_size=10),
    Parameters(prob_pred=0.02, max_group_size=15),
    Parameters(prob_pred=0.02, max_group_size=20),
    # Parameters(prob_pred=0.02, max_group_size=25),
    Parameters(prob_pred=0.02, max_group_size=30),
    Parameters(prob_pred=0.02, max_group_size=40),
    Parameters(prob_pred=0.02, max_group_size=50),
    Parameters(prob_pred=0.02, max_group_size=100),
    # varying max group size, prob pred = 0.02
    Parameters(prob_pred=0.2, max_group_size=5),
    Parameters(prob_pred=0.2, max_group_size=10),
    Parameters(prob_pred=0.2, max_group_size=15),
    Parameters(prob_pred=0.2, max_group_size=20),
    # Parameters(prob_pred=0.2, max_group_size=25),
    Parameters(prob_pred=0.2, max_group_size=30),
    Parameters(prob_pred=0.2, max_group_size=40),
    Parameters(prob_pred=0.2, max_group_size=50),
    Parameters(prob_pred=0.2, max_group_size=100),
    # varying max group size, prob pred = 0.02, Ni = 500
    Parameters(prob_pred=0.02, max_group_size=5, Ni=500),
    Parameters(prob_pred=0.02, max_group_size=10, Ni=500),
    Parameters(prob_pred=0.02, max_group_size=15, Ni=500),
    Parameters(prob_pred=0.02, max_group_size=20, Ni=500),
    Parameters(prob_pred=0.02, max_group_size=25, Ni=500),
    Parameters(prob_pred=0.02, max_group_size=30, Ni=500),
    Parameters(prob_pred=0.02, max_group_size=40, Ni=500),
    Parameters(prob_pred=0.02, max_group_size=50, Ni=500),
    Parameters(prob_pred=0.02, max_group_size=100, Ni=500),
    # varying max group size, prob pred = 0.2, Ni = 500
    Parameters(prob_pred=0.2, max_group_size=5, Ni=500),
    Parameters(prob_pred=0.2, max_group_size=10, Ni=500),
    Parameters(prob_pred=0.2, max_group_size=15, Ni=500),
    Parameters(prob_pred=0.2, max_group_size=20, Ni=500),
    Parameters(prob_pred=0.2, max_group_size=25, Ni=500),
    Parameters(prob_pred=0.2, max_group_size=30, Ni=500),
    Parameters(prob_pred=0.2, max_group_size=40, Ni=500),
    Parameters(prob_pred=0.2, max_group_size=50, Ni=500),
    Parameters(prob_pred=0.2, max_group_size=100, Ni=500),
    # varying prob pred and e gain
    Parameters(prob_pred=0.02, e_gain=0.5),
    Parameters(prob_pred=0.02, e_gain=1.5),
    Parameters(prob_pred=0.02, e_gain=2),
    Parameters(prob_pred=0.06, e_gain=0.5),
    Parameters(prob_pred=0.06, e_gain=1.5),
    Parameters(prob_pred=0.06, e_gain=2),
    Parameters(prob_pred=0.08, e_gain=0.5),
    Parameters(prob_pred=0.08, e_gain=1.5),
    Parameters(prob_pred=0.08, e_gain=2),
    Parameters(prob_pred=0.1, e_gain=0.5),
    Parameters(prob_pred=0.1, e_gain=1.5),
    Parameters(prob_pred=0.1, e_gain=2),
    Parameters(prob_pred=0.2, e_gain=0.5),
    Parameters(prob_pred=0.2, e_gain=1.5),
    Parameters(prob_pred=0.2, e_gain=2),
]
