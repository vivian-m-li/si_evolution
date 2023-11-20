import matplotlib.cm as cm
from si_types import *

out_file_path = "results"

DEFAULT_PARAMS = Parameters()

COLOR_MAP = [x for y in [cm.Set2.colors, cm.Set1.colors] for x in y]

PARAM_FUNCS = {
    "prob_pred": AnalysisParam(
        label="Prob of Predation", func=lambda p, r=None: p.prob_pred
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
        label="Average Group Size", func=lambda p, r: r.avg_group_size
    ),
}

prev_sim_params = [
    Parameters(prob_pred=0.2),
    Parameters(prob_pred=0.2, max_group_size=15),
    Parameters(prob_pred=0.2, max_group_size=50),
    Parameters(prob_pred=0.2, max_group_size=50, Ni=500),
    Parameters(prob_pred=0.002),
    Parameters(prob_pred=0.005),
    Parameters(prob_pred=0.01),
    Parameters(prob_pred=0.02),
    Parameters(prob_pred=0.1),
    Parameters(prob_pred=0.2, max_group_size=15, Ni=500),
    Parameters(prob_pred=0.2, max_group_size=25, Ni=500),
    Parameters(prob_pred=0.04),
    Parameters(prob_pred=0.06),
    Parameters(prob_pred=0.08),
    Parameters(prob_pred=0.2, max_group_size=15, Ni=500, maxf=1000),
    Parameters(prob_pred=0.2, max_group_size=25, Ni=500, maxf=1000),
    Parameters(prob_pred=0.2, max_group_size=50, Ni=500, maxf=1000),
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
    # From here on, number of deaths per group per timestep is capped to 1, default prob_pred = 0.02
    Parameters(prob_pred=0.002),
    Parameters(prob_pred=0.005),
    Parameters(prob_pred=0.01),
    Parameters(prob_pred=0.02),
    Parameters(prob_pred=0.04),
    Parameters(prob_pred=0.06),
    Parameters(prob_pred=0.08),
    Parameters(prob_pred=0.1),
    Parameters(prob_pred=0.2, max_group_size=15),
    Parameters(prob_pred=0.2, max_group_size=15, Ni=500),
    Parameters(prob_pred=0.2, max_group_size=15, Ni=500, maxf=1000),
    Parameters(prob_pred=0.2, max_group_size=25),
    Parameters(prob_pred=0.2, max_group_size=25, Ni=500),
    Parameters(prob_pred=0.2, max_group_size=25, Ni=500, maxf=1000),
    Parameters(prob_pred=0.2, max_group_size=50),
    Parameters(prob_pred=0.2, max_group_size=50, Ni=500),
    Parameters(prob_pred=0.2, max_group_size=50, Ni=500, maxf=1000),
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
