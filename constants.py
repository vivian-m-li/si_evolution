import numpy as np
import matplotlib.cm as cm
from si_types import *

OUT_FILE_DIR = "results"
PLOT_FILE_DIR = "plots"

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
        error_func=lambda p, r: (
            r.avg_group_size.confidence_interval[0],
            r.avg_group_size.confidence_interval[1],
        ),
    ),
    "pop_size": AnalysisParam(
        label="Average Group Size",
        func=lambda p, r: r.avg_group_size.mean,
        error_func=lambda p, r: (
            r.avg_group_size.confidence_interval[0],
            r.avg_group_size.confidence_interval[1],
        ),
        label_func=lambda p, r=None: f"Ni={p.Ni}",
    ),
    "lambda": AnalysisParam(
        label="Lambda",
        func=lambda p, r=None: p.prob_pred,
    ),
}
