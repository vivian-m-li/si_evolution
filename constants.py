import matplotlib.cm as cm
from si_types import *

out_file_path = "results"

DEFAULT_PARAMS = Parameters()

COLOR_MAP = [x for y in [cm.Set2.colors, cm.Set1.colors] for x in y]

PARAM_FUNCS = {
    "prob_pred": AnalysisParam(label="Prob of Predation", func=lambda x: x.prob_pred),
    "max_group_size": AnalysisParam(
        label="Max Group Size", func=lambda x: x.max_group_size
    ),
}
