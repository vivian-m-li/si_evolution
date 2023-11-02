import matplotlib.cm as cm
from si_types import *

out_file_path = "results"

DEFAULT_PARAMS = Parameters(
    Ni=100,  # number of individuals
    tf=30,  # time steps for a given generation
    e_gain=1,  # energy units gained each time step in which the ith individual does not flee
    coef_false=0.2,  # coefficient that determines the prob of a false alarm (a smaller value than f_pred)
    maxf=500,  # number of generations to run the model through
    prob_pred=0.2,
    max_group_size=25,
)

GROUP_BIN_SIZE = 10

COLOR_MAP = [x for y in [cm.Set2.colors, cm.Set1.colors] for x in y]
