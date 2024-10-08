import multiprocessing
import time
import datetime
from si_evolution import evo_fun
from constants import OUT_FILE_DIR
from si_types import Parameters

NUM_TRIALS = 100
SIM_PARAMS = [
    # varying lambda
    Parameters(prob_pred=0),
    Parameters(prob_pred=0.02),
    Parameters(prob_pred=0.05),
    Parameters(prob_pred=0.1),
    Parameters(prob_pred=0.2),
    Parameters(prob_pred=0.4),
    Parameters(prob_pred=0.6),
    Parameters(prob_pred=0.8),
    Parameters(prob_pred=1),
    Parameters(prob_pred=1.2),
    Parameters(prob_pred=1.4),
    Parameters(prob_pred=1.6),
    Parameters(prob_pred=1.8),
    Parameters(prob_pred=2),
    # varying max group size, lambda = 0.2
    Parameters(prob_pred=0.2, max_group_size=5),
    Parameters(prob_pred=0.2, max_group_size=10),
    Parameters(prob_pred=0.2, max_group_size=15),
    Parameters(prob_pred=0.2, max_group_size=20),
    Parameters(prob_pred=0.2, max_group_size=30),
    Parameters(prob_pred=0.2, max_group_size=40),
    Parameters(prob_pred=0.2, max_group_size=50),
    Parameters(prob_pred=0.2, max_group_size=100),
    # varying max group size, lambda = 2
    Parameters(prob_pred=2, max_group_size=5),
    Parameters(prob_pred=2, max_group_size=10),
    Parameters(prob_pred=2, max_group_size=15),
    Parameters(prob_pred=2, max_group_size=20),
    Parameters(prob_pred=2, max_group_size=30),
    Parameters(prob_pred=2, max_group_size=40),
    Parameters(prob_pred=2, max_group_size=50),
    Parameters(prob_pred=2, max_group_size=100),
    # varying lambda and e gain
    Parameters(prob_pred=0.2, e_gain=0.5),
    Parameters(prob_pred=0.2, e_gain=1.5),
    Parameters(prob_pred=0.2, e_gain=2),
    Parameters(prob_pred=2, e_gain=0.5),
    Parameters(prob_pred=2, e_gain=1.5),
    Parameters(prob_pred=2, e_gain=2),
    # increase maxf
    Parameters(prob_pred=0, maxf=1000),
    Parameters(prob_pred=0.4, maxf=1000),
    Parameters(prob_pred=0.8, maxf=1000),
    Parameters(prob_pred=1.2, maxf=1000),
    Parameters(prob_pred=1.6, maxf=1000),
    Parameters(prob_pred=2, maxf=1000),
    # parameter sensitivity testing, low e_gain
    Parameters(prob_pred=0, e_gain=0.5),
    Parameters(prob_pred=0.4, e_gain=0.5),
    Parameters(prob_pred=0.8, e_gain=0.5),
    Parameters(prob_pred=1.2, e_gain=0.5),
    Parameters(prob_pred=1.6, e_gain=0.5),
    Parameters(prob_pred=2, e_gain=0.5),
    # parameter sensitivity testing, high e_gain
    Parameters(prob_pred=0, e_gain=2.0),
    Parameters(prob_pred=0.4, e_gain=2.0),
    Parameters(prob_pred=0.8, e_gain=2.0),
    Parameters(prob_pred=1.2, e_gain=2.0),
    Parameters(prob_pred=1.6, e_gain=2.0),
    Parameters(prob_pred=2, e_gain=2.0),
    # parameter sensitivity testing, low max_group_size
    Parameters(prob_pred=0, max_group_size=10),
    Parameters(prob_pred=0.4, max_group_size=10),
    Parameters(prob_pred=0.8, max_group_size=10),
    Parameters(prob_pred=1.2, max_group_size=10),
    Parameters(prob_pred=1.6, max_group_size=10),
    Parameters(prob_pred=2, max_group_size=10),
    # parameter sensitivity testing, high max_group_size
    Parameters(prob_pred=0, max_group_size=50),
    Parameters(prob_pred=0.4, max_group_size=50),
    Parameters(prob_pred=0.8, max_group_size=50),
    Parameters(prob_pred=1.2, max_group_size=50),
    Parameters(prob_pred=1.6, max_group_size=50),
    Parameters(prob_pred=2, max_group_size=50),
    # parameter sensitivity testing, low coef_false
    Parameters(prob_pred=0, coef_false=0.1),
    Parameters(prob_pred=0.4, coef_false=0.1),
    Parameters(prob_pred=0.8, coef_false=0.1),
    Parameters(prob_pred=1.2, coef_false=0.1),
    Parameters(prob_pred=1.6, coef_false=0.1),
    Parameters(prob_pred=2, coef_false=0.1),
    # parameter sensitivity testing, high coef_false
    Parameters(prob_pred=0, coef_false=0.3),
    Parameters(prob_pred=0.4, coef_false=0.3),
    Parameters(prob_pred=0.8, coef_false=0.3),
    Parameters(prob_pred=1.2, coef_false=0.3),
    Parameters(prob_pred=1.6, coef_false=0.3),
    Parameters(prob_pred=2, coef_false=0.3),
]


def run_si_evolution_sims():
    start = time.time()
    with multiprocessing.Manager() as manager:
        p = multiprocessing.Pool(multiprocessing.cpu_count())
        args = []
        for params in SIM_PARAMS:
            for _ in range(NUM_TRIALS):
                args.append((OUT_FILE_DIR, params))

        p.starmap(evo_fun, args)
        p.close()
        p.join()

    end = time.time()
    print(
        f"Total amount of time to run {len(SIM_PARAMS) * NUM_TRIALS} simulations: {datetime.timedelta(seconds = end - start)}"
    )


if __name__ == "__main__":
    run_si_evolution_sims()
