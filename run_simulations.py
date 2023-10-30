import multiprocessing
import time
import datetime
import math
from multiprocessing import Pool
from si_evolution import evo_fun
from analyze import get_sim_id
from constants import DEFAULT_PARAMS, out_file_path


def run_sim_batch(start_sim_id: int, num_sims: int):
    params = DEFAULT_PARAMS
    for i in range(start_sim_id, start_sim_id + num_sims):
        evo_fun(out_file_path, params, sim_id=i)


def run_si_evolution_sims(num_simulations: int):
    cpu_count = multiprocessing.cpu_count()
    p = Pool(cpu_count)
    args = []

    start_sim_id = get_sim_id(f"{out_file_path}/all.csv")

    batch_size = math.ceil(num_simulations / cpu_count)
    num_simulations % cpu_count
    for i in range(cpu_count):
        if i + 1 > num_simulations:
            break

        curr_batch_size = (
            batch_size if i < num_simulations % cpu_count else batch_size - 1
        )
        args.append((start_sim_id, curr_batch_size))
        start_sim_id += curr_batch_size
    p.starmap(run_sim_batch, args)
    p.close()
    p.join()


if __name__ == "__main__":
    now = datetime.datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print(f"Starting simulations at {current_time}")
    start = time.time()

    num_simulations = 100
    run_si_evolution_sims(num_simulations)

    end = time.time()

    print(
        f"Total amount of time to run {num_simulations} simulations: {datetime.timedelta(seconds = end - start)}"
    )
