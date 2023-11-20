import multiprocessing
import time
import datetime
import math
from multiprocessing import Pool
from si_evolution import evo_fun
from analyze import get_sim_id
from constants import out_file_path
from si_types import Parameters
from typing import Optional, List


def run_sim_batch(
    start_sim_id: int, num_sims: int, params: Parameters, start_time: Optional[float]
):
    for i in range(start_sim_id, start_sim_id + num_sims):
        start = time.time()
        evo_fun(out_file_path, params, sim_id=i)
        end = time.time()
        if i == start_sim_id and start_time is not None:
            sim_runtime = end - start
            total_runtime = start_time + (sim_runtime * num_sims)
            print(
                f"Estimated completion time: {time.strftime('%H:%M:%S', time.localtime(total_runtime))}"
            )


def run_si_evolution_sims(start_time: float, params: Parameters, num_simulations: int):
    cpu_count = multiprocessing.cpu_count()
    p = Pool(cpu_count)
    args = []

    start_sim_id = get_sim_id(f"{out_file_path}/all.csv")

    batch_size = math.ceil(num_simulations / cpu_count)
    for i in range(cpu_count):
        if i + 1 > num_simulations:
            break

        curr_batch_size = (
            batch_size
            if (num_simulations % cpu_count == 0 or i < num_simulations % cpu_count)
            else batch_size - 1
        )
        print_estimate = start_time if i == 0 else None
        args.append((start_sim_id, curr_batch_size, params, print_estimate))
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
    sim_params = []
    for i, params in enumerate(sim_params):
        param_start = time.time()
        run_si_evolution_sims(param_start, params, num_simulations)
        now = datetime.datetime.now()
        print(
            f"Finished {num_simulations} simulations for {(i+1)}/{len(sim_params)} parameters at {now.strftime('%H:%M:%S')}"
        )

    end = time.time()

    print(
        f"Total amount of time to run {len(sim_params) * num_simulations} simulations: {datetime.timedelta(seconds = end - start)}"
    )
