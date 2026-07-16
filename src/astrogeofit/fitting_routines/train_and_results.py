"""
Copyright (C) 2025  CNRS
Lead author : J. Laskar
    
This file is part of the AstroGeoFit software.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""


import time
import random
import sys
import numpy as np  # type: ignore
import logging

from functools import partial
from typing import Optional, Union
from itertools import cycle
from tqdm import tqdm
from multiprocessing import Pool, Manager
from threading import Thread
from joblib import Parallel, delayed 


from pymoo.optimize import minimize  # type:ignore
from pymoo.operators.mutation.pm import PM  # type: ignore
from pymoo.algorithms.moo.nsga2 import NSGA2  # type: ignore
from pymoo.operators.crossover.sbx import SBX  # type: ignore
from pymoo.termination import get_termination  # type:ignore
from pymoo.operators.sampling.rnd import FloatRandomSampling  # type:ignore

from astrogeofit.fitting_routines.problems import (
    Callback_getF,
    ParameterInference
)

logger = logging.getLogger("ToolLogger")

# =================================================== #
# ============= TRAINING FUNCTIONS ================== #
# =================================================== #


def progress_listener(q, total):
    use_progress_bar = sys.stdout.isatty()
    if use_progress_bar:
        from tqdm import tqdm
        with tqdm(total=total, desc="Processing Fitting Solutions", unit="solution") as pbar:
            for _ in range(total):
                q.get()
                pbar.update(1)
    else:
        prev_time = time.time()
        for i in range(total):
            q.get()
            current_time = time.time()
            elapsed = current_time - prev_time
            prev_time = current_time
            logger.info(f"Completed {i+1}/{total} solutions - time for last solution: {elapsed:.2f} sec")


def creation_dictionary(
    interpolator,
    inverse_SR_lims: list,
    freqs_model,
    range_freq_limits: list,
    n_pieces: int,
    metric,
    population_size: int,
):
    dict_kwargs = {}
    dict_kwargs["interpolator"] = interpolator
    dict_kwargs["inverse_SR_lims"] = inverse_SR_lims
    dict_kwargs["freqs_model"] = freqs_model
    dict_kwargs["free_freqs_lims"] = range_freq_limits
    dict_kwargs["n_pieces"] = n_pieces
    dict_kwargs["metric"] = metric
    dict_kwargs["callback"] = Callback_getF
    def_algorithm = partial(
        NSGA2,
        pop_size=population_size,
        crossover=SBX(prob=0.9, eta=15),
        mutation=PM(eta=20),
        eliminate_duplicates=True,
    )
    dict_kwargs["def_algorithm"] = def_algorithm
    return dict_kwargs

def __calculate_significance_test(
    dict_kwargs:dict,
    num_solution:int = 0,
):
    dict_kwargs["num_sol"] = num_solution
    optimization_results = algorithm_training(**dict_kwargs)
        
    return result_sort(
        num_solutions=1,
        list_num_genes=dict_kwargs["list_num_genes"],
        num_generations=dict_kwargs["number_generations"],
        optimization_results=[optimization_results],
    )

def obtain_significance_results(
    data: list,
    num_parallel_jobs: int,
    list_num_genes: list,
    num_solutions: int,
    num_generations: int,
    seed: Optional[Union[int, list]],
    ar_series,
    dict_kwargs: dict,
    remote:bool
) -> tuple:
    """
    Processes and sorts optimization results based on the best values obtained.

    Args:
        data (list): Input data, where data[0] is depth and data[1] is the series.
        num_parallel_jobs (int): Number of parallel jobs to run.
        list_num_genes (list): Number of genes to optimize at each level.
        num_solutions (int): Number of solutions to compute.
        num_generations (int): Number of generations for the algorithm.
        seed (int or list): Seed or list of seeds for reproducibility.
        significance_test (bool): Whether to use series from ar_series or just use data.
        ar_series: List of series (if significance_test=True).
        dict_kwargs (dict): Additional arguments including 'def_algorithm', 'def_problem', etc.

    Returns:
        tuple: (best_values_per_job_sorted, metric_sorted_results)
    """
    # Build depth structure
    list_depth_genes = [
        np.linspace(data[0][0], data[0][-1], n_genes) for n_genes in list_num_genes
    ]

    # Generate seed list
    if isinstance(seed, int):
        seeds = [seed + 1234 * i for i in range(num_solutions)]
    elif isinstance(seed, list):
        if len(seed) != num_solutions:
            seeds = [s for s, _ in zip(cycle(seed), range(num_solutions))]
            
        else:
            seeds = seed
    else:
        seeds = [random.randint(0, 2**32 - 1) for _ in range(num_solutions)]
    # Add required argument
    dict_kwargs["def_problem"] = ParameterInference

    # Setup queue for progress tracking
    if not remote:
        manager = Manager()
        queue = manager.Queue()
    else:
        queue = None
    
    args = dict(
            list_num_genes=list_num_genes,
            list_depth_genes=list_depth_genes,
            number_generations=num_generations,
            seeds=seeds,
            data=data,
            queue=queue,
            **dict_kwargs
        )
    if ar_series:
        arguments = []
        for ar_serie in ar_series:
            data = [args["data"][0], ar_serie]
            args["data"] = data
            arguments.append(args)
        if not remote:
            # Run in parallel with Pool
            progress_thread = Thread(target=progress_listener, args=(queue, len(ar_series)))
            progress_thread.start()
            with Pool(processes=num_parallel_jobs) as pool:
                results_async = [pool.apply_async(__calculate_significance_test, kwds={"dict_kwargs":arguments,"num_solution": i}) for i,arguments in enumerate(arguments)]
                data_results = [r.get() for r in results_async]
            progress_thread.join()
        else:
            data_results = Parallel(n_jobs=num_parallel_jobs, verbose=100)\
                        (delayed(__calculate_significance_test)(dict_kwargs=argument, num_solution = i)
                        for i, argument in enumerate(arguments))

    else:
        if not remote:
            progress_thread = Thread(target=progress_listener, args=(queue, num_solutions))
            progress_thread.start()
            # Run in parallel with Pool
            with Pool(processes=num_parallel_jobs) as pool:
                results_async = [pool.apply_async(__calculate_significance_test, kwds={"dict_kwargs":args, "num_solution": i}) for i in range(num_solutions)]
                data_results = [r.get() for r in results_async]
            progress_thread.join()
        else:
            data_results = Parallel(n_jobs=num_parallel_jobs, verbose=100)\
                        (delayed(__calculate_significance_test)(args, i)
                        for i in range(num_solutions))
    print("")
    return data_results
    

def obtain_results(
    data: list,
    num_parallel_jobs: int,
    list_num_genes: list,
    num_solutions: int,
    num_generations: int,
    seed: Optional[Union[int, list]],
    dict_kwargs: dict,
    remote:bool
) -> tuple:
    """
    Processes and sorts optimization results based on the best values obtained.

    Args:
        data (list): Input data, where data[0] is depth and data[1] is the series.
        num_parallel_jobs (int): Number of parallel jobs to run.
        list_num_genes (list): Number of genes to optimize at each level.
        num_solutions (int): Number of solutions to compute.
        num_generations (int): Number of generations for the algorithm.
        seed (int or list): Seed or list of seeds for reproducibility.
        significance_test (bool): Whether to use series from ar_series or just use data.
        ar_series: List of series (if significance_test=True).
        dict_kwargs (dict): Additional arguments including 'def_algorithm', 'def_problem', etc.

    Returns:
        tuple: (best_values_per_job_sorted, metric_sorted_results)
    """

    # Build depth structure
    list_depth_genes = [
        np.linspace(data[0][0], data[0][-1], n_genes) for n_genes in list_num_genes
    ]

    # Generate seed list
    if isinstance(seed, int):
        seeds = [seed + 1234 * i for i in range(num_solutions)]
    elif isinstance(seed, list):
        if len(seed) != num_solutions:
            seeds = [s for s, _ in zip(cycle(seed), range(num_solutions))]
        else:
            seeds = seed
    else:
        seeds = [random.randint(0, 2**32 - 1) for _ in range(num_solutions)]
    # Add required argument
    dict_kwargs["def_problem"] = ParameterInference

    # Setup queue for progress tracking
    if not remote:
        manager = Manager()
        queue = manager.Queue()
    else:
        queue = None

    # Prepare args for multiprocessing
    args_list = []
    
    for sol in range(num_solutions):
        args = dict(
            list_num_genes=list_num_genes,
            list_depth_genes=list_depth_genes,
            num_sol=sol,
            number_generations=num_generations,
            seeds=seeds,
            data=data,
            queue=queue,
            **dict_kwargs
        )
        args_list.append(args)
    if not remote:
        progress_thread = Thread(target=progress_listener, args=(queue, num_solutions))
        progress_thread.start()

        # Run in parallel with Pool
        with Pool(processes=num_parallel_jobs) as pool:
            results_async = [pool.apply_async(algorithm_training, kwds=args) for args in args_list]
            optimization_results = [r.get() for r in results_async]
        
        progress_thread.join()
    else:
        optimization_results = Parallel(n_jobs=num_parallel_jobs, verbose=100)\
                        (delayed(algorithm_training)(**args)
                        for args in args_list)
    print("")
    return result_sort(
        num_solutions=num_solutions,
        list_num_genes=list_num_genes,
        num_generations=num_generations,
        optimization_results=optimization_results,
    )


def algorithm_training(
    list_num_genes: list,
    list_depth_genes: list,
    def_algorithm,
    num_sol: int,
    number_generations: int,
    seeds: list,
    def_problem,
    data,
    callback=Callback_getF(),
    queue=None,
    **kwargs,
):
    """
    like train_sequentially_v2, use list_depth_genes + algorithm + callback as input.
    """
    termination = get_termination("n_gen", number_generations)
    opitmization_result = []
    kwargs["data"] = data
    interpolator = kwargs["interpolator"]
    pop_size = def_algorithm.keywords["pop_size"]
    if seeds and seeds != []:
        seed = seeds[num_sol - 1]
        np.random.seed(seed)

    else:
        seed = None
    for i, depth_genes in enumerate(list_depth_genes):
        problem = def_problem(depth_genes, **kwargs)
        n_ff = problem.n_free_freqs

        if i == 0:
            X_init = FloatRandomSampling()(problem, pop_size).get("X")
        else:
            resp = opitmization_result[i - 1]
            args_best = [np.argmin(resp.F.mean(axis=1))] + [
                np.argmin(resp.F[:, i]) for i in range(resp.F.shape[1])
            ]
            X_init_prev = np.array(
                [
                    interpolator(
                        [resp.problem.depth_genes, resp.X[j, n_ff:]], depth_genes
                    )
                    for j in args_best
                ]
            )
            X_init_prev = np.hstack([resp.X[args_best, :n_ff], X_init_prev])
            X_init_rand = FloatRandomSampling()(
                problem, pop_size - len(X_init_prev)
            ).get("X")
            X_init = np.vstack([X_init_prev, X_init_rand])

        algorithm = def_algorithm(sampling=X_init)

        if seed:
            res = minimize(
                problem,
                algorithm,
                termination,
                seed=seed,
                save_history=False,
                verbose=False,
                callback=callback(),
            )
        else:
            res = minimize(
                problem,
                algorithm,
                termination,
                save_history=False,
                verbose=False,
                callback=callback(),
            )
        del res.problem.Data
        opitmization_result.append(res)
    if queue is not None:
        queue.put(1)
    return opitmization_result


def result_sort(
    num_solutions: int,
    list_num_genes: list,
    num_generations: int,
    optimization_results: list,
):
    metric_sorted_resuts = np.zeros(
        (num_solutions, len(list_num_genes), num_generations)
    )

    # Populate the array with negative best values from optimization results
    for i in range(num_solutions):
        for j in range(len(list_num_genes)):
            res = optimization_results[i][j]
            metric_sorted_resuts[i, j, :] = -np.array(
                res.algorithm.callback.data["best"]
            )

    # Sort indices based on the last generation's best values
    sorted_job_indices = np.argsort(metric_sorted_resuts[:, -1, -1])

    # Reorder optimization results based on sorted indices
    optimization_results = [optimization_results[i] for i in sorted_job_indices]

    # Return the sorted best values and sorted optimization results
    return metric_sorted_resuts[sorted_job_indices], optimization_results
