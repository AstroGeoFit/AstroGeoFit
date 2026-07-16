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

import sys
import time
import random
import numpy as np  # type:ignore
import emcee  # type:ignore
import logging

from joblib import Parallel, delayed  # type:ignore
from functools import partial
from itertools import cycle
from tqdm import tqdm
from multiprocessing import Manager
from threading import Thread


from astrogeofit.mcmc_routines import utils

logger = logging.getLogger("ToolLogger")

def progress_listener(q, total):
    use_progress_bar = sys.stdout.isatty()
    if use_progress_bar:
        from tqdm import tqdm
        with tqdm(total=total, desc="Processing MCMC Solutions", unit="solution") as pbar:
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


def setup_seed(seed: int):
    """Set the global seed for reproducibility."""
    random.seed(seed)
    np.random.seed(seed)


def sampler_wrapper(
    params_ini,
    log_PDF,
    n_step,
    execution,
    seeds,
    params_ini_lims=None,
    nwalkers=None,
    scale=1e-3,
    ndim_min=40,
    
):
    #logger.info(f"Starting job number {execution} of the MCMC")
    """
    define a emcee.EnsembleSampler to sample log_PDF from an ensemble centered
    epsilon-closely around  params_ini for n_step
    """
    setup_seed(seeds[execution])
    ndim = len(params_ini)
    scale_initial = np.ones(ndim) * scale
    if nwalkers is None:
        nwalkers = max(ndim * 2, ndim_min)
    if params_ini_lims is None:
        pos_ini = params_ini + scale_initial[None] * np.random.randn(nwalkers, ndim)
    else:
        pos_ini = np.zeros((nwalkers, ndim))
        for i in range(ndim):
            lo = max(params_ini[i] - scale_initial[i], params_ini_lims[i][0])
            up = min(params_ini[i] + scale_initial[i], params_ini_lims[i][1])
            pos_ini[:, i] = np.random.uniform(lo, up, nwalkers)
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_PDF)
    #logger.info(f"Running MCMC for job number {execution}")
    sampler.run_mcmc(pos_ini, n_step)
    #logger.info(f"Finished MCMC for job number {execution}")
    return sampler


def execution_mcmc(
    optimization_results: list,
    prior_frequencies: list,
    prior_distributions: list,
    num_solutions: int,
    list_mcmc_genes: list,
    list_num_genes: list,
    data: list,
    num_parallel_jobs: int,
    n_step,
    seed,
    remote:bool,
) -> dict:
    """
    Perform MCMC sampling with time tracking for each stage.
    """
    res = optimization_results[-1][-1]
    freqs_model = res.problem.freqs_model
    inverse_SR_lims = res.problem.inverse_SR_lims
    interpolator = res.problem.interpolator
    lags = 2

    seeds = []
    if isinstance(seed, int):
        seeds = [seed + 1234 * i for i in range(0, num_solutions)]
    elif isinstance(seed, list):
        if len(seed) < num_solutions:
            seeds = [seed for _, x in zip(range(num_solutions), cycle(seeds))]
        seeds = seed
    else:
        seeds = [random.randint(0, 2**32 - 1) for _ in range(num_solutions)]
    
    log_prior = partial(
        utils.get_log_prior,
        inverse_SR_lims=inverse_SR_lims,
        log_prior_freq=prior_distributions,
        log_prior_freq_params=prior_frequencies,
    )

    ## UNTIL HERE EVERYTHING IS THE SAME
    params_ini_par = []
    params_ini_lims_par = []
    scale_par = []
    log_PDF_par = []

    # Time tracking for parameter setup
    start_time_setup = time.time()

    # Number of genes. Create an option to choose the number of genes that u want to calcultate
    # the MCMC from.
    for j, num_genes_to_analyze in enumerate(list_mcmc_genes):
        index_num_gen_analyze = list_num_genes.index(num_genes_to_analyze)
        ## ADD: Should be added a error control in case that the number of genes does not exist.
        params_ini_par.append([])

        # We take the last num_solutions as they are the best (ordered before.)
        # The range from 1 to num_solutions +1 is like this bc we are taking the last solution: -1 to the -num_solutions
        for i in range(1, num_solutions + 1):
            res = optimization_results[-i][index_num_gen_analyze]
            arg_best = np.argmin(res.F.sum(axis=1))
            params_ini_par[j].append(res.X[arg_best])
        res = optimization_results[0][index_num_gen_analyze]
        num_free_freqs = res.problem.n_free_freqs
        depth_genes = res.problem.depth_genes
        invSR_lims_x = np.tile(np.array(inverse_SR_lims), (len(depth_genes), 1))
        params_ini_lims = np.vstack([res.problem.free_freqs_lims, invSR_lims_x])
        params_ini_lims_par.append(params_ini_lims)
        # ASK: WHEN DOES THE FACTOR 1e3 CHANGES
        scale = (
            params_ini_lims[:, 1] - params_ini_lims[:, 0]
        ) * 1e-3  # !!!!!! HARDCODED VALUE, NOT NECESSARY TO CHANGE IT BY THE MOMENT!
        # THIS IS TO GET A GOOD INITIAL VARIANCE OF THE g and s FREQUENCIES
        if num_free_freqs > 1:
            # ASK: WHEN DOES THE FACTOR 10 CHANGES
            scale[1:num_free_freqs] = (
                scale[1:num_free_freqs] * 10
            )  # !!!!!!  HARDCODED VALUE, NOT NECESSARY TO CHANGE IT BY THE MOMENT!
        scale_par.append(scale)

        logllh = partial(
            utils.get_log_likelihood,
            depth_of_invSR=depth_genes,
            freqs_model=freqs_model,
            inverse_SR_lims=inverse_SR_lims,
            data=data,
            interpolator=interpolator,
            lags=lags,
        )
        log_PDF = partial(
            utils.get_log_posterior, log_prior=log_prior, log_likelihood=logllh
        )
        log_PDF_par.append(log_PDF)

    setup_time = time.time() - start_time_setup
    logger.info(f"MCMC parameters setup completed in {setup_time:.2f} seconds")
    print("")

    sampler_par = []

    # Time tracking for parallel sampling

    # Number of genes. Create an option to choose the number of genes that u want to calcultate
    # the MCMC from.
    for j, _ in enumerate(list_mcmc_genes):
    #logger.info(f"Starting sampling for the solutions with {mcmc_genes} number of genes.")

    # Setup progress tracking
        if not remote:
            manager = Manager()
            queue = manager.Queue()
            
            progress_thread = Thread(target=progress_listener, args=(queue, num_solutions))
            progress_thread.start()
            
            # Parallel execution with queue signaling
            def wrapper_with_queue(*args, **kwargs):
                res = sampler_wrapper(*args, **kwargs)
                queue.put(1)
                return res

            sampler_par_j = Parallel(n_jobs=num_parallel_jobs, verbose=100)(
                delayed(wrapper_with_queue)(
                    params_ini=params_ini_par[j][i],
                    log_PDF=log_PDF_par[j],
                    n_step=n_step,
                    scale=scale_par[j],
                    params_ini_lims=params_ini_lims_par[j],
                    ndim_min=50,
                    execution=i,
                    seeds = seeds,
                )
                for i in range(num_solutions)
            )
            print("")
            progress_thread.join()
            sampler_par.append(sampler_par_j)

        else:
            sampler_par_j = Parallel(n_jobs=num_parallel_jobs, verbose=100)(
                delayed(sampler_wrapper)(
                    params_ini=params_ini_par[j][i],
                    log_PDF=log_PDF_par[j],
                    n_step=n_step,
                    scale=scale_par[j],
                    params_ini_lims=params_ini_lims_par[j],
                    ndim_min=50,
                    execution=i,
                    seeds = seeds,
                )
                for i in range(num_solutions)
            )
            print("")
            sampler_par.append(sampler_par_j)
        
    return sampler_par