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

import random
import numpy as np
import scipy as sp

from functools import partial

from astrogeofit.weight_calcula_routines import weight_calcula
from astrogeofit.mcmc_routines import utils

def setup_seed(seed: int):
    """Set the global seed for reproducibility."""
    random.seed(seed)
    np.random.seed(seed)

def execute_weight_computation(
    optimization_results,
    data,
    samples_MCMC,
    num_weight_executions: int,
    num_parallel_jobs_weights: int,
    prior_frequencies,
    prior_distributions,
    stability_factor: float,
    pareto_smoothing: bool,
    seed
):
    if isinstance(seed, int):
        setup_seed(seed)
    elif isinstance(seed, list):
        setup_seed(seed[0])
        
    logllh_seq = []
    for k in range(len(samples_MCMC)):
        problem = optimization_results[0][k].problem
        invSR_lims = problem.inverse_SR_lims
        depth_genes = problem.depth_genes
        freqs_model = problem.freqs_model  # !!!
        interpolator = problem.interpolator  # !!!
        lags = 2
        logllh = partial(
            utils.get_log_likelihood,
            depth_of_invSR=depth_genes,
            freqs_model=freqs_model,
            inverse_SR_lims=invSR_lims,
            data=data,
            interpolator=interpolator,
            lags=lags,
        )
        logllh_seq.append(logllh)

    ## GET THE PRIOR
    log_prior = partial(
        utils.get_log_prior,
        inverse_SR_lims=invSR_lims,
        log_prior_freq=prior_distributions,
        log_prior_freq_params=prior_frequencies,
    )
    evidences = []
    for k, sample_MCMC in enumerate(samples_MCMC):
        evidences.append([])
        N_seq = sample_MCMC.shape[0]
        for i in range(num_weight_executions):
            evidences[-1].append([])
            for i in range(N_seq):
                posterior_samples = sample_MCMC[i]
                posterior_m_samples = [
                    sp.stats.norm.rvs(
                        loc=posterior_samples[:, j].mean(),
                        scale=posterior_samples[:, j].std(axis=0, ddof=1),
                        size=posterior_samples.shape[0],
                    )
                    for j in range(posterior_samples.shape[1])
                ]
                posterior_m_samples = np.column_stack(posterior_m_samples)
                posterior_m_pdf = weight_calcula.get_m_gaussian_fit(
                    posterior_samples, stability_factor, posterior_m_samples
                )
                posterior_m_pdf = np.log(posterior_m_pdf).sum(axis=-1)
                log_evidence = weight_calcula.perakis(
                    posterior_m_samples,
                    posterior_m_pdf,
                    logllh_seq[k],
                    log_prior,
                    n_jobs=num_parallel_jobs_weights,
                    pareto_smoothing=pareto_smoothing,
                    show_progress=False,
                )
                evidences[-1][-1].append(log_evidence)
    return np.array(evidences)
