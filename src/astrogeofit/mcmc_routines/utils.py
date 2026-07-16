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
import numpy as np
import logging

from astrogeofit.utils import shared_functions

logger = logging.getLogger("ToolLogger")


def calculate_prior_frequencies_and_distributions(prior_frquencies_values: dict, prior_distributions_values:dict, frequency_values: dict, convers_factor:float):
    prior_frequencies = []
    prior_distributions = []
    if (
        prior_frquencies_values["p0_prior"] != []
        and prior_frquencies_values["p0_prior"] and 
        (frequency_values["use_precession"] or frequency_values["use_tilt"])
    ):
        prior_frequencies.append(prior_frquencies_values["p0_prior"])
        distribution = prior_distributions_values["p0_distribution"]
        prior_distributions.extend([distribution] * 1)

    if (
        prior_frquencies_values["gi_prior"] != []
        and prior_frquencies_values["gi_prior"] and 
        (frequency_values["use_eccentricity"] or frequency_values["use_precession"])
    ):
        prior_frequencies.extend(prior_frquencies_values["gi_prior"])
        distribution = prior_distributions_values["gi_distribution"]
        prior_distributions.extend([distribution] * 4)

    if (
        prior_frquencies_values["si_prior"] != []
        and prior_frquencies_values["si_prior"]
        and frequency_values["use_tilt"]
    ):
        prior_frequencies.extend(prior_frquencies_values["si_prior"])
        distribution = prior_distributions_values["si_distribution"]
        prior_distributions.extend([distribution] * 2)

    if prior_frequencies == []:
        prior_frequencies = None
    else:
        prior_frequencies = (
            np.stack(prior_frequencies) / convers_factor
        )
    return prior_frequencies, prior_distributions


def get_mcmc_samplers_and_logprob(sampler_par, discard=0, thin=1):
    """
    Retrieve MCMC samples and log probabilities with time tracking.
    """
    start_time_total = time.time()  # Track total time

    # Time tracking for retrieving samples
    start_time_samples = time.time()
    samples_combine_par = [
        np.stack(
            [
                sampler.get_chain(discard=discard, thin=thin, flat=True)
                for sampler in sampler_par[j]
            ]
        )
        for j in range(len(sampler_par))
    ]
    samples_time = time.time() - start_time_samples
    logger.info(f"Samples retrieval completed in {samples_time:.2f} seconds")

    # Time tracking for retrieving log probabilities
    start_time_logprob = time.time()
    logprob_par = [
        np.stack(
            [
                sampler.get_log_prob(discard=discard, thin=thin, flat=True)
                for sampler in sampler_par[j]
            ]
        )
        for j in range(len(sampler_par))
    ]
    logprob_time = time.time() - start_time_logprob
    logger.info(f"Log probabilities retrieval completed in {logprob_time:.2f} seconds")

    # Compute total time
    total_time = time.time() - start_time_total
    logger.info(
        f"Total time for MCMC samplers and log probability retrieval: {total_time:.2f} seconds"
    )

    return samples_combine_par, logprob_par


def get_log_likelihood(
    params, depth_of_invSR, freqs_model, inverse_SR_lims, interpolator, data, lags
):
    """wrapper of log_likelihood_empirical_bayes_arima"""
    number_of_genes = len(depth_of_invSR)
    n_free_freqs = len(params) - number_of_genes
    _, y_data = data
    freqs = params[:n_free_freqs]
    invSR = params[n_free_freqs:]
    depth_invSR = [depth_of_invSR, invSR]
    fs = freqs_model(freqs)
    y_pred = shared_functions.invSR_to_prediction(
        depth_invSR, data, inverse_SR_lims, fs, interpolator=interpolator
    )["y_pred"]
    noise = y_data - y_pred

    llh = shared_functions.calculate_log_likelihood_ar_piecewise(noise, lags)
    return llh


# ATENTION, I CHANGED THIS FUNCTIONS. MIGHT NOT WORK.
def get_log_prior(
    params, inverse_SR_lims: list, log_prior_freq, log_prior_freq_params: list
):
    # define log prior

    def log_uniform(x, x_lims):
        if np.all((x > x_lims[0]) & (x < x_lims[1])):
            return 0.0
        return -np.inf

    def log_gaussian(x, params):
        muy, sigma = params
        return -0.5 * (x - muy) ** 2 / sigma**2

    N_freq = len(log_prior_freq)
    freqs = params[:N_freq]
    invSR = params[N_freq:]
    lp_invSR = log_uniform(invSR, inverse_SR_lims)
    lp_freqs = []
    for i, log_freq in enumerate(log_prior_freq):
        if log_freq == "uniform":
            result = log_uniform(freqs[i], log_prior_freq_params[i])
        elif log_freq == "gaussian":
            result = log_gaussian(freqs[i], log_prior_freq_params[i])

        lp_freqs.append(result)
    lp_freqs_sum = sum(lp_freqs)
    return lp_invSR + lp_freqs_sum


def get_log_posterior(params, log_prior, log_likelihood):
    """
    define log posterior
    """
    lp = log_prior(params)
    if not np.isfinite(lp):
        return -np.inf
    llh = log_likelihood(params)
    return llh + lp
