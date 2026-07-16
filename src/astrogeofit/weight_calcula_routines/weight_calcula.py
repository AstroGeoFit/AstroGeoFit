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

import numpy as np
import scipy as sp

from tqdm import tqdm
from joblib import Parallel, delayed
from scipy.special import logsumexp


def get_m_posterior_kde(posterior_samples, kernel_scale=1.0):
    N, D = posterior_samples.shape
    stds = posterior_samples.std(axis=0, ddof=1)
    stds[stds == 0] = 1e-8
    scales = kernel_scale * N ** (-1 / 5) * stds

    kde = [
        sp.stats.gaussian_kde(posterior_samples[:, j], bw_method=scales[j])
        for j in range(D)
    ]
    pdf = np.column_stack([kde[j](posterior_samples[:, j]) for j in range(D)])

    return pdf


def get_m_gaussian_fit(posterior_samples, stability_factor: float, samples=None):
    N, D = posterior_samples.shape
    means = posterior_samples.mean(axis=0)
    stds = posterior_samples.std(axis=0, ddof=1)
    stds[stds == 0] = stability_factor

    if samples is None:
        samples = posterior_samples
    pdf = np.column_stack(
        [
            sp.stats.norm.pdf(samples[:, j], loc=means[j], scale=stds[j])
            for j in range(D)
        ]
    )

    return pdf


def shuffle(posterior_samples, posterior_m_pdf):
    shuffled_samples = posterior_samples.copy()
    shuffled_pdf = posterior_m_pdf.copy()

    for j in range(posterior_samples.shape[1]):
        ind = np.random.permutation(posterior_samples.shape[0])
        shuffled_samples[:, j] = posterior_samples[ind, j]
        shuffled_pdf[:, j] = posterior_m_pdf[ind, j]

    return shuffled_samples, shuffled_pdf


def pareto_loglike(x, theta):
    k = -np.mean(np.log(1 - theta * x))
    n = len(x)
    return n * (np.log(theta / k) + k - 1)


def parreto(x):
    n = len(x)
    x = np.sort(x)
    m = 20 + int(n**0.5)
    x_star = x[int(0.25 * n + 0.5)]
    theta = [1 / x[-1] + (1 - (m / (j - 0.5))) / (3 * x_star) for j in range(1, m + 1)]
    loglikes = np.array([pareto_loglike(x, theta_j) for theta_j in theta])
    w = [1 / np.exp(loglikes - loglikes[j]).sum() for j in range(m)]
    theta_hat = (np.array(theta) * np.array(w)).sum()
    k = -np.mean(np.log(1 - theta_hat * x))
    sigma = k / theta_hat
    return theta_hat, k, sigma


def transform_weights(r):
    M = int(min(0.2 * len(r), 3 * len(r) ** 0.5))
    r_sorted = np.sort(r)
    r_highest = r_sorted[-M:]
    u = r_highest[0]
    _, k, sigma = parreto(r_highest - u)
    new_highest = np.array(
        [u + (sigma / k) * (1 / (1 - (z - 0.5) / M) ** k - 1) for z in range(1, M + 1)]
    )
    max_r = r.max()
    new_highest[new_highest > max_r] = max_r
    # return np.concat((r_sorted[:-M], new_highest)), k
    return np.concatenate((r_sorted[:-M], new_highest)), k


def perakis(
    m_posterior_samples,
    m_posterior_log_pdf,
    loglikelihood_f,
    prior_f,
    pareto_smoothing: bool,
    show_progress=True,
    n_jobs=-1,
):
    l = Parallel(n_jobs=n_jobs)(
        delayed(loglikelihood_f)(sample)
        for sample in tqdm(m_posterior_samples, disable=not show_progress)
    )
    prior = Parallel(n_jobs=n_jobs)(
        delayed(prior_f)(sample) for sample in m_posterior_samples
    )

    l = np.array(l)
    prior = np.array(prior)

    r = np.exp(prior - m_posterior_log_pdf)
    new_r, k = transform_weights(r)
    if pareto_smoothing:
        r = new_r

    if k > 0.7:
        print(f"Warning, reporting k = {k}")

    normalize_factor = r.sum()
    log_weights = l + np.log(r) - np.log(normalize_factor)
    max_log = log_weights.max()

    return logsumexp(log_weights - max_log) + max_log
