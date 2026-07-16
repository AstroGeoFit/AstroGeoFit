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
import statsmodels.api as sm  # type: ignore

import astrogeofit.utils.shared_functions as shared_functions


def get_prediction(res, depth, y):
    data = [depth, y]
    arg_best = res.F.mean(axis=1).argmin()
    params = res.X[arg_best]
    problem = res.problem
    depth_genes = problem.depth_genes
    interpolator = problem.interpolator
    invSR_lims = problem.inverse_SR_lims
    n_ff = problem.n_free_freqs
    fs = problem.freqs_model(params[:n_ff])
    res_pred = shared_functions.invSR_to_prediction(
        [depth_genes, params[n_ff:]], data, invSR_lims, fs, interpolator
    )
    return res_pred


def sample(y, params, sigma):
    intercept, phi = params[0], np.array(params[1:])
    lag = len(phi)
    n_samples = len(y)
    samples = np.empty(n_samples)
    samples[:lag] = 0.0

    for i in range(lag, n_samples):
        noise = np.random.normal(0, sigma)
        new_value = intercept + np.dot(phi, samples[i - lag : i][::-1]) + noise
        samples[i] = new_value

    return samples


def get_r2(res, y, depth, L):
    eps = np.linalg.solve(L, y)
    pred = get_prediction(res, depth, y)
    eps_res = np.linalg.solve(L, y - pred["y_pred"])
    r2 = 1 - (eps_res.T @ eps_res) / (eps.T @ eps)
    return r2


def detrend(depth, y):
    depth_ori = depth
    y_ori = y  # These variables are not used. Should I eliminate them?

    y_lowess = sm.nonparametric.lowess(y, depth, frac=0.1)[:, 1]
    y_detr = y - y_lowess
    return y_detr
