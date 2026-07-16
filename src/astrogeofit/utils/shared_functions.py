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

import scipy as sp
import numpy as np
import pandas as pd
from scipy.linalg import lstsq  # type:ignore
from scipy.interpolate import interp1d  # type:ignore
from scipy.interpolate import (
    CubicSpline,
    splev,
    splrep,
    PchipInterpolator,
    Akima1DInterpolator,
    interp1d,
)
from sklearn.linear_model import LinearRegression  # type:ignore

from astrogeofit.main_routines import data_manager


def interpolate_CubicSpline(SR, depth):
    if len(SR[0]) == 1:
        SR_interpolate = np.ones_like(depth) * SR[1]
    else:
        spl_SR = CubicSpline(*SR)
        SR_interpolate = spl_SR(depth)
    return SR_interpolate


def interpolate_Pchip(xy, xs):
    spl_xy = PchipInterpolator(*xy)
    return spl_xy(xs)


def interpolate_BSpline(SR, depth):
    spl_SR = splrep(*SR)
    SR_interpolate = splev(depth, spl_SR)
    return SR_interpolate


def interpolate_Akima(SR, depth):
    if len(SR[0]) == 2:
        slope = (SR[1][1] - SR[1][0]) / (SR[0][1] - SR[0][0])
        SR_interpolate = slope * depth + SR[1][0]
    else:
        SR_interpolate = Akima1DInterpolator(SR[0], SR[1])(depth)
    return SR_interpolate


def interpolate_interp1d(SR, depth, *args, **kwargs):
    f = interp1d(*SR, *args, **kwargs)
    return f(depth)


def get_eccentricity_parameters(eccentricity_params) -> CubicSpline:
    e3_data, _ = data_manager.data_read_dataset_from_file(eccentricity_params, False, False, True)
    e3_data[0] = -e3_data[0]
    t0, tf = float(eccentricity_params["start_time"]), float(
        eccentricity_params["final_time"]
    )
    inds = (e3_data[0] >= t0) & (e3_data[0] <= tf)
    time_e3, e3 = e3_data[:, inds]
    e3 = (e3 - e3.mean()) / e3.std()
    return CubicSpline(time_e3 / 1e3, e3)  # CHANGE: NEED TO CHANGE THe 1/e3 VALUE


def generate_X_linReg(Amp, freq, times):
    ts_cos = Amp[None] * np.sin(freq[None] * times[:, None])
    ts_sin = Amp[None] * np.cos(freq[None] * times[:, None])
    return np.concatenate([ts_cos, ts_sin], axis=1)


def generate_sine_waves(Amp, freq, times, phi=0):
    Amp, freq, times = list(map(np.array, [Amp, freq, times]))
    ts_cos = Amp[None] * np.sin(freq[None] * times[:, None] + phi)
    ts_sin = Amp[None] * np.cos(freq[None] * times[:, None] + phi)
    return np.concatenate([ts_cos, ts_sin], axis=1)


def invSR_to_prediction(
    depth_invSR, data, inverse_SR_lims, fs, interpolator=interpolate_CubicSpline
):
    depth, y = data
    invSR_interpolate = interpolator(depth_invSR, depth)
    invSR_interpolate[invSR_interpolate < inverse_SR_lims[0]] = inverse_SR_lims[0]
    invSR_interpolate[invSR_interpolate > inverse_SR_lims[1]] = inverse_SR_lims[1]
    time = sp.integrate.cumulative_trapezoid(invSR_interpolate, depth, initial=0)
    X = generate_X_linReg(np.ones_like(fs), fs, time)
    reg_model = LinearRegression().fit(X, y)
    y_pred = reg_model.predict(X)
    return {"time": time, "y_pred": y_pred, "reg_model": reg_model}


# USED TO CALCULATE THE FUNCTIONS OF THE AGE DEPTH MODEL
def get_age_depth_model_funcs(depth_age_model_data):
    invSR_nominal = np.gradient(depth_age_model_data[1], depth_age_model_data[0])
    func_invSR_nominal = interp1d(
        depth_age_model_data[0], invSR_nominal, kind="nearest", fill_value="extrapolate"
    )
    func_time_nominal = interp1d(
        depth_age_model_data[0], depth_age_model_data[1], fill_value="extrapolate"
    )
    return invSR_nominal, func_invSR_nominal, func_time_nominal


def calculate_log_likelihood_ar_piecewise(time_series, p, k=1):
    def fit_ar_model(time_series, p):
        """
        Fits an autoregressive model of lag p to the given time series.

        Parameters:
        time_series (list or np.array): 1D time series data
        p (int): lag order of the autoregressive model

        Returns:
        tuple: coefficients of the AR model, intercept term, and residuals
        """
        time_series = np.asarray(time_series)
        n = len(time_series)

        if p >= n:
            raise ValueError(
                "The lag order p must be less than the length of the time series."
            )

        # Create the lagged matrix X and the response vector y
        X = np.lib.stride_tricks.sliding_window_view(time_series, p + 1)
        y = X[:, -1]
        X = X[:, :-1]

        # Add a column of ones for the intercept term
        X = np.column_stack([np.ones(X.shape[0]), X])

        # Solve the linear regression problem X * beta = y
        beta, _, _, _ = lstsq(X, y)

        # Separate the intercept term and AR coefficients
        intercept = beta[0]
        ar_coefficients = beta[1:]

        # Calculate the fitted values
        fitted_values = X @ np.concatenate(([intercept], ar_coefficients))

        # Calculate the residuals
        residuals = y - fitted_values

        return ar_coefficients, intercept, residuals

    """
    Calculates the log likelihood of the fitted AR model of lag p for the given time series.

    Parameters:
    time_series (list or np.array): 1D time series data
    p (int): lag order of the autoregressive model

    Returns:
    float: log likelihood of the fitted AR model
    """
    ar_coefficients, intercept, residuals = fit_ar_model(time_series, p)

    n = len(time_series)
    T = n - p

    # Calculate the variance of the residuals
    segment_size = n // k
    rss_k = np.zeros(k)
    for i in range(k):
        start_idx = i * segment_size
        end_idx = (i + 1) * segment_size if i < k - 1 else n
        rss_k[i] = np.sum(residuals[start_idx:end_idx] ** 2)

    rss = np.sum(residuals**2)
    sigma2 = rss / T

    C = -0.5 * T * (np.log(2 * np.pi) + np.log(sigma2))
    # Calculate the log likelihood
    log_likelihood = C / k - 0.5 / sigma2 * rss_k
    if k == 1:
        log_likelihood = log_likelihood[0]

    return log_likelihood


# USED IN PRE-PROCESS
def interpolate_data_to_max_points(depths, values, max_points, step=None):
    """
    Interpolate data to get a specified maximum number of evenly spaced points.

    Parameters:
    depths (list or numpy array): The original depths.
    values (list or numpy array): The original values associated with each depth.
    max_points (int): The maximum number of points to keep in the interpolated data.

    Returns:
    numpy array: Interpolated depths with a maximum of max_points.
    numpy array: Interpolated values corresponding to the interpolated depths.
    """
    # Convert depths and values to numpy arrays if they are not already
    depths = np.array(depths)
    values = np.array(values)

    # Calculate the smallest difference between consecutive depths
    if step == None:
        min_diff = np.min(np.diff(depths))
    else:
        min_diff = step

    # Generate evenly spaced depths based on min_diff
    even_depths = np.arange(depths.min(), depths.max() + min_diff, min_diff)

    # Create an interpolation function
    interp_func = interp1d(depths, values, kind="linear", fill_value="extrapolate")

    # Interpolate the values for the evenly spaced depths
    even_values = interp_func(even_depths)

    # Determine the total number of interpolated points
    total_points = len(even_depths)

    # Check if we need to reduce the number of points
    if total_points > max_points:
        # Calculate the step to select evenly spaced points
        step = total_points / max_points
        # Select indices for evenly spaced points
        indices = np.round(np.arange(0, total_points, step)).astype(int)
        # Use these indices to select evenly spaced points
        even_depths = even_depths[indices]
        even_values = even_values[indices]
    interp_data = (even_depths, even_values)
    data = {"Depth (interpolated)": even_depths, "proxy (interpolated)": even_values}
    return interp_data, pd.DataFrame(data)


# USED FOR FFT PLOTS. IN PRE-PROCESS AND TO VISULAIZE THE RESULTS
def count_points_in_meter_range(depths, start_meter=None):
    """
    Count the number of points in the depths array that fall between start_meter and start_meter + 1.

    Parameters:
    depths (list or numpy array): The original depths.
    start_meter (float): The starting meter value. If None, it will take the first value found in depths.

    Returns:
    int: The number of points between start_meter and start_meter + 1.
    """
    if start_meter == None:
        start_meter = depths[0]

    # Calculate the end meter value
    end_meter = start_meter + 1.0

    # Filter the depths to include only those within the specified range
    points_in_range = (depths >= start_meter) & (depths < end_meter)

    # Count the number of points that fall within the range
    count = np.sum(points_in_range)

    return count
