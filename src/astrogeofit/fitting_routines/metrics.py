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

from sklearn.metrics import r2_score  # type:ignore
from sklearn.linear_model import LinearRegression  # type:ignore

from astrogeofit.fitting_routines import misc
from astrogeofit.utils import shared_functions, exceptions


def taner_filter(x, dt, fc, fl, fh, c=10**1, npad=10, more_output=False):
    """
    %            Produces a bandpassed version of an input real-
    %            valued time series by FFT, multiplication by filter
    %            designed in frequency domain, and inverse FFT .
    %            Has zero phase response .
    %
    %    INPUTS:
    %           x = 1 column vector of data to filter
    %           dt = sample rate of x
    %           fc = center frequency of passband
    %           f1 = lower cut-off frequency
    %           fh = upper cut-off frequency
    %           c = roll-off/octave (default=1) (to avoid overflow error)
    %           tanerbandx = output filtered series

    """
    twopi = 2 * np.pi
    nlen = len(x)
    nptspow = 2 ** np.ceil(np.log2(nlen))
    newnlen = npad * int(nptspow)
    df = 1 / (newnlen * dt)

    # Filter setup
    dw = twopi / (newnlen * dt)
    wl = twopi * fl
    wc = twopi * fc
    wh = twopi * fh
    bw = wh - wl
    amp = 1 / np.sqrt(2)
    arg1 = 1 - (c * np.log(10)) / (20 * np.log(amp))
    arg1 = np.log(arg1)
    arg2 = ((bw + 2) / bw) ** 2
    arg2 = np.log(arg2)
    twod = 2 * arg1 / arg2
    # Array initialization for filter and frequency vector
    filter = np.zeros(newnlen)

    f = np.arange(newnlen // 2 + 1) * df

    w = np.arange(newnlen // 2 + 1) * dw
    # Compute filter coefficients for positive frequencies
    arg = (2 * np.abs(w - wc) / bw) ** twod
    filter[: len(w)] = amp * np.exp(-arg)

    filter[newnlen // 2 + 1 :] = filter[1 : newnlen // 2][
        ::-1
    ]  # Mirror the filter for negative frequencies
    filtout = filter * np.sqrt(2)  # Normalization for plotting

    # Forward FFT and filtering
    xftfor = sp.fft.fft(x, newnlen)
    xftfor *= filter
    xftinv = sp.fft.ifft(xftfor)

    tanerbandx = np.sqrt(2) * np.real(xftinv[:nlen])  # Remove padding

    if more_output:
        return (
            tanerbandx,
            filtout[: newnlen // 2 + 1],
            f,
        )  # Adjust array slice for f and filtout
    else:
        return tanerbandx


# ==================================== #
# ========= METRIC FUNCTIONS ========= #
# ==================================== #


def metric_piecewise_envelope_tanner(
    depth_invSR,
    data: list,
    fs_env,
    window_filter,
    interpolator,
    n_pieces=1,
    inverse_SR_lims=None,
    more_output=False,
):
    """_summary_:

        metric  (according to the metric type) for the linear model with predictors variable are fourier harmonics of frequencies fs
        of fitting the data, when used the age model derived from the invSR.
        age model: given points of inverse SR, interpolate with an interpolator and integrate along depth
        to obtain time.
    Args:
        depth_invSR (array[2,n]): inverse of sedimetation rates and their corresponding depth
        data (array[2,m]): the data and corresponding depth
        fs_env (_type_):[ecc] list of frequencies to model the envelope
        window_filter: window in the spectrum to define a bandpass filter of data.
        Returns:
        int: metric

    """

    depth, y_data = data

    time = misc.invSR_to_time(depth_invSR, depth, interpolator, inverse_SR_lims)

    time_grid = np.linspace(time[0], time[-1], len(y_data))
    dt = time_grid[1] - time_grid[0]
    f_sampling = 1 / dt
    # window_filter = window_filter.copy()
    fl, fh = window_filter
    if fh >= f_sampling / 2 - 1e-5:
        # Does not satisfy the nyquist condition
        return np.zeros(n_pieces)
    # fh = min(f_sampling/2-1e-5, window_filter[1])

    fc = (fl + fh) / 2
    y_filtered_tanner = taner_filter(y_data, dt, fc, fl, fh, c=10)
    y_complex = sp.signal.hilbert(y_filtered_tanner)
    y_env = abs(y_complex)

    X_env = shared_functions.generate_sine_waves(
        np.ones_like(fs_env), fs_env, time_grid
    )

    reg_model = LinearRegression().fit(X_env, y_env)
    y_pred_env = reg_model.predict(X_env)

    depth_pieces = np.linspace(depth[0], depth[-1], n_pieces + 1)
    metrix_env = np.zeros(n_pieces)
    for i in range(n_pieces):
        j1, j2 = np.searchsorted(depth, depth_pieces[i], "left"), np.searchsorted(
            depth, depth_pieces[i + 1], "right"
        )
        metrix_env[i] = r2_score(y_env[j1:j2], y_pred_env[j1:j2])
    if more_output:
        return metrix_env, y_env, y_pred_env
    else:
        return metrix_env


def metric_piecewise(
    depth_invSR,
    data: list,
    fs,
    interpolator,
    n_pieces: int = 1,
    inverse_SR_lims: list = None,
    lags: int = 2,
    metric_type: str = "r2",
):
    """_summary_:

        metric  (according to the metric type) for the linear model with predictors variable are fourier harmonics of frequencies fs
        of fitting the data, when used the age model derived from the invSR.
        age model: given points of inverse SR, interpolate with an interpolator and integrate along depth
        to obtain time.
    Args:
        depth_invSR (array[2,n]): inverse of sedimetation rates and their corresponding depth
        interpolator: interpolator
        data (array[2,m]): the data and corresponding depth
        fs (_type_): list of frequencies of the model
        metric_type (str, optional):"BIC", "AIC" or "RSS". Defaults to "BIC".
            Returns:
        int: metric

    """

    depth, y_data = data
    depth_genes, _ = depth_invSR
    invSR_interpolate = interpolator(depth_invSR, depth)

    if inverse_SR_lims is not None:
        invSR_lo = np.zeros(len(depth_genes)) + inverse_SR_lims[0]
        invSR_up = np.zeros(len(depth_genes)) + inverse_SR_lims[1]
        invSR_interpolate_lo = interpolator([depth_genes, invSR_lo], depth)
        invSR_interpolate_up = interpolator([depth_genes, invSR_up], depth)
        inds_lo, inds_up = (
            invSR_interpolate < invSR_interpolate_lo,
            invSR_interpolate > invSR_interpolate_up,
        )
        invSR_interpolate[inds_lo] = invSR_interpolate_lo[inds_lo]
        invSR_interpolate[inds_up] = invSR_interpolate_up[inds_up]

    time = sp.integrate.cumulative_trapezoid(invSR_interpolate, depth, initial=0)
    X = shared_functions.generate_X_linReg(np.ones_like(fs), fs, time)

    reg_model = LinearRegression().fit(X, y_data)
    y_pred = reg_model.predict(X)

    if metric_type == "r2":
        return r2_score_piecewise(y_data, y_pred, n_pieces)
    elif metric_type == "loglike":
        return shared_functions.calculate_log_likelihood_ar_piecewise(
            y_data - y_pred, lags, n_pieces
        )
    else:
        raise exceptions.WrongMetric(metric=metric_type)


def r2_score_piecewise(y_data, y_pred, k):
    n = len(y_data)
    segment_size = n // k
    r2s = np.zeros(k)
    for i in range(k):
        start_idx = i * segment_size
        end_idx = (i + 1) * segment_size if i < k - 1 else n
        r2s[i] = r2_score(y_data[start_idx:end_idx], y_pred[start_idx:end_idx])
    return r2s
