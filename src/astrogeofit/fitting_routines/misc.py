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

import scipy as sp  # type:ignore
import numpy as np  # type:ignore

from astrogeofit.utils import shared_functions


# ==================================== #
# ======= INV SED RATE TO PRED ======= #
# ==================================== #


## THIS FUNCTION IS MISSING THE depth_genes VARIABLE. WHY?


def invSR_to_time(depth_invSR, grid, interpolator, inverse_SR_lims=None):
    """integrate invSR to obtain an age model"""

    invSR_interpolate = interpolator(depth_invSR, grid)
    if inverse_SR_lims is not None:
        invSR_interpolate[invSR_interpolate < inverse_SR_lims[0]] = inverse_SR_lims[0]
        invSR_interpolate[invSR_interpolate > inverse_SR_lims[1]] = inverse_SR_lims[1]
    time = sp.integrate.cumulative_trapezoid(invSR_interpolate, grid, initial=0)
    return time


# ============================================= #
# ========= LOG LIKELY HOOD FUNCTIONS ========= #
# ============================================= #

# FUNCTIONS THAT ARE NOT USED SHOULD BE DELETED


def get_envelope_precession(
    params,
    depth_of_invSR,
    freqs_model,
    inverse_SR_lims,
    interpolator,
    data,
    ip,
    norm=True,
):
    n_genes = len(depth_of_invSR)
    n_ff = len(params) - n_genes
    invSR = params[n_ff:]
    freqs = params[:n_ff]
    fs = freqs_model(freqs)
    res = shared_functions.invSR_to_prediction(
        [depth_of_invSR, invSR], data, inverse_SR_lims, fs, interpolator
    )
    reg_model = res["reg_model"]
    times_inferred = res["time"]
    ip_b = ip + len(fs)

    AB_prec = np.hstack([reg_model.coef_[ip], reg_model.coef_[ip_b]])
    X_prec = shared_functions.generate_sine_waves(
        np.ones_like(fs[ip]), fs[ip], times_inferred
    )
    X_prec_90 = shared_functions.generate_sine_waves(
        np.ones_like(fs[ip]), fs[ip], times_inferred, phi=np.pi / 2
    )

    y_pred_prec = X_prec @ AB_prec
    y_pred_prec_90 = X_prec_90 @ AB_prec
    y_envelope = np.sqrt(y_pred_prec_90**2 + y_pred_prec**2)
    if norm:
        y_envelope = (y_envelope - y_envelope.mean()) / y_envelope.std()

    return times_inferred, y_envelope


def multi_pearsonr(x, y):
    xmean = x.mean(axis=1)
    ymean = y.mean()
    xm = x - xmean[:, None]
    ym = y - ymean
    normxm = np.linalg.norm(xm, axis=1)
    normym = np.linalg.norm(ym)
    return np.clip(np.dot(xm / normxm[:, None], ym / normym), -1.0, 1.0)


def get_t0_from_e3(
    t_e3_est, f_e3, t0_lims, n_res=1000, inds_omit=[0, 0], return_full=False
):
    times_est, e3_est = t_e3_est
    ib, di = inds_omit
    ie = len(e3_est) - di
    t0s = np.linspace(*t0_lims, n_res)
    e3s_true = np.array([f_e3(times_est + t0) for t0 in t0s])
    r_t0s = multi_pearsonr(e3s_true[:, ib:ie], e3_est[ib:ie])
    i_best = np.argmax(r_t0s)
    if return_full:
        return t0s, r_t0s
    else:
        return t0s[i_best], r_t0s[i_best]
