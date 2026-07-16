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
from functools import partial

from pymoo.core.problem import ElementwiseProblem  # type:ignore
from pymoo.core.callback import Callback  # type:ignore

from astrogeofit.fitting_routines import metrics


class ParameterInference(ElementwiseProblem):
    """
    freqs_model: a function to get the frequencies of the signal from free variables as input
    or it can be a constant array, when n_free_freqs =0
    """

    def __init__(
        self,
        depth_genes,
        inverse_SR_lims,
        interpolator,
        freqs_model,
        data,
        free_freqs_lims=None,
        n_pieces=1,
        metric=metrics.metric_piecewise,
    ):
        self.inverse_SR_lims = inverse_SR_lims
        self.depth_genes = depth_genes
        self.interpolator = interpolator
        self.Data = data
        self.n_free_freqs = len(free_freqs_lims) if free_freqs_lims is not None else 0
        self.freqs_model = freqs_model
        self.N_genes_invSR = len(depth_genes)
        self.n_pieces = n_pieces
        self.free_freqs_lims = free_freqs_lims

        if self.n_free_freqs == 0:
            #self.freqs_model = lambda x=None: freqs_model
            xl_free_freqs, xu_free_freqs = ([], [])
        else:
            self.freqs_model = freqs_model
            xl_free_freqs, xu_free_freqs = (
                free_freqs_lims[:, 0],
                free_freqs_lims[:, 1],
            )
        xl_invSR = np.zeros(self.N_genes_invSR) + inverse_SR_lims[0]
        xu_invSR = np.zeros(self.N_genes_invSR) + inverse_SR_lims[1]

        xl = np.hstack([xl_free_freqs, xl_invSR])
        xu = np.hstack([xu_free_freqs, xu_invSR])
        n_var = self.N_genes_invSR + self.n_free_freqs
        self.func_metric = partial(
            metric,
            interpolator=interpolator,
            n_pieces=n_pieces,
            inverse_SR_lims=inverse_SR_lims,
        )

        super().__init__(n_var=n_var, n_obj=n_pieces, xl=xl, xu=xu)

    def _evaluate(self, genes, out, *args, **kwargs):
        if self.n_free_freqs == 0:
            fs = self.freqs_model(None)
        else:
            fs = self.freqs_model(genes[: self.n_free_freqs])
        genes_invSR = genes[self.n_free_freqs :]
        out["F"] = -self.func_metric(
            [self.depth_genes, genes_invSR], fs=fs, data=self.Data
        )


class ParameterInference_env(ElementwiseProblem):
    """
    ParameterInference with envelope
    freqs_model: a function to get the frequencies of the signal from free variables as input
    or it can be a constant array, when n_free_freqs =0
    Note that in this version. first 5 frequencies must be 5 climatic precession frequencies (p0+g_i)
    """

    def __init__(
        self,
        depth_genes,
        inverse_SR_lims,
        interpolator,
        freqs_model,
        freqs_model_env,
        data,
        free_freqs_lims=None,
        n_pieces=1,
        metric=metrics.metric_piecewise,
        metric_env=metrics.metric_piecewise_envelope_tanner,
        df=5,
    ):
        self.inverse_SR_lims = inverse_SR_lims
        self.depth_genes = depth_genes
        self.interpolator = interpolator
        self.Data = data
        self.n_free_freqs = len(free_freqs_lims) if free_freqs_lims is not None else 0
        self.freqs_model = freqs_model
        self.number_of_genes = len(depth_genes)
        self.n_pieces = n_pieces
        self.free_freqs_lims = free_freqs_lims
        self.df = df

        if isinstance(freqs_model_env, (list, np.ndarray)):
            self.freqs_model_env = lambda x=None: np.array(freqs_model_env)

        # if self.n_free_freqs==0:
        if isinstance(freqs_model, (list, np.ndarray)):
            self.freqs_model = lambda x=None: np.array(freqs_model)
            xl_free_freqs, xu_free_freqs = ([], [])
        else:
            self.freqs_model = freqs_model
            xl_free_freqs, xu_free_freqs = (
                free_freqs_lims[:, 0],
                free_freqs_lims[:, 1],
            )  # type:ignore
        xl_invSR = np.zeros(self.number_of_genes) + inverse_SR_lims[0]
        xu_invSR = np.zeros(self.number_of_genes) + inverse_SR_lims[1]

        xl = np.hstack([xl_free_freqs, xl_invSR])
        xu = np.hstack([xu_free_freqs, xu_invSR])
        n_var = self.number_of_genes + self.n_free_freqs
        self.func_metric = partial(
            metric,
            interpolator=interpolator,
            n_pieces=n_pieces,
            inverse_SR_lims=inverse_SR_lims,
        )
        self.func_metric_env = partial(
            metric_env,
            interpolator=interpolator,
            n_pieces=n_pieces,
            inverse_SR_lims=inverse_SR_lims,
        )
        super().__init__(n_var=n_var, n_obj=n_pieces, xl=xl, xu=xu)

    def _evaluate(self, genes, out, *args, **kwargs):
        fs = self.freqs_model(genes[: self.n_free_freqs])
        fs_env = self.freqs_model_env(genes[: self.n_free_freqs])
        genes_invSR = genes[self.n_free_freqs :]
        window_filter = np.array(
            [min(fs[:5]) / 2 / np.pi - self.df, max(fs[:5]) / 2 / np.pi + self.df]
        )
        r2_env = self.func_metric_env(
            [self.depth_genes, genes_invSR],
            fs_env=fs_env,
            window_filter=window_filter,
            data=self.Data,
        )
        r2 = self.func_metric([self.depth_genes, genes_invSR], fs=fs, data=self.Data)
        out["F"] = -r2 * r2_env


class Callback_getF(Callback):
    def __init__(self) -> None:
        super().__init__()
        self.data["best"] = []

    def notify(self, algorithm):
        b = algorithm.pop.get("F").mean(axis=1).min()
        self.data["best"].append(b)
