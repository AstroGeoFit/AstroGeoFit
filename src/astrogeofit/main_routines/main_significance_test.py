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
import random

from statsmodels.tsa.ar_model import AutoReg
from statsmodels.tsa.arima_process import arma_acovf

from astrogeofit.main_routines.main_fitting import main_fitting_function
from astrogeofit.significance_test import utils as significance_test_utils

def setup_seed(seed: int):
    """Set the global seed for reproducibility."""
    random.seed(seed)
    np.random.seed(seed)

def get_shape(obj):
    print(type(obj))
    if isinstance(obj, (list, tuple)):
        if len(obj) == 0:
            return (0,)  # Handle empty lists/tuples
        return (len(obj),) + get_shape(obj[1])  # Check the first element recursively
    return ()  # Base case: if it's not a list/tuple, return empty tuple


def significance_test_execution(
    data_model_parameters: dict,
    significance_test_parameters: dict,
    data_set: list,
    remote:bool
) -> dict:
    seed = significance_test_parameters["seed"]
    if isinstance(seed, int):
        setup_seed(seed)
    elif isinstance(seed, list):
        setup_seed(seed[0])
        
    results = main_fitting_function(
        data_model_parameters=data_model_parameters,
        configuration_parameters=significance_test_parameters,
        data_set=data_set,
        significance_test=True,
        remote=remote
    )
    model = AutoReg(data_set[1], lags=2, old_names=False)
    fitted_model = model.fit()

    phi = fitted_model.params[1:]
    intercept = fitted_model.params[0]  # This is not used. Why is it here.
    sigma = np.std(fitted_model.resid)
    v = arma_acovf(
        ar=np.hstack([1, -phi]), ma=[1], nobs=len(data_set[1]), sigma2=sigma**2
    )
    V = sp.linalg.toeplitz(v)
    L = np.linalg.cholesky(V)
    best_result = max(results, key=lambda x: x[0][:,-1,-1])[1][0][-1]
    best_result_pred = significance_test_utils.get_prediction(best_result, data_set[0], data_set[1])['y_pred']
    residuals = data_set[1] - best_result_pred
    r2_data = []
    for _, res in enumerate(results):
        r2_data.append(
            [
                significance_test_utils.get_r2(res[1][0][j], data_set[1], data_set[0], L)
                for j in range(len(res[1][0]))
            ]
        )
    r2_data = np.array(r2_data)

    dict_to_save = {"r2_data": r2_data}

    sample_detrended = lambda: significance_test_utils.detrend(
        data_set[0],
        np.array(
            significance_test_utils.sample(
                data_set[1], fitted_model.params, np.std(fitted_model.resid, ddof=1)
            )
        )
    )
    print("We started detrending the noise data.")
    print("")
    ar_series = [
        sample_detrended()
        for _ in range(significance_test_parameters["number_noise_executions"])
    ]
    print("Finished detrending the noise data.")
    print("")
    results_ar = main_fitting_function(
        data_model_parameters=data_model_parameters,
        configuration_parameters=significance_test_parameters,
        data_set=data_set,
        significance_test=True,
        ar_series=ar_series,
        remote=remote
    )
    r2_ar = []
    for i, res in enumerate(results_ar):
        r2_ar.append(
            [
                significance_test_utils.get_r2(res[1][0][j], ar_series[i], data_set[0], L)
                for j in range(len(res[1][0]))
            ]
        )
    r2_ar = np.array(r2_ar)
    dict_to_save["r2_ar"] = r2_ar

    return dict_to_save
