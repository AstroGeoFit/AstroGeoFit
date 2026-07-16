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
import logging

from colorama import Fore, Style  # type: ignore
from functools import partial

from astrogeofit.utils import shared_functions
from astrogeofit.fitting_routines import metrics, frequency_model, train_and_results

logger = logging.getLogger("ToolLogger")

'''
def test_results(name_file_results:str, **kwargs):
    """
    Writes a given content to a text file. If the file already exists, it appends the new content.

    Parameters:
    - file_path (str): The path to the text file.
    - content (str): The content to write to the file.
    """
    with open(f"test_results_{name_file_results}.txt", 'a') as file:
        for key, value in kwargs.items():
            file.write(f"{key}: {value}\n")
        file.write('\n')  # Add a newline for separation between entries
'''


def chose_interpolator(chosen_interpolator: str):
    if chosen_interpolator == "default" or chose_interpolator == "CubicSpline":
        return shared_functions.interpolate_CubicSpline
    elif chosen_interpolator == "Pchip":
        return shared_functions.interpolate_Pchip
    elif chosen_interpolator == "BSpline":
        return shared_functions.interpolate_BSpline
    elif chosen_interpolator == "AKIMA":
        return shared_functions.interpolate_Akima
    else:
        return shared_functions.interpolate_interp1d


def main_fitting_function(
    data_model_parameters: dict,
    configuration_parameters: dict,
    data_set: list,
    significance_test: bool = False,
    ar_series=None,
    remote:bool = False
) -> dict:
    values_of_fitting: dict = {}
    convers_factor: float = (180 * 3600) / (np.pi * 1000000)
    values_of_fitting["convers_factor"] = convers_factor
    inverse_SR_lims: list = [
        0.1 / data_model_parameters["sedimentation_rate_max"],
        0.1 / data_model_parameters["sedimentation_rate_min"],
    ]
    values_of_fitting["inverse_SR_lims"] = inverse_SR_lims

    # ---------------- FREQUENCY MODEL --------------------- #
    logger.info(f"{Fore.GREEN} Calculating frequency model. {Style.RESET_ALL}")
    print("")
    default_values: dict = {
        "p0_values": 50.467718,
        "gi_values": [5.579378, 7.456665, 17.366595, 17.910194],
        "g5_value": 4.257564,
        "si_values": [-18.845166, -17.758310],
        "s6_value": -26.347880
    }
    range_freq_limits, freqs_model = frequency_model.frequency_calcula(
        frequency_parameters=data_model_parameters["frequency_values"],
        default_values=default_values,
        convers_factor=convers_factor,
    ) 
    # ---------------- GENETIC ALGORITHM VALUES RETRIEVAL -------------------- #
    if (
        not configuration_parameters["metric_type"]
        or configuration_parameters["metric_type"] == ""
    ):
        metric_type = "loglike"
    else:
        metric_type = configuration_parameters["metric_type"]
    interpolator = chose_interpolator(configuration_parameters["interpolator"])
    metric = partial(
        metrics.metric_piecewise, metric_type=metric_type, lags=2
    )  # METRIC TYPE SHOULD BE POSSIBLE TO BE CHANGED
    dict_kwargs = train_and_results.creation_dictionary(
        interpolator=interpolator,
        inverse_SR_lims=inverse_SR_lims,
        freqs_model=freqs_model,
        range_freq_limits=range_freq_limits,
        n_pieces=3,
        metric=metric,
        population_size=configuration_parameters["population_size"],
    )

    print("")
    if significance_test:
        if not ar_series:
            num_executions = configuration_parameters["number_algorithm_executions"]
            logger.info(
                f"{Fore.GREEN} All parameters obatined, starting algorithm significance test. {Style.RESET_ALL}"
            )
            print("")
        else:
            logger.info(
                f"{Fore.GREEN} All parameters obatined, starting significance test. {Style.RESET_ALL}"
            )
            print("")
            num_executions = configuration_parameters["number_noise_executions"]
        significance_test_results = train_and_results.obtain_significance_results(
            data=data_set,
            num_parallel_jobs=configuration_parameters["number_processors_used"],
            list_num_genes=configuration_parameters["list_number_genes"],
            num_solutions=num_executions,
            num_generations=configuration_parameters["number_generations"],
            seed=configuration_parameters["seed"],
            ar_series=ar_series,
            dict_kwargs=dict_kwargs,
            remote = remote
        )
        return significance_test_results
        
    else:
        logger.info(
            f"{Fore.GREEN} All parameters obatined, starting analysis. {Style.RESET_ALL}"
        )
        print("")
        num_solutions = configuration_parameters["number_algorithm_solutions"]
        metric_sorted_resuts, optimization_results = train_and_results.obtain_results(
            data=data_set,
            num_parallel_jobs=configuration_parameters["number_processors_used"],
            list_num_genes=configuration_parameters["list_number_genes"],
            num_solutions=num_solutions,
            num_generations=configuration_parameters["number_generations"],
            seed=configuration_parameters["seed"],
            dict_kwargs=dict_kwargs,
            remote = remote
        )
        # WE CHOOSE THE LAST SOLUTION WITH THE LAST NUMBER OF GENES (BEST RESULT)

        values_of_fitting["metric_sorted_resuts"] = metric_sorted_resuts
        values_of_fitting["optimization_results"] = optimization_results
        return values_of_fitting
