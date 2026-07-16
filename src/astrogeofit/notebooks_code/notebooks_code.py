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

import dill
import logging
import numpy as np

from typing import Optional
from scipy.special import logsumexp

from astrogeofit.main_routines import data_manager
from astrogeofit.plots.plots_post_mcmc import plots_post_mcmc
from astrogeofit.utils import shared_functions, setup_logger

logger = logging.getLogger("ToolLogger")
# Check if the logger already has handlers
if not logger.hasHandlers():
    logger = setup_logger.setup_logger(False)


def obtain_significance_test_file(
    config_file_path: Optional[str],
    output_folder: Optional[str],
    output_file: Optional[str],
):
    if config_file_path:
        return data_manager.read_yaml_file(config_file_path)
    else:
        return data_manager.read_yaml_file(
            f"{output_folder}/analysis_parameters_{output_file}.yaml"
        )


def compute_weights(log_posteriors):
    max_log_posterior = np.max(log_posteriors)
    scaled_log_posteriors = log_posteriors - max_log_posterior

    log_normalizer = logsumexp(scaled_log_posteriors)

    weights = np.exp(scaled_log_posteriors - log_normalizer)
    return weights


def data_obtain_analysis_file_and_result_files(
    obtain_mcmc_results: bool,
    config_file_path: Optional[str],
    output_folder: Optional[str],
    output_file: Optional[str],
    fitting_path_file:Optional[str]= None,
):
    """
    Load analysis configuration and result files required for further computations and plotting.

    This function reads in a configuration YAML file and corresponding `.pickle` result files
    (fitted results, optional MCMC results, and optional model weights). It supports either specifying
    the configuration file path directly or constructing paths from a given output folder and output file name.
    
    If a seed is defined in the configuration, the plotting module's random seed is initialized for reproducibility.

    Parameters:
        - *obtain_mcmc_results* : bool
            Whether to load the MCMC result file (`True`) or not (`False`).
        - *config_file_path* : str or None
            Path to the configuration YAML file. If `None`, both `output_folder` and `output_file` must be provided.
        - *output_folder* : str or None
            Folder containing the output `.pickle` files. Used if `config_file_path` is not provided.
        - *output_file* : str or None
            Base name of the output `.pickle` files. Used if `config_file_path` is not provided.

    Returns:
        - *results_fit* : dict or None
            Dictionary containing the fitted results loaded from `<output_file>_fit.pickle`.
        - *results_mcmc* : dict or None
            Dictionary with MCMC results loaded from `<output_file>_mcmc.pickle` if `obtain_mcmc_results=True`; otherwise `None`.
        - *weight_results* : list or np.ndarray or None
            List or array of computed weights (if `executed_weights=True` in MCMC results), otherwise empty list.
        - *config_analysis* : dict
            Dictionary with the full content of the analysis configuration file.
        - *data_set_values* : dict
            Dictionary of dataset-related parameters from the configuration file.
        - *age_depth_model_values* : dict or None
            Dictionary containing parameters related to age-depth modeling.
        - *genetic_algorithm_parameters* : dict or None
            Parameters for the genetic algorithm used during model fitting.
        - **data_model_parameters** : dict
            Parameters defining the structure of the data model.
        - *mcmc_parameters* : dict or None
            Parameters used for the MCMC procedure.
        - *eccentricity_solution_parameters* : dict or None
            Parameters related to astronomical eccentricity solution modeling.

    Raises:
        Exception
            If neither a configuration file path nor both `output_folder` and `output_file` are provided.
    """
    if not config_file_path and (not output_folder and not output_file):
        logger.error(
            "THERE IS NO VALUE INTRODUCED FOR THE CONFIGURATION FILE AND NO VALUES FOR THE OUTPUT FOLDER AND OUTPUT FILE NAMES. PLEASE INTRODUCE THE CONFIGURATION FILE OR THE FOLDER AND NAME OF THE RESULTS."
        )
        raise Exception
    
    if config_file_path:
        print(config_file_path)
        if config_file_path == "odp_926" or config_file_path == "ODP_926":
            config_analysis = "./examples/ODP_926/configuration_file/configuration_ODP_926.yml"
        elif config_file_path == "odp_1260" or config_file_path == "ODP_1260":
            config_analysis = "./examples/ODP_1260/configuration_file/configuration_ODP_1260.yml"
        elif config_file_path == "odp_1262" or config_file_path == "ODP_1262":
            config_analysis = "./examples/ODP_1262/configuration_file/configuration_ODP_1262.yml"
        elif config_file_path == "basic_example":
            config_file_path = "./examples/synthetic_data/configuration_file/configuration_synthetic_data.yml"
        
        config_analysis = data_manager.read_yaml_file(config_file_path)
    else:
        config_analysis = data_manager.read_yaml_file(
            f"{output_folder}/analysis_parameters_{output_file}.yaml"
        )
    if output_folder == "default" or output_folder == "" or not output_folder:
        output_folder = config_analysis["output_folder"]
    if output_file == "default"  or output_file == "" or not output_file:
        output_file = config_analysis["output_file_name"]
    if obtain_mcmc_results and fitting_path_file:
        path_results_file_fit = fitting_path_file
    else:
        path_results_file_fit = f"{output_folder}/{output_file}_fit.pickle"
    with open(path_results_file_fit, "rb") as f:
        results_fit = dill.load(f)
    if obtain_mcmc_results:
        path_results_file_mcmc = f"{output_folder}/{output_file}_mcmc.pickle"

        with open(path_results_file_mcmc, "rb") as f:
            results_mcmc = dill.load(f)
    else:
        results_mcmc = None

    data_set_values = config_analysis["data_set"]
    if "age_depth_model_data" in config_analysis:
        age_depth_model_values = config_analysis["age_depth_model_data"]
    else:
        age_depth_model_values = None
    if "genetic_algorithm_parameters" in config_analysis:
        genetic_algorithm_parameters = config_analysis["genetic_algorithm_parameters"]
    else:
        genetic_algorithm_parameters = None
    data_model_parameters = config_analysis["data_model_parameters"]
    if "mcmc_parameters" in config_analysis:
        mcmc_parameters = config_analysis["mcmc_parameters"]
    else:
        mcmc_parameters = None
    if "eccentricity_solution_data" in config_analysis:
        eccentricity_solution_parameters = config_analysis["eccentricity_solution_data"]
    else:
        eccentricity_solution_parameters = None
    if genetic_algorithm_parameters["seed"]:
        if isinstance(genetic_algorithm_parameters["seed"], int):
            plots_post_mcmc.setup_seed(genetic_algorithm_parameters["seed"])
        elif isinstance(genetic_algorithm_parameters["seed"], list):
            plots_post_mcmc.setup_seed(
                genetic_algorithm_parameters["seed"][0]
            )  # WE CHOOSE THE FIRST SEED OF THE LIST ARBITRARILY
    
    if obtain_mcmc_results and results_mcmc["executed_weights"]:
        weight_path_results = f"{output_folder}/{output_file}_weights.pickle"
        with open(weight_path_results, "rb") as f:
            weight_results = dill.load(f)
        weights_list = []
        evidences = weight_results["evidences"]
        for k in np.arange(len(evidences)):
            weights_list.append([])
            for i in np.arange(len(evidences[k])):
                evidence_slice = np.array(evidences)[k, i : i + 1]
                average_evidence = logsumexp(evidence_slice, axis=0) - np.log(
                    len(evidences)
                )
                weights_list[-1].append(compute_weights(average_evidence))

        weights_list = np.array(weights_list)
        weight_results = weights_list.mean(axis=1)
    else:
        weight_results = None

    return (
        results_fit,
        results_mcmc,
        weight_results,
        config_analysis,
        data_set_values,
        age_depth_model_values,
        genetic_algorithm_parameters,
        data_model_parameters,
        mcmc_parameters,
        eccentricity_solution_parameters,
    )


def data_obtain_age_depth_values(age_depth_model_values: dict):
    """
    Load and process age-depth model data from a file if a path is provided.

    This function reads age-depth model data from the path specified in `age_depth_model_values`
    and computes nominal sedimentation rate (inverse) and age-depth functions using the shared utility functions.
    If no data path is provided, all returned values are set to `None`, and no nominal model is added.

    Parameters:
        - **age_depth_model_values**:

    Returns:
        - **depth_age_model_data** (dict or None): Parsed age-depth model data from the file, or `None` if no path is provided.
        - **depth_age_model_df** (pandas.DataFrame or None): DataFrame containing age-depth model data, or `None` if no path is provided.
        - **invSR_nominal** (array-like or None): Inverse sedimentation rate profile (nominal model), or `None` if no path is provided.
        - **func_invSR_nominal** (callable or None): Function to evaluate inverse sedimentation rate over depth, or `None` if no path is provided.
        - **func_time_nominal** (callable or None): Function to convert depth to age using the nominal model, or `None` if no path is provided.
        - **add_nominal** (bool): Flag indicating whether the nominal model should be plotted. `True` if data was loaded, `False` otherwise.
    """
    if isinstance(age_depth_model_values,dict) and age_depth_model_values["data_path"]:
        depth_age_model_data, depth_age_model_df = data_manager.data_read_dataset_from_file(
            age_depth_model_values,
            False,
            True,
            False
        )
        (
            invSR_nominal,
            func_invSR_nominal,
            func_time_nominal,
        ) = shared_functions.get_age_depth_model_funcs(depth_age_model_data)
        add_nominal = (
            True  # CHANGE THIS LINE IF YOU DON'T WANT TO SEE THE NOMINAL DATA PLOTTED.
        )
    else:
        depth_age_model_data = None
        depth_age_model_df = None
        invSR_nominal = None
        func_invSR_nominal = None
        func_time_nominal = None
        add_nominal = False

    return (
        depth_age_model_data,
        depth_age_model_df,
        invSR_nominal,
        func_invSR_nominal,
        func_time_nominal,
        add_nominal,
    )