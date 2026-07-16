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

import os
import scipy as sp  # type: ignore
import numpy as np
import random

from typing import Optional
from datetime import datetime

from astrogeofit.utils import shared_functions
from astrogeofit.main_routines import data_manager
from astrogeofit.notebooks_code import notebooks_code

def setup_seed(seed: int):
    """Set the global seed for reproducibility."""
    random.seed(seed)
    np.random.seed(seed)

def print_and_save(message, print_time: bool = False, file_path="output.txt"):
    if print_time:
        timestamp = datetime.now().strftime(
            "[%Y-%m-%d %H:%M:%S]"
        )  # Format: [YYYY-MM-DD HH:MM:SS]
        full_message = f"{timestamp} {message}"
    else:
        full_message = message
    print(full_message)  # Print to screen
    with open(file_path, "a") as f:  # Append mode
        print(full_message, file=f)  # Write to file


def get_result_atributes(
    result,
) -> dict:
    atributes = {}
    atributes["arg_best"] = np.argmin(result.F.mean(axis=1))
    atributes["number_free_freqs"] = result.problem.n_free_freqs
    atributes["best_inverse_SR"] = result.X[
        atributes["arg_best"], atributes["number_free_freqs"] :
    ]  # TAKE CARE IF THIS IS CORRE
    atributes["depth_genes"] = result.problem.depth_genes
    atributes["frequencies"] = result.problem.freqs_model(
        result.X[atributes["arg_best"], : atributes["number_free_freqs"]]
    )
    atributes["interpolator"] = result.problem.interpolator
    atributes["inverse_SR_lims"] = result.problem.inverse_SR_lims
    atributes["frequency_model"] = result.problem.freqs_model
    return atributes


def get_ETP_components(
    params,
    depth_of_invSR,
    freqs_model,
    inverse_SR_lims,
    interpolator,
    data,
    ip=None,
    ie=None,
    it=None,
    env_prec=True,
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
    
    result = {}
    result["time"] = times_inferred
    
    if ip is not None: # ----YW, 2025-Apr-10
        ip_b = ip + len(fs)
    
        AB_prec = np.hstack([reg_model.coef_[ip], reg_model.coef_[ip_b]])
        X_prec = shared_functions.generate_sine_waves(np.ones_like(fs[ip]), fs[ip], times_inferred)
        y_pred_prec = X_prec @ AB_prec
        result["y_prec"] = y_pred_prec  # ----YW, 2025-Apr-10
        
        X_prec_90 = shared_functions.generate_sine_waves(
            np.ones_like(fs[ip]), fs[ip], times_inferred, phi=np.pi / 2
        )
        y_pred_prec_90 = X_prec_90 @ AB_prec
        y_envelope = np.sqrt(y_pred_prec_90**2 + y_pred_prec**2)
        result["y_env_prec"] = y_envelope
    if ie is not None:
        ie_b = ie + len(fs)
        AB_ecc = np.hstack([reg_model.coef_[ie], reg_model.coef_[ie_b]])
        X_ecc = shared_functions.generate_sine_waves(
            np.ones_like(fs[ie]), fs[ie], times_inferred
        )
        y_pred_ecc = X_ecc @ AB_ecc
        result["y_ecc"] = y_pred_ecc

    if it is not None:
        it_b = it + len(fs)
        AB_ecc = np.hstack([reg_model.coef_[it], reg_model.coef_[it_b]])
        X_ecc = shared_functions.generate_sine_waves(
            np.ones_like(fs[it]), fs[it], times_inferred
        )
        y_pred_ecc = X_ecc @ AB_ecc
        result["y_tilt"] = y_pred_ecc
    return result

def data_obtain_mcmc_results_per_a_number_of_knots(
    configuration_file_path:str,
    folder_with_the_results:Optional[str],
    file_with_the_results:Optional[str],
    selected_num_genes: Optional[int],
    positive_feedback: bool,
    use_prec_env: bool,
    ignore_weights: bool,
) -> None:
    """
    Extract and process MCMC results for a selected number of knots in an astrochronological model.

    This function loads analysis and MCMC results from disk, extracts the parameters corresponding to a 
    specified number of knots (if provided), and reconstructs the astrochronological time series using 
    the MAP estimate and bootstrapped parameter samples. It computes time series for eccentricity (either 
    based on precession envelope or gi-gj method), integrates the inverse sedimentation rate to obtain an 
    inferred timescale, and saves all relevant results to a temporary pickle file.

    Parameters:
        - **configuration_file_path** : str  
            Path to the configuration YAML file used for the analysis.

        - **folder_with_the_results** : str or None  
            Folder containing the result files. If `None`, the default path is used.

        - **file_with_the_results** : str or None  
            Base name of the results file. If `None`, the default name is used.

        - **selected_num_genes** : int or None  
            Number of knots (genes) to extract results for. If `None`, uses the number of depth genes from the selected solution.

        - **positive_feedback** : bool  
            Whether to multiply the eccentricity feedback signal by `+1` (`True`) or `-1` (`False`).

        - **use_prec_env** : bool  
            If `True`, uses the precession envelope for eccentricity. If `False`, uses the `gi - gj` formulation.

        - **ignore_weights** : bool  
            If `True`, ignores MCMC weights when resampling bootstrap parameters.


    Raises:
        ValueError: If the specified number of knots is not found in the MCMC results.

    Saves:
        A temporary pickle file named 'temp_selection_genes' inside a 'tmp' folder in the current working directory.
        This file contains a dictionary with the following keys:
            - 'number_of_knots_to_explore', 'selected_opt_result', 'selected_sampler_MCMC',
              'selected_logprob_MCMC', 'ind_MAP', 'params_MAP', 'invSR_MAP', 'freqs_MAP', 'fs',
              'res_pred', 'y_pred', 't_pred', 'params_bs', 'ip', 'ie', 'it', 'ecc_type',
              'times_e3_bs', 'times_e3_MAP', 'times_e3_ref', 'times_inferred', 'weight_results'.

    Example:
        >>> data_obtain_mcmc_results_per_a_number_of_knots(
        ...     configuration_file_path="configs/astro_config.yaml",
        ...     folder_with_the_results="results/",
        ...     file_with_the_results="astro_results",
        ...     selected_num_genes=5,
        ...     positive_feedback=True,
        ...     use_prec_env=False,
        ...     ignore_weights=False
        ... )
        # The function will save MCMC-based reconstructions and metadata to a temporary file for later use.
    """
    ########## LIBRARY MODE SECTION ##############
    
    results_fit, results_mcmc, weight_results, _, data_set_values, _, genetic_algorithm_parameters, data_model_parameters, mcmc_parameters, _ = notebooks_code.data_obtain_analysis_file_and_result_files(True, configuration_file_path, folder_with_the_results, file_with_the_results)

    data, _ = data_manager.data_read_dataset_from_file(data_set_values)
    tool_data_obtain_mcmc_results_per_a_number_of_knots(
        genetic_algorithm_parameters=genetic_algorithm_parameters,
        data_model_parameters=data_model_parameters,
        results_fit=results_fit,
        results_mcmc=results_mcmc,
        mcmc_parameters=mcmc_parameters,
        selected_num_genes=selected_num_genes,
        positive_feedback=positive_feedback,
        use_prec_env=use_prec_env,
        weight_results=weight_results,
        ignore_weights=ignore_weights
    )

def tool_data_obtain_mcmc_results_per_a_number_of_knots(
    genetic_algorithm_parameters:dict,
    data_model_parameters:dict,
    results_fit:dict,
    results_mcmc:dict,
    mcmc_parameters:dict,
    data: list,
    selected_num_genes: Optional[int],
    positive_feedback: bool,
    use_prec_env: bool,
    weight_results,
    ignore_weights: bool,
) -> None:
    
    list_number_genes_mcmc=mcmc_parameters["list_of_genes_mcmc"]
    original_list_number_genes= genetic_algorithm_parameters["list_number_genes"]
    frequency_values=data_model_parameters["frequency_values"]
    single_opt_result=results_fit["optimization_results"][-1]
    sampler_MCMC=results_mcmc["sampler_MCMC"]
    logprob_MCMC=results_mcmc["logprob_MCMC"]
    
    seed=genetic_algorithm_parameters["seed"]
    if seed:
        if isinstance(seed, int):
            setup_seed(seed)
        elif isinstance(seed, list):
            setup_seed(seed[0])
            
    dict_to_save = {}
    feedback_factor = 1 if positive_feedback else -1
    ecc_type = "(prec env)" if use_prec_env else "(gi-gj)"
    if selected_num_genes:
        try:
            selected_index = list_number_genes_mcmc.index(selected_num_genes)
            selected_index_opt_results = original_list_number_genes.index(
                selected_num_genes
            )
        except Exception as e:
            raise ValueError(
                f"The selected number of genes {selected_num_genes} does not exist."
            )
    else:
        selected_index = -1
        mcmc_genes = list_number_genes_mcmc[-1]
        selected_index_opt_results = original_list_number_genes.index(
                mcmc_genes
            )
    
    # WE SELECT THE LAST SOLUTION (-1) AND AFTER THE SELECTED NUMBER OF GENES (selected_index)
    selected_opt_result = single_opt_result[selected_index_opt_results]
    selected_sampler_MCMC = sampler_MCMC[selected_index]
    selected_logprob_MCMC = logprob_MCMC[selected_index]
    ind_MAP = np.unravel_index(
        np.argmax(selected_logprob_MCMC, axis=None),
        selected_logprob_MCMC.shape,
    )
    params_MAP = selected_sampler_MCMC[ind_MAP]
    atributes = get_result_atributes(selected_opt_result)
    invSR_MAP = params_MAP[atributes["number_free_freqs"] :]
    freqs_MAP = params_MAP[: atributes["number_free_freqs"]]
    fs = atributes["frequency_model"](freqs_MAP)
    res_pred = shared_functions.invSR_to_prediction(
        [atributes["depth_genes"], invSR_MAP],
        data,
        atributes["inverse_SR_lims"],
        fs,
        atributes["interpolator"],
    )
    y_pred = res_pred["y_pred"]
    t_pred = res_pred["time"]
    
    if not selected_num_genes:
        selected_num_genes = len(atributes["depth_genes"])

    n_bs = 1000  # what is this n_bs and why is it 2000?
    if not ignore_weights and isinstance(weight_results,np.ndarray):
        dict_to_save["weight_results"] = weight_results[selected_index]
        
        inds_1 = np.random.choice(
            range(len(selected_sampler_MCMC)),
            size=n_bs,
            replace=True,
            p=weight_results[selected_index],
        )  # AS, 2025-Feb-12
    else:
        dict_to_save["weight_results"] = None
        inds_1 = np.random.choice(
            range(len(selected_sampler_MCMC)), size=n_bs, replace=True
        )
    inds_2 = np.random.choice(
        range(len(selected_sampler_MCMC[0])), size=n_bs, replace=True
    )
    inds_1 = np.hstack([ind_MAP[0], inds_1])
    inds_2 = np.hstack([ind_MAP[1], inds_2])

    params_bs = np.stack([selected_sampler_MCMC[j, i] for j, i in zip(inds_1, inds_2)])

    ######################################################################################

    index_frequencies = 0
    if frequency_values["use_precession"]:
        ip = np.arange(5)
        index_frequencies = index_frequencies + 5
    else:
        ip = None
    if frequency_values["use_eccentricity"]:
        ie = np.arange(index_frequencies, index_frequencies + 5)
        index_frequencies = index_frequencies + 5
    else:
        ie = None
    if frequency_values["use_tilt"]:
        it = np.arange(index_frequencies, index_frequencies + 3)
    else:
        it = None

    ############################################

    times_e3_bs = []
    for params in params_bs:
        res_ETP = get_ETP_components(
            params,
            atributes["depth_genes"],
            atributes["frequency_model"],
            atributes["inverse_SR_lims"],
            atributes["interpolator"],
            data,
            ip,
            ie,
            it,
        )
        time = res_ETP["time"]
        if not use_prec_env:
            # then use gi-gj to build ecc
            y_pred_ecc = res_ETP["y_ecc"]
            y_pred_ecc = y_pred_ecc * feedback_factor
            ecc_est = y_pred_ecc
        else:
            # use precession envelope to build by default
            y_env_prec = res_ETP["y_env_prec"]
            ecc_est = y_env_prec
        ecc_est = (ecc_est - ecc_est.mean()) / ecc_est.std()
        times_e3_bs.append(np.array([time, ecc_est]))
    times_e3_bs = np.array(times_e3_bs)
    times_e3_MAP = times_e3_bs[0]
    times_e3_ref = times_e3_MAP

    ######################################
    invSR_MAP = params_MAP[atributes["number_free_freqs"] :]
    invSR_interpolate = atributes["interpolator"](
        [atributes["depth_genes"], invSR_MAP], data[0]
    )
    # CHANGE: THIS VARIABLE IS USED IN OTHER PLACES
    times_inferred = sp.integrate.cumulative_trapezoid(
        invSR_interpolate, data[0], initial=0
    )

    ### ADD SAVING OF THIS DATA.
    dict_to_save["number_of_knots_to_explore"] = selected_num_genes
    dict_to_save["selected_opt_result"] = selected_opt_result
    dict_to_save["selected_sampler_MCMC"] = selected_sampler_MCMC
    dict_to_save["selected_logprob_MCMC"] = selected_logprob_MCMC
    dict_to_save["ind_MAP"] = ind_MAP
    dict_to_save["params_MAP"] = params_MAP
    dict_to_save["invSR_MAP"] = invSR_MAP
    dict_to_save["freqs_MAP"] = freqs_MAP
    dict_to_save["fs"] = fs
    dict_to_save["res_pred"] = res_pred
    dict_to_save["y_pred"] = y_pred
    dict_to_save["t_pred"] = t_pred
    dict_to_save["params_bs"] = params_bs
    dict_to_save["ip"] = ip
    dict_to_save["ie"] = ie
    dict_to_save["it"] = it
    dict_to_save["ecc_type"] = ecc_type
    dict_to_save["times_e3_bs"] = times_e3_bs
    dict_to_save["times_e3_MAP"] = times_e3_MAP
    dict_to_save["times_e3_ref"] = times_e3_ref
    dict_to_save["times_inferred"] = times_inferred

    temp_path = f"{os.getcwd()}/tmp"
    temp_file = "temp_selection_genes"
    temp_file_path = f"{temp_path}/{temp_file}"
    if os.path.exists(temp_file_path):
        os.remove(temp_file_path)
    if not os.path.exists(temp_path):
        os.mkdir(temp_path)
    data_manager.save_data(
        output_file_path=temp_file_path, fit_mcmc_weights="mcmc", **dict_to_save
    )
