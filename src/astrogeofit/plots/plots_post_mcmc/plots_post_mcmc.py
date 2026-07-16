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
import random
import numpy as np  # type:ignore
import scipy as sp  # type:ignore
import logging
import matplotlib.pyplot as plt  # type:ignore
import matplotlib.transforms as mtransforms  # type:ignore

from typing import Optional
from functools import partial
from sklearn.metrics import r2_score  # type:ignore
from matplotlib.ticker import FuncFormatter  # type:ignore
from sklearn.linear_model import LinearRegression  # type:ignore
from scipy.interpolate import interp1d  # type:ignore

from astrogeofit.main_routines import data_manager
from astrogeofit.utils import shared_functions, setup_logger
from astrogeofit.plots.utils import shared_utils_plots as plots_utils
from astrogeofit.plots.plots_post_mcmc import utils_post_mcmc
from astrogeofit.mcmc_routines import utils as mcmc_utils
from astrogeofit.notebooks_code import notebooks_code

supported_formats = ["jpeg", "pdf", "eps", "tab"]

logger = logging.getLogger("ToolLogger")
# Check if the logger already has handlers
if not logger.hasHandlers():
    logger = setup_logger.setup_logger(False)


def setup_seed(seed: int):
    """Set the global seed for reproducibility."""
    random.seed(seed)
    np.random.seed(seed)

def plt_mcmc_SR_per_num_knots(
    configuration_file_path:str,
    folder_with_the_results:Optional[str],
    file_with_the_results:Optional[str],
    x_axis_plot_limits: Optional[list],
    y_axis_plot_limits: Optional[list],
    ignore_weights: bool,
    types_of_saving: list,
    output_path: str,
):
    """
    Plot the inverse sedimentation rate (SR⁻¹) derived from MCMC simulations 
    for different numbers of knots in the age-depth model.

    This function loads the results of an AstroGeoFit analysis, specifically the 
    MCMC sampling results for different age-depth model complexities (defined by 
    number of knots). It visualizes the median inverse sedimentation rate profiles, 
    confidence intervals, MAP estimates, and optionally the nominal model.

    Parameters:
        - **configuration_file_path** : str
            Path to the configuration file used for the AstroGeoFit analysis.

        - **folder_with_the_results** : str or None
            Path to the folder containing the result files. If None, the default folder 
            defined in the configuration file is used.

        - **file_with_the_results** : str or None
            Specific filename for the results. If None, a default filename is assumed.

        - **x_axis_plot_limits** : list or None
            Optional [xmin, xmax] for the depth axis. If None, limits are automatically inferred.

        - **y_axis_plot_limits** : list or None
            Optional [ymin, ymax] for the inverse SR axis. If None, limits are automatically inferred.

        - **ignore_weights** : bool
            If True, ignore the sampling weights when drawing MCMC samples.

        - **types_of_saving** : list of str
            List of file formats to save the resulting plot (e.g., ["png", "pdf"]). 
            If empty, the figure will not be saved.

        - **output_path** : str
            Output path where the figure will be saved. If "default" or empty, the 
            `figures` subfolder inside the output folder defined in the configuration 
            file is used.

    Notes:
        - Each subplot corresponds to a model with a different number of free knots 
        (i.e., free parameters in the age-depth model).
        - The plot shows the MAP solution, the median curve of the inverse SR samples, 
        and a shaded confidence interval (default 90% CI).
        - If `add_nominal` is True (from age-depth model settings), the nominal model 
        will also be plotted in red.
        - The output figure layout adjusts automatically depending on the number of 
        models being compared.

    Outputs:
        - A matplotlib figure with subplots, each showing inverse SR results per 
        number of knots used in the MCMC model.
        - Saved figure(s) in the specified format(s) and output location if requested.
    """
    ########## LIBRARY MODE SECTION ##############
    
    results_fit, results_mcmc, weight_results, parameters_analysis, data_set_values, age_depth_model_values, genetic_algorithm_parameters, _, mcmc_parameters, _ = notebooks_code.data_obtain_analysis_file_and_result_files(True, configuration_file_path, folder_with_the_results, file_with_the_results)

    _, _, _, func_invSR_nominal, _, add_nominal = notebooks_code.data_obtain_age_depth_values(age_depth_model_values)

    data, _ = data_manager.data_read_dataset_from_file(data_set_values)
    tool_plt_mcmc_SR_per_num_knots(
        genetic_algorithm_parameters=genetic_algorithm_parameters,
        results_fit=results_fit,
        results_mcmc=results_mcmc,
        mcmc_parameters=mcmc_parameters,
        data=data,
        func_invSR_nominal=func_invSR_nominal,
        parameters_analysis=parameters_analysis,
        x_lims=x_axis_plot_limits,
        y_lims=y_axis_plot_limits,
        weight_results=weight_results,
        ignore_weights=ignore_weights,
        types_of_saving=types_of_saving,
        output_path=output_path,
        add_nominal=add_nominal
    )

def tool_plt_mcmc_SR_per_num_knots(
    genetic_algorithm_parameters:dict,
    results_fit:dict,
    results_mcmc:dict,
    mcmc_parameters:dict,
    data: list,
    func_invSR_nominal,
    parameters_analysis: dict,
    x_lims: Optional[list],
    y_lims: Optional[list],
    weight_results: Optional[dict],
    ignore_weights: bool,
    types_of_saving: list,
    output_path: str,
    add_nominal: bool
):
    seed = genetic_algorithm_parameters["seed"]
    list_number_genes_mcmc=mcmc_parameters["list_of_genes_mcmc"]
    original_list_genes=genetic_algorithm_parameters["list_number_genes"]
    single_opt_result=results_fit["optimization_results"][-1]
    sampler_MCMC=results_mcmc["sampler_MCMC"]
    logprob_MCMC=results_mcmc["logprob_MCMC"]
    
    if seed:
        if isinstance(seed, int):
            setup_seed(seed)
        elif isinstance(seed, list):
            setup_seed(seed[0])
    
    n_plots = len(list_number_genes_mcmc)
    n_cols = int(np.ceil(np.sqrt(n_plots)))
    n_rows = int(np.ceil(n_plots / n_cols))

    fig, axs = plt.subplots(
        n_rows,
        n_cols,
        figsize=(12, 4),
        constrained_layout=True,
        sharex=True,
        sharey=True,
    )

    # Flatten axs if it's a 2D array
    if n_rows > 1 and n_cols > 1:
        axs = axs.flatten()
    elif n_rows == 1 or n_cols == 1:
        axs = np.array(axs).flatten()

    alpha = [0.05, 0.95]
    n = 1000

    for j, num_genes in enumerate(list_number_genes_mcmc):
        index_opt_result = original_list_genes.index(num_genes)
        result = single_opt_result[index_opt_result]
        atributes = plots_utils.get_result_atributes(result)
        depth_genes = atributes["depth_genes"]

        samples_combine = sampler_MCMC[j][:, :, atributes["number_free_freqs"] :]
        logprob_samples = logprob_MCMC[j][:]
        ind_MAP = np.unravel_index(
            np.argmax(logprob_samples, axis=None), logprob_samples.shape
        )
        if not ignore_weights and isinstance(weight_results, np.ndarray):
            inds_1 = np.random.choice(
                range(len(samples_combine)), size=n, replace=True, p=weight_results[j]
            )
        else:
            inds_1 = np.random.choice(range(len(samples_combine)), size=n, replace=True)

        inds_2 = np.random.choice(range(len(samples_combine[0])), size=n, replace=True)
        invSR_interpolate_par = np.array(
            [
                atributes["interpolator"]([depth_genes, samples_combine[j, i]], data[0])
                for j, i in zip(inds_1, inds_2)
            ]
        )
        mean_invSR = np.median(invSR_interpolate_par, axis=0)
        invSR_lo, invSR_up = np.quantile(invSR_interpolate_par, alpha, axis=0)
        axs[j].plot(data[0], mean_invSR, "k--", label="median", zorder=10)
        axs[j].fill_between(
            data[0],
            invSR_lo,
            invSR_up,
            alpha=0.3,
            color="black",
            label=f"{(alpha[1] - alpha[0]) * 100:.0f}% CI",
        )
        invSR_MAP = atributes["interpolator"](
            [depth_genes, samples_combine[ind_MAP]], data[0]
        )
        axs[j].plot(data[0], invSR_MAP, color="k", label="MAP")
        if add_nominal:
            axs[j].plot(
                data[0],
                func_invSR_nominal(data[0]),
                color="red",
                zorder=-1000,
                label="Nominal",
            )
        axs[j].set_xlim(data[0][0], data[0][-1])
        if y_lims:
            axs[j].set_ylim(y_lims[0], y_lims[1])
        if x_lims:
            axs[j].set_xlim(x_lims[0], x_lims[1])
        axs[j].set_title(f"{len(depth_genes)} knots")

    axs[0].legend(loc="upper left")

    # Set the y-label for the first column
    for i in range(n_rows):
        axs[i * n_cols].set_ylabel("SR^-1 [Myr/m]")

    # Set the x-label for the last row
    for j in range(n_cols):
        axs[(n_rows - 1) * n_cols + j].set_xlabel("Depth [m]")
    if types_of_saving and types_of_saving != []:
        if not output_path or output_path == "default" or output_path == "":
            output_path = parameters_analysis["output_folder"]
            output_path = f"{output_path}/figures"
            os.makedirs(output_path, exist_ok=True)
        figname = "Sedimentation_Rate_Per_Number_Of_Knots_MCMC"
        global supported_formats
        for fmt in types_of_saving:
            if fmt.lower() in supported_formats:
                if fmt.lower() == "tab":
                    logger.warning("tab format is not accepted for this plot.")
                else:
                    fig.savefig(
                        f"{output_path}/{figname}.{fmt.lower()}", format=fmt.lower()
                    )
                    logger.info(f"Saved {output_path}/{figname}.{fmt.lower()}")
            else:
                logger.warning(f"The format {fmt} is not accepted.")
    ##plt.show()

def plt_mcmc_aic_logprob_and_loglike_per_number_knots(
    configuration_file_path:str,
    folder_with_the_results:Optional[str],
    file_with_the_results:Optional[str],
    types_of_saving: list,
    output_path: str,
):
    """
    Plot AIC, Logprob and Loglike metrics for MCMC-inferred models with varying numbers of knots.

    This function evaluates the model fit for each number of knots used in an 
    age-depth model by calculating the Akaike Information Criterion (AIC) and the 
    coefficient of determination (R²) using MCMC posterior samples. The metrics 
    are plotted together using a dual y-axis plot.

    Parameters:
        - **configuration_file_path** : str
            Path to the configuration file used for the AstroGeoFit analysis.

        - **folder_with_the_results** : str or None
            Path to the folder containing the result files. If None, the default 
            folder specified in the configuration is used.

        - **file_with_the_results** : str or None
            Specific filename for the results. If None, the default name is assumed.

        - **types_of_saving** : list of str
            List of file formats in which to save the resulting plot (e.g., ["png", "pdf"]).

        - **output_path** : str
            Directory where the figure should be saved. If "default" or empty, it defaults 
            to a `figures` folder inside the output folder defined in the configuration.

    Notes:
        - For each model (defined by a different number of knots), the function calculates 
        the log-likelihood of a subset of MCMC samples.
        - AIC and BIC values are derived using the maximum log-likelihood value.
        - R² values are calculated from inverse SR-based predictions compared to the observed data.
        - A dual y-axis plot is generated with AIC on the left and R² on the right.
        - AIC values are normalized by subtracting the minimum AIC value (i.e., first entry) for visual clarity.

    Outputs:
        - A matplotlib figure showing AIC (black line) and R² (red line) as a function of 
        the number of knots in the model.
        - If `types_of_saving` is non-empty, the figure is saved to the specified `output_path`.

    Warnings:
        - The 'tab' file format is not supported and will be skipped if included.
    """
    results_fit, results_mcmc, _, parameters_analysis, data_set_values, _, genetic_algorithm_parameters, _, mcmc_parameters, _ = notebooks_code.data_obtain_analysis_file_and_result_files(True, configuration_file_path, folder_with_the_results, file_with_the_results)

    data, _ = data_manager.data_read_dataset_from_file(data_set_values)
    tool_plt_mcmc_aic_logprob_and_loglike_per_number_knots(
        genetic_algorithm_parameters=genetic_algorithm_parameters,
        results_fit=results_fit,
        results_mcmc=results_mcmc,
        mcmc_parameters=mcmc_parameters,
        data=data,
        parameters_analysis=parameters_analysis,
        types_of_saving=types_of_saving,
        output_path=output_path
    )

def tool_plt_mcmc_aic_logprob_and_loglike_per_number_knots(
    genetic_algorithm_parameters:dict,
    results_fit:dict,
    results_mcmc:dict,
    mcmc_parameters:dict,
    data: list,
    parameters_analysis:dict,
    types_of_saving:list,
    output_path:str
):
    
    def color_y_axis(ax, color):
        """Color your axes."""
        for t in ax.get_yticklabels():
            t.set_color(color)
        return None
    
    single_opt_result=results_fit["optimization_results"][-1]
    original_list_genes=genetic_algorithm_parameters["list_number_genes"]
    list_number_genes_mcmc=mcmc_parameters["list_of_genes_mcmc"]
    sampler_MCMC=results_mcmc["sampler_MCMC"]
    logprob_MCMC=results_mcmc["logprob_MCMC"]
    seed=genetic_algorithm_parameters["seed"]
    if seed:
        if isinstance(seed, int):
            setup_seed(seed)
        elif isinstance(seed, list):
            setup_seed(seed[0])
            
    logllh_seq = []
    max_llh = np.zeros(len(list_number_genes_mcmc)) # This is the length of the list of the number of genes as we have one result per each number of genes
    aics = np.zeros(len(list_number_genes_mcmc))
    bics = np.zeros(len(list_number_genes_mcmc))
    list_n_genes = np.zeros(len(list_number_genes_mcmc))
    r2s = np.zeros(len(list_number_genes_mcmc))
    logProb = np.zeros(len(list_number_genes_mcmc)) # YW
    loglike = np.zeros(len(list_number_genes_mcmc)) # YW
    for k, num_gen in enumerate(list_number_genes_mcmc):
        index_opt_result = original_list_genes.index(num_gen)
        result = single_opt_result[index_opt_result]
        atributes = plots_utils.get_result_atributes(result)
        lags = 2
        logllh = partial(
            mcmc_utils.get_log_likelihood,
            depth_of_invSR=atributes["depth_genes"],
            freqs_model=atributes["frequency_model"],
            inverse_SR_lims=atributes["inverse_SR_lims"],
            data=data,
            interpolator=atributes["interpolator"],
            lags=lags,
        )
        logllh_seq.append(logllh)
        selected_sample_MCMC = sampler_MCMC[k]
        # n = 100
        n = 1000   # ----YW, 2025-Apr-10
        inds_1 = np.random.choice(
            range(len(selected_sample_MCMC)), size=n, replace=True
        )
        inds_2 = np.random.choice(
            range(len(selected_sample_MCMC[0])), size=n, replace=True
        )
        llh_sp = [
            logllh_seq[k](selected_sample_MCMC[j1, j2])
            for j1, j2 in zip(inds_1, inds_2)
        ]
        max_llh[k] = max(llh_sp)
        # bics[k] = len(atributes["depth_genes"]) * np.log(len(data[1])) - 2 * max_llh[k]
        aics[k] = len(atributes["depth_genes"]) * 2 - 2 * max_llh[k]
        list_n_genes[k] = len(atributes["depth_genes"])

        # ----YW ----
        depth_genes = atributes["depth_genes"]
        samples_combine = sampler_MCMC[k][:, :, atributes["number_free_freqs"]:]
        logprob_samples = logprob_MCMC[k][:]
        ind_MAP = np.unravel_index(
            np.argmax(logprob_samples, axis=None), logprob_samples.shape
        )
        logProb[k]=logprob_samples[ind_MAP]
        # ----YW ----

        loglike[k]=logllh_seq[k](selected_sample_MCMC[ind_MAP])
        
        
    fig, ax = plt.subplots(1, 1)
    plt.subplots_adjust(left=0.15,right=0.85,bottom=0.15,
                    top=0.9)
    ax2 = ax.twinx()
    ax2.plot(list_n_genes, logProb, ".-", color="r",label='posterior') # YW
    ax2.plot(list_n_genes, loglike, ".-", color="b",label='likelihood') # YW
    ax2.legend(loc='best')
    ax2.set_ylabel('log probablity', color="r", rotation=-90, verticalalignment='bottom') # YW
    ax.plot(list_n_genes, aics - aics[0], ".-", color="k")
    ax.set_xlabel("Number of Knots")
    ax.set_ylabel("AIC")
    color_y_axis(ax, "k")
    color_y_axis(ax2, "r")
    
    if types_of_saving and types_of_saving != []:
        if not output_path or output_path == "default" or output_path == "":
            output_path = parameters_analysis["output_folder"]
            output_path = f"{output_path}/figures"
            os.makedirs(output_path, exist_ok=True)
        figname = "AIC_and_logprob"
        global supported_formats
        for fmt in types_of_saving:
            if fmt.lower() in supported_formats:
                if fmt.lower() == "tab":
                    logger.warning("tab format is not accepted for this plot.")
                else:
                    fig.savefig(f"{output_path}/{figname}.{fmt.lower()}", format=fmt.lower())
                    logger.info(f"Saved {output_path}/{figname}.{fmt.lower()}")
            else:
                logger.warning(f"The format {fmt} is not accepted.")
    #plt.show()

def plt_mcmc_aic_and_r2_per_number_knots(
    configuration_file_path:str,
    folder_with_the_results:Optional[str],
    file_with_the_results:Optional[str],
    types_of_saving: list,
    output_path: str,
):
    """
    Plot AIC and R² metrics for MCMC-inferred models with varying numbers of knots.

    This function evaluates the model fit for each number of knots used in an 
    age-depth model by calculating the Akaike Information Criterion (AIC) and the 
    coefficient of determination (R²) using MCMC posterior samples. The metrics 
    are plotted together using a dual y-axis plot.

    Parameters:
        - **configuration_file_path** : str
            Path to the configuration file used for the AstroGeoFit analysis.

        - **folder_with_the_results** : str or None
            Path to the folder containing the result files. If None, the default 
            folder specified in the configuration is used.

        - **file_with_the_results** : str or None
            Specific filename for the results. If None, the default name is assumed.

        - **types_of_saving** : list of str
            List of file formats in which to save the resulting plot (e.g., ["png", "pdf"]).

        - **output_path** : str
            Directory where the figure should be saved. If "default" or empty, it defaults 
            to a `figures` folder inside the output folder defined in the configuration.

    Notes:
        - For each model (defined by a different number of knots), the function calculates 
        the log-likelihood of a subset of MCMC samples.
        - AIC and BIC values are derived using the maximum log-likelihood value.
        - R² values are calculated from inverse SR-based predictions compared to the observed data.
        - A dual y-axis plot is generated with AIC on the left and R² on the right.
        - AIC values are normalized by subtracting the minimum AIC value (i.e., first entry) for visual clarity.

    Outputs:
        - A matplotlib figure showing AIC (black line) and R² (red line) as a function of 
        the number of knots in the model.
        - If `types_of_saving` is non-empty, the figure is saved to the specified `output_path`.

    Warnings:
        - The 'tab' file format is not supported and will be skipped if included.
    """
    results_fit, results_mcmc, _, parameters_analysis, data_set_values, _, genetic_algorithm_parameters, _, mcmc_parameters, _ = notebooks_code.data_obtain_analysis_file_and_result_files(True, configuration_file_path, folder_with_the_results, file_with_the_results)

    data, _ = data_manager.data_read_dataset_from_file(data_set_values)
    
    tool_plt_mcmc_aic_and_r2_per_number_knots(
        genetic_algorithm_parameters=genetic_algorithm_parameters,
        results_fit=results_fit,
        results_mcmc=results_mcmc,
        mcmc_parameters=mcmc_parameters,
        data=data,
        parameters_analysis=parameters_analysis,
        types_of_saving=types_of_saving,
        output_path=output_path
    )
    
def tool_plt_mcmc_aic_and_r2_per_number_knots(
    genetic_algorithm_parameters:dict,
    results_fit:dict,
    results_mcmc:dict,
    mcmc_parameters:dict,
    data:list,
    parameters_analysis: dict,
    types_of_saving: list,
    output_path: str,
):
    def color_y_axis(ax, color):
        """Color your axes."""
        for t in ax.get_yticklabels():
            t.set_color(color)
        return None
    
    single_opt_result=results_fit["optimization_results"][-1]
    original_list_genes=genetic_algorithm_parameters["list_number_genes"]
    list_number_genes_mcmc=mcmc_parameters["list_of_genes_mcmc"]
    sampler_MCMC=results_mcmc["sampler_MCMC"]
    
    seed=genetic_algorithm_parameters["seed"]
    if seed:
        if isinstance(seed, int):
            setup_seed(seed)
        elif isinstance(seed, list):
            setup_seed(seed[0])
    
    logllh_seq = []
    max_llh = np.zeros(
        len(list_number_genes_mcmc)
    )  # This is the length of the list of the number of genes as we have one result per each number of genes
    aics = np.zeros(len(list_number_genes_mcmc))
    bics = np.zeros(len(list_number_genes_mcmc))
    list_n_genes = np.zeros(len(list_number_genes_mcmc))
    r2s = np.zeros(len(list_number_genes_mcmc))
    for k, num_gen in enumerate(list_number_genes_mcmc):
        index_opt_result = original_list_genes.index(num_gen)
        result = single_opt_result[index_opt_result]
        atributes = plots_utils.get_result_atributes(result)
        lags = 2
        logllh = partial(
            mcmc_utils.get_log_likelihood,
            depth_of_invSR=atributes["depth_genes"],
            freqs_model=atributes["frequency_model"],
            inverse_SR_lims=atributes["inverse_SR_lims"],
            data=data,
            interpolator=atributes["interpolator"],
            lags=lags,
        )
        logllh_seq.append(logllh)
        selected_sample_MCMC = sampler_MCMC[k]
        n = 1000 # We should find a way to set this number so the plot does not change too much between executions
        inds_1 = np.random.choice(
            range(len(selected_sample_MCMC)), size=n, replace=True
        )
        inds_2 = np.random.choice(
            range(len(selected_sample_MCMC[0])), size=n, replace=True
        )
        llh_sp = [
            logllh_seq[k](selected_sample_MCMC[j1, j2])
            for j1, j2 in zip(inds_1, inds_2)
        ]
        max_llh[k] = max(llh_sp)
        bics[k] = len(atributes["depth_genes"]) * np.log(len(data[1])) - 2 * max_llh[k]
        aics[k] = len(atributes["depth_genes"]) * 2 - 2 * max_llh[k]
        list_n_genes[k] = len(atributes["depth_genes"])
        res_pred = shared_functions.invSR_to_prediction(
            [atributes["depth_genes"], atributes["best_inverse_SR"]],
            data,
            atributes["inverse_SR_lims"],
            atributes["frequencies"],
            atributes["interpolator"],
        )
        y_pred = res_pred["y_pred"]
        r2s[k] = r2_score(data[1], y_pred)

    fig, ax = plt.subplots(1, 1)
    plt.subplots_adjust(left=0.15,right=0.85,bottom=0.15,
                    top=0.9)
    ax2 = ax.twinx()
    ax2.plot(list_n_genes, r2s, ".-", color="r")
    ax2.set_ylabel("r$^2$", color="r", rotation=0)
    ax.plot(list_n_genes, aics - aics[0], ".-", color="k")
    ax.set_xlabel("Number of Knots")
    ax.set_ylabel("AIC")
    color_y_axis(ax, "k")
    color_y_axis(ax2, "r")
    if types_of_saving and types_of_saving != []:
        if not output_path or output_path == "default" or output_path == "":
            output_path = parameters_analysis["output_folder"]
            output_path = f"{output_path}/figures"
            os.makedirs(output_path, exist_ok=True)
        figname = "AIC_and_Metric_Scores_Per_Number_of_Knots"
        global supported_formats
        for fmt in types_of_saving:
            if fmt.lower() in supported_formats:
                if fmt.lower() == "tab":
                    logger.warning("tab format is not accepted for this plot.")
                else:
                    fig.savefig(
                        f"{output_path}/{figname}.{fmt.lower()}", format=fmt.lower()
                    )
                    logger.info(f"Saved {output_path}/{figname}.{fmt.lower()}")
            else:
                logger.warning(f"The format {fmt} is not accepted.")
    ##plt.show()

def plt_mcmc_eccentricity_correlation_with_solution(
    configuration_file_path:str,
    folder_with_the_results:Optional[str],
    file_with_the_results:Optional[str],
    t0_offset_range_plot: Optional[list],
    ecc_time_range: Optional[list],
    save_eccentricity_solutions: bool,
    save_eccentricity_solutions_comment: Optional[str],
    name_eccentricity_solution: Optional[str],
    types_of_saving: list,
    output_path: str,
):
    """
    Plot the results of eccentricity from MCMC runs and compare them with an astronomical solution.

    This function generates a multi-panel figure showing:
    - The mean square weighted deviation (MSWD) across different t0 offsets.
    - The aligned time offset between relative and absolute time.
    - The fitted eccentricity envelope overlaid with a reference astronomical solution.

    It also logs relevant metrics and saves the resulting figure in the specified formats.

    Parameters:
        - **configuration_file_path** : str
            Path to the configuration file used in the analysis.

        - **folder_with_the_results** : str or None
            Path to the folder containing the MCMC result files. Can be `None` if `file_with_the_results` is provided.

        - **file_with_the_results** : str or None
            Specific path to a result file. Can be `None` if `folder_with_the_results` is provided.

        - **t0_offset_range_plot** : list or None
            Range of t0 offsets to evaluate MSWD. Should be a list of two floats [min_offset, max_offset].
            If `None`, it is inferred from the nominal age model if available.

        - **ecc_time_range** : list or None
            Time range for plotting the eccentricity envelope. If `None`, it's automatically inferred from the dataset.

        - **save_eccentricity_solutions**: bool
            Flag indicating whether to save the reconstructed eccentricity solutions and their parameters to output files.
            If set to True, the function will generate:

                A .txt file containing the depth, age, proxy values, and statistical summaries (MAP, median, 5th, and 95th percentiles) of the eccentricity.

                A .sol file containing the frequency components, their amplitudes (beta coefficients), and phase information derived from the MAP (maximum a posteriori) solution.

            These files allow post-analysis use of the eccentricity time series, e.g., for replication, comparison, or modeling.

            If False, these outputs are skipped.
            
        - **save_eccentricity_solutions_comment**: str or None
            A user-defined comment or note to be included in the header of the saved eccentricity solution output files (e.g., .txt and .sol). 
            This allows users to annotate saved results with context, version information, data origin, or other metadata

        - **name_eccentricity_solution** : str or None
            Label for the astronomical solution being compared (e.g., "La2004", "Astronomical Solution"). Default is "Astronomical Solution".

        - **types_of_saving** : list
            List of formats in which to save the output figure (e.g., ['pdf', 'png']).

        - **output_path** : str
            Directory where the output figures and logs should be saved. If 'default', uses the default output folder
            specified in the analysis configuration.

    Returns:
        None
            This function saves the plots and logs to the file system but does not return any object.
    
    Notes:
        - The function expects a temporary pickle file (`temp_selection_genes_mcmc.pickle`) to be available in `tmp/`.
        - If an age-depth model is not available, `t0_offset_range_plot` must be provided.
        - Logs are saved in a text file named `logs_plots.txt` inside the output path.
        - The figure contains three subplots: MSWD vs t0, relative/absolute time offsets, and eccentricity comparison.
    
    Warnings:
        - The function will log a warning if unsupported save formats (like 'tab') are passed in `types_of_saving`.
    """
    _, _, _, parameters_analysis, data_set_values, age_depth_model_values, genetic_algorithm_parameters, _, _, _ = notebooks_code.data_obtain_analysis_file_and_result_files(True, configuration_file_path, folder_with_the_results, file_with_the_results)

    _, _, _, _, func_time_nominal, add_nominal = notebooks_code.data_obtain_age_depth_values(age_depth_model_values)

    data, _ = data_manager.data_read_dataset_from_file(data_set_values)
    
    tool_plt_mcmc_eccentricity_correlation_with_solution(
        genetic_algorithm_parameters=genetic_algorithm_parameters,
        func_time_nominal=func_time_nominal,
        data=data,
        parameters_analysis=parameters_analysis,
        save_eccentricity_solutions=save_eccentricity_solutions,
        save_eccentricity_solutions_comment=save_eccentricity_solutions_comment,
        t0_offset_range_plot=t0_offset_range_plot,
        ecc_time_range=ecc_time_range,
        name_eccentricity_solution=name_eccentricity_solution,
        types_of_saving=types_of_saving,
        output_path=output_path,
        add_nominal=add_nominal
    )
    
def tool_plt_mcmc_eccentricity_correlation_with_solution(
    genetic_algorithm_parameters:dict,
    func_time_nominal,
    data: list,
    parameters_analysis: dict,
    save_eccentricity_solutions: bool,
    save_eccentricity_solutions_comment: str,
    t0_offset_range_plot: Optional[list],
    ecc_time_range: Optional[list],
    name_eccentricity_solution: Optional[str],
    types_of_saving: list,
    output_path: str,
    add_nominal: bool
):
    def get_phase_beta(y, fs, times_inferred):
        """
        Input: y: the observed data
            fs: the orbital frequencies, such as eccentricity, obliquity, and climatic precession
            times_infferred: the time scale
        Output: phi, an array which contains the phases of orbital frequencies
        """
        X = shared_functions.generate_X_linReg(
            np.ones_like(fs), fs, times_inferred
        )  # Generate X
        reg_model = LinearRegression().fit(X, y)  # Get the regression funciton.
        beta_sin, beta_cos = (
            reg_model.coef_[: len(fs)],
            reg_model.coef_[len(fs) :],
        )  # Get coefficients
        phi = np.arctan2(-beta_sin, beta_cos)  # Calculate the phase
        return phi, beta_sin, beta_cos

    eccentricity_params=parameters_analysis["eccentricity_solution_data"]
    dict_results_gene_selection = data_manager.load_dictionary_from_pickle(
        f"{os.getcwd()}/tmp/temp_selection_genes_mcmc.pickle"
    )
    params_bs = dict_results_gene_selection["params_bs"]
    times_e3_bs = dict_results_gene_selection["times_e3_bs"]
    times_e3_MAP = dict_results_gene_selection["times_e3_MAP"]
    times_e3_ref = dict_results_gene_selection["times_e3_MAP"]
    ecc_type = dict_results_gene_selection["ecc_type"]  # ----YW， 2025-Mar-11
    number_of_knots_explored = dict_results_gene_selection["number_of_knots_to_explore"]

    selected_opt_result = dict_results_gene_selection["selected_opt_result"]
    selected_sampler_MCMC = dict_results_gene_selection["selected_sampler_MCMC"]
    selected_logprob_MCMC = dict_results_gene_selection["selected_logprob_MCMC"]
    freqs_MAP = dict_results_gene_selection["freqs_MAP"]
    fs = dict_results_gene_selection["fs"]
    atributes = plots_utils.get_result_atributes(selected_opt_result)
    depth_genes = selected_opt_result.problem.depth_genes
    freqs_model = selected_opt_result.problem.freqs_model
    interpolator = selected_opt_result.problem.interpolator
    number_free_freqs = atributes["number_free_freqs"]
    
    seed=genetic_algorithm_parameters["seed"]
    if seed:
        if isinstance(seed, int):
            setup_seed(seed)
        elif isinstance(seed, list):
            setup_seed(seed[0])

    if not name_eccentricity_solution or name_eccentricity_solution == "":
        name_eccentricity_solution = "Astronomical Solution"
    fig, axs = plt.subplot_mosaic(
        [["a)", "b)"], ["c)", "c)"]], constrained_layout=True, figsize=(12, 4)
    )

    for label, ax in axs.items():
        # label physical distance in and down:
        trans = mtransforms.ScaledTranslation(10 / 72, -5 / 72, fig.dpi_scale_trans)
        ax.text(
            0,
            1,
            label,
            transform=ax.transAxes + trans,
            fontsize="medium",
            verticalalignment="top",
            fontfamily="serif",
            bbox=dict(facecolor="0.7", edgecolor="none", pad=3.0),
        )
    ax_MSWD = axs["a)"]
    ax_time_offset = axs["b)"]
    ax_ecc_a = axs["c)"]

    """
    get local offset with respect to times_e3_ref
    """
    offset_range = [-0.1, 0.1]  # QUESTION: eccentricity plot, what is it?
    offset_bs = np.zeros(len(params_bs))
    correlation_bs = np.zeros(len(params_bs))
    for i, t_e3 in enumerate(times_e3_bs):
        offset_bs[i], correlation_bs[i] = utils_post_mcmc.get_offset(
            t_e3, times_e3_MAP, offset_range
        )

    """
    get aligned times_e3_bs, that is e3_bs interpolated (with timeoffset) in the times_ref 
    """
    aligned_e3_bs = []
    for (t_bs, e3_bs), best_offset in zip(times_e3_bs, offset_bs):
        t_corrected = t_bs + best_offset
        interp_y_corrected = interp1d(
            t_corrected, e3_bs, kind="linear", bounds_error=False, fill_value=np.nan
        )
        y_corrected_interp = interp_y_corrected(times_e3_ref[0])
        aligned_e3_bs.append(y_corrected_interp)

    aligned_e3_bs = np.array(aligned_e3_bs)

    e3_interp = shared_functions.get_eccentricity_parameters(
        eccentricity_params=eccentricity_params
    )

    """
    get absolute time according to the weighted sum of squared deviation (chi squared)
    """
    # Set the xlims/t0_range for the subfigure of "MSWD over t0"
    if t0_offset_range_plot == None or t0_offset_range_plot == []:
        if add_nominal == False:
            print(
                "IT IS NEEDED A RANGE FOR THE T0 OFFSET, AS THERE IS NO AGE DEPTH MODEL TO EXTRACT THE DATA. PLEASE INTRODUCE IT."
            )
        else:
            t0_offset_range_plot = [
                func_time_nominal(data[0][0]) - 1,
                func_time_nominal(data[0][0]) + 1,
            ]

    # Set the xlims/time_range for the subfigure of "eccentricity over time"
    if ecc_time_range == None or ecc_time_range == []:
        ecc_time_range = [
            func_time_nominal(data[0][0]) - 0.5,
            func_time_nominal(data[0][-1]) + 0.5,
        ]

    var_e3_est = np.nanvar(aligned_e3_bs, axis=0)
    # QUESTION: Here t0 is used as linespace of t0_offset_range, but under we just change t0 and it becomes the x limits of the plot. Should these values be the same? If no, what controls each one?
    t0s = np.linspace(t0_offset_range_plot[0], t0_offset_range_plot[1], 2000)
    chi_squared_t0s = np.sum(
        (e3_interp(times_e3_ref[0, None] + t0s[:, None]) - times_e3_ref[1]) ** 2
        / var_e3_est[None],
        axis=1,
    )

    ax = ax_MSWD

    chi_squared_t0s = np.mean(
        (e3_interp(times_e3_ref[0, None] + t0s[:, None]) - times_e3_ref[1]) ** 2
        / var_e3_est[None],
        axis=1,
    )
    ax.plot(t0s, chi_squared_t0s, color="k")
    t0 = t0s[np.argmin(chi_squared_t0s)]

    ax.plot(t0, chi_squared_t0s[np.argmin(chi_squared_t0s)], "x", color="r")
    if t0_offset_range_plot:
        ax.set_xlim(t0_offset_range_plot[0], t0_offset_range_plot[1])
    # ax.set_ylim([0,20]) #QUESTION: Should the user be able to put limits
    ax.set_xlim([t0s[0], t0s[-1]])
    ax.set_xlabel("offset t0 [Ma]")
    ax.set_ylabel("MSWD")
    if add_nominal:
        ax.axvline(func_time_nominal(data[0][0]), linestyle="--", color="r")

    #############################################################################
    ax = ax_time_offset
    alpha = [0.05, 0.5, 0.95]
    t_lo, t_mi, t_up = np.quantile(
        times_e3_bs[:, 0] + offset_bs[:, None], alpha, axis=0
    )
    ax.fill_between(
        data[0],
        t_lo,
        t_up,
        alpha=0.3,
        color="k",
        label=f"{(alpha[-1]-alpha[0])*100:.0f}% CI",
    )
    ax.plot(data[0], t_mi, "--", color="k", label="median")
    ax.plot(data[0], times_e3_bs[0, 0], "-", color="k", label="MAP")
    if add_nominal:
        ax.plot(
            data[0], func_time_nominal(data[0]) - t0, label="nominal", color="r", lw=1
        )
    ax.legend(loc="lower right")
    ax.set_xlabel("Depth [m]")
    ax.set_ylabel("Relative Time [Myr]")
    secax_y = ax.secondary_yaxis(
        "right", functions=(lambda x: x + t0, lambda x: x - t0)
    )
    secax_y.set_ylabel("Absolute Time [Ma]")

    # ----YW， 2025-Mar-11, start ------------------------------
    # LOGS OF THE PLOTS
    if not output_path or output_path == "default" or output_path == "":
        output_path_logs = parameters_analysis["output_folder"]
        os.makedirs(output_path_logs, exist_ok=True)
    else:
        output_path_logs = output_path
    file_path = f"{output_path_logs}/logs_plots.txt"
    plots_utils.print_and_save(
        message="###########################################", file_path=file_path
    )
    plots_utils.print_and_save(message="", file_path=file_path)
    message = f"Eccentricity Plot After MCMC {ecc_type}"
    plots_utils.print_and_save(message=message, print_time=True, file_path=file_path)
    plots_utils.print_and_save(message="", file_path=file_path)
    message = f"t0 offset range used: {t0_offset_range_plot}"
    plots_utils.print_and_save(message=message, file_path=file_path)
    message = f"Eccentricity time range used: {ecc_time_range}"
    plots_utils.print_and_save(message=message, file_path=file_path)
    message = f"Eccentricity solution used: {name_eccentricity_solution}"
    plots_utils.print_and_save(message=message, file_path=file_path)
    plots_utils.print_and_save(message="", file_path=file_path)
    message = "Time and duration with 90% confidence interval:"
    plots_utils.print_and_save(message=message, file_path=file_path)
    plots_utils.print_and_save(message="", file_path=file_path)
    # 90% CI of time scale at the beginning
    a, b, c = t0, t_lo[0] - t_mi[0], t_up[0] - t_mi[0]
    message = f"t_start(t0):    {a:.3f} ({b:.3f}, +{c:.3f})"
    plots_utils.print_and_save(message=message, file_path=file_path)
    plots_utils.print_and_save(message="", file_path=file_path)

    # 90% CI of time scale at the end
    a, b, c = times_e3_ref[0][-1] + t0, t_lo[-1] - t_mi[-1], t_up[-1] - t_mi[-1]
    message = f"t_end:          {a:.3f} ({b:.3f}, +{c:.3f})"
    plots_utils.print_and_save(message=message, file_path=file_path)
    plots_utils.print_and_save(message="", file_path=file_path)

    # 90% CI of duration
    ats = times_e3_bs[:, 0] + offset_bs[:, None]
    duration = np.array([ats[i][-1] - ats[i][0] for i in range(ats.shape[0])])
    dur_ci = np.quantile(duration, [0.05, 0.95]) - times_e3_ref[0][-1]
    a, b, c = duration[0], dur_ci[0], dur_ci[1]
    message = f"duration:       {a:.3f} ({b:.3f}, +{c:.3f})"
    plots_utils.print_and_save(message=message, file_path=file_path)
    plots_utils.print_and_save(message="", file_path=file_path)
    plots_utils.print_and_save(
        message="###########################################", file_path=file_path
    )
    ############################################################################
    ax = ax_ecc_a
    ax.plot(times_e3_ref[0] + t0, times_e3_ref[1], color="k", label="MAP")
    y_lo, y_me, y_up = np.nanquantile(aligned_e3_bs, [0.05, 0.5, 0.95], axis=0)
    ax.fill_between(
        times_e3_ref[0] + t0, y_lo, y_up, alpha=0.3, color="k", label="90% CI"
    )
    ax.plot(times_e3_ref[0] + t0, y_me, "k--", label="median")

    t_grid = np.linspace(
        *ecc_time_range, 1000
    )  # WHAT SHOULD I DO WITH THIS. IF WE DONT INPUT ecc_time_range, WHICH VALUES SHOULD I TAKE?
    ax.plot(
        t_grid,
        e3_interp(t_grid),
        color="r",
        zorder=-10,
        label=name_eccentricity_solution,
    )
    ax.set_xlim(*ecc_time_range)
    ax.legend(
        loc="upper left", bbox_to_anchor=(1, 1)
    )  # set the bbox to anchor to (1,0.5) in order to set it inside of the plot.
    r = sp.stats.pearsonr(e3_interp(times_e3_ref[0] + t0), times_e3_ref[1])[0]
    ax.text(
        0.05,
        0.95,
        f"$r = {r:.2f}$",
        transform=ax.transAxes,
        fontsize="x-large",
        verticalalignment="top",
        fontfamily="serif",
    )
    ax.set_xlabel("time [Ma]")
    if add_nominal:
        # What should we call this line?
        ax.axvline(func_time_nominal(data[0][0]), linestyle="--", color="r")
    ax.set_ylabel("eccentricity")

    if types_of_saving and types_of_saving != []:
        if not output_path or output_path == "default" or output_path == "":
            output_path = parameters_analysis["output_folder"]
            output_path = f"{output_path}/figures"
            os.makedirs(output_path, exist_ok=True)
        figname = f"Eccentricity_Plot_After_MCMC_{ecc_type}_{number_of_knots_explored}_Knots_{name_eccentricity_solution}"
        global supported_formats
        for fmt in types_of_saving:
            if fmt.lower() in supported_formats:
                if fmt.lower() == "tab":
                    logger.warning("tab format is not accepted for this plot.")
                else:
                    fig.savefig(
                        f"{output_path}/{figname}.{fmt.lower()}", format=fmt.lower()
                    )
                    logger.info(f"Saved {output_path}/{figname}.{fmt.lower()}")
            else:
                logger.warning(f"The format {fmt} is not accepted.")
    ##plt.show()
    # Please keep the commented print commands. They may be useful in the future.

    if save_eccentricity_solutions:
        column_len=26 # Don't change column_len. And don't change the header_str
        depth=data[0] # ----YW, 2025-Apr-04
        y=data[1]
        ecc_type_dict={'(prec env)':'envpre','(gi-gj)':'eccgg'}   # ----YW, 2025-Apr-04
        
        # Save reconstructed eccentricity solutions.
    
        header_str=f"{'depth(m)':^{column_len}}"+f"{'age(Ma)':^{column_len}}"+ f"{'proxy':^{column_len}}"+ f"{'ecc_MAP':^{column_len}}"+ f"{'ecc_median':^{column_len}}" + f"{'ecc_5th_percentile':^{column_len}}"+ f"{'ecc_95th_percentile':^{column_len}}" # ----YW, 2025-Apr-04
        filename_str= f"{depth[0]:.0f}m{depth[-1]:.0f}m_{ecc_type_dict[ecc_type]}_{number_of_knots_explored}k_{name_eccentricity_solution}" # ----YW, 2025-Apr-04
        comment_str=f"# Data info: eccentricity {ecc_type} solutions. "+save_eccentricity_solutions_comment+'\n'

        if aligned_e3_bs[0] is not None:
            e3_MAP = np.vstack([depth, times_e3_ref[0]+t0,y, aligned_e3_bs[0],y_me,y_lo,y_up]) # ----YW, 2025-Apr-04
        
        np.savetxt(f"{output_path}/{filename_str}.txt", e3_MAP.T, fmt='%25.18e', delimiter=' ', header=header_str, comments=comment_str)

        # Derive phase for MAP solution ====

        samples_combine = selected_sampler_MCMC.reshape(
            -1, selected_sampler_MCMC.shape[-1]
        )
        logprob_combine = selected_logprob_MCMC.reshape(
            -1, selected_logprob_MCMC.shape[-1]
        )
        index_MAP_1d = np.argmax(logprob_combine, axis=None)

        j = index_MAP_1d
        n_ff = number_free_freqs
        depth = data[0]
        coef_time = 1e6  # time unit, 1e6 yrs, namely Myr.
        arcsec_to_pi = 1 / 180 / 3600 * np.pi

        free_freqs = samples_combine[j, :n_ff]
        # print('free freqs (rad/Myr):',free_freqs)
        # print('free freqs (arcsec/yr):',free_freqs/arcsec_to_pi/coef_time)

        fs_MAP = freqs_model(free_freqs)  # Calculated the "fs"
        # print('fs_MAP (rad/Myr):',fs_MAP)
        # print('fs_MAP (arcsec/yr):',fs_MAP/arcsec_to_pi/coef_time)
        # print('period (kyr):',2*np.pi/fs_MAP*1e3)

        invSR_interpolate = interpolator(
            [depth_genes, samples_combine[j, n_ff:]], depth
        )
        times_inferred_MAP = sp.integrate.cumulative_trapezoid(
            invSR_interpolate, depth, initial=0
        )

        [phi_MAP, beta_sin_MAP, beta_cos_MAP] = get_phase_beta(
            y, fs, times_inferred_MAP
        )
        # print('phi_MAP:',phi_MAP)
        # print('beta_sin_MAP:',beta_sin_MAP)
        # print('beta_cos_MAP:',beta_cos_MAP)

        # Save parameters to reconstruct eccentricity solutions. Eq.(A2)
        header_str = (
            f"{'period_j(kyr)':^{column_len}}"
            + f"{'nu_j(arcsec/yr)':^{column_len}}"
            + f"{'nu_j(rad/Myr)':^{column_len}}"
            + f"{'beta_j':^{column_len}}"
            + f"{'beta_(j+k)':^{column_len}}"
            + f"{'phi_j':^{column_len}}"
        )
        comment_str = (
            f"# Data info: these amplitude (beta_j,beta_(j+k)) and frequency (nu_j) data are corresponding to the eccentricity {ecc_type} solutions. "
            + save_eccentricity_solutions_comment
            + "\n"
        )

        e3_sol_MAP = np.vstack(
            [
                2 * np.pi / fs_MAP * 1e3,
                fs_MAP / arcsec_to_pi / coef_time,
                fs_MAP,
                beta_sin_MAP,
                beta_cos_MAP,
                phi_MAP,
            ]
        )
        np.savetxt(
            f"{output_path}/{filename_str}.sol",
            e3_sol_MAP.T,
            fmt="%25.18e",
            delimiter=" ",
            header=header_str,
            comments=comment_str,
        )
        # This file type should be '.sol', don't change it.


############ MISC SECTION ###################
def plt_mcmc_summary_of_results(
    configuration_file_path:str,
    folder_with_the_results:Optional[str],
    file_with_the_results:Optional[str],
    x_axis_limits_peridogram: Optional[list],
    y_axis_limits_peridogram: Optional[list],
    peridogram_scale_x_axis: Optional[str],
    peridogram_scale_y_axis: Optional[str],
    x_axis_limits_depth: Optional[list],
    ignore_weights:bool,
    types_of_saving: list,
    output_path: str,
):
    """
    Generate a comprehensive summary plot from MCMC results in astrochronological analysis.

    This function plots five panels:
        1. Spectral power analysis of log-likelihood for each number of knots.
        2. Median time model with 95% confidence interval.
        3. Inverse sedimentation rate over time with credible intervals.
        4. Fit of data in the time domain with posterior predictive samples.
        5. Decomposition of the ETP (Eccentricity, Tilt, Precession) signal using fitted knots.

    The figure is saved in the specified formats and the output path is determined by the provided configuration.

    Parameters:
        - **configuration_file_path** : str  
            Path to the configuration file used in the analysis.

        - **folder_with_the_results** : str or None  
            Path to the folder containing the MCMC result files. Can be `None` if `file_with_the_results` is provided.

        - **file_with_the_results** : str or None  
            Specific path to a result file. Can be `None` if `folder_with_the_results` is provided.

        - **x_axis_limits_peridogram** : list or None  
            Limits for the x-axis in the peridogram plot. Should be [min, max] in frequency units. If `None`, defaults are used.

        - **y_axis_limits_peridogram** : list or None  
            Limits for the y-axis in the peridogram plot. Should be [min, max]. If `None`, auto-scaled from data.

        - **peridogram_scale_x_axis** : str or None  
            Scale for the x-axis of the peridogram plot. Can be `"linear"` or `"log"`. Defaults to `"linear"`.

        - **peridogram_scale_y_axis** : str or None  
            Scale for the y-axis of the peridogram plot. Can be `"linear"` or `"log"`. Defaults to `"linear"`.

        - **x_axis_limits_depth** : list or None  
            Limits for the x-axis (depth) used in multiple plots. Should be [min_depth, max_depth]. If `None`, it spans the dataset range.

        - **types_of_saving** : list of str  
            List of file formats to save the final plot. Supported formats include `"png"`, `"pdf"`, and `"svg"`. `"tab"` is ignored.

        - **output_path** : str  
            Path where the output plot should be saved. If `"default"`, it's inferred from the configuration file.


    Raises:
        FileNotFoundError: If the results directory or any expected file is missing.
        ValueError: If the configuration file is not formatted correctly or missing required fields.
        RuntimeError: If the data dimensions do not match expected shapes or plotting fails.

    Example:
        >>> plt_mcmc_summary_of_results("configs/my_run_config.yaml", types_of_saving=["png", "pdf"])
        # Saves the plot as "Misc_Plots_<K>_Knots.png" and ".pdf" in the output path.
    """
    ########## LIBRARY MODE SECTION ##############
    
    _, _, _, parameters_analysis, data_set_values, age_depth_model_values, genetic_algorithm_parameters, data_model_parameters, _, _ = notebooks_code.data_obtain_analysis_file_and_result_files(True, configuration_file_path, folder_with_the_results, file_with_the_results)

    _, _, _, func_invSR_nominal, func_time_nominal, add_nominal = notebooks_code.data_obtain_age_depth_values(age_depth_model_values)

    data, _ = data_manager.data_read_dataset_from_file(data_set_values)
    
    frequency_values=data_model_parameters["frequency_values"]
    colors_science = [
            "#4477AA",
            "#EE6677",
            "#228833",
            "#CCBB44",
            "#66CCEE",
            "#AA3377",
            "#BBBBBB",
        ]
    tool_plt_mcmc_summary_of_results(
        data_model_parameters=data_model_parameters,
        data=data,
        parameters_analysis=parameters_analysis,
        colors_science=colors_science,
        func_invSR_nominal=func_invSR_nominal,
        func_time_nominal=func_time_nominal,
        x_axis_limits_peridogram=x_axis_limits_peridogram,
        y_axis_limits_peridogram=y_axis_limits_peridogram,
        peridogram_scale_x_axis=peridogram_scale_x_axis,
        peridogram_scale_y_axis=peridogram_scale_y_axis,
        x_axis_limits_depth=x_axis_limits_depth,
        ignore_weights=ignore_weights,
        seed = genetic_algorithm_parameters["seed"],
        types_of_saving=types_of_saving,
        output_path=output_path,
        add_nominal=add_nominal
    )

def tool_plt_mcmc_summary_of_results(
    data_model_parameters:dict,
    data: list,
    parameters_analysis: dict,
    colors_science,
    func_invSR_nominal,
    func_time_nominal,
    x_axis_limits_peridogram: Optional[list],
    y_axis_limits_peridogram: Optional[list],
    peridogram_scale_x_axis: Optional[str],
    peridogram_scale_y_axis: Optional[str],
    x_axis_limits_depth: Optional[list],
    ignore_weights:bool,
    seed,
    types_of_saving: list,
    output_path: str,
    add_nominal: bool
):
    frequency_values=data_model_parameters["frequency_values"]
    
    dict_results_gene_selection = data_manager.load_dictionary_from_pickle(
        f"{os.getcwd()}/tmp/temp_selection_genes_mcmc.pickle"
    )

    selected_opt_result = dict_results_gene_selection["selected_opt_result"]
    selected_sampler_MCMC = dict_results_gene_selection["selected_sampler_MCMC"]
    selected_logprob_MCMC = dict_results_gene_selection["selected_logprob_MCMC"]
    ind_MAP = dict_results_gene_selection["ind_MAP"]
    params_MAP = dict_results_gene_selection["params_MAP"]
    invSR_MAP = dict_results_gene_selection["invSR_MAP"]
    calculated_frequencies_map = dict_results_gene_selection["fs"]
    y_pred = dict_results_gene_selection["y_pred"]
    times_e3_bs = dict_results_gene_selection["times_e3_bs"]
    ip = dict_results_gene_selection["ip"]
    ie = dict_results_gene_selection["ie"]
    it = dict_results_gene_selection["it"]
    times_inferred = dict_results_gene_selection["times_inferred"]
    number_of_knots_explored = dict_results_gene_selection["number_of_knots_to_explore"]
    weight_results = dict_results_gene_selection["weight_results"]

    atributes = plots_utils.get_result_atributes(selected_opt_result)
    
    if seed:
        if isinstance(seed, int):
            setup_seed(seed)
        elif isinstance(seed, list):
            setup_seed(seed[0])

    fig, axs = plt.subplot_mosaic(
        [["a)"], ["b)"], ["c)"], ["d)"], ["e)"]],
        constrained_layout=True,
        figsize=(12, 10),
    )

    ax_spectrum = axs["a)"]
    ## THESE FOUR SHARE X AXIS LABEL
    ax_time_relative = axs["b)"]
    ax_invSR = axs["c)"]
    ax_data = axs["d)"]
    ax_decompose = axs["e)"]
    
    for label, ax in axs.items():
        # label physical distance in and down:
        # xtext = 0
        trans = mtransforms.ScaledTranslation(10 / 72, -5 / 72, fig.dpi_scale_trans)
        ax.text(
            0,
            1.0,
            label,
            transform=ax.transAxes + trans,
            fontsize="medium",
            verticalalignment="top",
            fontfamily="serif",
            bbox=dict(facecolor="0.7", edgecolor="none", pad=3.0),
        )

    ###########################################################
    ################### SPECTRUM PLOT #########################
    ###########################################################

    if x_axis_limits_peridogram:
        ax_spectrum.set_xlim(x_axis_limits_peridogram)
    if y_axis_limits_peridogram:
        ax_spectrum.set_ylim(
            y_axis_limits_peridogram[0], y_axis_limits_peridogram[1] * 1.9
        )
        y_min, y_max = ax_spectrum.get_ylim()
        max_height = y_max - y_min

    if not peridogram_scale_x_axis or peridogram_scale_x_axis == "linear":
        ax_spectrum.set_xscale("linear")
    elif peridogram_scale_x_axis == "log":
        ax_spectrum.set_xscale("log")

    if not peridogram_scale_y_axis or peridogram_scale_y_axis == "linear":
        ax_spectrum.set_yscale("linear")
    elif peridogram_scale_y_axis == "log":
        ax_spectrum.set_yscale("log")

    ax_spectrum.set_xlabel("frequency (cycle/Myr)")

    times_grid = np.linspace(times_inferred[0], times_inferred[-1], len(times_inferred))
    ydata_sp = interp1d(times_inferred, data[1])
    y_equ_spa = ydata_sp(times_grid)

    dt = times_grid[1] - times_grid[0]

    f, Pxx_den = sp.signal.periodogram(y_equ_spa, 1 / dt, nfft=2**12)
    if not y_axis_limits_peridogram:
        max_height = np.max(Pxx_den)
        ax_spectrum.set_ylim(0, max_height * 2.4)
    adjusted_height_frequencies = max_height * 1.5
    adjusted_height_labels = max_height * 0.7
    ax_spectrum.plot(f, Pxx_den, color="k")
    secax = ax_spectrum.secondary_xaxis("top", functions=(lambda f: f, lambda p: p))
    formatter = FuncFormatter(lambda x_val, tick_pos: f"{1/x_val*1000:.1f}")
    secax.xaxis.set_major_formatter(formatter)
    secax.set_xlabel("period (kyr)")
    
    fs_labels = []
    prec_labels = [f"$g_{i+1}+k$" for i in range(5)]
    ecc_labels = [
            f"$g_2-g_5$",
            f"$g_4-g_5$",
            f"$g_4-g_2$",
            f"$g_3-g_5$",
            f"$g_3-g_2$",
        ]
    til_labels = ["$s_3+k$", "$s_4+k$", "$s_6+k$"]
    if frequency_values["use_precession"]:
        fs_labels.extend(prec_labels)
        ax_spectrum.annotate(
            "climatic \n precession",
            xy=(49, adjusted_height_labels),
            ha="center",
            va="bottom",
            rotation="vertical",
            color=colors_science[3],
        )
        
    if frequency_values["use_eccentricity"]:
        fs_labels.extend(ecc_labels)
        ax_spectrum.annotate(
            "eccentricity",
            xy=(6, adjusted_height_labels),
            ha="center",
            va="bottom",
            rotation="vertical",
            color=colors_science[2],
        )
        
    if frequency_values["use_tilt"]:
        fs_labels.extend(til_labels)
        ax_spectrum.annotate(
            "tilt",
            xy=(22, adjusted_height_labels),
            ha="center",
            va="bottom",
            rotation="vertical",
            color=colors_science[4],
        )
    print(fs_labels)
    for i, (f, f_label) in enumerate(zip(calculated_frequencies_map, fs_labels)):
        dx = -5 if f_label in [f"$g_3-g_2$", f"$g_3-g_5$", "$g_5+k$", "$g_3+k$", "$s_3+k$"] else 2
        if f_label == "$g_1+k$": dx=0
        if f_label in prec_labels:
            color = colors_science[3]
        elif f_label in ecc_labels:
            color = colors_science[2]
        elif f_label in til_labels:
            color = colors_science[4]
            
        if f_label == "$g_4+k$":
            ax_spectrum.annotate(
                f"{f_label}",
                xy=((f / 2 / np.pi) + 0.4, adjusted_height_frequencies),
                ha="center",
                va="bottom",
                xytext=(dx + 0.4, adjusted_height_frequencies),
                textcoords="offset points",
                rotation="vertical",
                backgroundcolor="white",
                color=color,
                zorder=-i,
            )
        else:
            ax_spectrum.annotate(
                f"{f_label}",
                xy=((f / 2 / np.pi), adjusted_height_frequencies),
                ha="center",
                va="bottom",
                xytext=(dx, adjusted_height_frequencies),
                textcoords="offset points",
                rotation="vertical",
                backgroundcolor="white",
                color=color,
                zorder=-i,
            )
        ax_spectrum.axvline(
            f / 2 / np.pi, linestyle="--", color=color, zorder=-100, alpha=0.5
        )
        
    ax_spectrum.tick_params(labelbottom=True)
    #####################################################################################

    ################################################################
    ################### TIME RELATIVE PLOT #########################
    ################################################################
    alpha = [0.05, 0.5, 0.95]
    t_lo, t_mi, t_up = np.quantile(times_e3_bs[:, 0], alpha, axis=0)
    ax_time_relative.fill_between(
        data[0], t_lo, t_up, alpha=0.3, color="k", label="90% CI"
    )
    ax_time_relative.plot(data[0], t_mi, "--", color="k", label="median")
    if add_nominal:
        ax_time_relative.plot(
            data[0],
            func_time_nominal(data[0]) - func_time_nominal(data[0][0]),
            label="nominal",
            color="r",
            lw=1,
        )
    ax_time_relative.plot(data[0], times_e3_bs[0, 0], "-", color="k", label="MAP")
    ax_time_relative.legend(ncol=2, loc="lower right")
    ax_time_relative.set_ylabel("Time [Myr]")
    ax_time_relative.tick_params(labelbottom=False)
    ax_time_relative.sharex(ax_decompose)

    #######################################################################

    #############################################################################
    ################### INVERSE SEDIMENTATION RATE PLOT #########################
    #############################################################################

    n = 1000  # CHANGE: No hardcoded values
    if not ignore_weights  and isinstance(weight_results,np.ndarray):
        inds_1 = np.random.choice(
            range(len(selected_sampler_MCMC)), size=n, replace=True, p=weight_results
        )  # AS, 2025-Feb-12
    else:
        inds_1 = np.random.choice(
            range(len(selected_sampler_MCMC)), size=n, replace=True
        )
    inds_2 = np.random.choice(
        range(len(selected_sampler_MCMC[0])), size=n, replace=True
    )

    invSR_interpolate_par = np.array(
        [
            selected_opt_result.problem.interpolator(
                [
                    selected_opt_result.problem.depth_genes,  # CHANGE: See how to import an interpolator
                    selected_sampler_MCMC[j, i, atributes["number_free_freqs"] :],
                ],
                data[0],
            )
            for j, i in zip(inds_1, inds_2)
        ]
    )

    mean_invSR = np.median(invSR_interpolate_par, axis=0)
    alpha = [0.05, 0.95]  # CHANGE: No hardcoded values
    invSR_lo, invSR_up = np.quantile(invSR_interpolate_par, alpha, axis=0)

    ax_invSR.plot(data[0], mean_invSR, "--", color="black", label="median", zorder=10)
    ax_invSR.fill_between(
        data[0],
        invSR_lo,
        invSR_up,
        alpha=0.3,
        color="black",
        label=f"{(alpha[1]-alpha[0])*100:.0f}% CI",
    )
    if add_nominal:
        ax_invSR.plot(
            data[0],
            func_invSR_nominal(data[0]),
            color="red",
            zorder=-1000,
            label="nominal",
        )

    ind_MAP = np.unravel_index(
        np.argmax(selected_logprob_MCMC, axis=None),
        selected_logprob_MCMC.shape,
    )
    invSR_MAP = atributes["interpolator"](
        [
            atributes["depth_genes"],
            selected_sampler_MCMC[ind_MAP][atributes["number_free_freqs"] :],
        ],
        data[0],
    )
    ax_invSR.plot(data[0], invSR_MAP, color="k", label="MAP")
    ax_invSR.tick_params(labelbottom=False)
    ax_invSR.legend(loc="upper right")
    ax_invSR.set_ylabel("SR$^{-1}$ [Myr/m]")
    ax_invSR.sharex(ax_decompose)

    ###########################################################################

    ###########################################################################
    ################### DATA VS RECOVERED SIGNAL PLOT #########################
    ###########################################################################

    ax_data.plot(data[0], data[1], color="k", label="Original Data")
    ax_data.plot(data[0], y_pred, color="r", label="Recovered Signal")

    ax_data.set_ylabel("Proxy Value")
    ax_data.legend(loc="upper right")
    R2 = r2_score(data[1], y_pred)
    ax_data.tick_params(labelbottom=False)
    ax_data.text(
        0.1,
        0.95,
        f"$R^2$ = {R2:.2f}",
        transform=ax_data.transAxes,
        fontsize="x-large",
        verticalalignment="top",
        fontfamily="serif",
    )
    ax_data.sharex(ax_decompose)

    ###########################################################################

    #######################################################################
    ################### DECOMPOSE ETP SIGNAL PLOT #########################
    #######################################################################
    res_ETP = plots_utils.get_ETP_components(
        params_MAP,
        atributes["depth_genes"],
        atributes["frequency_model"],
        atributes["inverse_SR_lims"],
        atributes["interpolator"],
        data,
        ip,
        ie,
        it,
    )
    all_y_values = []

    if frequency_values["use_precession"]:
        y_prec = res_ETP["y_prec"]
        y_env_prec = res_ETP["y_env_prec"]
        ax_decompose.plot(data[0], y_prec, label="precession", color=colors_science[3])
        ax_decompose.plot(data[0], y_env_prec, label="envelope", color=colors_science[0])
        all_y_values.extend(y_prec)
        all_y_values.extend(y_env_prec)

    if frequency_values["use_eccentricity"]:
        y_ecc = res_ETP["y_ecc"] + 3
        ax_decompose.plot(data[0], y_ecc, label="eccentricity", color=colors_science[2])
        all_y_values.extend(y_ecc)

    if frequency_values["use_tilt"]:
        y_tilt = res_ETP["y_tilt"] - 2
        ax_decompose.plot(data[0], y_tilt, label="tilt", color=colors_science[-3])
        all_y_values.extend(y_tilt)

    # Set x-axis limits
    ax_decompose.set_xlim([data[0][0], data[0][-1]])
    if x_axis_limits_depth:
        ax_decompose.set_xlim(x_axis_limits_depth)

    # Dynamically set y-limits
    if all_y_values:
        min_y = np.min(all_y_values)
        max_y = np.max(all_y_values)
        buffer = 0.4 * (max_y - min_y)  # 10% buffer for spacing
        ax_decompose.set_ylim(min_y - buffer, max_y + buffer)

    # Final labels and legend
    ax_decompose.set_xlabel("Depth [m]")
    ax_decompose.legend(loc="lower right")

    if types_of_saving and types_of_saving != []:
        if not output_path or output_path == "default" or output_path == "":
            output_path = parameters_analysis["output_folder"]
            output_path = f"{output_path}/figures"
            os.makedirs(output_path, exist_ok=True)
        figname = f"Misc_Plots_{number_of_knots_explored}_Knots"
        global supported_formats
        for fmt in types_of_saving:
            if fmt.lower() in supported_formats:
                if fmt.lower() == "tab":
                    logger.warning("tab format is not accepted for this plot.")
                else:
                    fig.savefig(
                        f"{output_path}/{figname}.{fmt.lower()}", format=fmt.lower()
                    )
                    logger.info(f"Saved {output_path}/{figname}.{fmt.lower()}")
            else:
                logger.warning(f"The format {fmt} is not accepted.")
    ##plt.show()
    
###########################################################################

###########################################################################
############################ FREQUENCIES PLOT #############################
def plt_mcmc_prior_frequencies_distributions(
    configuration_file_path:str,
    folder_with_the_results:Optional[str],
    file_with_the_results:Optional[str],
    ignore_weights:bool,
    types_of_saving: list,
    output_path: str,
):
    """
    Plot the prior and posterior distributions of the frequencies explored by MCMC.

    This function generates subplots of the marginal prior and posterior distributions 
    for the different orbital frequencies (precession, obliquity, etc.) used in the 
    astrochronological model. It overlays the MAP value, the current astronomical values, 
    and posterior statistics such as the weighted median and standard deviation.
    
    It also logs summary statistics (MAP, median, std) for each frequency and saves the plot 
    in one or more formats.

    Parameters:
        - **configuration_file_path** : str
            Path to the YAML configuration file used to run the analysis.
        
        - **folder_with_the_results** : str or None
            Path to the folder containing the result files. If `file_with_the_results` is given,
            this can be `None`.

        - **file_with_the_results** : str or None
            Path to a specific result file. If provided, `folder_with_the_results` can be `None`.

        - **types_of_saving** : list
            List of file formats to save the plot in. Supported formats are image types like 
            ['png', 'pdf', 'svg'], etc. 'tab' is ignored with a warning.

        - **output_path** : str
            Directory where the figure(s) and logs should be saved. If set to "default", uses the 
            default output folder specified in the configuration file.

    Returns:
        None
            The function saves the figures and logs to the specified file system path but does not 
            return any object.

    Notes:
        - The function expects a temporary file at `tmp/temp_selection_genes_mcmc.pickle` that 
        stores selected MCMC results and MAP parameters.
        - The number of subplots and the layout depend on the number of free frequencies being estimated.
        - Each subplot shows:
            - Prior distribution (if available)
            - Posterior histogram
            - MAP estimate (black line)
            - Astronomical (reference) frequency (red line)
            - Posterior median (dashed line) and weighted std if weights are available

        - Summary statistics for each frequency (MAP, median, std) are saved in a log file 
        `logs_plots.txt`.

    Warnings:
    - If `types_of_saving` includes unsupported formats (e.g., 'tab'), they will be ignored with a warning.
    - If the number of free frequencies is zero, the function will silently exit without plotting.

    """
    ########## LIBRARY MODE SECTION ##############
    
    results_fit, results_mcmc, weight_results, parameters_analysis, _, _, genetic_algorithm_parameters, _, mcmc_parameters, _ = notebooks_code.data_obtain_analysis_file_and_result_files(True, configuration_file_path, folder_with_the_results, file_with_the_results)
    tool_plt_mcmc_prior_frequencies_distributions(
        results_fit=results_fit,
        results_mcmc=results_mcmc,
        ignore_weights = ignore_weights,
        mcmc_parameters=mcmc_parameters,
        seed = genetic_algorithm_parameters["seed"],
        parameters_analysis=parameters_analysis,
        types_of_saving=types_of_saving,
        output_path=output_path
    )

def tool_plt_mcmc_prior_frequencies_distributions(
    results_fit:dict,
    results_mcmc:dict,
    mcmc_parameters:dict,
    seed,
    ignore_weights:bool,
    parameters_analysis: dict,
    types_of_saving: list,
    output_path: str,
    script_path: str = "",
):
    if seed:
        if isinstance(seed, int):
            setup_seed(seed)
        elif isinstance(seed, list):
            setup_seed(seed[0])
    prior_params_freq=results_mcmc["prior_frequencies"]
    prior_distributions=mcmc_parameters["prior_distributions"]
    conversion_factor=results_fit["convers_factor"]
    
    dict_results_gene_selection = data_manager.load_dictionary_from_pickle(
        f"{os.getcwd()}/tmp/temp_selection_genes_mcmc.pickle"
    )
    selected_opt_result = dict_results_gene_selection["selected_opt_result"]
    selected_sampler_MCMC = dict_results_gene_selection["selected_sampler_MCMC"]
    params_MAP = dict_results_gene_selection["params_MAP"]
    atributes = plots_utils.get_result_atributes(selected_opt_result)
    number_of_knots_explored = dict_results_gene_selection["number_of_knots_to_explore"]
    weight_results = dict_results_gene_selection["weight_results"]
    
    current_p0 = 50.467718
    current_gi = [5.579378, 7.456665, 17.366595, 17.910194]
    current_si = [-18.845166, -17.758310]
    current_freqs = np.array([current_p0] +current_gi + current_si)
    print(current_freqs)
    name_freqs = (
        ["$k$"] + [f"$g_{i+1}$" for i in range(4)] + [f"$s_{i}$" for i in [3, 4]]
    )
    letters_plots = ["a)", "b)", "c)", "d)", "e)", "f)", "g)"]
    if atributes["number_free_freqs"] == 0:
        return
    # FRQUENCIES SECTION
    distributions = []
    # ALL FREQUENCIES: P0, GI AND SI
    import math

    n = atributes["number_free_freqs"]
    ncols = 2
    nrows = math.ceil(n / ncols)

    # Create subplot labels dynamically, e.g., ["a)", "b)", ..., up to needed number]
    letters_plots = [f"{chr(97+i)})" for i in range(n)]

    # Build mosaic layout dynamically
    mosaic_layout = [letters_plots[i * ncols:(i + 1) * ncols] for i in range(nrows)]
    empty_sentinel = ""
    if len(mosaic_layout[-1]) < ncols:
        mosaic_layout[-1].append(empty_sentinel)

    figsize = (12, 3 * nrows)  # Adjust height dynamically based on number of rows
    fig, axs = plt.subplot_mosaic(
        mosaic_layout, constrained_layout=True, figsize=figsize
    )

    if atributes["number_free_freqs"] == 7:
        distributions.append(prior_distributions["p0_distribution"])
        distributions.extend([prior_distributions["gi_distribution"]] * 4)
        distributions.extend([prior_distributions["si_distribution"]] * 2)
        name_freqs = (
            ["$k$"] + [f"$g_{i+1}$" for i in range(4)] + [f"$s_{i}$" for i in [3, 4]]
        )
    # GI AND SI
    elif atributes["number_free_freqs"] == 6:
        distributions.extend([prior_distributions["gi_distribution"]] * 4)
        distributions.extend([prior_distributions["si_distribution"]] * 2)
        name_freqs = (
            [f"$g_{i+1}$" for i in range(4)] + [f"$s_{i}$" for i in [3, 4]]
        )
    # P0 AND GI
    elif atributes["number_free_freqs"] == 5:
        distributions.append(prior_distributions["p0_distribution"])
        distributions.extend([prior_distributions["gi_distribution"]] * 4)
        name_freqs = (
            ["$k$"] + [f"$g_{i+1}$" for i in range(4)]
        )
    # GI
    elif atributes["number_free_freqs"] == 4:
        distributions.extend([prior_distributions["gi_distribution"]] * 4)
        name_freqs = (
            [f"$g_{i+1}$" for i in range(4)]
        )
    # P0 AND SI
    elif atributes["number_free_freqs"] == 3:
        distributions.append(prior_distributions["p0_distribution"])
        distributions.extend([prior_distributions["si_distribution"]] * 2)
        name_freqs = (
            ["$k$"] + [f"$s_{i}$" for i in [3, 4]]
        )
    # SI
    elif atributes["number_free_freqs"] == 2:
        distributions.extend([prior_distributions["si_distribution"]] * 2)
        name_freqs = (
            [f"$s_{i}$" for i in [3, 4]]
        )
    # P0
    elif atributes["number_free_freqs"] == 1:
        distributions.append(prior_distributions["p0_distribution"])
        name_freqs = (["$k$"])
    # LOGS OF THE PLOTS
    if not output_path or output_path == "default" or output_path == "":
        output_path_logs = parameters_analysis["output_folder"]
        os.makedirs(output_path_logs, exist_ok=True)
    else:
        output_path_logs = output_path
    file_path = f"{output_path_logs}/logs_plots.txt"
    plots_utils.print_and_save(
        message="###########################################", file_path=file_path
    )
    plots_utils.print_and_save(message="", file_path=file_path)
    message = f"Frequency Plots"
    plots_utils.print_and_save(message=message, print_time=True, file_path=file_path)
    plots_utils.print_and_save(message="", file_path=file_path)
    message = " MAP" + "       " + "median" + "     " + "std"
    plots_utils.print_and_save(message=message, file_path=file_path)
    plots_utils.print_and_save(message="", file_path=file_path)
    for i in range(atributes["number_free_freqs"]):
        axs_freq = axs[letters_plots[i]]
        prior_params = prior_params_freq[i].copy() * conversion_factor
        sample_freq = selected_sampler_MCMC[:, :, i].flatten() * conversion_factor

        if distributions[i] == "gaussian":
            xplot_lims = (
                prior_params[0] - prior_params[1] * 5,
                prior_params[0] + prior_params[1] * 5,
            )
            prior_freq = sp.stats.norm
            xgrid = np.linspace(*xplot_lims, 1000)
            axs_freq.plot(xgrid, prior_freq.pdf(xgrid, *prior_params), label="prior")
            axs_freq.set_xlim(*xplot_lims)
        else:
            y_min, y_max = axs_freq.get_ylim()  # Get current y-axis limits
            y_pos = y_min + (y_max - y_min) * 0.25  # Compute 1/4th of the height
            axs_freq.axhspan(0, y_pos, alpha=0.5, label="prior")
        axs_freq.axvline(params_MAP[i] * conversion_factor, label="MAP", color="k")
        axs_freq.axvline(current_freqs[i], label="current", color="r")
        #######################################################
        if not ignore_weights and isinstance(weight_results, np.ndarray):
            axs_freq.axvline(
                utils_post_mcmc.quantile(
                    sample_freq,
                    0.5,
                    weights=np.repeat(weight_results, selected_sampler_MCMC.shape[1]),
                ),
                linestyle="--",
                color="k",
                label="median",
            )
            axs_freq.hist(
                sample_freq,
                density=True,
                bins=100,
                label="posterior",
                color="gray",
                alpha=0.7,
                weights=np.repeat(weight_results, selected_sampler_MCMC.shape[1]),
            )  # AS, 2025-Feb-12
            axs_freq.set_xlabel(f"{name_freqs[i]} [arcsec/yr]")
            a = params_MAP[i] * conversion_factor
            b = utils_post_mcmc.quantile(
                sample_freq,
                0.5,
                weights=np.repeat(weight_results, selected_sampler_MCMC.shape[1]),
            )
            c = utils_post_mcmc.weighted_std(
                sample_freq, np.repeat(weight_results, selected_sampler_MCMC.shape[1])
            )

        else:
            axs_freq.axvline(
                np.median(sample_freq), linestyle="--", color="k", label="median"
            )
            axs_freq.hist(
                sample_freq,
                density=True,
                bins=50,
                label="posterior",
                color="gray",
                alpha=0.7,
            )
            axs_freq.set_xlabel(f"{name_freqs[i]} [arcsec/yr]")

            a, b, c = (
                params_MAP[i] * conversion_factor,
                np.median(sample_freq),
                np.std(sample_freq),
            )

        message = f"{a:>7.3f}   {b:>7.3f}   {c:>7.3f}"
        plots_utils.print_and_save(message=message, file_path=file_path)
        plots_utils.print_and_save(message="", file_path=file_path)
    # Delete the last plot in case it is empty
    
    if atributes["number_free_freqs"] % 2 == 1:
        axs[""].remove()
    
    plots_utils.print_and_save(
        message="###########################################", file_path=file_path
    )
    plots_utils.print_and_save(message="", file_path=file_path)
    for label, ax in axs.items():
        # label physical distance in and down:
        # xtext = 0
        trans = mtransforms.ScaledTranslation(10 / 72, -5 / 72, fig.dpi_scale_trans)
        ax.text(
            0,
            1.0,
            label,
            transform=ax.transAxes + trans,
            fontsize="medium",
            verticalalignment="top",
            fontfamily="serif",
            bbox=dict(facecolor="0.7", edgecolor="none", pad=3.0),
        )

    # axs[0].set_xlim(50,52)
    axs["a)"].legend()
    if types_of_saving and types_of_saving != []:
        if not output_path or output_path == "default" or output_path == "":
            output_path = parameters_analysis["output_folder"]
            output_path = f"{output_path}/figures"
            os.makedirs(output_path, exist_ok=True)
        figname = f"Frequency_Plots_{number_of_knots_explored}_Knots"
        global supported_formats
        for fmt in types_of_saving:
            if fmt.lower() in supported_formats:
                if fmt.lower() == "tab":
                    logger.warning("tab format is not accepted for this plot.")
                else:
                    fig.savefig(
                        f"{output_path}/{figname}.{fmt.lower()}", format=fmt.lower()
                    )
                    logger.info(f"Saved {output_path}/{figname}.{fmt.lower()}")
            else:
                logger.warning(f"The format {fmt} is not accepted.")
    # plt.show()

def plt_mcmc_phase_of_frequencies(
    configuration_file_path:str,
    folder_with_the_results:Optional[str],
    file_with_the_results:Optional[str],
    types_of_saving: list, 
    output_path: str
):
    """
    Generate and save a multi-panel plot showing the phases of orbital frequencies 
    derived from MCMC solutions applied to astronomical time scale inference.

    This function loads the results of an MCMC analysis (using AstroGeoFit), computes
    the phase of orbital frequencies (such as eccentricity, obliquity, and precession),
    and plots histograms of the sampled phase distributions along with the MAP (maximum
    a posteriori) and median values. It also logs summary statistics to a file.

    Parameters:
        - **configuration_file_path** : str
            Path to the configuration file used for the AstroGeoFit model run.
        
        - **folder_with_the_results** : str or None
            Path to the folder containing the output results of the MCMC run.
            Can be None if `file_with_the_results` is specified instead.
        
        - **file_with_the_results** : str or None
            Path to a specific result file to load instead of an entire folder.
            Can be None if `folder_with_the_results` is provided.

        - **types_of_saving** : list of str
            A list of file formats to save the final phase plot (e.g., ['png', 'pdf']).
            If 'tab' is included, it is ignored for this plot type.

        - **output_path** : str
            Directory path to store the generated plot and the log file. If empty or set 
            to "default", the output directory defined in the configuration will be used.

    Returns:
        None
            The function generates and saves a matplotlib figure showing the phase histograms 
            of the orbital frequencies. It also creates a log file with MAP, median, mean, 
            and standard deviation of the phases in units of π.
    
    Notes:
        - The orbital frequencies and timescales are inferred from the MCMC samples.
        - The phase is calculated using a linear regression on the sinusoidal components 
        of the orbital forcing signal.
        - Histogram axes are adjusted to align the median to be centered in phase space.
        - The function expects a temporary file `temp_selection_genes_mcmc.pickle` to 
        exist in the working directory, containing the best MCMC samples and associated data.

    """
    _, _, _, parameters_analysis, data_set_values, _, _, _, _, _ = notebooks_code.data_obtain_analysis_file_and_result_files(True, configuration_file_path, folder_with_the_results, file_with_the_results)

    data, _ = data_manager.data_read_dataset_from_file(data_set_values)
    tool_plt_mcmc_phase_of_frequencies(
        data=data,
        parameters_analysis=parameters_analysis,
        types_of_saving=types_of_saving,
        output_path=output_path
    )
    
def tool_plt_mcmc_phase_of_frequencies(
    data: list, parameters_analysis: dict, types_of_saving: list, output_path: str
):
    def get_phase(y, fs, times_inferred):
        """
        Input: y: the observed data
            fs: the orbital frequencies, such as eccentricity, obliquity, and climatic precession
            times_infferred: the timeMAP scale
        Output: phi, an array which contains the phases of orbital frequencies
        """
        X = shared_functions.generate_X_linReg(
            np.ones_like(fs), fs, times_inferred
        )  # Generate X
        reg_model = LinearRegression().fit(X, y)  # Get the regression funciton.
        beta_sin, beta_cos = (
            reg_model.coef_[: len(fs)],
            reg_model.coef_[len(fs) :],
        )  # Get coefficients
        phi = np.arctan2(-beta_sin, beta_cos)  # Calculate the phase
        return phi

    def get_phase_beta(y, fs, times_inferred):
        """
        Input: y: the observed data
            fs: the orbital frequencies, such as eccentricity, obliquity, and climatic precession
            times_infferred: the time scale
        Output: phi, an array which contains the phases of orbital frequencies
        """
        X = shared_functions.generate_X_linReg(
            np.ones_like(fs), fs, times_inferred
        )  # Generate X
        reg_model = LinearRegression().fit(X, y)  # Get the regression funciton.
        beta_sin, beta_cos = (
            reg_model.coef_[: len(fs)],
            reg_model.coef_[len(fs) :],
        )  # Get coefficients
        phi = np.arctan2(-beta_sin, beta_cos)  # Calculate the phase
        return phi, beta_sin, beta_cos

    def get_phase_median(phi_samples, x_hist, y_hist):
        phi_samples_new = np.zeros(phi_samples.shape)
        index_cut = np.argmin(y_hist)
        phase_cut = x_hist[index_cut]

        for i in range(len(phi_samples)):
            if (phi_samples[i] >= phase_cut) and (phase_cut >= 0):
                phi_samples_new[i] = phi_samples[i] - 2 * np.pi
            elif (phi_samples[i] < phase_cut) and (phase_cut < 0):
                phi_samples_new[i] = phi_samples[i] + 2 * np.pi
            else:
                phi_samples_new[i] = phi_samples[i]

        return phi_samples_new

    def trim_axs(axs, N):
        """
        Reduce *axs* to *N* Axes. All further Axes are removed from the figure.
        """
        axs = axs.flat
        for ax in axs[N:]:
            ax.remove()
        return axs[:N]
    
    def check_median_2pi(m):
        if m>np.pi:
            m=m-2*np.pi
        elif m<=-np.pi:
            m=m+2*np.pi
        return m
            
    dict_results_gene_selection = data_manager.load_dictionary_from_pickle(
        f"{os.getcwd()}/tmp/temp_selection_genes_mcmc.pickle"
    )

    selected_opt_result = dict_results_gene_selection["selected_opt_result"]
    selected_sampler_MCMC = dict_results_gene_selection["selected_sampler_MCMC"]
    selected_logprob_MCMC = dict_results_gene_selection["selected_logprob_MCMC"]
    params_bs = dict_results_gene_selection["params_bs"]
    times_inferred = dict_results_gene_selection["times_inferred"]
    number_of_knots_explored = dict_results_gene_selection["number_of_knots_to_explore"]
    atributes = plots_utils.get_result_atributes(selected_opt_result)
    depth_genes = selected_opt_result.problem.depth_genes
    freqs_model = selected_opt_result.problem.freqs_model
    interpolator = selected_opt_result.problem.interpolator
    number_free_freqs = atributes["number_free_freqs"]
    frequencies = atributes["frequencies"]
    depth = data[0]
    proxy = data[1]
    frequencies_parameters = parameters_analysis["data_model_parameters"]["frequency_values"]

    samples_combine = selected_sampler_MCMC.reshape(-1, selected_sampler_MCMC.shape[-1])
    logprob_combine = selected_logprob_MCMC.reshape(-1, selected_logprob_MCMC.shape[-1])
    index_MAP_1d = np.argmax(logprob_combine, axis=None)

    j = index_MAP_1d
    free_freqs = samples_combine[j, :number_free_freqs]

    invSR_interpolate = selected_opt_result.problem.interpolator(
        [depth_genes, samples_combine[j, number_free_freqs:]], depth
    )
    times_inferred_MAP = sp.integrate.cumulative_trapezoid(
        invSR_interpolate, depth, initial=0
    )

    phi_MAP = get_phase(
        data[1], frequencies, times_inferred_MAP
    )  # !!! Note here, need to use times_inferred_MAP

    phi_1000 = np.zeros((params_bs.shape[0], len(frequencies)))
    beta_sin_1000 = np.zeros((params_bs.shape[0], len(frequencies)))
    beta_cos_1000 = np.zeros((params_bs.shape[0], len(frequencies)))
    times_inferred_1000 = np.zeros((params_bs.shape[0], len(depth)))
    fs_1000 = np.zeros((params_bs.shape[0], len(frequencies)))

    for i in range(params_bs.shape[0]):
        free_freqs = params_bs[i, :number_free_freqs]
        fs = freqs_model(free_freqs)
        fs_1000[i] = fs

        invSR_interpolate = interpolator(
            [depth_genes, params_bs[i, number_free_freqs:]], depth
        )
        times_inferred = sp.integrate.cumulative_trapezoid(
            invSR_interpolate, depth, initial=0
        )
        times_inferred_1000[i] = times_inferred

        phi, beta_sin, beta_cos = get_phase_beta(proxy, fs, times_inferred)
        phi_1000[i], beta_sin_1000[i], beta_cos_1000[i] = phi, beta_sin, beta_cos
    
    name_freqs = []
    if frequencies_parameters["use_precession"]:
        name_freqs.extend([f"k+g$_{i+1}$" for i in range(5)])
    if frequencies_parameters["use_eccentricity"]:
        name_freqs.extend([
            f"g$_{2}$-g$_{5}$",
            f"g$_{4}$-g$_{5}$",
            f"g$_{4}$-g$_{2}$",
            f"g$_{3}$-g$_{5}$",
            f"g$_{3}$-g$_{2}$",
        ])
    if frequencies_parameters["use_tilt"]:
        name_freqs.extend([f"k+s$_{i}$" for i in [3, 4, 6]])

    # Get the phases of 1000 random solutions, with phase-axis adjusted ====
    phi_1000_new = np.zeros(phi_1000.shape)
    phase_median_mean_std = np.zeros([len(frequencies), 3])

    fig, axs = plt.subplots(
        int(len(frequencies) / 5) + 1,
        5,
        constrained_layout=True,
        figsize=(12, 3.5 * int(len(frequencies) / 5) + 3),
    )
    trim_axs(axs, len(frequencies))
    last_subplot_first_row = min(len(frequencies) - 1, 4)
    for i in range(len(frequencies)):
        ax = axs[int(i / 5), np.mod(i, 5)]
        ax.set_xlabel(name_freqs[i], fontsize=12)

        n_bin = 25
        y_hist, x_hist0, temp = ax.hist(
            phi_1000[:, i].flatten(), density=True, bins=n_bin, alpha=0.7, color="gray"
        )
        x_hist = x_hist0[:-1] + np.pi * 2 / n_bin

        ax.set_xticks(np.arange(-np.pi - 0.01, np.pi + 0.01, np.pi / 2))
        labels = [r"-$\pi$", r"-$\pi/2$", "$0$", r"$\pi/2$", r"$\pi$"]
        ax.set_xticklabels(labels)
        ax.axvline(phi_MAP[i], label="MAP", color="k", zorder=1)

        # To derive the median, with phase-axis adjusted ====
        phi_samples = phi_1000[:, i]
        phi_samples_new = get_phase_median(phi_samples, x_hist, y_hist)
        phi_1000_new[:, i] = phi_samples_new

        #phase_median_mean_std[i][0] = np.median(phi_1000_new[:, i])
        phase_median_mean_std[i][0]=check_median_2pi(np.median(phi_1000_new[:,i]))
        phase_median_mean_std[i][1]=np.mean(phi_1000_new[:,i])
        phase_median_mean_std[i][2]=np.std(phi_1000_new[:,i])
        
        ax.set_ylim(min(y_hist), max(y_hist) * 1.4)
        ax.axvline(
            phase_median_mean_std[i][0], color="k", linestyle="--", label="median"
        )
        # Add legend only in the last subplot of the first row
        if i == last_subplot_first_row:
            ax.legend(loc="upper right", fontsize=10, frameon=True)

    if not output_path or output_path == "default" or output_path == "":
        output_path_logs = parameters_analysis["output_folder"]
        os.makedirs(output_path_logs, exist_ok=True)
    else:
        output_path_logs = output_path
    file_path = f"{output_path_logs}/logs_plots.txt"
    plots_utils.print_and_save(
        message="###########################################", file_path=file_path
    )
    plots_utils.print_and_save(message="", file_path=file_path)
    plots_utils.print_and_save(
        message="PHASE PLOT", print_time=True, file_path=file_path
    )
    plots_utils.print_and_save(message="", file_path=file_path)
    message = "MAP   median   mean   std"
    plots_utils.print_and_save(message=message, print_time=False, file_path=file_path)
    for i in range(len(fs)):
        a = phi_MAP[i] / np.pi
        b, c, d = phase_median_mean_std[i] / np.pi
        message = f"{a:>6.3f}   {b:>6.3f}   {c:>6.3f}   {d:>6.3f} pi"
        plots_utils.print_and_save(
            message=message, print_time=False, file_path=file_path
        )
    plots_utils.print_and_save(message="", file_path=file_path)
    plots_utils.print_and_save(
        message="###########################################", file_path=file_path
    )
    plots_utils.print_and_save(message="", file_path=file_path)

    if types_of_saving and types_of_saving != []:
        if not output_path or output_path == "default" or output_path == "":
            output_path = parameters_analysis["output_folder"]
            output_path = f"{output_path}/figures"
            os.makedirs(output_path, exist_ok=True)
        figname = f"Phase_Plot_{number_of_knots_explored}_Knots"
        global supported_formats
        for fmt in types_of_saving:
            if fmt.lower() in supported_formats:
                if fmt.lower() == "tab":
                    logger.warning("tab format is not accepted for this plot.")
                else:
                    fig.savefig(
                        f"{output_path}/{figname}.{fmt.lower()}", format=fmt.lower()
                    )
                    logger.info(f"Saved {output_path}/{figname}.{fmt.lower()}")
            else:
                logger.warning(f"The format {fmt} is not accepted.")
    # plt.show()
