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
import numpy as np
import pandas as pd
import scipy as sp  # type: ignore
import logging

from typing import Optional
from matplotlib import pyplot as plt  # type: ignore
from matplotlib.ticker import FuncFormatter  # type: ignore
from scipy.interpolate import interp1d  # type: ignore
from statsmodels.tsa.ar_model import ar_select_order  # type: ignore
from sklearn.metrics import r2_score  # type:ignore
import statsmodels.api as sm  # type: ignore
from scipy.stats import pearsonr  # type: ignore
from sklearn.linear_model import LinearRegression  # type:ignore

from astrogeofit.utils import shared_functions, setup_logger
from astrogeofit.main_routines import data_manager
from astrogeofit.plots.utils import shared_utils_plots as plot_utils
from astrogeofit.notebooks_code import notebooks_code

logger = logging.getLogger("ToolLogger")
# Check if the logger already has handlers
if not logger.hasHandlers():
    logger = setup_logger.setup_logger(False)

supported_formats = ["jpeg", "pdf", "eps", "tab"]

def  plt_fit_summary_of_the_results(
    configuration_file_path:str,
    folder_with_the_results:Optional[str],
    file_with_the_results:Optional[str],
    number_of_knots_to_explore: Optional[int],
    types_of_saving: list[str],
    output_path: Optional[str],
    x_axis_limits: Optional[list[float]] = None,
    invSR_plot_y_limits: Optional[list[float]] = None,
    age_plot_y_limits: Optional[list[float]] = None,
    fitted_data_plot_y_limits: Optional[list[float]] = None,
    correlation_plot_y_limits: Optional[list[float]] = None,
):
    """
    Generates and saves a summary figure of the AstroGeoFit fitting results.

    This function creates a 4-panel plot displaying:
        1. Inverse sedimentation rate (invSR) versus depth,
        2. Time-depth (age-depth) model,
        3. Fitted signal predictions compared to the original dataset,
        4. Rolling correlation between predicted and original signals.

    The function loads and processes results from a given configuration and result files,
    supports plotting solutions for a specific number of knots (genes), and optionally adds
    a nominal model as reference. Users can save the summary figure in multiple formats and
    export the data used in each subplot.

    Parameters:
        - **configuration_file_path** (str): Path to the configuration `.yml` file used in the analysis.
        - **folder_with_the_results** (Optional[str]): Folder containing the output results. If None, uses default from config.
        - **file_with_the_results** (Optional[str]): File name of the serialized results. If None, uses default from config.
        - **number_of_knots_to_explore** (Optional[int]): Number of knots (genes) to visualize from the list used during optimization.
            If None, the last available entry is used.
        - **types_of_saving** (list[str]): List of file formats to save outputs (e.g., ["png", "pdf", "tab"]).
        - **output_path** (Optional[str]): Path to the output directory for saving figures and `.tab` files.
            If None, uses the default output folder from the configuration.
        - **ws** (int): Window size for the rolling correlation calculation. Default is 200.
        - **x_axis_limits** (Optional[list[float]]): X-axis limits for all subplots (in depth units).
        - **invSR_plot_y_limits** (Optional[list[float]]): Y-axis limits for the inverse sedimentation rate subplot.
        - **age_plot_y_limits** (Optional[list[float]]): Y-axis limits for the age-depth (time-depth) subplot.
        - **fitted_data_plot_y_limits** (Optional[list[float]]): Y-axis limits for the fitted signal subplot.
        - **correlation_plot_y_limits** (Optional[list[float]]): Y-axis limits for the correlation subplot.

    Returns:
        None

    Saves:
        - A summary figure in the specified formats.
        - Tab-separated `.tab` files for each subplot (if "tab" is listed in types_of_saving), including:
            * invSR values per depth,
            * reconstructed age-depth values,
            * predicted signal values,
            * rolling correlation coefficients.

    Raises:
        ValueError: If the requested number_of_knots_to_explore is not found in the optimization results.
        Logs errors for unsupported save formats.
    """
    
    ################### LIBRARY MODE SECTION ####################
    results_fit, _, _, parameters_analysis, data_set_values, age_depth_model_values, genetic_algorithm_parameters, _, _, _ = notebooks_code.data_obtain_analysis_file_and_result_files(False, configuration_file_path, folder_with_the_results, file_with_the_results)
    
    _, _, _, func_invSR_nominal, func_time_nominal, add_nominal = notebooks_code.data_obtain_age_depth_values(age_depth_model_values)

    data, _ = data_manager.data_read_dataset_from_file(data_set_values)
    
    tool_plt_fit_summary_of_the_results(
        genetic_algorithm_parameters=genetic_algorithm_parameters,
        results_fit=results_fit,
        data=data,
        func_invSR_nominal = func_invSR_nominal,
        func_time_nominal=func_time_nominal,
        parameters_analysis=parameters_analysis,
        number_of_genes_to_explore=number_of_knots_to_explore,
        x_axis_limits=x_axis_limits,
        invSR_plot_y_limits=invSR_plot_y_limits,
        age_plot_y_limits=age_plot_y_limits,
        fitted_data_plot_y_limits=fitted_data_plot_y_limits,
        correlation_plot_y_limits=correlation_plot_y_limits,
        types_of_saving=types_of_saving,
        output_path=output_path,
        add_nominal=add_nominal
    )

def tool_plt_fit_summary_of_the_results(
    genetic_algorithm_parameters: dict,
    results_fit:dict,
    data: list,
    func_invSR_nominal,
    func_time_nominal,
    parameters_analysis: dict,
    number_of_genes_to_explore: Optional[int],
    x_axis_limits: Optional[list[float]] = None,
    invSR_plot_y_limits: Optional[list[float]] = None,
    age_plot_y_limits: Optional[list[float]] = None,
    fitted_data_plot_y_limits: Optional[list[float]] = None,
    correlation_plot_y_limits: Optional[list[float]] = None,
    types_of_saving: list[str] = [],
    output_path: Optional[str] = [],
    add_nominal: bool = False
):
    
    optimization_results = results_fit["optimization_results"]
    invSR_to_predx= results_fit["invSR_to_predx"]
    list_number_genes=genetic_algorithm_parameters["list_number_genes"]
    ws = 200
    
    fig, axs = plt.subplots(
        4, 1, figsize=(10, 10), constrained_layout=True, sharex=True, dpi=100
    )
    cm = plt.get_cmap("rainbow")

    # Handle gene selection
    if number_of_genes_to_explore:
        try:
            index_result_to_explore = list_number_genes.index(
                number_of_genes_to_explore
            )
        except ValueError:
            logger.error(
                f"The selected number of genes {number_of_genes_to_explore} does not exist in the genes used for fitting."
            )
            return
    else:
        index_result_to_explore = -1

    tab_data = {"invSR": [], "age_depth": [], "fitted_data": [], "correlation": []}

    for i, optimization_result in enumerate(optimization_results):
        atributes_result = plot_utils.get_result_atributes(
            optimization_result[index_result_to_explore]
        )

        res_prediction = invSR_to_predx(
            [atributes_result["depth_genes"], atributes_result["best_inverse_SR"]],
            fs=atributes_result["frequencies"],
        )
        time = res_prediction["time"]
        y_pred = res_prediction["y_pred"]
        if len(optimization_result) > 1:
            color = cm(i / (len(optimization_results) - 1))
        else:
            color = cm(i / (len(optimization_results)))
        # Plot inverse sedimentation rate (invSR)
        axs[0].plot(
            atributes_result["depth_genes"],
            atributes_result["best_inverse_SR"],
            ".",
            color=color,
            markersize=5,
        )
        axs[0].plot(
            data[0],
            atributes_result["interpolator"](
                [atributes_result["depth_genes"], atributes_result["best_inverse_SR"]],
                data[0],
            ),
            "-",
            color=color,
        )
        tab_data["invSR"].append(
            (atributes_result["depth_genes"], atributes_result["best_inverse_SR"])
        )

        # Plot age-depth model
        axs[1].plot(data[0], time, color=color, alpha=0.5)
        tab_data["age_depth"].append((data[0], time))

        # Plot fitted data
        axs[2].plot(data[0], y_pred, color="k")
        tab_data["fitted_data"].append((data[0], y_pred))

        # Plot correlation coefficient
        axs[3].plot(
            data[0][:] - (data[0][ws] - data[0][0]) / 2,
            pd.Series(data[1]).rolling(ws).corr(pd.Series(y_pred)),
            color=color,
        )
        tab_data["correlation"].append(
            (
                data[0][:] - (data[0][ws] - data[0][0]) / 2,
                pd.Series(data[1]).rolling(ws).corr(pd.Series(y_pred)),
            )
        )

    # Add nominal plots
    if add_nominal:
        axs[0].plot(
            data[0],
            func_invSR_nominal(data[0]),
            color="black",
            zorder=1000,
            label="Nominal",
        )
        axs[1].plot(
            data[0],
            func_time_nominal(data[0]) - func_time_nominal(data[0])[0],
            color="black",
            zorder=1000,
            lw=2,
        )
    if invSR_plot_y_limits:
        axs[0].set_ylim(*invSR_plot_y_limits)

    if age_plot_y_limits:
        axs[1].set_ylim(*age_plot_y_limits)

    axs[2].plot(data[0], data[1], zorder=-10, label="Original Data")
    axs[2].legend()
    if fitted_data_plot_y_limits:
        axs[2].set_ylim(*fitted_data_plot_y_limits)

    axs[3].set_ylabel("Correlation coefficient")
    if correlation_plot_y_limits:
        axs[3].set_ylim(*correlation_plot_y_limits)

    if x_axis_limits:
        axs[3].set_xlim(*x_axis_limits)
    else:
        axs[3].set_xlim(data[0][0], data[0][-1])

    fig.suptitle(f"Summary of the Results {number_of_genes_to_explore} Knots")
    axs[0].set_ylabel("Inverse Sedimentation Rate R [Myr/m]")
    axs[0].legend()
    axs[0].set_title("Inverse Sedimentation Rate Solutions in Function of Depth")
    axs[1].set_ylabel("Age (Ma)")
    axs[1].set_title("Reconstructed Functions Time-Depth")
    axs[2].set_ylabel("Model Fitted to Original Data")
    axs[2].set_title("Reconstructed Solutions Over Data")
    axs[2].legend()
    axs[-1].set_title("Correlation Coefficient Per Solution Over Depth")
    axs[-1].set_xlabel("Depth (m)")

    # Save plots and data
    if types_of_saving and types_of_saving != []:
        if not output_path or output_path == "default" or output_path == "":
            output_path = parameters_analysis["output_folder"]
            output_path = f"{output_path}/figures"
            os.makedirs(output_path, exist_ok=True)
        figname = f"Summary_of_the_Results_{number_of_genes_to_explore}_Knots"
        global supported_formats
        for fmt in types_of_saving:
            fmt_lower = fmt.lower()
            if fmt_lower not in supported_formats:
                logger.error(
                    f"The file type {fmt} is not supported. The figure will not be saved in this format."
                )
            if fmt_lower == "tab":
                for key, values in tab_data.items():
                    tab_filename = f"{output_path}/{figname}_{key}.tab"
                    with open(tab_filename, "w") as file:
                        for value_set in values:
                            file.write(
                                "\n".join(
                                    "\t".join(map(str, pair))
                                    for pair in zip(*value_set)
                                )
                                + "\n"
                            )
                    logger.info(f"Saved {tab_filename}")
            else:
                fig.savefig(
                    f"{output_path}/{figname}.{fmt_lower}", format=fmt_lower, dpi=300
                )
                logger.info(f"Saved {output_path}/{figname}.{fmt.lower()}")

def plt_fit_metric_per_number_of_knots(
    configuration_file_path:str,
    folder_with_the_results:Optional[str],
    file_with_the_results:Optional[str],
    types_of_saving: list,
    output_path: Optional[str],
):
    """
    Plots the evolution of a performance metric across generations and knot configurations.

    This function visualizes the chosen metric (e.g., r²) over the number of generations
    for different solutions obtained during the optimization process. It also adds vertical
    lines and labels corresponding to different knot (gene) configurations used during fitting.

    Parameters:
        - **configuration_file_path** (str): Path to the configuration `.yml` file used for analysis.
        - **folder_with_the_results** (Optional[str]): Folder containing the output results from the optimization process.
        - **file_with_the_results** (Optional[str]): File containing the results to load from.
        - **types_of_saving** (list[str]): List of file formats to save the plot in (e.g., ["png", "pdf"]).
        - **output_path** (Optional[str]): Directory where the output figure will be saved. If not provided,
            the default output folder from the configuration will be used.

    Returns:
        None

    Notes:
        - The figure includes vertical lines to indicate transitions between different knot values.
        - If "tab" is specified in `types_of_saving`, a warning is issued as tabular saving is not supported for this plot.
        - Supported formats are controlled by the global `supported_formats`.

    Saves:
        A figure showing the metric evolution over generations, with optional saving in specified formats.

    Raises:
        Logs warnings if unsupported formats are requested or if saving as "tab" is attempted.
    """
    
    ################### LIBRARY MODE SECTION ####################
    results_fit, _, _, parameters_analysis, _, _, genetic_algorithm_parameters, _, _, _ = notebooks_code.data_obtain_analysis_file_and_result_files(False, configuration_file_path, folder_with_the_results, file_with_the_results)
    
    tool_plt_fit_metric_per_number_of_knots(
        genetic_algorithm_parameters=genetic_algorithm_parameters,
        results_fit=results_fit,
        parameters_analysis=parameters_analysis,
        types_of_saving=types_of_saving,
        output_path=output_path
    )

def tool_plt_fit_metric_per_number_of_knots(
    genetic_algorithm_parameters:dict,
    results_fit:dict,
    parameters_analysis: dict,
    types_of_saving: list,
    output_path: Optional[str],
):
    number_algorithm_solutions = genetic_algorithm_parameters["number_algorithm_solutions"]
    population_size = genetic_algorithm_parameters["population_size"]
    number_generations = genetic_algorithm_parameters["number_generations"]
    optimization_results = results_fit["optimization_results"]
    metric_sorted_resuts = results_fit["metric_sorted_resuts"]
    list_number_genes = genetic_algorithm_parameters["list_number_genes"]
    metric_type = genetic_algorithm_parameters["metric_type"]
    fig, ax = plt.subplots(1, 1, figsize=(8, 3), dpi=300)
    plt.subplots_adjust(
        left=0.1, right=0.95, bottom=0.15, top=0.88, wspace=0.1, hspace=0.7
    )
    if metric_type == "r2":
        metric_type = "r$^2$"
    cm = plt.get_cmap("rainbow")
    for i, _ in enumerate(optimization_results):
        if len(optimization_results) > 1:
            color = cm(i / (len(optimization_results) - 1))
        else:
            color = cm(i / (len(optimization_results)))
        ax.plot(
            metric_sorted_resuts[i].flatten(),
            color=color,
            zorder=i,
        )
    ax.set_ylabel(metric_type, fontsize=10)
    ax.set_xlabel("Number of generation")

    for i in range(len(list_number_genes)):  # ---- YW
        ax.axvline(number_generations * (i + 1))
        x_text = (0.2 + i) * number_generations
        y_text = np.percentile(ax.get_ylim(), 10)
        ax.text(x_text, y_text, f"{list_number_genes[i]}", fontsize=7, zorder=1000) #----YW

    plt.title(
        f"{number_algorithm_solutions} number of solutions (Population Size={population_size}, Number of Generations={number_generations})"
    )
    # Save plots
    if types_of_saving and types_of_saving != []:
        if not output_path or output_path == "default" or output_path == "":
            output_path = parameters_analysis["output_folder"]
            output_path = f"{output_path}/figures"
            os.makedirs(output_path, exist_ok=True)
        figname = "Metric_Per_Number_of_Knots"

        for fmt in types_of_saving:
            fmt_lower = fmt.lower()
            if fmt_lower not in supported_formats:
                logger.warning(
                    f"The file type {fmt} is not supported. The figure will not be saved in this format."
                )
            elif fmt_lower == "tab":
                # NECESSARY TO IMPLEMENT TAB?
                logger.warning("This plot can not be saved as tab file.")
            else:
                fig.savefig(
                    f"{output_path}/{figname}.{fmt_lower}", format=fmt_lower, dpi=300
                )
                logger.info(f"Saved {output_path}/{figname}.{fmt.lower()}")
    # plt.show()

def plt_fit_aics_bics_per_number_knots(
    configuration_file_path:str,
    folder_with_the_results:Optional[str],
    file_with_the_results:Optional[str],
    types_of_saving: list[str],
    output_path: Optional[str] = None,
):
    """
    Plot AIC and BIC values across different knot configurations and algorithm solutions.

    This function computes and visualizes the Akaike Information Criterion (AIC) and Bayesian 
    Information Criterion (BIC) for each solution obtained during the optimization process. 
    For each number of knots (genes), the function fits an autoregressive model to the residuals 
    between the fitted and original data, and computes the corresponding AIC and BIC values. 
    The final figure contains two subplots, one for AIC and one for BIC, with color-coded lines 
    representing each solution and a black line for the minimum value per knot configuration.

    Parameters:
        - **configuration_file_path** (str): Path to the `.yml` configuration file used during the analysis.
        - **folder_with_the_results** (Optional[str]): Folder containing the results of the optimization.
        - **file_with_the_results** (Optional[str]): File from which the results should be loaded.
        - **types_of_saving** (list[str]): List of output file formats to save the figure in (e.g., ["png", "pdf"]).
        - **output_path** (Optional[str]): Directory to save the figure. If None, the path from the configuration
            file is used under the `figures` subdirectory.

    Returns:
        None

    Notes:
        - The function automatically computes the residuals between predicted and original data.
        - An autoregressive model is fit to the residuals for AIC/BIC calculation using `ar_select_order`.
        - Each subplot contains one line per algorithm solution, color-coded, and a black line for the minimum values.
        - If "tab" is included in `types_of_saving`, a warning is logged, as saving this plot in `.tab` format is not supported.
        - Only supported image formats (defined in `supported_formats`) will be saved. Unsupported ones will be ignored with a warning.

    Saves:
        - A figure with AIC and BIC plots in the specified formats in the output directory.

    Raises:
        - Logs warnings for unsupported formats and if tabular export is requested.
    """
    
    ################### LIBRARY MODE SECTION ####################
    results_fit, _, _, parameters_analysis, data_set_values, _, genetic_algorithm_parameters, _, _, _ = notebooks_code.data_obtain_analysis_file_and_result_files(False, configuration_file_path, folder_with_the_results, file_with_the_results)
    data, _ = data_manager.data_read_dataset_from_file(data_set_values)
    
    tool_plt_fit_aics_bics_per_number_knots(
        genetic_algorithm_parameters=genetic_algorithm_parameters,
        results_fit=results_fit,
        data=data,
        parameters_analysis=parameters_analysis,
        types_of_saving=types_of_saving,
        output_path=output_path
    )

def tool_plt_fit_aics_bics_per_number_knots(
    genetic_algorithm_parameters:dict,
    results_fit:dict,
    data: list,
    parameters_analysis: dict,
    types_of_saving: list[str],
    output_path: Optional[str] = None,
):
    """
    Plot AIC and BIC values for different algorithm solutions and numbers of genes.

    Parameters:
        number_algorithm_solutions (int): Total number of algorithm solutions.
        list_number_genes (List): List of numbers of genes.
        optimization_results (List): Optimization results for each solution.
        invSR_to_predx: Function to compute predictions from inverse sedimentation rates.
        data (List): Data used for computing residuals.
        save_plot (bool): Whether to save the plot.
        types_of_saving (List[str]): File formats for saving the plot (e.g., ["jpg", "eps"]).
        output_path (Optional[str]): Path to save the plots.

    Returns:
        None
    """
    number_algorithm_solutions = genetic_algorithm_parameters["number_algorithm_solutions"]
    list_number_genes = genetic_algorithm_parameters["list_number_genes"]
    optimization_results = results_fit["optimization_results"]
    invSR_to_predx= results_fit["invSR_to_predx"]
    
    cm = plt.get_cmap("rainbow")
    aics = np.zeros((number_algorithm_solutions, len(list_number_genes)))
    bics = np.zeros((number_algorithm_solutions, len(list_number_genes)))

    # Compute AIC and BIC values
    for i, optimization_result in enumerate(optimization_results):
        for j, result in enumerate(optimization_result):
            atributes = plot_utils.get_result_atributes(result=result)
            res_prediction = invSR_to_predx(
                [atributes["depth_genes"], atributes["best_inverse_SR"]],
                fs=atributes["frequencies"],
            )
            y_pred = res_prediction["y_pred"]
            residual = data[1] - y_pred
            sel = ar_select_order(residual, 13)
            mod = sel.model
            res_ar = mod.fit()
            aics[i, j] = res_ar.aic + len(atributes["depth_genes"]) * 2
            bics[i, j] = res_ar.bic + len(atributes["depth_genes"]) * np.log(
                len(data[1])
            )

    fig, axs = plt.subplots(1, 2, figsize=(6.4, 2.4), constrained_layout=True, dpi=300)
    fig.suptitle("AIC and BIC Metrics for the Fitting Results")
    # Plot AIC values
    for i in range(number_algorithm_solutions):
        if number_algorithm_solutions > 1:
            axs[0].plot(
                list_number_genes,
                aics[i],
                color=cm(i / (number_algorithm_solutions - 1)),
            )
        else:
            axs[0].plot(
                list_number_genes,
                aics[i],
                color=cm(i / (number_algorithm_solutions)),
            )
    axs[0].plot(list_number_genes, aics.min(axis=0), color="k")
    axs[0].set_ylabel("AIC")
    axs[0].set_xlabel("Number of Knots")

    # Plot BIC values
    for i in range(number_algorithm_solutions):
        axs[1].plot(
            list_number_genes,
            bics[i],
            color=cm(i / (number_algorithm_solutions - 1)),
        )
    axs[1].plot(list_number_genes, bics.min(axis=0), color="k")
    axs[1].set_ylabel("BIC")
    axs[1].set_xlabel("Number of Knots")

    # Save plots
    if types_of_saving and types_of_saving != []:
        if not output_path or output_path == "default" or output_path == "":
            output_path = parameters_analysis["output_folder"]
            output_path = f"{output_path}/figures"
            os.makedirs(output_path, exist_ok=True)
        figname = "AIC_and_BIC"
        for fmt in types_of_saving:
            fmt_lower = fmt.lower()
            if fmt_lower not in supported_formats:
                logger.warning(
                    f"The file type {fmt} is not supported. The figure will not be saved in this format."
                )
            elif fmt_lower == "tab":
                # NECESSARY TO IMPLEMENT TAB?
                logger.warning("This plot can not be saved as tab file.")
            else:
                plt.savefig(
                    f"{output_path}/{figname}.{fmt_lower}", format=fmt_lower, dpi=300
                )
                logger.info(f"Saved {output_path}/{figname}.{fmt.lower()}")

    # plt.show()

def plt_fit_aics_r2_per_number_knots(
    configuration_file_path:str,
    folder_with_the_results:Optional[str],
    file_with_the_results:Optional[str],
    types_of_saving: list,
    output_path: str,
):
    """
    Plot AIC and r² values for a single optimization solution across different knot configurations.

    This function generates a dual-axis plot comparing the Akaike Information Criterion (AIC) 
    and the coefficient of determination (r²) across varying numbers of knots (genes) for the 
    final optimization result. The AIC is calculated based on fitting an autoregressive model 
    to the residuals between the predicted and original data. The r² score evaluates the 
    goodness of fit of the model predictions.

    Parameters:
        - **configuration_file_path** (str): Path to the `.yml` configuration file used for the analysis.
        - **folder_with_the_results** (Optional[str]): Folder containing the optimization results.
        - **file_with_the_results** (Optional[str]): File from which results should be loaded.
        - **types_of_saving** (list): List of file formats to save the figure in (e.g., ["png", "pdf"]).
        - **output_path** (str): Directory to save the plot. If not provided or set to "default", it falls
            back to the path defined in the configuration under `output_folder/figures`.

    Returns:
        None

    Notes:
        - Only the last (final) solution in the optimization results is used for this plot.
        - AIC is adjusted by the number of model parameters (twice the number of knots).
        - r² is computed between the observed and predicted data series.
        - The left y-axis shows the AIC values (shifted to start at 0 for relative comparison), 
        and the right y-axis shows the r² score.
        - A warning is logged for unsupported formats or if "tab" is requested, which is not supported.

    Saves:
        - A figure comparing AIC and r² across knot counts, in the specified formats and output directory.

    Raises:
        - Logs warnings for unsupported formats or attempts to save as a "tab" file.
    """
    results_fit, _, _, parameters_analysis, data_set_values, _, genetic_algorithm_parameters, _, _, _ = notebooks_code.data_obtain_analysis_file_and_result_files(False, configuration_file_path, folder_with_the_results, file_with_the_results)
    data, _ = data_manager.data_read_dataset_from_file(data_set_values)
    tool_plt_fit_aics_r2_per_number_knots(
        genetic_algorithm_parameters=genetic_algorithm_parameters,
        results_fit=results_fit,
        data=data,
        parameters_analysis=parameters_analysis,
        types_of_saving=types_of_saving,
        output_path=output_path
    )

def tool_plt_fit_aics_r2_per_number_knots(
    genetic_algorithm_parameters:dict,
    results_fit:dict,
    data: list,
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
    list_number_genes=genetic_algorithm_parameters["list_number_genes"]
    invSR_to_predx= results_fit["invSR_to_predx"]

    aics = np.zeros(len(list_number_genes))
    r2s = np.zeros(len(list_number_genes))
    for k, result in enumerate(single_opt_result):
        atributes = plot_utils.get_result_atributes(result)
        res_prediction = invSR_to_predx(
            [atributes["depth_genes"], atributes["best_inverse_SR"]],
            fs=atributes["frequencies"],
        )
        y_pred = res_prediction["y_pred"]
        residual = data[1] - y_pred
        sel = ar_select_order(residual, 13)
        mod = sel.model
        res_ar = mod.fit()
        aics[k] = res_ar.aic + len(atributes["depth_genes"]) * 2
        r2s[k] = r2_score(data[1], y_pred)

    fig, ax = plt.subplots(1, 1)
    fig.suptitle("AICs and Metric Results per Each Number of Knots.")
    ax2 = ax.twinx()
    ax2.plot(list_number_genes, r2s, ".-", color="r")
    ax2.set_ylabel("r$^2$ ", color="r", rotation=0)
    ax.plot(list_number_genes, aics - aics[0], ".-", color="k")
    ax.set_xlabel("Number of Knots")
    ax.set_ylabel("AIC")
    color_y_axis(ax, "k")
    color_y_axis(ax2, "r")
    if types_of_saving and types_of_saving != []:
        if not output_path or output_path == "default" or output_path == "":
            output_path = parameters_analysis["output_folder"]
            output_path = f"{output_path}/figures"
            os.makedirs(output_path, exist_ok=True)
        figname = "AIC_And_r2"
        for fmt in types_of_saving:
            fmt_lower = fmt.lower()
            if fmt_lower not in supported_formats:
                logger.warning(
                    f"The file type {fmt} is not supported. The figure will not be saved in this format."
                )
            elif fmt_lower == "tab":
                # NECESSARY TO IMPLEMENT TAB?
                logger.warning("This plot can not be saved as tab file.")
            else:
                plot_filename = f"{output_path}/{figname}.{fmt_lower}"
                plt.savefig(plot_filename, format=fmt_lower, dpi=300)
                logger.info(f"Saved {plot_filename}")
    # plt.show()

def plt_fit_sedimentation_rate_per_depth(
    configuration_file_path:str,
    folder_with_the_results:Optional[str],
    file_with_the_results:Optional[str],
    y_axis_range: list,
    types_of_saving: list,
    output_path: Optional[str],
):
    """
    Plot inverse sedimentation rate (SR) as a function of depth for multiple knot configurations.

    This function generates subplots showing the inverse sedimentation rate profiles for each 
    number of knots (genes) used during the optimization. The plot helps to visualize how the 
    sedimentation model changes depending on the resolution (number of knots) and across multiple 
    optimization results. Optionally, a nominal sedimentation rate profile is overlaid for comparison.

    Parameters:
        - **configuration_file_path** (str): Path to the `.yml` configuration file used for the analysis.
        - **folder_with_the_results** (Optional[str]): Folder containing the result files.
        - **file_with_the_results** (Optional[str]): Specific file name with the results, if applicable.
        - **y_axis_range** (list): Optional range for the y-axis (e.g., [min, max]) to constrain plots vertically.
        - **types_of_saving** (list): List of file formats to save the figures in. Supported formats include
            image types (e.g., "png", "pdf") and "tab" for saving tabular data.
        - **output_path** (Optional[str]): Directory to save output plots. If set to "default" or left empty, 
            uses the configured output directory from the analysis configuration.

    Returns:
        None

    Saves:
        - A multi-panel figure showing inverse sedimentation rates vs. depth across different knot settings.
        - (If requested) A `.tab` file with depth, inverse sedimentation rate, and corresponding knot count.

    Notes:
        - The number of subplots is determined by the number of knot configurations.
        - For each configuration, the optimized profile is drawn using all available optimization results.
        - The nominal sedimentation rate profile (if available) is added in black for comparison.
        - When saving in `.tab` format, results for all knot configurations are appended in the same file.

    Raises:
        - Logs warnings for unsupported formats or when trying to save in an unsupported format.
    """
    
    ################### LIBRARY MODE SECTION ####################
    results_fit, _, _, parameters_analysis, data_set_values, age_depth_model_values, genetic_algorithm_parameters, _, _, _ = notebooks_code.data_obtain_analysis_file_and_result_files(False, configuration_file_path, folder_with_the_results, file_with_the_results)
    
    _, _, _, func_invSR_nominal, _, add_nominal = notebooks_code.data_obtain_age_depth_values(age_depth_model_values)

    data, _ = data_manager.data_read_dataset_from_file(data_set_values)
    tool_plt_fit_sedimentation_rate_per_depth(
        genetic_algorithm_parameters=genetic_algorithm_parameters,
        results_fit=results_fit,
        data=data,
        func_invSR_nominal=func_invSR_nominal,
        parameters_analysis=parameters_analysis,
        y_axis_range=y_axis_range,
        types_of_saving=types_of_saving,
        output_path=output_path,
        add_nominal=add_nominal
    )

def tool_plt_fit_sedimentation_rate_per_depth(
    genetic_algorithm_parameters:dict,
    results_fit:dict,
    data: list,
    func_invSR_nominal,
    parameters_analysis: dict,
    y_axis_range: list,
    types_of_saving: list,
    output_path: Optional[str],
    add_nominal: bool,
):
    optimization_results = results_fit["optimization_results"]
    list_number_genes=genetic_algorithm_parameters["list_number_genes"]
    
    n_plots = len(list_number_genes)
    n_cols = int(np.ceil(np.sqrt(n_plots)))
    n_rows = int(np.ceil(n_plots / n_cols))

    _, axs = plt.subplots(
        n_rows,
        n_cols,
        figsize=(10, 5),
        constrained_layout=True,
        sharex=True,
        sharey=True,
        dpi=100,
        squeeze=False,
    )
    ax = axs.flatten()
    cm = plt.get_cmap("rainbow")
    # Prepare to save tab-separated data if needed
    if "tab" in types_of_saving:
        if not output_path or output_path == "default" or output_path == "":
            output_path = parameters_analysis["output_folder"]
            output_path = f"{output_path}/figures"
            os.makedirs(output_path, exist_ok=True)
        figname = "Inverse_Sedimentation_Rate_with_Various_Number_of_Knots"
        tab_filename = os.path.join(output_path, figname, ".tab")
        with open(tab_filename, "w") as tab_file:
            tab_file.write("Depth   Inverse_Sedimentation_Rate  Number_of_Knots\n")

    for i in range(len(list_number_genes)):
        if add_nominal:
            ax[i].plot(
                data[0],
                func_invSR_nominal(data[0]),
                color="black",
                zorder=1000,
                label="Nominal",
            )

        for j, optimization_result in enumerate(optimization_results):
            atributes = plot_utils.get_result_atributes(result=optimization_result[i])
            color = cm(j / (len(optimization_results) - 1))
            depth_genes = atributes["depth_genes"]
            inverse_SR = atributes["best_inverse_SR"]
            ax[i].set_title(f"{len(depth_genes)} Knots", fontsize=8)

            # Plot the data
            ax[i].plot(depth_genes, inverse_SR, ".", color=color, markersize=5)
            ax[i].plot(
                data[0],
                atributes["interpolator"]([depth_genes, inverse_SR], data[0]),
                "-",
                color=color,
            )
            if y_axis_range:
                ax[i].set_ylim(y_axis_range[0], y_axis_range[1])
            # If saving data in .tab format, write it to the file
            if "tab" in types_of_saving:
                with open(tab_filename, "a") as tab_file:
                    for d, sr in zip(depth_genes, inverse_SR):
                        tab_file.write(f"{d}    {sr}    {len(depth_genes)}\n")

    # Add labels to the axes
    for i in range(n_rows):
        axs[i, 0].set_ylabel("Inverse SR[Myr/m]")
    for i in range(n_cols):
        axs[-1, i].set_xlabel("Depth (m)")

    # Save the plot in other formats if requested
    if types_of_saving and types_of_saving != []:
        global supported_formats
        if not output_path or output_path == "default" or output_path == "":
            output_path = parameters_analysis["output_folder"]
            output_path = f"{output_path}/figures"
            os.makedirs(output_path, exist_ok=True)
        figname = "Inverse_Sedimentation_Rate_with_Various_Number_of_Knots"
        for fmt in types_of_saving:
            fmt_lower = fmt.lower()
            if fmt_lower not in supported_formats:
                logger.warning(
                    f"The file type {fmt} is not supported. The figure will not be saved in this format."
                )
            elif "tab" in types_of_saving:
                logger.warning(f"Saved {output_path}/{figname}.{fmt.lower()}")
            else:
                plt.savefig(
                    f"{output_path}/{figname}.{fmt_lower}", format=fmt_lower, dpi=300
                )
                logger.info(f"Saved {output_path}/{figname}.{fmt.lower()}")
    # plt.show()

def plt_fit_sedimentation_rate_depth_with_uncertainty(
    configuration_file_path:str,
    folder_with_the_results:Optional[str],
    file_with_the_results:Optional[str],
    number_of_knots_to_explore: Optional[int],
    types_of_saving: list,
    output_path: str,
):
    """
    Plot the inverse sedimentation rate (SR) with uncertainty bands for a selected number of knots.

    This function computes the median inverse sedimentation rate and associated uncertainty intervals 
    (50% and 90% confidence intervals) across multiple optimization results for a given number of knots. 
    It visualizes these estimates as a function of depth, and optionally overlays the nominal sedimentation 
    rate for comparison. The resulting plot and optional tabular data can be saved in user-specified formats.

    Parameters:
        - **configuration_file_path** (str): Path to the `.yml` configuration file used for the analysis.
        - **folder_with_the_results** (Optional[str]): Folder containing result files.
        - **file_with_the_results** (Optional[str]): Name of the results file, if applicable.
        - **number_of_knots_to_explore** (Optional[int]): Number of knots (genes) to visualize the sedimentation 
            rate for. If not specified, the last configuration in the results will be used.
        - **types_of_saving** (list): List of file types to save the figure. Acceptable formats include image 
            formats (e.g., "png", "pdf") and "tab" for tab-separated data export.
        - **output_path** (str): Path to the directory where the outputs will be saved. If set to "default" 
            or left empty, the output path is taken from the configuration.

    Returns:
        None

    Saves:
        - A figure of the inverse sedimentation rate with uncertainty bounds (if requested).
        - A `.tab` file with depth, median SR, 50% and 90% confidence interval bounds (if requested).

    Notes:
        - The 50% confidence interval corresponds to the 25th and 75th percentiles.
        - The 90% confidence interval corresponds to the 5th and 95th percentiles.
        - The nominal inverse SR is plotted in black if available in the age-depth model configuration.

    Raises:
        ValueError: If the specified number of knots is not found in the available configurations.
        Logs warnings for unsupported file formats and records where files are saved.
    """
    
    ################### LIBRARY MODE SECTION ####################
    results_fit, _, _, parameters_analysis, data_set_values, age_depth_model_values, genetic_algorithm_parameters, _, _, _ = notebooks_code.data_obtain_analysis_file_and_result_files(False, configuration_file_path, folder_with_the_results, file_with_the_results)
    
    _, _, _, func_invSR_nominal, _, add_nominal = notebooks_code.data_obtain_age_depth_values(age_depth_model_values)

    data, _ = data_manager.data_read_dataset_from_file(data_set_values)
    tool_plt_fit_sedimentation_rate_depth_with_uncertainty(
        genetic_algorithm_parameters=genetic_algorithm_parameters,
        results_fit=results_fit,
        data=data,
        func_invSR_nominal=func_invSR_nominal,
        parameters_analysis=parameters_analysis,
        number_of_genes_to_explore=number_of_knots_to_explore,
        types_of_saving=types_of_saving,
        output_path=output_path,
        add_nominal=add_nominal
    )
    
def tool_plt_fit_sedimentation_rate_depth_with_uncertainty(
    genetic_algorithm_parameters:dict,
    results_fit:dict,
    data: list,
    func_invSR_nominal,
    parameters_analysis: dict,
    number_of_genes_to_explore: Optional[int],
    types_of_saving: list,
    output_path: str,
    add_nominal: bool,
):
    optimization_results = results_fit["optimization_results"]
    list_number_genes=genetic_algorithm_parameters["list_number_genes"]
    invSR_pop = []
    if number_of_genes_to_explore:
        try:
            index_result_to_explore = list_number_genes.index(
                number_of_genes_to_explore
            )
        except ValueError:
            logger.error(
                f"The selected number of genes {number_of_genes_to_explore} does not exist in the genes used for fitting."
            )
            raise ValueError
    else:
        index_result_to_explore = -1

    for _, optimization_result in enumerate(optimization_results):
        atributes = plot_utils.get_result_atributes(
            optimization_result[index_result_to_explore]
        )
        invSR_interpolate = atributes["interpolator"](
            [atributes["depth_genes"], atributes["best_inverse_SR"]], data[0]
        )
        invSR_pop.append(invSR_interpolate)

    fig, ax = plt.subplots(1, 1, figsize=(8, 3), constrained_layout=True)
    invSR_pop = np.array(invSR_pop)
    plt.fill_between(
        data[0],
        *np.quantile(invSR_pop, [0.25, 0.75], axis=0),
        alpha=0.2,
        label="50% CI",
    )
    plt.fill_between(
        data[0],
        *np.quantile(invSR_pop, [0.05, 0.95], axis=0),
        alpha=0.2,
        zorder=-100,
        label="90% CI",
    )

    plt.plot(data[0], np.median(invSR_pop, axis=0), label="Median")
    if add_nominal:
        plt.plot(
            data[0],
            func_invSR_nominal(data[0]),
            color="black",
            zorder=100,
            label="Nominal",
        )
    plt.legend()

    plt.xlim([data[0][0], data[0][-1]])
    plt.xlabel("Depth (m)")
    plt.ylabel("Inverse Sedimentation Rate (Myr/m)")

    # ====
    if types_of_saving and types_of_saving != []:
        global supported_formats
        if not output_path or output_path == "default" or output_path == "":
            output_path = parameters_analysis["output_folder"]
            output_path = f"{output_path}/figures"
            os.makedirs(output_path, exist_ok=True)
        figname = f"Best_Inverse_Sedimentation_Rate_{number_of_genes_to_explore}_Knots"
        for fmt in types_of_saving:
            fmt_lower = fmt.lower()
            if fmt_lower in supported_formats:
                if fmt_lower == "tab":
                    with open(f"{output_path}/{figname}.tab", "w") as tab_file:
                        tab_file.write(
                            "Depth\tMedian_Inverse_Sedimentation_Rate\t50%_CI_Lower\t50%_CI_Upper\t90%_CI_Lower\t90%_CI_Upper\n"
                        )
                        depth_values = data[0]
                        median_values = np.median(invSR_pop, axis=0)
                        ci_50 = np.quantile(invSR_pop, [0.25, 0.75], axis=0)
                        ci_90 = np.quantile(invSR_pop, [0.05, 0.95], axis=0)
                        for (
                            d,
                            med,
                            ci50_lower,
                            ci50_upper,
                            ci90_lower,
                            ci90_upper,
                        ) in zip(
                            depth_values,
                            median_values,
                            ci_50[0],
                            ci_50[1],
                            ci_90[0],
                            ci_90[1],
                        ):
                            tab_file.write(
                                f"{d}\t{med}\t{ci50_lower}\t{ci50_upper}\t{ci90_lower}\t{ci90_upper}\n"
                            )

                else:
                    fig.savefig(
                        f"{output_path}/{figname}.{fmt_lower}", format=fmt_lower
                    )
                logger.info(f"Saved {output_path}/{figname}.{fmt.lower()}")
            else:
                logger.warning(
                    f"The file type {fmt} is not supported. The figure will not be saved in this format."
                )

def plt_fit_frequency_spectrum(
    configuration_file_path:str,
    folder_with_the_results:Optional[str],
    file_with_the_results:Optional[str],
    x_axis_plot_limits: Optional[list],
    y_axis_plot_limits: Optional[list],
    peridogram_scale_x_axis: str,
    peridogram_scale_y_axis: str,
    number_of_knots_to_explore: Optional[int],
    types_of_saving: list,
    output_path: str,
):
    """
    Plot the power spectrum of the time series inferred from sedimentation rate fitting.

    This function estimates a time series from sedimentation rate fitting results and performs
    a Fast Fourier Transform (FFT) to compute its power spectrum. The resulting spectrum is
    plotted as a function of frequency (cycles/Myr), with an optional top x-axis showing the
    corresponding periods (Myr). The user can customize axis scaling, plot limits, and file 
    formats for saving the figure.

    Parameters:
        - **configuration_file_path** (str): Path to the `.ymal` configuration file used for the fitting.
        - **folder_with_the_results** (Optional[str]): Folder containing result files.
        - **file_with_the_results** (Optional[str]): Name of the result file, if applicable.
        - **x_axis_plot_limits** (Optional[list]): List of two floats [min, max] to limit the x-axis (frequency).
        - **y_axis_plot_limits** (Optional[list]): List of two floats [min, max] to limit the y-axis (amplitude).
        - **peridogram_scale_x_axis** (str): Scale type for the x-axis; "linear" (default) or "log".
        - **peridogram_scale_y_axis** (str): Scale type for the y-axis; "linear" (default) or "log".
        - **number_of_knots_to_explore** (Optional[int]): Number of knots (genes) to select a specific fitting result.
            If not specified, the last configuration will be used.
        - **types_of_saving** (list): List of file types to save the figure (e.g., "png", "pdf").
            The "tab" format is not supported for this plot and will raise a warning.
        - **output_path** (str): Output folder where figures will be saved. If set to "default" or empty,
            the path is obtained from the configuration file.

    Returns:
        None

    Saves:
        - A plot of the power spectrum with FFT-derived amplitudes.
        - Vertical dashed lines showing astronomical target frequencies from the fitting.

    Notes:
        - The y-axis displays scaled FFT amplitudes.
        - The top x-axis shows period values (1/frequency * 1000) in Myr.
        - Frequencies from the best-fit astronomical targets are marked with dashed lines.
        - Unsupported formats and saving issues are logged via `logger`.

    Raises:
        ValueError: If the selected number of knots is not found in the fitting results.
    """
    
    ################### LIBRARY MODE SECTION ####################
    results_fit, _, _, parameters_analysis, data_set_values, age_depth_model_values, genetic_algorithm_parameters, _, _, _ = notebooks_code.data_obtain_analysis_file_and_result_files(False, configuration_file_path, folder_with_the_results, file_with_the_results)

    data, _ = data_manager.data_read_dataset_from_file(data_set_values)
    tool_plt_fit_frequency_spectrum(
        genetic_algorithm_parameters=genetic_algorithm_parameters,
        results_fit=results_fit,
        data=data,
        number_of_knots_to_explore=number_of_knots_to_explore,
        parameters_analysis=parameters_analysis,
        x_limits=x_axis_plot_limits,
        y_limits=y_axis_plot_limits,
        peridogram_scale_x_axis=peridogram_scale_x_axis,
        peridogram_scale_y_axis= peridogram_scale_y_axis,
        types_of_saving=types_of_saving,
        output_path=output_path
    )
    
def tool_plt_fit_frequency_spectrum(
    genetic_algorithm_parameters:dict,
    results_fit:dict,
    data: list,
    number_of_knots_to_explore: Optional[int],
    parameters_analysis: dict,
    x_limits: Optional[list],
    y_limits: Optional[list],
    peridogram_scale_x_axis: str,
    peridogram_scale_y_axis: str, 
    types_of_saving: list,
    output_path: str,
):
    
    optimization_results = results_fit["optimization_results"]
    list_number_of_knots=genetic_algorithm_parameters["list_number_genes"]
    # This cell plot the spectrum of time series ===========================================

    # Get times_inferred from the fitting result ====
    if number_of_knots_to_explore:
        try:
            index_result_to_explore = list_number_of_knots.index(
                number_of_knots_to_explore
            )
        except ValueError:
            logger.error(
                f"The selected number of genes {number_of_knots_to_explore} does not exist in the knots used for fitting."
            )
            raise ValueError
    else:
        index_result_to_explore = -1
    atributes = plot_utils.get_result_atributes(
        optimization_results[-1][index_result_to_explore]
    )
    invSR_interpolate = atributes["interpolator"](
        [atributes["depth_genes"], atributes["best_inverse_SR"]], data[0]
    )
    times_inferred = sp.integrate.cumulative_trapezoid(
        invSR_interpolate, data[0], initial=0
    )

    # Get interpolated data[1] to plot FFT ====
    y_interp = interp1d(times_inferred, data[1])
    times_grid = np.linspace(
        times_inferred.min(), times_inferred.max(), len(times_inferred)
    )
    N = len(times_grid)
    y_equ_spa = y_interp(times_grid)

    # Perform FFT ====
    dt = times_grid[1] - times_grid[0]
    xf = sp.fft.fftfreq(N, dt)[: N // 2]
    yf = sp.fft.fft(y_equ_spa)

    # Plot ====
    fig, ax = plt.subplots(1, 1, figsize=(8, 4), constrained_layout=True)
    ax.plot(xf, 2.0 / N * np.abs(yf[0 : N // 2]))
    yf_magnitude = np.abs(yf[0 : N // 2])
    # Apply the scaling factor 2.0 / N
    scaled_y_values = 2.0 / N * yf_magnitude
    # Find the maximum value (this is the maximum height of the plot)
    max_height = np.max(scaled_y_values) * 1.4
    ax.vlines(
        atributes["frequencies"] / 2 / np.pi,
        0,
        max_height,
        linestyles="--",
        colors="black",
        alpha=0.5,
        zorder=-100,
    )
    if x_limits:
        ax.set_xlim(x_limits)
    if y_limits:
        ax.set_ylim(y_limits)
    if not peridogram_scale_x_axis or peridogram_scale_x_axis == "linear":
        ax.set_xscale("linear")
    elif peridogram_scale_x_axis == "log":
        ax.set_xscale("log")
    if not peridogram_scale_y_axis or peridogram_scale_y_axis == "linear":
        ax.set_xscale("linear")
    elif peridogram_scale_y_axis == "log":
        ax.set_xscale("log")
    ax.set_xlabel("Frequency (cycles/Myr)")
    secax = ax.secondary_xaxis("top", functions=(lambda f: f, lambda p: p))
    formatter = FuncFormatter(lambda x_val, tick_pos: f"{1/x_val*1000:.1f}")
    secax.xaxis.set_major_formatter(formatter)
    secax.set_xlabel("Period (kyr)")
    if types_of_saving and types_of_saving != []:
        if not output_path or output_path == "default" or output_path == "":
            output_path = parameters_analysis["output_folder"]
            output_path = f"{output_path}/figures"
            os.makedirs(output_path, exist_ok=True)
        global supported_formats
        figname = f"Spectrum_Plot_{number_of_knots_to_explore}_Knots"
        for fmt in types_of_saving:
            fmt_lower = fmt.lower()
            if fmt_lower in supported_formats:
                if fmt_lower == "tab":
                    logger.warning("This plot can not be saved as tab file.")
                else:
                    fig.savefig(
                        f"{output_path}/{figname}.{fmt_lower}", format=fmt_lower
                    )
                    logger.info(f"Saved {output_path}/{figname}.{fmt_lower}")
            else:
                logger.warning(
                    f"The file type {fmt} is not supported. The figure will not be saved in this format."
                )

def plt_fit_time_series(
    configuration_file_path:str,
    folder_with_the_results:Optional[str],
    file_with_the_results:Optional[str],
    number_of_knots_to_explore:Optional[int],
    types_of_saving: list,
    output_path: str,
):
    """
    Plot depth and time series of proxy data and model predictions from a genetic algorithm fitting,
    and optionally save the resulting figures and/or data in various formats.

    This function reads the dataset and optimization results from the AstroGeoFit analysis,
    retrieves the selected result based on the number of knots to explore, and then generates
    plots comparing the raw data and the model predictions in both the depth and time domains.

    The plots can be saved in several formats including EPS, PNG, PDF, SVG, and table format (.tab).

    Parameters:
        - **configuration_file_path** (str): Path to the configuration file used in the analysis.
        - **folder_with_the_results** (Optional[str]): Path to the folder containing the fitting results.
        - **file_with_the_results** (Optional[str]): Name of the file with the fitting results.
        - **types_of_saving** (list): List of formats to save the figures and/or data.
            Supported formats: "png", "pdf", "svg", "eps", "tab".
        - **output_path** (str): Path where the figures/tables should be saved. If empty or "default",
            it uses the default output folder from the analysis parameters.

    Raises:
        ValueError: If the specified number of knots (`number_of_knots_to_explore`) is not found
            in the list of explored values.
        FileNotFoundError: If any required file is missing.
        KeyError: If expected keys are missing in the results or parameter dictionaries.

    Notes:
        - The function uses `lowess` smoothing for visualization purposes.
    """
    ################### LIBRARY MODE SECTION ####################
    results_fit, _, _, parameters_analysis, data_set_values, _, genetic_algorithm_parameters, _, _, _ = notebooks_code.data_obtain_analysis_file_and_result_files(False, configuration_file_path, folder_with_the_results, file_with_the_results)
    data, _ = data_manager.data_read_dataset_from_file(data_set_values)
    tool_plt_fit_time_series(
        genetic_algorithm_parameters=genetic_algorithm_parameters,
        results_fit=results_fit,
        data=data,
        parameters_analysis=parameters_analysis,
        number_of_knots_to_explore=number_of_knots_to_explore,
        types_of_saving=types_of_saving,
        output_path=output_path
    )

def tool_plt_fit_time_series(
    genetic_algorithm_parameters:dict,
    results_fit:dict,
    data: list,
    parameters_analysis: dict,
    number_of_knots_to_explore: Optional[int],
    types_of_saving: list,
    output_path: str,
):
    optimization_results = results_fit["optimization_results"]
    list_number_genes=genetic_algorithm_parameters["list_number_genes"]
    invSR_to_predx=results_fit["invSR_to_predx"]
    if number_of_knots_to_explore:
        try:
            index_result_to_explore = list_number_genes.index(
                number_of_knots_to_explore
            )
        except ValueError:
            logger.error(
                f"The selected number of genes {number_of_knots_to_explore} does not exist in the knots used for fitting."
            )
            raise ValueError
    else:
        index_result_to_explore = -1
    atributes = plot_utils.get_result_atributes(
        optimization_results[-1][index_result_to_explore]
    )

    fig, ax = plt.subplots(
        2, 1, figsize=(8, 4), constrained_layout=False, sharex=False, dpi=300
    )
    plt.subplots_adjust(
        left=0.1, right=0.95, bottom=0.15, top=0.9, wspace=0.1, hspace=0.5
    )

    # Get time scale from fitting results ====
    res_prediction = invSR_to_predx(
        depth_invSR=[atributes["depth_genes"], atributes["best_inverse_SR"]],
        fs=atributes["frequencies"],
    )
    time = res_prediction["time"]
    y_pred = res_prediction["y_pred"]
    # Smoothen the data, to better compare the model and the data ====
    lowess_span = 20 / len(data[1])  # 10-points smoothed
    lowess = sm.nonparametric.lowess(data[1], data[0], frac=lowess_span).T
    y_smooth = lowess[1]
    y_smooth = sp.signal.detrend(y_smooth)

    # Plot depth series ====
    ax[0].plot(data[0], y_pred, color="k", label="Model")  # Plot fitted data
    ax[0].plot(
        data[0], data[1], zorder=-10, label="Raw data", alpha=0.8
    )  # Plot "raw" data

    ax[0].legend()
    ax[0].set_xlabel("Depth (m)")
    ax[0].set_ylabel("Proxy")

    # Plot time series ====
    ax[1].plot(time, y_pred, color="k", label="Model")  # Plot fitted data
    ax[1].plot(time, data[1], zorder=-10, label="Data", alpha=0.8)  # Plot "raw" data

    # plt.legend()
    ax[1].set_xlabel("Time (Myr)")
    ax[1].set_ylabel("Proxy")

    if types_of_saving and types_of_saving != []:
        if not output_path or output_path == "default" or output_path == "":
            output_path = parameters_analysis["output_folder"]
            output_path = f"{output_path}/figures"
            os.makedirs(output_path, exist_ok=True)
        figname = f"Time-Depth_Series_{number_of_knots_to_explore}_Knots"
        global supported_formats
        for fmt in types_of_saving:
            if fmt.lower() in supported_formats:
                if fmt.lower() == "tab":
                    # Save time series data as a table
                    time_table_filename = f"{output_path}/{figname}_TimeSeries.tab"
                    with open(time_table_filename, "w") as file:
                        file.write("Time (Myr)\tModel\tData\n")
                        for t, yp, yd in zip(time, y_pred, data[1]):
                            file.write(f"{t}\t{yp}\t{yd}\n")
                    logger.info(f"Saved {time_table_filename} in format {fmt}")

                    # Save depth series data as a table
                    depth_table_filename = f"{output_path}/{figname}_DepthSeries.tab"
                    with open(depth_table_filename, "w") as file:
                        file.write("Depth (m)\tModel\tData\n")
                        for d, yp, yd in zip(data[0], y_pred, data[1]):
                            file.write(f"{d}\t{yp}\t{yd}\n")
                    logger.info(f"Saved {depth_plot_filename} in format {fmt}")

                elif fmt.lower() == "eps":
                    # Save time series plot as EPS
                    fig_time, ax_time = plt.subplots(figsize=(10, 5), dpi=300)
                    ax_time.plot(time, y_pred, color="k", label="Model")
                    ax_time.plot(time, data[1], zorder=-10, label="Data", alpha=0.8)
                    ax_time.set_xlabel("Time (Myr)")
                    ax_time.set_ylabel("Proxy")
                    ax_time.legend()
                    time_plot_filename = f"{output_path}/{figname}_TimeSeries.eps"
                    fig_time.savefig(time_plot_filename, format="eps")
                    plt.close(fig_time)
                    logger.info(f"Saved {output_path}/{figname}.{fmt.lower()}")

                    # Save depth series plot as EPS
                    fig_depth, ax_depth = plt.subplots(figsize=(10, 5), dpi=300)
                    ax_depth.plot(data[0], y_pred, color="k", label="Model")
                    ax_depth.plot(data[0], data[1], zorder=-10, label="Data", alpha=0.8)
                    ax_depth.set_xlabel("Depth (m)")
                    ax_depth.set_ylabel("Proxy")
                    ax_depth.legend()
                    depth_plot_filename = f"{output_path}/{figname}_DepthSeries.eps"
                    fig_depth.savefig(depth_plot_filename, format="eps")
                    plt.close(fig_depth)
                    logger.info(f"Saved {depth_plot_filename} in format {fmt}")

                else:
                    # Save combined plot in other formats
                    fig.savefig(
                        f"{output_path}/{figname}.{fmt.lower()}", format=fmt.lower()
                    )
                    logger.info(f"Saved {output_path}/{figname}.{fmt.lower()}")
            else:
                logger.warning(
                    f"The file type {fmt} is not supported. The figure will not be saved in this format."
                )

#################################################
########### ECCENTRICITY PLOTS PART #############
#################################################

def data_calculate_eccentricity_parameters(
    configuration_file_path:str,
    folder_with_the_results:Optional[str],
    file_with_the_results:Optional[str],
    positive_feedback: bool,
    number_of_knots_to_explore: Optional[int],
    name_eccentricity_solution: Optional[str],
):
    """
    Calculates and saves the eccentricity-related parameters from a time-depth model fit using 
    astronomical components (eccentricity, precession, tilt).

    This function loads the results of a genetic algorithm fit, reconstructs the inferred time scale, 
    and computes the contribution of each astronomical component using linear regression. 
    The resulting parameters and derived time series are stored in a temporary file.

    Parameters:
        - **configuration_file_path** : str
            Path to the configuration file used to perform the time-depth model fitting.

        - **folder_with_the_results** : Optional[str]
            Path to the folder where the results of the fitting procedure are stored. If None, 
            the default from the configuration file is used.

        - **file_with_the_results** : Optional[str]
            Specific filename of the results to load, if available. If None, uses the default.

        - **positive_feedback** : bool
            Indicates whether the eccentricity signal should be positively or negatively correlated 
            with the data (used to adjust the sign of the resulting prediction).

        - **number_of_knots_to_explore** : Optional[int]
            Number of knots (genes) used in the fitting process to be explored. If None, the last 
            set of results is used by default.

        - **name_eccentricity_solution** : Optional[str]
            Custom name assigned to the eccentricity solution (e.g., "La2004", "Berger2020"). 
            Used for labeling and tracking. If None or empty, defaults to "Astronomical Solution".

    Returns:
        None
            Saves a dictionary of computed parameters and time series in a temporary `.pickle` file. 
            Keys in the dictionary may include:
            - 'feedback_factor'
            - 'name_eccentricity_solution'
            - 'number_of_knots_to_explore'
            - 'times_inferred'
            - 'y_envelope', 'y_envelope_norm'
            - 'y_pred_prec', 'y_pred_ecc', 'y_pred_ecc_norm', 'y_pred_tilt'
    
    Notes:
        - The time scale is reconstructed by integrating the interpolated inverse sedimentation rate (invSR).
        - Regression is performed using sine and cosine components of astronomical frequencies.
        - Only components enabled in `frequency_values` (precession, eccentricity, tilt) are processed.
        - The results are stored under the path: `<current working directory>/tmp/temp_eccentricity_params.pickle`. This file should not be eliminated
    """
    ################### LIBRARY MODE SECTION ####################
    results_fit, _, _, _, data_set_values, _, genetic_algorithm_parameters, data_model_parameters, _, _ = notebooks_code.data_obtain_analysis_file_and_result_files(False, configuration_file_path, folder_with_the_results, file_with_the_results)
    data, _ = data_manager.data_read_dataset_from_file(data_set_values)
    tool_data_calculate_eccentricity_parameters(
        genetic_algorithm_parameters=genetic_algorithm_parameters,
        data_model_parameters=data_model_parameters,
        results_fit=results_fit,
        data=data,
        positive_feedback=positive_feedback,
        number_of_knots_to_explore=number_of_knots_to_explore,
        name_eccentricity_solution=name_eccentricity_solution
    )

def tool_data_calculate_eccentricity_parameters(
    genetic_algorithm_parameters:dict,
    data_model_parameters:dict,
    results_fit:dict,
    data: list,
    positive_feedback: bool,
    number_of_knots_to_explore: Optional[int],
    name_eccentricity_solution: Optional[str],
    script_path: str = "",
):
    optimization_results=results_fit["optimization_results"]
    frequency_values=data_model_parameters["frequency_values"]
    list_number_of_knots=genetic_algorithm_parameters["list_number_genes"]
    
    dict_to_save = {}
    feedback_factor = 1 if positive_feedback else -1
    dict_to_save["feedback_factor"] = feedback_factor
    if not name_eccentricity_solution or name_eccentricity_solution == "":
        name_eccentricity_solution = "Astronomical Solution"
    dict_to_save["name_eccentricity_solution"] = name_eccentricity_solution
    #################################################################################
    if number_of_knots_to_explore:
        try:
            index_result_to_explore = list_number_of_knots.index(
                number_of_knots_to_explore
            )
        except ValueError:
            logger.error(
                f"The selected number of genes {number_of_knots_to_explore} does not exist in the genes used for fitting."
            )
            raise ValueError
    else:
        index_result_to_explore = -1
    atributes = plot_utils.get_result_atributes(
        optimization_results[-1][index_result_to_explore]
    )
    if not number_of_knots_to_explore:
        number_of_knots_to_explore = len(atributes["depth_genes"])
    dict_to_save["number_of_knots_to_explore"] = number_of_knots_to_explore
    invSR_interpolate = atributes["interpolator"](
        [atributes["depth_genes"], atributes["best_inverse_SR"]], data[0]
    )
    times_inferred = sp.integrate.cumulative_trapezoid(
        invSR_interpolate, data[0], initial=0
    )
    dict_to_save["times_inferred"] = times_inferred
    X = shared_functions.generate_X_linReg(
        np.ones_like(atributes["frequencies"]), atributes["frequencies"], times_inferred
    )  # Generate X
    reg_model = LinearRegression().fit(X, data[1])  # Get the regression funciton.
    y_pred = X @ reg_model.coef_

    # Get coefficients
    freq_index_final = 0
    if frequency_values["use_precession"]:
        freq_index_final = (
            freq_index_final + 5
        )  # THE ADDITION OF 5 IS BECAUSE EACH TIME THERE IS 5 FREQUENCIES OF EACH TYPE.
        AB_prec = np.hstack(
            [
                reg_model.coef_[:freq_index_final],
                reg_model.coef_[
                    len(atributes["frequencies"]) : len(atributes["frequencies"])
                    + freq_index_final
                ],
            ]
        )
        X_prec = shared_functions.generate_sine_waves(
            np.ones_like(atributes["frequencies"][:freq_index_final]),
            atributes["frequencies"][:freq_index_final],
            times_inferred,
        )
        X_prec_90 = shared_functions.generate_sine_waves(
            np.ones_like(atributes["frequencies"][:freq_index_final]),
            atributes["frequencies"][:freq_index_final],
            times_inferred,
            phi=np.pi / 2,
        )
        # Predict Y for ETP
        y_pred_prec = X_prec @ AB_prec
        y_pred_prec_90 = X_prec_90 @ AB_prec
        y_envelope = np.sqrt(y_pred_prec_90**2 + y_pred_prec**2)
        y_envelope_norm = (y_envelope - y_envelope.mean()) / y_envelope.std()
        dict_to_save["y_envelope"] = y_envelope
        dict_to_save["y_pred_prec"] = y_pred_prec
        dict_to_save["y_envelope_norm"] = y_envelope_norm

    if frequency_values["use_eccentricity"]:
        freq_index_ini = freq_index_final
        freq_index_final = freq_index_final + 5
        AB_ecc = np.hstack(
            [
                reg_model.coef_[freq_index_ini:freq_index_final],
                reg_model.coef_[
                    len(atributes["frequencies"])
                    + freq_index_ini : len(atributes["frequencies"])
                    + freq_index_final
                ],
            ]
        )
        X_ecc = shared_functions.generate_sine_waves(
            np.ones_like(atributes["frequencies"][freq_index_ini:freq_index_final]),
            atributes["frequencies"][freq_index_ini:freq_index_final],
            times_inferred,
        )
        y_pred_ecc = X_ecc @ AB_ecc
        y_pred_ecc = y_pred_ecc * feedback_factor
        y_pred_ecc_norm = (y_pred_ecc - y_pred_ecc.mean()) / y_pred_ecc.std()
        dict_to_save["y_pred_ecc"] = y_pred_ecc
        dict_to_save["y_pred_ecc_norm"] = y_pred_ecc_norm

    if frequency_values["use_tilt"]:
        freq_index_ini = freq_index_final
        freq_index_final = freq_index_ini + (
            len(atributes["frequencies"]) - freq_index_ini
        )  # THIS GENERATES A PROBLEM NOW BC THERE ARE ONLY TWO S FREQUENCIES BEING USED. TEST MORE AND TAKE CARE
        AB_tilt = np.hstack(
            [
                reg_model.coef_[freq_index_ini:freq_index_final],
                reg_model.coef_[
                    len(atributes["frequencies"])
                    + freq_index_ini : len(atributes["frequencies"])
                    + freq_index_final
                ],
            ]
        )
        X_tilt = shared_functions.generate_sine_waves(
            np.ones_like(atributes["frequencies"][freq_index_ini:freq_index_final]),
            atributes["frequencies"][freq_index_ini:freq_index_final],
            times_inferred,
        )
        y_pred_tilt = X_tilt @ AB_tilt
        dict_to_save["y_pred_tilt"] = y_pred_tilt

    # y_envelope = abs(sp.signal.hilbert(y_pred_prec))

    if script_path == "":
        temp_path = f"{os.getcwd()}/tmp"
    else:
        temp_path = f"{script_path}/tmp"
    temp_file = "temp_eccentricity_params"
    if os.path.exists(f"{temp_path}/{temp_file}.pickle"):
        os.remove(f"{temp_path}/{temp_file}.pickle")
    if not os.path.exists(temp_path):
        os.mkdir(temp_path)
    data_manager.save_data(
        output_file_path=f"{temp_path}/{temp_file}",
        fit_mcmc_weights=None,
        **dict_to_save,
    )

def plt_fit_ETP_signals_from_model(
    configuration_file_path:str,
    folder_with_the_results:Optional[str],
    file_with_the_results:Optional[str],
    types_of_saving: list,
    output_path: str,
):
    """
    Plots the astronomical components (eccentricity, precession, and tilt) from a previously 
    calculated model fit, and optionally saves the resulting figure in multiple formats.

    This function loads a set of inferred astronomical signals and their associated time scale, 
    and generates a plot showing:
    - The fitted precession signal and its envelope.
    - The fitted eccentricity signal.
    - The fitted tilt signal (if available).
    - The correlation (r and r²) between the normalized eccentricity and the normalized precession envelope.

    Parameters:
        - **configuration_file_path** : str
            Path to the configuration file used during the model fitting process.

        - **folder_with_the_results** : Optional[str]
            Path to the folder containing the results. If None, uses the path from the configuration file.

        - **file_with_the_results** : Optional[str]
            Name of the specific results file to load. If None, uses the default defined in the configuration.

        - **types_of_saving** : list
            List of file formats to save the figure in (e.g., ["png", "pdf"]). 
            The format must be one of the supported formats (excluding "tab").

        - **output_path** : str
            Output directory where the figure(s) should be saved. If "default", "", or None, 
            it uses the default output folder from the analysis parameters.

    Returns:
        None
            The function saves the figure(s) in the specified formats if `types_of_saving` is provided. 
            The figure shows the comparison between fitted astronomical components and the inferred time scale.

    Notes:
        - The function reads from a temporary `.pickle` file previously saved by `data_calculate_eccentricity_parameters`.
        - The correlation metrics (r, r²) are only shown if both eccentricity and precession are enabled.
        - Supported formats are defined in the global variable `supported_formats`.
        - If `tab` is included in the `types_of_saving`, a warning is issued (not a supported plot format).
        - The plot shows inferred signals over time (in millions of years).
    """

    ################### LIBRARY MODE SECTION ####################
    _, _, _, parameters_analysis, _, _, _, data_model_parameters, _, _ = notebooks_code.data_obtain_analysis_file_and_result_files(False, configuration_file_path, folder_with_the_results, file_with_the_results)
    
    tool_plt_fit_ETP_signals_from_model(
        data_model_parameters=data_model_parameters,
        parameters_analysis=parameters_analysis,
        types_of_saving=types_of_saving,
        output_path=output_path
    )

def tool_plt_fit_ETP_signals_from_model(
    data_model_parameters:dict,
    parameters_analysis: dict,
    types_of_saving: list,
    output_path: str,
    script_path: str = "",
):
    frequency_values=data_model_parameters["frequency_values"]
    if script_path == "":
        temp_path = f"{os.getcwd()}/tmp/temp_eccentricity_params"
    else:
        temp_path = f"{script_path}/tmp/temp_eccentricity_params"
    results_dict = data_manager.load_dictionary_from_pickle(
        file_path=f"{temp_path}.pickle"
    )
    times_inferred = results_dict["times_inferred"]
    # Plot ====
    fig, ax = plt.subplots(1, 1, figsize=(10, 4), constrained_layout=True, sharex=True)
    y_text = np.percentile(ax.get_ylim(), 90)
    all_ys = []
    if frequency_values["use_precession"]:
        y_pred_prec = results_dict["y_pred_prec"]
        y_envelope = results_dict["y_envelope"]
        ax.plot(times_inferred, y_pred_prec, label="prec")
        ax.plot(times_inferred, y_envelope, label="env prec")
        all_ys.extend(results_dict["y_pred_prec"])
        all_ys.extend(results_dict["y_envelope"])

    if frequency_values["use_eccentricity"]:
        y_pred_ecc = results_dict["y_pred_ecc"]
        ax.plot(times_inferred, y_pred_ecc, label="ecc", zorder=100)
        all_ys.extend(results_dict["y_pred_ecc"])
        if frequency_values["use_precession"]:
            y_envelope_norm = results_dict["y_envelope_norm"]
            y_pred_ecc_norm = results_dict["y_pred_ecc_norm"]
            r2 = pearsonr(y_envelope_norm, y_pred_ecc_norm)[0] ** 2
            r = pearsonr(y_envelope_norm, y_pred_ecc_norm)[0]

    if frequency_values["use_tilt"]:
        y_pred_tilt = results_dict["y_pred_tilt"]
        ax.plot(times_inferred, y_pred_tilt, label="tilt")
        all_ys.extend(results_dict["y_pred_tilt"])

    min_ys = np.min(all_ys)
    max_ys = np.max(all_ys)
    buffer = 0.5 * (max_ys - min_ys)
    ax.set_ylim(min_ys - buffer, max_ys+ buffer)
    
    if frequency_values["use_eccentricity"] and frequency_values["use_precession"]:
        ax.text(0, max_ys + buffer * 0.65, "Correlation between ecc and env prec:", fontsize=12, color="r")
        ax.text(0, max_ys + buffer * 0.3, "r$^2$ = " + str(np.round(r2, 2)), fontsize=12, color="r")
        ax.text(0, max_ys + buffer * 0.1, "r   = " + str(np.round(r, 2)), fontsize=12, color="r")
    
    plt.legend(loc="upper right")  # Updated by YW, Nov 12, 2024
    ax.set_xlabel("Time (Myr)")
    if types_of_saving and types_of_saving != []:
        number_of_knots_explored = results_dict["number_of_knots_to_explore"]
        if not output_path or output_path == "default" or output_path == "":
            output_path = parameters_analysis["output_folder"]
            output_path = f"{output_path}/figures"
            os.makedirs(output_path, exist_ok=True)
        figname = f"Correlation_Between_Eccentricity_and_Envelope_Precession_{number_of_knots_explored}_Knots"
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


def plt_fit_eccentricity_mswd_plot(
    configuration_file_path:str,
    folder_with_the_results:Optional[str],
    file_with_the_results:Optional[str],
    time_range: Optional[list],
    types_of_saving: list,
    output_path: str,
):
    """
    Plots the Mean Squared Weighted Deviation (MSWD) between the fitted eccentricity/precession envelope
    and a reference astronomical solution, as a function of t₀ (the starting time of the inferred time scale).

    This function evaluates how well the inferred normalized eccentricity and precession envelope align
    with the chosen astronomical solution over a range of t₀ values. It produces two subplots: one for
    eccentricity and another for the precession envelope. The plots highlight the best-fit t₀ values where
    the MSWD is minimized.

    Parameters:
        - **configuration_file_path** : str
            Path to the configuration file used during the model fitting process.

        - **folder_with_the_results** : Optional[str]
            Path to the folder containing the result files. If None, the folder is inferred from the configuration.

        - **file_with_the_results** : Optional[str]
            Name of the result file to be loaded. If None, it defaults to the one defined in the configuration.

        - **time_range** : Optional[list]
            A list containing two floats or integers [tmin, tmax] defining the time range to search for optimal t₀.
            If None, the function uses a default time range based on the nominal top age of the dataset.

        - **types_of_saving** : list
            A list of file formats to save the resulting plot (e.g., ["png", "pdf"]). Format must be in 
            the supported formats list (excluding "tab").

        - **output_path** : str
            Directory where the output figure(s) and log file will be saved. If "default", "", or None, the output 
            path is inferred from the configuration file's output folder.

    Returns:
        None
            The function saves the MSWD plot for eccentricity and precession envelope and logs the relevant 
            values (e.g., optimal t₀, top/bottom ages) into a text file. It also saves the t₀ values in a temporary
            pickle for further use.

    Notes:
        - The MSWD is computed as the mean of squared residuals normalized by signal variance for each t₀ value.
        - The eccentricity parameters are interpolated from the selected astronomical solution using the 
        `shared_functions.get_eccentricity_parameters`.
        - Two optimal t₀ values are extracted: one minimizing MSWD for eccentricity, the other for the precession envelope.
        - A text file with a full summary of the MSWD results and t₀s is saved in the output directory.
        - A dictionary with `t0_ecc_MSWD` and `t0_env_MSWD` is saved as a pickle file for downstream processes.
    """
    ################### LIBRARY MODE SECTION ####################
    _, _, _, parameters_analysis, data_set_values, age_depth_model_values, _, data_model_parameters, _, eccentricity_solution_parameters = notebooks_code.data_obtain_analysis_file_and_result_files(False, configuration_file_path, folder_with_the_results, file_with_the_results)
    
    _, _, _, _, func_time_nominal, _ = notebooks_code.data_obtain_age_depth_values(age_depth_model_values)

    data, _ = data_manager.data_read_dataset_from_file(data_set_values)
    
    tool_plt_fit_eccentricity_mswd_plot(
        eccentricity_params=eccentricity_solution_parameters,
        data_model_parameters=data_model_parameters,
        data=data,
        func_time_nominal=func_time_nominal,
        time_min_and_max=time_range,
        parameters_analysis=parameters_analysis,
        types_of_saving=types_of_saving,
        output_path=output_path
    )
    
def tool_plt_fit_eccentricity_mswd_plot(
    eccentricity_params: dict,
    data_model_parameters:dict,
    data: list,
    func_time_nominal,
    time_min_and_max: Optional[list],
    parameters_analysis: dict,
    types_of_saving: list,
    output_path: str,
    script_path: str ="",
):
    if not isinstance(eccentricity_params, dict):
        logger.error("No Eccentricity Solution was provided. The plot can't be processed.")
        return
    frequency_values=data_model_parameters["frequency_values"]
    if script_path == "":
        temp_path = f"{os.getcwd()}/tmp/temp_eccentricity_params"
    else:
        temp_path = f"{script_path}/tmp/temp_eccentricity_params"
    results_dict = data_manager.load_dictionary_from_pickle(
        file_path=f"{temp_path}.pickle"
    )
    name_eccentricity_solution = results_dict["name_eccentricity_solution"]
    # Set a range for t0. "t0" is the starting time of the time scale.  # Updated by YW, Nov 12, 2024
    eccentricity_interp = shared_functions.get_eccentricity_parameters(
        eccentricity_params
    )
    if time_min_and_max:
        if len(time_min_and_max) != 2 or not (
            all(isinstance(item, (int, float)) for item in time_min_and_max)
        ):
            print(
                f"The starting range {time_min_and_max} is incorrect. Give a list of two values that are floats or integers."
            )
            return
        else:
            tmin = time_min_and_max[0]
            tmax = time_min_and_max[1]
    else:
        tmin = round(func_time_nominal(data[0])[0], 1) - 1  # default tmin
        tmax = round(func_time_nominal(data[0])[0], 1) + 1  # default tmax

    t0s = np.linspace(tmin, tmax, 2001)
    times_inferred = results_dict["times_inferred"]
    # For directly-derived eccentricity
    if frequency_values["use_eccentricity"]:
        y_pred_ecc_norm = results_dict["y_pred_ecc_norm"]
        var_e3_est_ecc = np.nanvar(y_pred_ecc_norm, axis=0)
        chi_squared_t0s_ecc = np.mean(
            (eccentricity_interp(times_inferred + t0s[:, None]) - y_pred_ecc_norm) ** 2
            / var_e3_est_ecc[None],
            axis=1,
        )

    # For precession envelope
    if frequency_values["use_precession"]:
        y_envelope_norm = results_dict["y_envelope_norm"]
        var_e3_est_env = np.nanvar(y_envelope_norm, axis=0)
        chi_squared_t0s_env = np.mean(
            (eccentricity_interp(times_inferred + t0s[:, None]) - y_envelope_norm) ** 2
            / var_e3_est_env[None],
            axis=1,
        )

    # ---------------------------------------------------------------------------

    # Find optimal t0 in the given range
    bools = (t0s > tmin) & (t0s < tmax)
    t0s_n = t0s[bools]
    t0_ecc_MSWD = t0s_n[np.argmin(chi_squared_t0s_ecc[bools])]
    if frequency_values["use_precession"]: # ----YW, 2025-Apr-10
        t0_env_MSWD = t0s_n[np.argmin(chi_squared_t0s_env[bools])]
    #t0_ecc_MSWD, t0_env_MSWD

    # Plot
    fig, axs = plt.subplots(2, 1, figsize=(8, 4), constrained_layout=True)

    # For directly derived eccentricity
    ax = axs[0]
    ax.plot(
        t0s,
        chi_squared_t0s_ecc,
        label=f"Eccentricity vs. {name_eccentricity_solution}",
        color="#0072B2",
    )
    ax.axvline(t0_ecc_MSWD, linestyle="--", linewidth=1.5, color="r", zorder=-1)
    ax.set_ylabel("MSWD")
    ax.set_xlabel("t0 [Ma]")
    ax.legend(loc="upper right")

    # For precession envelope
    if frequency_values["use_precession"]:
        ax=axs[1]
        ax.plot(t0s, chi_squared_t0s_env,label=f'Precesion env vs. {name_eccentricity_solution}' ,color='#CC79A7')
        ax.axvline(t0_env_MSWD, linestyle="--", linewidth=1.5, color='r', zorder=-1)
        ax.set_ylabel('MSWD')
        ax.set_xlabel('t0 [Ma]')
        ax.legend(loc='upper right')
    else: axs[1].remove()

    # LOGS OF THE PLOTS
    if not output_path or output_path == "default" or output_path == "":
        output_path_logs = parameters_analysis["output_folder"]
        os.makedirs(output_path_logs, exist_ok=True)
    else:
        output_path_logs = output_path
    file_path = f"{output_path_logs}/logs_plots.txt"
    plot_utils.print_and_save(
        message="###########################################", file_path=file_path
    )
    plot_utils.print_and_save(message="", file_path=file_path)
    plot_utils.print_and_save(message="MSWD PLOT", print_time=True, file_path=file_path)
    plot_utils.print_and_save(message="", file_path=file_path)
    message = f"Time range used: {tmin}, {tmax}"
    plot_utils.print_and_save(message=message, file_path=file_path)
    plot_utils.print_and_save(message="", file_path=file_path)
    message = f"Astronomical Solution used: {name_eccentricity_solution}"
    plot_utils.print_and_save(message=message, file_path=file_path)
    plot_utils.print_and_save(message="", file_path=file_path)
    plot_utils.print_and_save(
        message="--------------- Nominal ---------------", file_path=file_path
    )
    plot_utils.print_and_save(message="", file_path=file_path)
    message = f"Top age: {func_time_nominal(data[0])[0]} ; Bottom age: {func_time_nominal(data[0])[-1]}"
    plot_utils.print_and_save(message=message, file_path=file_path)
    plot_utils.print_and_save(message="", file_path=file_path)
    message = "--------- Eccentricity (gi-gj) --------"
    plot_utils.print_and_save(message=message, file_path=file_path)
    plot_utils.print_and_save(message="", file_path=file_path)
    message = (
        f"Top age: t0 = {t0_ecc_MSWD} ; Bottom age: {times_inferred[-1]+t0_ecc_MSWD}"
    )
    plot_utils.print_and_save(message=message, file_path=file_path)
    plot_utils.print_and_save(message="", file_path=file_path)
    
    if frequency_values["use_precession"]: # ----YW, 2025-Apr-10
        message = '--------- Precession envelope ---------'
        plot_utils.print_and_save(message=message, file_path=file_path)
        plot_utils.print_and_save(message="", file_path=file_path)
        message = f'Top age: t0 = {t0_env_MSWD} ; Bottom age: {times_inferred[-1]+t0_env_MSWD}'
        plot_utils.print_and_save(message=message, file_path=file_path)
        plot_utils.print_and_save(message="", file_path=file_path)
        
    plot_utils.print_and_save(message="###########################################", file_path=file_path)
    dict_to_save = {"t0_ecc_MSWD":t0_ecc_MSWD}

    if frequency_values["use_precession"]: # ----YW, 2025-Apr-10
        dict_to_save["t0_env_MSWD"] = t0_env_MSWD
        
    data_manager.save_data(
        output_file_path=temp_path, fit_mcmc_weights=None, **dict_to_save
    )
    if types_of_saving and types_of_saving != []:
        number_of_knots_explored = results_dict["number_of_knots_to_explore"]
        if not output_path or output_path == "default" or output_path == "":
            output_path = parameters_analysis["output_folder"]
            output_path = f"{output_path}/figures"
            os.makedirs(output_path, exist_ok=True)
        figname = f"MSWD_Eccentricity_and_Precession_{number_of_knots_explored}_Knots_{name_eccentricity_solution}"
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

def plt_fit_correlation_eccentricity_and_solution(
    configuration_file_path:str,
    folder_with_the_results:Optional[str],
    file_with_the_results:Optional[str],
    x_axis_limits: Optional[list],
    types_of_saving: list,
    output_path: str,
):
    """
    Plot the correlation between the fitted eccentricity or precession envelope 
    and the astronomical solution used in the AstroGeoFit modeling.

    This function reads the analysis configuration, dataset, and results of the 
    eccentricity or precession envelope fitting. It then compares these model results 
    against the corresponding astronomical solution and plots both the fitted curves 
    and the astronomical signal. Correlation metrics (r and r²) are calculated and 
    displayed on the plots.

    Parameters:
        - **configuration_file_path** : str
            Path to the configuration file used for the AstroGeoFit analysis.

        - **folder_with_the_results** : str or None
            Path to the folder containing the result files. If None, the default folder 
            defined in the configuration file is used.

        - **file_with_the_results** : str or None
            Specific filename for the results. If None, a default filename is assumed.

        - **x_axis_limits** : list or None
            Optional list specifying the [xmin, xmax] limits for the x-axis (in Ma). 
            If None, limits are automatically inferred from the data.

        - **types_of_saving** : list of str
            List of output file formats for saving the figure (e.g., ["png", "pdf"]). 
            If empty, the figure is not saved.

        - **output_path** : str
            Path where the output figure will be saved. If "default" or empty, a default 
            subdirectory `figures` in the output folder from the configuration is used.

    Notes:
        - This function assumes that the eccentricity solution and inferred signals have 
        already been computed and stored in a temporary pickle file.
        - If `frequency_values["use_eccentricity"]` or 
        `frequency_values["use_precession"]` are set in the configuration, the function 
        will generate the corresponding subplots.
        - Pearson correlation coefficient `r` and coefficient of determination `r²` are 
        shown to quantify the fit between the modeled and astronomical signals.

    Outputs:
        - A matplotlib figure with one or two subplots showing the fitted signal(s) 
        and astronomical solution.
        - Saved figure(s) in the specified format(s) and location(s) if requested.

    """
    ################### LIBRARY MODE SECTION ####################
    _, _, _, parameters_analysis, data_set_values, age_depth_model_values, _, data_model_parameters, _, eccentricity_solution_parameters = notebooks_code.data_obtain_analysis_file_and_result_files(False, configuration_file_path, folder_with_the_results, file_with_the_results)
    
    _, _, _, _, func_time_nominal, _ = notebooks_code.data_obtain_age_depth_values(age_depth_model_values)

    data, _ = data_manager.data_read_dataset_from_file(data_set_values)
    tool_plt_fit_correlation_eccentricity_and_solution(
        eccentricity_params=eccentricity_solution_parameters,
        data_model_parameters=data_model_parameters,
        data=data,
        func_time_nominal=func_time_nominal,
        parameters_analysis=parameters_analysis,
        x_axis_limits=x_axis_limits,
        types_of_saving=types_of_saving,
        output_path=output_path
    )

def tool_plt_fit_correlation_eccentricity_and_solution(
    eccentricity_params: dict,
    data_model_parameters:dict,
    data: list,
    func_time_nominal,
    parameters_analysis: dict,
    x_axis_limits: Optional[list],
    types_of_saving: list,
    output_path: str,
    script_path: str = "",
):
    if not isinstance(eccentricity_params, dict):
        logger.error("No Eccentricity Solution was provided. The plot can't be processed.")
        return
    frequency_values=data_model_parameters["frequency_values"]
    if script_path == "":
        temp_path = f"{os.getcwd()}/tmp/temp_eccentricity_params"
    else:
        temp_path = f"{script_path}/tmp/temp_eccentricity_params"
    results_dict = data_manager.load_dictionary_from_pickle(
        file_path=f"{temp_path}.pickle"
    )
    name_eccentricity_solution = results_dict["name_eccentricity_solution"]
    # Plot the correlation between fitted ecc and astronomical solution ========

    if frequency_values["use_eccentricity"]:
        if frequency_values["use_precession"]:   
            t0_ecc,t0_env=results_dict["t0_ecc_MSWD"],results_dict["t0_env_MSWD"]
        
            if not x_axis_limits:
                x_min=round(min(t0_ecc,t0_env),1)-0.5                        # default xmin
                x_max=round(max(t0_ecc,t0_env)+results_dict["times_inferred"][-1],1)+0.5     # default xmax
            else:
                x_min = x_axis_limits[0]
                x_max = x_axis_limits[1]
        else:
            t0_ecc=results_dict["t0_ecc_MSWD"]
            
            if not x_axis_limits:
                x_min=round(t0_ecc,1)-0.5                        # default xmin
                x_max=round(t0_ecc+results_dict["times_inferred"][-1],1)+0.5     # default xmax
            else:
                x_min = x_axis_limits[0]
                x_max = x_axis_limits[1]
    # ---------------------------------------------------------------------------
    ib, di = 20, 20  # What we have to do here? Are these useful?
    e3_data, _ = data_manager.data_read_dataset_from_file(eccentricity_params, False, False, True)
    e3_data[0] = -e3_data[0]
    t0, tf = float(eccentricity_params["start_time"]), float(
        eccentricity_params["final_time"]
    )
    inds = (e3_data[0] >= t0) & (e3_data[0] <= tf)
    time_e3, e3 = e3_data[:, inds]
    e3 = (e3 - e3.mean()) / e3.std()
    fig, axs = plt.subplots(2, 1, figsize=(10, 4), constrained_layout=True, dpi=300)

    # For directly derived eccentricity ----------------------------------------
    if frequency_values["use_eccentricity"]:
        ax = axs[0]
        t0 = t0_ecc

        ax.plot(
            time_e3 / 1e3, e3, label=name_eccentricity_solution, color="k", linewidth=1
        )
        ax.plot(
            results_dict["times_inferred"] + t0,
            results_dict["y_pred_ecc_norm"],
            label="ecc",
            color="#0072B2",
        )

        e3_wd = shared_functions.interpolate_CubicSpline(
            [time_e3 / 1e3, e3], results_dict["times_inferred"][ib:-di] + t0
        )
        ax.axvline(
            func_time_nominal(data[0])[0], linestyle="--", color="r", linewidth=1
        )  # Updated by YW, Nov 12, 2024
        ax.axvline(
            func_time_nominal(data[0])[-1], linestyle="--", color="r", linewidth=1
        )  # Updated by YW, Nov 12, 2024

        ymin_ecc = min(
            min(e3),  # Min of astronomic solution
            min(results_dict["y_pred_ecc_norm"][ib:-di]),  # Min of fitted eccentricity
        )
        ymax_ecc = max(max(e3), max(results_dict["y_pred_ecc_norm"][ib:-di]))

        ax.set_xlabel("Age [Ma]")
        ax.set_ylabel("Eccentricity")
        ax.legend()
        ax.set_xlim([x_min, x_max])
        ax.legend(loc="upper right")
        ax.set_ylim(ymin_ecc, ymax_ecc * 1.5)

        r2 = pearsonr(e3_wd, results_dict["y_pred_ecc_norm"][ib:-di])[0] ** 2
        ax.text(
            x_min + 0.1,
            ymax_ecc * 1.2,
            "r$^2$ =" + str(np.round(r2, 2)),
            fontsize=12,
            color="r",
        )
        r = pearsonr(e3_wd, results_dict["y_pred_ecc_norm"][ib:-di])[0]
        ax.text(
            x_min + 0.1,
            ymax_ecc * 0.8,
            "r  =" + str(np.round(r, 2)),
            fontsize=12,
            color="r",
        )

    # For precession envelope --------------------------------------------------
    if frequency_values["use_precession"]:
        ax = axs[1]
        t0 = t0_env

        ax.plot(
            time_e3 / 1e3, e3, label=name_eccentricity_solution, color="k", linewidth=1
        )
        ax.plot(
            results_dict["times_inferred"] + t0,
            results_dict["y_envelope_norm"],
            label="prec env",
            color="#CC79A7",
        )

        e3_wd = shared_functions.interpolate_CubicSpline(
            [time_e3 / 1e3, e3], results_dict["times_inferred"][ib:-di] + t0
        )
        ax.axvline(
            func_time_nominal(data[0])[0], linestyle="--", color="r", linewidth=1
        )  # Updated by YW, Nov 12, 2024
        ax.axvline(
            func_time_nominal(data[0])[-1], linestyle="--", color="r", linewidth=1
        )  # Updated by YW, Nov 12, 2024

        ymin_prec = min(min(e3), min(results_dict["y_envelope_norm"][ib:-di]))
        ymax_prec = max(max(e3), max(results_dict["y_envelope_norm"][ib:-di]))

        ax.set_ylim(ymin_prec, ymax_prec * 1.5)
        ax.set_xlabel("Age [Ma]")
        ax.set_ylabel("Eccentricity")
        plt.legend()
        ax.set_xlim([x_min, x_max])
        ax.legend(loc="upper right")

        r2 = pearsonr(e3_wd, results_dict["y_envelope_norm"][ib:-di])[0] ** 2
        ax.text(
            x_min + 0.1,
            ymax_prec * 1.2,
            "r$^2$ =" + str(np.round(r2, 2)),
            fontsize=12,
            color="r",
        )
        r = pearsonr(e3_wd, results_dict["y_envelope_norm"][ib:-di])[0]
        ax.text(
            x_min + 0.1,
            ymax_prec * 0.8,
            "r  =" + str(np.round(r, 2)),
            fontsize=12,
            color="r",
        )
    else: axs[1].remove() # ----YW, 2025-Apr-10
    if types_of_saving and types_of_saving != []:
        number_of_knots_explored = results_dict["number_of_knots_to_explore"]
        if not output_path or output_path == "default" or output_path == "":
            output_path = parameters_analysis["output_folder"]
            output_path = f"{output_path}/figures"
            os.makedirs(output_path, exist_ok=True)
        figname = f"Eccentricity_Plot_{number_of_knots_explored}_Knots_{name_eccentricity_solution}"
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
