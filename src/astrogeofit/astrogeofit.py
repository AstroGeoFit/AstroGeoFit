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
import shutil
import numpy as np
import warnings

from colorama import init, Fore, Style  # type: ignore
from typing import Optional
from functools import partial

from astrogeofit.utils import shared_functions, exceptions, variable_check
from astrogeofit.utils.setup_logger import setup_logger
from astrogeofit.mcmc_routines import utils as utils_rmcmc
from astrogeofit.main_routines import (
    data_manager,
    main_fitting,
    mcmc_calcula,
    weight_computation,
    main_significance_test,
)
from astrogeofit.plots.plots_post_fitting import plots_post_fit as AGF_Plots_Post_Fit
from astrogeofit.notebooks_code import notebooks_code as AGF_Notebooks
from astrogeofit.plots.plots_post_mcmc import plots_post_mcmc as AGF_MCMC_Plots
from astrogeofit.plots.utils import shared_utils_plots as AGF_Shared_Functions
from astrogeofit.plots.significance_test_plot import significance_test_plot as AGF_Plot_Significance_Test

# Ignore all warnings
warnings.filterwarnings("ignore", message="Unable to import Axes3D.*")

# Set up the logger using the centralized config

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

logger = setup_logger(False)

def _info_beginning(
    config_file_path: str,
    significance: bool,
    fit: bool,
    mcmc: bool,
    weights: bool,
    execute_all: bool,
    remote: bool,
):
    init(autoreset=True)
    print("")
    print("##########################################")
    logger.info(f"{Fore.GREEN}EXECUTING ASTROGEOFIT TOOL ANALYSIS{Style.RESET_ALL}")
    print("##########################################")
    print("")
    config_analysis = data_manager.read_yaml_file(config_file_path)
    if config_analysis == None:
        raise exceptions.FileNotFound(
            message=f"The configuration file was not found in the path: {config_file_path}. Check that the introduced path is correct."
        )
    logger.info(
        f"{Fore.GREEN} Reading configuration file from {config_file_path} {Style.RESET_ALL}"
    )
    print("")

    output_folder = config_analysis["output_folder"]
    output_file = config_analysis["output_file_name"]
    if execute_all:
        # Check if the output folder exists
        if (
            os.path.exists(f"{output_folder}/analysis_parameters_{output_file}.yaml")
            and not remote
        ):
            # Folder exists, ask user if they want to continue
            user_input = (
                input(
                    f"{Fore.RED}WARNING: The folder {output_folder} with the results in {output_file} already exists. Do you want to continue and delete the folder? (y/n): {Style.RESET_ALL}"
                )
                .strip()
                .lower()
            )

            if user_input != "y":
                raise exceptions.ExecutionStoppedByUser()
            else:
                shutil.rmtree(config_analysis["output_folder"])
        elif (
            os.path.exists(f"{output_folder}/analysis_parameters_{output_file}.yaml")
            and remote
        ):
            logger.warning(
                f"{Fore.RED}WARNING: The folder {output_folder} already existed. In order to execute the new analysis it has been eliminated."
            )
            shutil.rmtree(config_analysis["output_folder"])
    if fit:
        # Check if the output folder exists
        if os.path.exists(f"{output_folder}/{output_file}_fit.pickle") and not remote:
            # Folder exists, ask user if they want to continue
            user_input = (
                input(
                    f"{Fore.RED}WARNING: The file {output_file}_fit.pickle already exists. Do you want to continue and delete overwrite the results? (y/n): {Style.RESET_ALL}"
                )
                .strip()
                .lower()
            )

            if user_input != "y":
                raise exceptions.ExecutionStoppedByUser()
            else:
                os.remove(f"{output_folder}/{output_file}_fit.pickle")
        elif os.path.exists(f"{output_folder}/{output_file}_fit.pickle") and remote:
            logger.warning(
                f"{Fore.RED}WARNING: The file {output_file}_fit.pickle already existed. In order to execute the new analysis it has been eliminated."
            )
            os.remove(f"{output_folder}/{output_file}_fit.pickle")

    if significance:
        # Check if the output folder exists
        if (
            os.path.exists(f"{output_folder}/{output_file}_significance.pickle")
            and not remote
        ):
            # Folder exists, ask user if they want to continue
            user_input = (
                input(
                    f"{Fore.RED}WARNING: The file {output_file}_significance.pickle already exists. Do you want to continue and delete overwrite the results? (y/n): {Style.RESET_ALL}"
                )
                .strip()
                .lower()
            )

            if user_input != "y":
                raise exceptions.ExecutionStoppedByUser()
            else:
                os.remove(f"{output_folder}/{output_file}_significance.pickle")
        elif (
            os.path.exists(f"{output_folder}/{output_file}_significance.pickle")
            and remote
        ):
            logger.warning(
                f"{Fore.RED}WARNING: The file {output_file}_significance.pickle already existed. In order to execute the new analysis it has been eliminated."
            )
            os.remove(f"{output_folder}/{output_file}_significance.pickle")

    if mcmc:
        # Check if the output folder exists
        if os.path.exists(f"{output_folder}/{output_file}_mcmc.pickle") and not remote:
            # Folder exists, ask user if they want to continue
            user_input = (
                input(
                    f"{Fore.RED}WARNING: The file {output_file}_mcmc.pickle already exists. Do you want to continue and delete overwrite the results? (y/n): {Style.RESET_ALL}"
                )
                .strip()
                .lower()
            )

            if user_input != "y":
                raise exceptions.ExecutionStoppedByUser()
            else:
                os.remove(f"{output_folder}/{output_file}_mcmc.pickle")
        elif os.path.exists(f"{output_folder}/{output_file}_mcmc.pickle") and remote:
            logger.warning(
                f"{Fore.RED}WARNING: The file {output_file}_mcmc.pickle already existed. In order to execute the new analysis it has been eliminated."
            )
            os.remove(f"{output_folder}/{output_file}_mcmc.pickle")

    if weights:
        # Check if the output folder exists
        if (
            os.path.exists(f"{output_folder}/{output_file}_weights.pickle")
            and not remote
        ):
            # Folder exists, ask user if they want to continue
            user_input = (
                input(
                    f"{Fore.RED}WARNING: The file {output_file}_weights.pickle already exists. Do you want to continue and delete overwrite the results? (y/n): {Style.RESET_ALL}"
                )
                .strip()
                .lower()
            )

            if user_input != "y":
                raise exceptions.ExecutionStoppedByUser()
            else:
                os.remove(f"{output_folder}/{output_file}_weights.pickle")
        elif (
            os.path.exists(f"{output_folder}/{output_file}_weights.pickle")
            and remote
        ):
            logger.warning(
                f"{Fore.RED}WARNING: The file {output_file}_weights.pickle already existed. In order to execute the new analysis it has been eliminated."
            )
            os.remove(f"{output_folder}/{output_file}_weights.pickle")

    return config_analysis


def _astrogeofit_main_function(
    config_file_path: str = "default",
    execute_all: bool = True,
    fit: bool = False,
    significance: bool = False,
    mcmc: bool = False,
    weights: bool = False,
    path_fit_results: Optional[list] = None,
    path_mcmc_results: Optional[str] = None,
    savefigs: bool = False,
    logs: bool = False,
    remote: bool = False,
    test: bool = False,  # ADD TESTS SECTION
):
    global logger
    logger = setup_logger(logs)
    # ---------------- OBTENTION OF THE CONFIGURATION FILE ---------------- #
    if config_file_path == "default":
        config_file_path = f"{os.getcwd()}/config/configuration_file/configuration.yml"

    # ---------------- INTRODUCTION PRINTS + CONFIG FILE ERROR HANDLER ---------------- #
    config_analysis = _info_beginning(
        config_file_path=config_file_path,
        significance=significance,
        fit=fit,
        mcmc=mcmc,
        weights=weights,
        execute_all=execute_all,
        remote=remote,
    )
    if execute_all:
        fit = True
        mcmc = True
        significance = True
        weights = True

    # ---------------- CONFIGURATION FILE CHECK ----------------- #
    variable_check.configuration_file_variables_check(
        config_file=config_analysis,
        significance = significance,
        fit = fit,
        mcmc = mcmc,
        weights = weights,
        remote=remote
    )

    # ---------------- DATA SET RETRIEVAL --------------------- #
    data_set_values = config_analysis["data_set"]
    data_set, _ = data_manager.data_read_dataset_from_file(data_set_values)

    # ---------------- DATA MODEL VALUES RETRIEVAL --------------------- #
    data_model_parameters = config_analysis["data_model_parameters"]
    if fit or mcmc:
        genetic_algorithm_parameters = config_analysis["genetic_algorithm_parameters"]

    # ---------------- OBTENTION OF THE OUTPUT FOLDER --------------------- #
    output_folder = str(config_analysis["output_folder"])
    name_file = str(
        config_analysis["output_file_name"]
    )  # Add this if name_file is part of the config
    output_file_path = os.path.join(output_folder, name_file)
    os.makedirs(name=output_folder, exist_ok=True)

    if execute_all or significance:
        logger.info(f"{Fore.GREEN} EXECUTING SIGNIFICANCE TEST{Style.RESET_ALL}")
        print("")
        significance_test_parameters = config_analysis["significance_test_parameters"]
        significance_test_results = main_significance_test.significance_test_execution(
            data_model_parameters=data_model_parameters,
            significance_test_parameters=significance_test_parameters,
            data_set=data_set,
            remote=remote
        )
        data_manager.save_analysis_info(**config_analysis)
        data_manager.save_data(
            output_file_path=output_file_path,
            fit_mcmc_weights="significance",
            r2_data=significance_test_results["r2_data"],
            r2_ar=significance_test_results["r2_ar"],
        )

    ######################################################################
    ###################### INCLUDE FITTING FUNCTION ######################
    ######################################################################
    if fit or execute_all:
        logger.info(
            f"{Fore.GREEN} EXECUTING GENETIC ALGORITHM FITTING {Style.RESET_ALL}"
        )
        print("")
        # THE PARAMETERS func_invSR_nominal, func_time_nominal AND output_folder ARE SET AS NONE AS THE INTENSIVE ANALYSIS PART IT'S BEING RE-DONE
        results_fitting = main_fitting.main_fitting_function(
            data_model_parameters=data_model_parameters,
            configuration_parameters=genetic_algorithm_parameters,
            data_set=data_set,
            remote = remote
        )

        best_result_fit = results_fitting["optimization_results"][-1][-1]
        invSR_to_predx_standard = partial(
            shared_functions.invSR_to_prediction,
            data=[data_set[0], data_set[1]],
            inverse_SR_lims=best_result_fit.problem.inverse_SR_lims,
            interpolator=best_result_fit.problem.interpolator,
        )

        # ------------- SAVE DATA FIT --------------- #
        data_manager.save_analysis_info(**config_analysis)
        data_manager.save_data(
            output_file_path=output_file_path,
            fit_mcmc_weights="fit",
            optimization_results=results_fitting["optimization_results"],
            inverse_SR_lims=results_fitting["inverse_SR_lims"],
            metric_sorted_resuts=results_fitting["metric_sorted_resuts"],
            convers_factor=results_fitting["convers_factor"],
            best_result_fit=best_result_fit,
            invSR_to_predx=invSR_to_predx_standard,
        )
        logger.info(
            f"{Fore.GREEN}We saved the data of the fitting process. {Style.RESET_ALL}"
        )
        print("")

    # -------- IF MCMC OR WEIGHTS WITHOUT FIT ------------ #
    elif mcmc or weights:
        ### LOAD FIT DATA
        if not path_fit_results:
            path_fit_results = f"{output_file_path}_fit.pickle"
        if not os.path.exists(path_fit_results):
            raise Exception(
                f"The file of the fitting results, does not exists. Please check the path."
            )
        else:
            results_fitting = data_manager.load_dictionary_from_pickle(path_fit_results)

    # ------------------------ MCMC PART -------------------------- #
    if mcmc or execute_all:
        mcmc_parameters = config_analysis["mcmc_parameters"]
        frequency_values = config_analysis["data_model_parameters"]["frequency_values"]
        prior_frequencies, prior_distributions = utils_rmcmc.calculate_prior_frequencies_and_distributions(
            prior_frquencies_values = mcmc_parameters["prior_frequencies"],
            prior_distributions_values = mcmc_parameters["prior_distributions"],
            frequency_values = frequency_values,
            convers_factor = results_fitting["convers_factor"]
        )
        logger.info(f"{Fore.GREEN} EXECUTING MCMC CALCULA {Style.RESET_ALL}")
        print("")

        ###################################################################################
        # WE SHOULD CREATE A SECTION AT THE BEGINNING THAT CHECKS ERRORS IN THE VARIABLES #
        if int(mcmc_parameters["number_mcmc_solutions"]) > (
            genetic_algorithm_parameters["number_algorithm_solutions"]
        ):
            raise Exception(
                "The number of MCMC solutions can not be higher thant the number of solutions of the algorithm. Please change the parameters."
            )
        sampler_par = mcmc_calcula.execution_mcmc(
            optimization_results=results_fitting["optimization_results"],
            prior_frequencies=prior_frequencies,
            prior_distributions=prior_distributions,
            num_solutions=mcmc_parameters["number_mcmc_solutions"],
            list_mcmc_genes=mcmc_parameters["list_of_genes_mcmc"],
            list_num_genes=genetic_algorithm_parameters["list_number_genes"],
            data=data_set,
            num_parallel_jobs=mcmc_parameters["number_processors_used_mcmc"],
            n_step=mcmc_parameters["length_mcmc_chains"],
            seed=genetic_algorithm_parameters["seed"],
            remote = remote
        )
        sampler_MCMC, logprob_par = utils_rmcmc.get_mcmc_samplers_and_logprob(
            sampler_par, mcmc_parameters["discard"], mcmc_parameters["thin"]
        )
        # ----------------------------- DATA SAVING MCMC ------------------------ #
        data_manager.save_data(
            output_file_path=output_file_path,
            fit_mcmc_weights="mcmc",
            sampler_MCMC=sampler_MCMC,  # [samples_combine[1:] for samples_combine in samples_combine_par], # WE HAVE CHANGED THESE TWO LINES IN ORDER TO BE ABLE TO USE 1 SINGLE SOLUTION
            logprob_MCMC=logprob_par,  # [logprob[1:] for logprob in logprob_par],              # I SHOULD ASK WHY THIS WAS AS BEFORE
            prior_distributions=prior_distributions,
            prior_frequencies=prior_frequencies,
            executed_weights=False,
        )
    # OBTENTION OF THE MCMC RESULTS IN CASE THE THAT THE MCMC IS NOT EXECUTED
    elif weights:
        if not path_mcmc_results:
            path_mcmc_results = f"{output_file_path}_mcmc.pickle"
        if not os.path.exists(path_mcmc_results):
            raise Exception(
                f"The file of the mcmc results , does not exists. Please check the path."
            )
        else:
            results_mcmc = data_manager.load_dictionary_from_pickle(path_mcmc_results)
            sampler_MCMC = results_mcmc["sampler_MCMC"]
            prior_frequencies = results_mcmc["prior_frequencies"]
            prior_distributions = results_mcmc["prior_distributions"]

    #####################################################################################################
    if execute_all or weights:
        logger.info(f"{Fore.GREEN} EXECUTING MCMC WEIGHT CALCULA {Style.RESET_ALL}")
        print("")
        weights_parameters = config_analysis["weight_calcula_configuration"]
        evidences = weight_computation.execute_weight_computation(
            optimization_results=results_fitting["optimization_results"],
            data=data_set,
            samples_MCMC=sampler_MCMC,
            num_weight_executions=weights_parameters[
                "number_weight_evaluation_per_chain"
            ],
            num_parallel_jobs_weights=weights_parameters[
                "number_processors_used_weights"
            ],
            prior_frequencies=prior_frequencies,
            prior_distributions=prior_distributions,
            stability_factor=weights_parameters["stability_factor"],
            pareto_smoothing=weights_parameters["pareto_smoothing"],
            seed = config_analysis["genetic_algorithm_parameters"]["seed"]
        )

        data_manager.save_data(
            output_file_path=output_file_path,
            fit_mcmc_weights="mcmc",
            executed_weights=True,
        )

        data_manager.save_data(
            output_file_path=output_file_path,
            fit_mcmc_weights="weights",
            evidences=evidences,
        )

    if savefigs:
        file_types_for_saving = ["jpeg", "pdf"]
        output_path = f"{output_folder}/plots"
        os.makedirs(os.path.dirname(f"{output_path}/"), exist_ok=True)
        if fit:
            AGF_Plots_Post_Fit.plt_fit_summary_of_the_results(
                config_file_path, 
                "default", 
                "default",
                types_of_saving=file_types_for_saving,
                output_path=output_path,
                x_axis_limits=None,
            )
            AGF_Plots_Post_Fit.plt_fit_metric_per_number_of_knots(
                config_file_path, 
                "default", 
                "default",
                types_of_saving=file_types_for_saving,
                output_path=output_path,
            )

            AGF_Plots_Post_Fit.plt_fit_aics_bics_per_number_knots(
                config_file_path, 
                "default", 
                "default",
                types_of_saving=file_types_for_saving,
                output_path=output_path,
            )
            AGF_Plots_Post_Fit.plt_fit_aics_r2_per_number_knots(
                config_file_path, 
                "default", 
                "default",
                types_of_saving=file_types_for_saving,
                output_path=output_path,
            )
            AGF_Plots_Post_Fit.plt_fit_sedimentation_rate_per_depth(
                config_file_path, 
                "default", 
                "default",
                y_axis_range=None,
                types_of_saving=file_types_for_saving,
                output_path=output_path,
            )
            AGF_Plots_Post_Fit.plt_fit_sedimentation_rate_depth_with_uncertainty(
                config_file_path, 
                "default", 
                "default",
                number_of_genes_to_explore=genetic_algorithm_parameters[
                    "list_number_genes"
                ][-1],
                types_of_saving=file_types_for_saving,
                output_path=output_path,
            )
            AGF_Plots_Post_Fit.plt_fit_frequency_spectrum(
                config_file_path, 
                "default", 
                "default",
                x_axis_plot_limits=None,
                y_axis_plot_limits=None,
                peridogram_scale_x_axis=None,
                peridogram_scale_y_axis=None,
                types_of_saving=file_types_for_saving,
                output_path=output_path,
            )
            AGF_Plots_Post_Fit.plt_fit_time_series(
                config_file_path, 
                "default", 
                "default",
                None,
                types_of_saving=file_types_for_saving,
                output_path=output_path,
            )
            _, _, _, _, _, _, _, _, _, eccentricity_solution_parameters = AGF_Notebooks.data_obtain_analysis_file_and_result_files(False, config_file_path, "default", "default")
            if isinstance(eccentricity_solution_parameters, dict):
                AGF_Plots_Post_Fit.data_calculate_eccentricity_parameters(
                    config_file_path, 
                    "default", 
                    "default",
                    positive_feedback=True,
                )
                AGF_Plots_Post_Fit.plt_fit_ETP_signals_from_model(
                    config_file_path, 
                    "default", 
                    "default",
                    types_of_saving=file_types_for_saving,
                    output_path=output_path,
                )
                AGF_Plots_Post_Fit.plt_fit_eccentricity_mswd_plot(
                    config_file_path, 
                    "default", 
                    "default",
                    time_range=None,
                    types_of_saving=file_types_for_saving,
                    output_path=output_path,
                )
                AGF_Plots_Post_Fit.plt_fit_correlation_eccentricity_and_solution(
                    config_file_path, 
                    "default", 
                    "default",
                    x_axis_limits=None,
                    types_of_saving=file_types_for_saving,
                    output_path=output_path,
                )
        if mcmc:    
            AGF_MCMC_Plots.plt_mcmc_SR_per_num_knots(
                config_file_path, 
                "default", 
                "default",
                x_axis_plot_limits=None,
                y_axis_plot_limits=None,
                ignore_weights=False,
                types_of_saving=file_types_for_saving,
                output_path=output_path,
            )
            AGF_MCMC_Plots.plt_mcmc_aic_and_r2_per_number_knots(
                config_file_path, 
                "default", 
                "default",
                types_of_saving=file_types_for_saving,
                output_path=output_path,
            )
            AGF_MCMC_Plots.plt_mcmc_aic_logprob_and_loglike_per_number_knots(
                config_file_path, 
                "default", 
                "default",
                types_of_saving=file_types_for_saving,
                output_path=output_path,
            )
            AGF_Shared_Functions.data_obtain_mcmc_results_per_a_number_of_knots(
                config_file_path, 
                "default", 
                "default",
                positive_feedback=True,
                use_prec_env=False,
                ignore_weights=False
            )
            ##########################
            ## THIS MIGHT NOT WORK ###
            ##########################
            AGF_MCMC_Plots.plt_mcmc_eccentricity_correlation_with_solution(
                config_file_path, 
                "default", 
                "default",
                t0_offset_range_plot=None,
                ecc_time_range=None,
                save_eccentricity_solutions = False,
                save_eccentricity_solutions_comment = "",
                name_eccentricity_solution="Eccentricity Solution",
                types_of_saving=file_types_for_saving,
                output_path=output_path,
            )
            AGF_MCMC_Plots.plt_mcmc_prior_frequencies_distributions(
                config_file_path, 
                "default", 
                "default",
                types_of_saving=file_types_for_saving,
                output_path=output_path,
            )
            AGF_MCMC_Plots.plt_mcmc_summary_of_results(
                config_file_path, 
                "default", 
                "default",
                x_axis_limits_peridogram=None,
                y_axis_limits_peridogram=None,
                peridogram_scale_x_axis="linear",
                peridogram_scale_y_axis="linear",
                x_axis_limits_depth=None,
                types_of_saving=file_types_for_saving,
                output_path=output_path,
            )
            AGF_MCMC_Plots.plt_mcmc_phase_of_frequencies(
                config_file_path, 
                "default", 
                "default", 
                types_of_saving=file_types_for_saving, 
                output_path=output_path
            )
        if significance:
            AGF_Plot_Significance_Test.plt_significance_test_results(
                config_file_path,
                "default", 
                "default", 
                "Custom Section",
                file_types_for_saving=file_types_for_saving, 
                output_path=output_path)


class AstroGeoFit_tool:
    """
    A high-level interface for executing AstroGeoFit model processes,
    including fitting, significance testing, MCMC analysis, and weighting.

    This tool wraps the core functionality provided by the internal
    `run` which  allows users to run different
    parts of the AstroGeoFit workflow based on configuration.

    Results from each executed process are saved in separate pickle files:
        - {name_given_in_the_config_file}_fit_results.pickle
        - {name_given_in_the_config_file}_significance_results.pickle
        - {name_given_in_the_config_file}_mcmc_results.pickle
        - {name_given_in_the_config_file}_weight_results.pickle

    Args:
        - **config_file_path** (str): Path to the configuration file.
        - **execute_all** (bool, optional): Whether to execute all available processes. Defaults to True.
        - **fit** (bool, optional): Execute the fitting process. Defaults to False.
        - **significance** (bool, optional): Execute the significance test. Defaults to False.
        - **mcmc** (bool, optional): Run the MCMC chain analysis. Defaults to False.
        - **weight** (bool, optional): Compute or apply weighting. Defaults to False.
        - **path_fit_results** (str, optional): Path to output fit results. Only required if the execution calculates the mcmc or the weights but not the fitting part.
        - **path_mcmc_results** (str, optional): Path to output MCMC results. Only required if the execution calculates the weights but not the mcmc part.
        - **savefigs** (bool, optional): Whether to save figures. Defaults to False.
        - **logs** (bool, optional): Whether to generate log files. Defaults to False.
        - **remote** (bool, optional): Set to True if executing remotely. This makes Defaults to False.
        - **test** (bool, optional): If True, run in test mode (dry run or minimal config). Defaults to False.

    Methods:
        run(): Executes the requested AstroGeoFit processes using the internal main function.
    """

    def __init__(
        self,
        config_file_path="default",
        execute_all: Optional[bool] = False,
        fit: Optional[bool] = False,
        significance: Optional[bool] = False,
        mcmc: Optional[bool] = False,
        weight: Optional[bool] = False,
        path_fit_results=None,
        path_mcmc_results=None,
        savefigs=False,
        logs=False,
        remote=False,
        test=False,
    ):
        self.config_file_path = config_file_path
        self.execute_all = execute_all
        self.fit = fit
        self.significance = significance
        self.mcmc = mcmc
        self.weight = weight
        self.path_fit_results = path_fit_results
        self.path_mcmc_results = path_mcmc_results
        self.savefigs = savefigs
        self.logs = logs
        self.remote = remote
        self.test = test

    def run(self):
        """
        Executes the requested AstroGeoFit processes using the internal main function.

        This method acts as a wrapper around the `_astrogeofit_main_function`, passing along all
        relevant configuration parameters required for different AstroGeoFit execution modes such as 
        fitting, significance testing, MCMC simulation, and plotting.

        The parameters used to execute the main astrogeofit function are the ones set when creating the class.
        
        Returns:
            None
        """

        _astrogeofit_main_function(
            config_file_path=self.config_file_path,
            execute_all=self.execute_all,
            fit=self.fit,
            significance=self.significance,
            mcmc=self.mcmc,
            weights=self.weight,
            path_fit_results=self.path_fit_results,
            path_mcmc_results=self.path_mcmc_results,
            savefigs=self.savefigs,
            logs=self.logs,
            remote=self.remote,
            test=self.test,
        )
