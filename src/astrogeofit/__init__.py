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

from .astrogeofit import AstroGeoFit_tool

from .utils.setup_logger import setup_logger

from .main_routines.data_manager import data_read_dataset_from_file

from .notebooks_code.notebooks_code import data_obtain_age_depth_values, data_obtain_analysis_file_and_result_files

from .plots.plots_post_fitting.plots_post_fit import (
    plt_fit_summary_of_the_results,
    plt_fit_metric_per_number_of_knots,
    plt_fit_aics_bics_per_number_knots,
    plt_fit_aics_r2_per_number_knots,
    plt_fit_sedimentation_rate_per_depth,
    plt_fit_sedimentation_rate_depth_with_uncertainty,
    plt_fit_frequency_spectrum,
    plt_fit_time_series,
    data_calculate_eccentricity_parameters,
    plt_fit_ETP_signals_from_model,
    plt_fit_eccentricity_mswd_plot,
    plt_fit_correlation_eccentricity_and_solution,
    tool_plt_fit_summary_of_the_results,
    tool_plt_fit_metric_per_number_of_knots,
    tool_plt_fit_aics_bics_per_number_knots,
    tool_plt_fit_aics_r2_per_number_knots,
    tool_plt_fit_sedimentation_rate_per_depth,
    tool_plt_fit_sedimentation_rate_depth_with_uncertainty,
    tool_plt_fit_frequency_spectrum,
    tool_plt_fit_time_series,
    tool_data_calculate_eccentricity_parameters,
    tool_plt_fit_ETP_signals_from_model,
    tool_plt_fit_eccentricity_mswd_plot,
    tool_plt_fit_correlation_eccentricity_and_solution
)
from .plots.plots_post_mcmc.plots_post_mcmc import (
    plt_mcmc_SR_per_num_knots,
    tool_plt_mcmc_SR_per_num_knots,
    plt_mcmc_aic_and_r2_per_number_knots,
    tool_plt_mcmc_aic_and_r2_per_number_knots,
    plt_mcmc_eccentricity_correlation_with_solution,
    tool_plt_mcmc_eccentricity_correlation_with_solution,
    plt_mcmc_summary_of_results,
    tool_plt_mcmc_summary_of_results,
    plt_mcmc_prior_frequencies_distributions,
    tool_plt_mcmc_prior_frequencies_distributions,
    plt_mcmc_phase_of_frequencies,
    tool_plt_mcmc_phase_of_frequencies,
    tool_plt_mcmc_aic_logprob_and_loglike_per_number_knots,
    plt_mcmc_aic_logprob_and_loglike_per_number_knots
)

from .plots.significance_test_plot.significance_test_plot import (
    plt_significance_test_results,
    tool_plt_significance_test_results
)

from .plots.utils.shared_utils_plots import (
    data_obtain_mcmc_results_per_a_number_of_knots,
    tool_data_obtain_mcmc_results_per_a_number_of_knots
)

__all__ = [
    "AstroGeoFit_tool",
    "plt_fit_summary_of_the_results",
    "plt_fit_metric_per_number_of_knots",
    "plt_fit_aics_bics_per_number_knots",
    "plt_fit_aics_r2_per_number_knots",
    "plt_fit_sedimentation_rate_per_depth",
    "plt_fit_sedimentation_rate_depth_with_uncertainty",
    "plt_fit_frequency_spectrum",
    "plt_fit_time_series",
    "data_calculate_eccentricity_parameters",
    "plt_fit_ETP_signals_from_model",
    "plt_fit_eccentricity_mswd_plot",
    "plt_fit_correlation_eccentricity_and_solution",
    "plt_mcmc_SR_per_num_knots",
    "plt_mcmc_aic_and_r2_per_number_knots",
    "plt_mcmc_eccentricity_correlation_with_solution",
    "plt_mcmc_summary_of_results",
    "plt_mcmc_prior_frequencies_distributions",
    "plt_mcmc_phase_of_frequencies",
    "data_obtain_mcmc_results_per_a_number_of_knots",
    "plt_significance_test_results",
    "setup_logger",
    "tool_plt_fit_summary_of_the_results",
    "tool_plt_fit_metric_per_number_of_knots",
    "tool_plt_fit_aics_bics_per_number_knots",
    "tool_plt_fit_aics_r2_per_number_knots",
    "tool_plt_fit_sedimentation_rate_per_depth",
    "tool_plt_fit_sedimentation_rate_depth_with_uncertainty",
    "tool_plt_fit_frequency_spectrum",
    "tool_plt_fit_time_series",
    "tool_data_calculate_eccentricity_parameters",
    "tool_plt_fit_ETP_signals_from_model",
    "tool_plt_fit_eccentricity_mswd_plot",
    "tool_plt_fit_correlation_eccentricity_and_solution",
    "tool_plt_mcmc_SR_per_num_knots",
    "tool_plt_mcmc_aic_and_r2_per_number_knots",
    "tool_plt_mcmc_eccentricity_correlation_with_solution",
    "tool_plt_mcmc_summary_of_results",
    "tool_plt_mcmc_prior_frequencies_distributions",
    "tool_plt_mcmc_phase_of_frequencies",
    "tool_data_obtain_mcmc_results_per_a_number_of_knots",
    "tool_plt_significance_test_results",
    "data_read_dataset_from_file",
    "data_obtain_age_depth_values",
    "data_obtain_analysis_file_and_result_files",
    "tool_plt_mcmc_aic_logprob_and_loglike_per_number_knots",
    "plt_mcmc_aic_logprob_and_loglike_per_number_knots"
]
