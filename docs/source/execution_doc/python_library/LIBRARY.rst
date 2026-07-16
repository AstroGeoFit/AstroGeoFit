.. _library_documentation:

AstroGeoFit Python Library
==========================

.. toctree::
    :maxdepth: 2
    :hidden:
    
    asterogeofit_class
    plt_fit_functions
    plt_mcmc_functions
    plt_significance_functions
    data_functions


.. image:: https://img.shields.io/badge/python-3.11%2B-blue
    :alt: Python Version

.. image:: https://img.shields.io/badge/license-GPLv3-blue
    :alt: License

Execution of the Library
------------------------

The **AstroGeoFit** library provides a wide range of functions designed to help users both obtain and visualize the results produced by the tool. 
To improve readability and usability, these functions are organized by functionality through the use of specific prefixes that indicate their purpose.

AstroGeoFit_tool class
^^^^^^^^^^^^^^^^^^^^^^

This class represents the main core of the tool, which will recieve all the necessary parameters
to obtain the results of the different sections (fitting, MCMC, Weight Computation and Significance Test).
This class contains the function `run()` which will execute the tool itself. Find the detailed guide in: :ref:`astrogeofit_class`.


*plt_fit_* functions
^^^^^^^^^^^^^^^^^^^^^

The functions with this prefix generate plots from the results obtained after fitting our genetic algorithm and predicting
our new data with it. Find the detailed guide in: :ref:`plt_fit_functions`.

The functions that we find on this category are:

    - `plt_fit_summary_of_the_results`
    - `plt_fit_metric_per_number_of_knots`
    - `plt_fit_aics_bics_per_number_knots`
    - `plt_fit_aics_r2_per_number_knots`
    - `plt_fit_sedimentation_rate_per_depth`
    - `plt_fit_sedimentation_rate_depth_with_uncertainty`
    - `plt_fit_spectrum_plot`
    - `plt_fit_time_series`
    - `plt_fit_ETP_signals_from_model`
    - `plt_fit_eccentricity_mswd_plot`
    - `plt_fit_correlation_eccentricity_and_solution`

**Note**: Some of these functions require another function from the *data_* category to be executed.

*plt_significance_* function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The function with this prefix generates a plot using the results of the significance test.
Find the detailed guide in: :ref:`plt_significance_functions`.

The function that we find on this category is:

    - `plt_significance_test_results`

*plt_mcmc_* functions
^^^^^^^^^^^^^^^^^^^^^

The functions that have this prefix use the results of the mcmc and the weight calcula (if executed) to obtain different types of plots.
Find the detailed guide in: :ref:`plt_mcmc_functions`.

The functions that we find on this category are:

    - `plt_mcmc_SR_per_num_knots`
    - `plt_mcmc_aic_and_r2_per_number_knots`
    - `plt_mcmc_eccentricity_correlation_with_solution`
    - `plt_mcmc_summary_of_results`
    - `plt_mcmc_prior_frequencies_distributions`
    - `plt_mcmc_phase_of_frequencies`

**Note**: Some of these functions require another function from the *data_* category to be executed.

*data_* functions
^^^^^^^^^^^^^^^^^^

Functions with this prefix are related to obtaining data required for specific plots. Plots that depend on a preceding function will be 
marked with the **Pre-Requisite**. Find the detailed guide in: :ref:`data_functions`.

The functions that we find on this category are:

    - `data_calculate_eccentricity_parameters`
    - `data_obtain_mcmc_results_per_a_number_of_knots`
    - `data_obtain_analysis_file_and_result_files`
    - `data_obtain_age_depth_values`

Configuration Handling
----------------------

The library functions can be run independently or configured through the `.yml` configuration file parser provided.

We recommend using the configuration loader to prepare the configuration file before calling plotting or execution functions.

Example Usage
-------------
This example shows how to obtain the results of the fitting part, and plot the summary of those results. The path of the fitting results file is taken from the configuration file by not introducing the parameters 
`output_folder` and `output_file` on the `data_obtain_analysis_file_and_result_files` function.

.. code-block:: python

    from astrogeofit import plt_fit_summary_of_the_results

    configuration_file_path = "path/to/your/configuration_file.yaml"
    number_of_knots_to_explore = 54 
    types_of_saving = ["pdf"]
    output_path = "output/results/path"
    x_axis_limits = None
    invSR_plot_y_limits = [y_min,y_max]
    age_plot_y_limits = [y_min,y_max]
    fitted_data_plot_y_limits = [y_min,y_max]
    correlation_plot_y_limits = [y_min,y_max]
    # Load optimization results, prediction function, and data...
    plt_fit_summary_of_the_results(
        configuration_file_path = configuration_file_path,
        number_of_knots_to_explore = number_of_knots_to_explore,
        types_of_saving = types_of_saving,
        output_path = output_path,
        x_axis_limits = x_axis_limits,
        invSR_plot_y_limits = invSR_plot_y_limits
        age_plot_y_limits = age_plot_y_limits
        fitted_data_plot_y_limits = fitted_data_plot_y_limits
        correlation_plot_y_limits = correlation_plot_y_limits
    )