.. _plt_fit_functions:

*plt_fit_* detailed function guide
==================================

.. image:: https://img.shields.io/badge/python-3.11%2B-blue
    :alt: Python Version

.. image:: https://img.shields.io/badge/license-GPLv3-blue
    :alt: License
    
.. _plt_fit_summary_of_the_results:

.. function:: astrogeofit.plt_fit_summary_of_the_results

.. code-block:: python

    plt_fit_summary_of_the_results(
        configuration_file_path:str,
        folder_with_the_results:Optional[str],
        file_with_the_results:Optional[str],
        number_of_knots_to_explore: Optional[int],
        types_of_saving: list[str],
        output_path: Optional[str],
        x_axis_limits: Optional[list[float]],
        invSR_plot_y_limits: Optional[list[float]],
        age_plot_y_limits: Optional[list[float]],
        fitted_data_plot_y_limits: Optional[list[float]],
        correlation_plot_y_limits: Optional[list[float]],
    )

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
    - `configuration_file_path`: str 
        Path to the configuration `.yml` file used in the analysis.

    - `folder_with_the_results`: str or None 
        Folder containing the output results. If None, uses default from config.

    - `file_with_the_results`: str or None 
        File name of the serialized results. If None, uses default from config.

    - `number_of_knots_to_explore`: int or None 
        Number of knots (genes) to visualize from the list used during optimization. If None, the last available entry is used.

    - `types_of_saving` : list of str or None 
        List of file formats to save outputs (e.g., ["png", "pdf", "tab"]).

    - `output_path`: str or None 
        Path to the output directory for saving figures and `.tab` files. If None, uses the default output folder from the configuration.
    
    - `x_axis_limits`: list of floats or None 
        X-axis limits for all subplots (in depth units).

    - `invSR_plot_y_limits`: list of floats or None 
        Y-axis limits for the inverse sedimentation rate subplot.

    - `age_plot_y_limits`: list of floats or None 
        Y-axis limits for the age-depth (time-depth) subplot.

    - `fitted_data_plot_y_limits`: list of floats or None 
        Y-axis limits for the fitted signal subplot.

    - `correlation_plot_y_limits`: list of floats or None 
        Y-axis limits for the rolling correlation subplot.
        
----

.. _plt_fit_metric_per_number_of_knots:

.. function:: astrogeofit.plt_fit_metric_per_number_of_knots

.. code-block:: python

    plt_fit_metric_per_number_of_knots(
        configuration_file_path:str,
        folder_with_the_results:Optional[str],
        file_with_the_results:Optional[str],
        types_of_saving: list,
        output_path: Optional[str],
    )

Plots the evolution of a performance metric across generations and knot configurations.

This function visualizes the chosen metric (e.g., r²) over the number of generations
for different solutions obtained during the optimization process. It also adds vertical
lines and labels corresponding to different knot (gene) configurations used during fitting.

Parameters:
    - `configuration_file_path`: str 
        Path to the configuration `.yml` file used for analysis.

    - `folder_with_the_results`: str or None 
        Folder containing the output results from the optimization process.

    - `file_with_the_results`: str or None 
        File containing the results to load from.

    - `types_of_saving` : list of str or None 
        List of file formats to save the plot in (e.g., ["png", "pdf"]).

    - `output_path`: str or None 
        Directory where the output figure will be saved. If not provided, the default output folder from the configuration will be used.

----

.. _plt_fit_aics_bics_per_number_knots:

.. function:: astrogeofit.plt_fit_aics_bics_per_number_knots

.. code-block:: python

    plt_fit_aics_bics_per_number_knots(
        configuration_file_path:str,
        folder_with_the_results:Optional[str],
        file_with_the_results:Optional[str],
        types_of_saving: list[str],
        output_path: Optional[str],
    )

Plot AIC and BIC values across different knot configurations and algorithm solutions.

This function computes and visualizes the Akaike Information Criterion (AIC) and Bayesian 
Information Criterion (BIC) for each solution obtained during the optimization process. 
For each number of knots (genes), the function fits an autoregressive model to the residuals 
between the fitted and original data, and computes the corresponding AIC and BIC values. 
The final figure contains two subplots, one for AIC and one for BIC, with color-coded lines 
representing each solution and a black line for the minimum value per knot configuration.

Parameters:
    - `configuration_file_path`: str  
        Path to the `.yml` configuration file used during the analysis.

    - `folder_with_the_results`: str or None 
        Folder containing the results of the optimization.

    - `file_with_the_results`: str or None 
        File from which the results should be loaded.

    - `types_of_saving` : list of str or None 
        List of output file formats to save the figure in (e.g., ["png", "pdf"]).

    - `output_path`: str or None 
        Directory to save the figure. If None, the path from the configuration file is used under the `figures` subdirectory.

----

.. _plt_fit_aics_r2_per_number_knots:

.. function:: astrogeofit.plt_fit_aics_r2_per_number_knots

.. code-block:: python

    plt_fit_aics_r2_per_number_knots(
        configuration_file_path:str,
        folder_with_the_results:Optional[str],
        file_with_the_results:Optional[str],
        types_of_saving: list,
        output_path: str,
    )

Plot AIC and r² values for a single optimization solution across different knot configurations.

This function generates a dual-axis plot comparing the Akaike Information Criterion (AIC) 
and the coefficient of determination (r²) across varying numbers of knots (genes) for the 
final optimization result. The AIC is calculated based on fitting an autoregressive model 
to the residuals between the predicted and original data. The r² score evaluates the 
goodness of fit of the model predictions.

Parameters:
    - `configuration_file_path`: str 
        Path to the `.yml` configuration file used for the analysis.

    - `folder_with_the_results`: str or None 
        Folder containing the optimization results.

    - `file_with_the_results`: str or None 
        File from which results should be loaded.

    - `types_of_saving` :  list of str or None 
        List of file formats to save the figure in (e.g., ["png", "pdf"]).

    - `output_path` : str or None 
        Directory to save the plot. If not provided or set to "default", it falls back to the path defined in the configuration under `output_folder/figures`.

----

.. _plt_fit_sedimentation_rate_per_depth:

.. function:: astrogeofit.plt_fit_sedimentation_rate_per_depth

.. code-block:: python

    plt_fit_sedimentation_rate_per_depth(
        configuration_file_path:str,
        folder_with_the_results:Optional[str],
        file_with_the_results:Optional[str],
        y_axis_range: list,
        types_of_saving: list,
        output_path: Optional[str],
    )

Plot inverse sedimentation rate (SR) as a function of depth for multiple knot configurations.

This function generates subplots showing the inverse sedimentation rate profiles for each 
number of knots (genes) used during the optimization. The plot helps to visualize how the 
sedimentation model changes depending on the resolution (number of knots) and across multiple 
optimization results. Optionally, a nominal sedimentation rate profile is overlaid for comparison.

Parameters:
    - `configuration_file_path`: str 
        Path to the `.yml` configuration file used for the analysis.

    - `folder_with_the_results`: str or None 
        Folder containing the result files.

    - `file_with_the_results`: str or None 
        Specific file name with the results, if applicable.

    - `y_axis_range`: list or None 
        Optional range for the y-axis (e.g., [min, max]) to constrain plots vertically.

    - `types_of_saving` : list of str or None 
        List of file formats to save the figure in (e.g., ["png", "pdf"]).

    - `output_path`: str or None 
        Directory to save output plots. If set to "default" or left empty, uses the configured output directory from the analysis configuration.

----

.. _plt_fit_sedimentation_rate_depth_with_uncertainty:

.. function:: astrogeofit.plt_fit_sedimentation_rate_depth_with_uncertainty

.. code-block:: python

    plt_fit_sedimentation_rate_depth_with_uncertainty(
        configuration_file_path:str,
        folder_with_the_results:Optional[str],
        file_with_the_results:Optional[str],
        number_of_genes_to_explore: Optional[int],
        types_of_saving: list,
        output_path: str,
    )

Plot the inverse sedimentation rate (SR) with uncertainty bands for a selected number of knots.

This function computes the median inverse sedimentation rate and associated uncertainty intervals 
(50% and 90% confidence intervals) across multiple optimization results for a given number of knots. 
It visualizes these estimates as a function of depth, and optionally overlays the nominal sedimentation 
rate for comparison. The resulting plot and optional tabular data can be saved in user-specified formats.

Parameters:
    - `configuration_file_path`: str 
        Path to the `.yml` configuration file used for the analysis.

    - `folder_with_the_results`: str or None 
        Folder containing result files.

    - `file_with_the_results`: str or None 
        Name of the results file, if applicable.

    - `number_of_genes_to_explore`: int or None 
        Number of knots (genes) to visualize the sedimentation rate for. If not specified, the last configuration in the results will be used.
    
    - `types_of_saving`: list of str or None 
        List of file types to save the figure. Acceptable formats include image formats (e.g., "png", "pdf") and "tab" for tab-separated data export.
    
    - `output_path`: str or None 
        Path to the directory where the outputs will be saved. If set to "default" or left empty, the output path is taken from the configuration.

----

.. _plt_fit_spectrum_plot:

.. function:: astrogeofit.plt_fit_spectrum_plot

.. code-block:: python

    plt_fit_spectrum_plot(
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
    )

Plot the power spectrum of the time series inferred from sedimentation rate fitting.

This function estimates a time series from sedimentation rate fitting results and performs
a Fast Fourier Transform (FFT) to compute its power spectrum. The resulting spectrum is
plotted as a function of frequency (cycles/Myr), with an optional top x-axis showing the
corresponding periods (Myr). The user can customize axis scaling, plot limits, and file 
formats for saving the figure.

Parameters:
    - `configuration_file_path`: str 
        Path to the `.ymal` configuration file used for the fitting.

    - `folder_with_the_results`: str or None 
        Folder containing result files.

    - `file_with_the_results`: str or None 
        Name of the result file, if applicable.

    - `x_axis_plot_limits`: list of int/float or None 
        List of two floats [min, max] to limit the x-axis (frequency).

    - `y_axis_plot_limits`: list of int/float or None 
        List of two floats [min, max] to limit the y-axis (amplitude).

    - `peridogram_scale_x_axis`: str or None
        Scale type for the x-axis; "linear" (default) or "log". If nothing is introduced, "linear" is selected.

    - `peridogram_scale_y_axis`: str or None 
        Scale type for the y-axis; "linear" (default) or "log". If nothing is introduced, "linear" is selected.

    - `number_of_knots_to_explore`: int or None 
        Number of knots (genes) to select a specific fitting result. If not specified, the last configuration will be used.

    - `types_of_saving`: list of str or None 
        List of file types to save the figure (e.g., "png", "pdf"). The "tab" format is not supported for this plot and will raise a warning.
    
    - `output_path`: str or None 
        Output folder where figures will be saved. If set to "default" or empty, the path is obtained from the configuration file.

----

.. _plt_fit_time_series:

.. function:: astrogeofit.plt_fit_time_series

.. code-block:: python

    plt_fit_time_series(
        configuration_file_path:str,
        folder_with_the_results:Optional[str],
        file_with_the_results:Optional[str],
        types_of_saving: list,
        output_path: str,
    )

Plot depth and time series of proxy data and model predictions from a genetic algorithm fitting,
and optionally save the resulting figures and/or data in various formats.

This function reads the dataset and optimization results from the AstroGeoFit analysis,
retrieves the selected result based on the number of knots to explore, and then generates
plots comparing the raw data and the model predictions in both the depth and time domains.

The plots can be saved in several formats including EPS, PNG, PDF, SVG, and table format (.tab).

Parameters:
    - `configuration_file_path`: str 
        Path to the configuration file used in the analysis.

    - `folder_with_the_results`: str or None 
        Path to the folder containing the fitting results.

    - `file_with_the_results`: str or None 
        Name of the file with the fitting results.

    - `types_of_saving`: list of str or None 
        List of formats to save the figures and/or data. Supported formats: "png", "pdf", "svg", "eps", "tab".
    
    - `output_path`: str or None
        Path where the figures/tables should be saved. If empty or "default", it uses the default output folder from the analysis parameters.

----

.. _plt_fit_ETP_signals_from_model:

.. function:: astrogeofit.plt_fit_ETP_signals_from_model

.. code-block:: python

    plt_fit_ETP_signals_from_model(
        configuration_file_path:str,
        folder_with_the_results:Optional[str],
        file_with_the_results:Optional[str],
        types_of_saving: list,
        output_path: str,
    )

Plots the astronomical components (eccentricity, precession, and tilt) from a previously 
calculated model fit, and optionally saves the resulting figure in multiple formats.

This function loads a set of inferred astronomical signals and their associated time scale, 
and generates a plot showing:

- The fitted precession signal and its envelope.
- The fitted eccentricity signal.
- The fitted tilt signal (if available).
- The correlation (r and r²) between the normalized eccentricity and the normalized precession envelope.

Parameters:
    - `configuration_file_path` : str
        Path to the configuration file used during the model fitting process.

    - `folder_with_the_results` : str or None
        Path to the folder containing the results. If None, uses the path from the configuration file.

    - `file_with_the_results` : str or None
        Name of the specific results file to load. If None, uses the default defined in the configuration.

    - `file_with_the_results` : list
        List of file formats to save the figure in (e.g., ["png", "pdf"]). 
        The format must be one of the supported formats (excluding "tab").

    - `output_path` : str
        Output directory where the figure(s) should be saved. If "default", "", or None, 
        it uses the default output folder from the analysis parameters.

**PRE-REQUISITES**: This function requires the results of `data_calculate_eccentricity_parameters`.

----

.. _plt_fit_eccentricity_mswd_plot:

.. function:: astrogeofit.plt_fit_eccentricity_mswd_plot

.. code-block:: python

    plt_fit_eccentricity_mswd_plot(
        configuration_file_path:str,
        folder_with_the_results:Optional[str],
        file_with_the_results:Optional[str],
        time_range: Optional[list[float]],
        types_of_saving: list,
        output_path: str,
    )

Plots the Mean Squared Weighted Deviation (MSWD) between the fitted eccentricity/precession envelope
and a reference astronomical solution, as a function of t₀ (the starting time of the inferred time scale).

This function evaluates how well the inferred normalized eccentricity and precession envelope align
with the chosen astronomical solution over a range of t₀ values. It produces two subplots: one for
eccentricity and another for the precession envelope. The plots highlight the best-fit t₀ values where
the MSWD is minimized.

Parameters:
    - `configuration_file_path` : str
        Path to the configuration file used during the model fitting process.

    - `folder_with_the_results` : str or None
        Path to the folder containing the result files. If None, the folder is inferred from the configuration.

    - `file_with_the_results` : str or None
        Name of the result file to be loaded. If None, it defaults to the one defined in the configuration.

    - `time_range` : list or None
        A list containing two floats or integers [tmin, tmax] defining the time range to search for optimal t₀.
        If None, the function uses a default time range based on the nominal top age of the dataset.

    - `types_of_saving` : list of str or None
        A list of file formats to save the resulting plot (e.g., ["png", "pdf"]). Format must be in 
        the supported formats list (excluding "tab").

    - `output_path` : str or None
        Directory where the output figure(s) and log file will be saved. If "default", "", or None, the output 
        path is inferred from the configuration file's output folder.

**PRE-REQUISITES**: This function requires the results of `data_calculate_eccentricity_parameters` and `plt_fit_ETP_signals_from_model`.

----

.. _plt_fit_correlation_eccentricity_and_solution:

.. function:: astrogeofit.plt_fit_correlation_eccentricity_and_solution

.. code-block:: python

    plt_fit_correlation_eccentricity_and_solution(
        configuration_file_path:str,
        folder_with_the_results:Optional[str],
        file_with_the_results:Optional[str],
        x_axis_limits: Optional[list],
        types_of_saving: list,
        output_path: str,
    )

Plot the correlation between the fitted eccentricity or precession envelope 
and the astronomical solution used in the AstroGeoFit modeling.

This function reads the analysis configuration, dataset, and results of the 
eccentricity or precession envelope fitting. It then compares these model results 
against the corresponding astronomical solution and plots both the fitted curves 
and the astronomical signal. Correlation metrics (r and r²) are calculated and 
displayed on the plots.

Parameters:
    - `configuration_file_path` : str
        Path to the configuration file used for the AstroGeoFit analysis.

    - `folder_with_the_results` : str or None
        Path to the folder containing the result files. If None, the default folder 
        defined in the configuration file is used.

    - `file_with_the_results` : str or None
        Specific filename for the results. If None, a default filename is assumed.

    - `x_axis_limits` : list or None
        Optional list specifying the [xmin, xmax] limits for the x-axis (in Ma). 
        If None, limits are automatically inferred from the data.

    - `types_of_saving` : list of str or None
        List of output file formats for saving the figure (e.g., ["png", "pdf"]). 
        If empty, the figure is not saved.

    - `output_path` : str
        Path where the output figure will be saved. If "default" or empty, a default 
        subdirectory `figures` in the output folder from the configuration is used.

**PRE-REQUISITES**: This function requires the results of `data_calculate_eccentricity_parameters`, `plt_fit_ETP_signals_from_model` and `plt_fit_correlation_eccentricity_and_solution`.
