.. _plt_mcmc_functions:

*plt_mcmc_* detailed function guide
===================================

.. image:: https://img.shields.io/badge/python-3.11%2B-blue
    :alt: Python Version

.. image:: https://img.shields.io/badge/license-GPLv3-blue
    :alt: License
    
.. function:: astrogeofit.plt_mcmc_SR_per_num_knots

.. code-block:: python

    plt_mcmc_SR_per_num_knots(
        configuration_file_path:str,
        folder_with_the_results:Optional[str],
        file_with_the_results:Optional[str],
        x_axis_plot_limits: Optional[list],
        y_axis_plot_limits: Optional[list],
        ignore_weights: bool,
        types_of_saving: list,
        output_path: str,
    )

Plot the inverse sedimentation rate (SR⁻¹) derived from MCMC simulations 
for different numbers of knots in the age-depth model.

This function loads the results of an AstroGeoFit analysis, specifically the 
MCMC sampling results for different age-depth model complexities (defined by 
number of knots). It visualizes the median inverse sedimentation rate profiles, 
confidence intervals, MAP estimates, and optionally the nominal model.

Parameters:
    - `configuration_file_path` : str
        Path to the configuration file used for the AstroGeoFit analysis.

    - `folder_with_the_results` : str or None
        Path to the folder containing the result files. If None, the default folder 
        defined in the configuration file is used.

    - `file_with_the_results` : str or None
        Specific filename for the results. If None, a default filename is assumed.

    - `x_axis_plot_limits` : list or None
        Optional [xmin, xmax] for the depth axis. If None, limits are automatically inferred.

    - `y_axis_plot_limits` : list or None
        Optional [ymin, ymax] for the inverse SR axis. If None, limits are automatically inferred.

    - `ignore_weights` : bool
        If True, ignore the sampling weights when drawing MCMC samples.

    - `types_of_saving` : list of str
        List of file formats to save the resulting plot (e.g., ["png", "pdf"]). 
        If empty, the figure will not be saved.

    - `output_path` : str
        Output path where the figure will be saved. If "default" or empty, the 
        `figures` subfolder inside the output folder defined in the configuration 
        file is used.

----

.. function:: astrogeofit.plt_mcmc_aic_logprob_and_loglike_per_number_knots

.. code-block:: python

    plt_mcmc_aic_logprob_and_loglike_per_number_knots(
        configuration_file_path:str,
        folder_with_the_results:Optional[str],
        file_with_the_results:Optional[str],
        types_of_saving: list,
        output_path: str,
    )

Plot AIC, Logprob and Loglike metrics for MCMC-inferred models with varying numbers of knots.

This function evaluates the model fit for each number of knots used in an 
age-depth model by calculating the Akaike Information Criterion (AIC) and the 
coefficient of determination (R²) using MCMC posterior samples. The metrics 
are plotted together using a dual y-axis plot.

Parameters:
    - `configuration_file_path` : str
        Path to the configuration file used for the AstroGeoFit analysis.

    - `folder_with_the_results` : str or None
        Path to the folder containing the result files. If None, the default 
        folder specified in the configuration is used.

    - `file_with_the_results` : str or None
        Specific filename for the results. If None, the default name is assumed.

    - `types_of_saving` : list of str
        List of file formats in which to save the resulting plot (e.g., ["png", "pdf"]).

    - `output_path` : str
        Directory where the figure should be saved. If "default" or empty, it defaults 
        to a `figures` folder inside the output folder defined in the configuration.

----

.. function:: astrogeofit.plt_mcmc_aic_and_r2_per_number_knots

.. code-block:: python

    plt_mcmc_aic_and_r2_per_number_knots(
        configuration_file_path:str,
        folder_with_the_results:Optional[str],
        file_with_the_results:Optional[str],
        types_of_saving: list,
        output_path: str,
    )

Plot AIC and R² metrics for MCMC-inferred models with varying numbers of knots.

This function evaluates the model fit for each number of knots used in an 
age-depth model by calculating the Akaike Information Criterion (AIC) and the 
coefficient of determination (R²) using MCMC posterior samples. The metrics 
are plotted together using a dual y-axis plot.

Parameters:
    - `configuration_file_path` : str
        Path to the configuration file used for the AstroGeoFit analysis.

    - `folder_with_the_results` : str or None
        Path to the folder containing the result files. If None, the default 
        folder specified in the configuration is used.

    - `file_with_the_results` : str or None
        Specific filename for the results. If None, the default name is assumed.

    - `types_of_saving` : list of str
        List of file formats in which to save the resulting plot (e.g., ["png", "pdf"]).

    - `output_path` : str
        Directory where the figure should be saved. If "default" or empty, it defaults 
        to a `figures` folder inside the output folder defined in the configuration.


----

.. function:: astrogeofit.plt_mcmc_eccentricity_correlation_with_solution

.. code-block:: python

    plt_mcmc_eccentricity_correlation_with_solution(
        configuration_file_path:str,
        folder_with_the_results:Optional[str],
        file_with_the_results:Optional[str],
        t0_offset_range_plot: Optional[list],
        ecc_time_range: Optional[list],
        name_eccentricity_solution: Optional[str],
        types_of_saving: list,
        output_path: str,
    )

Plot the results of eccentricity from MCMC runs and compare them with an astronomical solution.

This function generates a multi-panel figure showing:
- The mean square weighted deviation (MSWD) across different t0 offsets.
- The aligned time offset between relative and absolute time.
- The fitted eccentricity envelope overlaid with a reference astronomical solution.

It also logs relevant metrics and saves the resulting figure in the specified formats.

Parameters:
    - `configuration_file_path` : str
        Path to the configuration file used in the analysis.

    - `folder_with_the_results` : str or None
        Path to the folder containing the MCMC result files. Can be `None` if `file_with_the_results` is provided.

    - `file_with_the_results` : str or None
        Specific path to a result file. Can be `None` if `folder_with_the_results` is provided.

    - `t0_offset_range_plot` : list or None
        Range of t0 offsets to evaluate MSWD. Should be a list of two floats [min_offset, max_offset].
        If `None`, it is inferred from the nominal age model if available.

    - `ecc_time_range` : list or None
        Time range for plotting the eccentricity envelope. If `None`, it's automatically inferred from the dataset.

    - `save_eccentricity_solutions`: bool
        Flag indicating whether to save the reconstructed eccentricity solutions and their parameters to output files.
        If set to True, the function will generate:

            A .txt file containing the depth, age, proxy values, and statistical summaries (MAP, median, 5th, and 95th percentiles) of the eccentricity.

            A .sol file containing the frequency components, their amplitudes (beta coefficients), and phase information derived from the MAP (maximum a posteriori) solution.

        These files allow post-analysis use of the eccentricity time series, e.g., for replication, comparison, or modeling.

        If False, these outputs are skipped.
            
    - `save_eccentricity_solutions_comment`: str or None
        A user-defined comment or note to be included in the header of the saved eccentricity solution output files (e.g., .txt and .sol). 
        This allows users to annotate saved results with context, version information, data origin, or other metadata

    - `name_eccentricity_solution` : str or None
        Label for the astronomical solution being compared (e.g., "La2004", "Astronomical Solution"). Default is "Astronomical Solution".

    - `types_of_saving` : list
        List of formats in which to save the output figure (e.g., ['pdf', 'png']).

    - `output_path` : str
        Directory where the output figures and logs should be saved. If 'default', uses the default output folder
        specified in the analysis configuration.

**PRE-REQUISITES**: This function needs the results of `data_obtain_mcmc_results_per_a_number_of_knots`.

----

.. function:: astrogeofit.plt_mcmc_summary_of_results

.. code-block:: python

    plt_mcmc_summary_of_results(
        configuration_file_path:str,
        folder_with_the_results:Optional[str],
        file_with_the_results:Optional[str],
        x_axis_limits_peridogram: Optional[list],
        y_axis_limits_peridogram: Optional[list],
        peridogram_scale_x_axis: Optional[str],
        peridogram_scale_y_axis: Optional[str],
        x_axis_limits_depth: Optional[list],
        types_of_saving: list,
        output_path: str,
    )

Generate a comprehensive summary plot from MCMC results in astrochronological analysis.

This function plots five panels:
    1. Spectral power analysis of log-likelihood for each number of knots.
    2. Median time model with 95% confidence interval.
    3. Inverse sedimentation rate over time with credible intervals.
    4. Fit of data in the time domain with posterior predictive samples.
    5. Decomposition of the ETP (Eccentricity, Tilt, Precession) signal using fitted knots.

The figure is saved in the specified formats and the output path is determined by the provided configuration.

Parameters:
    - `configuration_file_path` : str  
        Path to the configuration file used in the analysis.

    - `folder_with_the_results` : str or None  
        Path to the folder containing the MCMC result files. Can be `None` if `file_with_the_results` is provided.

    - `file_with_the_results` : str or None  
        Specific path to a result file. Can be `None` if `folder_with_the_results` is provided.

    - **x_axis_limits_peridogram** : list or None  
        Limits for the x-axis in the peridogram plot. Should be [min, max] in frequency units. If `None`, defaults are used.

    - **y_axis_limits_peridogram** : list or None  
        Limits for the y-axis in the peridogram plot. Should be [min, max]. If `None`, auto-scaled from data.

    - `peridogram_scale_x_axis` : str or None  
        Scale for the x-axis of the peridogram plot. Can be `"linear"` or `"log"`. Defaults to `"linear"`.

    - `peridogram_scale_y_axis` : str or None  
        Scale for the y-axis of the peridogram plot. Can be `"linear"` or `"log"`. Defaults to `"linear"`.

    - **x_axis_limits_depth** : list or None  
        Limits for the x-axis (depth) used in multiple plots. Should be [min_depth, max_depth]. If `None`, it spans the dataset range.

    - `file_with_the_results` : list of str  
        List of file formats to save the final plot. Supported formats include `"png"`, `"pdf"`, and `"svg"`. `"tab"` is ignored.

    - `output_path` : str  
        Path where the output plot should be saved. If `"default"`, it's inferred from the configuration file.

**PRE-REQUISITES**: This function needs the results of `data_obtain_mcmc_results_per_a_number_of_knots`.

----

.. function:: astrogeofit.plt_mcmc_prior_frequencies_distributions

.. code-block:: python

    plt_mcmc_prior_frequencies_distributions(
        configuration_file_path:str,
        folder_with_the_results:Optional[str],
        file_with_the_results:Optional[str],
        types_of_saving: list,
        output_path: str,
    )

Plot the prior and posterior distributions of the frequencies explored by MCMC.

This function generates subplots of the marginal prior and posterior distributions 
for the different orbital frequencies (precession, obliquity, etc.) used in the 
astrochronological model. It overlays the MAP value, the current astronomical values, 
and posterior statistics such as the weighted median and standard deviation.
    
It also logs summary statistics (MAP, median, std) for each frequency and saves the plot 
in one or more formats.

Parameters:
    - `configuration_file_path` : str
        Path to the YAML configuration file used to run the analysis.
        
    - `folder_with_the_results` : str or None
        Path to the folder containing the result files. If `file_with_the_results` is given,
        this can be `None`.

    - `file_with_the_results` : str or None
        Path to a specific result file. If provided, `folder_with_the_results` can be `None`.

    - `file_with_the_results` : list
        List of file formats to save the plot in. Supported formats are image types like 
        ['png', 'pdf', 'svg'], etc. 'tab' is ignored with a warning.

    - `output_path` : str
        Directory where the figure(s) and logs should be saved. If set to "default", uses the 
        default output folder specified in the configuration file.

**PRE-REQUISITES**: This function needs the results of `data_obtain_mcmc_results_per_a_number_of_knots`.

----

.. function:: astrogeofit.plt_mcmc_phase_of_frequencies

.. code-block:: python

    plt_mcmc_phase_of_frequencies(
        configuration_file_path:str,
        folder_with_the_results:Optional[str],
        file_with_the_results:Optional[str],
        types_of_saving: list, 
        output_path: str
    )

Generate and save a multi-panel plot showing the phases of orbital frequencies 
derived from MCMC solutions applied to astronomical time scale inference.

This function loads the results of an MCMC analysis (using AstroGeoFit), computes
the phase of orbital frequencies (such as eccentricity, obliquity, and precession),
and plots histograms of the sampled phase distributions along with the MAP (maximum
a posteriori) and median values. It also logs summary statistics to a file.

Parameters:
    - `configuration_file_path` : str
        Path to the configuration file used for the AstroGeoFit model run.
        
    - `folder_with_the_results` : str or None
        Path to the folder containing the output results of the MCMC run.
        Can be None if `file_with_the_results` is specified instead.
        
    - `file_with_the_results` : str or None
        Path to a specific result file to load instead of an entire folder.
        Can be None if `folder_with_the_results` is provided.

    - `file_with_the_results` : list of str
        A list of file formats to save the final phase plot (e.g., ['png', 'pdf']).
        If 'tab' is included, it is ignored for this plot type.

    - `output_path` : str
        Directory path to store the generated plot and the log file. If empty or set 
        to "default", the output directory defined in the configuration will be used.

**PRE-REQUISITES**: This function needs the results of `data_obtain_mcmc_results_per_a_number_of_knots`.

