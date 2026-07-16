.. _data_functions:

*data_* detailed function guide
===============================

.. image:: https://img.shields.io/badge/python-3.11%2B-blue
    :alt: Python Version

.. image:: https://img.shields.io/badge/license-GPLv3-blue
    :alt: License
    
.. function:: astrogeofit.data_calculate_eccentricity_parameters

.. code-block:: python

    data_calculate_eccentricity_parameters(
        configuration_file_path:str,
        folder_with_the_results:Optional[str],
        file_with_the_results:Optional[str],
        positive_feedback: bool,
        number_of_knots_to_explore: Optional[int],
        name_eccentricity_solution: Optional[str],
    )

Calculates and saves the eccentricity-related parameters from a time-depth model fit using 
astronomical components (eccentricity, precession, tilt).

This function loads the results of a genetic algorithm fit, reconstructs the inferred time scale, 
and computes the contribution of each astronomical component using linear regression. 
The resulting parameters and derived time series are stored in a temporary file.

Parameters:
    - `configuration_file_path` : str
        Path to the configuration file used to perform the time-depth model fitting.

    - `folder_with_the_results` : str or None
        Path to the folder where the results of the fitting procedure are stored. If None, 
        the default from the configuration file is used.

    - `file_with_the_results` : str or None
        Specific filename of the results to load, if available. If None, uses the default.

    - `positive_feedback` : bool
        Indicates whether the eccentricity signal should be positively or negatively correlated 
        with the data (used to adjust the sign of the resulting prediction).

    - `number_of_knots_to_explore` : int or None
        Number of knots (genes) used in the fitting process to be explored. If None, the last 
        set of results is used by default.

    - `name_eccentricity_solution` : str or None
        Custom name assigned to the eccentricity solution (e.g., "La2004", "Berger2020"). 
        Used for labeling and tracking. If None or empty, defaults to "Astronomical Solution".

----

.. function:: astrogeofit.data_obtain_mcmc_results_per_a_number_of_knots

.. code-block:: python

    data_obtain_mcmc_results_per_a_number_of_knots(
        configuration_file_path:str,
        folder_with_the_results: Optional[str],
        file_with_the_results: Optional[str],
        selected_num_genes: Optional[int],
        positive_feedback: bool,
        use_prec_env: bool,
        ignore_weights: bool,
    )

Extract and process MCMC results for a selected number of knots in an astrochronological model.

This function loads analysis and MCMC results from disk, extracts the parameters corresponding to a 
specified number of knots (if provided), and reconstructs the astrochronological time series using 
the MAP estimate and bootstrapped parameter samples. It computes time series for eccentricity (either 
based on precession envelope or gi-gj method), integrates the inverse sedimentation rate to obtain an 
inferred timescale, and saves all relevant results to a temporary pickle file.

Parameters:
    - `configuration_file_path` : str  
        Path to the configuration YAML file used for the analysis.

    - `folder_with_the_results` : str or None  
        Folder containing the result files. If `None`, the default path is used.

    - `file_with_the_results` : str or None  
        Base name of the results file. If `None`, the default name is used.

    - `selected_num_genes` : int or None  
        Number of knots (genes) to extract results for. If `None`, uses the number of depth genes from the selected solution.

    - `positive_feedback` : bool  
        Whether to multiply the eccentricity feedback signal by `+1` (`True`) or `-1` (`False`).

    - `use_prec_env` : bool  
        If `True`, uses the precession envelope for eccentricity. If `False`, uses the `gi - gj` formulation.

    - `ignore_weights` : bool  
        If `True`, ignores MCMC weights when resampling bootstrap parameters.


----