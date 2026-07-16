.. _plt_significance_functions:

*plt_significance_* detailed function guide
===========================================

.. image:: https://img.shields.io/badge/python-3.11%2B-blue
    :alt: Python Version

.. image:: https://img.shields.io/badge/license-GPLv3-blue
    :alt: License
    
.. function:: astrogeofit.plt_significance_test_results

.. code-block:: python

    plt_significance_test_results(
        configuration_file_path:str,
        folder_with_the_results:Optional[str],
        file_with_the_results:Optional[str], 
        section_title: list,
        types_of_saving:Optional[list[str]],
        output_path:Optional[str]
    )

Plot the results of a significance test comparing data fits to autoregressive (AR) noise fits.

This function generates a multi-panel figure showing histograms of $R^2_C$ values for data and AR model
fits across different numbers of knots (model complexity). For each subplot:
    - Histograms for both data and noise fits are shown.
    - The mean $R^2_C$ value for the data is indicated.
    - A Gaussian fit is plotted over the AR model histogram.
    - The False Alarm Probability (FAP) is computed and displayed.

Parameters:
    - `configuration_file_path` : str
        Path to the configuration file used for the AstroGeoFit analysis.

    - `folder_with_the_results` : str or None
        Path to the folder containing the result files. If None, the default folder 
        defined in the configuration file is used.

    - `file_with_the_results` : str or None
        Specific filename for the results. If None, a default filename is assumed.

    - `section_title`: str or None 
        Title or header to display above the plot. If None, it will set "X".

    - `types_of_saving` : list of str or None
        List of output file formats for saving the figure (e.g., ["png", "pdf"]). 
        If empty, the figure is not saved.

    - `output_path` : str
        Path where the output figure will be saved. If "default" or empty, a default 
        subdirectory `figures` in the output folder from the configuration is used.

Example:
        
    configuration_file_path = "path/to/the/configuration_file.yaml"
    folder_with_the_results = None
    file_with_the_results = None
    plt_significance_test_results(configuration_file_path, folder_with_the_results, file_with_the_results, "Significance Test Results")
    # Displays the figure with subplots for each number of knots and FAP annotations.

