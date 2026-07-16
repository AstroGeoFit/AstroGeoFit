.. _astrogeofit_class:

AstroGeoFit_tool class
======================

.. image:: https://img.shields.io/badge/python-3.11%2B-blue
    :alt: Python Version

.. image:: https://img.shields.io/badge/license-GPLv3-blue
    :alt: License

This class represents the main core of the tool, which will recieve all the necessary parameters
to obtain the results of the different sections (significance test, fittig, mcmc and weight computation).
This class contains the function `run()` which will execute the tool itself.

AstroGeoFit_tool_class
----------------------

.. class:: astrogeofit.AstroGeoFit_tool

.. code-block:: python

    AstroGeoFit_tool(
        config_file_path:str,
        exec_all:Optional[bool],
        fit: Optional[bool],
        significance: Optional[bool],
        mcmc: Optional[bool],
        weight: Optional[bool],
        path_fit_results:Optional[str],
        path_mcmc_results:Optional[str],,
        savefigs:bool,
        logs:bool,
        remote:bool,
    )

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

Parameters:
    - `config_file_path`: str  
        Path to the configuration file.

    - `exec_all`: bool or None 
        Whether to execute all available processes. Defaults to True.

    - `fit`: bool or None 
        Execute the fitting process. Defaults to False.

    - `significance`: bool or None 
        Execute the significance test. Defaults to False.

    - `mcmc`: bool or None
        Run the MCMC chain analysis. Defaults to False.

    - `weight`: bool or None
        Compute or apply weighting. Defaults to False.

    - `path_fit_results`: str or None 
        Path to output fit results. Only required if the execution calculates the mcmc or the weights but not the fitting part.

    - `path_mcmc_results`: str or None 
        Path to output MCMC results. Only required if the execution calculates the weights but not the mcmc part.

    - `savefigs`: bool or None 
        Whether to save figures. Defaults to False.

    - `logs`: bool or None 
        Whether to generate log files. Defaults to False.

    - `remote`: bool or None  
        Set to True if executing remotely. This makes Defaults to False.

Methods:
    `run()`: Executes the requested AstroGeoFit processes using the internal main function.

----

.. method:: astrogeofit.AstroGeoFit_tool.run

.. code-block:: python

    run()
    
Executes the requested AstroGeoFit processes using the internal main function.

This method acts as a wrapper around the `_astrogeofit_main_function`, passing along all
relevant configuration parameters required for different AstroGeoFit execution modes such as 
fitting, significance testing, MCMC simulation, and plotting.

The parameters used to execute the main astrogeofit function are the ones set when creating the class.
        
Returns:
    None

----