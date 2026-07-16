.. _tutorial:

Installation and Execution Guide
================================

This guide provides detailed, step-by-step instructions to install and run the ``AstroGeoFit`` tool.

Pre-Requisites
--------------

Before proceeding, ensure all required dependencies are installed.  
The complete list of pre-requisites can be found in the section :ref:`prerequisites`.

----

Download the Package
--------------------

These instructions assume you have already downloaded the package from the official *AstroGeo* webpage.  
If you have not done so, please refer to the guide in the section :ref:`download_the_package`.

----

Installation of the Tool
------------------------

Once pre-requisites are met and the package is downloaded, proceed with installing ``AstroGeoFit``.  
Detailed installation steps are provided in the section :ref:`installation_updates`.

----

Execution of  ``AstroGeoFit``
----------------------------------

The simplest way to run ``AstroGeoFit``—without any Python programming knowledge—is through the ``AstroGeoFit_tool.py`` script.  
This script allows interaction with the software using shell (terminal) commands only.

The script is included within the downloaded package.

The Configuration File
^^^^^^^^^^^^^^^^^^^^^^

Before starting an analysis, a **configuration file** is required.  
This file specifies all parameters necessary to run the analysis and can be downloaded from the official source.

Instructions on how to prepare the configuration file are detailed in the section :ref:`configuration_guide`.

The configuration file is also included in the downloaded package.

Executing the Script
^^^^^^^^^^^^^^^^^^^^

Once the configuration file is ready, run the script following the instructions in :ref:`tool_documentation`.

Using the `AstroGeoFit` Python Library
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The core functionality of `AstroGeoFit` is also accessible as a **Python library**.  
Documentation describing all available functions and their parameters is provided in :ref:`library_documentation`.

----

Obtaining Results and Visualization
-----------------------------------

After completing the analysis, results can be visualized using the provided Jupyter notebooks:

- ``AstroGeoFit_Analysis_Model_Fitting.ipynb``  
- ``AstroGeoFit_Analysis_MCMC.ipynb``

These notebooks are located inside the `astrogeofit` package.  

Alternatively, plots can be generated directly through the Python library’s functions.  
Comprehensive information on visualizing results is available in :ref:`visualization-results`.

----

Execution of the Datasets Found on the Paper
--------------------------------------------

In the same `AstroGeoFit` package there is a folder called `examples` which contain the required files to obtain the results of the three datasets
that appear in the paper **AstroGeoFit: A Genetic Algorithm and Bayesian approach for the Astronomical Calibration of the Geological Timescale**,
*Paleoceanography and Paleoclimatology*, 2025.

The guide to execute these examples can be found in :ref:`execution_examples`.


