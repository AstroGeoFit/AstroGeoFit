.. _execution_page:

Execution of the Package
========================

.. toctree::
    :maxdepth: 2
    :hidden:
    
    TOOL_MODE
    python_library/index

.. image:: https://img.shields.io/badge/python-3.11%2B-blue
    :alt: Python Version

.. image:: https://img.shields.io/badge/license-GPLv3-blue
    :alt: License


The `AstroGeoFit` tool can be executed in two different ways: through the `AstroGeoFit_tool.py` script or by using the *Python functions* of the library.

.. note::
    - The execution of the *AstroGeoFit tool* bases its execution in all the values found in the *configuration_file*. Please check the :ref:`configuration_guide` before executing the software.

Executing the `AstroGeoFit_tool` Script
---------------------------------------

Using the Python script provided in the `astrogeofit` package is the easiest option, as it does not require any prior Python knowledge.
This script gives the user the possibility to use the software via terminal commands. To see how to use it, please refer to: :ref:`tool_documentation`.

Using the `AstroGeoFit` Library
-------------------------------

The functions of `AstroGeoFit` can be executed separately in a Python script, similar to any Python library. The description of the
available functions and the explanation of their parameters can be found in :ref:`library_documentation`.

How to Execute the Different Parts of the `AstroGeoFit Tool`
------------------------------------------------------------

`AstroGeoFit` can be divided in four separate parts: **Fitting of the GA**, **MCMC**, **MCMC Weights Calcula** and **Significance Test**. These parts are not independent of each other and there is 
a "standard" way to execute the full tool. 

- **Fitting of the GA**: This is the base step in order to obtain the results of the tool. The *fitting results* will be used to compute the *MCMC* and compared with the results of the *Significance Test*.
In resume, this section has to be always executed.

- **Significance Test**: We fit a noise model as the null hypothesis and generate synthetic datasets from it.
Each is fitted with the full model to compute weighted :math:`R^2`, forming a reference distribution.
Comparing this to :math:`R^2` values from real data allows us to test whether the null hypothesis holds. This method requires the **Fitting Results** to be able to compare with.

- **MCMC**: MCMC is used to sample from the posterior distribution of model parameters, capturing the uncertainty in their estimates. To explore multiple modes in the parameter space, we run MCMC chains from diverse high-ranking starting points identified by a genetic algorithm. 
The resulting samples are then combined using a weighted stacking approach to better reflect the relative importance of each mode. The *MCMC* process requires the results from the **Fitting of the GA**.


- **MCMC Weights Calcula**: As the name of the section suggests, the calculation of the weights is based on the **MCMC Results**. This process estimates the full posterior distribution by stacking multiple MCMC chains, each exploring different local modes, 
and weighting them using Bayesian evidence estimated via importance sampling. Therefore, we first need to execute the *MCMC* in order to obtain the *weight calculation*.

.. footer:: Requests/Inquiries: astrogeofit@astrogeo.eu