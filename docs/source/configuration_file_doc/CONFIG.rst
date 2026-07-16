.. _configuration_guide:

Configuration File Guide
========================

.. image:: https://img.shields.io/badge/python-3.11%2B-blue
   :alt: Python Version

.. image:: https://img.shields.io/badge/license-GPLv3-blue
   :alt: License

The **AstroGeoFit** tool uses a ``.yml`` configuration file to define
the parameters required for execution. This section
provides an overview of each section and its expected values.

This file is included file, in the package `astrogeofit_latest.zip`. To download the package check the section :ref:`download_the_package`

.. note::

   - All the variables marked as *Optional* are optional variables. If no value is introduced they will be treated as **Null** values or **default** values depending of the variable.

   - Not all the sections of the configuration file are necessary. 

Configuration of the Datasets
-----------------------------

A guide to set the parameters to configurate the obtention of the three main datasets (Primarly Data, Age-Depth Model and Eccentricity Solution) 
is given in the following section :ref:`config-dataset-parameters`.

Data Model Parameters
---------------------

This section defines the sedimentation rate range and astronomical frequency parameters (precession, eccentricity, tilt) used in the data model. 
See :ref:`data_model_parameters` for details.

Significance Test 
-----------------

Defines the parameters required to run the Significance Test, including configuration for the genetic algorithm and evaluation settings.
See :ref:`significance-test-parameters`.

Genetic Algorithm
-----------------

This section defines the parameters that control the genetic algorithm used to optimize model fitting, including population size, number of generations, interpolation method, and fitness metric.
:ref:`config-genetic-algorithm-parameters`.

Bayesian MCMC
-------------

The MCMC analysis settings are defined in the :ref:`mcmc-parameters`.  
These parameters control the behavior of the Markov Chain Monte Carlo algorithm, including chain length, thinning, processor usage, and prior distributions for the orbital frequencies.

MCMC Weights
------------

This section defines the configuration used for computing model weights based on MCMC output. It includes the number of evaluations, parallelization settings, numerical stability factor, and whether Pareto smoothing is applied.

For full details, see: :ref:`config-weigth-parameters`.

Output Configuration
--------------------

For details on how to configure where and how results are saved, see the :ref:`config-output-parameters`.

Example Configuration File
--------------------------

Some examples of what the configuration file might look like depending on which is the goal of the analysis that wantst to be done.

The values shown on these examples are the values used to obtain the results for the *ODP_926*. See the examples in :ref:`config-example`.
