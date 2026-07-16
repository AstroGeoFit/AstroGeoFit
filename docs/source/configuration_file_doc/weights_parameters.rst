.. _config-weigth-parameters:

Weight Parameters
-----------------

.. image:: https://img.shields.io/badge/python-3.11%2B-blue
    :alt: Python Version

.. image:: https://img.shields.io/badge/license-GPLv3-blue
    :alt: License
    
Defines the configuration for evaluating the model weights in the analysis.

.. code:: yaml

    weight_calcula_configuration:
        number_weight_evaluation_per_chain: 20
        number_processors_used_weights: 20
        stability_factor: 1e-8
        pareto_smoothing: True

Parameters
~~~~~~~~~~

- **number_weight_evaluation_per_chain**: `int`  
        Number of weight evaluations to perform per MCMC chain. A suggested value would be equal or larger than *20*.

- **number_processors_used_weights**: `int`  
        Number of processors used in parallel for weight computations. This value has to be lower or equal to the value set in **number_weight_evaluation_per_chain**.

- **stability_factor**: `float`  
        Small constant added for numerical stability during weight calculations. A suggested value would be *1e-8*.

- **pareto_smoothing**: `bool`  
        Whether to apply Pareto smoothing to the importance weights. Normally it is better to set this value as *True* as it is seen
        to give better results.
