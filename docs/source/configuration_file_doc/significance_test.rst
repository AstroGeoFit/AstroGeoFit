.. _significance-test-parameters:

Significance Test Parameters
----------------------------

.. image:: https://img.shields.io/badge/python-3.11%2B-blue
    :alt: Python Version

.. image:: https://img.shields.io/badge/license-GPLv3-blue
    :alt: License

Defines the parameters used to execute the Significance Test of the fitting of the Genetic Algorithm. These variables must be set only when the test is to be performed.

.. code:: yaml

    significance_test_parameters:
        seed: 
        interpolator: default
        metric_type: loglike
        population_size: 150
        number_generations: 150
        number_processors_used: 50
        number_algorithm_executions: 1000
        number_noise_executions: 1000
        list_number_genes: [2, 12, 24, 36]

Parameters
~~~~~~~~~~

- **seed**: `int`, `list[int]`, or `None`  
    Random seed(s) for reproducibility. Set to `None` for random behavior.

- **interpolator**: `str` or `"default"`  
    Interpolation method to use. Options include `CubicSpline`, `Pchip`, `BSpline`, etc.  
    If `"default"` or not specified, uses **"CubicSpline"**.

- **metric_type**: `str` or `"default"`  
    Fitness metric type. Accepted values are `r2` or `loglike`.  
    If not specified, defaults to **"loglike"**.

- **population_size**: `int`  
    The number of individuals in each generation, known as the *population size*. A typical value of 100 often yields good-quality solutions.
    For improved accuracy, consider increasing the population size. Larger populations require more computational resources, so only scale up if your hardware capabilities can support it.

- **number_generations**: `int`  
    Total number of generations to run during optimization.
    It’s recommended to start with around 100 generations and adjust based on how the fitness metric evolves over time. 
    The number of generations should be sufficient for the fitness metric to reach a high plateau. Once this plateau is 
    maintained for several generations, increasing the number further is unlikely to improve the fit significantly.

- **number_processors_used**: `int`  
    Number of processors used for parallel computing.
    This value should be adapted to the hardware of the machine being used. It directly impacts the 
    performance and efficiency of the software. 

- **number_algorithm_executions**: `int`  
    Number of times that the genetic algorithm is executed for the significance test.
    This parameter significantly impacts the total run time of the software.
    Typically, around 30 executions are sufficient to ensure good solution diversity.3

- **number_noise_executions**: `int`  
    Number of noise-based solutions to generate for significance evaluation.
    This parameter significantly impacts the total run time of the software.

- **list_number_genes**: `list[int]`  
    Number of genes refers to the number of depth points at which the inverse sedimentation rate (SR⁻¹) is observed or calculated.
    We recommend starting with 2 genes, then increasing to 4, and continuing in steps of 4 (e.g., 8, 12, 16…) up to a suitable value such as 36 or 60.
    While it is possible to perform the fitting using a single number of genes, using a list of values for sequential fitting generally yields better performance.

.. note::
    - Note: population_size and number_generations should be set according to the genetic algorithm (GA) configuration parameters.
    - Note: It is common practice to set number_algorithm_executions equal to number_noise_executions.
    Suggested value: ≥ 1000 for both.