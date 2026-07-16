.. _config-genetic-algorithm-parameters:

Genetic Algorithm Parameters
----------------------------

.. image:: https://img.shields.io/badge/python-3.11%2B-blue
    :alt: Python Version

.. image:: https://img.shields.io/badge/license-GPLv3-blue
    :alt: License

Controls the optimization algorithm used for model fitting.

.. code:: yaml

    genetic_algorithm_parameters:
        seed:
        interpolator: default
        metric_type: loglike
        population_size: 100
        number_generations: 100
        number_processors_used: 50
        number_algorithm_solutions: 50
        list_number_genes: [2, 4, 8, 12, 16, 20, 24, 30, 36, 42, 48, 56, 64, 72, 82, 96]

Parameters
~~~~~~~~~~

- **seed**: `int`, `list[int]`, or `None`  
        Random seed(s) for reproducibility. Use `None` or leave empty for random behavior.

- **interpolator**: `str` or `"default"`  
        Interpolation method used to get a continuous curve of the inverse sedimentation rate (SR^-1) over depth given some values of SR^-1 by the optimized knods.
        The ptions include `CubicSpline`, `Pchip`, `BSpline`, etc.  
        If `"default"` or not specified, defaults to **"CubicSpline"**.

- **metric_type**: `str` or `"default"`  
        Metric to evaluate the fitness of solutions. Accepted values: `loglike`, `r2`.  
        Defaults to **"loglike"** if not specified.

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

- **number_algorithm_solutions**: `int`
        Number of times the genetic algorithm is executed to obtain diverse solutions.
        This parameter significantly impacts the total run time of the software.
        Typically, around 30 executions are sufficient to ensure good solution diversity.

- **list_number_genes**: `list`
        List of number of genes that the algorithm will use to optimize the inverse sedimentation rate function.
        A *number of genes* refers to the number of depth points at which the inverse sedimentation rate (SR⁻¹) is observed or calculated.
        We recommend starting with 2 genes, then increasing to 4, and continuing in steps of 4 (e.g., 8, 12, 16…) up to a suitable value such as 36 or 60.
        While it is possible to perform the fitting using a single number of genes, using a list of values for sequential fitting generally yields better performance.

.. note::
        - The value of **number_processors_used** has to be equal or lower to the value of **number_algorithm_solutions**.