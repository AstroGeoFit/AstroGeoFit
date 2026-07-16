.. _mcmc-parameters:

MCMC Parameters
---------------

.. image:: https://img.shields.io/badge/python-3.11%2B-blue
    :alt: Python Version

.. image:: https://img.shields.io/badge/license-GPLv3-blue
    :alt: License

Defines settings for the Markov Chain Monte Carlo (MCMC) analysis used to evaluate uncertainty in frequency estimation.

.. code:: yaml

    mcmc_parameters:
        length_mcmc_chains: 10000
        discard: 5000
        thin: 100
        number_processors_used_mcmc: 8
        number_mcmc_solutions: 30
        list_of_genes_mcmc: [2, 12, 24, 36, 48, 64, 72, 82, 96]
        prior_distributions: 
            p0_distribution: gaussian
            gi_distribution: gaussian
            si_distribution: gaussian
        prior_frequencies:
            p0_prior: [50.6443465, 0.5]
            gi_prior: [[5.579378, 0.055], [7.456665, 0.004], [17.366595, 0.030], [17.910194, 0.032]]
            si_prior: [[-18.845166, 0.047], [-17.758310, 0.023]]

Parameters
~~~~~~~~~~

- **length_mcmc_chains**: `int`  
        Total number of MCMC iterations per chain.

- **discard**: `int`  
        Number of initial samples to discard (burn-in period). This helps to reduce the dependency of the chains on their starting points.

- **thin**: `int`  
        Interval for thinning the MCMC samples. Thinning is used to reduce autocorrelation in MCMC chains and obtain more independent samples that better represent the true posterior distribution. 

- **number_processors_used_mcmc**: `int`  
        Number of processors used for parallel computing.
        This value should be adapted to the hardware of the machine being used. It directly impacts the 
        performance and efficiency of the software.

- **number_mcmc_solutions**: `int`  
        Number of solutions selected from the genetic algorithm fitting solutions used to execute the MCMC. The GA solutions are used as starting points for the MCMC solutions.
        The tool takes the best X number of GA solutions.
        Must be equal or smaller than the value of `number_algorithm_solutions` (see :ref:`config-genetic-algorithm-parameters`).

- **list_of_genes_mcmc**: `list[int]`  
        Subset of the list of genes used in the GA section (see :ref:`config-genetic-algorithm-parameters`) used to calculate
        the MCMC solutions. It is important to choose the number of genes with the best results as the MCMC simulation is computationally expensive.

- **prior_distributions**: `dict`  
        Specifies the type of prior distribution used for each frequency group.  
        Options: `"gaussian"` or `"uniform"`.

  - **p0_distribution**: `str`  
        Defines the prior distribution for `p0`.

  - **gi_distribution**: `str`  
        Defines the prior distribution for `gi` (g1 to g4).

  - **si_distribution**: `str`  
        Defines the prior distribution for `si` (s3 and s5).

- **prior_frequencies**: `dict`  
        Mean and standard deviation (for Gaussian) or min/max (for uniform) for each frequency prior.

  - **p0_prior**: `list[float]`  
        Mean and std dev or min/max for `p0`.

  - **gi_prior**: `list[list[float]]`  
        List of mean and std (if it's gaussian) or min/max values (for uniform) for `gi` frequencies (g1 to g4).

  - **si_prior**: `list[list[float]]`  
        List of mean and std (if it's gaussian) or min/max values (for uniform) for `si` frequencies (s3 and s5).

.. note::
        - The value of **number_processors_used_mcmc** has to be equal or lower to the value of **number_mcmc_solutions**.
        - All the values of the **list_of_genes_mcmc** have to exist in the **list_number_genes** of the genetic algorithm section.
        - The values of **prior_frequencies** must be expressed in arcseconds per year (**arcsec/yr**). 