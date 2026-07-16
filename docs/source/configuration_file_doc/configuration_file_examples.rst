.. _config-example:

Example Configuration File
==========================
.. image:: https://img.shields.io/badge/python-3.11%2B-blue
    :alt: Python Version

.. image:: https://img.shields.io/badge/license-GPLv3-blue
    :alt: License

Table of Examples
-----------------

-  :ref:`fit_mcmc_age_ecc`
-  :ref:`fit_no_age_ecc`
-  :ref:`sig_test`
-  :ref:`mcmc_weight`
----

.. _fit_mcmc_age_ecc:

Execute Genetic Algorithm Fitting and MCMC With a Age-Depth Model and a Eccentricity Solution
---------------------------------------------------------------------------------------------

.. code:: yaml

    data_set:
        data_path:  path/dataset/ODP926.txt
        depth_column_name: depth                                
        proxy_column_name: proxy                                
        skiprows: 0                                             
        delimiter: " "
        header: True

    age_depth_model_data:
        data_path: path/model/U926_AgeModel_Wilkens_etal_2017.txt
        depth_column_name: Depth
        age_column_name: age
        skiprows: 0
        delimiter: ","
        header: True

    eccentricity_solution_data:
        data_path: path/solution/La2010d_ecc3L.dat
        age_column_name: 
        eccentricity_column_name: 
        skiprows: 0
        delimiter: \s+
        header: False
        start_time: 0e3
        final_time: 20e3

    data_model_parameters:
        sedimentation_rate_min: 0.4
        sedimentation_rate_max: 5

        #unit: arcsec/yr
        ### gi frequencies should be strictly in the order of g1, g2, g3, g4.
        frequency_values:
            use_precession: True
            use_eccentricity: True
            use_tilt: True

            p0_values: [50,55]
            gi_values: [[5.45, 5.75], [7.43, 7.48], [17.1, 17.4], [17.7, 18.0]]
            si_values: [[-19., -18.8], [-17.85, -17.7]] 

            g5_value: default
            s6_value: default

    genetic_algorithm_parameters:
        seed: 
        interpolator: linear
        metric_type: r2
        population_size: 300
        number_generations: 300
        number_processors_used: 10
        number_algorithm_solutions: 30
        list_number_genes: [2,  4,  8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96]

    mcmc_parameters:
        length_mcmc_chains: 5000
        discard: 2000
        thin: 50
        number_processors_used_mcmc: 15
        number_mcmc_solutions: 15
        list_of_genes_mcmc: [ 2, 12, 24, 36, 48, 60, 72, 84, 96]
        prior_distributions: 
            p0_distribution: gaussian
            gi_distribution: gaussian
            si_distribution: gaussian
        prior_frequencies:
            p0_prior: [50.6443465,0.15]
            gi_prior: [[5.579378, 0.055], [7.456665, 0.004], [17.366595, 0.030], [17.910194, 0.032]]
            si_prior: [[-18.845166, 0.047], [-17.758310, 0.023]]

    # ------------------- OUTPUT FILE TO GET THE RESULTS ------------------------ #
    output_folder:  path/to/the/output/folder
    output_file_name: name_of_the_result_files_ODP_926

----

.. _fit_no_age_ecc:

Execute Genetic Algorithm Fitting Without a Age-Depth Model and a Eccentricity Solution
---------------------------------------------------------------------------------------------

.. code:: yaml

    data_set:
        data_path:  path/dataset/ODP926.txt
        depth_column_name: depth                                
        proxy_column_name: proxy                                
        skiprows: 0                                             
        delimiter: " "
        header: True

    data_model_parameters:
        sedimentation_rate_min: 0.4
        sedimentation_rate_max: 5

        #unit: arcsec/yr
        ### gi frequencies should be strictly in the order of g1, g2, g3, g4.
        frequency_values:
            use_precession: True
            use_eccentricity: True
            use_tilt: True

            p0_values: [50,55]
            gi_values: [[5.45, 5.75], [7.43, 7.48], [17.1, 17.4], [17.7, 18.0]]
            si_values: [[-19., -18.8], [-17.85, -17.7]] 

            g5_value: default
            s6_value: default

    genetic_algorithm_parameters:
        seed: 
        interpolator: linear
        metric_type: r2
        population_size: 300
        number_generations: 300
        number_processors_used: 10
        number_algorithm_solutions: 30
        list_number_genes: [2,  4,  8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96]

    # ------------------- OUTPUT FILE TO GET THE RESULTS ------------------------ #
    output_folder:  path/to/the/output/folder
    output_file_name: name_of_the_result_files_ODP_926

----

.. _sig_test:

Execute Significance Test
-------------------------

.. code:: yaml

    data_set:
        data_path:  path/dataset/ODP926.txt
        depth_column_name: depth                                
        proxy_column_name: proxy                                
        skiprows: 0                                             
        delimiter: " "
        header: True

    data_model_parameters:
        sedimentation_rate_min: 0.4
        sedimentation_rate_max: 5

        #unit: arcsec/yr
        ### gi frequencies should be strictly in the order of g1, g2, g3, g4.
        frequency_values:
            use_precession: True
            use_eccentricity: True
            use_tilt: True

            p0_values: [50,55]
            gi_values: [[5.45, 5.75], [7.43, 7.48], [17.1, 17.4], [17.7, 18.0]]
            si_values: [[-19., -18.8], [-17.85, -17.7]] 

            g5_value: default
            s6_value: default

    significance_test_parameters:
        seed:
        interpolator: linear
        metric_type: loglike
        population_size: 300
        number_generations: 300
        number_processors_used: 10
        number_algorithm_executions: 50
        number_noise_executions: 50
        list_number_genes: [2, 12, 24, 36]

    # ------------------- OUTPUT FILE TO GET THE RESULTS ------------------------ #
    output_folder:  path/to/the/output/folder
    output_file_name: name_of_the_result_files_ODP_926

----

.. _mcmc_weight:

Example MCMC and Weights (Only using eccentricity and precession. Using uniform distributions for prior values )
----------------------------------------------------------------------------------------------------------------

.. code:: yaml

    data_set:
        data_path:  path/dataset/ODP926.txt
        depth_column_name: depth                                
        proxy_column_name: proxy                                
        skiprows: 0                                             
        delimiter: " "
        header: True

    data_model_parameters:
        sedimentation_rate_min: 0.4
        sedimentation_rate_max: 5

        #unit: arcsec/yr
        ### gi frequencies should be strictly in the order of g1, g2, g3, g4.
        frequency_values:
            use_precession: True
            use_eccentricity: True
            use_tilt: False

            p0_values: [50,55]
            gi_values: [[5.45, 5.75], [7.43, 7.48], [17.1, 17.4], [17.7, 18.0]]
            si_values: 

            g5_value: default
            s6_value: default

    mcmc_parameters:
        length_mcmc_chains: 5000
        discard: 2000
        thin: 50
        number_processors_used_mcmc: 15
        number_mcmc_solutions: 15
        list_of_genes_mcmc: [ 2, 12, 24, 36, 48, 60, 72, 84, 96]
        prior_distributions: 
            p0_distribution: uniform
            gi_distribution: uniform 
            si_distribution: uniform
        prior_frequencies:
            p0_prior: [50.6443465,0.15]
            gi_prior: [[5.579378, 0.055], [7.456665, 0.004], [17.366595, 0.030], [17.910194, 0.032]]
            si_prior: 

    weight_calcula_configuration:
        number_weight_evaluation_per_chain: 5
        number_processors_used_weights: 5
        stability_factor: 1e-8
        pareto_smoothing: True

    # ------------------- OUTPUT FILE TO GET THE RESULTS ------------------------ #
    output_folder:  path/to/the/output/folder
    output_file_name: name_of_the_result_files_ODP_926

----

.. _fit_fix_freqs:

Execute Genetic Algorithm Fitting with fixed frequencies (only using eccentricity and precession)
---------------------------------------------------------------------------------------------

.. code:: yaml

    data_set:
        data_path:  path/dataset/ODP926.txt
        depth_column_name: depth                                
        proxy_column_name: proxy                                
        skiprows: 0                                             
        delimiter: " "
        header: True

    data_model_parameters:
        sedimentation_rate_min: 0.4
        sedimentation_rate_max: 5

        #unit: arcsec/yr
        ### gi frequencies should be strictly in the order of g1, g2, g3, g4.
        frequency_values:
            use_precession: True
            use_eccentricity: True
            use_tilt: 

            p0_values: 55
            gi_values: [5.75, 7.48, 17.4, 18.0]
            si_values:  

            g5_value: default
            s6_value: default

    genetic_algorithm_parameters:
        seed: 
        interpolator: linear
        metric_type: r2
        population_size: 300
        number_generations: 300
        number_processors_used: 10
        number_algorithm_solutions: 30
        list_number_genes: [2,  4,  8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96]

    # ------------------- OUTPUT FILE TO GET THE RESULTS ------------------------ #
    output_folder:  path/to/the/output/folder
    output_file_name: name_of_the_result_files_ODP_926

----