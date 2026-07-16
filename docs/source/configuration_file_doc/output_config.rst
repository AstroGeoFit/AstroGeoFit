.. _config-output-parameters:

Output Parameters
-----------------

.. image:: https://img.shields.io/badge/python-3.11%2B-blue
    :alt: Python Version

.. image:: https://img.shields.io/badge/license-GPLv3-blue
    :alt: License
    
Defines the folder and file name for saving results.

.. warning::
    - Use absolute paths to avoid errors when saving output files.

.. code:: yaml

    output_folder: path/to/results/testing_results_ODP926
    output_file_name: testing_results_ODP926

Parameters
~~~~~~~~~~

- **output_folder**: `str`  
    Absolute path to the directory where results will be stored.

- **output_file_name**: `str`  
    Name of the file used to save the output results.
