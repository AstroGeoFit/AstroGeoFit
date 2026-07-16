.. _synthetic_example:

Basic Example with Synthetic Data
=================================

.. image:: https://img.shields.io/badge/python-3.11%2B-blue
    :alt: Python Version

.. image:: https://img.shields.io/badge/license-GPLv3-blue
    :alt: License

This basic example contains synthetic data created for the porpouse of testing the tool during its development

-----

Files in the *data* Folder
--------------------------

- **synthetic_data_ETP.txt**: Synthetic data.
- **age_model_ETP.txt**: Synthetic age-depth model that is adapted to the data.
- **La2010d_ecc3L.txt**: Astronomical solution data (eccentricity component).

-----

Execute the Example
-------------------

.. note::
    - The folder `examples` has to bee in the same directory as the `AstroGeoFit_tool.py` file. If not, this will not work.
    
To run the example, open a terminal in the `AstroGeoFit` folder. Once this is open execute the following

.. code-block:: bash

    python AstroGeoFit_tool.py --basic_example

This command will execute the significance test, the genetic algorithm fitting, the MCMC, and the MCMC weights calculation.

If you only want to run specific components (e.g., the genetic algorithm fitting and the MCMC), you can use:

.. code-block:: bash

    python AstroGeoFit_tool.py ----basic_example -fit -mcmc

The full list of execution options can be found in :ref:`tool_documentation`.

-----

Visualization of the Results
----------------------------

To visualize the results, you can use the two Jupyter notebooks provided. To obtain the results of the `basic_example` data, just modify the variable `configuration_file_path` in the following way:

.. code-block:: python

    configuration_file_path = "basic_example"

And then the visualization notebook will be ready to be executed.

For a detailed guide on using the notebooks and interpreting the figures, please refer to :ref:`visualization-results`.