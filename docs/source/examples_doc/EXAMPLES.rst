.. _execution_examples:

Obtention of the Paper Results
==============================

.. image:: https://img.shields.io/badge/python-3.11%2B-blue
    :alt: Python Version

.. image:: https://img.shields.io/badge/license-GPLv3-blue
    :alt: License

To help users get started with *AstroGeoFit*, three examples are provided. These examples correspond to the data analyzed in the original publication associated with the tool (ODP_926, ODP_1260 and ODP_1262). 
They demonstrate how to configure, run *AstroGeoFit* and obtain the results shown in the paper.

Each example includes:
    - A folder named *data*, which contains:
        - Data set used.
        - Age Depth Model used.
        - Eccentricity Solution.
    - A folder named *configuration_file*, which contains:
        - Configuration File template.

These examples will be usefule to replicate the analyses described in the paper, making it easier to understand the tool's functionality and workflow.

Required Changes on the Configuration Files
-------------------------------------------

.. warning::
    - **NO MODIFICATIONS ARE REQUIRED TO APPLY TO THE NOTEBOOKS OR TO THE CONFIGURATION FILES**
    All the parameters are set on the configuration files in order to be executed directly. 

If the folder distribution or the directory of the files are changed, some variables must be modified accordingly. 
The variables that must be changed to grant a correct execution are:

- In the section *data_set*:
    - The variable **data_path**

- In the section *age_depth_model_data*:
    - The variable **data_path**

- In the section *eccentricity_solution_data*:
    - The variable **data_path**

- The variable **output_folder**.

- The variable **output_file_name**.

Explanations of the Examples
----------------------------

A detailed explanation of each example and how to run them can be found in the following pages:

- :ref:`synthetic_example`

- :ref:`odp_926_example`

- :ref:`odp_1260_example`

- :ref:`odp_1262_example`