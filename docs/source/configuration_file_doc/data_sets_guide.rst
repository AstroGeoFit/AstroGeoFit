.. _config-dataset-parameters:

Configuration for the Datasets Used
===================================

.. image:: https://img.shields.io/badge/python-3.11%2B-blue
    :alt: Python Version

.. image:: https://img.shields.io/badge/license-GPLv3-blue
    :alt: License

Main Dataset 
------------

This section defines the primary dataset used for the analysis.

.. note::

   Always use **absolute paths** to avoid file access errors.

.. code:: yaml

    data_set:
        data_path: 
        depth_column_name: 
        proxy_column_name:
        skiprows: 
        delimiter:
        header: 

Parameters
~~~~~~~~~~

- **data_path**: `str`  
        Absolute path to the dataset file. Accepted file extensions are: `.txt`, `.csv`, `.xlsx`, `.tab`, `.dat`.

- **depth_column_name**: `str`, `int`, or `None`  
        Name or index of the column that represents depth. If `None`, the first column will be used.

- **proxy_column_name**: `str`, `int`, or `None`  
        Name or index of the column representing proxy data. If `None`, the second column will be used.

- **skiprows**: `int` or `None`  
        Number of rows to skip at the beginning of the file, excluding the header.

- **delimiter**: `str` or `None`  
        Character used to separate values in the file. If `None`, the following defaults apply based on file extension:  
        `.txt` → `\\s+`, `.csv` → `,`, `.dat` → `" "`.

- **header**: `bool`  
        Whether the dataset file includes a header (`True` or `False`).

.. note::

   - Depth values must be in **meters**.  
   - Proxy values must be in **cm/kyr**.

----

.. _config-age-model-parameters:

Age-Depth Model
---------------

This section defines the file that maps depth to geological age. This file is optional. If omitted, some plots may not be available or appear differently.

.. note::

   Always use **absolute paths** to avoid file access errors.

.. code:: yaml

    age_depth_model_data:
        data_path:
        depth_column_name:
        age_column_name:
        skiprows:
        delimiter:
        header:

Parameters
~~~~~~~~~~

- **data_path**: `str`  
        Absolute path to the age-depth model file. Accepted file extensions are the same as the main dataset.

- **depth_column_name**: `str`, `int`, or `None`  
        Name or index of the column that represents depth.

- **age_column_name**: `str`, `int`, or `None`  
        Name or index of the column representing geological age.

- **skiprows**: `int` or `None`  
        Number of rows to skip at the beginning of the file, excluding the header.

- **delimiter**: `str` or `None`  
        Character used to separate values in the file.

- **header**: `bool`  
        Whether the file includes a header (`True` or `False`).

.. note::

   - Age values must be in **millions of years (Myr)**.

----

.. _config-eccentricity-parameters:

Eccentricity Solution Dataset
-----------------------------

This section defines the astronomical solution dataset. This file is optional. If not provided, plots involving eccentricity correlation will be disabled or raise errors in notebooks.

.. note::

   Always use **absolute paths** to avoid file access errors.

.. code:: yaml

    eccentricity_solution_data:
        data_path:
        age_column_name:
        eccentricity_column_name:
        skiprows:
        delimiter:
        header:
        start_time:
        final_time:

Parameters
~~~~~~~~~~

- **data_path**: `str`  
        Absolute path to the eccentricity solution file.

- **age_column_name**: `str`, `int`, or `None`  
        Name or index of the column that represents age.

- **eccentricity_column_name**: `str`, `int`, or `None`  
        Name or index of the column that represents eccentricity.

- **skiprows**: `int` or `None`  
        Number of rows to skip at the beginning of the file, excluding the header.

- **delimiter**: `str` or `None`  
        Character used to separate values in the file.

- **header**: `bool`  
        Whether the file includes a header (`True` or `False`).

- **start_time**: `float`  
        Starting age (in **thousands of years**) for correlation.

- **final_time**: `float`  
        Ending age (in **thousands of years**) for correlation.
