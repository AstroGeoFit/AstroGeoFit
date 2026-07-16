.. _data_model_parameters:

Data Model Parameters
=====================

.. image:: https://img.shields.io/badge/python-3.11%2B-blue
    :alt: Python Version

.. image:: https://img.shields.io/badge/license-GPLv3-blue
    :alt: License

Defines the parameters for sedimentation rates and fundamental frequency
ranges or values.

.. code:: yaml

      data_model_parameters:
        sedimentation_rate_min: 0.4
        sedimentation_rate_max: 5

        frequency_values:
                use_precession: True
                use_eccentricity: True
                use_tilt: True

                p0_values: [50, 55]
                gi_values: [[5.45, 5.75], [7.43, 7.48], [17.1, 17.4], [17.7, 18.0]]
                si_values: [[-19.0, -18.8], [-17.85, -17.7]]

                g5_value: default
                s6_value: default

Parameters
----------

Sedimentation Rate Range
~~~~~~~~~~~~~~~~~~~~~~~~

.. note::
    - To set a range of sedimentation rate to be searched in the AstroGeoFit, please edit sedimentation_rate_min, sedimentation_rate_max.

- **sedimentation_rate_min**: `float`  
      Minimum sedimentation rate, in **cm/kyr**.

- **sedimentation_rate_max**: `float`  
      Maximum sedimentation rate, in **cm/kyr**.

Frequency Model 
~~~~~~~~~~~~~~~~~~~~~~~~

- **frequency_values**: `dict`  
      Dictionary of frequency-related parameters.

    .. note::
        - To choose what Milankovitch cycles to be included in the frequency model, please edit use_precession, use_eccentricity, use_tilt.

  - **use_precession**: `bool`  
      Whether to include precession in the model (`True` or `False`).

  - **use_eccentricity**: `bool`  
      Whether to include eccentricity in the model (`True` or `False`).

  - **use_tilt**: `bool`  
      Whether to include tilt in the model (`True` or `False`).

    .. note::
        - The Earth's axial precession frequency (p0) and fundamental frequencies (gi and si) could be set fixed or not in the frequency model. 
          If fixed, input a specific value. If not fixed, input a range like [lower_bound, upper_bound] per each frequency value. If set as "default", the frequencies would be fixed using the values shown.
            
        - All the frequencies are in unit of "arcsec/yr".

  - **p0_values**: `float`, `list[float]`, or `"default"`  
      Single float value or list of two floats defining a range. If the p0 value is given as range, the tool will optimize that frequency.
      If `"default"`, value is set to **50.467718**.

  - **gi_values**: `list[float]`, `list[list[float]]`, or `"default"`  
      Frequencies representing `g1`, `g2`, `g3`, and `g4`.  
      Accepts:
      
      - A list of four float values  
      - A list of four lists (each with two floats) representing ranges. If the gi values are given as range, the tool will optimize these frequencies.
      - `"default"` (uses: **[5.579378, 7.456665, 17.366595, 17.910194]**)

  - **si_values**: `list[float]`, `list[list[float]]`, or `"default"`  
      Frequencies representing `s3` and `s5`.  
      Accepts:

      - A list of two float values  
      - A list of two lists (each with two floats) representing ranges. If the si values are given as range, the tool will optimize these frequencies.
      - `"default"` (uses: **[-18.845166, -17.758310]**)

  - **g5_value**: `float` or `"default"`  
    Constant value for frequency `g5`. Default: **4.257564**.

  - **s6_value**: `float` or `"default"`  
    Constant value for frequency `s6`. Default: **-26.347880**.
