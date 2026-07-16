.. _download_the_package:

Download the `AstroGeoFit` Package
==================================

.. image:: https://img.shields.io/badge/python-3.11%2B-blue
    :alt: Python Version

.. image:: https://img.shields.io/badge/license-GPLv3-blue
    :alt: License

Once the package is downloaded from the link: https://www.astrogeo.eu/wp-content/astrogeofit/astrogeofit_latest.zip, the user will find a file named `astrogeofit_latest.zip`, which must be extracted.

After extraction, the folder will contain the following files:

File `astrogeofit-{version}-py3-none-any.whl`
---------------------------------------------

This file is the **Python package installer** and is used to install the software. To properly install the tool,
please refer to the guide in the section :ref:`installation_updates`.

File `configuration_file.yml`
-----------------------------

This file is used to define the configuration settings for the analysis the user intends to perform.

A detailed guide on how to set the values in this file can be found in the section :ref:`configuration_guide`.

File `AstroGeoFit_tool.py`
--------------------------

This script allows execution of the tool via shell (terminal) commands only. The necessary command descriptions
can be found in the section :ref:`tool_documentation`.

Folder *docs*
-------------

This folder contains the same documentation that is available on the **AstroGeo** webpage.
If offline access is required, the documentation can be viewed by opening the file `docs/index.html`.

File `AstroGeoFit_Analysis_Model_fitting.ipynb`
-----------------------------------------------

This file is a **Python notebook** used to visualize the results of the **significance test** and the **genetic algorithm**.

Descriptions of the generated figures can be found in:

- :ref:`figures_post_significance` for the significance test results.
- :ref:`figures_post_fitting` for the genetic algorithm results.

File `AstroGeoFit_Analysis_MCMC.ipynb`
--------------------------------------

This file is a **Python notebook** used to visualize the results of the **Bayesian MCMC** procedure.

A description of the figures included in the notebook is available in the section :ref:`figures_post_mcmc`.

Folder `Examples`
-----------------

On this folder, there can be found three subfolders that contain the required files to execute the analysis of the datasets analyzed on the paper
**AstroGeoFit: A Genetic Algorithm and Bayesian approach for the Astronomical Calibration of the Geological Timescale**, *Paleoceanography and Paleoclimatology*, 2025.

There is a detailed guide on how to execute and obtain the results of each of the datasets:

- :ref:`odp_926_example`
- :ref:`odp_1260_example`
- :ref:`odp_1262_example`

File `README.md`
----------------

This file contains documentation about the tool, including information on how to cite the related paper, a description of the software, and other relevant details.

File `LICENSE`
--------------

This file contains the license for the software.

.. footer:: Requests/Inquiries: astrogeofit@astrogeo.eu