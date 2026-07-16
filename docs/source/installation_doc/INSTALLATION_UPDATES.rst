.. _installation_updates:

Installation and Updates
=========================

.. image:: https://img.shields.io/badge/python-3.11%2B-blue
    :alt: Python Version

.. image:: https://img.shields.io/badge/license-GPLv3-blue
    :alt: License

.. toctree::
    :maxdepth: 1
    :hidden:
    
    updates

This page will guide you through setting up the software on your computer. The easiest and recommended way to obtain the software is by downloading the package directly from the official **AstroGeo** project webpage.

.. note::

   The software has been tested and can be installed on **Windows**, **Linux**, and **macOS** systems.

Installation from the ``AstroGeoFit`` Package
---------------------------------------------

.. note::
    Change the section `{version}` for the one found on the name of the file.

After downloading and extracting the archive from the official **AstroGeo** webpage, you will find a file named:

``astrogeofit-{version}-py3-none-any.whl``

This file is the core Python wheel used to install the library locally. To install it, open a terminal in the same directory as the wheel file and run the following command:

.. code:: bash

    pip install astrogeofit-{version}-py3-none-any.whl

Installation from PyPI (Python Package Index)
---------------------------------------------

Alternatively, the package can be installed from **PyPI** (the Python Package Index) using the following command:

.. code:: bash

    pip install astrogeofit

.. warning::

   If it is the first installation of the tool, this method is **not recommended** as it only installs the core library of the tool. This will be more useful for the actualizations of the code
   
   It does **not** include the following external files and folders:

    - ``AstroGeoFit_tool.py``
    - ``configuration_file.yml``
    - ``docs`` folder
    - ``AstroGeoFit_Analysis_Model_Fitting.ipynb``
    - ``AstroGeoFit_Analysis_MCMC.ipynb``
    - ``examples`` folder

   These files must be downloaded separately from the official **AstroGeo** webpage.

Update `AstroGeoFit`
--------------------

To update the code, follow the following guide :ref:`package_updates`.
