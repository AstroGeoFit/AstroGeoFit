.. _prerequisites:

Pre-Requisites of the AstroGeoFit Package
=========================================

.. image:: https://img.shields.io/badge/python-3.11%2B-blue
    :alt: Python Version

.. image:: https://img.shields.io/badge/license-GPLv3-blue
    :alt: License

In order to install and execute ``AstroGeoFit``, there are some pre-requisites that have to be installed by the user.

These pre-requisites are listed below.

Python (version 3.9 or higher)
------------------------------

``AstroGeoFit`` uses **Python** as base language. In order to execute the tool it will be required to previously install *Python*, especifically the version 3.9 or any version higher.

To install **Python** refer to the official website https://www.python.org/downloads/.

To check which version of **Python** there is in the PC, execute the next command in the terminal:

    .. code:: bash

        python --version

Setting Up a Virtual Environment (Recommended)
----------------------------------------------

It’s recommended to install ``AstroGeoFit`` inside a virtual environment to
avoid dependency conflicts. You can use either **venv** or **conda**.

Using ``conda``
^^^^^^^^^^^^^^^

To install conda, follow the guide found in the official website: `Conda Installation Guide <https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html>`_

1. **Create a new Conda environment**:

    .. code:: bash

        conda create --name astrogeofit python=3.11

2. **Activate the environment**:

    .. code:: bash

        conda activate astrogeofit

Using ``venv`` (Built-in Python Virtual Environment)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can install and use `venv` by following the instructions in the `official Python venv documentation <https://docs.python.org/3/library/venv.html>`_.

1. **Create a virtual environment**:

    .. code:: bash

        python -m venv astrogeofit

2. **Activate the virtual environment**:

   -  **On Linux/macOS**:

        .. code:: bash

            source astrogeofit/bin/activate

   -  **On Windows (cmd or PowerShell)**:

        .. code:: powershell

            astrogeofit\Scripts\activate

.. note::
    
    - For the virtual environments, the name *astrogeofit* can be changed for the name that the user prefers.

.. footer:: Requests/Inquiries: astrogeofit@astrogeo.eu