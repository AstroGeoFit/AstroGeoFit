.. _package_updates:

How to Upgrade the Software to the Latest Version
=================================================

.. image:: https://img.shields.io/badge/python-3.11%2B-blue
    :alt: Python Version

.. image:: https://img.shields.io/badge/license-GPLv3-blue
    :alt: License
    
To upgrade the software to the latest version, there are two available options:

Downloading the Package Again
-----------------------------

Re-downloading the package ensures that both the source code and all external files referenced in :ref:`installation_updates` are up to date.

After extracting the contents of the downloaded ZIP file, run the following command to upgrade the library:

.. code:: bash

    pip install --upgrade astrogeofit-{version}-py3-none-any.whl

.. note::

    Be sure to replace ``{version}`` with the actual version indicated in the filename.

Updating the Code Using PyPI
----------------------------

.. warning::

   This method only updates the core Python library. It does **not** update the additional required files and folders, such as configuration files, notebooks, or documentation.

To update the core code via PyPI, execute the following command:

.. code:: bash

    pip install --upgrade astrogeofit
