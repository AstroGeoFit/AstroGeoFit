.. _tool_documentation:

AstroGeoFit Script Execution
=============================

.. image:: https://img.shields.io/badge/python-3.11%2B-blue
    :alt: Python Version

.. image:: https://img.shields.io/badge/license-GPLv3-blue
    :alt: License

.. _prerequisits:

Script File
-----------

To use the simplified version of the software, it is necessary to download the `AstroGeoFit_tool.py` file, 
which can be found in the package `astrogeofit_latest.zip`. To download the package check the section :ref:`download_the_package`

.. _execution:

Execution Of The Script
------------------------

To execute the **AstroGeoFit** tool, you need to run the script with the
appropriate arguments from the command line. Below are the options you
can use to configure and control the execution of the tool:

Command-Line Arguments
~~~~~~~~~~~~~~~~~~~~~~

To execute the program, the following line has to be executed (or copied) in to a terminal (or cmd in Windows). The argument `[options]` has to be changed  

.. code:: bash

    python AstroGeoFit_tool.py [options]

.. _execution-options:

.. _types_of_execution:

Options of Execution
~~~~~~~~~~~~~~~~~~~~

There are options that can be added to the previous command line, which will control the sections that `AstroGeoFit` will execute (*fit*, *mcmc*, *weights* or *significance test*).

.. warning::
    - If no types of execuion are set, it will execute all of the parts of the tool: `significance test`, `fit`, `mcmc` and `weights`.

-  | ``-fit``, ``--fit``

    | **(Optional)**

    | Define the execution type. If ``fit`` is set, the program will execute the **model fit**. This section inlcudes the fitting of the genetic algorithm and it's results.

    | Example:

    .. code:: bash

        python AstroGeoFit_tool.py -fit

-  | ``-mcmc``, ``--mcmc``

    | **(Optional)**

    | Define the execution type. If ``mcmc`` is set, the program will execute the ``mcmc`` . This section inlcudes the mcmc calcula and it's results.
    | Note that if the ``mcmc`` is executed without executing the ``fit`` on the same run, the path of the ``fit`` results has to be introduced. Refer to: :ref:`add_path_results`.
    | Example:

    .. code:: bash

        python AstroGeoFit_tool.py -mcmc

-  | ``-weights``, ``--weights``

    | **(Optional)**

    | Define the execution type. If ``weights`` is set, the program will compute the **weights**. This section inlcudes the MCMC weights computation and it's results. 
    | Note that if the ``weights`` is executed without executing the ``mcmc`` on the same run, the path of the ``mcmc`` results has to be introduced. Refer to: :ref:`add_path_results`.

    | **The execution of this option will affect the form of the plots.**
    | Example:

    .. code:: bash

        python AstroGeoFit_tool.py -weights

-  | ``-sig``, ``--significance``

    | **(Optional)**

    | Define the execution type. If ``sig`` is set, the program will execute the **significance test**.

.. note::
    The test applies only to the **fitting of the Genetic Algorithm** not to the MCMC.

    | Example:

    .. code:: bash

        python AstroGeoFit_tool.py -sig

.. _other_options:

Other Options
~~~~~~~~~~~~~

-  | ``-p``, ``--path``

    | **(Optional)**

    | Specify the path to the configuration file. If no path is provided, the default configuration file at ``config/config_file`` will beused.

    | Example:

    .. code:: bash

        python AstroGeoFit_tool.py -p /path/to/config_file.yml

.. _add_path_results:

-  | ``-pr``, ``--pathresults``

    | **(Optional)**

    | Specify the path where the file of the results is located. This argument is needed if the execution computes the *mcmc* but not the *fit*, or if it computes the *weights* but not the *mcmc* . If not provided, the path defined in the configuration file will be used.
    | It is possible to add the path for the *fit* results and *mcmc* results in case we want to execute the **weights** section alone. To do this, we just have to add the `-pr` option and add first the path of the *fit* results and after the path of the *mcmc* results separated by a space. **The order of the paths is important**.

    | Example:

    .. code:: bash

        python AstroGeoFit_tool.py -pr /path/of/the/results.pickle

-  | ``-r``, ``--remote``

    | **(Optional)**

    | This option is used when the tool is executed in a server that the user can not interact with. When activated, the software will not give the user the option to chose in certain situations, which can make
    | the tool more restrictive.

    | Example:

    .. code:: bash

        python AstroGeoFit_tool.py -r

-  | ``-logs``, ``--logs``

    | **(Optional)**

    | If this option is set, the tool will print the outputs of the execution in a txt file. This file will be saved under the name `logs` and will be included in the results folder.

    | Example:

    .. code:: bash

        python AstroGeoFit_tool.py -logs

-  | ``-s``, ``--savefigures``

    | **(Optional)**

    | If set, the tool will save the generated figures as jpegs in the same folder where the results are stored.

    | Example:

    .. code:: bash

        python AstroGeoFit_tool.py -s

.. _execution-examples:

Commonly Used Execution Examples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Run the tool with the default configuration file and execute all steps**:

    .. code:: bash

        python AstroGeoFit_tool.py

2. **Run the tool with a custom configuration file and execute the fitting process**:

    .. code:: bash

        python AstroGeoFit_tool.py -fit -p /path/to/config_file.yml

3. **Run the tool with a custom configuration file and execute the mcmc process (the folder of the fitting results will be extracted from the configuration file)**:

    .. code:: bash

        python AstroGeoFit_tool.py -mcmc -p /path/to/config_file.yml

4. **Run the tool with a custom configuration file and execute the significance test process**:

    .. code:: bash

        python AstroGeoFit_tool.py -sig -p /path/to/config_file.yml

5. **Run the tool with a custom configuration file and execute the weights calcula process (the folder of the fitting results and mcmc results will be extracted from the configuration file)**:

    .. code:: bash

        python AstroGeoFit_tool.py -weights -p /path/to/config_file.yml

Other Examples
~~~~~~~~~~~~~~

1. **Run the tool using MCMC and specify the path to the MCMC results**:

    .. code:: bash

        python AstroGeoFit_tool.py -mcmc -pr /path/to/fit_results.pickle

2. **Run the tool and save the figures**:

    .. code:: bash

        python AstroGeoFit_tool.py -s

3. **Run the tool selecting a custom configuration file, using the MCMC, specifying the path to the MCMC results and saving the figures**:

    .. code:: bash

        python AstroGeoFit_tool.py -p /path/to/config_file.yml -mcmc -pr /path/to/mcmc_results -s

4. **Run the tool using Weights section, and adding the paths of the fitting and the mcmc results**:

    .. code:: bash

        python AstroGeoFit_tool.py -weights -pr /path/to/fit_results /path/to/mcmc_results


