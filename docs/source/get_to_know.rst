.. _get_started_agf:

Getting Started with `AstroGeoFit`
=================================

.. image:: https://img.shields.io/badge/python-3.11%2B-blue
    :alt: Python Version

.. image:: https://img.shields.io/badge/license-GPLv3-blue
    :alt: License

.. warning::
    - This section assumes that the previous sections :ref:`prerequisites`, :ref:`download_the_package` and :ref:`installation_updates` have been followed and the package is correctly installed in the computer of the user. If this has not been done, please refer to the attached sections.

In order to have an easy first contact with `AstroGeoFit`, we provide a small example which containes synthetic data, which we used to test the code while developing it.
Here you can find a step-by-step guide on how to execute this example and visualize the results.

Observe the Files Found on the Folder
-------------------------------------

The first important thing to do is to do a quick check on the files that the software uses. To do this please accede to the `examples/synthetic_data` folder.
There you will find three folders : `configuration_file`, `data` and `results`. 

The folder `configuration_file` contains all the needed variables that the tool uses as base configuration. For this example, nothing has to be modified, but if it can be used as 
"sandbox" to get a first idea of how the program uses the file. A full guide on what are these variables and how to set them can be found in :ref:`configuration_guide`.

The folder `data` containes the **data set** that will be used for this example as well as the **age-depth model** and the **eccentricity solution** used.

The folder `results` is empty and will be used to store the results of the execution.

Execute the Fitting and MCMC of the Example
-------------------------------------------

.. note::
    - The folder `examples` has to bee in the same directory as the `AstroGeoFit_tool.py` file. If not, this will not work.

To execute the basic example, first open a terminal (bash) in the main folder of the package, where you can find the `AstroGeoFit_tool.py` file.
Once it is open, execute the following command:

.. code-block:: bash

    python AstroGeoFit_tool.py --basic_example --fit --mcmc

This will execute the genetic algorithm fitting and the MCMC calcula using the synthetic data set. If you see the message:

.. code-block:: bash

    INFO -  All parameters obatined, starting analysis.

Means that the execution started correctly. Now a progress bar will appear.

Once the message

.. code-block:: bash

    INFO - Data updated and saved to ./examples/synthetic_data/results/syntehtic_data_results_mcmc.pickle

appears on the terminal, the execution of the fitting and the MCMC will be terminated.

Visualize the Results Obtained
------------------------------

To visualize the results obtained the notebooks provided can be used. The `AstroGeoFit_Analysis_Model_Fitting` notebook will display the results obtained by the genetic algorithm, while
the `AstroGeoFit_Analysis_MCMC` will display the results of the MCMC. To open these notebooks open a terminal in the main folder of `AstroGeoFit` and execute the folliwng command:

.. code-block:: bash

    jupyter notebook

When executing this, a navigator (chrome, mozilla, safari, etc.) should open. If not, in the terminal will appear a link similar to this: 

.. code-block:: bash
    
    http://localhost:8888/tree?token=bf9f95f8bb2d05ce5c7988287bd91a23e0e313d3ee2db8e2

Copy this link into a browser window and press enter.

Now you should see a folder with all the files that exist in the `AstroGeoFit` folder. Open the necessary notebook.

In order to obtain the correct results it is necessary to change the variable `configuration_file_path` found in the notebook in the following way:

.. code-block:: python

    configuration_file_path = "basic_example"

Once this is done, the notebook can be executed. In order to find a guide which explains the figures in detail, please refer to :ref:`visualization-results`

.. warning::
    - In the fitting notebook there is a plot that shows the **significance test** results. As we did not execute it this time, this will not work. It is best to not execute this cell.


Next Steps
----------

Now a first execution of the `AstroGeoFit` tool has been done it is time to explore in more detail the rest of the tool. The next step could be either **try to execute one 
of the examples found on the paper** which can be found in :ref:`execution_examples`, **learn how to adapt the configuration file to your data** by following the guide found
in :ref:`configuration_guide` or **learn the different execution options that AstroGeoFit proportionates** by following the guide found in :ref:`execution_page`.


.. footer:: Requests/Inquiries: astrogeofit@astrogeo.eu