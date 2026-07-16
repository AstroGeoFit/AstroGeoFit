.. _visualization-results:

Visualization of the Results Guide
==================================

.. toctree::
    :maxdepth: 2
    :hidden:
    
    figure_post_significance
    figures_post_fitting
    figures_post_mcmc

.. image:: https://img.shields.io/badge/python-3.11%2B-blue
    :alt: Python Version

.. image:: https://img.shields.io/badge/license-GPLv3-blue
    :alt: License

Methods of Visualization
------------------------

After running the study, there are two primary ways to visualize the results:

1. **Saving Figures Automatically**  
   Use the ``-s`` option when executing the tool (see :ref:`execution-options`). This saves all plots automatically; however, customization options (e.g., axis limits or number of genes) are not available.

2. **Using Python Notebooks**  
   Two dedicated notebooks allow interactive exploration of results from the genetic algorithm and MCMC analyses. These notebooks include plotting functions with editable parameters to fine-tune the visualizations.  
   A detailed guide for using and customizing these notebooks follows.

Installation and Running of the Notebooks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The notebooks can be launched from a terminal or used within an integrated development environment (IDE) such as Visual Studio Code. Instructions for different platforms are provided below.

Running Notebooks from the Terminal
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To run the Jupyter Notebooks from a terminal on **Linux, Windows, or macOS**, follow these steps:

1. **Ensure Jupyter is installed**  
   If Jupyter is not installed, install it with:

   .. code-block:: bash

      pip install notebook

2. **Navigate to the notebook directory**  

   .. code-block:: bash

      cd path/to/notebooks

3. **Launch Jupyter Notebook**

   .. code-block:: bash

      jupyter notebook

4. A web browser will open displaying the Jupyter Notebook interface. From here, select and open the desired notebook.

Using Notebooks in Visual Studio Code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Visual Studio Code has built-in support for Jupyter Notebooks through extensions. To use the notebooks within VS Code, follow these steps:

1. **Install VS Code** (if not already installed)  
   Download and install from `https://code.visualstudio.com/ <https://code.visualstudio.com/>`_.

2. **Install the Python extension**  
   - Open VS Code and access the Extensions Marketplace (``Ctrl+Shift+X`` or ``Cmd+Shift+X`` on macOS).  
   - Search for "Python" and install the official Microsoft extension.

3. **Ensure Jupyter is installed**  
   If not installed, run:

   .. code-block:: bash

      pip install jupyter

4. **Open the notebook in VS Code**  
   - Navigate to the folder containing the notebook (``.ipynb`` file).  
   - Open the notebook in VS Code.  
   - Use the "Run Cell" button to execute code interactively.

How to Use Jupyter Notebooks
----------------------------

Jupyter Notebooks provide an interactive Python environment. Basic operations include:

Executing a Cell
~~~~~~~~~~~~~~~~

- Select a code cell and press ``Ctrl + Enter`` to run it.

Modifying Variables and Re-executing Cells
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Change the value of a variable in a previously run cell, then rerun the cell with ``Ctrl + Enter``.

Adding a New Cell
~~~~~~~~~~~~~~~~~

- Click the ``+`` button in the toolbar to insert a new cell.

Restarting the Kernel
~~~~~~~~~~~~~~~~~~~~~

- To restart execution, click "Kernel" > "Restart Kernel" and rerun necessary cells.

By following these instructions, users on **Linux, Windows, or macOS** can efficiently utilize Jupyter Notebooks for visualization using this Python tool.

.. footer:: Requests/Inquiries: astrogeofit@astrogeo.eu