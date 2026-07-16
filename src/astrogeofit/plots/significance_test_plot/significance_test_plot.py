"""
Copyright (C) 2025  CNRS
Lead author : J. Laskar
    
This file is part of the AstroGeoFit software.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import os
import math
import numpy as np
import scipy as sp
from typing import Optional
import logging

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator

from astrogeofit.main_routines import data_manager
from astrogeofit.utils import setup_logger
from astrogeofit.notebooks_code.notebooks_code import obtain_significance_test_file

logger = logging.getLogger("ToolLogger")
# Check if the logger already has handlers
if not logger.hasHandlers():
    logger = setup_logger.setup_logger(False)

def plt_significance_test_results(
    configuration_file_path:str,
    folder_with_the_results:Optional[str],
    file_with_the_results:Optional[str], 
    section_title: list,
    file_types_for_saving:Optional[list],
    output_path:Optional[str]
):
    """
    Plot the results of a significance test comparing data fits to autoregressive (AR) noise fits.

    This function generates a multi-panel figure showing histograms of $R^2_C$ values for data and AR model
    fits across different numbers of knots (model complexity). For each subplot:
      - Histograms for both data and noise fits are shown.
      - The mean $R^2_C$ value for the data is indicated.
      - A Gaussian fit is plotted over the AR model histogram.
      - The False Alarm Probability (FAP) is computed and displayed.

    Parameters:
        - **configuration_file_path**: str 
            Path to the configuration file used for the AstroGeoFit analysis.
        
        - **folder_with_the_results**: str or None
            Path to the folder containing the result files. If None, the default folder defined in the configuration file is used.
        
        - **file_with_the_results**: str or None
            Specific filename for the results. If None, a default filename is assumed.
        
        - **section_title**: str
            Title or header to display above the plot.
            
        - **types_of_saving** : list of str
            List of file formats to save the resulting plot (e.g., ["png", "pdf"]). 
            If empty, the figure will not be saved.

        - **output_path** : str
            Output path where the figure will be saved. If "default" or empty, the 
            `figures` subfolder inside the output folder defined in the configuration 
            file is used.

    Raises:
        FileNotFoundError: If the expected significance test results file is missing.
        KeyError: If required keys are not present in the input dictionary.

    Example:
        >>> configuration_file_path = "path/to/the/configuration_file.yaml"
        >>> folder_with_the_results = None
        >>> file_with_the_results = None
        >>> plt_significance_test_results(configuration_file_path, folder_with_the_results, file_with_the_results, "Significance Test Results")
        # Displays the figure with subplots for each number of knots and FAP annotations.
    """
    
    ########## LIBRARY MODE SECTION ##############

    tool_plt_significance_test_results(
        configuration_file_path=configuration_file_path,
        folder_with_the_results=folder_with_the_results,
        file_with_the_results=file_with_the_results,
        section_title=section_title,
        file_types_for_saving=file_types_for_saving,
        output_path=output_path
    )

def tool_plt_significance_test_results(
    configuration_file_path: str,
    folder_with_the_results: str,
    file_with_the_results: str,
    section_title: list,
    file_types_for_saving:Optional[list],
    output_path:Optional[str]
):
    supported_formats = ["pdf", "jpeg", "png", "eps"]
    def statistic_v3(x, x_0):
        return 2 * sp.stats.norm.sf(x.mean(), loc=x_0.mean(), scale=np.std(x_0, ddof=1))

    parameter_analysis = obtain_significance_test_file(
        config_file_path=configuration_file_path,
        output_folder=folder_with_the_results,
        output_file=file_with_the_results,
    )
    output_folder = parameter_analysis["output_folder"]
    file_results = parameter_analysis["output_file_name"]
    results_path_file = f"{output_folder}/{file_results}_significance.pickle"
    significance_test_results = data_manager.load_dictionary_from_pickle(
        results_path_file
    )
    list_num_genes = parameter_analysis["significance_test_parameters"][
        "list_number_genes"
    ]

    n_plots = len(list_num_genes)
    ncols = 4
    nrows = math.ceil(n_plots / ncols)  # Round up to ensure enough rows for n plots
    # Create figure and GridSpec (8 rows: 4 for captions, 4 for subplots)

    fig = plt.figure(figsize=(4 * ncols, 5 * nrows))
    gs = gridspec.GridSpec(nrows, ncols)  # Create grid for subplots

    fig.text(0.5, 1, section_title, ha="center", fontsize=14)  # , fontweight="bold")

    # Loop to create subplots in **odd rows (1,3,5,7)**

    data_result = significance_test_results["r2_data"]
    ar_result = significance_test_results["r2_ar"]
    for i in range(len(list_num_genes)):
        row = i // ncols  # Determine which row the subplot should go in
        col = i % ncols  # Determine the column in that row

        # Create a subplot in the correct location
        ax = fig.add_subplot(gs[row, col])
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
        ax.tick_params(axis="both", labelsize=14, width=2, length=5)
        h1_vals, _, _ = ax.hist(
            ar_result[:, i],
            bins=50,
            color="cornflowerblue",
            alpha=0.8,
            label="Noise fit",
            density=True,
        )
        h2_vals, _, _ = ax.hist(
            data_result[:, i],
            bins=50,
            color="black",
            alpha=0.45,
            label="Data fit",
            density=True,
        )
        mean = data_result[:, i].mean()
        ax.axvline(mean, color="firebrick", label="Data mean", lw=2)

        xmin = np.min(np.hstack((ar_result[:, i], data_result[:, i])))
        xmax = np.max(np.hstack((ar_result[:, i], data_result[:, i])))
        t = np.linspace(xmin, xmax, 500)
        plt.plot(
            t,
            sp.stats.norm.pdf(
                t, loc=ar_result[:, i].mean(), scale=ar_result[:, i].std(ddof=1)
            ),
            label="Gaussian fit",
            color="navy",
            lw=2,
        )

        p_value = (ar_result[:, i] >= mean).sum() / ar_result.shape[
            0
        ]  # This variable is not used.
        ax.set_title(f"$n={list_num_genes[i]}$ knots", fontsize=14)
        fap = statistic_v3(data_result[:, i], ar_result[:, i])
        ax.set_ylim(0, 1.15 * max(np.max(h1_vals), np.max(h2_vals)))
        fap_string = f"$FAP={fap:.4f}$" if fap > 1e-4 else f"$FAP={fap:.2e}$"
        ax.text(
            0.05,
            0.96,
            fap_string,
            verticalalignment="top",
            horizontalalignment="left",
            transform=ax.transAxes,
            size=12,
        )

        ax.set_xlabel("$R^2_C$", fontsize=10)
        if row == 0 and (col == 3 or col == len(list_num_genes) - 1):
            ax.legend(fontsize=12, loc="upper right")

    # Adjust layout for better spacing
    plt.tight_layout()
    if file_types_for_saving and file_types_for_saving != []:
        if not output_path or output_path == "default" or output_path == "":
            output_path = output_folder
            output_path = f"{output_path}/figures"
            os.makedirs(output_path, exist_ok=True)
        figname = f"Significance_Test_{section_title}"
        for fmt in file_types_for_saving:
            if fmt.lower() in supported_formats:
                if fmt.lower() == "tab":
                    logger.warning("tab format is not accepted for this plot.")
                else:
                    fig.savefig(
                        f"{output_path}/{figname}.{fmt.lower()}", format=fmt.lower()
                    )
                    logger.info(f"Saved {output_path}/{figname}.{fmt.lower()}")
            else:
                logger.warning(f"The format {fmt} is not accepted.")
    # plt.show()
