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

import dill  # type:ignore
import os
import yaml  # type:ignore
import pickle
from typing import Optional, Tuple
import yaml  # type: ignore
import pandas as pd  # type:ignore
import numpy as np  # type:ignore
from colorama import Fore, Style  # type: ignore
import logging

from astrogeofit.utils import exceptions

logger = logging.getLogger("ToolLogger")


def save_data(output_file_path: str, fit_mcmc_weights: Optional[str], **kwargs):
    if fit_mcmc_weights:
        if fit_mcmc_weights == "fit":
            name_res = f"{output_file_path}_fit.pickle"
        elif fit_mcmc_weights == "mcmc":
            name_res = f"{output_file_path}_mcmc.pickle"
        elif fit_mcmc_weights == "significance":
            name_res = f"{output_file_path}_significance.pickle"
        else:
            name_res = f"{output_file_path}_weights.pickle"
        # Load existing data if the file exists
    else:
        name_res = f"{output_file_path}.pickle"
    if os.path.exists(name_res):
        with open(name_res, "rb") as f:
            existing_data = dill.load(f)
    else:
        existing_data = {}

    # Update the existing data with the new data
    existing_data.update(kwargs)

    # Save the updated data back to the file
    with open(name_res, "wb") as f:
        dill.dump(existing_data, f, dill.HIGHEST_PROTOCOL)

    logger.info(f"Data updated and saved to {name_res}")


def save_analysis_info(**kwargs):
    kwargs_dict = kwargs
    output_folder = kwargs_dict["output_folder"]
    output_file = kwargs_dict["output_file_name"]
    filename = os.path.join(output_folder, f"analysis_parameters_{output_file}.yaml")
    # Ensure the output folder exists
    if os.path.exists(filename):
        os.remove(filename)
    os.makedirs(output_folder, exist_ok=True)
    with open(filename, "w") as file:
        yaml.dump(kwargs_dict, file, default_flow_style=False)
    logger.info(f"Analysis information saved to {filename}")


def read_yaml_file(file_path: str) -> Optional[dict]:
    with open(file_path, "r") as file:
        try:
            # Load the YAML file
            return yaml.safe_load(file)
        except yaml.YAMLError as exc:
            logger.error(f"Error reading YAML file: {exc}")
            return None


def data_read_dataset_from_file(
    data_path_values: dict,
    data_set:bool = True,
    age_depth_dataset:bool = False,
    eccentricity_dataset:bool = False
) -> tuple:
    """
    ## ADD SUMMARY

    Args:
        data_path_values (dict): Dictionary containing all the values found in the config file.

    Raises:
        ValueError: In case that the data path found in *data_path_values* does not exist, rise an Error.

    Returns:
        tuple: The first element is an array with the selected columns (found in the configuration file).
        The second element is the complete pandas dataFrame.
    """

    ########################  GET DATA FUNCTION ################################
    def get_data(
        data_path: str,
        depth_column_name: Optional[str] = None,
        proxy_column_name: Optional[str] = None,
        skiprows: Optional[int] = None,
        delimiter: Optional[str] = None,
        header: bool = True,
    ) -> Tuple[np.ndarray, pd.DataFrame]:
        """
        Reads data from a specified file and returns selected columns as a NumPy array and the DataFrame.

        Parameters:
        - data_path: str : Path to the data file.
        - depth_column_name: Optional[str] : Name of the column to be used as X. Defaults to the first column if None.
        - proxy_column_name: Optional[str] : Name of the column to be used as a proxy. Defaults to the second column if None.
        - skiprows: Optional[int] : Number of rows to skip at the beginning of the file.
        - header: bool : Indicates if the file contains a header row. Defaults to True.

        Returns:
        - Tuple: A tuple containing the selected columns as a NumPy array (transposed) and the DataFrame.

        Raises:
        - ValueError: If the file extension is unsupported or if specified columns are not found.
        """

        # ----- SKIPROWS ----- #
        # Default skiprows to 0 if None
        if skiprows is None:
            skiprows = 0

        # ----- SET DELIMITER IF NONE FOUND  ----- #
        _, file_extension = os.path.splitext(data_path)

        if delimiter is None:
            # Determine delimiter based on file extension
            if file_extension == ".csv":
                delimiter = ","  # Default CSV delimiter
            elif file_extension == ".dat":
                delimiter = " "  # Default DAT delimiter
            elif file_extension in [".txt", ".tab"]:
                delimiter = "\s"  # Default TXT delimiter

        # ----- READ DATA  ----- #
        if file_extension == ".xlsx":
            if header:
                dataFrame = pd.read_excel(data_path, skiprows=skiprows)
            else:
                dataFrame = pd.read_excel(data_path, skiprows=skiprows, header=None)
        else:
            # Read CSV, DAT, or TXT files
            if header:
                dataFrame = pd.read_csv(
                    data_path, delimiter=delimiter, skiprows=skiprows, engine="python"
                )
            else:
                dataFrame = pd.read_csv(
                    data_path,
                    delimiter=delimiter,
                    skiprows=skiprows,
                    header=None,
                    engine="python",
                )

        # ----- CHECK DATA SET FORMAT  ----- #
        if dataFrame.shape[1] < 2:
            raise exceptions.WrongDataSetStructure(
                "The data file must have at least two columns. If there are two columns, check that the delimiter set is the correct one."
            )

        # ----- OBTAIN COLUMNS BY NAMES ----- #
        if header:
            if (
                depth_column_name not in dataFrame.columns
                or proxy_column_name not in dataFrame.columns
            ):
                raise exceptions.WrongColumnNames(
                    message=f"One or both specified columns not found in the file: {depth_column_name}, {proxy_column_name}"
                )
            selected_columns = dataFrame[[depth_column_name, proxy_column_name]]

        # ----- OBTAIN COLUMNS BY INDEXES ----- #
        else:
            if depth_column_name is not None and proxy_column_name is not None:
                x_column_index = int(depth_column_name)
                proxy_column_index = int(proxy_column_name)
                if (
                    x_column_index >= dataFrame.shape[1]
                    or proxy_column_index >= dataFrame.shape[1]
                ):
                    raise exceptions.WrongColumnIndexSet(
                        f"One or both indexes are higher than the data set shape: Indexes set: {x_column_index}, {proxy_column_index} // Number of colums of the data: {dataFrame.shape[1]}"
                    )
            else:
                x_column_index = 0
                proxy_column_index = 1
            selected_columns = dataFrame.iloc[:, [x_column_index, proxy_column_index]]

        # Return the selected columns as a NumPy array (transposed) and the DataFrame
        return selected_columns.to_numpy().T, dataFrame

    ########################################################################################

    data_path = data_path_values["data_path"]
    logger.info(f"{Fore.GREEN} Obtaining data from:{Style.RESET_ALL} {data_path}")
    print("")
    if data_set:
        depth_column_name = data_path_values["depth_column_name"]
        proxy_column_name = data_path_values["proxy_column_name"]
    if age_depth_dataset:
        depth_column_name = data_path_values["depth_column_name"]
        proxy_column_name = data_path_values["age_column_name"]
    if eccentricity_dataset:
        depth_column_name = data_path_values["age_column_name"]
        proxy_column_name = data_path_values["eccentricity_column_name"]
    skiprows = data_path_values["skiprows"]
    delimiter = data_path_values["delimiter"]
    header = data_path_values["header"]
    return get_data(
        data_path=data_path,
        depth_column_name=depth_column_name,
        proxy_column_name=proxy_column_name,
        skiprows=skiprows,
        delimiter=delimiter,
        header=header,
    )


def load_dictionary_from_pickle(file_path: str) -> dict:
    if not os.path.exists(file_path):
        logger.error(f"File {file_path} does not exist.")
        return None
    try:
        with open(file_path, "rb") as file:
            return pickle.load(file)
    except (pickle.PickleError, EOFError) as e:
        logger.error(f"Error loading pickle file: {e}")
        return None
