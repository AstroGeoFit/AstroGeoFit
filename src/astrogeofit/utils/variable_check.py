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
import logging

from colorama import Fore, Style
from typing import Optional, Union

from astrogeofit.utils import exceptions

logger = logging.getLogger("ToolLogger")


def _data_set_variables_check(data_set_variables: dict, accepted_files: list) -> dict:
    # ----- ERRORS ---- #

    if not data_set_variables["data_path"]:
        message = "DATA_SET SECTION: You have not introduced any path for the data file. Please introduce a path and execute again."
        raise exceptions.NoPathIntroduced(message=message)

    data_file_path = data_set_variables["data_path"]
    if not os.path.exists(data_file_path):
        message = f"DATA_SET SECTION: No file in the path {data_file_path} found. Please check the path."
        raise exceptions.FileNotFound(message=message)

    _, file_extension = os.path.splitext(data_set_variables["data_path"])
    if file_extension not in accepted_files:
        message = f"DATA_SET SECTION: Type of file {file_extension} is not supported. Please add a file of type xlsx, csv, txt, dat or tab."
        raise exceptions.WrongFileType(message=message)

    if not data_set_variables["header"]:
        if (
            data_set_variables["depth_column_name"]
            or data_set_variables["proxy_column_name"]
        ):
            if isinstance(data_set_variables["depth_column_name"], str) or isinstance(
                data_set_variables["proxy_column_name"], str
            ):
                message = "DATA_SET SECTION: Header variable is set to false but after there are column names set. If there are column names in your data, please set the variable header to true. If not, please delete the column names."
                raise exceptions.NoHeaderButColumnNames(message=message)
            if (
                data_set_variables["depth_column_name"]
                and not data_set_variables["proxy_column_name"]
            ) or (
                not data_set_variables["depth_column_name"]
                and data_set_variables["proxy_column_name"]
            ):
                message = "DATA_SET SECTION: One of the index of the columns is set but not the other one. Please set both index."
                raise exceptions.NotBothColumnNamesSet(message=message)
    else:
        if (
            not data_set_variables["proxy_column_name"]
            or not data_set_variables["depth_column_name"]
        ):
            message = "DATA_SET SECTION: The header variable is set but one or both column names are not. Please set both names."
            raise exceptions.NotBothColumnNamesSet(message=message)

    # ----- WARNINGS ---- #
    warnings_found: dict = {}
    if not data_set_variables["delimiter"] and file_extension != ".xlsx":
        if file_extension == ".csv":
            delimiter = ","  # Default CSV delimiter
        elif file_extension == ".dat":
            delimiter = " "  # Default DAT delimiter
        else:
            delimiter = "\s"  # Default TXT delimiter
        warnings_found[
            "DATA SET WARNING DELIMITER"
        ] = f"No delimiter have been introduced for the data set. Because the file type is {file_extension}, the delimiter used will be {delimiter}."

    if not data_set_variables["header"]:
        if (
            not data_set_variables["depth_column_name"]
            and not data_set_variables["proxy_column_name"]
        ):
            warnings_found[
                "DATA SET WARNING INDEXES"
            ] = "No indexes were introduced for the data set file. The columns taken will be the column 0 as depth column and the column 1 as proxy column."
            
    return warnings_found


def _age_depth_model_check(
    age_model_data_variables: dict, accepted_files: list
) -> dict:
    # ----- ERRORS ---- #

    if age_model_data_variables["data_path"]:
        data_file_path = age_model_data_variables["data_path"]
        if not os.path.exists(data_file_path):
            message = f"AGE_DEPTH_MODEL_DATA SECTION: No file in the path {data_file_path} found. Please check the path."
            raise exceptions.FileNotFound(message=message)

        _, file_extension = os.path.splitext(age_model_data_variables["data_path"])
        if file_extension not in accepted_files:
            message = f"AGE_DEPTH_MODEL_DATA: Type of file {file_extension} is not supported. Please add a file of type xlsx, csv, txt, dat or tab."
            raise exceptions.WrongFileType(message=message)

        if not age_model_data_variables["header"]:
            if (
                age_model_data_variables["depth_column_name"]
                or age_model_data_variables["age_column_name"]
            ):
                if isinstance(
                    age_model_data_variables["depth_column_name"], str
                ) or isinstance(age_model_data_variables["age_column_name"], str):
                    message = "AGE_DEPTH_MODEL_DATA SECTION: Header variable is set to false but after there are column names set. If there are column names in your data, please set the variable header to true. If not, please delete the column names."
                    raise exceptions.NoHeaderButColumnNames(message=message)
                if (
                    age_model_data_variables["depth_column_name"]
                    and not age_model_data_variables["age_column_name"]
                ) or (
                    not age_model_data_variables["depth_column_name"]
                    and age_model_data_variables["age_column_name"]
                ):
                    message = "AGE_DEPTH_MODEL_DATA SECTION: One of the index of the columns is set but not the other one. Please set both index."
                    raise exceptions.NotBothColumnNamesSet(message=message)
        else:
            if (
                not age_model_data_variables["age_column_name"]
                or not age_model_data_variables["depth_column_name"]
            ):
                message = "AGE_DEPTH_MODEL_DATA SECTION: The header variable is set but one or both column names are not. Please set both names."
                raise exceptions.NotBothColumnNamesSet(message=message)

    # ----- WARNINGS ---- #
    warnings_found: dict = {}
    if not age_model_data_variables["data_path"]:
        warnings_found[
            "AGE MODEL WARNING PATH"
        ] = "No path for the age depth model was introduced."
    else:
        if not age_model_data_variables["delimiter"] and file_extension != ".xlsx":
            if file_extension == ".csv":
                delimiter = ","  # Default CSV delimiter
            elif file_extension == ".dat":
                delimiter = " "  # Default DAT delimiter
            else:
                delimiter = "\s"  # Default TXT delimiter
            warnings_found[
                "AGE MODEL WARNING DELIMITER"
            ] = f"No delimiter have been introduced for the data set. Because the file type is {file_extension}, the delimiter used will be {delimiter}."

        if not age_model_data_variables["header"]:
            if (
                not age_model_data_variables["depth_column_name"]
                and not age_model_data_variables["age_column_name"]
            ):
                warnings_found[
                    "AGE MODEL WARNING INDEXES"
                ] = "No indexes were introduced for the data set file. The columns taken will be the column 0 as depth column and the column 1 as age column."
    return warnings_found


def _eccentricity_model_check(
    eccentricity_model_data_variables: dict, accepted_files: list
) -> dict:
    # ----- ERRORS ---- #

    if eccentricity_model_data_variables["data_path"]:
        data_file_path = eccentricity_model_data_variables["data_path"]
        if not os.path.exists(data_file_path):
            message = f"ECCENTRICITY_SOLUTION_DATA SECTION: No file in the path {data_file_path} found. Please check the path."
            raise exceptions.FileNotFound(message=message)

        _, file_extension = os.path.splitext(
            eccentricity_model_data_variables["data_path"]
        )
        if file_extension not in accepted_files:
            message = f"ECCENTRICITY_SOLUTION_DATA: Type of file {file_extension} is not supported. Please add a file of type xlsx, csv, txt, dat or tab."
            raise exceptions.WrongFileType(message=message)

        if not eccentricity_model_data_variables["header"]:
            if (
                eccentricity_model_data_variables["age_column_name"]
                or eccentricity_model_data_variables["eccentricity_column_name"]
            ):
                if isinstance(
                    eccentricity_model_data_variables["age_column_name"], str
                ) or isinstance(
                    eccentricity_model_data_variables["eccentricity_column_name"], str
                ):
                    message = "ECCENTRICITY_SOLUTION_DATA SECTION: Header variable is set to false but after there are column names set. If there are column names in your data, please set the variable header to true. If not, please delete the column names."
                    raise exceptions.NoHeaderButColumnNames(message=message)
                if (
                    eccentricity_model_data_variables["age_column_name"]
                    and not eccentricity_model_data_variables["eccentricity_column_name"]
                ) or (
                    not eccentricity_model_data_variables["age_column_name"]
                    and eccentricity_model_data_variables["eccentricity_column_name"]
                ):
                    message = "ECCENTRICITY_SOLUTION_DATA SECTION: One of the index of the columns is set but not the other one. Please set both index."
                    raise exceptions.NotBothColumnNamesSet(message=message)
        else:
            if (
                not eccentricity_model_data_variables["eccentricity_column_name"]
                or not eccentricity_model_data_variables["age_column_name"]
            ):
                message = "ECCENTRICITY_SOLUTION_DATA SECTION: The header variable is set but one or both column names are not. Please set both names."
                raise exceptions.NotBothColumnNamesSet(message=message)

        if (
            not eccentricity_model_data_variables["start_time"]
            or not eccentricity_model_data_variables["final_time"]
        ):
            message = "ECCENTRICITY_SOLUTION_DATA SECTION: The start time or the final time were not introduced. Please introduce both."
            raise exceptions.NoTimeForEccentricitySolution(message)

    # ----- WARNINGS ---- #
    warnings_found: dict = {}
    if not eccentricity_model_data_variables["data_path"]:
        warnings_found[
            "ECCENTRICITY SOLUTION WARNING PATH"
        ] = "No path for the eccentricity solution was introduced. The plots about the correlation of the solution found with a given eccentricity solution will not be available."
    else:
        if (
            not eccentricity_model_data_variables["delimiter"]
            and file_extension != ".xlsx"
        ):
            if file_extension == ".csv":
                delimiter = ","  # Default CSV delimiter
            elif file_extension == ".dat":
                delimiter = " "  # Default DAT delimiter
            else:
                delimiter = "\s"  # Default TXT delimiter
            warnings_found[
                "ECCENTRICITY SOLUTION WARNING DELIMITER"
            ] = f"No delimiter have been introduced for the data set. Because the file type is {file_extension}, the delimiter used will be {delimiter}."

        if not eccentricity_model_data_variables["header"]:
            if (
                not eccentricity_model_data_variables["age_column_name"]
                and not eccentricity_model_data_variables["eccentricity_column_name"]
            ):
                warnings_found[
                    "ECCENTRICITY SOLUTION WARNING INDEXES"
                ] = "No indexes were introduced for the data set file. The columns taken will be the column 0 as age column and the column 1 as eccentricity column."
    return warnings_found


def frequency_values_check(frequency_values: dict) -> dict:
    # ----- ERRORS ---- #

    # ---- P0 VALUE ---- #
    p0 = frequency_values["p0_values"]
    gi = frequency_values["gi_values"]
    si = frequency_values["si_values"]

    use_prec = frequency_values["use_precession"]
    use_ecc = frequency_values["use_eccentricity"]
    use_tilt = frequency_values["use_tilt"]

    if use_prec or use_tilt:
        if not p0:
            var = "use_precession" if use_prec else "use_tilt"
            raise exceptions.WrongFrequencyValue(
                f"p0_values VARIABLE PROBLEM: The variable {var} is set to true but there are no values in p0_values."
            )

        if isinstance(p0, str):
            if p0 != "default":
                raise exceptions.WrongFrequencyValue(
                    "p0_values VARIABLE PROBLEM: String value must be 'default'."
                )

        elif isinstance(p0, list):
            if len(p0) != 2:
                raise exceptions.WrongFrequencyValue(
                    "p0_values VARIABLE PROBLEM: You must introduce two values to define the range of p0."
                )
            if not all(isinstance(x, (int, float)) for x in p0):
                raise exceptions.WrongFrequencyValue(
                    "p0_values VARIABLE PROBLEM: The list must contain two numerical values."
                )

        elif not isinstance(p0, (int, float)):
            raise exceptions.WrongFrequencyValue(
                "p0_values VARIABLE PROBLEM: p0 can only be a list of two values, a single numerical value, or 'default'."
            )


    # ---- GI VALUES ---- #
    if use_ecc or use_prec:
        if not gi:
            raise exceptions.WrongFrequencyValue(
                "gi_values VARIABLE PROBLEM: The variable use_eccentricity or use_precession is set to true but there are no values in gi_values."
            )

        if isinstance(gi, str):
            if gi != "default":
                raise exceptions.WrongFrequencyValue(
                    "gi_values VARIABLE PROBLEM: String value must be 'default'."
                )

        elif isinstance(gi, list):
            if len(gi) != 4:
                raise exceptions.WrongFrequencyValue(
                    "gi_values VARIABLE PROBLEM: The gi values variable has to be a list of four elements."
                )
            if not (
                all(isinstance(x, (int, float)) for x in gi)
                or all(
                    isinstance(sub, list)
                    and len(sub) == 2
                    and all(isinstance(x, (int, float)) for x in sub)
                    for sub in gi
                )
            ):
                raise exceptions.WrongFrequencyValue(
                    "gi_values VARIABLE PROBLEM: Must be either 4 lists of 2 values (int or float), or 4 int/float values."
                )


    # ---- SI VALUES ---- #
    if use_tilt:
        if not si:
            raise exceptions.WrongFrequencyValue(
                "si_values VARIABLE PROBLEM: The variable use_tilt is set to true but there are no values in si_values."
            )

        if isinstance(si, str):
            if si != "default":
                raise exceptions.WrongFrequencyValue(
                    "si_values VARIABLE PROBLEM: String value must be 'default'."
                )

        elif isinstance(si, list):
            if len(si) != 2:
                raise exceptions.WrongFrequencyValue(
                    "si_values VARIABLE PROBLEM: The si values variable has to be a list of two elements."
                )
            if not (
                all(isinstance(x, (int, float)) for x in si)
                or all(
                    isinstance(sub, list)
                    and len(sub) == 2
                    and all(isinstance(x, (int, float)) for x in sub)
                    for sub in si
                )
            ):
                raise exceptions.WrongFrequencyValue(
                    "si_values VARIABLE PROBLEM: Must be either 2 lists of 2 values (int or float), or 2 int/float values."
                )


    # ---- G5 VALUE ---- #
    g5_value = frequency_values["g5_value"]
    if (
        isinstance(g5_value, str)
        and g5_value != "default"
    ):
        message = "g5_value VARIABLE PROBLEM: String value must be 'default'."
        raise exceptions.WrongFrequencyValue(message)

    elif (
        not isinstance(g5_value, str)
        and not isinstance(g5_value, int)
        and not isinstance(g5_value, float)
    ):
        message = "g5_value VARIABLE PROBLEM: String value must be 'default'."
        raise exceptions.WrongFrequencyValue(message)

    # ---- S6 VALUE ---- #
    s6_value = frequency_values["s6_value"]
    if (
        isinstance(s6_value, str)
        and  s6_value != "default"
    ):
        message = "s6_value VARIABLE PROBLEM: String value must be 'default'."
        raise exceptions.WrongFrequencyValue(message)

    elif (
        not isinstance(s6_value, str)
        and not isinstance(s6_value, int)
        and not isinstance(s6_value, float)
    ):
        message = "s6_value VARIABLE PROBLEM: String value must be 'default'."
        raise exceptions.WrongFrequencyValue(message)

    # ---- WARNINGS ---- #
    warnings_found: dict = {}

    if not use_prec and not use_tilt and p0:
        warnings_found[
            "p0_value WARNING "
        ] = "The variable use_precession is set to false, but there are p0_values set. Those values will be ignored and precession will not be used."

    if not use_ecc and not use_prec and gi:
        warnings_found[
            "gi_values WARNING "
        ] = "The variable use_eccentricity is set to false, but there are gi_values set. Those values will be ignored and eccentricity will not be used."

    if not use_tilt and si:
        warnings_found[
            "si_values WARNING "
        ] = "The variable use_tilt is set to false, but there are si_values set. Those values will be ignored and tilt will not be used."

    return warnings_found


def interpolator_check(interpolator: str) -> None:
    # ---- ERROR ---- #
    possibile_interpolators = [
        "default",
        "CubicSpline",
        "Pchip",
        "BSpline",
        "AKIMA",
        "linear",
    ]
    if not interpolator in possibile_interpolators:
        message = f"interpolator VARIABLE PROBLEM: The choosen {interpolator} does not exist, please change the value for a valid one."
        raise exceptions.WrongInterpolator(message=message)


def seed_check(
    seed_values: Optional[Union[list, int]], number_algorithm_solutions: int
) -> dict:
    # ---- ERROR ---- #
    if (
        seed_values
        and not isinstance(seed_values, list)
        and not isinstance(seed_values, int)
    ):
        message = "seed VARIABLE PROBLEM: Seed can only be a list of integers, a integer or None."
        exceptions.WrongSeedValue(message)
    if isinstance(seed_values, list) and not all(
        isinstance(value, int) for value in seed_values
    ):
        message = "seed VARIABLE PROBLEM: Seed can only be a list of integers, a integer or None."
        exceptions.WrongSeedValue(message)
    # ---- WARNINGS ---- #
    warnings_found: dict = {}
    if isinstance(seed_values, int):
        warnings_found[
            "SEED SINGLE VALUE WARNING"
        ] = f"Seed has been set as a integer. The same seed ({seed_values}) will be used for every random function."
    if isinstance(seed_values, list) and len(seed_values) != number_algorithm_solutions:
        warnings_found[
            "NO MATCHING SEED VALUE WARNING"
        ] = f"The number of seeds given is different to the number of solutions required, some seed may be used twice."

    return warnings_found


def metric_value_check(metric_value: Optional[str]) -> None:
    # ---- ERRORS ---- #
    accepted_metrics = ["r2", "loglike"]
    if metric_value and not isinstance(metric_value, str):
        message = "metric VALUE PROBLEM: The metric introduce is not correct. It can only be a string or None."
        raise exceptions.WrongMetricValue(message=message)

    if metric_value and not metric_value in accepted_metrics:
        message = "metric VALUE PROBLEM: The metric value introduce is not correct. It can only be r2, loglike or None."
        raise exceptions.WrongMetricValue(message=message)


def population_value_check(population_value: int) -> None:
    # ---- ERRORS ---- #
    if not isinstance(population_value, int):
        raise exceptions.WrongPopulationValue()


def generations_value_check(generations_value: int) -> None:
    # ---- ERRORS ---- #
    if not isinstance(generations_value, int):
        raise exceptions.WrongGenerationsValue()


def jobs_value_check(jobs_value: int) -> None:
    # ---- ERRORS ---- #
    if not isinstance(jobs_value, int):
        raise exceptions.WrongJobsValue()


def algorithm_solutions_value_check(algorithm_solutions: int) -> None:
    # ---- ERRORS ---- #
    if not isinstance(algorithm_solutions, int):
        raise exceptions.WrongSolutionsValue()


def list_number_genes_value_check(list_number_genes: int) -> None:
    # ---- ERRORS ---- #
    if not isinstance(list_number_genes, list) and not all(
        isinstance(gen, int) for gen in list_number_genes
    ):
        raise exceptions.WrongListGensValue()


def mcmc_length_discard_thin_check(length: int, discard: int, thin: int) -> None:
    # ---- ERRORS ---- #

    if not isinstance(length, int):
        raise exceptions.WrongLengthValue()
    if not isinstance(discard, int):
        message = f"The value discard can only be an integer."
        raise exceptions.WrongDiscardValue(message=message)
    if not isinstance(thin, int):
        raise exceptions.WrongThinValue()
    if discard >= length:
        message = f"The value discard can not be higher or equal to the value of the length of the chains."
        raise exceptions.WrongDiscardValue(message)


def mcmc_solutions_check(mcmc_solutions: int, algorithm_solutions: int) -> None:
    if not isinstance(mcmc_solutions, int):
        message = "number_mcmc_solutions VARIABLE PROBLEM: The value number_mcmc_solutions has to be a integer."
        raise exceptions.WrongMCMCSolutionsValue(message=message)
    if mcmc_solutions > algorithm_solutions:
        message = "number_mcmc_solutions VARIABLE PROBLEM: The value number_mcmc_solutions can not be higher than the value number_algorithm_solutions."
        raise exceptions.WrongMCMCSolutionsValue(message=message)


def mcmc_list_of_genes_check(list_of_genes_mcmc: list, list_of_genes: list) -> None:
    if not isinstance(list_of_genes_mcmc, list) or not all(
        isinstance(value, int) for value in list_of_genes_mcmc
    ):
        message = "list_of_genes_mcmc VARIABLE PROBLEM: The variable list_of_genes_mcmc has to be a list of integers"
        raise exceptions.WrongConfiguration()
    for mcmc_gen in list_of_genes_mcmc:
        if mcmc_gen not in list_of_genes:
            message = f"The number of gen {mcmc_gen} found in the list_of_genes_mcmc does not exist in the list_number_genes."
            raise exceptions.WrongMCMCGenNumberValue(message=message)


def prior_distribution_frequencies_check(prior_distributions: dict, frequency_values: dict) -> dict:
    warnings_found = {}

    # ---- VALID KEYS AND VALUES ---- #
    valid_keys = {"p0_distribution", "gi_distribution", "si_distribution"}
    valid_values = {"gaussian", "uniform", None}

    # ---- KEY CHECK ---- #
    for key in prior_distributions:
        if key not in valid_keys:
            raise exceptions.WrongConfiguration(
                f"prior_distributions VALUE ERROR: The key '{key}' should not exist. "
                f"Accepted keys are: {list(valid_keys)}."
            )

    # ---- VALUE CHECK ---- #
    for key, value in prior_distributions.items():
        if value not in valid_values:
            raise exceptions.WrongConfiguration(
                f"prior_distributions VALUE ERROR: The value '{value}' for key '{key}' is not accepted. "
                f"Accepted values are: {list(valid_values)}."
            )

    # ---- WARNINGS ---- #
    p0 = frequency_values.get("p0_values")
    gi = frequency_values.get("gi_values")
    si = frequency_values.get("si_values")

    if prior_distributions.get("p0_distribution") and not (isinstance(p0, list) and len(p0) == 2):
        warnings_found["p0_distribution ADDED"] = (
            "p0_prior VARIABLE PROBLEM: There is a distribution added for p0_distribution, but p0_values is not a valid range (list of two values)."
        )

    if prior_distributions.get("gi_distribution"):
        if not gi or not isinstance(gi, list):
            warnings_found["gi_distribution ADDED"] = (
                "gi_prior VARIABLE PROBLEM:There is a distribution added for gi_distribution, but gi_values is not a list."
            )
        elif not all(isinstance(g, list) and len(g) == 2 for g in gi):
            warnings_found["gi_distribution ADDED"] = (
                "gi_prior VARIABLE PROBLEM: There is a distribution added for gi_distribution, but gi_values is not a list of four [min, max] pairs."
            )

    if prior_distributions.get("si_distribution"):
        if not si or not isinstance(si, list):
            warnings_found["si_distribution ADDED"] = (
                "si_prior VARIABLE PROBLEM: There is a distribution added for si_distribution, but si_values is not a list."
            )
        elif not all(isinstance(s, list) and len(s) == 2 for s in si):
            warnings_found["si_distribution ADDED"] = (
                "si_prior VARIABLE PROBLEM: There is a distribution added for si_distribution, but si_values is not a list of two [min, max] pairs."
            )
    use_prec = frequency_values.get("use_precession")
    use_ecc = frequency_values.get("use_eccentricity")
    use_tilt = frequency_values.get("use_tilt")

    if not use_prec and not use_tilt and prior_distributions.get("p0_distribution"):
        warnings_found["p0_distribution ADDED"] = (
                "p0_distribution VARIABLE PROBLEM: Neither the precession nor the tilt are being used, but still, there is a p0_distribution added. This distribution will be ignored."
            )
    if not use_ecc and not use_prec and prior_distributions.get("gi_distribution"):
        warnings_found["gi_distribution ADDED"] = (
                "gi_distribution VARIABLE PROBLEM: Neither the precession nor the eccentricity are being used, but still, there is a gi_distribution added. This distribution will be ignored."
            )
    if not use_tilt and prior_distributions.get("si_distribution"):
        warnings_found["si_distributionADDED"] = (
                "si_distribution VARIABLE PROBLEM: The obliquity is not being used, but still, there is a si_distribution added. This distribution will be ignored."
            )
    return warnings_found

def prior_parameter_frequencies_check(prior_frequencies: dict) -> None:
    accepted_keys = ["p0_prior", "gi_prior", "si_prior"]
    keys = prior_frequencies.keys()
    for key in keys:
        if key not in accepted_keys:
            message = f"prior_frequencies VARIABLE ERROR: The key {key} is not accepted. The accepted keys are: {accepted_keys}"

    if isinstance(prior_frequencies["p0_prior"], list):
        if len(prior_frequencies["p0_prior"]) != 2:
            message = "p0_prior VARIABLE PROBLEM: You must introduce two values to define the range of p0_prior."
            raise exceptions.WrongFrequencyValue(message)

        elif not all(
            isinstance(item, (int, float)) for item in prior_frequencies["p0_prior"]
        ):
            message = (
                "p0_prior VARIABLE PROBLEM: The list must contain two numerical values."
            )
            raise exceptions.WrongFrequencyValue(message)
    elif prior_frequencies["p0_prior"]:
        message = "p0_prior VARIABLE PROBLEM: You must introduce a list of two values to define the range of p0_prior."
        raise exceptions.WrongFrequencyValue(message)

    if isinstance(prior_frequencies["gi_prior"], list):
        if len(prior_frequencies["gi_prior"]) != 4:
            message = "gi_prior VARIABLE PROBLEM: The gi values variable has to be a list of four elements."
            raise exceptions.WrongFrequencyValue(message)
        if not (
            all(
                isinstance(sublist, list)
                and len(sublist) == 2
                and all(isinstance(item, (int, float)) for item in sublist)
                for sublist in prior_frequencies["gi_prior"]
            )
        ) and not (
            all(
                isinstance(item, (int, float)) for item in prior_frequencies["gi_prior"]
            )
        ):
            message = "gi_prior VARIABLE PROBLEM: Must be either 4 lists of 2 values (int or float), or 4 int/float values."
            raise exceptions.WrongFrequencyValue(message=message)
    elif prior_frequencies["gi_prior"]:
        message = "gi_prior VARIABLE PROBLEM: You must introduce four lists of two values to define the ranges of gi_prior."
        raise exceptions.WrongFrequencyValue(message)

    if isinstance(prior_frequencies["si_prior"], list):
        if len(prior_frequencies["si_prior"]) != 2:
            message = "si_prior VARIABLE PROBLEM: The si values variable has to be a list of two elements."
            raise exceptions.WrongFrequencyValue(message)
        if not (
            all(
                isinstance(sublist, list)
                and len(sublist) == 2
                and all(isinstance(item, (int, float)) for item in sublist)
                for sublist in prior_frequencies["si_prior"]
            )
        ) and not (
            all(
                isinstance(item, (int, float)) for item in prior_frequencies["si_prior"]
            )
        ):
            message = (
                "si_prior VARIABLE PROBLEM: Must be 2 lists of 2 values (int or float)."
            )
            raise exceptions.WrongFrequencyValue(message=message)
    elif prior_frequencies["si_prior"]:
        message = "si_prior VARIABLE PROBLEM: You must introduce two lists of two values to define the ranges of si_prior."
        raise exceptions.WrongFrequencyValue(message)



def configuration_file_variables_check(
    config_file: dict, 
    significance:bool,
    fit:bool, 
    mcmc:bool,
    weights:bool,
    remote: bool):
    warnings_list = []
    accepted_files = [".csv", ".dat", ".xlsx", ".tab", ".txt"]

    # ------ DATA SET VARIABLES CHECK ------ #
    data_set_variables = config_file["data_set"]
    data_set_warnings = _data_set_variables_check(
        data_set_variables=data_set_variables, accepted_files=accepted_files
    )
    if data_set_warnings != {}:
        warnings_list.append(data_set_warnings)

    # ------ AGE DEPTH MODEL VARIABLES CHECK ------ #
    if "age_depth_model_data" in config_file:
        age_model_data_variables = config_file["age_depth_model_data"]
        age_model_warnings = _age_depth_model_check(
            age_model_data_variables=age_model_data_variables, accepted_files=accepted_files
        )
        if age_model_warnings != {}:
            warnings_list.append(age_model_warnings)
    else:
        age_model_warnings = {"AGE_DEPTH_MODEL_DATA SECTION":"No Age-Depth model added."}
        warnings_list.append(age_model_warnings)

    # ------ ECCENTRICITY MODEL VARIABLES CHECK ------ #
    if "eccentricity_solution_data" in config_file:
        eccentricity_model_variables = config_file["eccentricity_solution_data"]
        eccentricity_model_warnings = _eccentricity_model_check(
            eccentricity_model_data_variables=eccentricity_model_variables,
            accepted_files=accepted_files,
        )
        if eccentricity_model_warnings != {}:
            warnings_list.append(eccentricity_model_warnings)
    else:
        eccentricity_model_warnings = {"ECCENTRICITY_SOLUTION_DATA SECTION":"eccentricity solution added."}
        warnings_list.append(eccentricity_model_warnings)
    # ---- FREQUENCY VALUES CHECK ---- #
    frequency_values = config_file["data_model_parameters"]["frequency_values"]
    frequency_warnings = frequency_values_check(frequency_values=frequency_values)
    if frequency_warnings != {}:
        warnings_list.append(frequency_warnings)

    if significance:
            # ---- INTERPOLATOR VALUE CHECK ---- #
        interpolator_check(config_file["significance_test_parameters"]["interpolator"])

        # ---- SEED VALUE CHECK ---- #
        seed_warning = seed_check(
            seed_values=config_file["significance_test_parameters"]["seed"],
            number_algorithm_solutions=config_file["significance_test_parameters"][
                "number_algorithm_executions"
            ],
        )
        if frequency_warnings != {}:
            warnings_list.append(seed_warning)

        # ---- METRIC VALUE CHECK ---- #
        metric_value_check(config_file["significance_test_parameters"]["metric_type"])

        # ---- POPULATION VALUE CHECK ---- #
        population_value_check(
            config_file["significance_test_parameters"]["population_size"]
        )

        # ---- GENERATIONS VALUE CHECK ---- #
        generations_value_check(
            config_file["significance_test_parameters"]["number_generations"]
        )

        # ---- JOBS VALUE CHECK ---- #
        jobs_value_check(
            config_file["significance_test_parameters"]["number_processors_used"]
        )

        # ---- SOLUTIONS VALUE CHECK ---- #
        algorithm_solutions_value_check(
            config_file["significance_test_parameters"]["number_algorithm_executions"]
        )

        # ---- LIST NUM GENES VALUE CHECK ---- #
        list_number_genes_value_check(
            config_file["significance_test_parameters"]["list_number_genes"]
        )

    if fit:
        # ---- INTERPOLATOR VALUE CHECK ---- #
        interpolator_check(config_file["genetic_algorithm_parameters"]["interpolator"])

        # ---- SEED VALUE CHECK ---- #
        seed_warning = seed_check(
            seed_values=config_file["genetic_algorithm_parameters"]["seed"],
            number_algorithm_solutions=config_file["genetic_algorithm_parameters"][
                "number_algorithm_solutions"
            ],
        )
        if frequency_warnings != {}:
            warnings_list.append(seed_warning)

        # ---- METRIC VALUE CHECK ---- #
        metric_value_check(config_file["genetic_algorithm_parameters"]["metric_type"])

        # ---- POPULATION VALUE CHECK ---- #
        population_value_check(
            config_file["genetic_algorithm_parameters"]["population_size"]
        )

        # ---- GENERATIONS VALUE CHECK ---- #
        generations_value_check(
            config_file["genetic_algorithm_parameters"]["number_generations"]
        )

        # ---- JOBS VALUE CHECK ---- #
        jobs_value_check(
            config_file["genetic_algorithm_parameters"]["number_processors_used"]
        )

        # ---- SOLUTIONS VALUE CHECK ---- #
        algorithm_solutions_value_check(
            config_file["genetic_algorithm_parameters"]["number_algorithm_solutions"]
        )

        # ---- LIST NUM GENES VALUE CHECK ---- #
        list_number_genes_value_check(
            config_file["genetic_algorithm_parameters"]["list_number_genes"]
        )
    if mcmc:
        # ---- MCMC PARAMETERS CHECK ---- #
        mcmc_length_discard_thin_check(
            length=config_file["mcmc_parameters"]["length_mcmc_chains"],
            discard=config_file["mcmc_parameters"]["discard"],
            thin=config_file["mcmc_parameters"]["thin"],
        )

        # ---- NUMBER MCMC SOLUTIONS CHECK ---- #
        mcmc_solutions_check(
            mcmc_solutions=config_file["mcmc_parameters"]["number_mcmc_solutions"],
            algorithm_solutions=config_file["genetic_algorithm_parameters"][
                "number_algorithm_solutions"
            ],
        )

        # ---- MCMC LIST OF GENS CHECK ---- #
        mcmc_list_of_genes_check(
            list_of_genes_mcmc=config_file["mcmc_parameters"]["list_of_genes_mcmc"],
            list_of_genes=config_file["genetic_algorithm_parameters"]["list_number_genes"],
        )

        # ---- PRIOR DISTRIBUTION CHECK ---- #
        prior_parameter_frequencies_check(
            prior_frequencies=config_file["mcmc_parameters"]["prior_frequencies"],
        )
        
        prior_distributions_warnings = prior_distribution_frequencies_check(
            prior_distributions=config_file["mcmc_parameters"]["prior_distributions"],
            frequency_values = frequency_values
        )
        
        if prior_distributions_warnings != {}:
            warnings_list.append(prior_distributions_warnings)

    print("")
    logger.info("NO ERRORS FOUND.")
    print("")
    if warnings_list == []:
        logger.info("NO WARNINGS FOUND")
        print("")
    else:
        logger.info("WARNINGS!")
        print("")
        for warning_dict in warnings_list:
            for key, value in warning_dict.items():
                logger.warning(f"{Fore.RED}{key}: {value}{Style.RESET_ALL}")
                print("")
        if remote:
            return
        continue_execution = input(
            f"{len(warnings_list)} WARNINGS FOUND. DO YOU WISH TO CONTINUE?(y/n) "
        )
        print("")
        if continue_execution != "y":
            raise exceptions.ExecutionStoppedByUser()
