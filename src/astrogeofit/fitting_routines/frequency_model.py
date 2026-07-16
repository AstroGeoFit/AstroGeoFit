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

import numpy as np  # type:ignore
from typing import Optional
from functools import partial


def frequency_calcula(
    frequency_parameters: dict, default_values: dict, convers_factor: float
):
    use_prec = frequency_parameters["use_precession"]
    use_ecc = frequency_parameters["use_eccentricity"]
    use_tilt = frequency_parameters["use_tilt"]

    p0_values = frequency_parameters["p0_values"]
    eccentricity_values = frequency_parameters["gi_values"]
    tilt_values = frequency_parameters["si_values"]

    g5_value = frequency_parameters["g5_value"]
    s6_value = frequency_parameters["s6_value"]

    frequency_values_type, frequencies_values = frequency_types_and_values(
        use_prec=use_prec,
        use_ecc=use_ecc,
        use_tilt=use_tilt,
        p0_values=p0_values,
        eccentricity_values=eccentricity_values,
        tilt_values=tilt_values,
        g5_value=g5_value,
        s6_value=s6_value,
        convers_factor=convers_factor,
    )
    default_frequencies = get_default_frequencies(
        convers_factor=convers_factor, default_values=default_values
    )
    # WE CHANGE THE DEFAULT STRING VALUES FOR THE REAL DEFAULT VALUES.
    for key, value in frequencies_values.items():
        if isinstance(value, str):
            if value == "default":
                frequencies_values[key] = default_frequencies[key]
                
    range_free_frequencies = []
    arguments_for_model = {
        "use_precession": use_prec,
        "use_eccentricity": use_ecc,
        "use_tilt": use_tilt,
    }

    for key, value in frequency_values_type.items():
        if value == "range":
            range_free_frequencies.append(frequencies_values[key])
            arguments_for_model[key] = None
        elif value == "fixed":
            arguments_for_model[key] = frequencies_values[key]

    arguments_for_model["g5_value"] = frequencies_values["g5_value"]
    arguments_for_model["s6_value"] = frequencies_values["s6_value"]
    if range_free_frequencies == []:  
        return (None, partial(get_freqs_model, **arguments_for_model))
    else:
        return (
            np.vstack(range_free_frequencies),
            partial(get_freqs_model, **arguments_for_model),
        )


def frequency_types_and_values(
    use_prec: bool,
    use_ecc: bool,
    use_tilt: bool,
    p0_values,
    eccentricity_values,
    tilt_values,
    g5_value: float,
    s6_value: float,
    convers_factor: float,
):
    def obtain_precession_values(
        p0_values, frequencies_values_type, frequencies_user_values
    ):
        if isinstance(p0_values, str):
            frequencies_values_type["p0_values"] = "fixed"
            frequencies_user_values["p0_values"] = "default"
        elif isinstance(p0_values, list):
            frequencies_values_type["p0_values"] = "range"
            frequencies_user_values["p0_values"] = np.array(p0_values) / convers_factor
        elif isinstance(p0_values, int) or isinstance(p0_values, float):
            frequencies_values_type["p0_values"] = "fixed"
            frequencies_user_values["p0_values"] = p0_values / convers_factor

    def obtain_eccentricity_values(
        eccentricity_values, frequencies_values_type, frequencies_user_values
    ):
        if isinstance(eccentricity_values, str):
            frequencies_values_type["gi_values"] = "fixed"
            frequencies_user_values["gi_values"] = "default"
            return

        # Check if it's a list of 4 lists, each containing exactly 2 values of int or float
        elif isinstance(eccentricity_values, list):
            # Check if all elements are lists of two numbers (int or float)
            if all(isinstance(sublist, list) for sublist in eccentricity_values):
                frequencies_values_type["gi_values"] = "range"
                frequencies_user_values["gi_values"] = (
                    np.array(eccentricity_values) / convers_factor
                )
                return
            # Check if all elements are int or float
            elif all(isinstance(item, (int, float)) for item in eccentricity_values):
                frequencies_values_type["gi_values"] = "fixed"
                frequencies_user_values["gi_values"] = (
                    np.array(eccentricity_values) / convers_factor
                )
                return

    def obtain_tilt_values(
        tilt_values, frequencies_values_type, frequencies_user_values
    ):
        if isinstance(tilt_values, str):
            frequencies_values_type["si_values"] = "fixed"
            frequencies_user_values["si_values"] = "default"
            return

        # Check if it's a list of 4 lists, each containing exactly 2 values of int or float
        elif isinstance(tilt_values, list):
            # Check if all elements are lists of two numbers (int or float)
            if all(
                isinstance(sublist, list)
                and len(sublist) == 2
                and all(isinstance(item, (int, float)) for item in sublist)
                for sublist in tilt_values
            ):
                frequencies_values_type["si_values"] = "range"
                frequencies_user_values["si_values"] = (
                    np.array(tilt_values) / convers_factor
                )
                return
            # Check if all elements are int or float
            elif all(isinstance(item, (int, float)) for item in tilt_values):
                frequencies_values_type["si_values"] = "fixed"
                frequencies_user_values["si_values"] = (
                    np.array(tilt_values) / convers_factor
                )
                return

    frequencies_values_type = {}
    frequencies_user_values = {}
    # CHECK IF THE VARIABLES OF PRECESSION, ECCENTRICITY AND OBLIQUITY HAVE A CORRECT VALUE.
    if use_prec:
        obtain_precession_values(
            p0_values, frequencies_values_type, frequencies_user_values
        )
        if not use_ecc:
            obtain_eccentricity_values(
                eccentricity_values, frequencies_values_type, frequencies_user_values
            )

    if use_ecc:
        obtain_eccentricity_values(
            eccentricity_values, frequencies_values_type, frequencies_user_values
        )
    if use_tilt:
        if not use_prec:
            obtain_precession_values(
                p0_values, frequencies_values_type, frequencies_user_values
            )
        obtain_tilt_values(
            tilt_values, frequencies_values_type, frequencies_user_values
        )

    if isinstance(g5_value, str):
        if g5_value == "default":
            frequencies_user_values["g5_value"] = "default"
    elif isinstance(g5_value, int) or isinstance(g5_value, float):
        frequencies_user_values["g5_value"] = g5_value / convers_factor

    if isinstance(s6_value, str):
        if s6_value == "default":
            frequencies_user_values["s6_value"] = "default"
    elif isinstance(s6_value, int) or isinstance(s6_value, float):
        frequencies_user_values["s6_value"] = s6_value / convers_factor
    return frequencies_values_type, frequencies_user_values


def get_default_frequencies(convers_factor: float, default_values: dict) -> dict:
    default_variables_dict = {}
    for key, value in default_values.items():
        if key == "gi_values" or key == "si_values":
            default_variables_dict[key] = np.array(value) / convers_factor
        elif key == "p0_values":
            if isinstance(value, list):
                default_variables_dict[key] = np.array(value) / convers_factor
            else:
                default_variables_dict[key] = value / convers_factor
        else:
            default_variables_dict[key] = value / convers_factor
    return default_variables_dict


def get_freqs_model(
    free_frequencies: Optional[list],  # these are the free freq ranges
    use_precession: bool,
    use_eccentricity: bool,
    use_tilt: bool,
    p0_values=None,
    gi_values=None,
    si_values=None,
    g5_value=None,
    s6_value=None,
):
    i = 0
    if (p0_values is None) and (
        use_precession or use_tilt
    ):  # In case only ecc and til are needed.
        p0_values = free_frequencies[i]
        i += 1
    if (gi_values is None) and (use_precession or use_eccentricity):
        gi_values = free_frequencies[i : i + 4]
        i += 4
    if si_values is None and use_tilt:
        si_values = free_frequencies[i:]  # In case more si terms are needed.
    fs = []
    if use_precession:
        fs.append(np.hstack([gi_values + p0_values, g5_value + p0_values]))
    if use_eccentricity:
        fs.append(
            np.hstack(
                [
                    gi_values[1] - g5_value,
                    gi_values[3] - g5_value,
                    gi_values[3] - gi_values[1],
                    gi_values[2] - g5_value,
                    gi_values[2] - gi_values[1],
                ]
            )
        )
    if use_tilt:
        fs.append(np.hstack([si_values + p0_values, s6_value + p0_values]))
    fs = np.hstack(fs)
    return fs
