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

import argparse
import astrogeofit
import warnings

# Ignore all warnings
warnings.filterwarnings("ignore", message="Unable to import Axes3D.*")

# Set up the logger using the centralized config

logger = astrogeofit.setup_logger(False)


def path_results_choice(path_results: list):
    path_results_fit = None
    path_results_mcmc = None
    for path in path_results:
        if "fit.pickle" in path:
            path_results_fit = path
        elif "mcmc.pickle" in path:
            path_results_mcmc = path
    if path_results_fit == None and path_results_mcmc == None:
        if len(path_results) > 0:
            path_results_fit = path_results[0]
        if len(path_results) > 1:
            path_results_mcmc = path_results[1]
    return path_results_fit, path_results_mcmc


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ASTROGEOFIT TOOL")
    parser.add_argument(
        "-p",
        "--path",
        type=str,
        required=False,
        help="Path of configuration file. If nothing is introduced it will take the file found in config/config_file.",
    )
    parser.add_argument(
        "-fit",
        "--fit",
        action="store_true",
        required=False,
        help="Execution type. If *fit* is set, the tool executes *fit*",
    )
    parser.add_argument(
        "-sig",
        "--significance",
        action="store_true",
        required=False,
        help="Execution type. If *sig* is set, the tool executes *significance test*",
    )
    parser.add_argument(
        "-mcmc",
        "--mcmc",
        action="store_true",
        required=False,
        help="Execution type. If *mcmc* is set, the tool executes *mcmc*",
    )
    parser.add_argument(
        "-weights",
        "--weights",
        action="store_true",
        required=False,
        help="Execution type. If *weights* is set, the tool executes *weights*",
    )
    parser.add_argument(
        "-pr",
        "--pathresults",
        type=str,
        nargs="+",
        required=False,
        help="Path from where we obtain the fitting results. This argument is only needed if we execute the execution type: mcmc. Otherwise, it will not be used. In case nothing is introduced and the execution type is mcmc, it will be used the output path of the config file.",
    )
    parser.add_argument(
        "-s",
        "--savefigures",
        action="store_true",
        required=False,
        help="If set, the figures will be saved in the same folder as the results.",
    )
    parser.add_argument(
        "-r",
        "--remote",
        action="store_true",
        required=False,
        help="Set this option if the execution is done in a remote machine and you can't interact with it during the execution.",
    )
    
    parser.add_argument(
        "-logs",
        "--logs",
        action="store_true",
        required=False,
        help="Obtain the logs of the tool in a file.",
    )
    parser.add_argument(
        "-926",
        "--odp_926",
        action="store_true",
        required=False,
        help="Execute the ODP_926 example.",
    )
    parser.add_argument(
        "-1260",
        "--odp_1260",
        action="store_true",
        required=False,
        help="Execute the ODP_1260 example.",
    )
    parser.add_argument(
        "-1262",
        "--odp_1262",
        action="store_true",
        required=False,
        help="Execute the ODP_1262 example.",
    )
    parser.add_argument(
        "-basic",
        "--basic_example",
        action="store_true",
        required=False,
        help="Execute the basic example which uses synthetic data.",
    )
    args = parser.parse_args()

    if args.path == None:
        path = "default"
    else:
        path = args.path

    fit = True if args.fit else False
    significance = True if args.significance else False
    mcmc = True if args.mcmc else False
    weights = True if args.weights else False

    if (not fit and not significance and not mcmc and not weights) or (
        fit and significance and mcmc and weights
    ):
        exec_all = True
    else:
        exec_all = False
    if args.remote:
        remote = True
    else:
        remote = False

    # IF NOT PATHS ARE ENTERED, IT WILL USE THE SAME PATH AS IN THE configuration.yaml FILE.
    path_fit_results = None
    path_mcmc_results = None
    if not args.pathresults:
        path_results = None

    elif args.pathresults and not mcmc and weights:
        path_results = args.pathresults
        if isinstance(path_results, list) and len(path_results) != 2:
            logger.warning(
                f"{len(path_results)} have been introduced. The path or paths that have not been introduced will be trated as default paths (path of the configuration file)."
            )
            if not remote:
                default_paths = input("Do you wish to continue? (y/n) ")
                if default_paths == "n":
                    raise Exception
        logger.warning(
            "Make sure that the first path introduced is the fit results and the second is the mcmc results"
        )
        path_fit_results, path_mcmc_results = path_results_choice(path_results)

    elif args.pathresults and mcmc:
        path_results = args.pathresults
        if isinstance(path_results, list) and len(path_results) != 1:
            logger.warning(
                f"{len(path_results)} have been introduced. Only the path of the fitting results if needed. If ot has not been introduced, it will be trated as default path (path of the configuration file)."
            )
            if not remote:
                default_paths = input("Do you wish to continue? (y/n) ")
                if default_paths == "n":
                    raise Exception
        logger.warning(
            "Make sure that the first path introduced is the fit results and the second is the mcmc results"
        )
        path_fit_results, path_mcmc_results = path_results_choice(path_results)

    if args.savefigures:
        savefigs = True
    else:
        savefigs = False
    if args.logs:
        logs = True
    else:
        logs = False
    
    if args.odp_926:
        path = "./examples/ODP_926/configuration_file/configuration_ODP_926.yml"
    if args.odp_1260:
        path = "./examples/ODP_1260/configuration_file/configuration_ODP_1260.yml"
    if args.odp_1262:
        path = "./examples/ODP_1262/configuration_file/configuration_ODP_1262.yml"
    if args.basic_example:
        path = "./examples/synthetic_data/configuration_file/configuration_synthetic_data.yml"

    astrogeofit_class = astrogeofit.AstroGeoFit_tool(
        config_file_path= path,
        execute_all= exec_all,
        fit = fit,
        significance=significance,
        mcmc=mcmc,
        weight=weights,
        path_fit_results=path_fit_results,
        path_mcmc_results=path_mcmc_results,
        savefigs=savefigs,
        logs=logs,
        remote=remote,
        test=False
    )
    astrogeofit_class.run()
