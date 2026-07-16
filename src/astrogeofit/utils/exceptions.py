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

import logging
from colorama import Fore, Style

logger = logging.getLogger("ToolLogger")


class FileNotFound(Exception):
    def __init__(
        self,
        message: str = "",
    ) -> None:
        logger.error(f"{Fore.RED}{message}{Style.RESET_ALL}")
        super().__init__()


class NoPathIntroduced(Exception):
    def __init__(
        self,
        message: str = "",
    ) -> None:
        logger.error(f"{Fore.RED}{message}{Style.RESET_ALL}")
        super().__init__()


class WrongFileType(Exception):
    def __init__(self, message="") -> None:
        logger.error(f"{Fore.RED}{message}{Style.RESET_ALL}")
        super().__init__()


class NoHeaderButColumnNames:
    def __init__(
        self,
        message: str = "",
    ) -> None:
        logger.error(f"{Fore.RED}{message}{Style.RESET_ALL}")
        super().__init__()


class NotBothColumnNamesSet(Exception):
    def __init__(
        self,
        message: str = "",
    ) -> None:
        logger.error(f"{Fore.RED}{message}{Style.RESET_ALL}")
        super().__init__()


class WrongColumnIndexSet(Exception):
    def __init__(
        self,
        message: str = "",
    ) -> None:
        logger.error(f"{Fore.RED}{message}{Style.RESET_ALL}")
        super().__init__()


class WrongDataSetStructure(Exception):
    def __init__(self, message="") -> None:
        logger.error(f"{Fore.RED}{message}{Style.RESET_ALL}")
        super().__init__()


class WrongColumnNames(Exception):
    def __init__(self, message="") -> None:
        logger.error(f"{Fore.RED}{message}{Style.RESET_ALL}")
        super().__init__()


class ExecutionStoppedByUser(Exception):
    def __init__(self) -> None:
        logger.error(f"{Fore.RED}Execution Stopped By the User.{Style.RESET_ALL}")
        super().__init__()


class NoTimeForEccentricitySolution(Exception):
    def __init__(self, message="") -> None:
        logger.error(f"{Fore.RED}{message}{Style.RESET_ALL}")
        super().__init__()


class WrongFrequencyValue(Exception):
    def __init__(self, message="") -> None:
        logger.error(f"{Fore.RED}{message}{Style.RESET_ALL}")
        super().__init__()


class WrongInterpolator(Exception):
    def __init__(
        self,
        message="",
    ) -> None:
        logger.error(f"{Fore.RED}{message}{Style.RESET_ALL}")
        super().__init__()


class WrongSeedValue(Exception):
    def __init__(
        self,
        message="",
    ) -> None:
        logger.error(f"{Fore.RED}{message}{Style.RESET_ALL}")
        super().__init__()


class WrongMetricValue(Exception):
    def __init__(
        self,
        message="",
    ) -> None:
        logger.error(f"{Fore.RED}{message}{Style.RESET_ALL}")
        super().__init__(f"{Fore.RED}{message}{Style.RESET_ALL}")


class WrongPopulationValue(Exception):
    def __init__(
        self,
    ) -> None:
        message = f"population_size VARIABLE PROBLEM: The value population_size can only be an integer."
        logger.error(f"{Fore.RED}{message}{Style.RESET_ALL}")
        super().__init__()


class WrongGenerationsValue(Exception):
    def __init__(
        self,
    ) -> None:
        message = f"number_generations VARIABLE PROBLEM: The value number_generations can only be an integer."
        logger.error(f"{Fore.RED}{message}{Style.RESET_ALL}")
        super().__init__()


class WrongJobsValue(Exception):
    def __init__(
        self,
    ) -> None:
        message = f"number_processors_used VARIABLE PROBLEM: The value number_processors_used can only be an integer."
        logger.error(f"{Fore.RED}{message}{Style.RESET_ALL}")
        super().__init__()


class WrongSolutionsValue(Exception):
    def __init__(
        self,
    ) -> None:
        message = f"population_size VARIABLE PROBLEM: The value number_algorithm_solutions can only be ."
        logger.error(f"{Fore.RED}{message}{Style.RESET_ALL}")
        super().__init__()


class WrongListGensValue(Exception):
    def __init__(
        self,
    ) -> None:
        message = f"number_algorithm_solutions VARIABLE PROBLEM: The value list_number_genes can only be an list of integers."
        logger.error(f"{Fore.RED}{message}{Style.RESET_ALL}")
        super().__init__()


class WrongLengthValue(Exception):
    def __init__(
        self,
    ) -> None:
        message = f"length_mcmc_chains VARIABLE PROBLEM: The value length_mcmc_chains can only be an integer."
        logger.error(f"{Fore.RED}{message}{Style.RESET_ALL}")
        super().__init__()


class WrongDiscardValue(Exception):
    def __init__(
        self,
        message="",
    ) -> None:
        logger.error(f"{Fore.RED}{message}{Style.RESET_ALL}")
        super().__init__()


class WrongThinValue(Exception):
    def __init__(
        self,
    ) -> None:
        message = f"thin VARIABLE PROBLEM: The value thin can only be an integer."
        logger.error(f"{Fore.RED}{message}{Style.RESET_ALL}")
        super().__init__()


class WrongMCMCSolutionsValue(Exception):
    def __init__(
        self,
        message="",
    ) -> None:
        logger.error(f"{Fore.RED}{message}{Style.RESET_ALL}")
        super().__init__()


class WrongMCMCGenNumberValue(Exception):
    def __init__(
        self,
        message="",
    ) -> None:
        super().__init__(f"{Fore.RED}{message}{Style.RESET_ALL}")


class WrongPriorConfiguration(Exception):
    def __init__(
        self,
        message="",
    ) -> None:
        logger.error(f"{Fore.RED}{message}{Style.RESET_ALL}")
        super().__init__()


class WrongConfiguration(Exception):
    def __init__(
        self,
        message="",
    ) -> None:
        logger.error(f"{Fore.RED}{message}{Style.RESET_ALL}")
        super().__init__()


class WrongTimeValue(Exception):
    def __init__(self, timeValue) -> None:
        message = f"The value of time: {timeValue} is not correct. Please set a value greater than 0. "
        logger.error(f"{Fore.RED}{message}{Style.RESET_ALL}")
        super().__init__()


class WrongFrequencyConfig(Exception):
    def __init__(self) -> None:
        message = f"TThe frequency configuration is not correct. Please check the documentation and set a correct one."
        logger.error(f"{Fore.RED}{message}{Style.RESET_ALL}")
        super().__init__()
