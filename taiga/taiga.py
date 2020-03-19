# TaIGa - Taxonomic Information Gatherer
# Author: Maycon Douglas de Oliveira
# Year: 2020
# Version: 2.0
# License: MIT, check the LICENSE file in the root directory

from .common import helpers
from .core import taxonomy


def run():
    """Wrapper for the main execution of the program, performs a python version checking first

    Parameters:
    None

    Returns:
    None

    """

    helpers.sanitize_version()
    taxonomy.run_taiga()
