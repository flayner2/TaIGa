# TaIGa - Taxonomic Information Gatherer
# Author: Maycon Douglas de Oliveira
# Year: 2020
# Version: 2.0
# License: MIT, check the LICENSE file in the root directory

import helpers
import taxonomy


def run() -> None:
    """Wrapper for the main execution of the program, performs a python version checking first

    Parameters:
    None

    Returns:
    None

    """

    helpers.sanitize_version()
    taxonomy.run_taiga()


if __name__ == "__main__":
    run()
