import pandas as pd
import logging as log
import sys
import os
from typing import List, Set


def create_df(taxon_list: List) -> pd.DataFrame:
    """Creates a Pandas DataFrame with the information for each input taxon

    Parameters:
    taxon_list (list): A list with all the Taxa the will provide the input information for the
                       DataFrame

    Returns:
    frame (DataFrame): A Pandas DataFrame object with all the taxonomic information organized for
                       each taxon

    """

    log.info("\n> Generating result DataFrame")

    # Headers set, will store all possible ranks from the input taxa
    raw_ranks: Set = set()

    # The header set is a combination of the rank sets from all taxa
    for taxon in taxon_list:
        if not (taxon.missing_name or taxon.missing_taxon_id or taxon.missing_classification):
            raw_ranks |= taxon.list_ranks()

    # To preserve the order, a list with all possible ranks from NCBI taxonomy is constructed
    ordered_ranks: List = [
        "no rank", "superkingdom", "kingdom", "subkingdom", "phylum", "subphylum",
        "superclass", "class", "subclass", "infraclass", "cohort", "subcohort",
        "superorder", "order", "suborder", "infraorder", "parvorder", "superfamily",
        "family", "subfamily", "tribe", "subtribe", "genus", "subgenus", "section",
        "subsection", "series", "species-group", "species", "subspecies", "forma"
    ]

    # To assure only ranks existing in the input taxa are part of the header, filter the list
    final_ranks: List = [rank for rank in ordered_ranks if rank in raw_ranks]

    # Add final header variables
    final_ranks.insert(0, "genome_id")
    final_ranks.insert(0, "taxon_id")

    # Create lists from the names and classifications of each taxon, in order
    taxon_names: List = [taxon.name for taxon in taxon_list if not (taxon.missing_name or
                                                                    taxon.missing_taxon_id or
                                                                    taxon.missing_classification)]
    taxon_classification: List = [taxon.classification for taxon in taxon_list
                                  if not (taxon.missing_name or taxon.missing_taxon_id
                                          or taxon.missing_classification)]

    # Create a dataframe from the lists of classifications, names and ranks
    frame: pd.DataFrame = pd.DataFrame(taxon_classification,
                                       index=taxon_names, columns=final_ranks)

    # Add the values for Taxon ID and Genome ID for each taxon
    for taxon in taxon_list:
        if not (taxon.missing_name or taxon.missing_taxon_id or taxon.missing_classification):
            frame.at[taxon.name, "taxon_id"] = taxon.taxon_id
            frame.at[taxon.name, "genome_id"] = taxon.genome_id

    # Convert the Taxon IDs to integers
    frame.taxon_id = frame.taxon_id.astype(int)

    # Change all 'NaN' occurences for 'N/A'
    frame.fillna("N/A", inplace=True)

    log.info("\n> Done generating result DataFrame")

    return frame


def create_output(output_path: str, frame: pd.DataFrame, taxon_list: List) -> None:
    """Creates the output directories (if they don't exist) and files for TaIGa

    Parameters:
    output_path (string): The path to the output folder as a string. It doesn't need to exist
    frame (DataFrame): The DataFrame containing all the information for all Taxa
    taxon_list (list): A list of Taxa objects used to create the file with the missing information

    Returns:
    None

    """
    log.info("\n> Checking if output folder exists")

    if not os.path.exists(output_path):
        log.info("\n> Creating output folder")

        try:
            os.makedirs(output_path)
        except (Exception):
            log.error(
                "\nERROR: Path to output folder may not be valid. Try again\n")
            sys.exit()

    log.info(
        "\n> Creating output file. You'll find it inside the provided output folder")

    # Export the DataFrame to the resulting .csv file
    frame.to_csv(output_path + "TaIGa_result.csv")

    log.info(
        "\n> Creating a file for the organisms with missing information. It might be empty")

    with open(output_path + "TaIGa_missing.txt", "w") as missing_file:
        missing_file.write("Missing corrected name: \n")
        for taxon in taxon_list:
            if taxon.missing_corrected:
                missing_file.write(f"\t{taxon.name}\n")

        missing_file.write("Missing Taxon ID: \n")
        for taxon in taxon_list:
            if taxon.missing_taxon_id:
                missing_file.write(f"\t{taxon.name}\n")

        missing_file.write("Missing name: \n")
        for taxon in taxon_list:
            if taxon.missing_name:
                missing_file.write(f"\t{taxon.taxon_id}\n")

        missing_file.write("Missing classification: \n")
        for taxon in taxon_list:
            if taxon.missing_classification:
                if not taxon.missing_name:
                    missing_file.write(f"\t{taxon.name}\n")
                else:
                    missing_file.write(f"\t{taxon.taxon_id}\n")
