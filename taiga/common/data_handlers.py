import pandas as pd
import logging as log
import sys
import os


def create_df(taxon_list):
    # Done finishing the tax_info dictionary, print a message informing the next step
    log.info(">> Generating result table")

    # Headers set, will store all possible ranks from the input taxa
    raw_ranks = set()

    # The header set is a combination of the rank sets from all taxa
    for taxon in taxon_list:
        raw_ranks |= taxon.list_ranks()

    # To preserve the order, a list with all possible ranks from NCBI taxonomy is constructed
    ordered_ranks = [
        'no rank', 'superkingdom', 'kingdom', 'subkingdom', 'phylum', 'subphylum',
        'superclass', 'class', 'subclass', 'infraclass', 'cohort', 'subcohort',
        'superorder', 'order', 'suborder', 'infraorder', 'parvorder', 'superfamily',
        'family', 'subfamily', 'tribe', 'subtribe', 'genus', 'subgenus', 'section',
        'subsection', 'series', 'species-group', 'species', 'subspecies', 'forma'
    ]

    # To assure only ranks existing in the input taxa are part of the header, filter the list
    final_ranks = [rank for rank in ordered_ranks if rank in raw_ranks]

    # Add final header variables
    final_ranks.insert(0, "genome_id")
    final_ranks.insert(0, "taxon_id")

    # Create lists from the names and classifications of each taxon, in order
    taxon_names = [taxon.name for taxon in taxon_list]
    taxon_classification = [taxon.classification for taxon in taxon_list]

    # Create a dataframe from the lists of classifications, names and ranks
    frame = pd.DataFrame(taxon_classification, index=taxon_names, columns=final_ranks)

    # Add the values for taxon id and genome id for each taxon
    for taxon in taxon_list:
        frame.at[taxon.name, 'taxon_id'] = taxon.taxon_id
        frame.at[taxon.name, 'genome_id'] = taxon.genome_id

    # Convert the taxon id and genome id to integers
    frame.taxon_id = frame.taxon_id.astype(int)
    frame.genome_id = frame.genome_id.astype(int)

    # Change all 'NaN' occurences for 'N/A'
    frame.fillna("N/A", inplace=True)

    log.info("\n>> Done generating result table")

    return frame


def create_output(output_path, frame, taxon_list):
    log.info(
        "\n>> Checking if output folder exists")
    # Check if the output directory exists. If not, create it.
    if not os.path.exists(output_path):
        log.info("\n>> Creating output folder")
        try:  # Checking if the path for the output folder is valid
            os.makedirs(output_path)
        except (Exception):
            log.error(
                "\nERROR: Path to output folder may not be valid. Try again.")
            sys.exit()

    log.info(
        "\n>> Creating output file. You'll find it inside the provided output folder, named 'TaIGa_result.csv'"
    )
    # Export the DataFrame to the resulting .csv file
    frame.to_csv(output_path + 'TaIGa_result.csv')

    # Checking if there are missing correct names or TaxIDs. If there are, generating log files for those.
    log.info(
        "\n>> Creating a file for the organisms with missing information. You'll find it inside the provided output folder, named 'TaIGa_missing.txt'"
    )
    with open(output_path + 'TaIGa_missing.txt', 'w') as missing_file:
        missing_file.write("Missing corrected names: \n")
        for taxon in taxon_list:
            if taxon.missing_corrected:
                missing_file.write("\t\t\t{}\n".format(taxon.name))

        missing_file.write("Missing TaxID: \n")
        for taxon in taxon_list:
            if taxon.missing_taxon_id:
                missing_file.write("\t\t\t{}\n".format(taxon.taxon_id))

        missing_file.write("TaxIDs with missing names: \n")
        for taxon in taxon_list:
            if taxon.missing_name:
                missing_file.write("\t\t\t{}\n".format(taxon.taxon_id))
