import sys
import logging as log
from taiga.common import fetchers


def retrieve_from_taxid(taxon_list, user_email, retries):
    """Wrapper for the functions that get the name and Genome ID for each taxon from their Taxon ID.
    It writes the fetched information directly in the Taxon objects provided

    Parameters:
    taxon_list (list): A list containing all Taxon objects parsed from the input file
    user_email (string): A valid email address provided for the user and used for Entrez.email
    retries (int): The maximum number of retries after an unsuccsessful fetch attempt

    Returns:
    None

    """

    for each_taxon in taxon_list:
        try:
            log.info("\n> Searching for corresponding organism name and Genome ID of '{}'"
                     .format(each_taxon.taxon_id))

            fetchers.fetch_name_from_taxon_id(user_email, each_taxon, retries)

            log.info(" >>>> Name for '{}': '{}'".format(each_taxon.taxon_id, each_taxon.name))

            fetchers.fetch_id_from_name(user_email, "genome", each_taxon, retries)

            log.info(" >>>> Genome ID for '{}': '{}'".format(each_taxon.taxon_id,
                                                             each_taxon.genome_id))
        except (KeyboardInterrupt):
            log.warning("\nQUIT: TaIGa was stopped by the user\n")
            sys.exit()


def retrieve_from_names(taxon_list, user_email, correction, retries):
    """Wrapper for the functions that get the Taxon ID and Genome ID for each taxon from their name.
    It writes the fetched information directly in the Taxon objects provided

    Parameters:
    taxon_list (list): A list containing all Taxon objects parsed from the input file
    user_email (string): A valid email address provided for the user and used for Entrez.email
    retries (int): The maximum number of retries after an unsuccsessful fetch attempt
    correction (bool): Tells the function if it should run the spell correction function on the
                       taxon names before trying to fetch the other information

    Returns:
    None

    """

    for each_taxon in taxon_list:
        try:
            if correction:
                log.info("\n> Correcting organism name of '{}'".format(each_taxon.name))

                fetchers.fetch_correct_spelling(user_email, each_taxon, retries)

            log.info("\n> Searching Taxon ID and Genome ID of organism '{}'"
                     .format(each_taxon.name))

            fetchers.fetch_id_from_name(user_email, "taxonomy", each_taxon, retries)

            log.info(" >>>> Taxon ID for '{}': '{}'".format(each_taxon.name, each_taxon.taxon_id))

            fetchers.fetch_id_from_name(user_email, "genome", each_taxon, retries)

            log.info(" >>>> Genome ID for '{}': '{}'".format(each_taxon.name,
                                                              each_taxon.genome_id))
        except (KeyboardInterrupt):
            log.warning("\nQUIT: TaIGa was stopped by the user\n")
            sys.exit()


def retrieve_taxonomy(taxon_list, user_email, retries):
    """Wrapper for the function that fetches a taxon's taxonomic classification information

    Parameters:
    taxon_list (list): A list containing all Taxon objects parsed from the input file
    user_email (string): A valid email address provided for the user and used for Entrez.email
    retries (int): The maximum number of retries after an unsuccsessful fetch attempt

    Returns:
    None

    """

    log.info("\n> Gathering taxonomic information from NCBI Taxonomy. This might take a while...")

    for each_taxon in taxon_list:
        if not (each_taxon.missing_name or each_taxon.missing_taxon_id):
            try:
                fetchers.fetch_taxonomic_info(user_email, each_taxon, retries)
            except (KeyboardInterrupt):
                log.warning("\nQUIT: TaIGa was stopped by the user\n")
                sys.exit()
        elif each_taxon.missing_name:
            log.warning("\n>WARNING: Taxon {} is missing critical information. Skipping it..."
                        .format(each_taxon.taxon_id))
        elif each_taxon.missing_taxon_id:
            log.warning("\n>WARNING: Taxon {} is missing critical information. Skipping it..."
                        .format(each_taxon.name))

    log.info("\n> Done gathering taxonomic information")
