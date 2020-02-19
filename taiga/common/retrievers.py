import sys
import logging as log
from taiga.common import fetchers


def retrieve_from_taxid(taxon_list, user_email, retries):
    """Wrapper for the fetching functions that gets the name and genome id for each taxon from their
    TaxonID. It writes the fetched information directly in the Taxon objects provided.

    Parameters:
    taxon_list (list): A list containing all Taxon objects parsed from the input file.
    user_email (string): A valid email address provided for the user and used for Entrez.email.
    retries (int): The maximum number of retries after an unsuccsessful fetch attempt.

    Returns:
    None

    """

    for each_taxon in taxon_list:
        try:
            log.info('> Searching for corresponding organism name and genome '
                     'ID of "{}"'.format(each_taxon.taxon_id))

            fetchers.fetch_name_from_taxon_id(user_email, each_taxon, retries)

            fetchers.fetch_id_from_name(user_email, "genome", each_taxon)

            log.info(' >>>> Name and genome ID for "{}" : "{}" | "{}"\n'
                     .format(each_taxon.taxon_id, each_taxon.name, each_taxon.genome_id))
        except (KeyboardInterrupt):
            log.warning('\nQUIT: TaIGa was stopped by the user.\n')
            sys.exit()


def retrieve_from_names(taxon_list, user_email, correction, retries):
    """Wrapper for the fetching functions that gets the taxon id and genome id for each taxon from
    their name. It writes the fetched information directly in the Taxon objects provided.

    Parameters:
    taxon_list (list): A list containing all Taxon objects parsed from the input file.
    user_email (string): A valid email address provided for the user and used for Entrez.email.
    retries (int): The maximum number of retries after an unsuccsessful fetch attempt.
    correction (bool): Tells the function if it should run the spell correction function on the
                       taxon names before trying to fetch the other information.

    Returns:
    None

    """

    for each_taxon in taxon_list:
        try:
            if correction:
                log.info("> Correcting organism name of '{}'".format(each_taxon.name))

                fetchers.fetch_correct_spelling(user_email, each_taxon, retries)

            log.info("> Searching TaxID and genome ID of organism '{}'".format(each_taxon.name))

            fetchers.fetch_id_from_name(user_email, "taxonomy", each_taxon)

            log.info(" >>>> TaxID for '{}' : '{}'".format(each_taxon.name, each_taxon.taxon_id))

            fetchers.fetch_id_from_name(user_email, "genome", each_taxon)

            log.info(" >>>> Genome ID for '{}' : '{}'\n".format(each_taxon.name,
                                                                each_taxon.genome_id))
        except (KeyboardInterrupt):
            log.warning("\nQUIT: TaIGa was stopped by the user.\n")
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

    log.info(">> Gathering taxonomic information from NCBI Taxonomy. This might take a while.\n")

    for each_taxon in taxon_list:
        try:
            fetchers.fetch_taxonomic_info(user_email, each_taxon, retries)
        except (KeyboardInterrupt):
            log.warning("\nQUIT: TaIGa was stopped by the user.\n")
            sys.exit()

    log.info(">> Done gathering taxonomic information\n")
