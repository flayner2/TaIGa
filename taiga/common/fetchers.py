from Bio import Entrez
import logging as log
import sys


def fetch_taxonomic_info(user_email, taxon, retries):
    """Receives a Taxon object and tries to fetch its full taxonomic classification information from
    NCBI's Taxonomy

    Parameters:
    user_email (string): A valid email provided by the user and used for Entrez.email
    taxon (object): A Taxon object that will hold the fetched information and provide the input
                    information
    retries (int): The maximum number of retries after an unsuccsessful fetch attempt

    Returns:
    None

    """

    Entrez.email = user_email
    Entrez.max_tries = retries
    Entrez.sleep_between_tries = 15

    try:
        query = Entrez.efetch(db="taxonomy", id=taxon.taxon_id, retmode="xml")
        parsed = Entrez.read(query)

        taxonomic_info = parsed[0]["LineageEx"]

        for taxon_level in taxonomic_info:
            taxon.classification[taxon_level["Rank"]] = taxon_level["ScientificName"]
    except (KeyboardInterrupt):
        log.warning("\nQUIT: TaIGa was stopped by the user.\n")
        sys.exit()
    except (IndexError):
        log.error("\nERROR: Couldn't fetch taxonomic information for organism '{}'"
                  .format(taxon.name))
    except (Exception):
        log.info("\n\t>> Unknown error occurred while trying to fetch the taxonomic information "
                 "for organism '{}'. It could be due to TaIGa reaching the maximum number "
                 "of retries or issues with the NCBI servers. Maybe wait and try again a bit "
                 "later.\n".format(taxon.name))


def fetch_correct_spelling(user_email, taxon, retries):
    """Receives a Taxon object and tries to fetch a correct name for it from NCBI.

    Parameters:
    user_email (string): A valid email provided by the user and used for Entrez.email
    taxon (object): A Taxon object that will hold the fetched information and provide the input
                    information
    retries (int): The maximum number of retries after an unsuccsessful fetch attempt

    Returns:
    None

    """

    Entrez.email = user_email
    Entrez.max_tries = retries
    Entrez.sleep_between_tries = 15

    try:
        query = Entrez.espell(db="taxonomy", term=taxon.name)
        parsed = Entrez.read(query)
        corrected_name = parsed["CorrectedQuery"]

        if (len(corrected_name) == 0):
            log.warning("\n\t>> Couldn't find the correct organism name for '{}'\n"
                        .format(taxon.name))
            taxon.missing_corrected = True
        else:
            taxon.name = corrected_name
    except (KeyboardInterrupt):
        log.warning("\nQUIT: TaIGa was stopped by the user.\n")
        sys.exit()
    except (RuntimeError):
        log.warning("\n\t>> Couldn't find the correct organism name for '{}'\n"
                    .format(taxon.name))
        taxon.missing_corrected = True
    except (Exception):
        log.info("\n\t>> Unknown error occurred while trying to correct the spelling for organism "
                 "'{}'.\n".format(taxon.name))
        taxon.missing_corrected = True


def fetch_id_from_name(user_email, db, taxon, retries):
    """Fetches either the taxon id or the genome id for a Taxon object, using the taxon's name

    Parameters:
    user_email (string): A valid email provided by the user and used for Entrez.email
    db (string): The database the function should try to fetch the information from. If the value
                 is 'genome', it will search Genome for a genome id. If it is 'taxonomy', it will
                 search Taxonomy for a taxon id
    taxon (object): A Taxon object that will hold the fetched information and provide the input
                    information
    retries (int): The maximum number of retries after an unsuccsessful fetch attempt

    Returns:
    None

    """

    Entrez.email = user_email
    Entrez.max_tries = retries
    Entrez.sleep_between_tries = 15

    query = Entrez.esearch(db=db, term=taxon.name, retmode="xml")
    parsed = Entrez.read(query)

    if db == "taxonomy":
        try:
            taxon_id = parsed["IdList"][0]
            taxon.taxon_id = taxon_id
        except (KeyboardInterrupt):
            log.warning("\nQUIT: TaIGa was stopped by the user.\n")
            sys.exit()
        except (IndexError):
            log.warning("\n\t>> Couldn't find a valid TaxID for the organism '{}'.\n"
                        .format(taxon.name))
            taxon.missing_taxon_id = True
        except (Exception):
            log.warning("\n\t>> Unknown error occurred while trying to find a valid TaxID for "
                        "organism '{}'.\n".format(taxon.name))
            taxon.missing_taxon_id = True
    elif db == "genome":
        try:
            genome_id = parsed["IdList"][-1]
            taxon.genome_id = genome_id
        except (KeyboardInterrupt):
            log.warning("\nQUIT: TaIGa was stopped by the user.\n")
            sys.exit()
        except (IndexError):
            log.warning('\nWARNING: Couldn\'t find a genome ID for oragnism "{}". It probably '
                        'doesn\'t have one available on NCBI.\n'.format(taxon.name))
        except (NameError):
            log.warning('\nERROR: No genome ID found for organism "{}". TaIGa probably '
                        'didn\'t find the organism name for this TaxID.\n'
                        .format(taxon.taxon_id))
        except (Exception):
            log.warning('\nERROR: An unknown error occurred while searching Taxonomy '
                        'for the genome id of "{}".\n'.format(taxon.taxon_id))


def fetch_name_from_taxon_id(user_email, taxon, retries):
    """Receives a Taxon object and tries to fetch a name for it from its taxon id

    Parameters:
    user_email (string): A valid email provided by the user and used for Entrez.email
    taxon (object): A Taxon object that will hold the fetched information and provide the input
                    information
    retries (int): The maximum number of retries after an unsuccsessful fetch attempt

    Returns:
    None

    """

    Entrez.email = user_email
    Entrez.max_tries = retries
    Entrez.sleep_between_tries = 15

    try:
        query = Entrez.efetch(db="taxonomy", id=taxon.taxon_id, retmode="xml")
        parsed = Entrez.read(query)
        name = parsed[0]["ScientificName"]

        taxon.name = name
    except (KeyboardInterrupt):
        log.warning('\nQUIT: TaIGa was stopped by the user.\n')
        sys.exit()
    except (IndexError):
        log.warning('\n\t>> Couldn\'t find an organism name for "{}"\n'
                    .format(taxon.taxon_id))
        taxon.missing_name = True
    except (Exception):
        log.warning('\nERROR: An unknown error occurred while searching Taxonomy for "{}".\n'
                    .format(taxon.taxon_id))
        taxon.missing_name = True
