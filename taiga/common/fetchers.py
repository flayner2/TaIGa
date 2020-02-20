import sys
import logging as log
from Bio import Entrez
from urllib.error import HTTPError


def fetch_taxonomic_info(user_email, taxon, retries):
    """Receives a Taxon object and tries to fetch its full taxonomic classification information

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
    taxon.classification = dict()  # To avoid bugs, this has to be initialized as an empty dict

    try:
        query = Entrez.efetch(db="taxonomy", id=taxon.taxon_id, retmode="xml")
        parsed = Entrez.read(query)

        taxonomic_info = parsed[0]["LineageEx"]

        for taxon_level in taxonomic_info:
            if taxon_level["Rank"] == "no rank":
                if "no rank" not in list(taxon.classification):
                    taxon.classification[taxon_level["Rank"]] = []
                taxon.classification[taxon_level["Rank"]].append(taxon_level["ScientificName"])
            else:
                taxon.classification[taxon_level["Rank"]] = taxon_level["ScientificName"]
    except (KeyboardInterrupt):
        log.warning("\nQUIT: TaIGa was stopped by the user\n")
        sys.exit()
    except (HTTPError):
        log.error("\nERROR: Connection error. It could be due to a lack of internet connection "
                  "or a broken response from the NCBI servers. Try again a bit later\n")
        sys.exit()
    except (IndexError):
        log.error("\nERROR: Couldn't fetch taxonomic information for organism '{}'\n"
                  .format(taxon.name))
        sys.exit()
    except (Exception):
        log.error("\nERROR: Unknown error occurred while trying to fetch the taxonomic information "
                  "for organism '{}'. It could be due to TaIGa reaching the maximum number of "
                  "retries or issues with the NCBI servers. Maybe wait and try again a bit later\n"
                  .format(taxon.name))
        sys.exit()


def fetch_correct_spelling(user_email, taxon, retries):
    """Receives a Taxon object and tries to fetch a correct name for it from NCBI

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
            log.warning("\nWARNING: Couldn't find the correct organism name for '{}'"
                        .format(taxon.name))
            taxon.missing_corrected = True
        else:
            taxon.name = corrected_name
    except (KeyboardInterrupt):
        log.warning("\nQUIT: TaIGa was stopped by the user\n")
        sys.exit()
    except (RuntimeError):
        log.warning("\nWARNING: Couldn't find the correct organism name for '{}'"
                    .format(taxon.name))
        taxon.missing_corrected = True
    except (Exception):
        log.warning("\nWARNING: Unknown error occurred while trying to correct the spelling for "
                    "organism '{}'".format(taxon.name))
        taxon.missing_corrected = True


def fetch_id_from_name(user_email, db, taxon, retries):
    """Fetches either the Taxon ID or the Genome ID for a Taxon object, using the taxon's name

    Parameters:
    user_email (string): A valid email provided by the user and used for Entrez.email
    db (string): The database the function should try to fetch the information from. If the value
                 is 'genome', it will search Genome for a Genome ID. If it is 'taxonomy', it will
                 search Taxonomy for a Taxon ID
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
            log.warning("\nQUIT: TaIGa was stopped by the user\n")
            sys.exit()
        except (IndexError):
            log.warning("\nWARNING: Couldn't find a valid TaxID for the organism '{}'"
                        .format(taxon.name))
            taxon.missing_taxon_id = True
        except (Exception):
            log.warning("\nWARNING: Unknown error occurred while trying to find a valid TaxID for "
                        "organism '{}'".format(taxon.name))
            taxon.missing_taxon_id = True
    elif db == "genome":
        try:
            genome_id = parsed["IdList"][-1]
            taxon.genome_id = genome_id
        except (KeyboardInterrupt):
            log.warning("\nQUIT: TaIGa was stopped by the user\n")
            sys.exit()
        except (IndexError):
            log.warning("\nWARNING: Couldn't find a Genome ID for oragnism '{}'. It probably "
                        "doesn't have one available on NCBI".format(taxon.name))
        except (NameError):
            log.warning("\nWARNING: No Genome ID found for organism '{}'. TaIGa probably "
                        "didn't find the organism name for this Taxon ID"
                        .format(taxon.taxon_id))
        except (Exception):
            log.warning("\nWARNING: An unknown error occurred while searching Taxonomy for the "
                        "Genome ID of '{}'".format(taxon.taxon_id))


def fetch_name_from_taxon_id(user_email, taxon, retries):
    """Receives a Taxon object and tries to fetch a name for it from its Taxon ID

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
        log.warning("\nQUIT: TaIGa was stopped by the user\n")
        sys.exit()
    except (IndexError):
        log.warning("\nWARNING: Couldn't find an organism name for '{}'"
                    .format(taxon.taxon_id))
        taxon.missing_name = True
    except (Exception):
        log.warning("\nWARNING: An unknown error occurred while searching Taxonomy for '{}'"
                    .format(taxon.taxon_id))
        taxon.missing_name = True
