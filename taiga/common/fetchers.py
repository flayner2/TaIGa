from Bio import Entrez
from time import sleep
from random import randint
import logging as log


def organize_tax_info(user_email, orgs_dict, retries):
    """ Creates a list of dictionaries for every organism from the input file,
    each dictionary containing all relevant output information for that
    organism. The 'taxonomy' key refers to a list o dictionaries, each
    containing 'rank': 'name' key:value pairs. Returns the full list. """
    taxonomic_info = []
    for org in list(orgs_dict.keys()):
        org_info = retrive_taxonomy(user_email, orgs_dict[org], retries)
        taxonomic_info.append({"name": org, "taxonomy": []})
        for taxon in org_info:
            taxons = {taxon["Rank"]: taxon["ScientificName"]}
            for item in taxonomic_info:
                if item["name"] == org:
                    item["taxonomy"].append(taxons)

    return taxonomic_info


def correct_spell(user_email, names):
    """ Uses Entrez.espell utility to correct the spelling of organism names.
    Returns the corrected name. """
    Entrez.email = user_email
    query = Entrez.espell(db="taxonomy", term=names)
    corrected_names = Entrez.read(query)

    return corrected_names["CorrectedQuery"]


def search(user_email, db, names):
    """ Uses Entrez.esearch to either search for oranisms TaxIDs or Genome IDs.
    Returns either one of those based on the 'db' argument. """
    Entrez.email = user_email
    query = Entrez.esearch(db=db, term=names, retmode="xml")
    parsed = Entrez.read(query)
    if db == "taxonomy":
        return parsed["IdList"][0]
    elif db == "genome":
        return parsed["IdList"][-1]


def retrive_taxonomy(user_email, tax_id, retries):
    """ Uses Entrez.efetch to fetch taxonomical data for the input taxon.
    Returns all taxonomical information for that taxon. """
    for i in range(retries):
        try:
            Entrez.email = user_email
            query = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
            tax_info = Entrez.read(query)

            return tax_info[0]["LineageEx"]
        except (IndexError):
            log.error("\nERROR: fetching for TaxID '{}' | retry: {} ".format(
                tax_id, i + 1))
            sleep(randint(0, 5))
            continue

    raise IndexError


def get_name_from_id(user_email, tax_id):
    Entrez.email = user_email
    query = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
    tax_info = Entrez.read(query)

    return tax_info[0]["ScientificName"]
