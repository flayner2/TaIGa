import sys
import logging as log
from time import sleep
from random import randint
from taiga.common import fetchers

# TODO: refactor all these functions. They can probably be smaller and more efficient.
# TODO: also docstring them all and remove possibly redundant ones.
def retrieve_from_taxid(taxon_list, user_email, retries):
    i = 0

    for each_taxon in taxon_list:
        # Searching for each organism name associated to each TaxID
        try:
            log.info("> Searching for corresponding organism name and genome\
                     ID of '{}'".format(each_taxon.taxon_id))
            while i < (retries + 1):
                try:
                    # Getting the corresponding name for each TaxID
                    each_taxon.name = fetchers.get_name_from_id(user_email, each_taxon.taxon_id)
                    log.info(" >>>> Name for '{}' : '{}'\n"
                             .format(each_taxon.taxon_id, each_taxon.name))
                    break
                except (KeyboardInterrupt):
                    log.warning("\nQUIT: TaIGa was stopped by the user.\n")
                    sys.exit()
                except (IndexError):
                    # If the name isn't found, ignore this organism for now
                    log.warning(
                            "\n\t>> Couldn't find an organism name for '{}'\n"
                            .format(each_taxon.taxon_id))
                    each_taxon.missing_name = True
                    continue
                except (Exception):
                    log.warning(" * Error, trying again * ")
                    i += 1
                    sleep(randint(3, 8))
                    continue
        except (KeyboardInterrupt):
            log.warning("\nQUIT: TaIGa was stopped by the user.\n")
            sys.exit()
        except (Exception):
            log.warning("\nERROR: An unknown error occurred while searching Taxonomy for '{}'.\n"
                        .format(each_taxon.taxon_id))
        # Searching for each GenomeID associated to each name
        try:
            each_taxon.genome_id = search(user_email, "genome", each_taxon.name)
        except (KeyboardInterrupt):
            log.warning("\nQUIT: TaIGa was stopped by the user.\n")
            sys.exit()
        except (IndexError):
            each_taxon.genome_id = "N/A"
        except (NameError):
            log.warning("\nERROR: No genome ID found for organism '{}'. TaIGa probably didn't find\
                        the organism name for this TaxID.\n".format(each_taxon.taxon_id))
            continue
        except (Exception):
            log.warning("\nERROR: An unknown error occurred while searching Taxonomy for '{}'.\n"
                        .format(each_taxon.taxon_id))



def retrieve_from_multiple_names(input_names, user_email, c, retries):
    missing_corrected = []
    missing_taxid = []
    taxon_ids_collection = {}
    genome_ids_collection = {}

    i = 0

    # Using Python's collections module to uniquefy all names in list
    input_names = list(dict.fromkeys(input_names))

    for name in input_names:
        # Use Entrez.espell to correct the spelling of organism names
        if c:
            try:
                log.info("> Correcting organism name of '{}'".format(name))
                correct_name = fetchers.correct_spell(user_email, name)
                if len(correct_name) == 0:
                    log.warning(
                        "\n\t>> Couldn't find the correct organism name for '{}'\n"
                        .format(name))
                    missing_corrected.append(name)
            except (RuntimeError):
                log.warning(
                    "\n\t>> Couldn't find the correct organism name for '{}'\n"
                    .format(name))
                missing_corrected.append(name)
            except (KeyboardInterrupt):
                log.warning("\nQUIT: TaIGa was stopped by the user.\n")
                sys.exit()
            except (Exception):
                log.info(
                    "\n\t>> Unknown error occurred while trying to correct the spelling for organism '{}'.\n"
                    .format(name))
                missing_corrected.append(name)
        # If name correction is disabled, just use the input organism name itself
        else:
            correct_name = name

        # Search Taxonomy for the TaxID of the input organism names and Genome for the GenomeID
        log.info("> Searching TaxID and genome ID of organism '{}'"
                 .format(name))

        while i < (retries + 1):
            try:
                t_id = fetchers.search(user_email, "taxonomy", correct_name)
                log.info(" >>>> TaxID for '{}' : '{}'".format(correct_name,
                                                              t_id))
                try:
                    g_id = fetchers.search(user_email, "genome", correct_name)
                    log.info(" >>>> Genome ID for '{}' : '{}'\n".format(correct_name, g_id))
                except (IndexError):
                    g_id = '-'
                    log.warning(" >>>> Organism '{}' doesn't seem to have a valid genome ID, so assigning '-' as a placeholder.\n".format(correct_name))
                break
            except (KeyboardInterrupt):
                log.warning("\nQUIT: TaIGa was stopped by the user.\n")
                sys.exit()
            except (IndexError):
                log.warning(
                    "\n\t>> Couldn't find a valid TaxID for the organism '{}'.\n".
                    format(correct_name))
                missing_taxid.append(correct_name)
                break
            except (Exception):
                log.warning(" * Error, trying again * ")
                t_id = None
                i += 1
                sleep(randint(3, 8))
                continue

        if t_id is None:
            if len(correct_name) == 0:
                log.warning(
                    "\n\t>> Organism '{}' lacks a corrected name. Unable to search for TaxID.\n"
                    .format(name))
            else:
                log.warning(
                    "\n\t>> Unknown error occurred while trying to find a valid TaxID for organism '{}'.\n"
                    .format(correct_name))
            missing_taxid.append(correct_name)

        # Add the name, taxid and genomeid to the dictionary
        try:
            if correct_name not in missing_corrected and correct_name not in missing_taxid:
                taxon_ids_collection[correct_name] = t_id
                genome_ids_collection[correct_name] = g_id
        except (NameError):
            log.warning(
                "\n\t>> Will ignore organism '{}' for now. Try to handle it manually later.\n"
                .format(correct_name))
        except (KeyboardInterrupt):
            log.warning("\nQUIT: TaIGa was stopped by the user.\n")
            sys.exit()
        except (Exception):
            log.warning(
                "\n\t>> Unknown error occurred while trying to save the TaxID for organism '{}'.\n"
                .format(correct_name))

    return missing_corrected, missing_taxid, taxon_ids_collection, genome_ids_collection


def retrieve_from_single_name(input_names, user_email, c, retries):
    missing_corrected = []
    missing_taxid = []
    taxon_ids_collection = {}
    genome_ids_collection = {}

    i = 0

    # Use Entrez.espell to correct the spelling of organism names
    if c:
        try:
            log.info("> Correcting organism name of '{}'".format(input_names))
            correct_name = fetchers.correct_spell(user_email, input_names)
            if len(correct_name) == 0:
                log.warning(
                    "\n\t>> Couldn't find the correct organism name for '{}'\n"
                    .format(input_names))
                missing_corrected.append(input_names)
        except (RuntimeError):
            log.warning(
                "\n\t>> Couldn't find the correct organism name for '{}'\n".
                format(input_names))
            missing_corrected.append(input_names)
        except (KeyboardInterrupt):
            log.warning("\nQUIT: TaIGa was stopped by the user.\n")
            sys.exit()
        except (Exception):
            log.warning(
                "\n\t>> Unknown error occurred while trying to correct the spelling for organism '{}'\n"
                .format(input_names))
            missing_corrected.append(input_names)
    else: # If name correction is disabled, just use the input organism name itself
        correct_name = input_names

    # Search Taxonomy for the TaxID of the input organism names and Genome for the GenomeID
    log.info("> Searching TaxID and genome ID of organism '{}'".format(correct_name))

    while i < (retries + 1):
        try:
            t_id = fetchers.search(user_email, "taxonomy", correct_name)
            log.info(" >>>> TaxID for '{}' : '{}'".format(correct_name, t_id))
            try:
                g_id = fetchers.search(user_email, "genome", correct_name)
                log.info(" >>>> Genome ID for '{}' : '{}'\n".format(correct_name, g_id))
            except (IndexError):
                g_id = '-'
                log.warning(" >>>> Organism '{}' doesn't seem to have a valid genome ID, so assigning '-' as a placeholder.\n".format(correct_name))
            break
        except (KeyboardInterrupt):
            log.warning("\nQUIT: TaIGa was stopped by the user.\n")
            sys.exit()
        except (IndexError):
            log.warning(
                "\n\t>> Couldn't find a valid TaxID for the organism {}\n".format(
                    correct_name))
            missing_taxid.append(correct_name)
        except (Exception):
            log.warning(" * Error, trying again * ")
            t_id = None
            i += 1
            sleep(randint(3, 8))
            continue

    if t_id is None:
        if len(correct_name) == 0:
            log.warning(
                "\n\t>> Organism '{}' lacks a corrected name. Unable to search for TaxID.\n"
                .format(input_names))
        else:
            log.warning(
                "\n\t>> Unknown error occurred while trying to find a valid TaxID for organism '{}'.\n"
                .format(correct_name))
        missing_taxid.append(correct_name)
    # Add the name, taxid and genomeid to the dictionary
    try:
        if correct_name not in missing_corrected and correct_name not in missing_taxid:
            taxon_ids_collection[correct_name] = t_id
            genome_ids_collection[correct_name] = g_id
    except (NameError):
        log.warning(
            "\n\t>> Will ignore organism '{}' for now. Try to handle it manually later.\n"
            .format(correct_name))
    except (KeyboardInterrupt):
        log.error("\nQUIT: TaIGa was stopped by the user.\n")
        sys.exit()
    except (Exception):
        log.warning(
            "\n\t>> Unknown error occurred while trying to save the TaxID for organism '{}'.\n"
            .format(correct_name))

    return (missing_corrected, missing_taxid, taxon_ids_collection, genome_ids_collection)


def get_tax_info(user_email, taxon_ids_collection, retries):
    log.info(">> Gathering taxonomic information from NCBI Taxonomy. \
        This might take a while.\n")

    # Creating the list of dictionaries to store the taxonomic information for the organisms.
    tax_info = fetchers.organize_tax_info(user_email, taxon_ids_collection, retries)

    log.info(">> Done gathering taxonomic information\n")

    return tax_info
