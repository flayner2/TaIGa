# TaIGa - Taxonomy Information Gatherer
# Author: Maycon D. Oliveira

import os
import sys
import argparse
import pandas as pd
from Bio import Entrez
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def correct_spell(user_email, names):
    Entrez.email = user_email
    query = Entrez.espell(db="taxonomy", term=names)
    corrected_names = Entrez.read(query)

    return corrected_names["CorrectedQuery"]


def search(user_email, names):
    Entrez.email = user_email
    query = Entrez.esearch(db="taxonomy", term=names, retmode="xml")
    tax_ids = Entrez.read(query)

    return tax_ids["IdList"][0]


def retrive_taxonomy(user_email, tax_id):
    Entrez.email = user_email
    query = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
    tax_info = Entrez.read(query)

    return tax_info[0]["LineageEx"]


# Taking the user input from the command line. The path to the input file and the path to the output folder. Both are taken from stdin. Aditional options let TaIGa understand the format of the input.
taiga = argparse.ArgumentParser(
    description=
    "TaIGa retrieves metadata of an organism (or a collection of organisms) from a Genbank format genome file."
)
taiga.add_argument(
    "input",
    help=
    "Full path to input file (or only the file name if it is on the same folder)"
)
taiga.add_argument(
    "outdir",
    help=
    "Full path to output folder (TaIGa will create subfolders and output files automatically)"
)
taiga.add_argument(
    "email",
    help=
    "Any valid email of yours. It is a common practice when using E-utils, which TaIGa does use"
)
taiga.add_argument(
    "--single",
    help="Use this if the input genome file contains a single record.",
    action="store_true")
taiga.add_argument(
    "--same",
    help=
    "Use this if the input genome file contains multiple records but all for the same organism (eg. multiple scaffolds for the same 'Apis mellifera' genome). Don't use this if your input file has multiple records for different organisms, even if one or another is repeated (eg. Genbank file with two 'Apis mellifera' records, one 'Bombus impatiens' record and one 'Homo sapiens' record).",
    action="store_true")

# Parsing the program arguments
args = taiga.parse_args()

print("""*********************************************
*                                           *
*   TaIGa - Taxonomy Information Gatherer   *
*                                           *
*********************************************""")

# Checking the type of input Genbank file
if args.single:  # Single record on Genbank file
    try:
        print("\n> Parsing input as a Genbank file with a single record.\n")
        try:
            records = SeqIO.read(args.input, "genbank")
        except (KeyboardInterrupt):
            print("\nQUIT: TaIGa was stopped by the user.\n")
        except:
            print(
                "\nERROR: Couldn't parse input genome file. Check your file and try running TaIGa again.\n"
            )
        names = records.annotations[
            "organism"]  # Take the name of the organism on the input file
        if names:
            print("> '{}' ---> All OK".format(names))
        else:
            print(
                "\nERROR: Something went wrong while trying to parse the input file. Check the file and the program execution commands and try again.\n"
            )
    except (KeyboardInterrupt):
        print("\nQUIT: TaIGa was stopped by the user.\n")
    # Catch an error if the file contains more than one record
    except (ValueError):
        print(
            "\nERROR: The file contains more than one record. Try running TaIGa again without the '--single' option.\n"
        )
elif args.same:  # Multiple records from the same organism
    print(
        "\n> Parsing input as a Genbank file with multiple records for the same organism.\n"
    )
    try:
        records = list(SeqIO.parse(args.input, "genbank"))
    except (KeyboardInterrupt):
        print("\nQUIT: TaIGa was stopped by the user.\n")
    except:
        print(
            "\nERROR: Couldn't parse input genome file. Check your file and try running TaIGa again.\n"
        )
    names = records[0].annotations[
        "organism"]  # Only take the first taxon name retrieved, as all are from the same organism
    if names:
        print("> '{}' ---> All OK".format(names))
    else:
        print(
            "\nERROR: Something went wrong while trying to parse the input file. Check the file and the program execution commands and try again.\n"
        )
else:  # Multiple records from multiple organisms
    name_counter = 0
    print(
        "\n> Parsing input as a Genbank file with multiple records for multiple organisms.\n"
    )
    try:
        records = SeqIO.parse(args.input, "genbank")
    except (KeyboardInterrupt):
        print("\nQUIT: TaIGa was stopped by the user.\n")
    except:
        print(
            "\nERROR: Couldn't parse input genome file. Check your file and try running TaIGa again.\n"
        )
    names = [seq.annotations["organism"]
             for seq in records]  # List the names of all organisms
    if names:
        print(
            "\n> All OK with parsing the input file. Checking the records...\n"
        )
        for name in names:
            name_counter += 1
            if name:
                print("> '{}' ---> All OK".format(name))
            else:
                print(
                    "> Something went wrong while trying to retrieve the organism name of record number '{}'"
                    .format(name_counter))
                continue
    else:
        print(
            "\nERROR: Something went wrong while trying to parse the input file. Check the file and the program execution commands and try again.\n"
        )

# Define the user email for Entrez.email, as it is a good practice of E-utils
user_email = args.email
# A dictionary to hold {'organism name': taxid} key:value pairs
taxon_ids_collection = {}
# A list to hold all organisms for which a correct name wasn't found
missing_corrected = []
# A list to hold all organisms for which the TaxID wasn't found
missing_taxid = []
# A dictionary of {'orignal name': correct name} key:value pairs
original_names = {}
# A list of dictionaries to hold the taxonomic information for each organism. The keys are the corrected organism names, and the values are lists of dictionaries, each containing the {taxonomic level: taxon} key:value pair
tax_info = []

print("\n> Searching for taxonomic information...\n")

# Checking if there are multiple records on the file (for different organisms)
if type(names) == list:
    for name in names:
        # Use Entrez.espell to correct the spelling of organism names
        try:
            correct_name = correct_spell(user_email, name)
            print("> Correcting organism name of '{}' to: '{}'".format(
                name, correct_name))
            original_names[name] = correct_name
        except (RuntimeError):
            pass
            print("> Couldn't find the correct organism name for '{}'".format(
                name))
            missing_corrected.append(name)
        except (KeyboardInterrupt):
            print("\nQUIT: TaIGa was stopped by the user.\n")
        except:
            pass
            print(
                "> Unknown error occurred while trying to correct the spelling for organism '{}'."
                .format(name))
            missing_corrected.append(name)
        # Search Taxonomy for the TaxID of the input organism names
        try:
            t_id = search(user_email, correct_name)
            print("> TaxID of organism '{}': '{}'".format(name, t_id))
        except (IndexError):
            pass
            print(
                "> Couldn't find a valid TaxID for the organism '{}'.".format(
                    correct_name))
            missing_taxid.append(correct_name)
        except (KeyboardInterrupt):
            print("\nQUIT: TaIGa was stopped by the user.\n")
        except:
            pass
            print(
                "> Unknown error occurred while trying to find a valid TaxID for organism '{}'."
                .format(correct_name))
            missing_taxid.append(correct_name)
        # Add the name and taxid to the dictionary
        try:
            taxon_ids_collection[correct_name] = t_id
        except (NameError):
            pass
            print(
                "> Will ignore organism '{}' for now. Try to handle it manually later."
                .format(correct_name))
        except (KeyboardInterrupt):
            print("\nQUIT: TaIGa was stopped by the user.\n")
        except:
            pass
            print(
                "> Unknown error occurred while trying to save the TaxID for organism '{}'."
                .format(correct_name))
    for org in list(taxon_ids_collection.keys()):
        org_info = retrive_taxonomy(user_email, taxon_ids_collection[org])
        tax_info.append({org: []})
        for taxon in org_info:
            taxons = {taxon["Rank"]: taxon["ScientificName"]}
            for item in tax_info:
                item[org].append(taxons)

else:  # If there's only one record, or only one organism, a loop isn't needed
    # Use Entrez.espell to correct the spelling of organism names
    try:
        correct_name = correct_spell(user_email, names)
        print("> Correcting organism name of '{}' to: '{}'".format(
            names, correct_name))
        original_names[names] = correct_name
    except (RuntimeError):
        pass
        print(
            "> Couldn't find the correct organism name for '{}'".format(names))
        missing_corrected.append(names)
    except (KeyboardInterrupt):
        print("\nQUIT: TaIGa was stopped by the user.\n")
    except:
        pass
        print(
            "> Unknown error occurred while trying to correct the spelling for organism '{}'"
            .format(names))
        missing_corrected.append(names)
    # Search Taxonomy for the TaxID of the input organism names
    try:
        t_id = search(user_email, correct_name)
        print("> TaxID of organism '{}': '{}'".format(names, t_id))
    except (IndexError):
        pass
        print("> Couldn't find a valid TaxID for the organism {}".format(
            correct_name))
        missing_taxid.append(correct_name)
    except (KeyboardInterrupt):
        print("\nQUIT: TaIGa was stopped by the user.\n")
    except:
        pass
        print(
            "> Unknown error occurred while trying to find a valid TaxID for organism '{}'."
            .format(correct_name))
        missing_taxid.append(correct_name)
    # Add the name and taxid to the dictionary
    try:
        taxon_ids_collection[correct_name] = t_id
    except (NameError):
        pass
        print(
            "> Will ignore organism '{}' for now. Try to handle it manually later."
            .format(correct_name))
    except (KeyboardInterrupt):
        print("\nQUIT: TaIGa was stopped by the user.\n")
    except:
        pass
        print(
            "> Unknown error occurred while trying to save the TaxID for organism '{}'."
            .format(correct_name))
    for org in list(taxon_ids_collection.keys()):
        org_info = retrive_taxonomy(user_email, taxon_ids_collection[org])
        tax_info.append({org: []})
        for taxon in org_info:
            taxons = {taxon["Rank"]: taxon["ScientificName"]}
            for item in tax_info:
                item[org].append(taxons)

# print(SeqIO.read(genome, "genbank"))
# for seq in records:
# print(seq.annotations["organism"])
# print(tax_names)

# Context manager used to retrieve all the needed information from the genome file. Is also used to generate the name list for the request to the "Taxonomy name/id Status Report Page"
# with open(out_folder + "/tax_info", "w") as f:
#     if key == "length":
#         for seq in records:
#             print(seq.annotations['organism'])
#             f.write(str(len(seq.seq)) + "\n")
#     elif key == "id" or key == "name":
#         for seq in records:
#             print(seq.annotations['organism'])
#             f.write(seq.id + "\n")
#     else:
#         for seq in records:
#             print(seq.annotations['organism'])
#             f.write(seq.annotations[key] + "\n")
#     print("Done!")
