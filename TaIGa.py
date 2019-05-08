# TaIGa - Taxonomy Information Gatherer
# Author: Maycon D. Oliveira

import os
import sys
import argparse
import pandas as pd
from Bio import Entrez
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def organize_tax_info(user_email, orgs_dict):
    taxonomic_info = []
    for org in list(orgs_dict.keys()):
        org_info = retrive_taxonomy(user_email, orgs_dict[org])
        taxonomic_info.append({"name": org, "taxonomy": []})
        for taxon in org_info:
            taxons = {taxon["Rank"]: taxon["ScientificName"]}
            for item in taxonomic_info:
                if item["name"] == org:
                    item["taxonomy"].append(taxons)

    return taxonomic_info


def correct_spell(user_email, names):
    Entrez.email = user_email
    query = Entrez.espell(db="taxonomy", term=names)
    corrected_names = Entrez.read(query)

    return corrected_names["CorrectedQuery"]


def search(user_email, db, names):
    Entrez.email = user_email
    query = Entrez.esearch(db=db, term=names, retmode="xml")
    parsed = Entrez.read(query)
    if db == "taxonomy":
        return parsed["IdList"][0]
    elif db == "genome":
        return parsed["IdList"][-1]


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
taiga.add_argument(
    "--name",
    help=
    "Use this option if you want to give TaIGa a list of species names. This way, TaIGa will retrieve all relevant taxonomic information for the list of names, but nothing more (no genome length, Genbank GenomeID). When using this option, give TaIGa a list of species names in a text file, separated by line (linebreaks).",
    action="store_true")

# Parsing the program arguments
args = taiga.parse_args()
# Define the path to the input file
input_path = args.input
# Define the path to the output file
output_path = args.outdir + '/' if args.outdir[-1] != '/' else args.outdir
# Define the user email for Entrez.email, as it is a good practice of E-utils
user_email = args.email
# A list to hold all organism names
names = []
# A dictionary to hold {'organism name': taxid} key:value pairs
taxon_ids_collection = {}
# A dictionary to hold {'organism name': genomeid} key:value pairs
genome_ids_collection = {}
# A list to hold all organisms for which a correct name wasn't found
missing_corrected = []
# A list to hold all organisms for which the TaxID wasn't found
missing_taxid = []
# A list of dictionaries to hold the taxonomic information for each organism. The keys are the corrected organism names, and the values are lists of dictionaries, each containing the {taxonomic level: taxon} key:value pair
tax_info = []

print("""*********************************************
*                                           *
*   TaIGa - Taxonomy Information Gatherer   *
*                                           *
*********************************************""")

# First checking if input file is a simple list of names
if args.name:
    try:
        print("\n>> Parsing input file a simple list of species names.\n")
        with open(input_path, "r") as list_of_names:
            names = list_of_names.readlines()
            for name in range(len(names)):
                names[name] = names[name].replace("\n", "")
                print("{} ---> All OK".format(names[name]))
    except (KeyboardInterrupt):
        print("\nQUIT: TaIGa was stopped by the user.\n")
        raise
    except:
        print(
            "\nERROR: Couldn't parse name list. Check your file an try running TaIGa again.\n"
        )
# Else, checking the type of input Genbank file
elif args.single:  # Single record on Genbank file
    try:
        print("\n>> Parsing input as a Genbank file with a single record.\n")
        try:
            records = SeqIO.read(input_path, "genbank")
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
        "\n>> Parsing input as a Genbank file with multiple records for the same organism.\n"
    )
    try:
        records = list(SeqIO.parse(input_path, "genbank"))
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
        "\n>> Parsing input as a Genbank file with multiple records for multiple organisms.\n"
    )
    try:
        records = SeqIO.parse(input_path, "genbank")
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
            "\n>> All OK with parsing the input file. Checking the records...\n"
        )
        for name in names:
            name_counter += 1
            if name:
                print("> '{}' ---> All OK".format(name))
            else:
                print(
                    ">> Something went wrong while trying to retrieve the organism name of record number '{}'"
                    .format(name_counter))
                continue
    else:
        print(
            "\nERROR: Something went wrong while trying to parse the input file. Check the file and the program execution commands and try again.\n"
        )

print("\n>> Searching for taxonomic information...\n")

# Checking if there are multiple records on the file (for different organisms)
if type(names) == list:
    for name in names:
        # Use Entrez.espell to correct the spelling of organism names
        try:
            print("> Correcting organism name of '{}'".format(name))
            correct_name = correct_spell(user_email, name)
        except (RuntimeError):
            pass
            print(">> Couldn't find the correct organism name for '{}'".format(
                name))
            missing_corrected.append(name)
        except (KeyboardInterrupt):
            print("\nQUIT: TaIGa was stopped by the user.\n")
            break
        except:
            pass
            print(
                ">> Unknown error occurred while trying to correct the spelling for organism '{}'."
                .format(name))
            missing_corrected.append(name)
        # Search Taxonomy for the TaxID of the input organism names
        try:
            print("> Searching TaxID of organism '{}'".format(name))
            t_id = search(user_email, "taxonomy", correct_name)
            g_id = search(user_email, "genome", correct_name)
            print(" >>>> TaxID for '{}' : '{}'\n".format(correct_name, t_id))
        except (IndexError):
            pass
            print(
                ">> Couldn't find a valid TaxID for the organism '{}'.".format(
                    correct_name))
            missing_taxid.append(correct_name)
        except (KeyboardInterrupt):
            print("\nQUIT: TaIGa was stopped by the user.\n")
            break
        except:
            pass
            print(
                ">> Unknown error occurred while trying to find a valid TaxID for organism '{}'."
                .format(correct_name))
            missing_taxid.append(correct_name)
        # Add the name and taxid to the dictionary
        try:
            taxon_ids_collection[name] = t_id
            genome_ids_collection[name] = g_id
        except (NameError):
            pass
            print(
                ">> Will ignore organism '{}' for now. Try to handle it manually later."
                .format(correct_name))
        except (KeyboardInterrupt):
            print("\nQUIT: TaIGa was stopped by the user.\n")
            break
        except:
            pass
            print(
                ">> Unknown error occurred while trying to save the TaxID for organism '{}'."
                .format(correct_name))
else:  # If there's only one record, or only one organism, a loop isn't needed
    # Use Entrez.espell to correct the spelling of organism names
    try:
        print("> Correcting organism name of '{}'".format(names))
        correct_name = correct_spell(user_email, names)
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
        print("> Searching TaxID of organism '{}'".format(names))
        t_id = search(user_email, "taxonomy", correct_name)
        g_id = search(user_email, "genome", correct_name)
        print(" >>>> TaxID for '{}' : '{}'\n".format(correct_name, t_id))
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
        taxon_ids_collection[names] = t_id
        genome_ids_collection[names] = g_id
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

print(
    ">> Gathering taxonomic information from NCBI Taxonomy. This might take a while.\n"
)

tax_info = organize_tax_info(user_email, taxon_ids_collection)

# for i in tax_info:
#     if "no rank" in i["taxonomy"]:
#         print(i["name"], i["taxonomy"]["no rank"])

print(">> Done gathering taxonomic information\n")

for i in list(taxon_ids_collection.keys()):
    for j in tax_info:
        if j["name"] == i:
            j["tax_id"] = taxon_ids_collection[i]
            j["genome_id"] = genome_ids_collection[i]

print(">> Generating output file\n")

ranks = [
    'tax_id', 'genome_id', 'no rank', 'superkingdom', 'kingdom', 'subkingdom',
    'phylum', 'subphylum', 'superclass', 'class', 'subclass', 'infraclass',
    'cohort', 'subcohort', 'superorder', 'order', 'suborder', 'infraorder',
    'parvorder', 'superfamily', 'family', 'subfamily', 'tribe', 'subtribe',
    'genus', 'subgenus', 'section', 'subsection', 'series', 'species-group',
    'species', 'subspecies', 'forma'
]

to_keep = {
    rank
    for org in tax_info for dic in org["taxonomy"] for rank in dic.keys()
}

to_keep.add("tax_id")
to_keep.add("genome_id")

tmp_remove = set([i for i in ranks if i not in to_keep])
for i in tmp_remove:
    ranks.remove(i)

final_info = []

for orgn in tax_info:
    tmp_orgn = {}
    tmp_orgn['no rank'] = []
    tmp_orgn['tax_id'] = orgn['tax_id']
    tmp_orgn['genome_id'] = orgn['genome_id']
    for tax_rank in orgn["taxonomy"]:
        if list(tax_rank.keys())[0] == 'no rank':
            tmp_orgn['no rank'].append(tax_rank['no rank'])
        else:
            tmp_orgn[list(tax_rank.keys())[0]] = tax_rank[list(
                tax_rank.keys())[0]]
    tmp_orgn['no rank'] = ", ".join(tmp_orgn['no rank'])
    final_info.append(tmp_orgn)

frame = pd.DataFrame(final_info, index=names, columns=ranks)
frame.fillna('-', inplace=True)

frame.to_csv(output_path + 'tax_report.csv')
# print(">> Done generating table.\n")

# print(">> Showing a preview of the table.\n")

# print(
#     ">> TaIGa was run successfully! You can check your results on the informed output folder. If there's any missing data, check the 'TaIGa_missing.txt' file on the same folder.\n"
# )