################################################################################################################################
"""                                          TaIGa - TAxonomy Information GAtherer
This is a simple script that interacts with various utilities from the NCBI's Entrez api in order to retrieve relevant taxonomic
information for a collection of organisms. As of now, TaIGa is able to handle multiple types of Genbank format genome files, as
well as a text file format list of organism names, separated by lines. TaIGa recieves a file as input, an output folder path, a valid user e-mail and one optional argument to identify the type of input file. TaIGa then uses Entrez to retrieve the TaxID, 
Genome ID and all taxonomic information of all taxa up to the organism name provided. Then, it builds a DataFrame and outputs it
to a .csv file, so the user can visualize it as a table. TaIGa's goal is to make easier for researchers to gather mass 
taxonomical metadata for their projects. Therefore, TaIGa is best used when you have a big list of organisms or a big collection 
of genomes in a file. TaIGa is also a very cute anime character from the japanese romance animation ToraDora. You should watch 
it. 
                        TaIGa was developed and is maintained by Maycon Douglas de Oliveira - 2019                          """
################################################################################################################################

import os
import sys
import argparse
import pandas as pd
import logging as log
from time import sleep
from random import randint
from collections import OrderedDict
from Bio import Entrez
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def organize_tax_info(user_email, orgs_dict, retries):
    """ Creates a list of dictionaries for every organism from the input file, each dictionary containing all relevant output information for that organism. The 'taxonomy' key refers to a list o dictionaries, each containing 'rank': 'name' key:value pairs. Returns the full list. """
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
    """ Uses Entrez.espell utility to correct the spelling of organism names. Returns the corrected name. """
    Entrez.email = user_email
    query = Entrez.espell(db="taxonomy", term=names)
    corrected_names = Entrez.read(query)

    return corrected_names["CorrectedQuery"]


def search(user_email, db, names):
    """ Uses Entrez.esearch to either search for oranisms TaxIDs or Genome IDs. Returns either one of those based on the 'db' argument. """
    Entrez.email = user_email
    query = Entrez.esearch(db=db, term=names, retmode="xml")
    parsed = Entrez.read(query)
    if db == "taxonomy":
        return parsed["IdList"][0]
    elif db == "genome":
        return parsed["IdList"][-1]


def retrive_taxonomy(user_email, tax_id, retries):
    """ Uses Entrez.efetch to fetch taxonomical data for the input taxon. Returns all taxonomical information for that taxon. """
    for i in range(retries):
        try:
            Entrez.email = user_email
            query = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
            tax_info = Entrez.read(query)

            return tax_info[0]["LineageEx"]
        except (IndexError):
            print("\nERROR: fetching for TaxID '{}' | retry: {} ".format(
                tax_id, i + 1))
            sleep(randint(0, 5))
            continue

    raise IndexError


#Taking the user input from the command line. The path to the input file and the path to the output folder. Both are taken from stdin. Aditional options let TaIGa understand the format of the input.
taiga = argparse.ArgumentParser(
    description=
    "TaIGa retrieves metadata of an organism (or a collection of organisms) from a text file format list of names or a Genbank format genome file."
)
taiga.add_argument(
    "input",
    help=
    "Full path to input file (or only the file name if it is on the same folder)"
)
taiga.add_argument(
    "outdir",
    help=
    "Full path to output folder (folder doensn't need to be pre-existent. TaIGa will create output files automatically)"
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
    "--multi",
    help=
    "Use this option if you want to run TaIGa from a Genbank format genome file with multiple records from multiple different organisms. TaIGa does check for duplicate names and automatically removes them.",
    action="store_true")
taiga.add_argument(
    "-c",
    help=
    "Use this to disable TaIGa's name correcting function. Sometimes, this function will alter organism names unecessarily and thus result in missing information returned (which could be misleading).",
    action="store_true")
taiga.add_argument(
    "-t",
    help=
    "Set the maximum ammount of retries for TaIGa's requests. By default, this number is 5.",
    nargs=1)
taiga.add_argument(
    "-v",
    help=
    "Turn off TaIGa's standard verbose mode. This way, all prints that would usually go to stdout will be logged to a file on TaIGa's current folder and she will print to your screen no more.",
    action="store_true")

# Parsing the program arguments
args = taiga.parse_args()
# Define the path to the input file
input_path = args.input
# Define the path to the output file
output_path = args.outdir + '/' if args.outdir[-1] != '/' else args.outdir
# Define the user email for Entrez.email, as it is a good practice of E-utils
user_email = args.email
# Define maximum number of retries
retries = int(args.t[0]) if args.t else 5
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

if args.v:
    log.basicConfig(
        filename='TaIGa_run.log', format="%(message)s", level=log.DEBUG)
else:
    log.basicConfig(format="%(message)s", level=log.DEBUG)

log.info("""*********************************************
*                                           *
*   TaIGa - Taxonomy Information Gatherer   *
*                                           *
*********************************************""")

# Checking if only one optional argument was passed to TaIGa
if (args.multi and args.same) or (args.multi
                                  and args.single) or (args.single
                                                       and args.same):
    log.error(
        "\nERROR: Please run TaIGa with only one of the possible input type optional arguments."
    )
    sys.exit()

# First checking if input file is a genome file with multiple records from multiple organisms
if args.multi:
    name_counter = 0
    log.info(
        "\n>> Parsing input as a Genbank file with multiple records for multiple organisms.\n"
    )
    try:
        records = SeqIO.parse(input_path, "genbank")
    except (KeyboardInterrupt):
        log.warning("\nQUIT: TaIGa was stopped by the user.\n")
        sys.exit()
    except:
        log.error(
            "\nERROR: Couldn't parse input genome file. Check your file and try running TaIGa again.\n"
        )
        sys.exit()
    names = [seq.annotations["organism"]
             for seq in records]  # List the names of all organisms
    if names:
        log.info(
            "\n>> All OK with parsing the input file. Checking the records...\n"
        )
        for name in names:
            name_counter += 1
            if name:
                log.info("> '{}' ---> All OK".format(name))
            else:
                log.info(
                    ">> Something went wrong while trying to retrieve the organism name of record number '{}'"
                    .format(name_counter))
                continue
    else:
        log.error(
            "\nERROR: Something went wrong while trying to parse the input file. Check the file and the program execution commands and try again.\n"
        )
        sys.exit()

# Else, checking the type of input Genbank file
elif args.single:  # Single record on Genbank file
    try:
        log.info(
            "\n>> Parsing input as a Genbank file with a single record.\n")
        try:
            records = SeqIO.read(input_path, "genbank")
        except (KeyboardInterrupt):
            log.warning("\nQUIT: TaIGa was stopped by the user.\n")
            sys.exit()
        except:
            log.error(
                "\nERROR: Couldn't parse input genome file. Check your file and try running TaIGa again.\n"
            )
            sys.exit()
        names = records.annotations[
            "organism"]  # Take the name of the organism on the input file
        if names:
            log.info("> '{}' ---> All OK".format(names))
        else:
            log.error(
                "\nERROR: Something went wrong while trying to parse the input file. Check the file and the program execution commands and try again.\n"
            )
            sys.exit()
    except (KeyboardInterrupt):
        log.warning("\nQUIT: TaIGa was stopped by the user.\n")
        sys.exit()
    # Catch an error if the file contains more than one record
    except (ValueError):
        log.error(
            "\nERROR: The file contains more than one record. Try running TaIGa again without the '--single' option.\n"
        )
        sys.exit()
elif args.same:  # Multiple records from the same organism
    log.info(
        "\n>> Parsing input as a Genbank file with multiple records for the same organism.\n"
    )
    try:
        records = list(SeqIO.parse(input_path, "genbank"))
    except (KeyboardInterrupt):
        log.warning("\nQUIT: TaIGa was stopped by the user.\n")
        sys.exit()
    except:
        log.error(
            "\nERROR: Couldn't parse input genome file. Check your file and try running TaIGa again.\n"
        )
        sys.exit()
    names = records[0].annotations[
        "organism"]  # Only take the first taxon name retrieved, as all are from the same organism
    if names:
        log.info("> '{}' ---> All OK".format(names))
    else:
        log.error(
            "\nERROR: Something went wrong while trying to parse the input file. Check the file and the program execution commands and try again.\n"
        )
        sys.exit()
else:  # Simple text file with a list of organism names
    try:
        log.info("\n>> Parsing input file a simple list of species names.\n")
        with open(input_path, "r") as list_of_names:
            names = list_of_names.readlines(
            )  # Take the names of the organisms
            for name in range(len(names)):
                names[name] = names[name].replace(
                    "\n", "")  # Correct the formatting of each name
                log.info("{} ---> All OK".format(names[name]))
    except (KeyboardInterrupt):
        log.warning("\nQUIT: TaIGa was stopped by the user.\n")
        sys.exit()
    except:
        log.error(
            "\nERROR: Couldn't parse name list. Check your file an try running TaIGa again.\n"
        )
        sys.exit()

# Done getting all names. Checking for duplicates.
log.info(
    "\n>> Ignore if any duplicate name was printed. Checking and removing duplicate names."
)

names = list(dict.fromkeys(
    names))  # Using Python's collections module to uniquefy all names in list

# Done collecting all organism names from the input file. Print a message indicating the next step.
log.info("\n>> Searching for taxonomic information...\n")

# Checking if there are multiple records on the file (for different organisms)
if type(names) == list:
    for name in names:
        # Use Entrez.espell to correct the spelling of organism names
        if not args.c:
            try:
                log.info("> Correcting organism name of '{}'".format(name))
                correct_name = correct_spell(user_email, name)
                if len(correct_name) == 0:
                    log.warning(
                        "\n\t>> Couldn't find the correct organism name for '{}'\n"
                        .format(name))
                    missing_corrected.append(name)
            except (RuntimeError):
                pass
                log.warning(
                    "\n\t>> Couldn't find the correct organism name for '{}'\n"
                    .format(name))
                missing_corrected.append(name)
            except (KeyboardInterrupt):
                log.warning("\nQUIT: TaIGa was stopped by the user.\n")
                sys.exit()
            except:
                pass
                log.info(
                    "\n\t>> Unknown error occurred while trying to correct the spelling for organism '{}'.\n"
                    .format(name))
                missing_corrected.append(name)
        # Search Taxonomy for the TaxID of the input organism names and Genome for the GenomeID
        try:
            if args.c:
                correct_name = name
            log.info("> Searching TaxID of organism '{}'".format(name))
            t_id = search(user_email, "taxonomy", correct_name)
            try:
                g_id = search(user_email, "genome", correct_name)
            except (IndexError):
                pass
                g_id = '-'
            log.info(" >>>> TaxID for '{}' : '{}'\n".format(
                correct_name, t_id))
        except (IndexError):
            pass
            log.warning(
                "\n\t>> Couldn't find a valid TaxID for the organism '{}'.\n".
                format(name))
            missing_taxid.append(name)
        except (KeyboardInterrupt):
            log.warning("\nQUIT: TaIGa was stopped by the user.\n")
            sys.exit()
        except:
            pass
            if len(correct_name) == 0:
                log.warning(
                    "\n\t>> Organism '{}' lacks a corrected name. Unable to search for TaxID.\n"
                    .format(name))
            else:
                log.warning(
                    "\n\t>> Unknown error occurred while trying to find a valid TaxID for organism '{}'.\n"
                    .format(name))
            missing_taxid.append(name)
        # Add the name, taxid and genomeid to the dictionary
        try:
            if name not in missing_corrected and name not in missing_taxid:
                taxon_ids_collection[name] = t_id
                genome_ids_collection[name] = g_id
        except (NameError):
            pass
            log.warning(
                "\n\t>> Will ignore organism '{}' for now. Try to handle it manually later.\n"
                .format(name))
        except (KeyboardInterrupt):
            log.warning("\nQUIT: TaIGa was stopped by the user.\n")
            sys.exit()
        except:
            pass
            log.warning(
                "\n\t>> Unknown error occurred while trying to save the TaxID for organism '{}'.\n"
                .format(name))
else:  # If there's only one record, or only one organism, a loop isn't needed
    # Use Entrez.espell to correct the spelling of organism names
    if not args.c:
        try:
            log.info("> Correcting organism name of '{}'".format(names))
            correct_name = correct_spell(user_email, names)
            if len(correct_name) == 0:
                log.warning(
                    "\n\t>> Couldn't find the correct organism name for '{}'\n"
                    .format(names))
                missing_corrected.append(names)
        except (RuntimeError):
            pass
            log.warning(
                "\n\t>> Couldn't find the correct organism name for '{}'\n".
                format(names))
            missing_corrected.append(names)
        except (KeyboardInterrupt):
            log.warning("\nQUIT: TaIGa was stopped by the user.\n")
            sys.exit()
        except:
            pass
            log.warning(
                "\n\t>> Unknown error occurred while trying to correct the spelling for organism '{}'\n"
                .format(names))
            missing_corrected.append(names)
    # Search Taxonomy for the TaxID of the input organism names and Genome for the GenomeID
    try:
        if args.c:
            correct_name = names
        log.info("> Searching TaxID of organism '{}'".format(names))
        t_id = search(user_email, "taxonomy", correct_name)
        try:
            g_id = search(user_email, "genome", correct_name)
        except (IndexError):
            g_id = '-'
        log.info(" >>>> TaxID for '{}' : '{}'\n".format(correct_name, t_id))
    except (IndexError):
        pass
        log.warning(
            "\n\t>> Couldn't find a valid TaxID for the organism {}\n".format(
                names))
        missing_taxid.append(names)
    except (KeyboardInterrupt):
        log.warning("\nQUIT: TaIGa was stopped by the user.\n")
        sys.exit()
    except:
        pass
        if len(correct_name) == 0:
            log.warning(
                "\n\t>> Organism '{}' lacks a corrected name. Unable to search for TaxID.\n"
                .format(name))
        else:
            log.warning(
                "\n\t>> Unknown error occurred while trying to find a valid TaxID for organism '{}'.\n"
                .format(names))
        missing_taxid.append(names)
    # Add the name, taxid and genomeid to the dictionary
    try:
        if names not in missing_corrected and names not in missing_taxid:
            taxon_ids_collection[names] = t_id
            genome_ids_collection[names] = g_id
    except (NameError):
        pass
        log.warning(
            "\n\t>> Will ignore organism '{}' for now. Try to handle it manually later.\n"
            .format(names))
    except (KeyboardInterrupt):
        log.error("\nQUIT: TaIGa was stopped by the user.\n")
        sys.exit()
    except:
        pass
        log.warning(
            "\n\t>> Unknown error occurred while trying to save the TaxID for organism '{}'.\n"
            .format(names))

try:
    if not args.c:
        # Check for the organisms with missing corrected name or taxid and remove their name form the names list
        for missed in missing_corrected:
            if missed in names:
                names.remove(missed)
    for missed in missing_taxid:
        if missed in names:
            names.remove(missed)

    # Done getting all preliminary information, print a message informing the next step
    log.info(
        ">> Gathering taxonomic information from NCBI Taxonomy. This might take a while.\n"
    )

    # Creating the list of dictionaries to store the taxonomic information for the organisms.
    tax_info = organize_tax_info(user_email, taxon_ids_collection, retries)

    # Done getting and organizing all taxonomic information, print a message informing that
    log.info(">> Done gathering taxonomic information\n")

    # Append to the tax_info list of dictionaries the corresponding tax_id and genome_id for each organism
    for i in list(taxon_ids_collection.keys()):
        for j in tax_info:
            if j["name"] == i:
                j["tax_id"] = taxon_ids_collection[i]
                j["genome_id"] = genome_ids_collection[i]

    # Done finishing the tax_info dictionary, print a message informing the next step
    log.info(">> Generating result table")

    # Creating a list to store all possible fields of information on the output file
    ranks = [
        'tax_id', 'genome_id', 'no rank', 'superkingdom', 'kingdom',
        'subkingdom', 'phylum', 'subphylum', 'superclass', 'class', 'subclass',
        'infraclass', 'cohort', 'subcohort', 'superorder', 'order', 'suborder',
        'infraorder', 'parvorder', 'superfamily', 'family', 'subfamily',
        'tribe', 'subtribe', 'genus', 'subgenus', 'section', 'subsection',
        'series', 'species-group', 'species', 'subspecies', 'forma'
    ]

    # Creating a set of unique ranks from all organisms with retrieved taxonomic information
    to_keep = {
        rank
        for org in tax_info for dic in org["taxonomy"] for rank in dic.keys()
    }

    # Adding the fields tax_id and genome_id to the set to avoid them from being removed later
    to_keep.add("tax_id")
    to_keep.add("genome_id")

    # Creating a set of rank fields that need to be removed from the final ranks
    tmp_remove = set([i for i in ranks if i not in to_keep])

    # Removing the unecessary ranks to avoid cluttering the output table
    for i in tmp_remove:
        ranks.remove(i)

    # Creating a list of dictionaries to hold the information for each organism. This makes it easier to create DataFrames with Pandas
    final_info = []

    # Populating the final_info list with the information for each organism, except the names. Using names as a key isn't needed, as lists are oredered and TaIGa keeps the order.
    for orgn in tax_info:
        tmp_orgn = {
        }  # Creating a temporary dict to hold info for each organism

        # Initializing 'no rank' key with the value as an empty array. 'no rank' will hold multiple values per cell
        tmp_orgn['no rank'] = []
        tmp_orgn['tax_id'] = orgn[
            'tax_id']  # Adding the tax_id for the organism
        tmp_orgn['genome_id'] = orgn[
            'genome_id']  # Adding the genome_id for the organism

        # Looping through each dictionary inside the 'taxonomy' key for each organism on the original tax_info dictionary
        for tax_rank in orgn["taxonomy"]:
            # Checking if the key for this particular dictionary is 'no rank'. If it is, append that value to the final 'no rank' key
            if list(tax_rank.keys())[0] == 'no rank':
                tmp_orgn['no rank'].append(tax_rank['no rank'])
            # Else, create a key for that particular rank and add its value on the new dictionary
            else:
                tmp_orgn[list(tax_rank.keys())[0]] = tax_rank[list(
                    tax_rank.keys())[0]]
        tmp_orgn['no rank'] = ", ".join(
            tmp_orgn['no rank']
        )  # Transform the 'no rank' list in a string separated by comma
        final_info.append(
            tmp_orgn)  # Append that organism dictionary to the final list

    # Create a DataFrame with the results. The values come from final_info, the names of the organisms are the indexes for the rows and the taxonomic ranks in ranks are the labels for the columns
    frame = pd.DataFrame(final_info, index=names, columns=ranks)
    # Transform the missing data to a more visual indicator for the user
    frame.fillna('-', inplace=True)

    log.info(
        "\n>> Done generating result table. Checking if output folder exists")

    # Check if the output directory exists. If not, create it.
    if not os.path.exists(output_path):
        log.info("\n>> Creating output folder")
        try:  # Checking if the path for the output folder is valid
            os.makedirs(output_path)
        except:
            log.error(
                "\nERROR: Path to output folder may not be valid. Try again.")
            sys.exit()

    log.info(
        "\n>> Creating output file. You'll find it inside the provided output folder, named 'TaIGa_result.csv'"
    )
    # Export the DataFrame to the resulting .csv file
    frame.to_csv(output_path + 'TaIGa_result.csv')

    # Checking if there are missing correct names or TaxIDs. If there are, generating log files for those.
    if missing_corrected or missing_taxid:
        log.info(
            "\n>> Creating a file for the organisms with missing information. You'll find it inside the provided output folder, named 'TaIGa_missing.txt'"
        )
        with open(output_path + 'TaIGa_missing.txt', 'w') as missing_file:
            missing_file.write("Missing corrected names: \n")
            if (missing_corrected):
                for name in missing_corrected:
                    missing_file.write("\t\t\t{}\n".format(name))
            missing_file.write("Missing TaxID: \n")
            if (missing_taxid):
                for taxid in missing_taxid:
                    missing_file.write("\t\t\t{}\n".format(taxid))
except (KeyboardInterrupt):
    log.warning("\nQUIT: TaIGa was stopped by the user.\n")
    sys.exit()
except (IndexError):
    log.error("\nQUIT: Too many broken responses.")
    sys.exit()
except:
    log.error("\nQUIT: Unknown error occured while generating output files.")
    sys.exit()

log.info(
    "\n>> TaIGa was run successfully! You can check your results on the informed output folder. If there's any missing data, check the 'TaIGa_missing.txt' file on the same folder.\n"
)
