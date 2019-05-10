################################################################################################################################
"""                                          TaIGa - TAxonomy Information GAtherer
This is a simple script that interacts with various utilities from the NCBI's Entrez api in order to retrieve relevant taxonomic information for a collection of organisms. As of now, TaIGa is able to handle multiple types of Genbank format genome files, as well as a text file format list of organism names, separated by lines. TaIGa recieves a file as input, an output folder path, a valid user e-mail and one optional argument to identify the type of input file. TaIGa then uses Entrez to retrieve the TaxID, Genome ID and all taxonomic information of all taxa up to the organism name provided. Then, it builds a DataFrame and outputs it to a .csv file, so the user can visualize it as a table. TaIGa's goal is to make easier for researchers to gather mass taxonomical metadata for their projects. Therefore, TaIGa is best used when you have a big list of organisms or a big collection of genomes in a file. TaIGa is also a very cute anime character from the japanese romance animation ToraDora. You should watch it. 
                        TaIGa was developed and is maintained by Maycon Douglas de Oliveira - 2019                          """

import os
import sys
import argparse
import pandas as pd
from Bio import Entrez
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def organize_tax_info(user_email, orgs_dict):
    """ Creates a list of dictionaries for every organism from the input file, each dictionary containing all relevant output information for that organism. The 'taxonomy' key refers to a list o dictionaries, each containing 'rank': 'name' key:value pairs. Returns the full list. """
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


def retrive_taxonomy(user_email, tax_id):
    """ Uses Entrez.efetch to fetch taxonomical data for the input taxon. Returns all taxonomical information for that taxon. """
    Entrez.email = user_email
    query = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
    tax_info = Entrez.read(query)

    return tax_info[0]["LineageEx"]


#Taking the user input from the command line. The path to the input file and the path to the output folder. Both are taken from stdin. Aditional options let TaIGa understand the format of the input.
taiga = argparse.ArgumentParser(
    description=
    "TaIGa retrieves metadata of an organism (or a collection of organisms) from a Genbank format genome file or a list of names."
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
    "--name",
    help=
    "Use this option if you want to give TaIGa a list of species names. When using this option, give TaIGa a list of species names in a text file, separated by line (linebreaks).",
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

if not os.path.exists(
        output_path):  # Checking if provided output path is valid
    sys.exit(
        "\nERROR: The provided output folder is not a valid path. It doesn't need to be created, but the path must be valid. Check it and run TaIGa again.."
    )

# Checking if only one optional argument was passed to TaIGa
if (args.name and args.same) or (args.name and args.single) or (args.single
                                                                and args.same):
    sys.exit(
        "\nERROR: Please run TaIGa with only one of the possible optional arguments."
    )

# First checking if input file is a simple list of names
if args.name:
    try:
        print("\n>> Parsing input file a simple list of species names.\n")
        with open(input_path, "r") as list_of_names:
            names = list_of_names.readlines(
            )  # Take the names of the organisms
            for name in range(len(names)):
                names[name] = names[name].replace(
                    "\n", "")  # Correct the formatting of each name
                print("{} ---> All OK".format(names[name]))
    except (KeyboardInterrupt):
        sys.exit("\nQUIT: TaIGa was stopped by the user.\n")
    except:
        sys.exit(
            "\nERROR: Couldn't parse name list. Check your file an try running TaIGa again.\n"
        )
# Else, checking the type of input Genbank file
elif args.single:  # Single record on Genbank file
    try:
        print("\n>> Parsing input as a Genbank file with a single record.\n")
        try:
            records = SeqIO.read(input_path, "genbank")
        except (KeyboardInterrupt):
            sys.exit("\nQUIT: TaIGa was stopped by the user.\n")
        except:
            sys.exit(
                "\nERROR: Couldn't parse input genome file. Check your file and try running TaIGa again.\n"
            )
        names = records.annotations[
            "organism"]  # Take the name of the organism on the input file
        if names:
            print("> '{}' ---> All OK".format(names))
        else:
            sys.exit(
                "\nERROR: Something went wrong while trying to parse the input file. Check the file and the program execution commands and try again.\n"
            )
    except (KeyboardInterrupt):
        sys.exit("\nQUIT: TaIGa was stopped by the user.\n")
    # Catch an error if the file contains more than one record
    except (ValueError):
        sys.exit(
            "\nERROR: The file contains more than one record. Try running TaIGa again without the '--single' option.\n"
        )
elif args.same:  # Multiple records from the same organism
    print(
        "\n>> Parsing input as a Genbank file with multiple records for the same organism.\n"
    )
    try:
        records = list(SeqIO.parse(input_path, "genbank"))
    except (KeyboardInterrupt):
        sys.exit("\nQUIT: TaIGa was stopped by the user.\n")
    except:
        sys.exit(
            "\nERROR: Couldn't parse input genome file. Check your file and try running TaIGa again.\n"
        )
    names = records[0].annotations[
        "organism"]  # Only take the first taxon name retrieved, as all are from the same organism
    if names:
        print("> '{}' ---> All OK".format(names))
    else:
        sys.exit(
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
        sys.exit("\nQUIT: TaIGa was stopped by the user.\n")
    except:
        sys.exit(
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
        sys.exit(
            "\nERROR: Something went wrong while trying to parse the input file. Check the file and the program execution commands and try again.\n"
        )

# Done collecting all organism names from the input file. Print a message indicating the next step.
print("\n>> Searching for taxonomic information...\n")

# Checking if there are multiple records on the file (for different organisms)
if type(names) == list:
    for name in names:
        # Use Entrez.espell to correct the spelling of organism names
        try:
            print("> Correcting organism name of '{}'".format(name))
            correct_name = correct_spell(user_email, name)
            if len(correct_name) == 0:
                print(
                    "\n\t>> Couldn't find the correct organism name for '{}'\n"
                    .format(name))
                missing_corrected.append(name)
        except (RuntimeError):
            pass
            print("\n\t>> Couldn't find the correct organism name for '{}'\n".
                  format(name))
            missing_corrected.append(name)
        except (KeyboardInterrupt):
            sys.exit("\nQUIT: TaIGa was stopped by the user.\n")
        except:
            pass
            print(
                "\n\t>> Unknown error occurred while trying to correct the spelling for organism '{}'.\n"
                .format(name))
            missing_corrected.append(name)
        # Search Taxonomy for the TaxID of the input organism names and Genome for the GenomeID
        try:
            print("> Searching TaxID of organism '{}'".format(name))
            t_id = search(user_email, "taxonomy", correct_name)
            g_id = search(user_email, "genome", correct_name)
            print(" >>>> TaxID for '{}' : '{}'\n".format(correct_name, t_id))
        except (IndexError):
            pass
            print(
                "\n\t>> Couldn't find a valid TaxID for the organism '{}'.\n".
                format(name))
            missing_taxid.append(name)
        except (KeyboardInterrupt):
            sys.exit("\nQUIT: TaIGa was stopped by the user.\n")
        except:
            pass
            if len(correct_name) == 0:
                print(
                    "\n\t>> Organism '{}' lacks a corrected name. Unable to search for TaxID.\n"
                    .format(name))
            else:
                print(
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
            print(
                "\n\t>> Will ignore organism '{}' for now. Try to handle it manually later.\n"
                .format(name))
        except (KeyboardInterrupt):
            sys.exit("\nQUIT: TaIGa was stopped by the user.\n")
        except:
            pass
            print(
                "\n\t>> Unknown error occurred while trying to save the TaxID for organism '{}'.\n"
                .format(name))
else:  # If there's only one record, or only one organism, a loop isn't needed
    # Use Entrez.espell to correct the spelling of organism names
    try:
        print("> Correcting organism name of '{}'".format(names))
        correct_name = correct_spell(user_email, names)
        if len(correct_name) == 0:
            print("\n\t>> Couldn't find the correct organism name for '{}'\n".
                  format(names))
            missing_corrected.append(names)
    except (RuntimeError):
        pass
        print(
            "\n\t>> Couldn't find the correct organism name for '{}'\n".format(
                names))
        missing_corrected.append(names)
    except (KeyboardInterrupt):
        sys.exit("\nQUIT: TaIGa was stopped by the user.\n")
    except:
        pass
        print(
            "\n\t>> Unknown error occurred while trying to correct the spelling for organism '{}'\n"
            .format(names))
        missing_corrected.append(names)
    # Search Taxonomy for the TaxID of the input organism names and Genome for the GenomeID
    try:
        print("> Searching TaxID of organism '{}'".format(names))
        t_id = search(user_email, "taxonomy", correct_name)
        g_id = search(user_email, "genome", correct_name)
        print(" >>>> TaxID for '{}' : '{}'\n".format(correct_name, t_id))
    except (IndexError):
        pass
        print(
            "\n\t>> Couldn't find a valid TaxID for the organism {}\n".format(
                names))
        missing_taxid.append(names)
    except (KeyboardInterrupt):
        sys.exit("\nQUIT: TaIGa was stopped by the user.\n")
    except:
        pass
        if len(correct_name) == 0:
            print(
                "\n\t>> Organism '{}' lacks a corrected name. Unable to search for TaxID.\n"
                .format(name))
        else:
            print(
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
        print(
            "\n\t>> Will ignore organism '{}' for now. Try to handle it manually later.\n"
            .format(names))
    except (KeyboardInterrupt):
        sys.exit("\nQUIT: TaIGa was stopped by the user.\n")
    except:
        pass
        print(
            "\n\t>> Unknown error occurred while trying to save the TaxID for organism '{}'.\n"
            .format(names))

try:
    # Check for the organisms with missing corrected name or taxid and remove their name form the names list
    for missed in missing_corrected:
        if missed in names:
            names.remove(missed)
    for missed in missing_taxid:
        if missed in names:
            names.remove(missed)

    # Done getting all preliminary information, print a message informing the next step
    print(
        ">> Gathering taxonomic information from NCBI Taxonomy. This might take a while.\n"
    )

    # Creating the list of dictionaries to store the taxonomic information for the organisms
    tax_info = organize_tax_info(user_email, taxon_ids_collection)

    # Done getting and organizing all taxonomic information, print a message informing that
    print(">> Done gathering taxonomic information\n")

    # Append to the tax_info list of dictionaries the corresponding tax_id and genome_id for each organism
    for i in list(taxon_ids_collection.keys()):
        for j in tax_info:
            if j["name"] == i:
                j["tax_id"] = taxon_ids_collection[i]
                j["genome_id"] = genome_ids_collection[i]

    # Done finishing the tax_info dictionary, print a message informing the next step
    print(">> Generating result table")

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

    print(
        "\n>> Done generating result table. Checking if output folder exists")

    # Check if the output directory exists. If not, create it.
    if not os.path.exists(output_path):
        print("\n>> Creating output folder")
        os.makedirs(output_path)

    print(
        "\n>> Creating output file. You'll find it inside the provided output folder, named 'TaIGa_result.csv'"
    )
    # Export the DataFrame to the resulting .csv file
    frame.to_csv(output_path + 'TaIGa_result.csv')

    # Checking if there are missing correct names or TaxIDs. If there are, generating log files for those.
    if missing_corrected or missing_taxid:
        print(
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
    sys.exit("\nQUIT: TaIGa was stopped by the user.\n")
except:
    sys.exit(
        "\nQUIT: Unknown error occured while generating output files. Try checking your inputs and running TaIGa again."
    )

print(
    "\n>> TaIGa was run successfully! You can check your results on the informed output folder. If there's any missing data, check the 'TaIGa_missing.txt' file on the same folder.\n"
)
