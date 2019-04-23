# TaIGa - Taxonomy Information Gatherer
# Author: Maycon D. Oliveira

import argparse
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

    return tax_ids


# Taking the user input from the command line. The path to the input file and the path to the output folder. Both are taken from stdin. Aditional options let TaIGa understand the format of the input.
taiga = argparse.ArgumentParser(
    description=
    "TaIGa retrieves metadata of an organism (or a collection of organisms) from a Genbank format genome file."
)
taiga.add_argument(
    "input_file",
    help=
    "Full path to input file (or only the file name if it is on the same folder)"
)
taiga.add_argument(
    "output_folder",
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

# Checking the type of input Genbank file
if args.single:  # Single record on Genbank file
    try:
        records = SeqIO.read(args.input_file, "genbank")
        names = records.annotations[
            "organism"]  # Take the name of the organism on the input file
    except (ValueError
            ):  # Catch an error if the file contains more than one record
        print(
            "ERROR: The file contains more than one record. Try again without the '--single' option."
        )
elif args.same:  # Multiple records from the same organism
    records = list(SeqIO.parse(args.input_file, "genbank"))
    names = records[0].annotations[
        "organism"]  # Only take the first taxon name retrieved, as all are from the same organism
    print(names)
else:  # Multiple records from multiple organisms
    records = SeqIO.parse(args.iput_file, "genbank")
    names = [seq.annotations["organism"]
             for seq in records]  # List the names of all organisms

user_email = args.email  # Define the user email for Entrez.email, as it is a good practice of E-utils
correct_names = correct_spell(
    user_email,
    names)  # Use Entrez.espell to correct the spelling of organism names
t_ids = search(
    user_email,
    correct_names)  # Search Taxonomy for the TaxID of the input organism names

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
