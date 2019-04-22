# TaIGa - Taxonomy Information Gatherer
# Author: Maycon D. Oliveira

import requests
import sys
import argparse
from Bio import Entrez
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# Taking the user input from the command line. The path to the input file and the path to the output folder. Both are taken from stdin
if len(sys.argv) != 4:
    sys.exit(
        "Run TaIGa as: python3 TaIGa.py <path to input file> <path to output folder> <your email>"
    )
else:
    genome = sys.argv[1]
    out_folder = sys.argv[2]
    user_email = sys.argv[3]

# Storing all the records of the genome file on the 'records' variable
records = SeqIO.parse(genome, "genbank")
print(SeqIO.read(genome, "genbank"))
for seq in records:
    print(seq.annotations["organism"])
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
