import argparse

taiga = argparse.ArgumentParser(
        description="""TaIGa retrieves metadata of an organism (or a collection of organisms)\
        from a list of names, taxon IDs or a Genbank format genome file.""")

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
taiga.add_argument("--tid", help="Use this option if you want to run TaIGa from a text file with a list of valid Taxonomy TaxIDs.", action="store_true")
taiga.add_argument(
    "-c",
    help=
    "Use this to enable TaIGa's name correcting function. Beware: sometimes, this function will alter organism names unecessarily and thus result in missing information returned (which could be misleading).",
    action="store_true")
taiga.add_argument(
    "-t",
    help=
    "Set the maximum ammount of retries for TaIGa's requests. By default, this number is 5.",
    nargs=1, default="5")
taiga.add_argument(
    "-v",
    help=
    "Turn off TaIGa's standard verbose mode. This way, all prints that would usually go to stdout will be logged to a file on TaIGa's current folder and she will print to your screen no more.",
    action="store_true")

args = taiga.parse_args()

# Define the path to the input file
input_path = args.input
# Define the path to the output file, adding an ending forward slash if needed
output_path = args.outdir + '/' if args.outdir[-1] != '/' else args.outdir
# Define the user email for Entrez.email, as it is a good practice of E-utils
user_email = args.email
# Define maximum number of retries
retries = int(args.t[0])
# True or false for verbose
verbose = args.v
# True or false for correction
correction = args.c
# True or false for taxon id list input
tid = args.tid
# True or false for multiple gbff file input
multi = args.multi
# True or false for same gbff file input
same = args.same
# True or false for single gbff file input
single = args.single
