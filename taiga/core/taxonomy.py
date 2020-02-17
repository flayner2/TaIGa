import sys
import logging as log
import argparse
from taiga.common import parsers, helpers, retrievers, data_handlers


def run_taiga():
    """ Wrapper for all of TaIGa's main functionalities

    Gets the user input, parses the input files and run all the needed
    logic and functions, producing the output files at the specified
    folder.
    """

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


    # A list to hold all organism names
    input_names = []
    # A dictionary to hold {'organism name': taxid} key:value pairs
    taxon_ids_collection = {}
    # A dictionary to hold {'organism name': genomeid} key:value pairs
    genome_ids_collection = {}
    # A list to hold all organisms for which a correct name wasn't found
    missing_corrected = []
    # A list to hold all organisms for which the TaxID wasn't found
    missing_taxid = []
    # A list to hold all organisms for which the name wasn't found
    missing_name = []
    # A list of dictionaries to hold the taxonomic information for each organism. The keys are the corrected organism names, and the values are lists of dictionaries, each containing the {taxonomic level: taxon} key:value pair
    tax_info = []

    helpers.config_log(verbose)

    log.info("""
    *********************************************
    *                                           *
    *   TaIGa - Taxonomy Information Gatherer   *
    *                                           *
    *********************************************""")

    # Checking if only one input format argument was passed to TaIGa
    helpers.sanitize_args(args)

    # Checking if TaIGa is being run on TaxID mode with the '-c' argument.
    # This is needed because, when run with '--tid', TaIGa never actually calls 'correct_spell'.
    # The retrieved name is assumed to be correct.
    if tid and correction:
        log.error("\nERROR: Please, when running TaIGa with the '--tid' option, don't use '-c', as TaIGa skips the name correction.")
        sys.exit()

    # First checking if input file is a text file with multiple TaxIDs
    if tid:
        input_taxid = parsers.parse_taxid(input_path)
    # Then checking if input file is a genome file with multiple records from multiple organisms
    elif multi:
        input_names = parsers.parse_gb_multi(input_path)
    # Else, checking the type of input Genbank file
    elif single:  # Single record on Genbank file
        input_names = parsers.parse_gb_single(input_path)
    elif same:  # Multiple records from the same organism
        input_names = parsers.parse_gb_same(input_path)
    else:  # Simple text file with a list of organism names
        input_names = parsers.parse_names(input_path)

    log.info("\n>> Ignore if any duplicate name was printed, duplicates are removed.")

    log.info("\n>> Searching for taxonomic information...\n")

    # Checking if the input is a list of TaxIDs instead of any type of name input
    if tid:
        taxon_ids_collection, \
            genome_ids_collection, \
            missing_name = retrievers.retrieve_from_taxid(
                input_taxid,
                user_email,
                retries)
    # Checking if there are multiple records on the file (for different organisms)
    elif type(input_names) == list:
        missing_corrected, \
            missing_taxid, \
            taxon_ids_collection, \
            genome_ids_collection = retrievers.retrieve_from_multiple_names(
                input_names,
                user_email,
                correction,
                retries)
        print(taxon_ids_collection)
    else:  # If there's only one record, or only one organism, a loop isn't needed
        missing_corrected, \
            missing_taxid, \
            taxon_ids_collection, \
            genome_ids_collection = retrievers.retrieve_from_single_name(
                input_names,
                user_email,
                correction,
                retries)

    try:
        helpers.sanitize_name_list(correction, tid,
                                   missing_corrected, missing_taxid,
                                   input_names)

        tax_info = retrievers.get_tax_info(
            user_email,
            taxon_ids_collection,
            retries)

        frame = data_handlers.create_df(
            taxon_ids_collection,
            genome_ids_collection,
            tax_info,
            tid,
            input_names)

        data_handlers.create_output(
            output_path,
            frame,
            missing_corrected,
            missing_taxid,
            missing_name)
    except (KeyboardInterrupt):
        log.warning("\nQUIT: TaIGa was stopped by the user.\n")
        sys.exit()
    except (IndexError):
        log.error("\nQUIT: Too many broken responses.")
        sys.exit()
    except (Exception) as e:
        #print(e)
        log.error("\nQUIT: Unknown error occured while generating output files.")
        sys.exit()

    log.info(
        "\n>> TaIGa was run successfully! You can check your results on the informed output folder. If there's any missing data, check the 'TaIGa_missing.txt' file on the same folder.\n"
    )
