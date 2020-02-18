import sys
import logging as log
import argparse
from taiga.common import parsers, helpers, retrievers, data_handlers


def run_taiga():
    """ Wrapper for all of TaIGa's main functionalities

    Parameters:
    None

    Returns:
    None

    """

    taiga = argparse.ArgumentParser(
        description="TaIGa retrieves metadata of an organism (or a collection of organisms) from a\
        list of names, taxon IDs or a Genbank format genome file.")

    taiga.add_argument(
        "infile",
        help="Full path to input file (or only the file name if it is on the same folder)")
    taiga.add_argument(
        "outdir",
        help="Full path to output folder (folder doensn't need to be pre-existent. TaIGa will\
        create the output files automatically)")
    taiga.add_argument(
        "email",
        help="Any valid email of yours. It is a common practice when using E-utils, which TaIGa\
        does use")
    taiga.add_argument(
        "--gb-mode",
        help="Change the expected Genbank file format input. Only use if your input file is in\
        Genbank format. Refer to the documentation for an explanation on the options.",
        type=int,
        nargs=1,
        default=0,
        choices=[0, 1, 2, 3],
        dest='mode'
    )
    taiga.add_argument(
        "--tid", help="Use this option if you want to run TaIGa from a text file with a list of\
        valid Taxonomy TaxIDs.", action="store_true")
    taiga.add_argument(
        "-c",
        help="Use this to enable TaIGa's name correcting function. Beware: sometimes, this function\
        will alter organism names unecessarily and thus result in missing information returned\
        (which could be misleading).",
        action="store_true")
    taiga.add_argument(
        "-t",
        help="Set the maximum ammount of retries for TaIGa's requests. By default, this number is\
        5.",
        nargs=1, default="5", type=int)
    taiga.add_argument(
        "-v",
        help="Turn off TaIGa's standard verbose mode. This way, all prints that would usually go to\
        stdout will be logged to a file on TaIGa's current folder and she will print to your screen\
        no more.",
        action="store_true")

    args = taiga.parse_args()

    # Define the path to the input file
    input_path = args.infile
    # Define the path to the output file, adding an ending forward slash if needed
    output_path = args.outdir + '/' if args.outdir[-1] != '/' else args.outdir
    # Define the user email for Entrez.email, as it is a good practice of E-utils
    user_email = args.email
    # Define maximum number of retries. TODO: maybe use the biopython native approach instead.
    retries = args.t[0]
    # True or false for verbose
    verbose = args.v
    # True or false for correction
    correction = args.c
    # True or false for taxon id list input
    tid = args.tid
    # Input mode for genbank format files
    mode = args.mode[0]

    # A list to hold all organism names
    taxon_list = []

    helpers.config_log(verbose)

    log.info("""
    *********************************************
    *                                           *
    *   TaIGa - Taxonomy Information Gatherer   *
    *                                           *
    *********************************************""")

    # Checking if TaIGa is being run on TaxID mode with the '-c' argument.
    # This is needed because, when run with '--tid', TaIGa never actually calls 'correct_spell'.
    # The retrieved name is assumed to be correct.
    if tid and correction:
        log.error("\nERROR: Please, when running TaIGa with the '--tid' option, don't use '-c', as\
                  TaIGa skips the name correction.")
        sys.exit()

    # Check if input mode is for a genbank format file or a text file and parse the file
    if not (mode == 0):
        taxon_list = parsers.parse_gb(input_path, mode)
    else:
        taxon_list = parsers.parse_txt(input_path, tid)

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
    elif type(taxon_list) == list:
        missing_corrected, \
            missing_taxid, \
            taxon_ids_collection, \
            genome_ids_collection = retrievers.retrieve_from_multiple_names(
                taxon_list,
                user_email,
                correction,
                retries)
        print(taxon_ids_collection)
    else:  # If there's only one record, or only one organism, a loop isn't needed
        missing_corrected, \
            missing_taxid, \
            taxon_ids_collection, \
            genome_ids_collection = retrievers.retrieve_from_single_name(
                taxon_list,
                user_email,
                correction,
                retries)

    try:
        helpers.sanitize_name_list(correction, tid,
                                   missing_corrected, missing_taxid,
                                   taxon_list)

        tax_info = retrievers.get_tax_info(
            user_email,
            taxon_ids_collection,
            retries)

        frame = data_handlers.create_df(
            taxon_ids_collection,
            genome_ids_collection,
            tax_info,
            tid,
            taxon_list)

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
    except (Exception):
        log.error("\nQUIT: Unknown error occured while generating output files.")
        sys.exit()

    log.info(
        "\n>> TaIGa was run successfully! You can check your results on the informed output folder.\
        If there's any missing data, check the 'TaIGa_missing.txt' file on the same folder.\n"
    )
