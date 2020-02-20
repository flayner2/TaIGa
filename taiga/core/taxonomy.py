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
        nargs=1, default=5, type=int)
    taiga.add_argument(
        "-v",
        help="Turn off TaIGa's standard verbose mode. This way, all prints that would usually go to\
        stdout will be logged to a file on TaIGa's current folder and she will print to your screen\
        no more.",
        action="store_true")

    args = taiga.parse_args()
    
    input_path = args.infile
    # Add an ending forward slash to output path if needed
    output_path = args.outdir + '/' if args.outdir[-1] != '/' else args.outdir
    # Providing the email when doing requests through E-Utils is recommended
    user_email = args.email
    # This will be used to set the value for Entrez.max_tries
    retries = args.t
    verbose = args.v
    correction = args.c
    tid = args.tid
    mode = args.mode[-1]

    # A list to hold Taxon objects
    taxon_list = []

    # Inital configuration for the logging module.
    # At this point, the output may be set to verbose or not.
    helpers.config_log(verbose)

    log.info("""
    *********************************************
    *                                           *
    *   TaIGa - Taxonomy Information Gatherer   *
    *                                           *
    *********************************************""")

    # Checking if TaIGa is being run on TaxID mode with the '-c' argument.
    # This is needed because, when run with '--tid', TaIGa never actually calls 'correct_spell' as
    # the retrieved name is assumed to be correct.
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

    # Checking the type of input (taxon id or names) and fetching the rest of the information
    if tid:
        retrievers.retrieve_from_taxid(taxon_list, user_email, retries)
    else:
        retrievers.retrieve_from_names(taxon_list, user_email, correction, retries)

    retrievers.retrieve_taxonomy(taxon_list, user_email, retries)

    frame = data_handlers.create_df(taxon_list)

    data_handlers.create_output(output_path, frame, taxon_list)

    log.info(
        "\n>> TaIGa was run successfully! You can check your results on the informed output folder.\
        If there's any missing data, check the 'TaIGa_missing.txt' file on the same folder.\n"
    )
