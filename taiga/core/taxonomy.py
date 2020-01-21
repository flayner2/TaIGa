import sys
import logging as log
from taiga.common import parsers, helpers, global_vars, retrievers, data_handlers


def run_taiga():
    """ Wrapper for all of TaIGa's main functionalities

    Gets the user input, parses the input files and run all the needed
    logic and functions, producing the output files at the specified
    folder.
    """

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

    helpers.config_log(global_vars.verbose)

    log.info("""
    *********************************************
    *                                           *
    *   TaIGa - Taxonomy Information Gatherer   *
    *                                           *
    *********************************************""")

    # Checking if only one input format argument was passed to TaIGa
    helpers.sanitize_args(global_vars.args)

    # Checking if TaIGa is being run on TaxID mode with the '-c' argument.
    # This is needed because, when run with '--tid', TaIGa never actually calls 'correct_spell'.  
    # The retrieved name is assumed to be correct.
    if global_vars.tid and global_vars.correction:
        log.error("\nERROR: Please, when running TaIGa with the '--tid' option, don't use '-c', as TaIGa skips the name correction.") 
        sys.exit()

    # First checking if input file is a text file with multiple TaxIDs
    if global_vars.tid:
        input_taxid = parsers.parse_taxid(global_vars.input_path)
    # Then checking if input file is a genome file with multiple records from multiple organisms
    elif global_vars.multi:
        input_names = parsers.parse_gb_multi(global_vars.input_path)
    # Else, checking the type of input Genbank file
    elif global_vars.single:  # Single record on Genbank file
        input_names = parsers.parse_gb_single(global_vars.input_path)
    elif global_vars.same:  # Multiple records from the same organism
        input_names = parsers.parse_gb_same(global_vars.input_path)
    else:  # Simple text file with a list of organism names
        input_names = parsers.parse_names(global_vars.input_path)

    log.info("\n>> Ignore if any duplicate name was printed, duplicates are removed.")

    log.info("\n>> Searching for taxonomic information...\n")

    # Checking if the input is a list of TaxIDs instead of any type of name input
    if global_vars.tid:
        taxon_ids_collection, \
            genome_ids_collection, \
            missing_name = retrievers.retrieve_from_taxid(
                input_taxid,
                global_vars.user_email,
                global_vars.retries)
    # Checking if there are multiple records on the file (for different organisms)
    elif type(input_names) == list:
        missing_corrected, \
            missing_taxid, \
            taxon_ids_collection, \
            genome_ids_collection = retrievers.retrieve_from_multiple_names(
                input_names,
                global_vars.user_email,
                global_vars.correction,
                global_vars.retries)
        print(taxon_ids_collection)
    else:  # If there's only one record, or only one organism, a loop isn't needed
        missing_corrected, \
            missing_taxid, \
            taxon_ids_collection, \
            genome_ids_collection = retrievers.retrieve_from_single_name(
                input_names,
                global_vars.user_email,
                global_vars.correction,
                global_vars.retries)

    try:
        helpers.sanitize_name_list(global_vars.correction, global_vars.tid,
                                   missing_corrected, missing_taxid,
                                   input_names)

        tax_info = retrievers.get_tax_info(
            global_vars.user_email,
            taxon_ids_collection,
            global_vars.retries)

        frame = data_handlers.create_df(
            taxon_ids_collection,
            genome_ids_collection,
            tax_info,
            global_vars.tid,
            input_names)

        data_handlers.create_output(
            global_vars.output_path,
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
