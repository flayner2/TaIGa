import os
import sys
import pandas as pd
import logging as log
from taiga.common import fetchers, parsers, helpers, global_vars, retrievers


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
        taxon_ids_collection,
        genome_ids_collection,
        missing_name = retrievers.retrieve_from_taxid(
            input_taxid,
            global_vars.user_email,
            global_vars.retries)
    # Checking if there are multiple records on the file (for different organisms)
    elif type(input_names) == list:
        missing_corrected,
        missing_taxid,
        taxon_ids_collection,
        genome_ids_collection = retrievers.retrieve_from_multiple_names(
            input_names,
            global_vars.user_email,
            global_vars.correction,
            global_vars.retries)
    else:  # If there's only one record, or only one organism, a loop isn't needed
        missing_corrected,
        missing_taxid,
        taxon_ids_collection,
        genome_ids_collection = retrievers.retrieve_from_single_name(
            input_names,
            global_vars.user_email,
            global_vars.correction,
            global_vars.retries)

    try:
        if (global_vars.correction) and (not global_vars.tid):
            # Check for the organisms with missing corrected name or taxid and remove their name form the names list
            for missed in missing_corrected:
                if missed in input_names:
                    input_names.remove(missed)
        if not global_vars.tid:    
            for missed in missing_taxid:
                if missed in input_names:
                    input_names.remove(missed)

        # Done getting all preliminary information, print a message informing the next step
        log.info(
            ">> Gathering taxonomic information from NCBI Taxonomy. This might take a while.\n"
        )

        # Creating the list of dictionaries to store the taxonomic information for the organisms.
        tax_info = fetchers.organize_tax_info(global_vars.user_email, taxon_ids_collection, global_vars.retries)

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
        if global_vars.tid:
            names_list = [i for i in taxon_ids_collection.keys()]
            frame = pd.DataFrame(final_info, index=names_list, columns=ranks)
        elif type(names) == list:  # Checking if names is a list
            frame = pd.DataFrame(final_info, index=names, columns=ranks)
        else:  # If it's not, make it a list so Pandas can use it as indexes for the rows of the DataFrame
            names_list = []
            names_list.append(names)
            frame = pd.DataFrame(final_info, index=names_list, columns=ranks)

        # Transform the missing data to a more visual indicator for the user
        frame.fillna('-', inplace=True)

        log.info(
            "\n>> Done generating result table. Checking if output folder exists")

        # Check if the output directory exists. If not, create it.
        if not os.path.exists(global_vars.output_path):
            log.info("\n>> Creating output folder")
            try:  # Checking if the path for the output folder is valid
                os.makedirs(global_vars.output_path)
            except:
                log.error(
                    "\nERROR: Path to output folder may not be valid. Try again.")
                sys.exit()

        log.info(
            "\n>> Creating output file. You'll find it inside the provided output folder, named 'TaIGa_result.csv'"
        )
        # Export the DataFrame to the resulting .csv file
        frame.to_csv(global_vars.output_path + 'TaIGa_result.csv')

        # Checking if there are missing correct names or TaxIDs. If there are, generating log files for those.
        if missing_corrected or missing_taxid or missing_name:
            log.info(
                "\n>> Creating a file for the organisms with missing information. You'll find it inside the provided output folder, named 'TaIGa_missing.txt'"
            )
            with open(global_vars.output_path + 'TaIGa_missing.txt', 'w') as missing_file:
                missing_file.write("Missing corrected names: \n")
                if (missing_corrected):
                    for name in missing_corrected:
                        missing_file.write("\t\t\t{}\n".format(name))
                missing_file.write("Missing TaxID: \n")
                if (missing_taxid):
                    for taxid in missing_taxid:
                        missing_file.write("\t\t\t{}\n".format(taxid))
                missing_file.write("TaxIDs with missing names: \n")
                if (missing_name):
                    for tax_id in missing_name:
                        missing_file.write("\t\t\t{}\n".format(tax_id))
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
