import pandas as pd
import logging as log
import sys
import os


def create_df(taxon_ids_collection, genome_ids_collection, tax_info, tid, input_names):
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
    if tid:
        names_list = [i for i in taxon_ids_collection.keys()]
        frame = pd.DataFrame(final_info, index=names_list, columns=ranks)
    elif type(input_names) == list:  # Checking if names is a list
        frame = pd.DataFrame(final_info, index=input_names, columns=ranks)
    else:  # If it's not, make it a list so Pandas can use it as indexes for the rows of the DataFrame
        names_list = []
        names_list.append(input_names)
        frame = pd.DataFrame(final_info, index=names_list, columns=ranks)

    # Transform the missing data to a more visual indicator for the user
    frame.fillna('-', inplace=True)

    log.info(
        "\n>> Done generating result table")

    return frame


def create_output(output_path, frame, missing_corrected, missing_taxid, missing_name):
    log.info(
        "\n>> Checking if output folder exists")
    # Check if the output directory exists. If not, create it.
    if not os.path.exists(output_path):
        log.info("\n>> Creating output folder")
        try:  # Checking if the path for the output folder is valid
            os.makedirs(output_path)
        except (Exception):
            log.error(
                "\nERROR: Path to output folder may not be valid. Try again.")
            sys.exit()

    log.info(
        "\n>> Creating output file. You'll find it inside the provided output folder, named 'TaIGa_result.csv'"
    )
    # Export the DataFrame to the resulting .csv file
    frame.to_csv(output_path + 'TaIGa_result.csv')

    # Checking if there are missing correct names or TaxIDs. If there are, generating log files for those.
    if missing_corrected or missing_taxid or missing_name:
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
            missing_file.write("TaxIDs with missing names: \n")
            if (missing_name):
                for tax_id in missing_name:
                    missing_file.write("\t\t\t{}\n".format(tax_id))
