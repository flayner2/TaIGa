import sys
import logging as log
from Bio import SeqIO


def parse_names(input_path):
    try:
        log.info("\n>> Parsing input file a simple list of species names.\n")
        with open(input_path, "r") as input_names:
            names = input_names.readlines()  # Take the names of the organisms
            for each_name in range(len(names)):
                names[each_name] = names[each_name].replace("\n", "")  # Correct the formatting of each name
                log.info("{} ---> All OK".format(names[each_name]))
    except (KeyboardInterrupt):
        log.warning("\nQUIT: TaIGa was stopped by the user.\n")
        sys.exit()
    except (Exception):
        log.error(
            "\nERROR: Couldn't parse name list. Check your file an try running TaIGa again.\n"
        )
        sys.exit()

    return names


def parse_taxid(input_path):
    try:
        log.info("\n>> Parsing input as a text file with multiple TaxIDs.\n")
        with open(input_path, "r") as input_ids:
            ids = input_ids.readlines()
            for each_id in range(len(ids)):
                ids[each_id] = ids[each_id].replace("\n", "")
                log.info("{} ---> All OK".format(ids[each_id]))
    except (KeyboardInterrupt):
        log.warning("\nQUIT: TaIGa was stopped by the user.\n")
        sys.exit()
    except (Exception):
        log.error(
            """\nERROR: Couldn't parse name list.\
            Check your file an try running TaIGa again.\n"""
        )
        sys.exit()

    return ids


def parse_gb_multi(input_path):
    try:
        log.info("\n>> Parsing input as a Genbank file with multiple records for multiple organisms.\n")
        input_multi_records = SeqIO.parse(input_path, "genbank")
    except (KeyboardInterrupt):
        log.warning("\nQUIT: TaIGa was stopped by the user.\n")
        sys.exit()
    except(Exception):
        log.error(
            "\nERROR: Couldn't parse input genome file. Check your file and try running TaIGa again.\n"
        )
        sys.exit()

    multi_records = [seq.annotations["organism"] for seq in input_multi_records]  # List the names of all organisms

    if multi_records:
        log.info(
            "\n>> All OK with parsing the input file. Checking the records...\n"
        )
        for record in multi_records:
            log.info("> '{}' ---> All OK".format(record))   
    else:
        log.error("\nERROR: Something went wrong while trying to parse the input file. Check the file and the program execution commands and try again.\n")
        sys.exit()

    return multi_records


def parse_gb_single(input_path):
    try:
        log.info("\n>> Parsing input as a Genbank file with a single record.\n")
        
        try:
            input_single_record = SeqIO.read(input_path, "genbank")
        except (KeyboardInterrupt):
            log.warning("\nQUIT: TaIGa was stopped by the user.\n")
            sys.exit()
        except (Exception):
            log.error("\nERROR: Couldn't parse input genome file. Check your file and try running TaIGa again.\n")
            sys.exit()
        
        single_record = input_single_record.annotations["organism"]  # Take the name of the organism on the input file
        
        if single_record:
            log.info("> '{}' ---> All OK".format(single_record))
        else:
            log.error("\nERROR: Something went wrong while trying to parse the input file. Check the file and the program execution commands and try again.\n")
            sys.exit()
    except (KeyboardInterrupt):
        log.warning("\nQUIT: TaIGa was stopped by the user.\n")
        sys.exit()
    # Catch an error if the file contains more than one record
    except (ValueError):
        log.error("\nERROR: The file contains more than one record. Try running TaIGa again without the '--single' option.\n")
        sys.exit()

    return single_record


def parse_gb_same(input_path):
    log.info("\n>> Parsing input as a Genbank file with multiple records for the same organism.\n")
    
    try:
        input_same_records = list(SeqIO.parse(input_path, "genbank"))
    except (KeyboardInterrupt):
        log.warning("\nQUIT: TaIGa was stopped by the user.\n")
        sys.exit()
    except(Exception):
        log.error(
            "\nERROR: Couldn't parse input genome file. Check your file and try running TaIGa again.\n"
        )
        sys.exit()
    same_records = input_same_records[0].annotations["organism"]  # Only take the first taxon name retrieved, as all are from the same organism
    if same_records:
        log.info("> '{}' ---> All OK".format(same_records))
    else:
        log.error("\nERROR: Something went wrong while trying to parse the input file. Check the file and the program execution commands and try again.\n")
        sys.exit()

    return same_records
