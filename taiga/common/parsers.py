import sys
import logging as log
from Bio import SeqIO
from .data_models import Taxon
from typing import List, Set


def parse_txt(input_path: str, tid: bool) -> List:
    """Parses a file containing multiple strings and return a list of objects

    Parameters:
    input_path (string): The path to the input file as a string
    tid (boolean): A variable informing if the input file contains Taxon IDs

    Returns:
    list_of_taxa (list): A list of Taxon objects, each containing a self.name value for each name
                         in the input file or a self.tid value for each Taxon ID in the input file

    """

    list_of_taxa: List = []

    try:

        with open(input_path, "r") as infile:
            # Using set to uniquefy possible duplicate names or IDs from the input
            all_inputs: Set = set(infile.readlines())

            # Checks first if input file is an ids or names file
            if tid:
                log.info("\n> Parsing input file a simple list of Taxon IDs\n")

                for each_id in all_inputs:
                    # Each taxon might contain a '\n' character at the end, so remove it
                    new_id = Taxon(taxon_id=int(each_id.replace("\n", "")))

                    list_of_taxa.append(new_id)

                    log.info(f"{new_id.taxon_id} ---> All OK")
            else:
                log.info("\n> Parsing input file a simple list of species names\n")

                for each_taxon in all_inputs:
                    # Each taxon might contain a '\n' character at the end, so remove it
                    new_taxon = Taxon(name=each_taxon.replace("\n", ""))

                    list_of_taxa.append(new_taxon)

                    log.info(f"{new_taxon.name} ---> All OK")
    except (KeyboardInterrupt):
        log.warning("\nQUIT: TaIGa was stopped by the user\n")
        sys.exit()
    except (Exception):
        log.error(
            "\nERROR: Couldn't parse input text file. Check your file an try running TaIGa again\n"
        )
        sys.exit()

    return list_of_taxa


def parse_gb(input_path: str, mode: int) -> List:
    """Parses a Genbank format file containing sequence records and organism information

    Parameters:
    input_path (string): The path to the input file as a string
    mode (int): An integer representing the reading mode for the input file. Options are:
                0: a plain text file, which is not handled by this function
                1: Genbank file with multiple records from different organisms
                2: Genbank file with a single record from a single organism
                3: Genbank file with multiple records from the same organism

    Returns:
    list_of_taxa (list): A list of Taxon objects, each containing a self.name value for each record
                  in the input file

    """

    list_of_taxa: List = []

    try:
        if mode == 1:
            log.info(
                "\n> Parsing input as a Genbank file with multiple records and organisms")

            input_records = SeqIO.parse(input_path, "genbank")

            # Using a set to uniquefy possible duplicate names or IDs from the input
            list_of_taxa = list(
                {Taxon(name=seq.annotations["organism"]) for seq in input_records})
        elif mode == 2:
            log.info("\n> Parsing input as a Genbank file with a single record")
            try:
                input_records = SeqIO.read(input_path, "genbank")

                list_of_taxa.append(
                    Taxon(name=input_records.annotations["organism"]))
            except (ValueError):
                # Catch an error if there's more than one record in the input file.
                log.error("\nERROR: The file contains more than one record. Try running TaIGa "
                          "again with a different Genbank file mode\n")
                sys.exit()
        elif mode == 3:
            log.info(
                "\n> Parsing input as a Genbank file with multiple records for one organism")

            input_records = SeqIO.parse(input_path, "genbank")

            for record in input_records:
                # Only take the first item of the generator. Is faster than conveting to a list.
                list_of_taxa.append(Taxon(name=record.annotations["organism"]))
                break
    except (KeyboardInterrupt):
        log.warning("\nQUIT: TaIGa was stopped by the user.\n")
        sys.exit()
    except(Exception):
        log.error(
            "\nERROR: Couldn't parse input genome file. Check your file and try running TaIGa "
            "again.\n")
        sys.exit()

    if list_of_taxa:
        log.info("\n> All OK with parsing the input file. Checking the records...\n")
        for taxon in list_of_taxa:
            log.info(f"> '{taxon.name}' ---> All OK")
    else:
        log.error("\nERROR: Something went wrong while trying to parse the input file. Check the "
                  "file and the program execution commands and try again\n")
        sys.exit()

    return list_of_taxa
