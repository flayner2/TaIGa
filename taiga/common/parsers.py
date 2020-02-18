import sys
import logging as log
from Bio import SeqIO
from taiga.common.data_models import Taxon


def parse_txt(input_path, tid):
    """Parses a file containing multiple strings and return a list of objects.

    Parameters:
    input_path (string): The path to the input file as a string.
    tid (boolean): A variable informing if the input file contains TaxonIDs.

    Returns:
    list_of_taxa (list): A list of Taxon objects, each containing a self.name value for each name
    in the input file or a self.tid value for each TaxonID in the input file.

    """

    list_of_taxa = []

    try:
        log.info("\n>> Parsing input file a simple list of species names.\n")

        with open(input_path, "r") as infile:
            all_inputs = infile.readlines()

            # Checks first if input file is an ids or names file
            if tid:
                for each_id in all_inputs:
                    list_of_taxa.append(Taxon(t_id=each_id))

                    log.info("{} ---> All OK".format(each_id))
            else:
                for each_taxon in all_inputs:
                    # Each taxon might contain a '\n' character at the end, so remove it
                    new_taxon = Taxon(name=each_taxon.replace('\n', ''))

                    list_of_taxa.append(new_taxon)

                    log.info("{} ---> All OK".format(each_taxon))
    except (KeyboardInterrupt):
        log.warning("\nQUIT: TaIGa was stopped by the user.\n")
        sys.exit()
    except (Exception): # TODO: improve generic exception handling
        log.error(
            "\nERROR: Couldn't parse input text file. Check your file an try running TaIGa again.\n"
        )
        sys.exit()

    return list_of_taxa


def parse_gb(input_path, mode):
    """Parses a genbank format file containing sequence records and organism information.

    Parameters:
    input_path (string): The path to the input file as a string.
    mode (int): An integer representing the reading mode for the input file. Options are:
        1: genbank file with multiple records from different organisms
        2: genbank file with a single record from a single organism
        3: genbank file with multiple records from a single organism

    Returns:
    list_of_taxa: A list of Taxon objects, each containing a self.name value for each record
    in the input file.

    """

    list_of_taxa = []

    try:
        if mode == 1:
            log.info("\n>> Parsing input as a Genbank file with multiple records for multiple organisms.\n")

            input_records = SeqIO.parse(input_path, "genbank")

            list_of_taxa = [Taxon(name=seq.annotations["organism"]) for seq in input_records]
        elif mode == 2:
            log.info("\n>> Parsing input as a Genbank file with a single record.\n")
            try:
                input_records = SeqIO.read(input_path, "genbank")

                list_of_taxa.append(Taxon(name=input_records.annotations['organism']))
            except (ValueError):
                # Catch an error if there's more than one record in the input file.
                log.error("\nERROR: The file contains more than one record. Try running TaIGa again without the '--single' option.\n")
                sys.exit()
        elif mode == 3:
            log.info("\n>> Parsing input as a Genbank file with multiple records for the same organism.\n")

            input_records = SeqIO.parse(input_path, "genbank")

            for record in input_records:
                # Only take the first item of the generator. Is faster than conveting to a list.
                list_of_taxa.append(Taxon(name=record.annotations['organism']))
                break
    except (KeyboardInterrupt):
        log.warning("\nQUIT: TaIGa was stopped by the user.\n")
        sys.exit()
    except(Exception):
        log.error(
            "\nERROR: Couldn't parse input genome file. Check your file and try running TaIGa again.\n"
        )
        sys.exit()


    if list_of_taxa:
        log.info(
            "\n>> All OK with parsing the input file. Checking the records...\n"
        )
        for taxon in list_of_taxa:
            log.info("> '{}' ---> All OK".format(taxon.name))
    else:
        log.error("\nERROR: Something went wrong while trying to parse the input file. Check the file and the program execution commands and try again.\n")
        sys.exit()

    return list_of_taxa
