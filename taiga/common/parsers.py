import sys
import logging as log
from Bio import SeqIO
from .taxon import Taxon
from .errors import EmptyFileError


def parse_txt(input_path: str, tid: bool) -> list[Taxon]:
    """Parses a file containing multiple strings and return a list of objects

    Parameters:
    input_path (string): The path to the input file
    tid (boolean): A variable informing if the input file contains Taxon IDs

    Returns:
    parsed_taxa (list): A list of Taxon objects, each containing a self.name value for each name
                        in the input file or a self.tid value for each Taxon ID in the input file
    """

    parsed_taxa = []

    try:
        with open(input_path, "r") as infile:
            # Using set to uniquefy possible duplicate names or IDs from the input,
            # with a list comprehension to remove newline characters and convert
            # all names to lowercase
            all_inputs = set(
                [each_line.strip().lower() for each_line in infile.readlines()]
            )

            if len(all_inputs) <= 0:
                raise EmptyFileError(input_path)

            if tid:
                log.info("\n> Parsing input file as a list of Taxon IDs...\n")

                for each_id in all_inputs:
                    try:
                        new_id = int(each_id)
                        new_taxon = Taxon(taxon_id=new_id)
                        parsed_taxa.append(new_taxon)
                    except ValueError:
                        log.error(
                            f"\nERROR: Taxon ID '{each_id}' is possibly not a valid id."
                        )
                        sys.exit(1)
            else:
                log.info("\n> Parsing input file as a list of taxon names...\n")

                for each_name in all_inputs:
                    new_taxon = Taxon(name=each_name)
                    parsed_taxa.append(new_taxon)
    except KeyboardInterrupt:
        log.warning("\nQUIT: TaIGa was stopped by the user\n")
        sys.exit(1)
    except FileNotFoundError:
        log.error(f"\nERROR: File at '{input_path}' doensn't exist.")
        sys.exit(1)
    except EmptyFileError:
        log.error(f"\nERROR: File at '{input_path}' is empty.")
        sys.exit(1)
    except Exception as e:
        log.error(f"\nERROR: {e}\n")
        sys.exit(1)

    log.info("> Done parsing input file.\n")

    return parsed_taxa


def parse_gb(input_path: str, mode: int) -> List[Taxon]:
    """Parses a Genbank format file containing sequence records and organism information

    Parameters:
    input_path (string): The path to the input file as a string
    mode (int): An integer representing the reading mode for the input file. Options are:
                0: a plain text file, which is not handled by this function
                1: Genbank file with multiple records from different organisms
                2: Genbank file with a single record from a single organism
                3: Genbank file with multiple records from the same organism

    Returns:
    parsed_taxa (list): A list of Taxon objects, each containing a self.name value for each record
                  in the input file

    """

    parsed_taxa: List[Taxon] = []

    try:
        if mode == 1:
            log.info(
                "\n> Parsing input as a Genbank file with multiple records and organisms"
            )

            input_records = SeqIO.parse(input_path, "genbank")

            # Using a set to uniquefy possible duplicate names or IDs from the input
            parsed_taxa = list(
                {Taxon(name=seq.annotations["organism"]) for seq in input_records}
            )
        elif mode == 2:
            log.info("\n> Parsing input as a Genbank file with a single record")
            try:
                input_records = SeqIO.read(input_path, "genbank")

                parsed_taxa.append(Taxon(name=input_records.annotations["organism"]))
            except (ValueError):
                # Catch an error if there's more than one record in the input file.
                log.error(
                    "\nERROR: The file contains more than one record. Try running TaIGa "
                    "again with a different Genbank file mode\n"
                )
                sys.exit()
        elif mode == 3:
            log.info(
                "\n> Parsing input as a Genbank file with multiple records for one organism"
            )

            input_records = SeqIO.parse(input_path, "genbank")

            for record in input_records:
                # Only take the first item of the generator. Is faster than conveting to a list.
                parsed_taxa.append(Taxon(name=record.annotations["organism"]))
                break
    except (KeyboardInterrupt):
        log.warning("\nQUIT: TaIGa was stopped by the user.\n")
        sys.exit()
    except (Exception):
        log.error(
            "\nERROR: Couldn't parse input genome file. Check your file and try running TaIGa "
            "again.\n"
        )
        sys.exit()

    if parsed_taxa:
        log.info("\n> All OK with parsing the input file. Checking the records...\n")
        for taxon in parsed_taxa:
            log.info(f"> '{taxon.name}' ---> All OK")
    else:
        log.error(
            "\nERROR: Something went wrong while trying to parse the input file. Check the "
            "file and the program execution commands and try again\n"
        )
        sys.exit()

    return parsed_taxa
