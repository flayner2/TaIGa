"""
A Taxon holds all the information parsed an fetched by TaIGa for one organism
"""
from typing import Dict, List, Set


class Taxon:
    def __init__(self, name: str = None, genome_id: int = None, taxon_id: int = None, classification: Dict = dict(),
                 missing_name: bool = False, missing_taxon_id: bool = False, missing_corrected: bool = False,
                 missing_classification: bool = False) -> None:
        self.name = name
        self.genome_id = genome_id
        self.taxon_id = taxon_id
        self.classification = classification
        self.missing_name = missing_name
        self.missing_taxon_id = missing_taxon_id
        self.missing_corrected = missing_corrected
        self.missing_classification = missing_classification

    def list_taxa(self) -> List:
        """Returns a list of all the names of the taxon ranks for an organism

        Parameters:
        self (object): A reference to the object itself

        Returns:
        (list): A list of all the values from the 'self.classification' dict

        """

        return list(self.classification.values())

    def list_ranks(self) -> Set:
        """Returns a set with all the ranks for an organism

        Parameters:
        self (object): A reference to the object itself

        Returns:
        (set): A set of all the unique keys from the 'self.classification' dict

        """

        return set(self.classification)
