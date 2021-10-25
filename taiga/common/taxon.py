from typing import Dict, List, Set


class Taxon:
    """
    A Taxon is a representation of the taxonomic information of an organism.
    """

    name: str
    genome_id: int
    taxon_id: int
    classification: Dict  # This may become a `defaultdict`
    missing_name: bool
    missing_taxon_id: bool
    missing_corrected: bool
    missing_classification: bool

    def __init__(
        self,
        name: str = "",
        genome_id: int = 0,
        taxon_id: int = 0,
        classification: Dict = {},
        missing_name: bool = False,
        missing_taxon_id: bool = False,
        missing_corrected: bool = False,
        missing_classification: bool = False,
    ) -> None:
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

        if self.classification is not None:
            return list(self.classification.values())
        else:
            return []

    def list_ranks(self) -> Set:
        """Returns a set with all the ranks for an organism

        Parameters:
        self (object): A reference to the object itself

        Returns:
        (set): A set of all the unique keys from the 'self.classification' dict

        """

        if self.classification is not None:
            return set(self.classification)
        else:
            return set()
