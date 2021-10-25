from typing import Dict, List, Set


class Taxon:
    """
    A Taxon is a representation of the taxonomic information of an organism.
    """

    # Class variable declarations

    name: str
    genome_id: int
    taxon_id: int
    classification: Dict  # This may become a `defaultdict`

    ###########################################################

    # Dunder methods

    def __init__(
        self,
        name: str = "",
        genome_id: int = 0,
        taxon_id: int = 0,
        classification: Dict = {},
    ):
        self.name = name
        self.genome_id = genome_id
        self.taxon_id = taxon_id
        self.classification = classification

    def __lt__(self, other) -> bool:
        pass
        # if self.taxon_id > 0 and other.taxon_id > 0:
        # return self.taxon_id < other.taxon_id
        # elif self.name and other.name:
        # return self.name < other.name
        # elif self.genome_id > 0 and

    ###########################################################

    # Getters and setters

    def get_taxon_names_list(self) -> List:
        """
        @getter

        Returns a list of all the names of the taxon ranks for an organism.

        Returns:
            (list): A list of all the values from the `self.classification` dict.
        """

        if self.classification:
            return list(self.classification.values())

        return []

    def get_taxon_ranks_list(self) -> Set:
        """
        @getter

        Returns a set with all the ranks for an organism.

        Parameters:
            self (object): A reference to the object itself.

        Returns:
            (set): A set of all the unique keys from the `self.classification` dict.

        """

        if self.classification is not None:
            return set(self.classification)
        else:
            return set()

    ###########################################################

    # Utility methods

    def has_name(self) -> bool:
        """
        Checks if the Taxon has an organism name.

        Returns:
            (bool): True if the Taxon has a name, else False.
        """

        return bool(self.name)

    def has_taxon_id(self) -> bool:
        """
        Checks if the Taxon has a TaxonID (meaning it is > 0).

        Returns:
            (bool): True if the Taxon has a TaxonID, else False.
        """

        return bool(self.taxon_id)

    def has_genome_id(self) -> bool:
        """
        Checks if the Taxon has an genome id (meaning it is > 0).

        Returns:
            (bool): True if the Taxon has a genome id, else False.
        """

        return bool(self.genome_id)

    def has_classification(self) -> bool:
        """
        Checks if the Taxon has a classification dictionary.

        Returns:
            (bool): True if the Taxon has a classification, else False.
        """

        return bool(self.classification)
