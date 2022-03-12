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
        if self.has_taxon_id() and other.has_taxon_id():
            return self.taxon_id < other.taxon_id
        elif self.has_name() and other.has_name():
            return self.name < other.name
        elif self.has_genome_id() and other.has_genome_id():
            return self.genome_id < other.genome_id

        return False

    def __le__(self, other) -> bool:
        if self.has_taxon_id() and other.has_taxon_id():
            return self.taxon_id <= other.taxon_id
        elif self.has_name() and other.has_name():
            return self.name <= other.name
        elif self.has_genome_id() and other.has_genome_id():
            return self.genome_id <= other.genome_id

        return False

    def __gt__(self, other) -> bool:
        if self.has_taxon_id() and other.has_taxon_id():
            return self.taxon_id > other.taxon_id
        elif self.has_name() and other.has_name():
            return self.name > other.name
        elif self.has_genome_id() and other.has_genome_id():
            return self.genome_id > other.genome_id

        return False

    def __ge__(self, other) -> bool:
        if self.has_taxon_id() and other.has_taxon_id():
            return self.taxon_id >= other.taxon_id
        elif self.has_name() and other.has_name():
            return self.name >= other.name
        elif self.has_genome_id() and other.has_genome_id():
            return self.genome_id >= other.genome_id

        return False

    def __eq__(self, other) -> bool:
        if self.has_taxon_id() and other.has_taxon_id():
            return self.taxon_id == other.taxon_id
        elif self.has_name() and other.has_name():
            return self.name == other.name
        elif self.has_genome_id() and other.has_genome_id():
            return self.genome_id == other.genome_id

        return False

    def __ne__(self, other) -> bool:
        if self.has_taxon_id() and other.has_taxon_id():
            return self.taxon_id != other.taxon_id
        elif self.has_name() and other.has_name():
            return self.name != other.name
        elif self.has_genome_id() and other.has_genome_id():
            return self.genome_id != other.genome_id

        return False

    def __str__(self) -> str:
        return (
            f"Name: {self.name}\n"
            f"TaxonID: {self.taxon_id if self.has_taxon_id() else 'N/A'}\n"
            f"GenomeID: {self.genome_id if self.has_genome_id() else 'N/A'}\n"
            f"Classification: {self.get_taxon_names_list()}"
        )

    ###########################################################

    # Getters and setters - only for compound data like dicts and lists

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
