# A Taxon holds all the information parsed an fetched by TaIGa for one organism
class Taxon:
    def __init__(self, name="", genome_id="N/A", taxon_id="N/A", classification=dict(),
                 missing_name=False, missing_taxon_id=False, missing_corrected=False,
                 missing_classification=False):
        self.name = name
        self.genome_id = genome_id
        self.taxon_id = taxon_id
        self.classification = classification
        self.missing_name = missing_name
        self.missing_taxon_id = missing_taxon_id
        self.missing_corrected = missing_corrected
        self.missing_classification = missing_classification

    def list_taxa(self):
        """Returns a list of all the names of the taxon ranks for an organism

        Parameters:
        self (object): A reference to the object itself

        Returns:
        (list): A list of all the values from the 'self.classification' dict

        """

        return list(self.classification.values())

    def list_ranks(self):
        """Returns a set with all the ranks for an organism

        Parameters:
        self (object): A reference to the object itself

        Returns:
        (set): A set of all the unique keys from the 'self.classification' dict

        """

        return set(self.classification)
