class Taxon:
    def __init__(self, name='', g_id='N/A', t_id='N/A', classification={},
                 missing_name=False, missing_taxid=False, missing_corrected=False):
        self.name = name
        self.g_id = g_id
        self.t_id = t_id
        self.classification = classification
        self.missing_name = missing_name
        self.missing_taxid = missing_taxid
        self.missing_corrected = missing_corrected

        # A list of all the ranks present in the classification of an organism
        self.ranks = set(self.classification)

    def list_taxa(self):
        return list(self.classification.values())
