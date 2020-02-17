class Taxon:
    def __init__(self, name, g_id='N/A', t_id='N/A', classification={}):
        name = self.name
        g_id = self.g_id
        t_id = self.t_id
        classification = self.classification

        # A list of all the ranks present in the classification of an organism
        ranks = set(classification)

    def list_taxa(self):
        return list(self.classification.values())
