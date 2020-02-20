class Taxon:
    def __init__(self, name='', genome_id='N/A', taxon_id='N/A', classification=dict(),
                 missing_name=False, missing_taxon_id=False, missing_corrected=False):
        self.name = name
        self.genome_id = genome_id
        self.taxon_id = taxon_id
        self.classification = classification
        self.missing_name = missing_name
        self.missing_taxon_id = missing_taxon_id
        self.missing_corrected = missing_corrected


    def list_taxa(self):
        return list(self.classification.values())


    def list_ranks(self):
        return set(self.classification)
