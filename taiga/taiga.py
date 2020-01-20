"""                                          TaIGa - TAxonomy Information GAtherer
This is a simple script that interacts with various utilities from the NCBI's Entrez api in order to retrieve relevant taxonomic
information for a collection of organisms. As of now, TaIGa is able to handle multiple types of Genbank format genome files, as
well as a text file format list of organism names or taxon IDs, separated by lines. TaIGa recieves a file as input, an output 
folder path, a valid user e-mail and one optional argument to identify the type of input file. TaIGa then uses Entrez to 
retrieve the TaxID, Genome ID and all taxonomic information of all taxa up to the organism name provided. Then, it builds a 
DataFrame and outputs it to a .csv file, so the user can visualize it as a table. TaIGa's goal is to make easier for researchers 
to gather mass taxonomical metadata for their projects. Therefore, TaIGa is best used when you have a big list of organisms or a 
big collection of genomes in a file. TaIGa is also a very cute anime character from the japanese romance animation ToraDora. You 
should watch it. """ 

from taiga.common import helpers
from taiga.core import taxonomy


def run():
    helpers.sanitize_version()
    taxonomy.run_taiga()
