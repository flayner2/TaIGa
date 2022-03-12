import os
import sys

import helpers

# Fix to allow relative imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from taiga.common import parsers
from taiga.common.taxon import Taxon


def test_parse_names_text_file():
    result = parsers.parse_txt("../examples/input_names.txt", False)

    # Order doesn't matter since `parse_txt` messes the original ordering by
    # converting the list to a set first to remove duplicates.
    expected = [
        Taxon("Homo sapiens"),
        Taxon("Apis mellifera"),
        Taxon("Bombus impatiens"),
        Taxon("Bombus terrestris"),
        Taxon("Gorilla gorilla"),
        Taxon("Escherichia coli"),
        Taxon("Caenorhabditis elegans"),
        Taxon("Mus musculus"),
        Taxon("Equus caballus"),
        Taxon("Bos taurus"),
    ]

    # Assert letting the helper function force the same order for both lists.
    assert helpers.batch_assert_objects(result, expected, "name", False)


def test_parse_ids_text_file():
    result = parsers.parse_txt("../examples/input_taxon_ids.txt", True)

    # Order doesn't matter since `parse_txt` messes the original ordering by
    # converting the list to a set first to remove duplicates.
    expected = [
        Taxon(taxon_id=9606),
        Taxon(taxon_id=7460),
        Taxon(taxon_id=132113),
        Taxon(taxon_id=30195),
        Taxon(taxon_id=9593),
        Taxon(taxon_id=562),
        Taxon(taxon_id=6239),
        Taxon(taxon_id=10090),
        Taxon(taxon_id=9796),
        Taxon(taxon_id=9913),
    ]

    # Assert letting the helper function force the same order for both lists.
    assert helpers.batch_assert_objects(result, expected, "taxon_id", False)


if __name__ == "__main__":
    # Debugging
    test_parse_names_text_file()
