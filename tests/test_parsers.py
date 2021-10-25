import sys
import os
import helpers

# Fix to allow relative imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from taiga.common import parsers
from taiga.common.data_models import Taxon


def test_parse_names_text_file():
    result = parsers.parse_txt("examples/input_names.txt", False)

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


if __name__ == "__main__":
    test_parse_names_text_file()
