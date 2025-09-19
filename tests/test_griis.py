# tests/test_griis.py
import pytest
from b3alien import griis
from b3alien import utils
import pandas as pd
from b3alien.griis import CheckList

def test_griis_checklist():
    checklist = CheckList("tests/data/dwca-griis-portugal-v1.9/merged_distr.txt")
    assert hasattr(checklist, "species")
    assert isinstance(checklist.species, list)
    assert len(checklist.species) > 0

def test_matching():
    species = griis.get_speciesKey("Solanum nigrum")
    assert isinstance(species, list)
    assert len(species) > 0

def test_get_species_under_genus():
    species = griis.get_species_under_genus(9818805)
    assert isinstance(species, list)
    assert len(species) > 0

def test_split_event_date():
    date = "2001/2020"
    intro, outro = griis.split_event_date(date)
    assert intro == "2001"
    assert outro == "2020"

def test_split_event_date_invalid():
    date = "invalid_date"
    intro, outro = griis.split_event_date(date)
    assert intro is None or intro != intro  # Check for NaN
    assert outro is None or outro != outro  # Check for NaN

def test_split_event_date_non_string():
    date = 20200101
    intro, outro = griis.split_event_date(date)
    assert intro is None or intro != intro  # Check for NaN
    assert outro is None or outro != outro  # Check for NaN

def test_get_speciesKey_unresolvable():
    species = griis.get_speciesKey("Unresolvable species name")
    assert species == ["Uncertain"]

def test_get_speciesKey_genus():
    species = griis.get_speciesKey("Solanum")
    assert isinstance(species, list)
    assert len(species) > 0

def test_runtime_detection():
    assert utils.in_ipython() is True or utils.in_ipython() is False  
    assert utils.in_jupyter() is False or utils.in_jupyter() is True  
    assert utils.in_script() is False or utils.in_script is True  # Just ensure it runs without error

def make_distribution_file(tmp_path, rows):
    """
    Helper to write a tab-separated file with a 'speciesKey' column.
    `rows` is an iterable of values to put in the speciesKey column.
    Values can be:
      - Python lists (we'll stringify so they round-trip through literal_eval)
      - ints
      - None / float('nan')
    """
    # Ensure lists become proper Python-list strings (e.g., "[1, 2, 'Uncertain']")
    def to_cell(v):
        if isinstance(v, list):
            # Keep strings quoted so literal_eval parses correctly
            return str(v)
        return v

    df = pd.DataFrame({"speciesKey": [to_cell(v) for v in rows]})
    path = tmp_path / "distribution.txt"
    df.to_csv(path, sep="\t", index=False)
    return path


def test_load_griis_parses_lists_ints_and_filters_uncertain_and_nulls(tmp_path):
    # Rows include:
    # - list with ints
    # - list with repeated value + 'Uncertain' + None
    # - a single int (as string "789" in file -> literal_eval -> 789)
    # - a NaN cell
    rows = [
        [123, 456],
        [456, "Uncertain", None],
        789,
        float("nan"),
    ]
    file_path = make_distribution_file(tmp_path, rows)

    cl = CheckList(str(file_path))

    # Should keep unique integer species, drop 'Uncertain' and nulls, and deduplicate
    expected = {123, 456, 789}
    assert set(cl.species) == expected
    # Ensure everything is int
    assert all(isinstance(x, int) for x in cl.species)


def test_add_species_deduplicates_and_preserves_existing(tmp_path):
    # Minimal file to start with a couple of species
    file_path = make_distribution_file(tmp_path, [[1, 2]])
    cl = CheckList(str(file_path))

    cl._add_species([2, 3, 4, 3])
    assert cl.species == [1, 2, 3, 4]


def test_remove_species_filters_multiple_entries(tmp_path):
    file_path = make_distribution_file(tmp_path, [[1, 2, 3, 4]])
    cl = CheckList(str(file_path))

    cl._remove_species([2, 4, 999])  # 999 not present; should be no-op for that id
    assert cl.species == [1, 3]


def test_mixed_encoded_values_round_trip(tmp_path):
    # Include:
    # - list encoded as string with an 'Uncertain'
    # - a single int cell
    # - None in list (becomes None after literal_eval)
    # - duplicate values across rows
    rows = [
        [10, "Uncertain", 11],
        12,
        [11, None, 12, 13],
        [13, 13],  # duplicates
    ]
    file_path = make_distribution_file(tmp_path, rows)

    cl = CheckList(str(file_path))
    assert set(cl.species) == {10, 11, 12, 13}


def test_empty_file_yields_empty_species_list(tmp_path):
    file_path = make_distribution_file(tmp_path, [])
    # Pandas writes headers only; reading back gives empty dataframe
    cl = CheckList(str(file_path))
    assert cl.species == []
