# tests/test_griis.py
import pytest
from b3alien import griis
from b3alien import utils

def test_griis_checklist():
    checklist = griis.CheckList("tests/data/dwca-griis-portugal-v1.9/merged_distr.txt")
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
    assert utils.in_script() is False  or utils.in_script is True # Just ensure it runs without error