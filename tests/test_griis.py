# tests/test_griis.py
import pytest
import pandas as pd
import numpy as np
from unittest.mock import patch
from b3alien import griis, utils
from b3alien.griis import CheckList, filter_by_speciesprofile


# ===========================================================================
# Shared helpers
# ===========================================================================

def make_distribution_file(tmp_path, rows):
    """
    Write a tab-separated file with a 'speciesKey' column.
    We name it 'merged_distr.txt' because CheckList often appends this name.
    """
    def to_cell(v):
        if isinstance(v, list):
            return str(v)
        return v

    df = pd.DataFrame({"speciesKey": [to_cell(v) for v in rows]})
    # Create a directory and put the file inside to mimic expected structure
    data_dir = tmp_path / "data"
    data_dir.mkdir(exist_ok=True)
    path = data_dir / "merged_distr.txt"
    df.to_csv(path, sep="\t", index=False)
    return path


def make_species_profile(**kwargs):
    """Return a minimal speciesprofile DataFrame."""
    defaults = {
        "id": ["sp1", "sp2", "sp3", "sp4"],
        "isInvasive": ["Invasive", "Invasive", "Non-invasive", "Non-invasive"],
        "habitat": ["Terrestrial", "Marine", "Terrestrial", "Freshwater"],
    }
    defaults.update(kwargs)
    return pd.DataFrame(defaults)


def make_merged_distr(**kwargs):
    """Return a minimal merged_distr DataFrame (speciesKey already as list strings)."""
    defaults = {
        "id": ["sp1", "sp2", "sp3", "sp4"],
        "speciesKey": ["[100]", "[200]", "[300]", "[400]"],
    }
    defaults.update(kwargs)
    return pd.DataFrame(defaults)


def make_read_csv_side_effect(profile_df, merged_df):
    """
    Returns a side_effect for pd.read_csv that dispatches on filename,
    mimicking the two files CheckList reads when use_profile=True.
    """
    def side_effect(path, *args, **kwargs):
        path_str = str(path)
        if "speciesprofile" in path_str:
            return profile_df.copy()
        elif "merged_distr" in path_str:
            return merged_df.copy()
        # Fallback for simple tests that don't differentiate
        return merged_df.copy()
    return side_effect


# ===========================================================================
# CheckList — basic loading
# ===========================================================================

class TestCheckListBasic:

    @patch("pandas.read_csv")
    def test_griis_checklist(self, mock_read):
        # Mocking to avoid dependency on local file system paths
        mock_read.return_value = make_merged_distr()
        checklist = CheckList("fake/path/merged_distr.txt")
        assert hasattr(checklist, "species")
        assert isinstance(checklist.species, list)
        assert len(checklist.species) > 0

    def test_load_griis_parses_lists_ints_and_filters_uncertain_and_nulls(self, tmp_path):
        rows = [
            [123, 456],
            [456, "Uncertain", None],
            789,
            float("nan"),
        ]
        file_path = make_distribution_file(tmp_path, rows)
        cl = CheckList(str(file_path))
        assert set(cl.species) == {123, 456, 789}
        assert all(isinstance(x, int) for x in cl.species)

    def test_mixed_encoded_values_round_trip(self, tmp_path):
        rows = [
            [10, "Uncertain", 11],
            12,
            [11, None, 12, 13],
            [13, 13],
        ]
        file_path = make_distribution_file(tmp_path, rows)
        cl = CheckList(str(file_path))
        assert set(cl.species) == {10, 11, 12, 13}

    def test_empty_file_yields_empty_species_list(self, tmp_path):
        file_path = make_distribution_file(tmp_path, [])
        cl = CheckList(str(file_path))
        assert cl.species == []


# ===========================================================================
# CheckList — mutations
# ===========================================================================

class TestCheckListMutations:

    def test_add_species_deduplicates_and_preserves_existing(self, tmp_path):
        file_path = make_distribution_file(tmp_path, [[1, 2]])
        cl = CheckList(str(file_path))
        cl._add_species([2, 3, 4, 3])
        assert cl.species == [1, 2, 3, 4]

    def test_remove_species_filters_multiple_entries(self, tmp_path):
        file_path = make_distribution_file(tmp_path, [[1, 2, 3, 4]])
        cl = CheckList(str(file_path))
        cl._remove_species([2, 4, 999])  # 999 not present
        assert cl.species == [1, 3]


# ===========================================================================
# CheckList — use_profile parameter
# ===========================================================================

class TestCheckListUseProfile:

    @patch("b3alien.griis.filter_by_speciesprofile")
    @patch("pandas.read_csv")
    def test_use_profile_calls_filter_function(self, mock_read_csv, mock_filter):
        profile_df = make_species_profile()
        merged_df = make_merged_distr()
        mock_read_csv.side_effect = make_read_csv_side_effect(profile_df, merged_df)
        mock_filter.return_value = merged_df.copy()

        CheckList("/fake/path/", use_profile=True,
                  invasive_val="Invasive", habitat_val="Terrestrial")

        mock_filter.assert_called_once()

    @patch("b3alien.griis.filter_by_speciesprofile")
    @patch("pandas.read_csv")
    def test_use_profile_species_extracted_correctly(self, mock_read_csv, mock_filter):
        profile_df = make_species_profile()
        merged_df = make_merged_distr()
        mock_read_csv.side_effect = make_read_csv_side_effect(profile_df, merged_df)
        mock_filter.return_value = merged_df[merged_df["id"].isin(["sp1", "sp2"])].copy()

        cl = CheckList("/fake/path/", use_profile=True, invasive_val="Invasive")

        assert set(cl.species) == {100, 200}

    @patch("b3alien.griis.filter_by_speciesprofile")
    @patch("pandas.read_csv")
    def test_use_profile_false_skips_filter(self, mock_read_csv, mock_filter):
        merged_df = make_merged_distr()
        mock_read_csv.return_value = merged_df.copy()

        cl = CheckList("/fake/path/", use_profile=False)

        mock_filter.assert_not_called()
        assert set(cl.species) == {100, 200, 300, 400}

    @patch("b3alien.griis.filter_by_speciesprofile")
    @patch("pandas.read_csv")
    def test_uncertain_values_excluded_with_profile(self, mock_read_csv, mock_filter):
        profile_df = make_species_profile()
        merged_df = make_merged_distr(
            id=["sp1", "sp2"],
            speciesKey=["[100]", "['Uncertain']"],
        )
        mock_read_csv.side_effect = make_read_csv_side_effect(profile_df, merged_df)
        mock_filter.return_value = merged_df.copy()

        cl = CheckList("/fake/path/", use_profile=True, invasive_val="Invasive")

        assert "Uncertain" not in cl.species
        assert 100 in cl.species

    @patch("b3alien.griis.filter_by_speciesprofile")
    @patch("pandas.read_csv")
    def test_checklist_species_are_ints_with_profile(self, mock_read_csv, mock_filter):
        profile_df = make_species_profile()
        merged_df = make_merged_distr()
        mock_read_csv.side_effect = make_read_csv_side_effect(profile_df, merged_df)
        mock_filter.return_value = merged_df.copy()

        cl = CheckList("/fake/path/", use_profile=True)

        assert all(isinstance(s, int) for s in cl.species)

    @patch("b3alien.griis.filter_by_speciesprofile")
    @patch("pandas.read_csv")
    def test_init_params_stored_on_instance(self, mock_read_csv, mock_filter):
        profile_df = make_species_profile()
        merged_df = make_merged_distr()
        # Fix: Use side_effect here to ensure 'isInvasive' column exists when read
        mock_read_csv.side_effect = make_read_csv_side_effect(profile_df, merged_df)
        mock_filter.return_value = merged_df

        cl = CheckList("/fake/path/", use_profile=True,
                       invasive_val="Invasive", habitat_val="Marine")

        assert cl.use_profile is True
        assert cl.invasive_val == "Invasive"
        assert cl.habitat_val == "Marine"


# ===========================================================================
# filter_by_speciesprofile
# ===========================================================================

class TestFilterBySpeciesprofile:

    def setup_method(self):
        self.df_source = make_species_profile()
        self.df_target = make_merged_distr()

    def test_no_filters_returns_full_target(self):
        result = filter_by_speciesprofile(self.df_source, self.df_target)
        pd.testing.assert_frame_equal(result.reset_index(drop=True),
                                      self.df_target.reset_index(drop=True))

    def test_filter_invasive_only(self):
        result = filter_by_speciesprofile(
            self.df_source, self.df_target, invasive_val="Invasive")
        assert set(result["id"]) == {"sp1", "sp2"}

    def test_filter_invasive_and_terrestrial(self):
        result = filter_by_speciesprofile(
            self.df_source, self.df_target,
            invasive_val="Invasive", habitat_val="Terrestrial")
        assert set(result["id"]) == {"sp1"}


# ===========================================================================
# Utility / API functions
# ===========================================================================

class TestApiAndUtils:

    def test_matching(self):
        # This test relies on an external API (GBIF) usually
        # If this fails in CI, consider mocking get_speciesKey
        species = griis.get_speciesKey("Solanum nigrum")
        assert isinstance(species, list)

    def test_split_event_date(self):
        intro, outro = griis.split_event_date("2001/2020")
        assert intro == "2001"
        assert outro == "2020"

    def test_runtime_detection(self):
        assert isinstance(utils.in_ipython(), bool)