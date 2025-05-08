# tests/test_cube.py
import pytest
import numpy as np
import pandas as pd
from b3alien import b3cube
from b3alien import griis
from b3alien import simulation

def test_checklist():

    checklist = griis.CheckList("tests/data/dwca-griis-portugal-v1.9/merged_distr.txt")

    assert len(checklist.species) > 0

def test_solow_costello():

    cube = b3cube.OccurrenceCube("tests/data/data_PT-30.parquet")
    checklist = griis.CheckList("tests/data/dwca-griis-portugal-v1.9/merged_distr.txt")

    df_sparse, df_cumulative = b3cube.cumulative_species(cube, checklist.species)

    annual_time_gbif, annual_rate_gbif = b3cube.calculate_rate(df_cumulative)

    c1 = simulation.simulate_solow_costello(annual_time_gbif, annual_rate_gbif, vis=False)
    
    assert len(c1) > 0
