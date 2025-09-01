# tests/test_simulation.py
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

def test_simulate_solow_costello_basic():
    # Example input data
    annual_time = np.array([2000, 2001, 2002, 2003])
    annual_rate = np.array([5, 10, 15, 20])
    C1, vec1 = simulation.simulate_solow_costello(annual_time, annual_rate, vis=False)
    assert isinstance(C1, np.ndarray) or isinstance(C1, list)
    assert len(C1) == len(annual_time)

def test_simulate_solow_costello_empty():
    annual_time = np.array([])
    annual_rate = np.array([])
    result, _ = simulation.simulate_solow_costello_scipy(annual_time, annual_rate, vis=False)
    assert len(result) == 0

def test_simulate_scipy():
    annual_time = np.array([2000, 2001, 2002])
    annual_rate = np.array([2, 4, 6])
    c1, vec1 = simulation.simulate_solow_costello_scipy(annual_time, annual_rate)
    assert len(c1) == len(annual_time)

def test_simulate_bootstrap():
    # If your module raises exceptions for bad input
    annual_time = np.array([2000, 2001, 2002])
    annual_rate = np.array([2, 4, 6])
    results = simulation.parallel_bootstrap_solow_costello(annual_time, annual_rate, n_iterations=10)
    assert len(results) == 6