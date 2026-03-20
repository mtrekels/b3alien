import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import pytest
matplotlib.use("Agg")

# ---- import your module under test ----
# If your file is named differently, change this import:
from b3alien.simulation.simulation import *


# ---------- Basic numeric helpers ----------

def test_count_m_simple_trend():
    # m(t) = exp(m0 + m1*t)
    params = np.array([0.0, 0.1])  # only first two used
    t = np.array([0, 1, 2, 3])
    out = count_m(t, params)
    expected = np.exp(0.0 + 0.1 * t)
    np.testing.assert_allclose(out, expected, rtol=1e-12)


def test_count_pi_stability_and_range():
    # count_pi should be in [0,1], stable even for large exponents (via clipping)
    # Only elements 2,3,4 are used in count_pi
    params = np.array([0, 0, 5.0, 2.0, 3.0])  # (pi0, pi1, pi2) large-ish but clipped inside
    s = np.array([[1, 2], [3, 4]])
    t = np.array([[2, 3], [4, 5]])

    pi = count_pi(s, t, params)
    assert pi.shape == s.shape
    assert np.all(np.isfinite(pi))
    assert np.all((pi >= 0) & (pi <= 1))

    # Monotonic-ish sanity: if exponent is huge, pi -> 1
    big_params = np.array([0, 0, 1000.0, 1000.0, 1000.0])
    pi_big = count_pi(s, t, big_params)
    assert np.allclose(pi_big, 1.0, atol=1e-6)


def test_count_p_shape_and_bounds():
    # t is scalar in the function, output length should be t
    t = 5
    # Provide params with 5 entries for count_pi usage
    params = np.array([-1.0, 0.0, -0.5, 0.2, 0.1])
    p = count_p(t, params)
    assert p.shape == (t,)
    assert np.all(np.isfinite(p))
    assert np.all((p >= 0) & (p <= 1))


# ---------- Lambda & log-like ----------

def test_count_lambda_offset_and_nonnegativity():
    # count_lambda expects last element of params as offset (per your docstring)
    # First 5 drive m and p; last one is offset.
    params = np.array([-2.0, 0.0, -1.0, 0.1, 0.1, 0.5])  # offset=0.5
    N = 4
    lam = count_lambda(params, N)
    assert lam.shape == (N,)
    # offset applied and clamped non-negative
    assert np.all(lam >= 0.0)
    # crude check that offset visibly affects:
    params2 = params.copy(); params2[-1] = 1.5
    lam2 = count_lambda(params2, N)
    assert np.all(lam2 >= lam + 0.9)  # loose lower bound


def test_count_log_like_restrictions_and_shapes():
    # num_discov length defines N
    num_discov = pd.Series([2.0, 0.0, 1.0, 3.0])
    # params of length 6 since count_lambda uses the last as offset
    params = np.array([-1.0, 0.05, -1.2, 0.1, 0.1, 0.2])

    # Restrict vector: any value not 99 is a fixed value
    restrict = np.array([99, 99, -1.2, 0.1, 0.1, 99])  # fix indices 2,3,4
    ll, lambdas = count_log_like(params, restrict, num_discov)

    assert np.isfinite(ll)
    assert lambdas.shape == (len(num_discov),)

    # If we change a restricted param in "params", objective should NOT change
    params_changed = params.copy()
    params_changed[2] = -9.0  # restricted to -1.2 anyway
    ll2, lambdas2 = count_log_like(params_changed, restrict, num_discov)
    np.testing.assert_allclose(ll, ll2, rtol=1e-10)
    np.testing.assert_allclose(lambdas, lambdas2, rtol=1e-10)


# ---------- Simulators ----------

@pytest.fixture
def tiny_series():
    # small synthetic “discoveries per year”
    T = pd.Series([2000, 2001, 2002, 2003, 2004])
    y = pd.Series([0.5, 1.0, 1.5, 1.2, 0.8])
    return T, y


def test_simulate_solow_costello_runs_and_shapes(tiny_series):
    T, y = tiny_series
    C1, vec = simulate_solow_costello(T, y, vis=False)
    # vec is length 5 (Nelder–Mead via fmin with 5D guess), C1 length = len(y)
    assert len(C1) == len(y)
    assert len(vec) == 5
    assert np.all(np.isfinite(C1))
    assert np.all(C1 >= 0)  # from count_lambda clamp


def test_simulate_solow_costello_scipy_runs_and_shapes(tiny_series):
    T, y = tiny_series
    C1, vec = simulate_solow_costello_scipy(T, y, vis=False)
    # vec is length 6 (offset included), C1 length = len(y)
    assert len(C1) == len(y)
    assert len(vec) == 6
    assert np.all(np.isfinite(C1))
    assert np.all(C1 >= 0)


# ---------- Plotting ----------

def test_plot_with_confidence_no_gui(monkeypatch, tiny_series):
    T, y = tiny_series
    # Fabricate bootstrap results in expected shape
    n = len(T)
    fake_results = {
        "beta1_samples": np.array([0.01, 0.02, 0.03]),
        "beta1_ci": (0.01, 0.03),
        "c1_mean": np.linspace(0.5, 2.0, n),
        "c1_lower": np.linspace(0.4, 1.8, n),
        "c1_upper": np.linspace(0.6, 2.2, n),
        "c1_all": np.vstack([np.linspace(0.5, 2.0, n)]*3),
    }

    # Do not actually show the plot
    monkeypatch.setattr(plt, "show", lambda *a, **k: None)

    # Should not raise
    plot_with_confidence(T, y, fake_results)

# ------- test real simulation results from SC paper ------

def test_SC_original():
    inputf = 'tests/data/sc_original/NumDis.csv'
    n_discover = pd.read_csv(inputf, sep=';')
    numdis = np.array(n_discover["numdis"])
    T = np.array(n_discover["time"])

    C1, vec = simulate_solow_costello_scipy(T, numdis)

    assert vec[1] > 0.0134 and vec[1] < 0.0136


# This fake worker must be a top-level function to be pickleable by multiprocessing.
def _fake_bootstrap_worker_for_test(i, time_list, rate_list):
    """
    A predictable, fake bootstrap worker that returns known values for testing.
    """
    n_timesteps = len(time_list)
    # Return a fake parameter and a fake cumulative curve.
    # The parameter is just the iteration number `i`.
    # The curve is a simple linear array based on `i`.
    fake_param = float(i)
    fake_cumulative_curve = np.arange(n_timesteps) * (i + 1)
    return fake_param, fake_cumulative_curve


def test_run_bootstrap_analysis(monkeypatch, tiny_series):
    """
    Tests the aggregation and DataFrame creation logic of run_bootstrap_analysis
    by mocking the parallel worker.
    """
    T, y = tiny_series
    n_iterations = 10
    n_timesteps = len(T)

    # Replace the real worker with our fake one
    monkeypatch.setattr('b3alien.simulation.simulation.bootstrap_worker', _fake_bootstrap_worker_for_test)

    # --- Run the function under test ---
    results_df = run_bootstrap_analysis(list(T), list(y), n_iterations=n_iterations)

    # --- Assertions ---
    assert isinstance(results_df, pd.DataFrame)
    assert len(results_df) == n_timesteps
    assert 'Annual_Rate' in results_df.columns
    assert 'Cumulative_Value' in results_df.columns
    assert 'Lower_CI_95' in results_df.columns
    assert not results_df.isnull().values.any()


def test_get_bootstrap_errors(monkeypatch, tiny_series):
    """
    Tests the statistics calculation in get_bootstrap_errors by mocking the
    underlying simulation function.
    """
    T, y = tiny_series
    n_iterations = 5
    n_timesteps = len(T)

    # --- Mock the simulation function ---
    # This list will hold the predictable results our fake function will return.
    fake_results = []
    for i in range(n_iterations):
        # Create a unique but predictable result for each iteration
        fake_vec = np.array([0.1, 0.2, i, i, i, i])  # vec1_boot
        fake_c1 = np.full(n_timesteps, i + 1.0)      # C1_boot
        fake_results.append((fake_c1, fake_vec))

    # This fake function will pop one of the pre-made results on each call.
    def fake_simulate_scipy(annual_time, annual_rate, vis=False):
        return fake_results.pop(0)

    monkeypatch.setattr('b3alien.simulation.simulation.simulate_solow_costello_scipy', fake_simulate_scipy)

    # --- Run the function under test ---
    vec1_mean, vec1_std, C1_mean, C1_std = get_bootstrap_errors(T, y, iterations=n_iterations)

    # --- Assertions ---
    # Expected C1 values are [1.0, 2.0, 3.0, 4.0, 5.0] for each timestep
    assert C1_mean.shape == (n_timesteps,)
    np.testing.assert_allclose(C1_mean, np.full(n_timesteps, 3.0)) # Mean of 1,2,3,4,5 is 3
    np.testing.assert_allclose(C1_std, np.full(n_timesteps, np.std([1,2,3,4,5])))