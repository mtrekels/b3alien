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