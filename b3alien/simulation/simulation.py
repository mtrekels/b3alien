import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import fmin
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import warnings

def count_m(t, params):
    """Calculates the mean, mu, from Solow and Costello (2004)."""
    m0 = params[0]
    m1 = params[1]
    m = np.exp(m0 + m1 * t)
    return m

def count_pi(s, t, params):
    """Calculates the variable pi from Solow and Costello (2004)."""
    pi0 = params[2]
    pi1 = params[3]
    pi2 = params[4]
    exponent = np.clip(pi0 + pi1 * t + pi2 * np.exp(t - s), -700, 700)
    num = np.exp(exponent)
    pi = np.divide(num, (1 + num), out=np.zeros_like(num), where=(1 + num) != 0)
    pi = np.where(np.isinf(num), 1, pi)
    return pi

def count_p(t, params):
    """Calculates the value p from Solow and Costello (2004).
    It uses matrix coding for efficiency.
    """
    S = np.tile(np.arange(1, t + 1), (t, 1))
    thing = 1 - count_pi(S, S.T, params)
    thing[t - 1, :] = 1
    up = np.triu(np.ones_like(thing), 1)
    thing2 = np.tril(thing) + up
    product = np.prod(thing2, axis=0)
    pst = product * count_pi(np.arange(1, t + 1), t, params)
    return pst

def count_lambda(params, N):
    """
    This function calculates lambda from Solow and Costello, 2004.
    params is a vector of parameters
    N is the number of time points
    Note: an additional offset parameter is included to allow for a non-zero baseline.
    """
    lambda_result = np.zeros(N)
    for t in range(1, N + 1):
        S = np.arange(1, t + 1)
        Am = count_m(S, params)
        Ap = count_p(t, params)
        lambda_result[t - 1] = np.dot(Am, Ap)

    # apply the offset
    offset = params[-1]
    lambda_result += offset

    # Ensure no negative expected values
    lambda_result = np.clip(lambda_result, 0, None)

    return lambda_result

def count_log_like(params, restrict, num_discov):
    """
    Calculates the negative log likelihood from Solow and Costello (2004), 
    supporting optional parameter restrictions.
    """

    f = np.where(restrict != 99)[0]
    g = np.where(restrict == 99)[0]
    new_params = params.copy()
    new_params[g] = params[g]
    new_params[f] = restrict[f]

    # Use count_lambda (which now includes the offset)
    lambda_values = count_lambda(new_params, len(num_discov))

    # Compute log-likelihood components safely
    summand2 = np.where(
        lambda_values > 0,
        num_discov * np.log(lambda_values) - lambda_values,
        -lambda_values
    )

    LL = -np.sum(summand2)
    return LL, lambda_values


def simulate_solow_costello(annual_time_gbif, annual_rate_gbif, vis=False): 
    """
        Solow-Costello simulation of the rate of establishment.

        Parameters
        ----------
        annual_time_gbif : pandas.Series
            Time series of the rate of establishment.
        annual_rate_gbif : pandas.Series
            Rates corresponding to the time series.
        vis : bool, optional
            Create a plot of the simulation. Default is False.
            
        Returns
        -------
        C1: numpy.Series
            Result of the simulation.
        val1: numpy.Series
            Parameters of the fitting.
    """

    #  global num_discov;  #  No need for global, pass as argument
    num_discov = pd.Series(annual_rate_gbif).T   #  Load and transpose
    T = pd.Series(annual_time_gbif) #np.arange(1851, 1996)  #  Create the time period
    #  options = optimset('TolFun',.01,'TolX',.01);  #  Tolerance is handled differently in scipy

    guess = np.array([-1.1106, 0.0135, -1.4534, 0.1, 0.1])  #  Initial guess
    constr = 99 * np.ones_like(guess)  #  Constraint vector

    vec1 = fmin(
        lambda x: count_log_like(x, constr, num_discov)[0],
        guess,
        xtol=0.01,
        ftol=0.01,
        disp=0  # disables all output
    )
    
    val1 = count_log_like(vec1, constr, num_discov)[0]  #  Get the function value at the minimum


    C1 = count_lambda(vec1, len(num_discov))  #  Calculate the mean of Y

    if vis:
        #  Create the plot
        plt.plot(T, np.cumsum(num_discov), 'k-', T, np.cumsum(C1), 'k--')
        plt.legend(['Discoveries', 'Unrestricted'])
        plt.xlabel('Time')
        plt.ylabel('Cumulative Discovery')
        plt.show()

    return C1, vec1

def simulate_solow_costello_scipy(annual_time_gbif, annual_rate_gbif, vis=False): 
    """
        Solow-Costello simulation of the rate of establishment. Uses scipy's minimize for optimization.

        Parameters
        ----------
        annual_time_gbif : pandas.Series
            Time series of the rate of establishment.
        annual_rate_gbif : pandas.Series
            Rates corresponding to the time series.
        vis : bool, optional
            Create a plot of the simulation. Default is False.
            
        Returns
        -------
        C1: numpy.Series
            Result of the simulation.
        val1: numpy.Series
            Parameters of the fitting.
    """

    #  global num_discov;  #  No need for global, pass as argument
    num_discov = pd.Series(annual_rate_gbif).T   #  Load and transpose
    T = pd.Series(annual_time_gbif) #np.arange(1851, 1996)  #  Create the time period
    #  options = optimset('TolFun',.01,'TolX',.01);  #  Tolerance is handled differently in scipy
    guess = np.array([-1.1106, 0.0135, -1.4534, 0.1, 0.1, 0.0])  #  Initial guess
    constr = 99 * np.ones_like(guess) 

    # Objective function for minimize (returns scalar log-likelihood)
    def objective(x):
        return count_log_like(x, constr, num_discov)[0]  # still log-likelihood

    # Define bounds for each parameter
    # These must match the size and meaning of `guess`
    bounds = [
        (-5, 0),     # e.g., parameter 1: negative decay
        (0, 1),      # e.g., parameter 2: rate between 0 and 1
        (-5, 0),     # e.g., parameter 3: another decay
        (0.01, 2),   # e.g., parameter 4: noise scale
        (0.01, 2),   # e.g., parameter 5: another scale
        (-1, 1),  # adding an additional offset
    ]

    # Run bounded optimization
    result = minimize(
        objective,
        guess,
        method="Nelder-Mead",    
        options={"xatol": 0.01, "fatol": 0.01, "disp": False, "maxiter": 1000}
    )

    vec1 = result.x


    C1 = count_lambda(vec1, len(num_discov))  #  Calculate the mean of Y

    if vis:
        #  Create the plot
        plt.plot(T, np.cumsum(num_discov), 'k-', T, np.cumsum(C1), 'k--')
        plt.legend(['Discoveries', 'Unrestricted'])
        plt.xlabel('Time')
        plt.ylabel('Cumulative Discovery')
        plt.show()

    return C1, vec1

def bootstrap_worker(i, time_list, rate_list):
    '''
    Bootstrap on the residuals
    Returns=
    - fitting parameters (vec1)
    - C1_sim cumulative prediction from refit
    '''
    time_series = pd.Series(time_list)
    rate_series = pd.Series(rate_list)

    try:
        # Fit once to get baseline model
         with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)

            # Fit once to get baseline model
            C1_fit, vec1 = simulate_solow_costello_scipy(time_series, rate_series)
            residuals = rate_series.reset_index(drop=True) - C1_fit

            # Bootstrap residuals and create new synthetic data
            resampled_residuals = residuals.sample(frac=1, replace=True).reset_index(drop=True)
            simulated_rate = C1_fit + resampled_residuals

            # Fit again on simulated data
            C1_sim, vec1_boot = simulate_solow_costello_scipy(time_series, simulated_rate)
            
            return vec1_boot[1], np.cumsum(C1_sim)
    
    except Exception:
        return None

def parallel_bootstrap_solow_costello(annual_time_gbif, annual_rate_gbif, n_iterations=1000, ci=95):
    """
        Perform parallel bootstrapping of the Solow-Costello model 
        to estimate confidence intervals.
        Parameters
        ----------
        annual_time_gbif : pandas.Series
            Time series of the rate of establishment.
        annual_rate_gbif : pandas.Series
            Rates corresponding to the time series.
        n_iterations : int, optional
            Number of bootstrap iterations. Default is 1000.
        ci : float, optional
            Confidence interval percentage. Default is 95.
        Returns
        -------
        dict
            A dictionary containing bootstrap results and confidence intervals.
    """
    time_list = list(annual_time_gbif)
    rate_list = list(annual_rate_gbif)
    n_cores = max(1, multiprocessing.cpu_count() - 1)

    beta1_samples = []
    c1_curves = []

    with ProcessPoolExecutor(max_workers=n_cores) as executor:
        futures = [
            executor.submit(bootstrap_worker, i, time_list, rate_list)
            for i in range(n_iterations)
        ]

        for f in tqdm(as_completed(futures), total=n_iterations, desc="Bootstrapping"):
            try:
                result = f.result()
                if result is not None:
                    beta1, c1_curve = result
                    beta1_samples.append(beta1)
                    c1_curves.append(c1_curve)
            except Exception as e:
                print(f"Unhandled error in future: {e}")

    if not beta1_samples:
        raise RuntimeError("All bootstrap iterations failed. No valid samples.")

    beta1_samples = np.array(beta1_samples)
    c1_curves = np.array(c1_curves)  # shape: (B, T)

    ci_lower = np.percentile(beta1_samples, (100 - ci) / 2)
    ci_upper = np.percentile(beta1_samples, 100 - (100 - ci) / 2)
    ci_beta1 = (ci_lower, ci_upper)

    lower_band = np.percentile(c1_curves, (100 - ci) / 2, axis=0)
    upper_band = np.percentile(c1_curves, 100 - (100 - ci) / 2, axis=0)
    mean_cumsum = np.mean(c1_curves, axis=0)

    return {
        "beta1_samples": beta1_samples,
        "beta1_ci": ci_beta1,
        "c1_mean": mean_cumsum,
        "c1_lower": lower_band,
        "c1_upper": upper_band,
        "c1_all": c1_curves
    }

def plot_with_confidence(T, observed, results):
    """
        Plot the observed cumulative discoveries 
        with bootstrap confidence intervals.

        Parameters
        ----------
        T : pandas.Series
            Time series of the rate of establishment.
        observed : pandas.Series
            Observed cumulative discoveries.
        results : dict
            Dictionary containing bootstrap results and confidence intervals.
    """
    plt.figure(figsize=(10, 6))
    plt.plot(T, np.cumsum(observed), 'k-', label='Observed Discoveries')
    plt.plot(T, results["c1_mean"], 'b--', label='Bootstrap Mean C1')
    plt.fill_between(T, results["c1_lower"], results["c1_upper"], color='blue', alpha=0.3, label='95% CI')
    plt.xlabel('Time')
    plt.ylabel('Cumulative Discoveries')
    plt.title(f'Solow-Costello Fit with CI (β₁ 95% CI: {results["beta1_ci"][0]:.4f} – {results["beta1_ci"][1]:.4f})')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()