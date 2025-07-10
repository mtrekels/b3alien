import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import fmin
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

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
    """
    lambda_result = np.zeros(N)
    for t in range(1, N + 1):
        S = np.arange(1, t + 1)
        Am = count_m(S, params)
        Ap = count_p(t, params)
        lambda_result[t - 1] = np.dot(Am, Ap)
    return lambda_result

def count_log_like(params, restrict, num_discov):
    """
    This function file calculates the log likelihood function for Solow and
    Costello (2004).  It takes into account any possible restrictions (See
    below)

    params is a vector of parameters
    restrict is a vector (same size as params) that places restrictions on the
    parameters. If restrict[i]=99, then there is no restriction for the ith
    parameter. If restrict[i]=0 (for example) then the restriction is exactly
    that.
    """

    f = np.where(restrict != 99)[0]
    g = np.where(restrict == 99)[0]
    new_params = params.copy()
    new_params[g] = params[g]
    new_params[f] = restrict[f]

    summand2 = np.zeros_like(num_discov, dtype=float)
    lambda_values = np.zeros_like(num_discov, dtype=float)

    for t in range(1, len(num_discov) + 1):
        S = np.arange(1, t + 1)
        Am = count_m(S, new_params)
        Ap = count_p(t, new_params)
        lambda_t = np.dot(Am, Ap)
        lambda_values[t - 1] = lambda_t
        summand2[t - 1] = num_discov[t - 1] * np.log(lambda_t) - lambda_t if lambda_t > 0 else -lambda_t

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
    ]

    # Run bounded optimization
    result = minimize(
        objective,
        guess,
        method="L-BFGS-B",     # supports bounds
        bounds=bounds,
        options={"ftol": 0.01, "gtol": 0.01, "disp": False}
    )

    vec1 = result.x
    val1 = result.fun


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
    '''
    time_series = pd.Series(time_list)
    rate_series = pd.Series(rate_list)

    # Fit once to get baseline model
    C1_fit, vec1 = simulate_solow_costello_scipy(time_series, rate_series)
    residuals = rate_series.reset_index(drop=True) - C1_fit

    # Bootstrap residuals and create new synthetic data
    resampled_residuals = residuals.sample(frac=1, replace=True).reset_index(drop=True)
    simulated_rate = C1_fit + resampled_residuals

    # Fit again on simulated data
    C1_sim, _ = simulate_solow_costello_scipy(time_series, simulated_rate)
    return np.cumsum(C1_sim)

def parallel_bootstrap_solow_costello(annual_time_gbif, annual_rate_gbif, n_iterations=1000, ci=95):
    time_list = list(annual_time_gbif)
    rate_list = list(annual_rate_gbif)
    n_cores = max(1, multiprocessing.cpu_count() - 1)

    C1_samples = []

    with ProcessPoolExecutor(max_workers=n_cores) as executor:
        futures = [
            executor.submit(bootstrap_worker, i, time_list, rate_list)
            for i in range(n_iterations)
        ]

        for f in tqdm(as_completed(futures), total=n_iterations, desc="Bootstrapping"):
            try:
                C1_samples.append(f.result())
            except Exception as e:
                print(f"Error in worker: {e}")

    C1_samples = np.array(C1_samples)
    lower_bound = np.percentile(C1_samples, (100 - ci) / 2, axis=0)
    upper_bound = np.percentile(C1_samples, 100 - (100 - ci) / 2, axis=0)
    mean_cumsum = np.mean(C1_samples, axis=0)

    return mean_cumsum, lower_bound, upper_bound

def plot_with_confidence(T, observed, mean_cumsum, lower_bound, upper_bound):
    plt.figure(figsize=(10, 6))
    plt.plot(T, np.cumsum(observed), 'k-', label='Observed Discoveries')
    plt.plot(T, mean_cumsum, 'b--', label='Model Prediction')
    plt.fill_between(T, lower_bound, upper_bound, color='blue', alpha=0.3, label='95% CI')
    plt.xlabel('Time')
    plt.ylabel('Cumulative Discoveries')
    plt.legend()
    plt.title('Solow-Costello Model with 95% Confidence Interval')
    plt.show()