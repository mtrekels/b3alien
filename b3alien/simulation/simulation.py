import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import fmin

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
    """

    #  global num_discov;  #  No need for global, pass as argument
    num_discov = pd.Series(annual_rate_gbif).T   #  Load and transpose
    T = pd.Series(annual_time_gbif) #np.arange(1851, 1996)  #  Create the time period
    #  options = optimset('TolFun',.01,'TolX',.01);  #  Tolerance is handled differently in scipy

    guess = np.array([-1.1106, 0.0135, -1.4534, 0.1, 0.1])  #  Initial guess
    constr = 99 * np.ones_like(guess)  #  Constraint vector

    vec1 = fmin(lambda x: count_log_like(x, constr, num_discov)[0], guess, xtol=0.01, ftol=0.01)
    val1 = count_log_like(vec1, constr, num_discov)[0]  #  Get the function value at the minimum


    C1 = count_lambda(vec1, len(num_discov))  #  Calculate the mean of Y

    if vis:
        #  Create the plot
        plt.plot(T, np.cumsum(num_discov), 'k-', T, np.cumsum(C1), 'k--')
        plt.legend(['Discoveries', 'Unrestricted'])
        plt.xlabel('Time')
        plt.ylabel('Cumulative Discovery')
        plt.show()

    return C1