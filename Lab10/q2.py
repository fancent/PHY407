"""
This file aims to 
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def model(x,a,b):
    return a*x + b
    
def get_tau_step():
    """
    Calculate how far a photon travels before it gets scattered.
    Input : tau - optical depth of the atmosphere
    Output: optical depth traveled
    """
    delta_tau = -np.log(np.random.random())
    return delta_tau

def emit_photon(tau_max):
    """
    Emit a photon from the stellar core.
    Input: tau max - max optical depth
    Output: tau: optical depth at which the photon is created
            mu: directional cosine of the photon emitted
    """
    tau = tau_max
    delta_tau = get_tau_step()
    mu = np.random.random()
    return tau-delta_tau*mu, mu

def scatter_photon(tau):
    """
    Scatter a photon.
    Input : tau âĹŠ optical depth of the atmosphere
    Output: tau: new optical depth
            mu: directional cosine of the photon scattered
    """
    delta_tau = get_tau_step()
    # sample mu uniformly from -1 to 1
    mu = 2*np.random.random()-1
    tau = tau + delta_tau*mu
    return tau, mu

def photon_path(tau_max):
    """
    This function tracks the path of photon exiting the core.
    uncomment to get the mu and number of scattering.
    """
    tau, mu = emit_photon(tau_max)
    numOfScattered = 0
    while tau >= 0:
        tau, mu = scatter_photon(tau)
        if tau > tau_max:
            tau, mu = emit_photon(tau_max)
            numOfScattered = 0
        numOfScattered += 1
    # print("mu = {} scattered {} times".format(mu, numOfScattered))
    return mu

def simulation(N, tau_max):
    """
    This function simulates N amount of photons and calculate 
    the normalised frequency of each mu range.
    Then it checks if it matches the proposed equation 1.
    """
    mu_results = []
    intensity_results = []
    for _ in range(N):
        mu = photon_path(tau_max)
        mu_results.append(mu)
    mu_results = np.abs(np.array(mu_results))
    mu_avg_results = []
    count_results = []
    hist_results, bins, patches = plt.hist(mu_results, bins=20)
    for i in range(len(hist_results)):
        mu_avg = (bins[i]+bins[i+1])/2
        normalised_count = hist_results[i]/hist_results[-1]
        intensity = normalised_count/mu_avg
        mu_avg_results.append(mu_avg)
        intensity_results.append(intensity)
        count_results.append(normalised_count)
    mu_avg_results = np.array(mu_avg_results)

    # plotting histogram
    fig1, graph1 = plt.subplots()
    graph1.bar(mu_avg_results, count_results,width=0.05)
    graph1.set(xlabel='mu range', ylabel='frequency',
        title="Histogram")
    fig1.savefig("q2b_bar_tau={}.png".format(tau_max))
    plt.show()

    popt, pcov = curve_fit(model, mu_avg_results, intensity_results)
    print("proposed parameters:", popt)
    fig2, graph2 = plt.subplots()
    graph2.plot(mu_avg_results, intensity_results, label="our results")
    graph2.plot(mu_avg_results, model(mu_avg_results, popt[0], popt[1]), label="curve fit line")
    graph2.set(xlabel='mu range', ylabel='intensity',
        title="best fit")
    graph2.legend()
    fig2.savefig("q2b_graph_tau={}.png".format(tau_max))
    plt.show()
    
# simulation(10**5, 10)
simulation(10**5, 0.0001)
    

