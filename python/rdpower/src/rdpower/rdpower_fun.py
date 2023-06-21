#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import numpy as np
import scipy.stats as stats

#################################################################
# Power function
#################################################################

def rdpower_powerfun(n, tau, stilde, z):
    x = 1 - stats.norm.cdf(np.sqrt(n)*tau/stilde + z) + stats.norm.cdf(np.sqrt(n)*tau/stilde - z)
    return x

#################################################################
# Power function derivative w.r.t. n
#################################################################

def rdpower_powerfun_dot(n, tau, stilde, z):
    x = (stats.norm.pdf(np.sqrt(n)*tau/stilde - z) - stats.norm.pdf(np.sqrt(n)*tau/stilde + z)) * tau / (2 * stilde * np.sqrt(n))
    return x

#################################################################
# Power function derivative w.r.t. tau
#################################################################

def rdpower_powerfun_dot_tau(n, tau, stilde, z):
    x = (stats.norm.pdf(np.sqrt(n)*tau/stilde - z) - stats.norm.pdf(np.sqrt(n)*tau/stilde + z)) * np.sqrt(n) / stilde
    return x

#################################################################
# Newton-Raphson for sample size calculation
#################################################################

def rdpower_powerNR(x0, tau, stilde, z, beta):
    tol = 1
    iter = 0

    while tol > sys.float_info.epsilon:
        iter += 1
        k = 1

        # Check if derivative at x0 is too small
        while rdpower_powerfun_dot(x0, tau, stilde, z) < 0.00001:
            x0 = 1.2 * x0 * (rdpower_powerfun(x0, tau, stilde, z) <= beta) + 0.8 * x0 * (rdpower_powerfun(x0, tau, stilde, z) > beta)
            iter += 1

        x1 = x0 - (rdpower_powerfun(x0, tau, stilde, z) - beta) / rdpower_powerfun_dot(x0, tau, stilde, z)

        # Check if x1 is negative or too small
        while x1 < 2:
            x1 = x0 - k * (rdpower_powerfun(x0, tau, stilde, z) - beta) / rdpower_powerfun_dot(x0, tau, stilde, z)
            k /= 2
            iter += 1

        tol = abs(rdpower_powerfun(x1, tau, stilde, z) - beta)
        x0 = x1

    b = rdpower_powerfun(x1, tau, stilde, z)
    m = np.ceil(x1)

    output = {
        'm': m,
        'iter': iter,
        'powercheck': b
    }

    return output

#################################################################
# Newton-Raphson for MDE calculation
#################################################################

def rdpower_powerNR_mde(n, tau0, stilde, z, beta):
    tol = 1
    iter = 0

    while tol > sys.float_info.epsilon:
        iter += 1

        # Check if derivative at x0 is too small
        while rdpower_powerfun_dot_tau(n, tau0, stilde, z) < 0.00001:
            tau0 = 1.5 * tau0 * (rdpower_powerfun(n, tau0, stilde, z) <= beta) + 0.5 * tau0 * (rdpower_powerfun(n, tau0, stilde, z) > beta)
            iter += 1

        tau1 = tau0 - (rdpower_powerfun(n, tau0, stilde, z) - beta) / rdpower_powerfun_dot_tau(n, tau0, stilde, z)

        tol = abs(rdpower_powerfun(n, tau1, stilde, z) - beta)
        tau0 = tau1

    b = rdpower_powerfun(n, tau1, stilde, z)

    output = {
        'mde': tau1,
        'iter': iter,
        'powercheck': b
    }

    return output