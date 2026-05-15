#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.stats import norm

_EPS = np.finfo(float).eps
_MIN_DERIVATIVE = 1e-5

#################################################################
# Power function
#################################################################

def rdpower_powerfun(n, tau, stilde, z):
    x = 1 - norm.cdf(np.sqrt(n)*tau/stilde + z) + norm.cdf(np.sqrt(n)*tau/stilde - z)
    return x

#################################################################
# Power function derivative w.r.t. n
#################################################################

def rdpower_powerfun_dot(n, tau, stilde, z):
    x = (norm.pdf(np.sqrt(n)*tau/stilde - z) - norm.pdf(np.sqrt(n)*tau/stilde + z)) * tau / (2 * stilde * np.sqrt(n))
    return x

#################################################################
# Power function derivative w.r.t. tau
#################################################################

def rdpower_powerfun_dot_tau(n, tau, stilde, z):
    x = (norm.pdf(np.sqrt(n)*tau/stilde - z) - norm.pdf(np.sqrt(n)*tau/stilde + z)) * np.sqrt(n) / stilde
    return x

#################################################################
# Newton-Raphson for sample size calculation
#################################################################

def rdpower_powerNR(x0, tau, stilde, z, beta):
    tol = 1
    iter = 0

    while tol > _EPS:
        iter += 1
        k = 1

        # Check if derivative at x0 is too small
        dot0 = rdpower_powerfun_dot(x0, tau, stilde, z)
        power0 = rdpower_powerfun(x0, tau, stilde, z)
        while dot0 < _MIN_DERIVATIVE:
            x0 = 1.2 * x0 * (power0 <= beta) + 0.8 * x0 * (power0 > beta)
            iter += 1
            dot0 = rdpower_powerfun_dot(x0, tau, stilde, z)
            power0 = rdpower_powerfun(x0, tau, stilde, z)

        x1 = x0 - (power0 - beta) / dot0

        # Check if x1 is negative or too small
        while x1 < 2:
            x1 = x0 - k * (power0 - beta) / dot0
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

    while tol > _EPS:
        iter += 1

        # Check if derivative at x0 is too small
        dot0 = rdpower_powerfun_dot_tau(n, tau0, stilde, z)
        power0 = rdpower_powerfun(n, tau0, stilde, z)
        while dot0 < _MIN_DERIVATIVE:
            tau0 = 1.5 * tau0 * (power0 <= beta) + 0.5 * tau0 * (power0 > beta)
            iter += 1
            dot0 = rdpower_powerfun_dot_tau(n, tau0, stilde, z)
            power0 = rdpower_powerfun(n, tau0, stilde, z)

        tau1 = tau0 - (power0 - beta) / dot0

        tol = abs(rdpower_powerfun(n, tau1, stilde, z) - beta)
        tau0 = tau1

    b = rdpower_powerfun(n, tau1, stilde, z)

    output = {
        'mde': tau1,
        'iter': iter,
        'powercheck': b
    }

    return output
