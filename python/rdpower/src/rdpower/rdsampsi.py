#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import math
from rdrobust import rdrobust
from rdpower.rdpower_fun import rdpower_powerNR

def rdsampsi(data=None,
             cutoff=0,
             tau=None,
             alpha=0.05,
             beta=0.8,
             samph=None,
             nsamples=None,
             all=False,
             bias=None,
             variance=None,
             nratio=None,
             init_cond=None,
             plot=False,
             graph_range=None,
             covs=None,
             covs_drop=True,
             deriv=0,
             p=1,
             q=None,
             h=None,
             b=None,
             rho=None,
             kernel='triangular',
             bwselect='mserd',
             vce='nn',
             cluster=None,
             scalepar=1,
             scaleregul=1,
             fuzzy=None,
             level=95,
             weights=None,
             masspoints='adjust',
             bwcheck=None,
             bwrestrict=True,
             stdvars=False):
    
    """
    rdsampsi() performs sample size calculations for RD designs.
    
    Parameters:
        data: a matrix (Y, R) containing the outcome variable and the running variable (as column vectors).
        cutoff: the RD cutoff (default is 0).
        tau: specifies the treatment effect under the alternative at which the power function is evaluated. The default is half the standard deviation of the outcome for the untreated group.
        alpha: specifies the significance level for the power function. Default is 0.05.
        beta: specifies the desired power. Default is 0.8.
        nsamples: sets the total sample size to the left, sample size to the left inside the bandwidth, total sample size to the right, and sample size to the right of the cutoff inside the bandwidth to calculate the variance when the running variable is not specified. When not specified, the values are calculated using the running variable.
        samph: sets the bandwidths at each side of the cutoff for power calculation. The first number is the bandwidth to the left of the cutoff, and the second number is the bandwidth to the right. Default values are the bandwidths used by rdrobust.
        all: displays the power using the conventional variance estimator, in addition to the robust bias-corrected one.
        bias: set bias to the left and right of the cutoff. If not specified, the biases are estimated using rdrobust.
        variance: set variance to the left and right of the cutoff. If not specified, the variances are estimated using rdrobust.
        nratio: specifies the proportion of treated units in the window. Default is the ratio of the standard deviation of the treated to the sum of the standard deviations for treated and controls.
        init_cond: sets the initial condition for the Newton-Raphson algorithm that finds the sample size. Default is the number of observations in the sample with non-missing values of the outcome and running variable.
        plot: plots the power function using the conventional and robust bias-corrected standard errors from rdrobust.
        graph_range: range of the plot.
        covs: option for rdrobust(): specifies additional covariates to be used for estimation and inference.
        covs_drop: option for rdrobust(): if True, it checks for collinear additional covariates and drops them. Default is True.
        deriv: option for rdrobust(): specifies the order of the derivative of the regression functions to be estimated.
        p: option for rdrobust(): specifies the order of the local-polynomial used to construct the point-estimator.
        q: option for rdrobust(): specifies the order of the local-polynomial used to construct the bias-correction.
        h: option for rdrobust(): specifies the values of the main bandwidth to be used on the left and on the right of the cutoff, respectively.
        b: option for rdrobust(): specifies the values of the bias bandwidth b to be used on the left and on the right of the cutoff, respectively.
        rho: option for rdrobust(): specifies the value of rho so that the bias bandwidth b equals b = h / rho.
        kernel: option for rdrobust(): kernel function used to construct the local-polynomial estimators.
        bwselect: option for rdrobust(): specifies the bandwidth selection procedure to be used.
        vce: option for rdrobust(): specifies the procedure used to compute the variance-covariance matrix estimator.
        cluster: option for rdrobust(): indicates the cluster ID variable used for the cluster-robust variance estimation with degrees-of-freedom weights.
        scalepar: option for rdrobust(): specifies scaling factor for RD parameter of interest.
        scaleregul: option for rdrobust(): specifies scaling factor for the regularization terms of bandwidth selectors.
        fuzzy: option for rdrobust(): specifies the treatment status variable used to implement fuzzy RD estimation.
        level: option for rdrobust(): sets the confidence level for confidence intervals.
        weights: option for rdrobust(): is the variable used for optional weighting of the estimation procedure. The unit-specific weights multiply the kernel function.
        masspoints: option for rdrobust(): checks and controls for repeated observations in the running variable.
        bwcheck: option for rdrobust(): if a positive integer is provided, the preliminary bandwidth used in the calculations is enlarged so that at least bwcheck unique observations are used.
        bwrestrict: option for rdrobust(): if True, computed bandwidths are restricted to lie within the range of x. Default is True.
        stdvars: option for rdrobust(): if True, x and y are standardized before computing the bandwidths. Default is True.
    
    Returns:
        alpha: significance level
        beta: desired power
        tau: treatment effect under alternative hypothesis
        sampsi_h_tot: total number of observations inside the window
        sampsi_h_r: number of observations inside the window to the right of the cutoff
        sampsi_h_l: number of observations inside the window to the left of the cutoff
        N_r: Total sample size to the right of the cutoff
        N_l: Total sample size to the left of the cutoff
        samph_r: bandwidth to the right of the cutoff
        samph_l: bandwidth to the left of the cutoff
        var_r: Robust bias-corrected variance to the right of the cutoff
        Var_l: Robust bias-corrected variance to the left of the cutoff
        sampsi_h_tot_cl: implied total number of observations inside the window using conventional s.e.
        sampsi_h_r_cl: number of observations inside the window to the right of the cutoff using conventional s.e.
        sampsi_h_l_cl: number of observations inside the window to the left of the cutoff using conventional s.e.
        no_iter: number of iterations until convergence of the Newton-Raphson algorithm
        init_cond: initial condition of the Newton-Raphson algorithm

    Example:
        X = np.random.normal(size=(1000, 2))
        R = X[:, 0] + X[:, 1] + np.random.normal(size=1000)
        Y = 1 + R - 0.5 * R**2 + 0.3 * R**3 + (R >= 0) + np.random.normal(size=1000)

        # Sample size to achieve power of 0.8 against tau = 1
        tmp = rdsampsi(data=np.column_stack((Y, R)), tau=1)

        # Sample size against tau = 1 including covariates
        tmp = rdsampsi(data=np.column_stack((Y, R)), tau=1, covs=X)
    """
    
    #################################################################
    # Options, default values, and error checking
    #################################################################

    if data is not None:
        if data.shape[1] < 2:
            raise ValueError('Too few variables specified in data')
        else:
            Y = data.iloc[:, 0].values
            R = data.iloc[:, 1].values

    if nsamples is not None:
        if data is not None:
            print('nsamples ignored when data is specified')
        else:
            if len(nsamples) != 4:
                raise ValueError('Insufficient arguments in nsamples')
            else:
                if any([ns % 1 != 0 for ns in nsamples]):
                    raise ValueError('Sample sizes in nsamples have to be integers')
                if min(nsamples) <= 0:
                    raise ValueError('Sample sizes in nsamples have to be > 0')
                nminus = nsamples[0]
                n_hnew_l = nsamples[1]
                nplus = nsamples[2]
                n_hnew_r = nsamples[3]

    if samph is not None:
        if len(samph) == 1:
            hnew_l = samph
            hnew_r = samph
        elif len(samph) == 2:
            hnew_l = samph[0]
            hnew_r = samph[1]
        else:
            raise ValueError('samph incorrectly specified')

    if bias is not None:
        if len(bias) == 1:
            raise ValueError('Need to specify both Bl and Br')
        elif len(bias) == 2:
            bias_l = bias[0]
            bias_r = bias[1]
        else:
            raise ValueError('bias incorrectly specified')

    if variance is not None:
        if min(variance) <= 0:
            raise ValueError('Variances have to be > 0')
        if len(variance) == 1:
            raise ValueError('Need to specify both Vl and Vr')
        elif len(variance) == 2:
            vl = variance[0]
            vr = variance[1]
        else:
            raise ValueError('variance incorrectly specified')
        if data is None:
            vl_cl = vl
            vr_cl = vr

    if variance is not None and samph is None:
        raise ValueError('Need to set samph when variance is specified')

    if nratio is not None:
        nratio_cl = nratio
        if nratio >= 1 or nratio <= 0:
            raise ValueError('nratio has to be in (0, 1)')

    if q is None:
        q = p + 1
    
    #################################################################
    # Bias and variance
    #################################################################

    if data is not None:
        if bias is None or variance is None:
            aux = rdrobust(Y, R, c=cutoff, all=True, covs=covs, covs_drop=covs_drop, deriv=deriv, p=p, q=q, h=h, b=b,
                                    rho=rho, kernel=kernel, bwselect=bwselect, vce=vce, cluster=cluster,
                                    scalepar=scalepar, scaleregul=scaleregul, fuzzy=fuzzy, level=level,
                                    weights=weights, masspoints=masspoints, bwcheck=bwcheck, bwrestrict=bwrestrict,
                                    stdvars=stdvars)
            
            h_aux = aux.bws.values
            h_l = h_aux[0, 0]
            h_r = h_aux[0, 1]

            if bias is None:
                bias = aux.bias
                bias_l = bias[0] / (h_l ** (1 + p - deriv))
                bias_r = bias[1] / (h_r ** (1 + p - deriv))
            
            if variance is None:
                VL_CL = aux.V_cl_l
                VR_CL = aux.V_cl_r
                VL_RB = aux.V_rb_l
                VR_RB = aux.V_rb_r

                if cluster is not None:
                    N = len(np.unique(cluster[~np.isnan(Y) & ~np.isnan(R)]))
                else:
                    N = np.sum(~np.isnan(Y) & ~np.isnan(R))

                pos = deriv
                vl = N * (h_l ** (1 + 2 * deriv)) * VL_RB[pos, pos]
                vr = N * (h_r ** (1 + 2 * deriv)) * VR_RB[pos, pos]
                vl_cl = N * (h_l ** (1 + 2 * deriv)) * VL_CL[pos, pos]
                vr_cl = N * (h_r ** (1 + 2 * deriv)) * VR_CL[pos, pos]

        if samph is None:
            hnew_l = h_l
            hnew_r = h_r

        if vl_cl is None or vr_cl is None:
            vl_cl = vl
            vr_cl = vr

    # Bias adjustment
    bias = bias_r * hnew_r ** (1 + p - deriv) + bias_l * hnew_l ** (1 + p - deriv)

    # Variance adjustment
    V_rbc = vl / (hnew_l ** (1 + 2 * deriv)) + vr / (hnew_r ** (1 + 2 * deriv))
    stilde = np.sqrt(V_rbc)
    V_cl = vl_cl / (hnew_l ** (1 + 2 * deriv)) + vr_cl / (hnew_r ** (1 + 2 * deriv))
    stilde_cl = np.sqrt(V_cl)

    #################################################################
    # Sample size calculation
    #################################################################

    # Critical value

    z = norm.ppf(1-alpha/2)

    # Set default value of tau

    if tau is None:
        sd0 = np.nanstd(Y[(R >= cutoff - hnew_l) & (R < cutoff)])
        tau = 0.5 * sd0

     # Set initial value for Newton-Raphson

    if init_cond is None:
        init_cond = np.sum(~np.isnan(Y) & ~np.isnan(R))

    # Find m

    print('Calculating sample size...')

    maux = rdpower_powerNR(init_cond, tau, stilde, z, beta)
    m = maux['m']

    if all:
        maux1 = rdpower_powerNR(init_cond, tau + bias, stilde_cl, z, beta)
        m_cl = maux1['m']

    print('Sample size obtained.')

    if nratio is None:
        nratio = math.sqrt(vr) / (math.sqrt(vr) + math.sqrt(vl))
        nratio_cl = math.sqrt(vr_cl) / (math.sqrt(vr_cl) + math.sqrt(vl_cl))

    if data is not None:
        if cluster is not None:
            N = len(np.unique(cluster[~np.isnan(Y) & ~np.isnan(R)]))
            nplus = len(np.unique(cluster[(R >= cutoff) & (~np.isnan(Y)) & (~np.isnan(R))]))
            nminus = len(np.unique(cluster[(R < cutoff) & (~np.isnan(Y)) & (~np.isnan(R))]))
            n_hnew_r = len(np.unique(cluster[(R >= cutoff) & (R <= cutoff + hnew_r) & (~np.isnan(Y)) & (~np.isnan(R))]))
            n_hnew_l = len(np.unique(cluster[(R < cutoff) & (R >= cutoff - hnew_l) & (~np.isnan(Y)) & (~np.isnan(R))]))
        else:
            N = np.sum(~np.isnan(Y) & ~np.isnan(R))
            nplus = np.sum((R >= cutoff) & (~np.isnan(Y)) & (~np.isnan(R)))
            nminus = np.sum((R < cutoff) & (~np.isnan(Y)) & (~np.isnan(R)))
            n_hnew_r = np.sum((R >= cutoff) & (R <= cutoff + hnew_r) & (~np.isnan(Y)) & (~np.isnan(R)))
            n_hnew_l = np.sum((R < cutoff) & (R >= cutoff - hnew_l) & (~np.isnan(Y)) & (~np.isnan(R)))

    denom = nratio * nplus / n_hnew_r + (1 - nratio) * nminus / n_hnew_l
    denom_cl = nratio_cl * nplus / n_hnew_r + (1 - nratio_cl) * nminus / n_hnew_l

    M = m / denom
    Mr = np.ceil(M * nratio)
    Ml = np.ceil(M * (1 - nratio))
    M = Ml + Mr

    if all:
        M_cl = m_cl / denom_cl
        Mr_cl = np.ceil(M_cl * nratio_cl)
        Ml_cl = np.ceil(M_cl * (1 - nratio_cl))
        M_cl = Ml_cl + Mr_cl

    #################################################################
    # Descriptive statistics for display
    #################################################################

    if data is not None:

        # Left panel

        nplus_disp = np.sum((R >= cutoff) & (~np.isnan(Y)) & (~np.isnan(R)))
        nminus_disp = np.sum((R < cutoff) & (~np.isnan(Y)) & (~np.isnan(R)))

        if h is None:
            hl = hnew_l
            hr = hnew_r
        else:
            hl = h_l
            hr = h_r

        n_hnew_r_disp = np.sum((R >= cutoff) & (R <= cutoff + hr) & (~np.isnan(Y)) & (~np.isnan(R)))
        n_hnew_l_disp = np.sum((R < cutoff) & (R >= cutoff - hl) & (~np.isnan(Y)) & (~np.isnan(R)))

        if cluster is not None:
            gplus = len(np.unique(cluster[(R >= cutoff) & (~np.isnan(Y)) & (~np.isnan(R))]))
            gminus = len(np.unique(cluster[(R < cutoff) & (~np.isnan(Y)) & (~np.isnan(R))]))
            gplus_h_r = len(np.unique(cluster[(R >= cutoff) & (R <= cutoff + hr) & (~np.isnan(Y)) & (~np.isnan(R))]))
            gminus_h_l = len(np.unique(cluster[(R < cutoff) & (R >= cutoff - hl) & (~np.isnan(Y)) & (~np.isnan(R))]))

        # Right panel
        
        N_disp = np.sum(~np.isnan(Y) & ~np.isnan(R))

        if bias is None or variance is None:
            bwselect = aux.bwselect
            kernel_type = aux.kernel
            vce_type = aux.vce
        else:
            bwselect = None
            kernel_type = None
            vce_type = None
    else:

        # Left panel

        if cluster is not None:
            gplus = None
            gminus = None
        nplus_disp = None
        nminus_disp = None
        n_hnew_l_disp = None
        n_hnew_r_disp = None

        hl = None
        hr = None
        p = None

        # Right panel

        N_disp = None
        bwselect = None
        kernel_type = None
        vce_type = None

    # Size distortion

    if all:
        se_cl_aux = stilde_cl / np.sqrt(m_cl)
        size_dist = 1 - norm.cdf(bias / se_cl_aux + z) + norm.cdf(bias / se_cl_aux - z)

    #################################################################
    # Output
    #################################################################

    print('\n')
    print(f'Number of obs = {N_disp}')
    print(f'BW type       = {bwselect}')
    print(f'Kernel type   = {kernel_type}')
    print(f'VCE method    = {vce_type}')
    print(f'Derivative    = {deriv}')
    print(f'HA:       tau = {round(tau, 3)}')
    print(f'Power         = {round(beta, 3)}')
    if all:
        print(f'Size dist.    = {round(size_dist, 3)}')
    print('\n')

    print(f'{format("Cutoff c = " + str(round(cutoff, 3)), "22")} '
        f'{format("Left of c", "16")} '
        f'{format("Right of c", "16")}')
    
    if data is not None:
        print(f'{format("Number of obs", "22")} '
            f'{format(str(nminus_disp), "16")} '
            f'{format(str(nplus_disp), "16")}')
        print(f'{format("Eff. number of obs", "22")} '
            f'{format(str(n_hnew_l_disp), "16")} '
            f'{format(str(n_hnew_r_disp), "16")}')
        print(f'{format("BW loc. poly.", "22")} '
            f'{format(str(round(hl, 3)), "16")} '
            f'{format(str(round(hr, 3)), "16")}')
        print(f'{format("Order loc. poly.", "22")} '
            f'{format(str(p), "16")} '
            f'{format(str(p), "16")}')

    text_aux = "New sample"
    if cluster is not None:
        print(f'{format("Number of clusters", "22")} '
            f'{format(str(gminus), "16")} '
            f'{format(str(gplus), "16")}')
        print(f'{format("Eff. num. of clusters", "22")} '
            f'{format(str(gminus_h_l), "16")} '
            f'{format(str(gplus_h_r), "16")}')
        text_aux = "New cluster sample"
    print(f'{format("Sampling BW", "22")} '
        f'{format(str(round(hnew_l, 3)), "16")} '
        f'{format(str(round(hnew_r, 3)), "16")}')
    
    print('\n')

    print('='*89)
    if cluster is not None:
        print(f'{format("", "28")} '
            f'{format("Number of clusters in window", "35")} '
            f'{format("Proportion", "15")}')
    else:
        print(f'{format("", "33")} '
            f'{format("Number of obs in window", "35")} '
            f'{format("Proportion", "15")}')
    print(f'{format("", "25")} '
        f'{format("[c-h,c)", "15")} '
        f'{format("[c,c+h]", "15")} '
        f'{format("Total", "15")} '
        f'{format("[c,c+h]", "13")}')
    print('-'*89)
    print(f'{format("Robust bias-corrected", "25")} '
        f'{format(str(round(Ml,0)), "15")} '
        f'{format(str(Mr), "15")} '
        f'{format(str(M), "15")} '
        f'{format(str(round(nratio, 3)), "13")}')
    if all:
        print(f'{format("Conventional", "25")} '
            f'{format(str(round(Ml_cl,0)), "15")} '
            f'{format(str(Mr_cl), "15")} '
            f'{format(str(M_cl), "15")} '
            f'{format(str(round(nratio_cl, 3)), "13")}')
        print('='*89)
    else:
        print('='*89)

    #################################################################
    # Power function plot
    #################################################################

    if plot:
        N_plot = np.sum(~np.isnan(Y) & ~np.isnan(R))
        left = 0
        right = N_plot
        if graph_range is not None:
            left = graph_range[0]
            right = graph_range[1]

        def power_func(x):
            return 1 - norm.cdf(np.sqrt(x*denom)*tau/stilde + norm.ppf(1-alpha/2)) + norm.cdf(np.sqrt(x*denom)*tau/stilde - norm.ppf(1-alpha/2))

        plt.plot(np.arange(left, right+1), power_func(np.arange(left, right+1)))
        plt.xlabel('total sample size in window')
        plt.ylabel('power')
        plt.xlim(left, right)
        plt.title('Power function')
        plt.axvline(x=M, linestyle='dashed', color='gray')
        plt.axhline(y=beta, linestyle='dashed', color='gray')
        if all:
            plt.plot(np.arange(left, right+1), power_func(np.arange(left, right+1)), label='conventional')
            plt.legend(loc='lower left')
        plt.show()

    #################################################################
    # Return values
    #################################################################

    output = {
        'sampsi.h.tot': M,
        'sampsi.h.r': Mr,
        'sampsi.h.l': Ml,
        'sampsi.tot': m,
        'N.r': nplus,
        'N.l': nminus,
        'Nh.r': n_hnew_r,
        'Nh.l': n_hnew_l,
        'bias.r': bias_r,
        'bias.l': bias_r,
        'var.r': vr,
        'var.l': vl,
        'samph.r': hnew_r,
        'samph.l': hnew_l,
        'tau': tau,
        'beta': beta,
        'alpha': alpha,
        'init.cond': init_cond,
        'no.iter': maux['iter']
    }
    if all:
        output.update({
            'sampsi.h.tot.cl': M_cl,
            'sampsi.h.r.cl': Mr_cl,
            'sampsi.h.l.cl': Ml_cl,
            'sampsi.tot.cl': m_cl,
            'var.r.cl': vr_cl,
            'var.l.cl': vl_cl
        })

    return output