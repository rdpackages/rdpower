#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
from rdrobust import rdrobust

def rdpower(data=None, cutoff=0, tau=None, alpha=0.05, nsamples=None, sampsi=None, samph=None, all=False,
            bias=None, variance=None, plot=False, graph_range=None, covs=None, covs_drop=True, deriv=0,
            p=1, q=None, h=None, b=None, rho=None, kernel='triangular', bwselect='mserd', vce='nn',
            cluster=None, scalepar=1, scaleregul=1, fuzzy=None, level=95, weights=None, masspoints='adjust',
            bwcheck=None, bwrestrict=True, stdvars=False):
    
    """
    Power Calculations for RD Designs

    rdpower performs power calculations for RD designs.

    Author:
    Matias Cattaneo, Princeton University. Email: cattaneo@princeton.edu
    Rocio Titiunik, Princeton University. Email: titiunik@princeton.edu
    Gonzalo Vazquez-Bare, UC Santa Barbara. Email: gvazquez@econ.ucsb.edu

    References:
    Cattaneo, M. D., R. Titiunik, and G. Vazquez-Bare. (2019).
    Power Calculations for Regression Discontinuity Designs. Stata Journal, 19(1): 210-245.
    URL: https://rdpackages.github.io/references/Cattaneo-Titiunik-VazquezBare_2019_Stata.pdf

    Parameters:
    - data: A matrix (Y, R) containing the outcome variable and the running variable (as column vectors).
    - cutoff: The RD cutoff (default is 0).
    - tau: Specifies the treatment effect under the alternative at which the power function is evaluated. Default is half the standard deviation of the outcome for the untreated group.
    - alpha: Specifies the significance level for the power function. Default is 0.05.
    - nsamples: Sets the total sample size to the left, sample size to the left inside the bandwidth, total sample size to the right, and sample size to the right of the cutoff inside the bandwidth to calculate the variance when the running variable is not specified. When not specified, the values are calculated using the running variable.
    - sampsi: Sets the sample size at each side of the cutoff for power calculation. The first number is the sample size to the left of the cutoff, and the second number is the sample size to the right. Default values are the sample sizes inside the chosen bandwidth.
    - samph: Sets the bandwidths at each side of the cutoff for power calculation. The first number is the bandwidth to the left of the cutoff, and the second number is the bandwidth to the right. Default values are the bandwidths used by rdrobust.
    - all: Displays the power using the conventional variance estimator, in addition to the robust bias-corrected one.
    - bias: Sets bias to the left and right of the cutoff. If not specified, the biases are estimated using rdrobust.
    - variance: Sets variance to the left and right of the cutoff. If not specified, the variances are estimated using rdrobust.
    - plot: Plots the power function using the conventional and robust bias-corrected standard errors from rdrobust.
    - graph.range: Range of the plot.
    - covs: Option for rdrobust: specifies additional covariates to be used for estimation and inference.
    - covs_drop: Option for rdrobust: if True, it checks for collinear additional covariates and drops them. Default is True.
    - deriv: Option for rdrobust: specifies the order of the derivative of the regression functions to be estimated.
    - p: Option for rdrobust: specifies the order of the local polynomial used to construct the point estimator.
    - q: Option for rdrobust: specifies the order of the local polynomial used to construct the bias-correction.
    - h: Option for rdrobust: specifies the values of the main bandwidth to be used on the left and on the right of the cutoff, respectively.
    - b: Option for rdrobust: specifies the values of the bias bandwidth b to be used on the left and on the right of the cutoff, respectively.
    - rho: Option for rdrobust: specifies the value of rho so that the bias bandwidth b equals b=h/rho.
    - kernel: Option for rdrobust: kernel function used to construct the local-polynomial estimators.
    - bwselect: Option for rdrobust: specifies the bandwidth selection procedure to be used.
    - vce: Option for rdrobust: specifies the procedure used to compute the variance-covariance matrix estimator.
    - cluster: Option for rdrobust: indicates the cluster ID variable used for the cluster-robust variance estimation with degrees-of-freedom weights.
    - scalepar: Option for rdrobust: specifies scaling factor for RD parameter of interest.
    - scaleregul: Option for rdrobust: specifies scaling factor for the regularization terms of bandwidth selectors.
    - fuzzy: Option for rdrobust: specifies the treatment status variable used to implement fuzzy RD estimation.
    - level: Option for rdrobust: sets the confidence level for confidence intervals.
    - weights: Option for rdrobust: is the variable used for optional weighting of the estimation procedure. The unit-specific weights multiply the kernel function.
    - masspoints: Option for rdrobust: checks and controls for repeated observations in the running variable.
    - bwcheck: Option for rdrobust: if a positive integer is provided, the preliminary bandwidth used in the calculations is enlarged so that at least bwcheck unique observations are used.
    - bwrestrict: Option for rdrobust: if True, computed bandwidths are restricted to lie within the range of x. Default is bwrestrict=True.
    - stdvars: Option for rdrobust: if True, x and y are standardized before computing the bandwidths. Default is stdvars=True.

    Returns:
    - power_rbc: Power against tau using robust bias-corrected standard error.
    - se_rbc: Robust bias-corrected standard error.
    - sampsi_r: Number of observations inside the window to the right of the cutoff.
    - sampsi_l: Number of observations inside the window to the left of the cutoff.
    - samph_r: Bandwidth to the right of the cutoff.
    - samph_l: Bandwidth to the left of the cutoff.
    - alpha: Significance level used in power function.
    - tau: Treatment effect under the alternative hypothesis.
    - bias_r: Bias to the right of the cutoff.
    - bias_l: Bias to the left of the cutoff.
    - Vr_rb: Robust bias-corrected variance to the right of the cutoff.
    - Vl_rb: Robust bias-corrected variance to the left of the cutoff.
    - N_r: Total sample size to the right of the cutoff.
    - N_l: Total sample size to the left of the cutoff.
    - power_conv: Power against tau using conventional inference.
    - se_conv: Conventional standard error.

    Examples:
    # Toy dataset
    import numpy as np

    X = np.random.randn(1000, 2)
    R = X[:, 0] + X[:, 1] + np.random.randn(1000)
    Y = 1 + R - 0.5 * R**2 + 0.3 * R**3 + (R >= 0) + np.random.randn(1000)

    # Power against tau = 1
    tmp = rdpower(data=np.column_stack((Y, R)), tau=1)

    # Power against tau = 1 including covariates
    tmp = rdpower(data=np.column_stack((Y, R)), tau=1, covs=X)
    """
    
    #################################################################
    # Options, default values, and error checking
    #################################################################

    if data is not None:
        if data.shape[1] < 2:
            raise ValueError('Too few variables specified in data')
        elif data.shape[1] > 2:
            raise ValueError('Too many variables specified in data')
        else:
            Y = data.iloc[:, 0].values
            R = data.iloc[:, 1].values

    if nsamples is not None:
        if data is not None:
            print('nsamples ignored when data is specified')
        else:
            if bias is None or variance is None or samph is None or sampsi is None or tau is None:
                raise ValueError('Not enough information to calculate power without data')
            if len(nsamples) != 4:
                raise ValueError('Incorrect number of arguments in nsamples')
            else:
                if any(n % 1 != 0 for n in nsamples):
                    raise ValueError('Sample sizes in nsamples have to be integers')
                if any(n <= 0 for n in nsamples):
                    raise ValueError('Sample sizes in nsamples have to be >0')
                nminus, n_hnew_l, nplus, n_hnew_r = nsamples
    
    if sampsi is not None:
        sampsi = np.array(sampsi)
        if np.sum(sampsi % 1) != 0:
            raise ValueError('Sample sizes in sampsi have to be integers')
        if np.min(sampsi) <= 0:
            raise ValueError('Sample sizes in sampsi have to be > 0')
        if np.isscalar(sampsi):
            ntilde_l = sampsi
            ntilde_r = sampsi
        elif len(sampsi) == 2:
            ntilde_l, ntilde_r = sampsi
        else:
            raise ValueError('sampsi incorrectly specified')

    if samph is not None:
        samph = np.array(samph)
        if len(samph) == 1:
            hnew_l = samph
            hnew_r = samph
        elif len(samph) == 2:
            hnew_l = samph[0]
            hnew_r = samph[1]
        else:
            raise ValueError('samph incorrectly specified')

    if bias is not None:
        bias = np.array(bias)
        if np.isscalar(bias):
            raise ValueError('Need to specify both Bl and Br')
        elif len(bias) == 2:
            bias_l, bias_r = bias
        else:
            raise ValueError('bias incorrectly specified')

    if variance is not None:
        variance = np.array(variance)
        if min(variance) <= 0:
            raise ValueError('Variances have to be > 0')
        if np.isscalar(variance):
            raise ValueError('Need to specify both Vl and Vr')
        elif len(variance) == 2:
            Vl_rb, Vr_rb = variance
        else:
            raise ValueError('variance incorrectly specified')

    if q is None: q = p + 1

    #################################################################
    # Definition of bias, variance, sample sizes and bandwidths
    #################################################################

    if data is not None:
        if bias is None or variance is None:
            aux = rdrobust(Y, R, c=cutoff, all=True, covs=covs, covs_drop=covs_drop, deriv=deriv, p=p, q=q, h=h, b=b, rho=rho, cluster=cluster,
                            kernel=kernel, bwselect=bwselect, vce=vce, scalepar=scalepar, scaleregul=scaleregul, fuzzy=fuzzy, level=level,
                            weights=weights, masspoints=masspoints, bwcheck=bwcheck, bwrestrict=bwrestrict, stdvars=stdvars)
            
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
                    N = np.unique(cluster[~np.isnan(Y) & ~np.isnan(R)]).shape[0]
                else:
                    N = np.sum(~np.isnan(Y) & ~np.isnan(R))
                
                pos = deriv
                Vl_cl = N * (h_l ** (1 + 2 * deriv)) * VL_CL[pos, pos]
                Vr_cl = N * (h_r ** (1 + 2 * deriv)) * VR_CL[pos, pos]
                Vl_rb = N * (h_l ** (1 + 2 * deriv)) * VL_RB[pos, pos]
                Vr_rb = N * (h_r ** (1 + 2 * deriv)) * VR_RB[pos, pos]
        
        # Set default new bandwidth
        if samph is None:
            hnew_l = h_l
            hnew_r = h_r
        
        # Set default value of tau
        if tau is None:
            sd0 = np.nanstd(Y[(R >= cutoff - hnew_l) & (R < cutoff)],ddof =0)
            tau = 0.5 * sd0
        
        # Calculate sample sizes
        if cluster is not None:
            N = np.unique(cluster[~np.isnan(Y) & ~np.isnan(R)]).shape[0]
            nplus = np.unique(cluster[(R >= cutoff) & (~np.isnan(Y)) & (~np.isnan(R))]).shape[0]
            nminus = np.unique(cluster[(R < cutoff) & (~np.isnan(Y)) & (~np.isnan(R))]).shape[0]
            n_hnew_r = np.unique(cluster[(R >= cutoff) & (R <= cutoff + hnew_r) & (~np.isnan(Y)) & (~np.isnan(R))]).shape[0]
            n_hnew_l = np.unique(cluster[(R < cutoff) & (R >= cutoff - hnew_l) & (~np.isnan(Y)) & (~np.isnan(R))]).shape[0]
        else:
            N = np.sum(~np.isnan(Y) & ~np.isnan(R))
            nplus = np.sum((R >= cutoff) & (~np.isnan(Y)) & (~np.isnan(R)))
            nminus = np.sum((R < cutoff) & (~np.isnan(Y)) & (~np.isnan(R)))
            n_hnew_r = np.sum((R >= cutoff) & (R <= cutoff + hnew_r) & (~np.isnan(Y)) & (~np.isnan(R)))
            n_hnew_l = np.sum((R < cutoff) & (R >= cutoff - hnew_l) & (~np.isnan(Y)) & (~np.isnan(R)))
        
        if sampsi is None:
            ntilde_l = n_hnew_l
            ntilde_r = n_hnew_r

    ntilde = nplus * (ntilde_r / n_hnew_r) + nminus * (ntilde_l / n_hnew_l)

    # Variance and bias adjustment
    V_rbc = Vl_rb / (ntilde * (hnew_l ** (1 + 2 * deriv))) + Vr_rb / (ntilde * (hnew_r ** (1 + 2 * deriv)))
    se_rbc = np.sqrt(V_rbc)

    if all:
        V_conv = Vl_cl / (ntilde * (hnew_l ** (1 + 2 * deriv))) + Vr_cl / (ntilde * (hnew_r ** (1 + 2 * deriv)))
        se_conv = np.sqrt(V_conv)
    else:
        V_conv = V_rbc
        se_conv = se_rbc

    bias = bias_r * (hnew_r ** (1 + p - deriv)) + bias_l * (hnew_l ** (1 + p - deriv))

    # Power calculation
    power_rbc = 1 - norm.cdf((tau) / se_rbc + norm.ppf(1 - alpha / 2)) + norm.cdf((tau) / se_rbc - norm.ppf(1 - alpha / 2))
    power_conv = 1 - norm.cdf((tau + bias) / se_conv + norm.ppf(1 - alpha / 2)) + norm.cdf((tau + bias) / se_conv - norm.ppf(1 - alpha / 2))

    power_rbc_list = []
    power_conv_list = []

    for r in [0, 2, 5, 8]:
        te = tau * r / 10

        power_rbc_aux = 1 - norm.cdf((te) / se_rbc + norm.ppf(1 - alpha / 2)) + norm.cdf((te) / se_rbc - norm.ppf(1 - alpha / 2))
        power_conv_aux = 1 - norm.cdf((te + bias) / se_conv + norm.ppf(1 - alpha / 2)) + norm.cdf((te + bias) / se_conv - norm.ppf(1 - alpha / 2))

        power_rbc_list.append(power_rbc_aux)
        power_conv_list.append(power_conv_aux)

    #################################################################
    # Descriptive statistics for display
    #################################################################
    
    if data is not None:

        # Left panel

        nplus_disp = sum((R >= cutoff) & (~np.isnan(Y)) & (~np.isnan(R)))
        nminus_disp = sum((R < cutoff) & (~np.isnan(Y)) & (~np.isnan(R)))

        if h is None:
            hl = hnew_l
            hr = hnew_r
        else:
            hl = h_l
            hr = h_r

        nhr_disp = sum((R >= cutoff) & (R <= cutoff + hr) & (~np.isnan(Y)) & (~np.isnan(R)))
        nhl_disp = sum((R < cutoff) & (R >= cutoff - hl) & (~np.isnan(Y)) & (~np.isnan(R)))

        if cluster is not None:
            gplus = len(set(cluster[(R >= cutoff) & (~np.isnan(Y)) & (~np.isnan(R))]))
            gminus = len(set(cluster[(R < cutoff) & (~np.isnan(Y)) & (~np.isnan(R))]))
            gplus_h_r = len(set(cluster[(R >= cutoff) & (R <= cutoff + hr) & (~np.isnan(Y)) & (~np.isnan(R))]))
            gminus_h_l = len(set(cluster[(R < cutoff) & (R >= cutoff - hl) & (~np.isnan(Y)) & (~np.isnan(R))]))

        # Right panel

        N_disp = sum((~np.isnan(Y)) & (~np.isnan(R)))

        if bias is not None and variance is not None:
            bwselect = None
            kernel_type = None
            vce_type = None
        else:
            bwselect = aux.bwselect
            kernel_type = aux.kernel
            vce_type = aux.vce

    else:

        # Left panel

        if cluster is not None:
            gplus = None
            gminus = None
        nplus_disp = None
        nminus_disp = None
        nhr_disp = None
        nhl_disp = None

        hl = None
        hr = None
        p = None

        # Right panel

        N_disp = None
        bwselect = None
        kernel_type = None
        vce_type = None


    #################################################################
    # Output
    #################################################################

    print('\n')
    print(f"Number of obs  =  {N_disp}")
    print(f"BW type        =  {bwselect}")
    print(f"Kernel type    =  {kernel_type}")
    print(f"VCE method     =  {vce_type}")
    print(f"Derivative     =  {deriv}")
    print(f"HA:       tau  =  {round(tau, 3)}")

    if all == True:
        print(f"Size dist.     =  {round(power_conv_list[3] - alpha, 3)}")

    print('\n')

    print(f"{'Cutoff c = ':10}{round(cutoff, 3):<10}{'Left of c':>16}{'Right of c':>16}")
    if data is not None:
        print(f"{'Number of obs':20}{nminus_disp:16}{nplus_disp:16}")
        print(f"{'Eff. # of obs':20}{nhl_disp:16}{nhr_disp:16}")
        print(f"{'BW loc. poly.':20}{round(hl, 3):16}{round(hr, 3):16}")
        print(f"{'Order loc. poly.':20}{p:16}{p:16}")

    text_aux = "New sample"
    if cluster is not None:
        print(f"{'Number of clusters':20}{gminus:16}{gplus:16}")
        print(f"{'Eff. # of clusters':20}{gminus_h_l:16}{gplus_h_r:16}")
        text_aux = "New cluster sample"

    print(f"{'Sampling BW':20}{round(hnew_l, 3):16}{round(hnew_r, 3):16}")
    print(f"{text_aux:20}{ntilde_l:16}{ntilde_r:16}")

    print("\n")

    print('='*100)
    print(f"{'Power against:':25}{'H0: tau = ':>15}{'0.2*tau = ':>15}{'0.5*tau = ':>15}{'0.8*tau = ':>15}{'tau = ':>15}")
    print(f"{'              ':22}{round(0, 3):15}{round(0.2*tau, 3):15}{round(0.5*tau, 3):15}{round(0.8*tau, 3):15}{round(tau, 3):15}")
    print('-'*100)
    print(f"{'Robust bias-corrected':22}{round(power_rbc_list[0], 3):15}{round(power_rbc_list[1], 3):15}{round(power_rbc_list[2], 3):15}{round(power_rbc_list[3], 3):15}{round(power_rbc, 3):15}")

    if all == True:
        print(f"{'Conventional':22}{round(power_conv_list[0], 3):15}{round(power_conv_list[1], 3):15}{round(power_conv_list[2], 3):15}{round(power_conv_list[3], 3):15}{round(power_conv, 3):15}")

    print('='*100)

    #################################################################
    # Power function plot
    #################################################################

    if plot:
        if graph_range is None:
            left = round(-1.5 * tau)
            right = round(1.5 * tau)
        else:
            left = graph_range[0]
            right = graph_range[1]

        def power_func(x):
            return 1 - norm.cdf((x) / se_rbc + norm.ppf(1 - alpha / 2)) + norm.cdf((x) / se_rbc - norm.ppf(1 - alpha / 2))

        plt.plot(np.linspace(left, right, num=100), [power_func(x) for x in np.linspace(left, right, num=100)])
        plt.xlabel('tau')
        plt.ylabel('power')
        plt.xlim(left, right)
        plt.title('Power function')
        plt.axvline(0, linestyle='dashed', color='gray')
        plt.axhline(alpha, linestyle='dashed', color='gray')
        if all:
            plt.plot(np.linspace(left, right, num=100), [(power_func(x + bias) for x in np.linspace(left, right, num=100))], linestyle='dashed')
            plt.legend(['robust bc', 'conventional'], loc='lower left')

        plt.show()


    #################################################################
    # Return values
    #################################################################

    output = {
        'power_rbc': power_rbc,
        'se_rbc': se_rbc,
        'sampsi_r': ntilde_r,
        'sampsi_l': ntilde_l,
        'samph_r': hnew_r,
        'samph_l': hnew_l,
        'N_r': nplus,
        'N_l': nminus,
        'Nh_l': n_hnew_l,
        'Nh_r': n_hnew_r,
        'tau': tau,
        'bias_r': bias_r,
        'bias_l': bias_l,
        'Vr_rb': Vr_rb,
        'Vl_rb': Vl_rb,
        'alpha': alpha
    }

    if all:
        output.update({
            'power_conv': power_conv,
            'se_conv': se_conv
        })

    return output