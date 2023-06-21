#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.stats import norm
from rdrobust import rdrobust
from rdpower.rdpower_fun import rdpower_powerNR_mde

def rdmde(data=None,
          cutoff=0,
          alpha=0.05,
          beta=0.8,
          nsamples=None,
          sampsi=None,
          samph=None,
          all=False,
          bias=None,
          variance=None,
          init_cond=None,
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
    MDE Calculations for RD Designs
    
    Perform MDE calculations for RD designs.
    
    Parameters:
        data: numpy array or pandas DataFrame
            A matrix (Y, R) containing the outcome variable and the running variable (as column vectors).
        cutoff: float, optional (default=0)
            The RD cutoff.
        alpha: float, optional (default=0.05)
            The significance level for the power function.
        beta: float, optional (default=0.8)
            The desired power.
        nsamples: tuple, optional
            Sets the total sample size to the left, sample size to the left inside the bandwidth, total sample size to the right, 
            and sample size to the right of the cutoff inside the bandwidth to calculate the variance when the running variable is not specified.
            When not specified, the values are calculated using the running variable.
        sampsi: tuple, optional
            Sets the sample size at each side of the cutoff for power calculation. The first number is the sample size to the left of the cutoff,
            and the second number is the sample size to the right. Default values are the sample sizes inside the chosen bandwidth.
        samph: tuple, optional
            Sets the bandwidths at each side of the cutoff for power calculation. The first number is the bandwidth to the left of the cutoff,
            and the second number is the bandwidth to the right. Default values are the bandwidths used by rdrobust.
        all: bool, optional (default=False)
            Displays the power using the conventional variance estimator, in addition to the robust bias-corrected one.
        bias: tuple, optional
            Sets bias to the left and right of the cutoff. If not specified, the biases are estimated using rdrobust.
        variance: tuple, optional
            Sets variance to the left and right of the cutoff. If not specified, the variances are estimated using rdrobust.
        init_cond: float, optional
            Sets the initial condition for the Newton-Raphson algorithm that finds the MDE. Default is 0.2 times the standard deviation of the outcome below the cutoff.
        covs: numpy array or pandas DataFrame, optional
            Specifies additional covariates to be used for estimation and inference.
        covs_drop: bool, optional (default=True)
            If True, checks for collinear additional covariates and drops them.
        deriv: int, optional
            Specifies the order of the derivative of the regression functions to be estimated.
        p: int, optional
            Specifies the order of the local-polynomial used to construct the point estimator.
        q: int, optional
            Specifies the order of the local-polynomial used to construct the bias-correction.
        h: tuple, optional
            Specifies the values of the main bandwidth to be used on the left and on the right of the cutoff, respectively.
        b: tuple, optional
            Specifies the values of the bias bandwidth b to be used on the left and on the right of the cutoff, respectively.
        rho: float, optional
            Specifies the value of rho so that the bias bandwidth b equals b=h/rho.
        kernel: str, optional
            Kernel function used to construct the local-polynomial estimators.
        bwselect: str, optional
            Bandwidth selection procedure to be used.
        vce: str, optional
            Procedure used to compute the variance-covariance matrix estimator.
        cluster: str, optional
            Indicates the cluster ID variable used for the cluster-robust variance estimation with degrees-of-freedom weights.
        scalepar: float, optional
            Scaling factor for RD parameter of interest.
        scaleregul: float, optional
            Scaling factor for the regularization terms of bandwidth selectors.
        fuzzy: str, optional
            Treatment status variable used to implement fuzzy RD estimation.
        level: float, optional
            Sets the confidence level for confidence intervals.
        weights: numpy array or pandas Series, optional
            The variable used for optional weighting of the estimation procedure. The unit-specific weights multiply the kernel function.
        masspoints: bool, optional
            Checks and controls for repeated observations in the running variable.
        bwcheck: int, optional
            If a positive integer is provided, the preliminary bandwidth used in the calculations is enlarged so that at least bwcheck unique observations are used.
        bwrestrict: bool, optional (default=True)
            If True, computed bandwidths are restricted to lie within the range of the running variable.
        stdvars: bool, optional (default=True)
            If True, x and y are standardized before computing the bandwidths.
    
    Returns:
        mde: float
            MDE using robust bias-corrected standard error.
        se_rbc: float
            Robust bias-corrected standard error.
        sampsi_r: int
            Number of observations inside the window to the right of the cutoff.
        sampsi_l: int
            Number of observations inside the window to the left of the cutoff.
        samph_r: float
            Bandwidth to the right of the cutoff.
        samph_l: float
            Bandwidth to the left of the cutoff.
        alpha: float
            Significance level used in the power function.
        bias_r: float
            Bias to the right of the cutoff.
        bias_l: float
            Bias to the left of the cutoff.
        Vr_rb: float
            Robust bias-corrected variance to the right of the cutoff.
        Vl_rb: float
            Robust bias-corrected variance to the left of the cutoff.
        N_r: int
            Total sample size to the right of the cutoff.
        N_l: int
            Total sample size to the left of the cutoff.
        mde_conv: float
            MDE using conventional inference.
        se_conv: float
            Conventional standard error.

    Example:
        # Toy dataset
        X = np.random.randn(2000).reshape((1000, 2))
        R = X[:, 0] + X[:, 1] + np.random.randn(1000)
        Y = 1 + R - 0.5 * R ** 2 + 0.3 * R ** 3 + (R >= 0) + np.random.randn(1000)

        # MDE calculation
        tmp = rdmde(np.column_stack((Y, R)), init_cond=0.5)
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
            if bias is None or variance is None or samph is None or init_cond is None:
                raise ValueError('Not enough information to calculate power without data')
            if len(nsamples) != 4:
                raise ValueError('Incorrect number of arguments in nsamples')
            else:
                if sum(np.mod(nsamples, 1)) != 0:
                    raise ValueError('Sample sizes in nsamples have to be integers')
                if min(nsamples) <= 0:
                    raise ValueError('Sample sizes in nsamples have to be >0')
                nminus = nsamples[0]
                n_hnew_l = nsamples[1]
                nplus = nsamples[2]
                n_hnew_r = nsamples[3]
    
    if sampsi is not None:
        sampsi = np.array(sampsi)
        if np.sum(sampsi % 1) != 0:


            raise ValueError('Sample sizes in sampsi have to be integers')
        if np.min(sampsi) <= 0:
            raise ValueError('Sample sizes in sampsi have to be >0')
        if len(sampsi) == 1:
            ntilde_l = sampsi
            ntilde_r = sampsi
        elif len(sampsi) == 2:
            ntilde_l = sampsi[0]
            ntilde_r = sampsi[1]
        else:
            raise ValueError('sampsi incorrectly specified')
    
    if samph is not None:
        if np.isscalar(samph):
            hnew_l = samph
            hnew_r = samph
        elif len(samph) == 2:
            hnew_l, hnew_r = samph
        else:
            raise ValueError('samph incorrectly specified')
    
    if bias is not None:
        if np.isscalar(bias):
            raise ValueError('Need to specify both Bl and Br')
        elif len(bias) == 2:
            bias_l, bias_r  = bias
        else:
            raise ValueError('bias incorrectly specified')
    
    if variance is not None:
        if np.min(variance) <= 0:
            raise ValueError('Variances have to be >0')
        if np.isscalar(variance):
            raise ValueError('Need to specify both Vl and Vr')
        elif len(variance) == 2:
            Vl_rb, Vr_rb = variance
        else:
            raise ValueError('variance incorrectly specified')
    
    if q is None:
        q = p + 1
    
    #################################################################
    # Definition of bias, variance, sample sizes, and bandwidths
    #################################################################
    
    if data is not None:
        if bias is None or variance is None:
            aux = rdrobust(Y, R, c=cutoff, all=True, covs=covs, covs_drop=covs_drop, deriv=deriv, p=p, q=q, h=h, b=b,
                           rho=rho, cluster=cluster, kernel=kernel, bwselect=bwselect, vce=vce, scalepar=scalepar,
                           scaleregul=scaleregul, fuzzy=fuzzy, level=level, weights=weights, masspoints=masspoints,
                           bwcheck=bwcheck, bwrestrict=bwrestrict, stdvars=stdvars)
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
                Vl_cl = N * (h_l ** (1 + 2 * deriv)) * VL_CL[pos, pos]
                Vr_cl = N * (h_r ** (1 + 2 * deriv)) * VR_CL[pos, pos]
                Vl_rb = N * (h_l ** (1 + 2 * deriv)) * VL_RB[pos, pos]
                Vr_rb = N * (h_r ** (1 + 2 * deriv)) * VR_RB[pos, pos]
    
        ## Set default new bandwidth
    
        if samph is None:
            hnew_l = h_l
            hnew_r = h_r
    
        ## Calculate sample sizes
    
        if cluster is not None:
            N = len(np.unique(cluster[(~np.isnan(Y)) & (~np.isnan(R))]))
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
    
    if sampsi is None:
        ntilde_l = n_hnew_l
        ntilde_r = n_hnew_r
    ntilde = nplus * (ntilde_r / n_hnew_r) + nminus * (ntilde_l / n_hnew_l)

    #################################################################
    # Variance and bias adjustment
    #################################################################

    ## Variance adjustment

    V_rbc = Vl_rb / (hnew_l ** (1 + 2 * deriv)) + Vr_rb / (hnew_r ** (1 + 2 * deriv))
    se_rbc = np.sqrt(V_rbc)

    if all:
        V_conv = Vl_cl / (hnew_l ** (1 + 2 * deriv)) + Vr_cl / (hnew_r ** (1 + 2 * deriv))
        se_conv = np.sqrt(V_conv)
    else:
        V_conv = V_rbc
        se_conv = se_rbc

    ## Bias adjustment

    bias = bias_r * hnew_r ** (1 + p - deriv) + bias_l * hnew_l ** (1 + p - deriv)

    #################################################################
    # MDE calculation
    #################################################################

    ## Set initial condition for Newton-Raphson

    if init_cond is None:
        sd0 = np.nanstd(Y[R < cutoff])
        tau0 = 0.2 * sd0
    else:
        tau0 = init_cond

    print('Calculating MDE...', end="")

    mde_aux = rdpower_powerNR_mde(ntilde, tau0, se_rbc, norm.ppf(1 - alpha / 2), beta)
    mde = mde_aux['mde']

    beta_list = np.zeros(4)
    mde_rbc_list = np.zeros(4)

    count = 0
    for r in [-0.125, -0.0625, 0.0625, 0.125]:
        baux = beta * (1 + r)
        beta_list[count] = baux

        if baux < 1:
            rdpow_aux = rdpower_powerNR_mde(ntilde, tau0, se_rbc, norm.ppf(1 - alpha / 2), baux)
            mde_rbc_list[count] = rdpow_aux['mde']

        count += 1

    if all:
        mde_conv_aux = rdpower_powerNR_mde(ntilde, tau0 + bias, se_conv, norm.ppf(1 - alpha / 2), beta)
        mde_conv = mde_conv_aux['mde']

        mde_conv_list = np.zeros(4)
        count = 0
        for r in [-0.125, -0.0625, 0.0625, 0.125]:
            baux = beta * (1 + r)

            if baux < 1:
                rdpow_conv_aux = rdpower_powerNR_mde(ntilde, tau0 + bias, se_conv, norm.ppf(1 - alpha / 2), baux)
                mde_conv_list[count] = rdpow_conv_aux['mde']

            count += 1

    print('MDE obtained')

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

        nhr_disp = np.sum((R >= cutoff) & (R <= cutoff + hr) & (~np.isnan(Y)) & (~np.isnan(R)))
        nhl_disp = np.sum((R < cutoff) & (R >= cutoff - hl) & (~np.isnan(Y)) & (~np.isnan(R)))

        if cluster is not None:
            gplus = len(np.unique(cluster[(R >= cutoff) & (~np.isnan(Y)) & (~np.isnan(R))]))
            gminus = len(np.unique(cluster[(R < cutoff) & (~np.isnan(Y)) & (~np.isnan(R))]))
            gplus_h_r = len(np.unique(cluster[(R >= cutoff) & (R <= cutoff + hr) & (~np.isnan(Y)) & (~np.isnan(R))]))
            gminus_h_l = len(np.unique(cluster[(R < cutoff) & (R >= cutoff - hl) & (~np.isnan(Y)) & (~np.isnan(R))]))

        # Right panel

        N_disp = np.sum(~np.isnan(Y) & ~np.isnan(R))

        if bias is not None and variance is not None:
            bwselect = np.nan
            kernel_type = np.nan
            vce_type = np.nan
        else:
            bwselect = aux.bwselect
            kernel_type = aux.kernel
            vce_type = aux.vce

    else:

        # Left panel

        if cluster is not None:
            gplus = np.nan
            gminus = np.nan
        nplus_disp = np.nan
        nminus_disp = np.nan
        nhr_disp = np.nan
        nhl_disp = np.nan

        hl = np.nan
        hr = np.nan
        p = np.nan

        # Right panel

        N_disp = np.nan
        bwselect = np.nan
        kernel_type = np.nan
        vce_type = np.nan

    #################################################################
    # Output
    #################################################################

    print('')
    print(f'Number of obs = {N_disp}')
    print(f'BW type       = {bwselect}')
    print(f'Kernel type   = {kernel_type}')
    print(f'VCE method    = {vce_type}')
    print(f'Derivative    = {deriv}')
    
    print('\n')

    print(f"{'Cutoff c = ':10}{round(cutoff, 3):<10}{'Left of c':>16}{'Right of c':>16}")
    if data is not None:
        print(f"{'Number of obs':20}{nminus_disp:16}{nplus_disp:16}")
        print(f"{'Eff. number of obs':20}{nhl_disp:16}{nhr_disp:16}")
        print(f"{'BW loc. poly.':20}{round(hl, 3):16}{round(hr, 3):16}")
        print(f"{'Order loc. poly.':20}{p:16}{p:16}")

    text_aux = "New sample"
    if cluster is not None:
        print(f"{'Number of clusters':20}{gminus:16}{gplus:16}")
        print(f"{'Eff. # of clusters':20}{gminus_h_l:16}{gplus_h_r:16}")
        text_aux = "New cluster sample"

    print(f"{'Sampling BW':20}{round(hnew_l, 3):16}{round(hnew_r, 3):16}")
    print(f"{text_aux:20}{ntilde_l:16}{ntilde_r:16}")

    print('\n')

    print('=' * 100)
    print(f"{'MDE for power = ':20}{'beta =':>15}{'beta =':>15}{'beta =':>15}{'beta =':>15}{'beta =':>15}")
    print(f"{'              ':20}{round(beta_list[0],3):15}{round(beta_list[1],3):15}{round(beta,3):15}{round(beta_list[2],3):15}{round(beta_list[3],3):15}")
    print('-' * 100)
    print(f"{'Robust bias-corrected':20}{round(mde_rbc_list[0], 3):15}{round(mde_rbc_list[1], 3):15}{round(mde, 3):15}{round(mde_rbc_list[2], 3):15}{round(mde_rbc_list[3], 3):15}")


    if all:
        print('')
        print(f'Conventional', f'{round(mde_conv_list[0], 3)}', f'{round(mde_conv_list[1], 3)}', f'{round(mde_conv, 3)}', f'{round(mde_conv_list[2], 3)}', f'{round(mde_conv_list[3], 3)}')
    print('=' * 100)

    #################################################################
    # Return values
    #################################################################

    output = {
        'mde': mde,
        'se_rbc': se_rbc,
        'sampsi_r': ntilde_r,
        'sampsi_l': ntilde_l,
        'samph_r': hnew_r,
        'samph_l': hnew_l,
        'N_r': nplus,
        'N_l': nminus,
        'Nh_l': n_hnew_l,
        'Nh_r': n_hnew_r,
        'bias_r': bias_r,
        'bias_l': bias_l,
        'Vr_rb': Vr_rb,
        'Vl_rb': Vl_rb,
        'alpha': alpha,
        'beta': beta
    }

    if all:
        output.update({
            'mde_conv': mde_conv,
            'se_conv': se_conv
        })

    return output

