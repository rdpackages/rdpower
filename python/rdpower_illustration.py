###################################################################
# rdpower: power calculations for RD designs
# Illustration file
# Authors: Matias Cattaneo, Rocio Titiunik, Gonzalo Vazquez-Bare
###################################################################
## NOTE: if you are using rdrobust version 2020 or newer, the option
## masspoints="off" and stdvars="on" may be needed in order to replicate the
## results in the Stata journal article.
## For example, line 39:
##    aux = rdpower(data=Z,tau=5)
## should be replaced by:
##    aux = rdpower(data=Z,tau=5,masspoints="off",stdvars="on")
###################################################################

#################################################################
# Setup
#################################################################

import pandas as pd
import numpy as np
from rdrobust import rdrobust
from rdpower import rdpower, rdsampsi, rdmde

########## To run from the source ##################################################
# (in terminal) pip uninstall rdpower
# import sys
# sys.path.insert(0, '/')
# from rdpower import rdpower
# from rdsampsi import rdsampsi
# from rdmde import rdmde
# from rdrobust import rdrobust
###################################################################################

data = pd.read_csv('rdpower_senate.csv')
Y = data['demvoteshfor2']
R = data['demmv']
Z = data[['demvoteshfor2', 'demmv']]
covs = data[['population', 'dopen', 'dmidterm']]
clustvar = data['state']

################################################################
# RDPOWER
################################################################

#  rdpower against tau = 5
aux = rdpower(data=Z, tau=5)

#  rdpower with covariates
aux = rdpower(data=Z, tau=5, covs=covs)

#  rdpower with plot
aux = rdpower(data=Z, tau=5, plot=True)

#  rdpower with rdrobust options
aux = rdpower(data=Z, tau=5, h=[16, 18], b=[18, 20])
aux = rdpower(data=Z, tau=5, kernel='uniform', cluster=clustvar, vce='hc0')
aux = rdpower(data=Z, tau=5, bwselect='certwo', vce='hc3', scaleregul=0, rho=1)

#  rdpower with conventional inference

aux = rdpower(data=Z, tau=5, all=True)

#  rdpower with user-specified bias and variance

rdr = rdrobust(Y, R)
samph = rdr.bws.values[0]
sampsi = rdr.N_h
bias = np.array(rdr.bias)/(samph ** 2)

N = sum(rdr.N)

Vl_rb = N * samph[0] * rdr.V_rb_l[0, 0]
Vr_rb = N * samph[1] * rdr.V_rb_r[0, 0]
variance = [Vl_rb, Vr_rb]

aux = rdpower(data=Z, tau=5, bias=bias, variance=variance, samph=samph, sampsi=sampsi)

#  rdpower manually increasing variance by 20%

rdr = rdrobust(Y, R)
samph = rdr.bws.values[0]

N = sum(rdr.N)

Vl_rb = N * samph[0] * rdr.V_rb_l[0, 0] * 1.2
Vr_rb = N * samph[1] * rdr.V_rb_r[0, 0] * 1.2
variance = [Vl_rb, Vr_rb]

aux = rdpower(data=Z, tau=5, variance=variance)

# rdpower without data

aux1 = rdpower(data=Z, tau=5)

aux = rdpower(tau=5, nsamples=[aux1['N_l'], aux1['Nh_l'], aux1['N_r'], aux1['Nh_r']],
              bias=[aux1['bias_l'], aux1['bias_r']], variance=[aux1['Vl_rb'], aux1['Vr_rb']],
              sampsi=[aux1['sampsi_l'], aux1['sampsi_r']], samph=[aux1['samph_l'], aux1['samph_r']])

# comparing exp-post power across specifications

aux1 = rdpower(data=Z, tau=5, p=1, h=20, plot=True)
aux2 = rdpower(data=Z, tau=5, p=2, h=20, plot=True)
aux3 = rdpower(data=Z, tau=5, p=1, plot=True)
aux4 = rdpower(data=Z, tau=5, p=2, plot=True)

# rdpower with clustering
aux = rdpower(data=Z, tau=5, cluster=clustvar)

#################################################################
## RDSAMPSI
#################################################################

# rdsampsi with tau = 5
aux = rdsampsi(data=Z, tau=5)

# rdsampsi setting bandwidth and nratio with plot
aux = rdsampsi(data=Z, tau=5, beta=0.9, samph=[18, 19], nratio=0.5, plot=True)

# rdsampsi with conventional inference
aux = rdsampsi(data=Z, tau=5, all=True)

# rdsampsi vs rdpower
aux1 = rdsampsi(data=Z, tau=5)
aux2 = rdpower(data=Z, tau=5, sampsi=[aux1['sampsi.h.l'], aux1['sampsi.h.r']])

# rdsampsi without data
aux1 = rdsampsi(data=Z, tau=5)
aux2 = rdsampsi(tau=5, nsamples=[aux1['N.l'], aux1['Nh.l'], aux1['N.r'], aux1['Nh.r']],
                bias=[aux1['bias.l'], aux1['bias.r']], variance=[aux1['var.l'], aux1['var.r']],
                samph=[aux1['samph.l'], aux1['samph.r']], init_cond=aux1['init.cond'])

# comparing sample sizes across designs

aux1 = rdsampsi(data=Z, tau=5, p=0, h=20, plot=True)
aux2 = rdsampsi(data=Z, tau=5, p=1, h=20, plot=True)
aux3 = rdsampsi(data=Z, tau=5, p=0, plot=True)
aux4 = rdsampsi(data=Z, tau=5, p=1, plot=True)

# rdsampsi with clustering
aux = rdsampsi(data=Z, tau=5, cluster=clustvar)

#################################################################
# RDMDE
#################################################################

# rdmde with default options
aux = rdmde(data=Z)

# rdmde for a power of 0.75
aux = rdmde(data=Z, beta=0.75)

# rdmde manually setting the bandwidths
aux = rdmde(data=Z, samph=[12, 13])

# rdmde manually setting the sample size inside the bandwidths
aux = rdmde(data=Z, sampsi=[350, 320])

# rdmde without data
rdr = rdrobust(Y, R)
samph = rdr.bws.values[0]
bias = np.array(rdr.bias)/(samph ** 2)
N = sum(rdr.N)
Vl_rb = N * samph[0] * rdr.V_rb_l[0, 0]
Vr_rb = Vl_rb*1.1
variance = [Vl_rb, Vr_rb]

aux = rdmde(nsamples=[600, 350, 700, 320], bias=bias, variance=variance, samph=17, init_cond=4)

# compare rdmde and rdpower
aux1 = rdmde(data=Z)
mde = aux1['mde']
aux2 = rdpower(data=Z, tau=mde)

# compare rdmde and rdsampsi
aux1 = rdmde(data=Z)
mde = aux1['mde']
nhl = aux1['Nh_l']
nhr = aux1['Nh_r']
aux2 = rdsampsi(data=Z, tau=mde, nratio=nhr / (nhl + nhr))

# rdmde with clustering
aux = rdmde(data=Z, cluster=clustvar)