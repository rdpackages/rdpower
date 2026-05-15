#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.stats import norm


def as_y_r(data):
    """Return outcome and running-variable columns from a DataFrame or array."""
    if hasattr(data, "iloc"):
        ncol = data.shape[1]
    else:
        data = np.asarray(data)
        if data.ndim != 2:
            raise ValueError('data must be a two-column matrix')
        ncol = data.shape[1]

    if ncol < 2:
        raise ValueError('Too few variables specified in data')
    if ncol > 2:
        raise ValueError('Too many variables specified in data')

    if hasattr(data, "iloc"):
        return data.iloc[:, 0].to_numpy(), data.iloc[:, 1].to_numpy()

    return data[:, 0], data[:, 1]


def as_1d_array(value):
    return np.atleast_1d(np.asarray(value))


def valid_mask(Y, R):
    return ~np.isnan(Y) & ~np.isnan(R)


def unique_count(values):
    return np.unique(values).size


def current_vce(vce):
    if vce is None or not isinstance(vce, str):
        return vce

    vce = vce.lower()
    if vce in {"cluster", "nncluster"}:
        raise ValueError(
            "vce='cluster' and vce='nncluster' are no longer RDROBUST "
            "Python/R options. Use cluster=... with vce='cr1', 'cr2', "
            "or 'cr3'."
        )
    return vce


def subset_sample(Y, R, cluster=None, subset=None):
    Y = np.asarray(Y)
    R = np.asarray(R)
    cluster = None if cluster is None else np.asarray(cluster).reshape(-1)

    if cluster is not None and cluster.size != Y.size:
        raise ValueError("cluster must have the same length as data")

    if subset is None:
        return Y, R, cluster

    subset = np.asarray(subset)
    if subset.dtype == bool:
        if subset.size != Y.size:
            raise ValueError("Boolean subset must have the same length as data")
        index = subset
    elif np.issubdtype(subset.dtype, np.integer) or np.issubdtype(subset.dtype, np.floating):
        if (
            not np.all(np.isfinite(subset))
            or np.any(subset < 0)
            or np.any(subset >= Y.size)
            or not np.all(subset == subset.astype(int))
        ):
            raise ValueError(f"Numeric subset must contain integer indices in 0..{Y.size - 1}")
        index = subset.astype(int)
    else:
        raise ValueError("subset must be boolean or integer")

    if cluster is None:
        return Y[index], R[index], None
    return Y[index], R[index], cluster[index]


def two_sided_power(effect, se, z):
    effect = np.asarray(effect)
    out = 1 - norm.cdf(effect / se + z) + norm.cdf(effect / se - z)
    return out.item() if np.ndim(out) == 0 else out
