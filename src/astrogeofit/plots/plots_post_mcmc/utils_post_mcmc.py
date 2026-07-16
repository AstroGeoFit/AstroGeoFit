"""
MIT License

Copyright (c) 2025 [CNRS-Observatoire de Paris]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction...
(Shortened for brevity)
"""

import numpy as np
from scipy.interpolate import interp1d  # type:ignore


def quantile(x, q, weights):
    if not isinstance(weights, np.ndarray):
        weights = np.array(weights)
    weights = weights / np.sum(weights)
    indices = np.argsort(x)
    p = 0
    # for i, w in enumerate(weights[indices]):
    for i in range(weights.shape[0]):
        w = weights[indices[i]]
        p += w
        if p >= q:
            return x[indices[i]]


def weighted_mean(x, w):
    return np.sum(w * x) / np.sum(w)


def weighted_variance(x, w):
    m = weighted_mean(x, w)
    W = np.sum(w)
    sq_dev = np.sum(w * (x - m) ** 2)
    norm_factor = W * (1 - np.sum(w**2) / W**2)
    return sq_dev / norm_factor


def weighted_std(x, w):
    return weighted_variance(x, w) ** 0.5


def get_offset(t_y, t_y_reference, offset_range, num_offsets=1000, return_full=False):
    t, y = t_y
    interp_ref = interp1d(
        *t_y_reference, kind="linear", bounds_error=False, fill_value=np.nan
    )
    min_offset, max_offset = offset_range

    offsets = np.linspace(min_offset, max_offset, num_offsets)
    t_shifted_matrix = t[:, np.newaxis] + offsets
    yref_interp_matrix = interp_ref(t_shifted_matrix)
    valid_mask = ~np.isnan(yref_interp_matrix)
    valid_yref_interp = np.where(valid_mask, yref_interp_matrix, 0)
    valid_y = np.where(valid_mask, y[:, np.newaxis], 0)

    correlations = multi_pearsonr_v2(valid_yref_interp, valid_y)
    best_idx = np.argmax(correlations)
    best_offset = offsets[best_idx]
    best_correlation = correlations[best_idx]

    if return_full:
        return offsets, correlations
    return best_offset, best_correlation


def multi_pearsonr_v2(x, y):
    """
    Compute the Pearson correlation coefficient for each column of y with x.

    Parameters:
    x: array-like, shape (n, m)
    y: array-like, shape (n, m)

    Returns:
    corrs: array, shape (m,)
    """
    mx = x.mean(axis=0)
    my = y.mean(axis=0)
    xm, ym = x - mx, y - my
    r_num = np.add.reduce(xm * ym, axis=0)
    r_den = np.sqrt(
        np.multiply(np.add.reduce(xm * xm, axis=0), np.add.reduce(ym * ym, axis=0))
    )
    return r_num / r_den
