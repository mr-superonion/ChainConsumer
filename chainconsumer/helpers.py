# -*- coding: utf-8 -*-
import numpy as np


def mean_weight(x, w):
    return (x * w).sum() / w.sum()


def median_weight(x, w):
    a = np.argsort(x)
    w = w[a]
    x = x[a]
    wc = np.cumsum(w)
    wc /= wc[-1]
    return np.interp(0.5, wc, x)


def percentile_weight(x, w, p):
    a = np.argsort(x)
    w = w[a]
    x = x[a]
    wc = np.cumsum(w)
    wc /= wc[-1]
    return np.interp(p / 100.0, wc, x)


def std_weight(x, w):
    """estimate std with weight"""
    mu = mean_weight(x, w)
    r = x - mu
    return np.sqrt((w * r**2).sum() / w.sum())


def get_extents(data, weight, plot=False, wide_extents=True, tiny=False, pad=False):
    """estimate boundary of the data. The two ends of the boundary is defined
    as 1e-4 or 1e-5 of the cdf and 1-cdf, respectively.
    Parameters:
        data (ndarray):     chain
        weight (ndarray):   wegihts
    """
    hist, be = np.histogram(data, weights=weight, bins=2000)
    bc = 0.5 * (be[1:] + be[:-1])
    cdf = hist.cumsum()
    cdf = cdf / cdf.max()
    icdf = (1 - cdf)[::-1]
    icdf = icdf / icdf.max()
    cdf = 1 - icdf[::-1]
    threshold = 1e-4 if plot else 1e-5
    if plot and not wide_extents:
        threshold = 0.05
    if tiny:
        threshold = 0.3
    i1 = np.where(cdf > threshold)[0][0]
    i2 = np.where(icdf > threshold)[0][0]
    lower = bc[i1]
    upper = bc[-i2]
    if pad:
        width = upper - lower
        lower -= 0.1 * width
        upper += 0.1 * width
    lower = max(lower, np.min(data))
    upper = min(upper, np.max(data))
    return lower, upper


def get_bins(chains):
    proposal = [
        max(
            35,
            np.floor(1.0 * np.power(chain.chain.shape[0] / chain.chain.shape[1], 0.25)),
        )
        for chain in chains
    ]
    return proposal


def get_smoothed_bins(
    smooth, bins, data, weight, marginalised=True, plot=False, pad=False
):
    minv, maxv = get_extents(data, weight, plot=plot, pad=pad)
    if smooth is None or not smooth or smooth == 0:
        return np.linspace(minv, maxv, int(bins)), 0
    else:
        return (
            np.linspace(minv, maxv, int((2 if marginalised else 2) * smooth * bins)),
            smooth,
        )


def get_grid_bins(data):
    bin_c = np.sort(np.unique(data))
    delta = 0.5 * (bin_c[1] - bin_c[0])
    bins = np.concatenate((bin_c - delta, [bin_c[-1] + delta]))
    return bins


def get_latex_table_frame(caption, label):  # pragma: no cover
    base_string = r"""\begin{table}
    \centering
    \caption{%s}
    \label{%s}
    \begin{tabular}{%s}
        %s    \end{tabular}
\end{table}"""
    return base_string % (caption, label, "%s", "%s")
