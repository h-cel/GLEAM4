import re

import numpy as np
import pandas as pd


def extract_site_code(filename):
    """
    Function to extract the official eddy covariance site name (XX-YYY)
    form modified versions

    Paramaters
    ---------
    filename: string
        String containing the filename

    Returns
    -------
    site: string or None
        If the filename contains a valid site name, it is returned.
        Returns None otherwise
    """
    # Regular expression to match the site code pattern
    match = re.search(r"_(\w{2}-\w{3})_", filename)  # e.g. BE-Bra
    match2 = re.search(r"_(\w{2}-\w{3}_V\w{1})_", filename)  # e.g. BE-Bra_V2
    if match:
        if match2:
            site_long = match2.group(1)
            match_long = re.search(r"^(\w{2}-\w{3})(?:_V\w{1})?$", site_long)
            site = match_long.group(1)
        else:
            site = match.group(1)
    else:
        site = None
    return site


def kge(y_true, y_pred, modified=False):
    """
    Calculate the Kling-Gupta Efficiency (KGE) between two time series.

    Parameters
    ----------
    y_true : np.array
        The true (often observed) time series
    y_pred : np.array
        The predicted (often simulated) time series
    modified : bool, optional
        Default is False, calculating the classic KGE, e.g. see:
        https://en.wikipedia.org/wiki/Kling%E2%80%93Gupta_efficiency
        If True, the modified KGE is calculated, see:
        https://doi.org/10.1016/j.jhydrol.2012.01.011

    Returns
    -------
    kge_acc : float
        KGE
    r : float
        Pearson correlation coefficient
    beta : float
        Ratio of the mean of the predicted to the mean of the true time series
    alpha_or_gamma: float
        Ratio of the standard deviation of the predicted to the standard deviation
        of the true time series (alpha) if modified = False. If modified = True,
        the ratio of the coefficient of variations is used instead (gamma)
    """
    # Select where not Nan
    nan_bool_true = np.isnan(y_true)
    nan_bool_pred = np.isnan(y_pred)
    nan_bool = np.logical_or(nan_bool_true, nan_bool_pred)
    y_true, y_pred = y_true[~nan_bool], y_pred[~nan_bool]

    # calculate
    if y_true.size > 0:
        r = np.corrcoef(y_true, y_pred)[0, 1]
        mu_true, mu_pred = np.mean(y_true), np.mean(y_pred)
        sigma_true, sigma_pred = np.std(y_true), np.std(y_pred)
        beta = mu_pred / mu_true
        if modified:
            gamma_or_alpha = (sigma_pred / mu_pred) / (sigma_true / mu_true)
        else:
            gamma_or_alpha = sigma_pred / sigma_true
        kge_acc = 1 - np.sqrt(
            (r - 1) ** 2 + (beta - 1) ** 2 + (gamma_or_alpha - 1) ** 2
        )
    else:
        kge_acc, r, beta, gamma_or_alpha = None, None, None, None
    return kge_acc, r, beta, gamma_or_alpha


def kge_agg_func(x, name_true, name_pred, modified=False):
    """
    Wrapper function around `kge` to apply to a pandas DataFrame.
    For example usage, see
    https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.apply.html

    See also
    --------
    kge: the function this one wraps around
    """
    tuple_out = kge(x[name_true], x[name_pred], modified)
    kge_acc, r, beta, gamma_or_alpha = tuple_out
    result = {"KGE": kge_acc, "correlation": r, "beta": beta}

    if modified:
        result["gamma"] = gamma_or_alpha
    else:
        result["alpha"] = gamma_or_alpha

    return pd.Series(result)
