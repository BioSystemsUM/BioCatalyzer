from typing import List

import numpy as np
import pandas as pd


def match_value(v: float, values: List[float], tol: float = 0.1) -> tuple:
    """
    Match a value with a list of values.

    Parameters
    ----------
    v: float
        The value to match.
    values: List[float]
        The list of values to match.
    tol: float
        The tolerance to match the values.

    Returns
    -------
    tuple:
        The matched value and the index of the matched value.
    """
    for i, value in enumerate(values):
        if value - tol <= v <= value + tol:
            return True, i
    return False, None


def _empty_dfs(dfs: List[pd.DataFrame]):
    """
    Check if at least one dataframe is not empty.

    Parameters
    ----------
    dfs: List[pd.DataFrame]
        The list of dataframes to check.

    Returns
    -------
    bool:
        True if at least one dataframe is not empty. Otherwise, False.
    """
    for r in dfs:
        if not r.empty:
            return False
    return True


def _merge_ec_numbers(x):
    """
    Merge multiple EC numbers.

    Parameters
    ----------
    x: str
        The EC numbers to merge.

    Returns
    -------
    str:
        The merged EC numbers.
    """
    if x == '':
        return np.NaN
    seen = set()
    seen_add = seen.add
    return ';'.join([y for y in x.split(';') if not (y in seen or seen_add(y) or y == '')])


def _merge_fields(value):
    """
    Merge multiple fields.

    Parameters
    ----------
    value: str
        The fields to merge.

    Returns
    -------
    str:
        The merged fields.
    """
    if len(value.split(';')) == 1:
        return value
    seen = set()
    seen_add = seen.add
    return ';'.join([x for x in value.split(';') if not (x in seen or seen_add(x))])
