from typing import List


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
