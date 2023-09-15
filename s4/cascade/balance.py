"""Utility to quickly balance reactions."""
from functools import lru_cache
from typing import List, Tuple

import numpy
import sympy
from numpy.linalg import LinAlgError
from pymatgen.core import Composition as _C

from s4.types import Comp

__author__ = 'Haoyan Huo'
__email__ = 'haoyan.huo@lbl.gov'
__maintainer__ = 'Haoyan Huo'

__all__ = [
    'CannotBalance',
    'quick_balance',
]


class CannotBalance(ValueError):
    """A reaction cannot be balanced."""


@lru_cache(maxsize=128)
def _quick_balance(precursors: Tuple[Comp], open_comp: Tuple[Comp], target: Comp) \
        -> Tuple[List[float], List[float]]:
    precursors = [_C(x) if not isinstance(x, _C) else x for x in precursors]
    open_comp = [_C(x) if not isinstance(x, _C) else x for x in open_comp]
    target = _C(target) if not isinstance(target, _C) else target

    if target in precursors or target in open_comp:
        raise CannotBalance('Target is included in precursor/open_comp list.')

    elements = set(target)
    for i in precursors + open_comp:
        elements |= set(i)
    elements = list(elements)

    target_coef = numpy.array([round(target[x], 4) for x in elements])
    coefs = []
    for precursor in precursors:
        coefs.append([round(precursor[x], 4) for x in elements])
    for gas in open_comp:
        coefs.append([round(gas[x], 4) for x in elements])

    coefs = numpy.array(coefs).T
    try:
        # Numpy is faster than sympy, but cannot deal with overdetermined equations.
        if coefs.shape[0] != coefs.shape[1]:
            raise LinAlgError('Invalid shape for solve!')

        solution = numpy.linalg.solve(coefs, target_coef)
    except LinAlgError as exception:
        solution, params = sympy.Matrix(coefs).gauss_jordan_solve(sympy.Matrix(target_coef))[:2]
        if params.shape[0] != 0:
            raise CannotBalance(f'Found undetermined variables {params}.') from exception

    solution = [float(x) for x in solution]
    return solution[:len(precursors)], solution[len(precursors):]


def quick_balance(precursors: List[Comp], open_comp: List[Comp], target: Comp) \
        -> Tuple[List[float], List[float]]:
    """
    Balance a chemical reaction.

    :param precursors: List of precursors.
    :param open_comp: List of open compositions.
    :param target: The target material.
    :return: Coefficients for the precursor and open materials.
    """
    precursors = tuple(precursors)
    open_comp = tuple(open_comp)
    return _quick_balance(precursors, open_comp, target)
