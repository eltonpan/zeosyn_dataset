"""Utility functions for thermodynamic calculations."""
import scipy.stats

__author__ = 'Haoyan Huo'
__email__ = 'haoyan.huo@lbl.gov'
__maintainer__ = 'Haoyan Huo'

from pymatgen.core import Composition

from s4.types import Comp

_energy_convert_table = {
    'cal': {
        'ev': 0.0000433641,
        'kcal': 1e-3,
        'kj': 0.004184,
        'j': 4.184,
        'cal': 1,
    },
    'ev': {
        'ev': 1,
        'kcal': 23.0609,
        'kj': 96.4869,
        'j': 96486.9,
        'cal': 23060.9,
    },
    'j': {
        'ev': 0.00001036410,
        'kcal': 0.000239001,
        'kj': 0.001,
        'j': 1,
        'cal': 0.239001
    },
    'kj': {
        'ev': 0.01036410,
        'kcal': 0.239001,
        'kj': 1,
        'j': 1000,
        'cal': 239.001
    },
    'kcal': {
        'ev': 0.0433641,
        'kcal': 1,
        'kj': 4.184,
        'j': 4184,
        'cal': 1000,
    }
}


def as_composition(comp: Comp) -> Composition:
    """Return a composition if not."""
    if not isinstance(comp, Composition):
        return Composition(comp)
    return comp


def fit_and_predict(data_x, data_y, x_hat: float) -> float:
    """
    Quickly fit a linear model from data_y = a data_x + b
    and use it to predict for x_hat
    """
    regressor = scipy.stats.linregress(data_x, data_y)

    return x_hat * regressor.slope + regressor.intercept


def convert_unit(value, unit, fromunit='cal'):
    """
    Convert energy from one unit to another unit.

    :param value: Energy value.
    :param unit: The source unit.
    :param fromunit: The target unit.
    :return: Converted energy value.
    """
    if fromunit not in _energy_convert_table:
        raise TypeError(f'Invalid from unit. Possible values '
                        f'are {list(_energy_convert_table.keys())}')
    if unit not in _energy_convert_table[fromunit]:
        raise TypeError(f'Invalid to unit. Possible values '
                        f'are {list(_energy_convert_table[fromunit].keys())}')

    return _energy_convert_table[fromunit][unit] * value
