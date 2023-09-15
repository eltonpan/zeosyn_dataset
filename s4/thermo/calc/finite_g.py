"""Finite Gibbs energy correction by Chris Bartel, implemented by Haoyan Huo."""
import json
from itertools import combinations
from math import log

import scipy.interpolate
from pymatgen.entries.computed_entries import ComputedEntry

from s4.data import open_data

__author__ = 'Haoyan Huo'
__email__ = 'haoyan.huo@lbl.gov'
__maintainer__ = 'Haoyan Huo'

__all__ = [
    'finite_dg_correction',
]

with open_data('Element_mass.json') as _f:
    element_masses = json.load(_f)

with open_data('Element_G.json') as _f:
    element_g = json.load(_f)
    _interp_x, _interp_y = zip(*element_g.items())
    _interp_x = list(map(float, _interp_x))
    element_g_interp = {
        el: scipy.interpolate.interp1d(_interp_x, [_t[el] for _t in _interp_y], kind='quadratic')
        for el in _interp_y[0]
    }


def finite_dg_correction(mp_entry: ComputedEntry, temperature: float, dhf: float) -> float:
    """
    Compute finite-temperature :math:`dG(T)` correction using Chris Bartel's method,
    see [Chris2018]_.


    :param mp_entry: The entry to a Materials Project entry, which must contain the
        volume of the structure.
    :param temperature: Finite temperature for which :math:`dG(T)` is approximated.
    :param dhf: Zero-temperature formation enthalpy.
    :returns: Interpolated gibbs energy of formation at finite temperature.

    .. [Chris2018] Bartel, Christopher J., et al. "Physical descriptor for the Gibbs energy
        of inorganic crystalline solids and temperature-dependent materials chemistry."
        Nature communications 9.1 (2018): 1-10.
    """
    comp = mp_entry.composition
    natom = sum(comp.values())

    reduced_mass_sum = 0
    for element_a, element_b in combinations(comp.keys(), 2):
        element_a, element_b = element_a.symbol, element_b.symbol
        reduced_mass = element_masses[element_a] * element_masses[element_b] / \
                       (element_masses[element_a] + element_masses[element_b])
        weight = comp[element_a] + comp[element_b]
        reduced_mass_sum += weight * reduced_mass
    reduced_mass_sum /= (len(comp) - 1) * natom
    vol = mp_entry.data['volume'] / natom

    gdelta = (
            (-2.48e-4 * log(vol) - 8.94e-5 * reduced_mass_sum / vol) * temperature
            + 0.181 * log(temperature) - 0.882
    )

    refs = 0
    for element, fraction in comp.items():
        refs += element_g_interp[element.symbol](temperature) * fraction / natom

    return dhf + gdelta - refs
