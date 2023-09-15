"""Compute finite-temperature thermodynamic quantities using materials project data."""
import logging
import os
import pickle
import time
from functools import lru_cache
from typing import List, Union

import numpy as np
from filelock import FileLock
from pymatgen.core import Composition as C
from pymatgen.core import Element
from pymatgen.ext.matproj import MPRester
from pymatgen.core.periodic_table import _pt_data
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.ext.matproj import MPRestError

from s4.cache import CACHE_DIR
from s4.thermo.calc.finite_g import finite_dg_correction
from s4.thermo.constants import CO3_CORRECTION, MP_OXYGEN_CORRECTION
from s4.thermo.utils import convert_unit
from s4.types import Comp

__author__ = 'Haoyan Huo'
__email__ = 'haoyan.huo@lbl.gov'
__maintainer__ = 'Haoyan Huo'

__protocol_version__ = 20210316

__all__ = [
    'compute_corrected_dhf',
    'compute_corrected_dgf',
    'query_system',
    'MPThermoDatabase',
    'database',
]

rester = MPRester()

Carbon = Element('C')
Oxygen = Element('O')


def compute_corrected_dhf(mp_entry: ComputedEntry) -> float:
    """
    Compute corrected :math:`\\Delta H_f(0K)` in ev/atom for a material in
    the Materials Project.

    We apply two corrections as follows:

    1. Elemental contributions as used in the Materials Project.
    2. For carbonates, we apply a correction for CO32- anion (fitted by Haoyan).

    :param mp_entry: A computed entry retrieved from the Materials Project.
    :return: The corrected formation enthalpy at zero Kelvin.
    """
    comp = mp_entry.composition

    dhf = mp_entry.energy

    # Fix CO3
    if Carbon in comp and np.isclose(comp.get(Oxygen, 0) / comp[Carbon], 3):
        dhf -= MP_OXYGEN_CORRECTION * comp[Oxygen]
        dhf += CO3_CORRECTION * comp[Carbon]

    for element, fraction in comp.items():
        pure = query_system(element.symbol)[0]
        dhf -= fraction * pure.energy / sum(pure.composition.values())

    return dhf / sum(comp.values())


def compute_corrected_dgf(mp_entry: ComputedEntry, temperature: float = 300) -> float:
    """
    Compute corrected :math:`\\Delta G_f(T)` in ev/atom for a material in
    the Materials Project.

    We apply three corrections as follows:

    1. Elemental contributions as used in the Materials Project.
    2. For carbonates, we apply a correction for CO32- anion (fitted by Haoyan).
    3. Chris Bartels' finite :math:`dG_f(T)` correction.

    :param mp_entry: A computed entry retrieved from the Materials Project.
    :param temperature: Temperature for which :math:`dG_f(T)` will be computed.
    :return: The corrected gibbs energy of formation at the temperature specified.
    """
    if len(mp_entry.composition) == 1:
        return 0
    dhf = compute_corrected_dhf(mp_entry)
    dgf = finite_dg_correction(mp_entry, temperature, dhf)
    return dgf


def _cache_valid(content):
    if not isinstance(content, dict):
        return False
    if content.get('protocol_version', 0) < __protocol_version__:
        return False

    return True


@lru_cache(maxsize=128)
def _read_or_download_cache(system):
    cache_dir = os.path.join(
        CACHE_DIR,
        'mp_cache')

    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir)

    cache_fn = os.path.join(cache_dir, f'{system}.pickle')

    def _return_cache():
        if os.path.exists(cache_fn):
            with open(cache_fn, 'rb') as cache_f:
                content = pickle.load(cache_f)
                if _cache_valid(content):
                    return content['data']
        raise ValueError('Cache missing or invalid.')

    try:
        return _return_cache()
    except (ValueError, pickle.UnpicklingError, EOFError):
        pass

    lock = FileLock(os.path.join(cache_dir, 'write.lock'))

    max_retry = 5
    for i in range(max_retry):
        with lock:
            # Someone else who locked the file may have already written this
            # file. So try to read again.
            try:
                return _return_cache()
            except (ValueError, pickle.UnpicklingError, EOFError):
                pass

        try:
            start = time.time()
            result = rester.get_entries(
                system,
                inc_structure='final',
                sort_by_e_above_hull=True,
                property_data=['volume', 'icsd_ids'])
            logging.debug('Querying MP for %s took %.3f seconds.', system, time.time() - start)

            with lock:
                with open(cache_fn, 'wb') as cache_write_f:
                    pickle.dump({
                        'data': result,
                        'protocol_version': __protocol_version__
                    }, cache_write_f)

            return result
        except MPRestError:
            if i == max_retry - 1:
                raise
            time.sleep(i + 1)


def query_system(system: Union[List[str], str]) -> List[ComputedEntry]:
    """
    Query the Materials Project to get a list of compositions within a system or a specific
    composition. There are multiple ways of specifying a system:

    1. List of chemical elements, e.g., `Ba-Ti-O`.
    2. Exact composition, e.g., `BaCO3`.

    :param system: Chemical system identifier.
    :returns: List of matched entries in the Materials Project.
    """
    if isinstance(system, list):
        assert len([Element(x) for x in system]) > 0
        system_s = '-'.join(sorted(system))
    elif system in _pt_data:
        system_s = system
    elif isinstance(system, str):
        composition = C(system)
        assert len(composition) > 0
        system_s = '-'.join(sorted(x.symbol for x in composition))
        data = _read_or_download_cache(system_s)

        target_comp = composition.reduced_formula
        entries = list(filter(
            lambda x: x.composition.reduced_formula == target_comp,
            data))
        return entries
    else:
        raise TypeError(f'Invalid material specifier, '
                        f'must be a list of elements or a '
                        f'chemical formula, got {system}')

    return _read_or_download_cache(system_s)


class MPThermoDatabase:
    """
    Thermo database by using MP data.
    Just to keep a same interface as ExpThermoDatabase.
    """

    def __str__(self):
        return '<Calculated thermo database using MP data>'

    __repr__ = __str__

    @staticmethod
    def dgf(composition: Comp, temperature: float, unit='ev') -> float:
        """Compute dGf formation gibbs free energy at finite temperature"""
        composition = C(composition)
        entry = query_system(composition.reduced_formula)[0]
        value = compute_corrected_dgf(entry, temperature=temperature)

        if unit == 'ev/atom':
            return value

        value *= sum(composition.values())
        return convert_unit(value, unit, fromunit='ev')

    @staticmethod
    def dhf(composition: Comp,
            # Keep compatible with other APIs
            temperature: float,  # pylint: disable=unused-argument
            unit='ev') -> float:
        """Compute dHf formation enthalpy at finite temperature"""
        composition = C(composition)
        entry = query_system(composition.reduced_formula)[0]
        value = compute_corrected_dhf(entry)

        if unit == 'ev/atom':
            return value

        value *= sum(composition.values())
        return convert_unit(value, unit, fromunit='ev')


database = MPThermoDatabase()
