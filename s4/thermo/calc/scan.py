"""Compute finite-temperature thermodynamic quantities using SCAN data from MP."""
import os
import pickle
from typing import List, Union

from filelock import FileLock
from maggma.stores import MongoStore
from pymatgen.core import Composition as C
from pymatgen.core import Element
from pymatgen.entries.computed_entries import ComputedEntry

from s4.thermo.calc.finite_g import finite_dg_correction

__author__ = 'Haoyan Huo'
__email__ = 'haoyan.huo@lbl.gov'
__maintainer__ = 'Haoyan Huo'

__all__ = [
    'compute_corrected_dhf',
    'compute_corrected_dgf',
    'query_system_scan',
]

scan_db = MongoStore(database="mp_rk_calculations",
                     collection_name="materials",
                     host="***REDACTED***",
                     port=27017,
                     username="***REDACTED***",
                     password="***REDACTED***",
                     last_updated_field="last_updated",
                     key="task_id")


def compute_corrected_dhf(mp_entry: ComputedEntry) -> float:
    """
    Compute corrected :math:`\\Delta H_f(0K)` in ev/atom for a material in
    the SCAN database.

    :param mp_entry: A computed entry retrieved from the SCAN database.
    :return: The corrected formation enthalpy at zero Kelvin.
    """
    comp = mp_entry.composition
    dhf = mp_entry.energy

    for element, fraction in comp.items():
        pure = query_system_scan(element.symbol)[0]
        dhf -= fraction * pure.energy / sum(pure.composition.values())

    return dhf / sum(comp.values())


def compute_corrected_dgf(mp_entry: ComputedEntry, temperature: float = 300) -> float:
    """
    Compute corrected :math:`\\Delta G_f(T)` in ev/atom for a material in
    the SCAN database.

    We apply only one corrections as follows:

    1. Chris Bartels' finite :math:`dG_f(T)` correction.

    :param mp_entry: A computed entry retrieved from the SCAN database.
    :param temperature: Temperature for which :math:`dG_f(T)` will be computed.
    :return: The corrected gibbs energy of formation at the temperature specified.
    """
    dhf = compute_corrected_dhf(mp_entry)
    dgf = finite_dg_correction(mp_entry, temperature, dhf)
    return dgf


def get_entries(query, sort_by_energy=True, use_gga='r2scan'):
    """Get entries from SCAN database."""
    assert use_gga in {'gga', 'scan', 'r2scan'}

    query.update({
        'entries.%s' % use_gga: {'$exists': True}
    })

    computed_entries = []

    scan_db.connect()
    for entry in scan_db.query(query):
        data = entry["entries"][use_gga]
        repeat = sum(entry['composition'].values()) / sum(entry['composition_reduced'].values())
        data.update({
            'data': {
                'volume': entry["r2scan_structure"]["lattice"]["volume"] / repeat,
                'symmetry': entry["symmetry"]
            }
        })
        computed_entries.append(ComputedEntry.from_dict(data))

    if sort_by_energy:
        computed_entries = sorted(computed_entries, key=lambda x: x.energy)

    return computed_entries


def query_system_scan(system: Union[List[str], str]) -> List[ComputedEntry]:
    """
    Query Ryan's SCAN database to get a list of compositions within a system or a specific
    composition. There are multiple ways of specifying a system:

    1. List of chemical elements, e.g., `Ba-Ti-O`.
    2. Exact composition, e.g., `BaCO3`.

    :param system: Chemical system identifier.
    :returns: List of matched entries in the SCAN database.
    """
    if isinstance(system, list):
        assert len([Element(x) for x in system]) > 0

        system_s = '-'.join(sorted(system))
        query = {'elements': {'$all': system}}
    elif isinstance(system, str):
        composition = C(system)
        assert len(composition) > 0

        system_s = composition.reduced_formula
        query = {'formula_pretty': composition.reduced_formula}
    else:
        raise TypeError(f'Invalid material specifier, '
                        f'must be a list of elements or a '
                        f'chemical formula, got {system}')

    cache_dir = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        '_ryan_db_cache')

    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir)

    cache_fn = os.path.join(cache_dir, f'{system_s}.pickle')

    if os.path.exists(cache_fn):
        with open(cache_fn, 'rb') as cache_f:
            return pickle.load(cache_f)

    lock = FileLock(os.path.join(cache_dir, 'write.lock'))
    with lock:
        # Someone else who locked the file may have already written this
        # file. So try to read again.
        if os.path.exists(cache_fn):
            with open(cache_fn, 'rb') as cache_f:
                return pickle.load(cache_f)

        result = get_entries(query, use_gga='r2scan')

        with open(cache_fn, 'wb') as cache_f:
            pickle.dump(result, cache_f)

        return result
