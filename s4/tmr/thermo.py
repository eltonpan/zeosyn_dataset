"""Thermodynamic calculations for text-mined dataset entries."""
import hashlib
import json
import os
import sqlite3
from dataclasses import dataclass
from math import log
from threading import Lock
from typing import List

from pymatgen.core import Composition as C
from pymatgen.entries.computed_entries import ComputedEntry

from s4.cache import CACHE_DIR
from s4.cascade.thermo import get_gas_mu
from s4.structural.polytype import GeometryAnalyzer
from s4.thermo.calc.mp import query_system, compute_corrected_dgf, compute_corrected_dhf
from s4.thermo.constants import RT
from s4.thermo.exp.freed import database
from s4.thermo.utils import convert_unit
from s4.tmr.interp import MPUniverseInterpolation

__author__ = 'Haoyan Huo'
__email__ = 'haoyan.huo@lbl.gov'
__maintainer__ = 'Haoyan Huo'

__all__ = [
    'GasMaterial',
    'MPInterpolatedMaterial',
    'ExpDeterminedMaterial',
    'InvalidMPIdError',
]

sqlite_fn = os.path.join(CACHE_DIR, 'structure_cache.sqlite')
structure_db = sqlite3.connect(sqlite_fn, check_same_thread=False)
structure_db_cur = structure_db.cursor()
structure_db_lock = Lock()

with structure_db_lock:
    structure_db_cur.execute(
        "CREATE TABLE IF NOT EXISTS mp_structures ("
        "id INTEGER PRIMARY KEY AUTOINCREMENT, "
        "mp_id VARCHAR(255) NOT NULL, "
        "structure_hash CHAR(32) NOT NULL, "
        "data TEXT);")

mp_interp = MPUniverseInterpolation()


class InvalidMPIdError(KeyError):
    """The material cannot be associated with a Materials Project ID."""


@dataclass
class GasMaterial:
    """Gaseous material, such as O2, CO2, etc."""
    #: Composition of the gas.
    composition: C
    #: Effective partial pressure of the gas.
    fugacity: float

    def dhf(self, temperature: float) -> float:
        """
        Returns the formation enthalpy of this material in ev/atom.

        :param temperature: Temperature for which formation enthalpy is calculated.
        :return: Formation enthalpy at the specified temperature.
        """
        return database.dhf(self.composition, temperature, unit='ev/atom')

    def dgf(self, temperature: float) -> float:
        """
        Returns the Gibbs energy of formation of this material in ev/atom.

        :param temperature: Temperature for which Gibbs energy of formation is calculated..
        :return: Gibbs energy of formation at the specified temperature.
        """
        if len(self.composition) == 1:
            dgf0 = 0
        else:
            dgf0 = database.dgf(self.composition, temperature, unit='ev/atom')
        return (dgf0 +
                convert_unit(RT, fromunit='j', unit='ev') * temperature * log(self.fugacity) / sum(
                    self.composition.values()))

    def mu(self, temperature: float) -> float:
        """
        Returns the chemical potential of this gas in ev/atom

        :param temperature: Temperature for which chemical potential is calculated..
        :return: Chemical potential at the specified temperature.
        """

        return get_gas_mu(
            self.composition, temperature, self.fugacity
        ) / sum(self.composition.values())


@dataclass
class ExpDeterminedMaterial:
    """Material with experimentally measured (from FREED, see
    :class:`s4.thermo.exp.ExpThermoDatabase`) thermodynamic quantities."""

    #: The composition of the material.
    composition: C

    def dhf(self, temperature: float) -> float:
        """
        Returns the formation enthalpy of this material in ev/atom.

        :param temperature: Temperature for which formation enthalpy is calculated.
        :return: Formation enthalpy at the specified temperature.
        """
        return database.dhf(self.composition, temperature, unit='ev/atom', allow_extrapolate=True)

    def dgf(self, temperature: float) -> float:
        """
        Returns the Gibbs energy of formation of this material in ev/atom.

        :param temperature: Temperature for which Gibbs energy of formation is calculated..
        :return: Gibbs energy of formation at the specified temperature.
        """
        return database.dgf(self.composition, temperature, unit='ev/atom', allow_extrapolate=True)

    @staticmethod
    def from_material_dict(material, allow_many=False):
        """
        Generate from material data. The parameter `material` should contain a field
        named `"thermo"`. The following is an example:

        .. code-block:: python

            material = {
                "thermo": [
                {
                    "interpolation": "BaCO3",
                    "formula": "BaCO3",
                }
            ]}
        """
        assert material.get('thermo', None) is not None, "The material does not have thermo data."
        if len(material['thermo']) > 1 and not allow_many:
            raise ValueError(f'This material has {len(material["thermo"])} thermo entries, '
                             f'please specify you need a list by setting allow_many=True.')

        results = []
        for thermo in material['thermo']:
            if '+' in thermo['interpolation']:
                continue

            composition = C(thermo['formula'])
            if composition not in database.compositions:
                continue

            results.append(ExpDeterminedMaterial(composition=composition))

        if not allow_many:
            return results[0]

        return results


@dataclass
class MPInterpolatedMaterial:
    """Material with thermo quantities interpolated using MP."""
    compositions: List[C]
    amounts: List[float]
    mp_entries: List[ComputedEntry]

    @property
    def polytypes(self):
        """List of polytypes in the structure as predicted by geometry analyzer."""
        analyzer = GeometryAnalyzer()

        polytype_by_element = {}
        for i, (amt, entry) in enumerate(zip(self.amounts, self.mp_entries)):
            structure_json = entry.structure.as_dict()
            md5 = hashlib.md5()
            md5.update(json.dumps(structure_json).encode())
            json_hash = md5.hexdigest()

            with structure_db_lock:
                result = structure_db_cur.execute(
                    'SELECT data FROM mp_structures '
                    'WHERE mp_id = ? AND structure_hash = ?', (entry.entry_id, json_hash)
                ).fetchone()
            if result is None:
                result = analyzer.get_geometry_all_atoms(entry.structure)
                with structure_db_lock:
                    structure_db_cur.execute(
                        'INSERT INTO mp_structures (mp_id, structure_hash, data) '
                        'VALUES (?, ?, ?)', (entry.entry_id, json_hash, json.dumps({
                            'result': result,
                            'structure': structure_json
                        })))
                    structure_db.commit()
            else:
                result = json.loads(result[0])['result']

            for site, data in result.items():
                element = data['element']

                if element not in polytype_by_element:
                    polytype_by_element[element] = []
                data['entry_ind'] = i
                data['site_ind'] = site
                data['%atoms'] = amt
                polytype_by_element[element].append(data)
        for element, data in polytype_by_element.items():
            total = sum(x['%atoms'] for x in data)
            for i in data:
                i['%atoms'] /= total

        return polytype_by_element

    @property
    def composition(self):
        """Total composition."""
        comp = self.compositions[0].copy() * self.amounts[0]
        for amt, i in zip(self.amounts[1:], self.compositions[1:]):
            comp += i * amt
        return comp

    def dhf(self, temperature) -> float:  # pylint: disable=unused-argument
        """
        Returns the dHf of this interpolated material in ev/atom.

        :return:
        """
        value = 0
        weights = 0
        for amt, entry, comp in zip(self.amounts, self.mp_entries, self.compositions):
            value += amt * compute_corrected_dhf(entry) * sum(comp.values())
            weights += amt * sum(comp.values())

        return value / weights

    def dgf(self, temperature: float) -> float:
        """
        Returns the dGf of this interpolated material in ev/atom.

        :param temperature: Temperature of dGf.
        :return:
        """
        value = 0
        weights = 0
        for amt, entry, comp in zip(self.amounts, self.mp_entries, self.compositions):
            value += amt * compute_corrected_dgf(entry, temperature) * sum(comp.values())
            weights += amt * sum(comp.values())

        return value / weights

    @staticmethod
    def from_material_dict(material, allow_many=False):
        """Generate an instance from materials data."""
        assert material.get('thermo', None) is not None, "The material does not have thermo data."

        thermo = material['thermo']
        if len(thermo) > 1 and not allow_many:
            raise ValueError(f'This material has {len(thermo)} thermo entries, '
                             f'please specify you need a list by setting allow_many=True.')

        entries = []
        for entry in thermo:
            interp_result = mp_interp.interpolate(entry['formula'])
            comps = list(interp_result)
            amt = [interp_result[x]['amt'] for x in comps]
            comp_mp_ids = [query_system(x.reduced_formula)[0].entry_id for x in comps]

            # interp = list(map(lambda x: x.strip().split(' '), entry['interpolation'].split('+')))
            # amt, comp = zip(*interp)
            # amt = [float(x) for x in amt]
            # comp = [C(x) for x in comp]
            # comp_mp_ids = entry['mp_ids']
            mp_entries = []

            for comp, mp_id in zip(comps, comp_mp_ids):
                system = [x.symbol for x in comp]
                _mp_entries = query_system(system)
                try:
                    mp_entries.append(next(x for x in _mp_entries if x.entry_id == mp_id))
                except StopIteration as exception:
                    raise InvalidMPIdError(f'No matching entry in MP with id '
                                           f'{mp_id} could be found.') from exception

            entries.append(MPInterpolatedMaterial(
                compositions=comps,
                amounts=amt,
                mp_entries=mp_entries,
            ))

        if not allow_many:
            return entries[0]

        return entries
