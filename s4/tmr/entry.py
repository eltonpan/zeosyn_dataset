"""Entries in text-mined dataset"""
import os
import sys
import traceback
from collections import defaultdict
from dataclasses import dataclass
from multiprocessing import Pool
from typing import Dict, List, Optional, Union, Tuple

import numpy as np
from pymatgen.core import Composition as _C
from pymatgen.core.composition import CompositionError
from tqdm import tqdm

from s4.cascade.balance import quick_balance
from s4.thermo.constants import ATM_GAS_PARTIALS
from s4.tmr.thermo import (
    MPInterpolatedMaterial, GasMaterial, ExpDeterminedMaterial, InvalidMPIdError)

__author__ = 'Haoyan Huo'
__email__ = 'haoyan.huo@lbl.gov'
__maintainer__ = 'Haoyan Huo'

__all__ = [
    'MaterialWithEnergy',
    'ReactionEnergies',
    'ReactionEntry',
    'from_reactions',
    'from_reactions_multiprocessing',
]

possible_gases = {
    _C('CO2'),
    _C('O2'),
    # _C('NH4'),
    # _C('H2O'),
}


@dataclass
class MaterialWithEnergy:
    """
    A material with thermodynamics quantities resolved.
    """

    #: The resolved thermodynamic entity.
    thermo: Union[MPInterpolatedMaterial, GasMaterial]
    #: The composition of this material.
    composition: _C
    #: Whether this material is a target material.
    is_target: bool
    #: Which side of the reaction this material belongs to.
    side: str
    #: The amount of the material.
    amt: float

    def __repr__(self):
        return f'<Material {self.composition.reduced_formula.ljust(13)} ' \
               f'as {"%r" % round(self.amt, 3)} {self.side.rjust(8)}, ' \
               f'dGf_300K={"%.3f" % self.thermo.dgf(300)} ev/at, ' \
               f'dGf_600K={"%.3f" % self.thermo.dgf(600)} ev/at, ' \
               f'dGf_900K={"%.3f" % self.thermo.dgf(900)} ev/at>'

    def __str__(self):
        return self.__repr__()


@dataclass
class ReactionEnergies:
    """
    Reaction with thermodynamics quantities.
    """
    #: The target composition.
    target: _C
    #: Species involved in this reaction.
    species: List[MaterialWithEnergy]
    #: Variable substitutions
    vars_sub: Dict[str, float]

    def dgrxn(self, temperature: float) -> float:
        """
        Compute the Gibbs energy of reaction at the specified temperature.

        :param temperature: Temperature for which the reaction Gibbs energy
            should be calculated.
        :return: Calculated Gibbs energy of reaction at the specified temperature.
        """
        delta_g = 0
        target = next(x.composition for x in self.species if x.is_target)
        target_atom = sum(target.values())
        target_amt = next(x.amt for x in self.species if x.is_target)

        for item in self.species:
            if item.side == 'product':
                delta_g += item.thermo.dgf(temperature) * item.amt * sum(item.composition.values())
            elif item.side == 'reactant':
                delta_g -= item.thermo.dgf(temperature) * item.amt * sum(item.composition.values())
        delta_g /= target_amt * target_atom
        return round(delta_g, 4)

    def dgcrxn(self, temperature: float) -> float:
        """
        Compute the grand canonical energy of reaction at the specified temperature.

        :param temperature: Temperature for which the reaction grand canonical energy
            should be calculated.
        :return: Calculated grand canonical energy of reaction at the specified temperature.
        """
        gcrxm = 0
        target = next(comp.composition for comp in self.species if comp.is_target)
        target_atom = sum(target.values())
        target_amt = next(x.amt for x in self.species if x.is_target)

        for item in self.species:
            if isinstance(item.thermo, GasMaterial):
                value = item.thermo.mu(temperature)
            else:
                value = item.thermo.dgf(temperature)

            if item.side == 'product':
                gcrxm += value * item.amt * sum(item.composition.values())
            elif item.side == 'reactant':
                gcrxm -= value * item.amt * sum(item.composition.values())
        gcrxm /= target_amt * target_atom
        return round(gcrxm, 4)

    @property
    def reaction_string(self) -> str:
        """String representation of the reaction."""
        reactants = []
        products = []
        for item in self.species:
            part = '%r %s' % (round(item.amt, 3), item.composition.reduced_formula)
            if item.side == 'reactant':
                reactants.append(part)
            else:
                products.append(part)
        return ' + '.join(reactants) + ' == ' + ' + '.join(products)

    def __repr__(self):
        return f'<Reaction {self.reaction_string}, Vars: {self.vars_sub}, \n' \
               f'  dGrxn_300K={"%.3f" % self.dgrxn(300)} ev/at, ' \
               f'dGrxn_1000K={"%.3f" % self.dgrxn(1000)} ev/at, ' \
               f'dGrxn_1500K={"%.3f" % self.dgrxn(1500)} ev/at,\n' \
               f'  dGCrxn_300K={"%.3f" % self.dgcrxn(300)} ev/at, ' \
               f'dGCrxn_1000K={"%.3f" % self.dgcrxn(1000)} ev/at, ' \
               f'dGCrxn_1500K={"%.3f" % self.dgcrxn(1500)} ev/at,\n' \
               f'  Materials:\n    ' + '\n    '.join([
            str(x).replace("\n", "\n    ") for x in self.species
        ]) + '>'

    def __str__(self):
        return self.__repr__()


@dataclass
class ReactionEntry:
    """
    An entry derived from a text-mined recipe.
    """
    #: Identifier or key of the reaction.
    k: str  # pylint: disable=invalid-name
    #: Reaction string.
    reaction_string: str
    #: Experimental temperature.
    exp_t: float
    #: Experimental time.
    exp_time: float
    #: List of reactions contained in this recipe.
    reactions: List[ReactionEnergies]

    def experimental_dgrxn(self) -> Dict[_C, float]:
        """Gibbs free energy of the reaction."""
        return {
            reaction.target: reaction.dgrxn(self.exp_t + 273.15)
            for reaction in self.reactions
        }

    def experimental_dgcrxn(self) -> Dict[_C, float]:
        """Grand canonical free energy of the reaction."""
        return {
            reaction.target: reaction.dgcrxn(self.exp_t + 273.15)
            for reaction in self.reactions
        }

    def __repr__(self):
        return f'<k={self.k}, {len(self.reactions)} Reactions from {self.reaction_string}, \n' \
               f' Experimental temperature: {self.exp_t}, time: {self.exp_time}, \n ' + \
               ',\n'.join(str(x) for x in self.reactions).replace('\n', '\n ') + '>'

    def __str__(self):
        return self.__repr__()

    @staticmethod
    def from_reaction_entry(data, key, override_fugacity=None, use_database='mp'):
        """
        Compute ReactionEntry from text-mined data, identified by keys.

        :param data: Dictionary consisting of reaction entries.
        :param key: Key of entry.
        :param override_fugacity: If not None, set the fugacity (effective partial pressure
            of certain gases.
        :param use_database: Which thermodynamic database to use. Choices are {'mp', 'freed'}.
        """

        # pylint: disable=too-many-locals,too-many-statements
        recipe = data[key]

        precursors = recipe['precursors']
        precursor_comp_strings = [x['material_formula'] for x in precursors]
        open_comp = list(possible_gases)

        override_fugacity = override_fugacity or {}
        fugacity = ATM_GAS_PARTIALS.copy()
        fugacity.update({_C(c): v for c, v in override_fugacity.items()})

        if use_database == 'mp':
            cls = MPInterpolatedMaterial
        elif use_database == 'freed':
            cls = ExpDeterminedMaterial
        else:
            raise ValueError(f'Invalid use_database {use_database}, must be from ["mp", "freed"]')

        target_thermos = cls.from_material_dict(recipe['target'], allow_many=True)
        reactions = []
        for target_thermo, thermo_data in zip(target_thermos, recipe['target']['thermo']):
            # Try to balance
            precursor_amts, open_comp_amts = quick_balance(precursor_comp_strings, open_comp,
                                                           target_thermo.composition)

            species = [MaterialWithEnergy(
                thermo=target_thermo, composition=target_thermo.composition,
                is_target=True, side='product', amt=1)]

            for amt, precursor_d in zip(precursor_amts, precursors):
                if np.isclose(amt, 0.0):
                    continue
                p_thermo = cls.from_material_dict(precursor_d)
                assert amt > 0, \
                    f"Amount of precursor {p_thermo.composition.reduced_formula} is negative!"

                species.append(MaterialWithEnergy(
                    thermo=p_thermo, composition=p_thermo.composition,
                    is_target=False, side='reactant', amt=amt))
            for amt, gas in zip(open_comp_amts, open_comp):
                if np.isclose(amt, 0.0):
                    continue

                side = 'product' if amt < 0 else 'reactant'

                if gas not in fugacity:
                    raise ValueError('{%s} is a gas, but there is no fugacity data for it!' % gas)

                species.append(MaterialWithEnergy(
                    thermo=GasMaterial(gas.copy(), fugacity[gas]), composition=gas.copy(),
                    is_target=False, side=side, amt=abs(amt)))

            vars_sub = (thermo_data['amts_vars'] or {}).copy()
            vars_sub.update(thermo_data.get('el_var', {}))
            reactions.append(ReactionEnergies(
                target=target_thermo.composition, species=species, vars_sub=vars_sub
            ))

        def max_exp_t(item):
            temp_val = float(-1)
            for operation in item.get('operations', {}):
                if operation is None:
                    continue
                cond = operation.get('conditions', {})
                if cond is None:
                    continue
                for temp in cond.get('heating_temperature', []):
                    if temp['max_value'] is not None:
                        temp_val = max(temp['max_value'], temp_val)
                    elif temp['values']:
                        temp_val = max(max(temp['values']), temp_val)
            return temp_val

        def max_exp_time(item):
            time_val = float(-1)
            for operation in item.get('operations', {}):
                if operation is None:
                    continue
                cond = operation.get('conditions', {})
                if cond is None:
                    continue
                for time in cond.get('heating_time', []):
                    if time['max_value'] is not None:
                        time_val = max(time['max_value'], time_val)
                    elif time['values']:
                        time_val = max(max(time['values']), time_val)
            return time_val

        if not reactions:
            raise ValueError('No reactions in the list')

        return ReactionEntry(
            k=key,
            reaction_string=data[key]['reaction_string'],
            exp_t=max_exp_t(data[key]),
            exp_time=max_exp_time(data[key]),
            reactions=reactions)


_mp_storage = {}


def from_reactions(data: Dict[str, dict], keys: List[str], *,
                   override_fugacity=None, use_database='mp'):
    """
    Generate list of ReactionEntry.

    :param data: Dictionary containing solid-state synthesis dataset.
    :param keys: List of keys for which reactions are generate.
    :param override_fugacity: Override some of the fugacities.
    :param use_database: Which database to use.
    :return:
    """

    result = []
    for key in keys:
        try:
            entry = ReactionEntry.from_reaction_entry(
                data, key, override_fugacity=override_fugacity,
                use_database=use_database)
            result.append(entry)
        except ValueError:
            pass

    return result


def _mp_initializer(data, use_database, override_fugacity=None, suppress_output=True):
    _mp_storage.update({
        'data': data,
        'use_database': use_database,
        'override_fugacity': override_fugacity,
        'suppress_output': suppress_output
    })
    if suppress_output:
        sys.stdout = open(os.devnull, 'w')  # pylint: disable=consider-using-with
        sys.stderr = open(os.devnull, 'w')  # pylint: disable=consider-using-with


def _mp_worker(key):
    data = _mp_storage['data']
    override_fugacity = _mp_storage['override_fugacity']
    use_database = _mp_storage['use_database']
    try:
        entry = ReactionEntry.from_reaction_entry(
            data, key, override_fugacity=override_fugacity,
            use_database=use_database)
        return {
            'error': False,
            'data': entry
        }
    except (ValueError, TypeError, IndexError, InvalidMPIdError,
            CompositionError, AssertionError) as exp:
        trace = traceback.format_exc()
        return {
            'error': True,
            'traceback': trace,
            'exception': type(exp),
            'exception_desc': str(exp)
        }


def from_reactions_multiprocessing(  # pylint: disable=too-many-arguments
        data: Dict[str, dict], keys: List[str], *,
        processes=16, suppress_output=True, use_database='mp',
        override_fugacity: Optional[Dict[_C, float]] = None,
        return_errors: bool = False) -> Union[Dict[str, ReactionEntry], Tuple]:
    """
    Generate list of ReactionEntry using multiprocessing.

    :param data: Dictionary containing solid-state synthesis dataset.
    :param keys: List of keys for which reactions are generate.
    :param processes: Number of processes to use.
    :param suppress_output: Suppress the output of workers.
    :param use_database: Which database to use.
    :param override_fugacity: Override some of the fugacities.
    :param return_errors: Whether or not to return errors.
    :return:
    """
    reaction_data = {}
    errors = defaultdict(list)
    with Pool(processes=processes,
              initializer=_mp_initializer,
              initargs=(data, use_database, override_fugacity, suppress_output)) as pool:
        for result in tqdm(pool.imap_unordered(_mp_worker, keys), total=len(keys)):
            if result['error']:
                errors[result['exception']].append(result)
            else:
                reaction_data[result['data'].k] = result['data']

    if return_errors:
        return reaction_data, dict(errors)

    return reaction_data
