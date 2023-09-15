"""Thermodynamic calculations of the cascade model."""
import logging
from functools import lru_cache, reduce
from math import log
from operator import add
from typing import List, Dict, Set, Mapping, Tuple, Optional

import dataclasses as dataclasses
import numpy
from monty.fractions import gcd
from pymatgen.core import Element, Composition
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.util.string import formula_double_format
from scipy.interpolate import interp1d

from s4.cascade.balance import Comp, quick_balance
from s4.thermo.calc.mp import query_system, compute_corrected_dgf
from s4.thermo.calc.scan import (
    query_system_scan, compute_corrected_dgf as compute_corrected_dgf_scan
)
from s4.thermo.constants import ATM_GAS_ENTROPIES, RT
from s4.thermo.exp.freed import database
from s4.thermo.utils import convert_unit, as_composition as C

__author__ = 'Haoyan Huo'
__email__ = 'haoyan.huo@lbl.gov'
__maintainer__ = 'Haoyan Huo'

__all__ = [
    'reduce_formula',
    'get_dgf_fu',
    'get_gas_mu',
    'ReactionDrivingForce',
    'reaction_driving_force',
]

atom_gas_entropy_interp = {
    x: interp1d(
        data['x'],
        convert_unit(numpy.array(data['y']), fromunit='j', unit='ev'))
    for x, data in ATM_GAS_ENTROPIES.items()
}


def reduce_formula(sym_amt: Mapping[Element, float], iupac_ordering: bool = False) -> \
        Tuple[str, int]:
    """
    Faster implementation of pymatgen.periodic_table.reduce_formula.

    The original pymatgen implementation is too slow. For example,
    some conversions between string and Element are not necessary. Since we
    will call this function for significantly many times, we need to optimize
    this function as much as possible.

    :param sym_amt: Dictionary that contains {Elements: amount}
    :param iupac_ordering: Whether to use IUPAC ordering.
    :return: The reduced composition as string and the factor of reduction.
    """
    elems = [(x.X, x) for x in sym_amt.keys()]
    syms = [x[1] for x in sorted(elems)]
    syms = list(filter(lambda x: abs(sym_amt[x]) > Composition.amount_tolerance, syms))

    factor = 1
    # Enforce integers for doing gcd.
    if all((int(i) == i for i in sym_amt.values())):
        factor = abs(gcd(*(int(i) for i in sym_amt.values())))

    polyanion = []
    # if the composition contains a poly anion
    if len(syms) >= 3 and syms[-1].X - syms[-2].X < 1.65:
        poly_sym_amt = {syms[i]: sym_amt[syms[i]] / factor for i in [-2, -1]}
        (poly_form, poly_factor) = reduce_formula(
            poly_sym_amt, iupac_ordering=iupac_ordering
        )

        if poly_factor != 1:
            polyanion.append("({}){}".format(poly_form, int(poly_factor)))

    syms = syms[: len(syms) - 2 if polyanion else len(syms)]

    if iupac_ordering:
        syms = sorted(syms, key=lambda x: [x.iupac_ordering, x])

    reduced_form = []
    for s in syms:
        normamt = sym_amt[s] * 1.0 / factor
        reduced_form.append(s.symbol)
        reduced_form.append(formula_double_format(normamt))

    reduced_form = "".join(reduced_form + polyanion)
    return reduced_form, factor


@lru_cache(maxsize=512)
def get_reduced_formula(composition: Composition, iupac_ordering: bool = False):
    """
    Faster implementation of Composition.get_reduced_formula.

    :param composition: Composition to reduce.
    :param iupac_ordering: Whether to use IUPAC ordering.
    """
    all_int = all(
        abs(x - round(x)) < Composition.amount_tolerance for x in composition.values()
    )
    if not all_int:
        return composition.formula.replace(" ", ""), 1
    d = {k: int(round(v)) for k, v in composition.items()}
    formula, _ = reduce_formula(d, iupac_ordering=iupac_ordering)

    if formula in Composition.special_formulas:
        formula = Composition.special_formulas[formula]

    return formula


@lru_cache(maxsize=512)
def _get_mp_entry(composition: Composition) -> ComputedEntry:
    """
    Fetch the first Materials Project entry matching the composition.

    :param composition: Composition to match.
    """
    system = [x.symbol for x in composition]
    entries = query_system(system)

    def same_comp(a, b):
        if len(a) != len(b) or any(x not in b for x in a):
            return False
        if a == b:
            return True

        return get_reduced_formula(a) == get_reduced_formula(b)

    entries = list(filter(lambda x: same_comp(x.composition, composition), entries))

    # Return the first one
    if len(entries) == 0:
        raise ValueError(f'No such composition {composition.reduced_formula} in MP!')
    return entries[0]


@lru_cache(maxsize=512)
def _get_scan_entry(composition):
    system = [x.symbol for x in composition]
    entries = query_system_scan(system)
    entries = list(filter(
        lambda x: x.composition.get_reduced_formula_and_factor()[0]
                  == composition.get_reduced_formula_and_factor()[0], entries))

    # Return the first one
    if len(entries) == 0:
        raise ValueError(f'No such composition {composition.reduced_formula} in SCAN db!')
    return entries[0]


def get_dgf_fu(composition: Comp, temperature: float,
               use_mp=False, use_scan=False) -> float:
    """
    Get the gibbs formation energy for a material. If `use_mp=True`, the value is from
    Materials Project with finite gibbs formation energy interpolation. Otherwise, it's
    obtained from FREED (experimental) database.

    :param composition: Composition of the material.
    :param temperature: Temperature of the material.
    :param use_mp: Whether to use MP data.
    :param use_scan: Whether to use SCAN data.
    :return: Gibbs formation energy of the compound.
    """
    composition = C(composition)

    # Pure elements has no dGf
    if len(composition) == 1:
        return 0.

    if use_mp:
        return compute_corrected_dgf(
            _get_mp_entry(composition), temperature
        ) * sum(composition.values())

    if use_scan:
        return compute_corrected_dgf_scan(
            _get_scan_entry(composition), temperature
        ) * sum(composition.values())

    return database.dgf(composition, temperature, unit='ev', allow_extrapolate=True)


def get_gas_mu(composition: Comp, temperature: float, fugacity: float = 1.0) -> float:
    """
    Compute chemical potential of gas. Enthalpy values are from the FREED experimental database.
    Entropy values are from NIST data files.

    :param composition: Composition of the gas.
    :param temperature: Temperature of the gas.
    :param fugacity: Fugacity (partial pressure) of the gas.
    :return: Gas molecule chemical potential.
    """
    composition = C(composition)
    # No extrapolate since the database covers most of the temperatures
    enthalpy = database.h(composition, temperature, allow_extrapolate=False)
    entropy = atom_gas_entropy_interp[composition](temperature)

    return (
            enthalpy - temperature * entropy +
            convert_unit(RT, fromunit='j', unit='ev') * temperature * log(fugacity)
    )


@dataclasses.dataclass
class ReactionDrivingForce:
    """
    Calculated driving force of a inorganic synthesis reaction.
    """
    reaction_string: str
    driving_force: float

    reactants: List[Tuple[Comp, float]]
    gases: List[Tuple[Comp, float]]


def reaction_driving_force(  # pylint: disable=too-many-arguments
        precursors: List[Comp], open_comp: List[Comp],
        target: Comp, target_mixture: Dict[C, float],
        atom_set: Set[Element], temperature: float, gas_partials: Dict[C, float],
        use_scan: bool) -> Optional[ReactionDrivingForce]:
    """
    Balance reaction and compute grand canonical driving force at once. The computed
    driving forces are in ev/metal_atom. Note that the reaction equation is normalized
    to have 1 metal atom per target composition.

    :param precursors: List of precursors of the reaction.
    :param open_comp: List of open compositions of the reaction.
    :param target: Target material composition.
    :param target_mixture: Dictionary containing the mixtures of target material.
    :param atom_set: Set of "non-volatile" atoms.
    :param temperature: Temperature of the reaction.
    :param gas_partials: Dictionary containing gas partial pressures.
    :param use_scan: Whether to use SCAN database.
    :return: Calculated reaction driving force, or None if no reaction can be balanced.
    """

    # pylint: disable=too-many-locals
    try:
        p_amt, o_amt = quick_balance(precursors, open_comp, target)
    except ValueError:
        logging.debug('Warning: skip target %s because I cannot balance %r, %r ==> %s',
                      target.reduced_formula, precursors, open_comp, target.reduced_formula)
        return None

    # Compute thermodynamics
    target_dgf = 0
    target_comp = []
    for comp, amount in target_mixture.items():
        comp = C(comp)
        target_dgf += get_dgf_fu(comp, temperature, use_mp=not use_scan,
                                 use_scan=use_scan) * amount
        target_comp.append(comp * amount)
    delta_g = target_dgf
    target_comp = reduce(add, target_comp)
    if target != target_comp:
        raise ValueError(f'Target composition {target} does not match mixture {target_comp}.')

    precursor_dg_contrib = []
    for amt, precursor in zip(p_amt, precursors):
        contrib = get_dgf_fu(precursor, temperature, use_mp=not use_scan, use_scan=use_scan)
        precursor_dg_contrib.append(contrib)
        delta_g -= amt * contrib

    # We maximize the driving force of grand canonical potential, instead of
    # gibbs energy of reaction. This is because the system is open to gas molecules.
    open_comp_dg_contrib = []
    for amt, gas in zip(o_amt, open_comp):
        contrib = get_gas_mu(gas, temperature, gas_partials[C(gas)])
        # contrib = try_get_dgf_fu(o, temperature)
        # contrib += convert_unit(rt, fromunit='j', unit='ev'
        # ) * temperature * log(_gas_partials[C(o)])
        open_comp_dg_contrib.append(contrib)
        delta_g -= amt * contrib

    tgt_atoms = sum([target[x] for x in target if x in atom_set])
    factor = tgt_atoms
    # factor = 1
    delta_g /= factor

    reaction = []
    for amt, precursor, contrib in zip(p_amt, precursors, precursor_dg_contrib):
        if abs(amt) > 1e-3:
            reaction.append(
                '%.2f %s (%.3f ev)' % (round(amt / factor, 2), precursor.reduced_formula, contrib))
    for amt, gas, contrib in zip(o_amt, open_comp, open_comp_dg_contrib):
        if abs(amt) > 1e-3:
            reaction.append(
                '%.2f %s (%.3f ev)' % (round(amt / factor, 2), gas.reduced_formula, contrib))

    reaction_string = '%s == %.2f %s (%.3f ev) (DF_rxn=%.3f ev/atom)' % (
        ' + '.join(reaction), 1 / factor, target.reduced_formula, target_dgf, delta_g)

    return ReactionDrivingForce(
        reaction_string=reaction_string,
        driving_force=delta_g,
        reactants=[(precursor, amt / factor) for precursor, amt in zip(precursors, p_amt)],
        gases=[(gas, amt / factor) for gas, amt in zip(open_comp, o_amt)],
    )
