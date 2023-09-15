"""Greedy finding competing pairwise reactions."""
import logging
from typing import List, Optional, Dict

import numpy
from pymatgen.core import Composition as C

from s4.cascade.thermo import reaction_driving_force
from s4.thermo.constants import ATM_GAS_PARTIALS

__author__ = 'Haoyan Huo'
__email__ = 'haoyan.huo@lbl.gov'
__maintainer__ = 'Haoyan Huo'

__all__ = [
    'find_competing_reactions',
]


def find_competing_reactions(  # pylint: disable=too-many-arguments
        precursors: List[C],
        open_comp: List[C],
        all_targets: List[C],
        target_mixture: Dict[C, Dict[C, float]],
        atom_set,
        temperature: float,
        gas_partials: Optional[Dict[C, float]],
        use_scan=False):
    """
    Find competing pairwise reactions at a interface by finding largest driving force.

    :param precursors: List of precursors to consider.
    :param open_comp: List of open compositions such as gas.
    :param all_targets: List of possible target compounds.
    :param target_mixture: Composition mixture.
    :param atom_set: Set of "non-volatile" atoms.
    :param temperature: Temperature of the reaction.
    :param gas_partials: Dictionary containing gas partial pressures.
    :param use_scan: Whether to use SCAN database.
    """

    # pylint: disable=too-many-locals
    _gas_partials = ATM_GAS_PARTIALS
    _gas_partials.update(gas_partials or {})

    candidate_targets = []
    reaction_strings = []
    dgrxn = []
    reaction_eqs = []

    targets_pool = target_mixture.copy()
    for target in all_targets:
        targets_pool[target] = {target: 1}

    for target, mixture in targets_pool.items():
        thermo = reaction_driving_force(
            precursors, open_comp,
            target, mixture, atom_set,
            temperature, _gas_partials, use_scan)
        if not thermo:
            continue

        reaction_strings.append(thermo.reaction_string)
        logging.info('dGrxn: %.3f ev/atom ==> %s',
                     thermo.driving_force, thermo.reaction_string)
        dgrxn.append(thermo.driving_force)
        candidate_targets.append(target.formula.replace(' ', ''))
        reaction_eqs.append((
            [x[1] for x in thermo.reactants],
            [x[1] for x in thermo.gases],
        ))

    if not dgrxn:
        msg = 'None of the targets can be balanced. %r + %r ==> %r' % (
            precursors, open_comp, all_targets)
        logging.debug(msg)
        raise ValueError(msg)

    sort_index = numpy.array(dgrxn).argsort()
    candidate_targets = [candidate_targets[i] for i in sort_index]
    dgrxn = [dgrxn[i] for i in sort_index]
    reaction_strings = [reaction_strings[i] for i in sort_index]
    reaction_eqs = [reaction_eqs[i] for i in sort_index]

    info = {
        'candidates': candidate_targets,
        'dgrxn (ev/m.a.)': [round(x, 3) for x in dgrxn],
        '_dgrxn (ev/m.a.)': dgrxn,
        'reactions': reaction_strings,
        '_reaction_eqs': reaction_eqs,
    }

    # Find the first reaction with all positive precursor amounts
    rxn_idx = 0
    for i in range(len(sort_index)):
        if all(x > 0 for x in reaction_eqs[i][0]):
            rxn_idx = i
            break
        logging.debug("Skipping reaction %d because precursor amount is negative", i)
    return candidate_targets[rxn_idx], dgrxn[rxn_idx], info
