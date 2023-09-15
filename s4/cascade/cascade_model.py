"""The cascade model implementation."""
import logging
from functools import reduce
from itertools import combinations
from operator import or_
from typing import Dict, List, Optional, NamedTuple, Tuple, Any

import numpy
from pymatgen.core import Composition as C

from s4.cascade.balance import quick_balance
from s4.cascade.greedy_resolver import find_competing_reactions
from s4.thermo.calc.mp import query_system
from s4.thermo.calc.scan import query_system_scan
from s4.types import Vessel

__author__ = 'Haoyan Huo'
__email__ = 'haoyan.huo@lbl.gov'
__maintainer__ = 'Haoyan Huo'

__all__ = [
    'find_competing_reactions',
    'find_cascade_reaction',
]


def unique(items):
    """Return unique items in a list."""
    return list(set(list(items)))


def n_metal_atoms(precursors: List[C], open_comps: List[C], target: C):
    """Number of metal cations."""
    open_comp_elements = reduce(or_, map(set, open_comps))
    system_elements = reduce(or_, map(set, precursors)) - open_comp_elements
    return sum([target[x] for x in target if x in system_elements])


def find_cascade_reaction(  # pylint: disable=too-many-arguments
        precursors: List[C], open_comp: List[C], temperature: float,
        target_mixture: Dict[C, Dict[C, float]],
        gas_partials: Optional[Dict[C, float]] = None,
        use_scan=False, only_icsd=True):
    """Find pairwise reaction by comparing energies at the interface."""

    # pylint: disable=too-many-locals
    open_comp_elements = reduce(or_, map(set, open_comp))
    system_elements = reduce(or_, map(set, precursors)) - open_comp_elements

    if use_scan:
        completing_entries = list(
            query_system_scan([str(i) for i in system_elements] + ['O']))
    else:
        completing_entries = [
            x for x in query_system([str(i) for i in system_elements] + ['O'])
            if x.data['e_above_hull'] < 0.05]  # < 50 mev/atom
    if only_icsd:
        completing_entries = list(filter(
            lambda x: len(x.data['icsd_ids']) > 0,
            completing_entries
        ))
    all_targets = unique([x.composition.reduced_formula for x in completing_entries])
    all_targets = [C(x) for x in all_targets]

    logging.info('Identify the reaction with lowest energy from %r ==> %r',
                 precursors, all_targets)
    try:
        best, dgrxn, info = find_competing_reactions(
            precursors=precursors,
            open_comp=open_comp,
            all_targets=all_targets,
            target_mixture=target_mixture,
            atom_set=system_elements,
            temperature=temperature,
            gas_partials=gas_partials,
            use_scan=use_scan)
        logging.info('The best reaction is %s with dGrxn = %.3f ev/atom', best, dgrxn)
        return best, dgrxn, info
    except ValueError:
        logging.info('No balanced reactions found.')
        raise


class CascadeResult(NamedTuple):
    """Cascade model result."""
    vessel: Vessel
    driving_force: float
    stop: bool
    reaction_info: Dict[Tuple[str], Dict[str, Any]]


def cascade_step(  # pylint: disable=too-many-arguments
        precursors_dict: Dict[str, float],
        open_comp: List[str],
        temperature: float,
        target_mixture: Dict[str, Dict[str, float]] = None,
        gas_partials: Optional[Dict[str, float]] = None,
        use_scan=False,
        only_icsd=True,
        return_k=None):
    """
    Perform a cascade model step. Note: here the IO materials are all strings.

    Returned driving forces are in ev/step.

    :param precursors_dict: A dictionary containing amounts of precursors.
    :param open_comp: List of open compositions.
    :param temperature: Temperature of the vessel
    :param target_mixture: If present, specify a list of possible mixtures of
        targets that are MP composition interpolations. Note: target_mixture
        must be self-consistent.
    :param gas_partials: Dictionary containing partial pressures of gases phases.
    :param use_scan: Whether to use SCAN database.
    :param only_icsd: Whether to only use matched ICSD structures.
    :param return_k: If not None, return top_k results sorted by reaction energy.
    :return:
    """

    # pylint: disable=too-many-locals
    intermediate_energies = {}
    reactions = {}
    all_info = {}

    open_comp = [C(x) for x in open_comp]
    target_mixture = target_mixture or {}
    target_mixture = {C(t): {C(s): a for s, a in v.items()} for t, v in target_mixture.items()}
    gas_partials = {C(x): a for x, a in (gas_partials or {}).items()}

    precursors_dict = Vessel(precursors_dict)
    for pairwise in combinations(precursors_dict, 2):
        try:
            intermediate, energy, info = find_cascade_reaction(
                precursors=[C(x) for x in pairwise],
                open_comp=open_comp,
                temperature=temperature,
                target_mixture=target_mixture,
                gas_partials=gas_partials,
                use_scan=use_scan,
                only_icsd=only_icsd)
        except ValueError:
            continue

        # p_amt, o_amt = info['_reaction_eqs'][0]
        # amt_reaction = min([
        #     precursors_dict[x] / max(1e-3, y)
        #     for x, y in zip(pairwise, p_amt)])

        all_info[tuple(pairwise)] = info
        if intermediate is None:
            continue

        # The reaction does not know the overall composition, so sort by ev/m.a.
        intermediate_energies[intermediate] = energy  # * amt_reaction
        reactions[intermediate] = pairwise

    if len(intermediate_energies) == 0:
        logging.info('No reactions could be found, stopping')
        return CascadeResult(
            vessel=precursors_dict,
            driving_force=float('nan'),
            stop=True,
            reaction_info=all_info)

    candidates = sorted(intermediate_energies, key=lambda x: intermediate_energies[x])
    top_k_results = []

    while candidates:
        candidate = candidates.pop(0)
        precursor_amounts, _ = quick_balance(reactions[candidate], open_comp, candidate)

        if any(x < 0 for x in precursor_amounts):
            continue

        logging.info('The cascade step chooses %s (%.3f ev/metal_atom) from %r',
                     candidate, intermediate_energies[candidate], intermediate_energies)

        new_precursors_dict = Vessel(**precursors_dict)
        amt_reaction = min([new_precursors_dict[x] / max(1e-3, y) for x, y in
                            zip(reactions[candidate], precursor_amounts)])

        for precursor, amt in zip(reactions[candidate], precursor_amounts):
            amt_consumed = amt_reaction * amt
            new_precursors_dict[precursor] -= amt_consumed
            if numpy.isclose(new_precursors_dict[precursor], 0):
                del new_precursors_dict[precursor]

        new_precursors_dict[candidate] = new_precursors_dict.get(candidate, 0) + amt_reaction
        logging.info('Current precursors: %r', new_precursors_dict)

        # open_comp_elements = reduce(or_, map(set, open_comp))
        # system_elements = reduce(or_, map(set, [C(x) for x in reactions[candidate]])
        # ) - open_comp_elements
        # target = C(candidate)
        tgt_atom = n_metal_atoms(
            precursors=[C(x) for x in reactions[candidate]],
            open_comps=open_comp,
            target=C(candidate)
        )

        # print(intermediate_energies[candidate], amt_reaction, precursors_dict, )
        result = CascadeResult(
            vessel=new_precursors_dict,
            # reaction energy is in ev/m.a.
            driving_force=intermediate_energies[candidate] * tgt_atom * amt_reaction,
            stop=len(new_precursors_dict) == 1,
            reaction_info=all_info)
        if return_k is None:
            return result

        top_k_results.append(result)
        if len(top_k_results) >= return_k:
            break

    if not top_k_results:
        logging.info('All reaction precursor amounts are negative, stopping')
        return CascadeResult(
            vessel=precursors_dict,
            driving_force=float('nan'),
            stop=True,
            reaction_info=all_info)

    return top_k_results
