"""Stuff for analyzing cascades."""
import math
from functools import reduce
from operator import or_, add
from typing import List

from pymatgen.core import Composition as C

from s4.cascade.balance import quick_balance
from s4.cascade.cascade_model import cascade_step
from s4.tmr.entry import ReactionEnergies
from s4.types import Vessel

__author__ = 'Haoyan Huo'
__email__ = 'haoyan.huo@lbl.gov'
__maintainer__ = 'Haoyan Huo'

__all__ = [
    'compute_cascade',
]


def compute_cascade(reaction: ReactionEnergies, step_temperatures: List[float],
                    only_icsd: bool = True, use_scan: bool = False):
    """
    Compute cascade for a reaction.

    :param reaction: The reaction to calculate cascades for.
    :param step_temperatures: List of temperatures for the cascade reactions.
    :param only_icsd: Only use ICSD universe.
    :param use_scan: Use SCAN energies.
    :return: List of pairwise reactions.
    """
    # pylint: disable=too-many-locals

    # Set up the vessel
    vessel = Vessel()
    for specie in reaction.species:
        if not specie.is_target and 'GasMaterial' not in str(specie.thermo.__class__):
            vessel[specie.composition.reduced_formula] = specie.amt

    # Set up the target species
    tgt_specie = next(x for x in reaction.species if x.is_target)

    # Make sure the reaction balances.
    all_p = list(vessel)
    all_p_amt, _ = quick_balance(all_p, ['CO2', 'O2'], tgt_specie.composition)
    amt_tgt = vessel[all_p[0]] / all_p_amt[0]

    tgt_formula = tgt_specie.composition.formula.replace(' ', '')
    if len(tgt_specie.thermo.compositions) == 1:
        target_mixture = None
    else:
        mixture = {}
        overall_comps = []
        for amt, comp in zip(tgt_specie.thermo.amounts, tgt_specie.thermo.compositions):
            mixture[comp.formula.replace(' ', '')] = amt
            overall_comps.append(comp * amt)
        target_mixture = {tgt_formula: mixture}

        if C(list(target_mixture)[0]) != reduce(add, overall_comps):
            raise ValueError('Target mixture is not balanced, check your reaction.')

    open_comp_elements = reduce(or_, map(set, [C('CO2'), C('O2')]))
    system_elements = reduce(or_, map(lambda x: set(C(x)), vessel)) - open_comp_elements
    tgt_natoms = sum(tgt_specie.composition[x] for x in system_elements)

    # Run the vessel reaction.
    steps_passed = 0
    step_information = []
    for temperature in step_temperatures:
        reason = 'cascade: determined by minimizing dG/m.a'
        new_vessel, dgrxn, stop, _ = cascade_step(
            vessel, ['CO2', 'O2'], temperature=temperature,
            target_mixture=target_mixture, only_icsd=only_icsd, use_scan=use_scan)

        # If the target material is ever visited, compute dGrxn directly.
        # if len(vessel) == 2:
        #     p1, p2 = list(info)[0]
        #     for i, (target, (p_amt, o_amt), _dgrxn) in enumerate(zip(
        #             info[(p1, p2)]['candidates'],
        #             info[(p1, p2)]['_reaction_eqs'],
        #             info[(p1, p2)]['_dgrxn (ev/m.a.)'])):
        #         if numpy.isclose(p_amt[0] / p_amt[1], vessel[p1] / vessel[p2]) and \
        #                 C(target).reduced_formula == tgt_specie.composition.reduced_formula:
        #             logging.info('Already visited target material, '
        #             'using this as the final reaction.')
        #             tgt_comp = info[(p1, p2)]['candidates'][i]
        #
        #             p_amt, o_amt = quick_balance([p1, p2], ['CO2', 'O2'], tgt_comp)
        #
        #             new_vessel = vessel.copy()
        #             amt_reaction = new_vessel[p1] / p_amt[0]
        #             new_vessel[tgt_comp] = amt_reaction
        #             del new_vessel[p1], new_vessel[p2]
        #             # print(dgrxn, amt_reaction, p_amt)
        #             target = C(target)
        #             tgt_atoms = sum(target[a] for a in system_elements)
        #             # dgrxn is already in ev/m.a., which is ev/tgt.m.a
        #             # To get ev/step, we multiply amt*tgt_atoms
        #             dgrxn  = _dgrxn * amt_reaction * tgt_atoms
        #             stop = True
        #             reason = 'target matched: one of the visited reactions produce target phase'
        #             break

        if not math.isnan(dgrxn):
            step_information.append({
                # dgrxn is in ev/step, here, we normalize by number of target atoms.
                # the final unit is ev/tgt.m.a.
                'driving_force': dgrxn / tgt_natoms / amt_tgt,
                'temperature': temperature,
                'previous_vessel': vessel,
                'current_vessel': new_vessel,
                'reason': reason
            })
            # print(f'dG = {dgrxn:+.3f}ev/atom @ {temperature}K, {vessel} ==> {new_vessel}')
            steps_passed += 1

        vessel = new_vessel

        if stop:
            break

    # if len(vessel) > 1 and steps_passed < len(step_temperatures):
    #     # If any remaining exists, compute the driving force of the final reaction
    #     targets_pool = target_mixture if target_mixture is not None else {
    #         tgt_formula: {tgt_formula: 1}}
    #     precursors = list(vessel)
    #     target_comp = list(targets_pool)[0]
    #     temperature = step_temperatures[steps_passed]
    #     steps_passed += 1
    #
    #     final_result = driving_force(
    #         [C(x) for x in precursors], [C('O2'), C('CO2')], C(target_comp),
    #         {C(x): y for x, y in targets_pool[target_comp].items()},
    #         atom_set=set(x for x in C(target_comp).keys() if x.symbol != 'O'),
    #         temperature=temperature, gas_partials=atm_gas_partials, use_scan=False)
    #
    #     if final_result is None:
    #         raise ValueError('Cannot complete the last step of the reaction.')
    #
    #     p_amt, o_amt = quick_balance(precursors, ['CO2', 'O2'], target_comp)
    #
    #     amt_reaction = vessel[precursors[0]] / p_amt[0]
    #     # Check the reaction actually balances
    #     for p, _amt in zip(precursors, p_amt):
    #         amt_r_prime = vessel[p] / _amt
    #         if not numpy.isclose(amt_r_prime, amt_reaction, atol=1e-3):
    #             raise ValueError('The reaction vessel does not balance to the target phase')
    #
    #     new_vessel = Vessel({target_comp: amt_reaction})
    #
    #     step_information.append({
    #         'driving_force': final_result["driving_force"] * amt_reaction / amt_tgt,
    #         'temperature': temperature,
    #         'previous_vessel': vessel,
    #         'current_vessel': new_vessel,
    #         'reason': 'force formation of target phase'
    #     })
    #     # print(f'dG = {final_result["driving_force"]:+.3f}ev/atom @ {temperature}K,'
    #     f' {vessel} ==> {new_vessel}')

    return step_information
