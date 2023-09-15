"""Chris' interpolation method to break down an unknown material into entries in MP."""
import json
from typing import Union, Dict, Set, Tuple, List

import numpy as np
import scipy.optimize
from pymatgen.core import Composition as C
from scipy.optimize import minimize

from s4.data import open_data

__author__ = 'Haoyan Huo'
__email__ = 'haoyan.huo@lbl.gov'
__maintainer__ = 'Haoyan Huo'

__all__ = ['MPUniverseInterpolation']


class MPUniverseInterpolation:
    """
    Interpolate a unknown material to the Materials Project (MP)
    using all compounds in the MP.

    This code is adapted from the original version by Christopher J. Bartel.

    The interpolation is done by optimizing the geometry energy,
    calculated by :math:`-\\exp(-D)` where :math:`D=sqrt(\\sum fractional\\_comp\\_diff^2)`,
    under the constraint that all atoms must conserve during optimization.
    """

    def __init__(self, mp_space=None, mp_data_by_temp=None):
        """
        Constructor.

        :param mp_space: Dictionary of {formula, set of atoms} describing
            the MP universe.
        :param mp_data_by_temp: Dictionary of {temperature: {formula, data}} describing
            the MP data of each entry in mp_space.
        """
        self._mp_space = mp_space
        self._mp_data_by_temp = mp_data_by_temp

    @property
    def mp_space(self):
        """Dictionary describing the MP universe."""
        if self._mp_space is None:
            with open_data('mp_spaces.json') as spaces_f:
                mp_space = json.load(spaces_f)
            self._mp_space = {formula: set(elements) for formula, elements in mp_space.items()}
        return self._mp_space

    @property
    def mp_data_by_temp(self):
        """MP thermodynamic data by temperature."""
        if self._mp_data_by_temp is None:
            with open_data('mp_stability.json') as stability_f:
                mp_data = json.load(stability_f)
            self._mp_data_by_temp = mp_data

        return self._mp_data_by_temp

    def neighbors(self, composition: C, target_space):
        """Compute the neighboring compounds for a given compound."""
        neighbors = [f for f, atom_sets in self.mp_space.items() if atom_sets <= target_space]
        neighbors = [f for f in neighbors if len(f) > 1]
        neighbors = [f for f in neighbors if self.mp_data_by_temp['0'][f]['Ed'] <= 0.3]
        neighbors = [C(x) for x in neighbors]
        neighbors += [C(x.symbol) for x in composition]

        return neighbors

    @staticmethod
    def geometry_energy(composition: C, target_space, neighbors):
        """
        Compute geometry energy given a target compound and a list of neighboring compounds.

        """

        def distance(comp_a: C, comp_b: C):
            total_a, total_b = sum(comp_a.values()), sum(comp_b.values())
            return np.sqrt(sum(
                (comp_a[x] / total_a - comp_b[x] / total_b) ** 2
                for x in target_space
            ))

        # Use geometry distance as energy
        distances = np.array([distance(x, composition) for x in neighbors])
        energy_distance = -np.exp(-distances)
        energies = dict(zip(neighbors, energy_distance))
        return energies

    @staticmethod
    def competing_phases(
            composition: C, target_space: Set[str], neighbors,  # pylint: disable=unused-argument
            neighbor_energies
    ) -> Tuple[Dict, List[C]]:
        """
        Find all competing phases for a target compound given list of neighboring
        compounds.

        :param composition: Composition to search for.
        :param target_space: Set of chemical elements in the target composition.
        :param neighbors: Neighboring compounds.
        :param neighbor_energies: Energies of the neighboring compounds.
        :returns: Energies of the neighboring compounds, and the energy hull.
        """

        # # List hull spaces
        # chem_systems = set(['-'.join(sorted(x.symbol for x in c)) for c in neighbors])
        # chem_spaces = [set(x.split('-')) for x in chem_systems]
        # chem_subspaces = [x for x in chem_spaces if any(x < y for y in chem_spaces)]
        # hull_spaces = list(filter(lambda x: x not in chem_subspaces, chem_spaces))

        # # Sort hull data
        # hull_data = {}
        # for space in hull_spaces:
        #     relevant_compounds = [x for x in neighbor_energies
        #     if set(y.symbol for y in x) <= space]
        #     data = {}
        #     for c in relevant_compounds:
        #         data[c] = {'E': neighbor_energies[c]}
        #
        #     hull_data['_'.join(sorted(space))] = data
        relevant_compounds = [x for x in neighbor_energies
                              if set(y.symbol for y in x) <= target_space]
        hull = {c: {'E': neighbor_energies[c]} for c in relevant_compounds}

        # Find competing phases
        # hull = hull_data['_'.join(sorted(target_space))]

        # Note here we **SHOULD NOT** exclude the phases that are
        # equal to the target phase. This happens when we try to
        # interpolate a known compound, which will lead to the trivial
        # but meaningful solution of {target_comp: {amt: 1}}.
        competing_compounds = list(filter(
            lambda x: set(y.symbol for y in x) <= target_space, hull))
        return hull, competing_compounds

    @staticmethod
    def optimize_energy(composition: C, target_space: Set[str],
                        hull: Dict[C, Dict[str, float]], competing_phases: List[C]
                        ) -> scipy.optimize.OptimizeResult:
        """
        Optimize geometry energy and find the combination of competing phases that
        generate the lowest geometry energy.

        :param composition: Composition of the target material.
        :param target_space: List of chemical elements in the target material.
        :param hull: Dictionary whose keys are phases and values are {'E': energy}.
        :param competing_phases: List of compositions as competing phases.
        :returns: Optimization result that contains the solution.
        """
        coefs = np.array([[phase[x] for x in target_space] for phase in competing_phases]).T
        target = np.array([composition[x] for x in target_space])
        energy = np.array([hull[x]['E'] * sum(x.values()) for x in competing_phases])
        initial_sol = np.full((len(competing_phases),), 0.01)
        max_bound = sum(composition.values())
        bounds = [(0, max_bound) for _ in competing_phases]

        def competing_formation_energy(sol):
            return np.dot(sol, energy)

        constraints = [
            {'type': 'eq',
             'fun': lambda x: np.dot(coefs, x) - target}]
        # Try different precisions until no solution can be found.
        for tol in [1e-5, 1e-4, 1e-3, 5e-3, 1e-2]:
            solution = minimize(
                competing_formation_energy,
                initial_sol,
                method='SLSQP',
                bounds=bounds,
                constraints=constraints,
                tol=tol,
                options={'maxiter': 1000, 'disp': False})
            if solution.success:
                break
        return solution

    def interpolate(self, composition: Union[str, C]) -> Dict[C, Dict[str, float]]:
        """
        Interpolate a composition using compounds in the Materials Project.

        The returned data looks like the following:

        .. code-block:: python

            result = {
                Composition("BaO"): {'amt': 1.0, 'E': 0.0}
            }

        :param composition: The composition to interpolate.
        :returns: Dictionary that contains the decomposed compositions and their information.
        """
        if not isinstance(composition, C):
            composition = C(composition)
        target_space = {x.symbol for x in composition}

        neighbors = self.neighbors(composition, target_space)
        energies = self.geometry_energy(composition, target_space, neighbors)
        hull, competing_phases = self.competing_phases(composition, target_space, neighbors,
                                                       energies)

        trivial_solution = {C(el.symbol): {'amt': amt, 'E': 0} for el, amt in composition.items()}

        if not competing_phases or all(len(x) <= 1 for x in competing_phases) == 1:
            return trivial_solution

        solution = self.optimize_energy(composition, target_space, hull, competing_phases)

        if solution.success:
            eps = 1e-4
            mixture = {
                formula: {'amt': amt, 'E': hull[formula]['E']}
                for amt, formula in zip(solution.x, competing_phases) if amt > eps
            }
            return mixture

        if composition in hull and hull[composition]['E'] > 0:
            # Decompose into elemental compositions
            return trivial_solution

        raise ValueError('Failed to decompose %r into competing phases %r' %
                         (composition, competing_phases))
