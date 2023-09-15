from unittest import TestCase

from pymatgen.core import Composition as C

from s4.tmr.interp import MPUniverseInterpolation


class TestCarbonates(TestCase):
    def setUp(self):
        self.interp = MPUniverseInterpolation()

    def test_known_compounds(self):
        compounds = [
            ('Fe', [('Fe', 1.)]),  # Pure element
            ('TiO2', [('TiO2', 1.)]),  # Known compounds binary
            ('Ti2O4', [('TiO2', 2.)]),
            ('Ti0.5O', [('TiO2', 0.5)]),
            ('BaCO3', [('BaCO3', 1.)]),  # Known compounds ternary

            ('Ba0.2Sr0.8CO3', [('BaCO3', 0.2), ('SrCO3', 0.8)])
        ]
        for compound, solution in compounds:
            x = self.interp.interpolate(compound)
            self.assertEqual(len(x), len(solution), "Mixture count not equal")
            for mix, amt in solution:
                self.assertAlmostEqual(
                    amt, x[C(mix)]['amt'],
                    msg="Mixture amount not equal")
