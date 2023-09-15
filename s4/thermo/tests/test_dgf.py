from unittest import TestCase

from s4.thermo.calc.mp import compute_corrected_dgf, query_system
from s4.thermo.exp.freed import database


class TestCarbonates(TestCase):
    def test_alkali(self):
        for temperature in range(300, 1000, 100):
            for material in ['Li2CO3', 'Na2CO3', 'K2CO3', 'Cs2CO3']:
                dgf_exp = database.dgf(material, temperature, unit='ev/atom', allow_extrapolate=True)

                entry = query_system(material)[0]
                dgf_calc = compute_corrected_dgf(entry, temperature)

                self.assertAlmostEqual(
                    dgf_exp, dgf_calc, delta=0.1, msg=f'{material} @ {temperature} diverges')

    def test_earth_alkali(self):
        for temperature in range(300, 1000, 100):
            for material in ['MgCO3', 'CaCO3', 'SrCO3', 'BaCO3']:
                dgf_exp = database.dgf(material, temperature, unit='ev/atom', allow_extrapolate=True)

                entry = query_system(material)[0]
                dgf_calc = compute_corrected_dgf(entry, temperature)

                self.assertAlmostEqual(
                    dgf_exp, dgf_calc, delta=0.1, msg=f'{material} @ {temperature} diverges')
