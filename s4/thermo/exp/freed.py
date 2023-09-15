"""FREED database."""
import pickle
from dataclasses import dataclass
from math import log, sqrt
from typing import List, Optional, Dict

import numpy
import pandas
from pandas import Series
from pymatgen.core import Composition as _C

from s4.data import open_data
from s4.thermo.utils import convert_unit, fit_and_predict
from s4.types import Comp

__author__ = 'Haoyan Huo'
__email__ = 'haoyan.huo@lbl.gov'
__maintainer__ = 'Haoyan Huo'

__all__ = [
    'EnthalpyEq', 'GibbsFormationEq', 'FREEDEntry', 'ExpThermoDatabase',
    'database'
]


@dataclass
class EnthalpyEq:
    """
    Represents an enthalpy equation in the form of:

    :math:`H(T) - H_0(298K) = AT + BT^2 + \\frac{C}{T} + DT^{\\frac{1}{2}} + ET^3 + F`
    """

    # Keep same name as in FREED database.
    # pylint: disable=too-many-instance-attributes

    #: Coefficient :math:`A`.
    A: float  # pylint: disable=invalid-name
    #: Coefficient :math:`B`.
    B: float  # pylint: disable=invalid-name
    #: Coefficient :math:`C`.
    C: float  # pylint: disable=invalid-name
    #: Coefficient :math:`D`.
    D: float  # pylint: disable=invalid-name
    #: Coefficient :math:`E`.
    E: float  # pylint: disable=invalid-name
    #: Coefficient :math:`F`.
    F: float  # pylint: disable=invalid-name
    #: Transition temperature (upper bound of this equation).
    T_transition: float  # pylint: disable=invalid-name
    #: Transition enthalpy value.
    H_transition: Optional[float]  # pylint: disable=invalid-name
    #: Comments on this enthalpy equation.
    comment: Optional[str]

    def __str__(self):
        return f'<H(T)-H(298K) (below {self.T_transition}K) = {self.A}*T + {self.B}*T^2 + ' \
               f'{self.C}/T + {self.D}*T^0.5 + {self.E}*T^3 + {self.F} cal>'

    def __call__(self, temperature, unit='ev'):
        """Get the enthalpy at temperature."""
        if temperature > self.T_transition:
            raise ValueError(f'Incompatible equation used. '
                             f'Desired {temperature}K, max {self.T_transition}K.')

        dgf = (
                self.A * temperature
                + self.B * temperature ** 2
                + self.C / temperature
                + self.D * sqrt(temperature)
                + self.E * temperature ** 3
                + self.F
        )

        return convert_unit(dgf, unit)

    __repr__ = __str__


@dataclass
class GibbsFormationEq:
    """
    Represents an gibbs formation energy equation in the form of:

    :math:`\\Delta G_f(T) = AT\\ln(T) + BT + CT^2 + \\frac{D}{T} + ET^{\\frac{1}{2}} + FT^3 + G`
    """

    # pylint: disable=too-many-instance-attributes
    #: Coefficient :math:`A`.
    A: float  # pylint: disable=invalid-name
    #: Coefficient :math:`B`.
    B: float  # pylint: disable=invalid-name
    #: Coefficient :math:`C`.
    C: float  # pylint: disable=invalid-name
    #: Coefficient :math:`D`.
    D: float  # pylint: disable=invalid-name
    #: Coefficient :math:`E`.
    E: float  # pylint: disable=invalid-name
    #: Coefficient :math:`F`.
    F: float  # pylint: disable=invalid-name
    #: Coefficient :math:`G`.
    G: float  # pylint: disable=invalid-name
    #: Transition temperature (upper bound of this equation).
    T_transition: float  # pylint: disable=invalid-name

    def __str__(self):
        return f'<dGf (below {self.T_transition}K) = {self.A}*Tln(T) + {self.B}*T + ' \
               f'{self.C}*T^2 + {self.D}*T^-1 + {self.E}*T^0.5 + {self.F}*T^3 + {self.G} cal>'

    def __call__(self, temperature, unit='ev'):
        """Get the dGf at temperature."""
        if temperature > self.T_transition:
            raise ValueError(f'Incompatible equation used. Desired {temperature}K, '
                             f'max {self.T_transition}K.')

        dgf = (
                self.A * temperature * log(temperature)
                + self.B * temperature
                + self.C * temperature ** 2
                + self.D / temperature
                + self.E * sqrt(temperature)
                + self.F * temperature ** 3
                + self.G
        )

        return convert_unit(dgf, unit)

    __repr__ = __str__


@dataclass
class FREEDEntry:
    """A verbatim entry in FREED database."""

    # pylint: disable=too-many-instance-attributes
    #: Formula of the compound.
    formula: str
    #: Description of the material.
    desc: str
    #: Name of the material.
    name: str
    #: Mineral name of the material.
    mineral_name: str

    #: Molar mass of the material.
    molar_mass: float
    #: Density of the material.
    density: float

    #: Standard enthalpy of formation at 298K (unit in cal).
    h0: float  # pylint: disable=invalid-name
    #: Entropy of formation at 298K (unit in cal/kelvin).
    ent0: float
    #: Maximal characterized temperature.
    tmax: float

    #: List of enthalpy equations.
    dh_eqs: List[EnthalpyEq]
    #: List of Gibbs energy equations.
    dgf_eqs: List[GibbsFormationEq]

    def __str__(self):
        return f'<FREEDEntry for {self.formula} {self.desc}. Material name: {self.name}\n' \
               f'  Molar mass: {self.molar_mass}, density: {self.density}, ' \
               f'max thermo temperature {self.tmax}K\n' \
               f'  Enthalpy (298K): {self.h0} cal, Entropy (298K): {self.ent0} cal/kelvin\n' \
               f'  Enthalpy equations:\n    ' + \
               '\n    '.join([str(x) for x in self.dh_eqs]) + \
               '\nGibbs formation energy equations:\n    ' + \
               '\n    '.join([str(x) for x in self.dgf_eqs])

    __repr__ = __str__

    @staticmethod
    def from_row(row: Series):
        """Construct a <FREEDEntry> from a pandas row."""

        # pylint: disable=too-many-locals
        cols = list(row.keys())

        dh_eqs = []
        start = 13
        for i in range(int(row['# Ht Eq'])):
            if i < int(row['# Ht Eq']) - 1:
                (
                    A, B, C, D, E, F,  # pylint: disable=invalid-name
                    T_transition, H_transition, comment  # pylint: disable=invalid-name
                ) = row[cols[start:start + 9]]
                start += 9
            else:
                (
                    A, B, C, D, E, F, T_transition  # pylint: disable=invalid-name
                ) = row[cols[start:start + 7]]
                H_transition, comment = None, None  # pylint: disable=invalid-name
                start += 7
            equation = EnthalpyEq(
                A=A, B=B, C=C, D=D, E=E, F=F,
                T_transition=T_transition, H_transition=H_transition, comment=comment)
            dh_eqs.append(equation)

        dgf_eqs = []
        for i in range(int(row['# Gf Eq'])):
            (
                A, B, C, D, E, F, G, T_transition  # pylint: disable=invalid-name
            ) = row[cols[start:start + 8]]
            start += 8
            equation = GibbsFormationEq(
                A=A, B=B, C=C, D=D, E=E, F=F, G=G,
                T_transition=T_transition)
            dgf_eqs.append(equation)

        return FREEDEntry(
            formula=row['Formula'],
            desc=row['Descr.'],
            name=row['Name'],
            mineral_name=row['Mineral Name'],
            molar_mass=row['GFW'],
            density=row['Dens. (298K)'],
            h0=row['DHf(298K)'],
            ent0=row['S(298 K)'],
            tmax=row['Tmax'],
            dh_eqs=dh_eqs,
            dgf_eqs=dgf_eqs
        )

    def plot_h(self, temp_range=None, unit='ev/atom', ax=None, show=True) -> None:
        """
        Plot enthalpy at various temperatures

        :param temp_range: The range of temperatures. If None, defaults to 298K-TMax.
        :param unit: Unit of the plot
        :param ax: If not None, plot to an existing axis.
        :param show: If True, display this plot after function call.
        """
        import matplotlib.pyplot as plt  # pylint: disable=import-outside-toplevel

        if temp_range is None:
            temp_range = numpy.linspace(298, self.tmax, 10)
        if ax is None:
            _, ax = plt.subplots()

        ax.plot(temp_range, [self.h(x, unit=unit) for x in temp_range], 'x-')
        ax.set_xlabel('Temperature')
        ax.set_ylabel(f'Enthalpy ({unit})')
        if show:
            plt.show()

    def plot_dgf(self, temp_range=None, unit='ev/atom', ax=None, show=True) -> None:
        """
        Plot Gibbs energy of formation at various temperatures

        :param temp_range: The range of temperatures. If None, defaults to 298K-TMax.
        :param unit: Unit of the plot
        :param ax: If not None, plot to an existing axis.
        :param show: If True, display this plot after function call.
        """
        import matplotlib.pyplot as plt  # pylint: disable=import-outside-toplevel

        if temp_range is None:
            temp_range = numpy.linspace(298, self.tmax, 10)
        if ax is None:
            _, ax = plt.subplots()

        ax.plot(temp_range, [self.dgf(x, unit=unit) for x in temp_range], 'x-')
        ax.set_xlabel('Temperature')
        ax.set_ylabel(f'Gibbs formation energy ({unit})')
        if show:
            plt.show()

    def h(self,  # pylint:disable=invalid-name
          temperature: float, unit='ev', allow_extrapolate=False) -> float:
        """
        Compute enthalpy of formation at a temperature.
        If you specify `allow_extrapolate=True` and the temperature falls outside
        of the range defined by FREED database, you will get an extrapolated value.

        :param temperature: Temperature for the enthalpy calculation.
        :param unit: Can be one of {'ev', 'ev/atom', 'kcal', 'kj', 'cal'}
        :param allow_extrapolate: Whether extrapolation is allowed.
        :returns: Enthalpy of formation at the specified temperature.
        """
        if unit == 'ev/atom':
            dhf0 = convert_unit(self.h0, 'ev')
        else:
            dhf0 = convert_unit(self.h0, unit)

        if temperature < 298 or temperature > self.tmax:
            if not allow_extrapolate:
                raise ValueError(f'Temperature {temperature} beyond range (298K, {self.tmax}K) '
                                 f'defined by FREED equations. If an interpolated enthalpy is '
                                 f'acceptable, please specify allow_extrapolate=True')

            temps = numpy.linspace(298, self.tmax, 10)
            return fit_and_predict(
                temps,
                [self.h(i, unit=unit, allow_extrapolate=False) for i in temps],
                temperature)

        for equation in self.dh_eqs:
            try:
                if unit == 'ev/atom':
                    energy = dhf0 + equation(temperature, unit='ev')
                    composition = _C(self.formula)
                    return energy / sum(composition.values())

                return dhf0 + equation(temperature, unit=unit)
            except ValueError:
                continue

        raise ValueError('Something happened')

    def dgf(self, temperature: float, unit='ev', allow_extrapolate=False) -> float:
        """
        Compute gibbs energy of formation at a temperature.
        If you specify `allow_extrapolate=True` and the temperature falls outside
        of the range defined by FREED database, you will get an extrapolated value.

        :param temperature: Temperature for the energy calculation.
        :param unit: Can be one of {'ev', 'ev/atom', 'kcal', 'kj', 'cal'}
        :param allow_extrapolate: Whether extrapolation is allowed.
        :returns: Gibbs energy of formation at the specified temperature
        """
        if temperature < 298 or temperature > self.tmax:
            if not allow_extrapolate:
                raise ValueError(f'Temperature {temperature} beyond range (298K, {self.tmax}K) '
                                 f'defined by FREED equations. If an interpolated gibbs '
                                 f'free energy of formation is acceptable, '
                                 f'please specify allow_extrapolate=True')

            temps = numpy.linspace(298, self.tmax, 10)
            return fit_and_predict(
                temps, [self.dgf(i, unit=unit, allow_extrapolate=False) for i in temps],
                temperature)

        for equation in self.dgf_eqs:
            try:
                if unit == 'ev/atom':
                    composition = _C(self.formula)
                    return equation(temperature, unit='ev') / sum(composition.values())

                return equation(temperature, unit=unit)
            except ValueError:
                continue

        raise ValueError('Something happened')


class ExpThermoDatabase:
    """A database containing all entries in FREED database."""

    def __init__(self, freed_db_fn='./FREED 11.0.xlsm'):
        self.compositions: Dict[str, List[FREEDEntry]] = {}

        data = pandas.read_excel(freed_db_fn, sheet_name='Database', header=1)
        data = data.dropna(subset=('Formula',))
        for _, row in data.iterrows():
            comp = _C(row['Formula'])
            if comp not in self.compositions:
                self.compositions[comp] = []

            self.compositions[comp].append(FREEDEntry.from_row(row))

    def __str__(self):
        return f'<Experimental thermo database containing {len(self.compositions)} compositions>'

    __repr__ = __str__

    def dgf(self, composition: Comp, temperature: float, unit='ev',
            allow_extrapolate=False) -> float:
        """Get the formation gibbs free energy of a compound at desired temperature."""
        return self[composition].dgf(temperature, unit=unit, allow_extrapolate=allow_extrapolate)

    def h(self,  # pylint:disable=invalid-name
          composition: Comp, temperature: float, unit='ev', allow_extrapolate=False) -> float:
        """
        Get the enthalpy of a compound at desired temperature.

        :param composition: Composition of the compound.
        :param temperature: Desired temperature.
        :param unit: Return unit.
        :param allow_extrapolate: Whether an extrapolated value is accepted.
        :return:
        """
        return self[composition].h(temperature, unit=unit, allow_extrapolate=allow_extrapolate)

    def dhf(self, composition: Comp, temperature: float, unit='ev',
            allow_extrapolate=False) -> float:
        """Get the formation enthalpy of a compound at desired temperature."""
        if unit == 'ev/atom':
            per_atom = True
            unit = 'ev'
        else:
            per_atom = False

        if not isinstance(composition, _C):
            composition = _C(composition)

        enthalpy = self.h(composition, temperature, unit=unit, allow_extrapolate=allow_extrapolate)
        for element, amount in composition.items():
            enthalpy -= amount * self.h(element.symbol, temperature, unit=unit,
                                        allow_extrapolate=allow_extrapolate)

        if per_atom:
            enthalpy /= sum(composition.values())

        return enthalpy

    def __getitem__(self, composition: Comp) -> FREEDEntry:
        """Return the polymorph with lowest formation energy (does not always make sense)"""
        if not isinstance(composition, _C):
            composition = _C(composition)

        return sorted(self.compositions[composition], key=lambda x: x.h0)[0]


try:
    with open_data('FREED_parsed.pickle') as f:
        database: Optional[ExpThermoDatabase] = pickle.load(f)
except FileNotFoundError:
    database = None  # pylint: disable=invalid-name
