"""
Module for computing thermodynamics using calculated databases (Materials
Project).

The general application of this module is to compute thermodynamics quantities
for calculated material entries. The following is the step to do this:

1. Obtain the system of interest using :meth:`query_system_mp` or
    :meth:`query_system_scan`.
2. Obtain the zero-temperature formation enthalpy using
    :meth:`compute_corrected_dhf_mp` or :meth:`compute_corrected_dhf_scan`.
3. Interpolate the finite temperature Gibbs energy of formation using
    :meth:`compute_corrected_dgf_mp` or :meth:`compute_corrected_dgf_scan`.

"""
from .finite_g import finite_dg_correction
from .mp import (
    query_system as query_system_mp, MPThermoDatabase,
    compute_corrected_dgf as compute_corrected_dgf_mp,
    compute_corrected_dhf as compute_corrected_dhf_mp)
from .scan import (
    query_system_scan,
    compute_corrected_dhf as compute_corrected_dhf_scan,
    compute_corrected_dgf as compute_corrected_dgf_scan,
)

__all__ = [
    'finite_dg_correction',

    'query_system_mp',
    'compute_corrected_dgf_mp',
    'compute_corrected_dhf_mp',

    'query_system_scan',
    'compute_corrected_dgf_scan',
    'compute_corrected_dhf_scan',
]
