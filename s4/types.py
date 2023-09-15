"""Generic types used in this package."""
from typing import Union
from pymatgen.core import Composition as C

__author__ = 'Haoyan Huo'
__email__ = 'haoyan.huo@lbl.gov'
__maintainer__ = 'Haoyan Huo'

__all__ = [
    'Comp',
    'Vessel',
]

Comp = Union[C, str]


class Vessel(dict):
    """A reaction vessel."""

    def __repr__(self):
        items = ['%r %s' % (round(v, 4), k) for (k, v) in self.items()]
        return '{%s}' % ', '.join(items)

    def copy(self):
        """Make a copy of itself."""
        return Vessel(self)
