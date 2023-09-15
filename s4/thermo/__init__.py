"""
Patch pymatgen hash function in order to get it work better than lru_cache().
"""

from pymatgen.core import Composition


def pymatgen_composition_hash(self):
    """
    Better hash function that pymatgen.Composition.__hash__
    """
    hashcode = 0
    for element, amt in self.items():
        hashcode += element.Z ** 3
        hashcode += int(round(amt, 4) * 1000)

    return hash(hashcode)


Composition.__hash__ = pymatgen_composition_hash
