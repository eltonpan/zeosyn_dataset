"""Structural features of reactions."""
from pymatgen.analysis.local_env import CrystalNN
# from robocrys.condense.site import SiteAnalyzer


class GeometryAnalyzer:
    """Analyze a material's structure."""
    near_neighbors = CrystalNN()
    transition_metals = {
        'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
        'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
        'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
        'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
        'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No',
        'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn',
    }
    alkali_metals = {
        'Li', 'Na', 'K', 'Rb', 'Cs', 'Fr'
    }
    earth_alkali_metals = {
        'Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra',
    }

    def __init__(self, analyzer_params=None):
        self.analyzer_params = {
            'symprec': 0.01,
            'use_symmetry_equivalent_sites': False
        }
        self.analyzer_params.update(analyzer_params or {})

    def _get_geometry_for_atoms(self, structure, atoms=None):
        bonded_s = self.near_neighbors.get_bonded_structure(structure)
        site_analyzer = SiteAnalyzer(bonded_s, **self.analyzer_params)

        sites_info = {}
        for eq_site in site_analyzer.equivalent_sites:
            if atoms is None or bonded_s.structure[eq_site].specie.symbol in atoms:
                element = str(bonded_s.structure[eq_site].specie)
                geometry = site_analyzer.get_site_geometry(eq_site)

                nn_sites = site_analyzer.get_nearest_neighbors(
                    eq_site, inc_inequivalent_site_index=True)
                nn_indices = [nn_site['inequiv_index'] for nn_site in nn_sites]

                sites_info[eq_site] = {
                    'element': element,
                    'geometry': geometry,
                    'nn': nn_indices,
                }
                # site_analyzer.get_site_summary(eq_site)
        return sites_info

    def get_geometry_all_atoms(self, structure):
        """Get geometry of all atoms."""
        return self._get_geometry_for_atoms(structure)

    def get_geometry_all_metals(self, structure):
        """Get geometry of all metal cations."""
        return self._get_geometry_for_atoms(
            structure, self.transition_metals | self.alkali_metals | self.earth_alkali_metals)

    def get_geometry_transition_metals(self, structure):
        """Get geometry of all transition metal cations."""
        return self._get_geometry_for_atoms(
            structure, self.transition_metals)
