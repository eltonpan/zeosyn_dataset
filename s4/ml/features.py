"""Computes features for each reaction."""
import json
import re
from collections import defaultdict
from typing import Dict, Tuple

import numpy
from pymatgen.core import Composition as C
from pymatgen.core.composition import CompositionError

from s4.cascade.analysis import compute_cascade
from s4.data import open_data
from s4.tmr.entry import ReactionEnergies
from s4.tmr.thermo import GasMaterial
from s4.thermo.exp.freed import database as exp_database

__all__ = [
    'Featurizer',
]

wp_mp = {}
with open_data('Materials_Melting_Pt.json') as f:
    item = None
    for item in json.load(f):
        if item['value'] is None:
            continue
        try:
            wp_mp[C(item['composition'])] = item['value']
        except (ValueError, CompositionError):
            continue
del item, f
wp_mp_median = numpy.median(list(wp_mp.values()))


class Featurizer:
    """
    This featurizer computes all 133 features used in this work.

    The features are divided into four types:
    1. Precursor properties (melting points, Gibbs formation energy, formation enthalpy, etc.)
    2. Target compositions.
    3. Experiment-adjacent features.
    4. Thermodynamic driving forces.
    """

    #: Chemical elements that qualify as a compositional feature.
    elements = ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'K',
                'Ca', 'Sc', 'Ti', 'V',
                'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Rb', 'Sr',
                'Y', 'Zr', 'Nb',
                'Mo', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'Cs', 'Ba', 'La', 'Ce',
                'Pr', 'Nd', 'Sm',
                'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os',
                'Ir', 'Hg', 'Tl',
                'Pb', 'Bi', 'Th', 'U']

    #: Temperatures for which thermodynamic driving forces are calculated.
    thermo_temps = [x + 273.15 for x in [800, 900, 1000, 1100, 1200, 1300]]

    def __init__(self, thermo_temps=None):
        if thermo_temps is not None:
            self.thermo_temps = thermo_temps

    def _feature_compositional(self, reaction: ReactionEnergies):
        feature_dict = {}
        target_comp = reaction.species[0].composition

        for element in self.elements:
            feature_dict[f'feature_ele_{element}'] = element in target_comp

        return feature_dict

    @staticmethod
    def _feature_textual(tms_data: Dict):
        # pylint: disable=too-many-locals
        text = ''.join(tms_data['ext_paragraph'])

        is_ball_milling = re.search(r'ball[-\s]{0,5}mill', text) is not None
        is_multi_heating = len(re.findall(r'\d+\s*(°C|K)', text)) > 1
        is_regrinding = re.search(
            r'(intermediate|interim|repeated)\s+grind', text,
            flags=re.IGNORECASE) is not None
        is_binder = re.search(r'binder', text, flags=re.IGNORECASE) is not None
        # is_Japan = re.search(r'Japan', text, flags=re.IGNORECASE) is not None
        is_distilled = re.search(r'distill(ed)?', text, flags=re.IGNORECASE) is not None
        is_zirconia = re.search(r'zirconia', text, flags=re.IGNORECASE) is not None
        is_diameter = re.search(r'diameter', text, flags=re.IGNORECASE) is not None
        is_sintered = re.search(r'sinter(ed)?', text, flags=re.IGNORECASE) is not None
        is_wet = re.search(r'wet', text, flags=re.IGNORECASE) is not None
        is_PVA = re.search(r'PVA', text, flags=re.IGNORECASE) is not None
        is_lead = re.search(r'lead', text, flags=re.IGNORECASE) is not None
        is_polycrystal = re.search(r'polycrystal', text, flags=re.IGNORECASE) is not None
        is_phosphors = re.search(r'phosphors', text, flags=re.IGNORECASE) is not None
        is_homoge = re.search(r'homoge', text, flags=re.IGNORECASE) is not None
        is_ground = re.search(r'ground', text, flags=re.IGNORECASE) is not None
        is_again = re.search(r'again', text, flags=re.IGNORECASE) is not None

        feature_dict = {
            'feature_text_is_ball_milling': is_ball_milling,
            'feature_text_is_multi_heating': is_multi_heating,
            'feature_text_is_regrinding': is_regrinding,
            'feature_text_is_binder': is_binder,
            # 'feature_text_is_Japan': is_Japan,
            'feature_text_is_distilled': is_distilled,
            'feature_text_is_zirconia': is_zirconia,
            'feature_text_is_diameter': is_diameter,
            'feature_text_is_sintered': is_sintered,
            'feature_text_is_wet': is_wet,
            'feature_text_is_PVA': is_PVA,
            # 'feature_text_is_lead': is_lead,
            'feature_text_is_polycrystal': is_polycrystal,
            'feature_text_is_phosphors': is_phosphors,
            'feature_text_is_homoge': is_homoge,
            'feature_text_is_ground': is_ground,
            'feature_text_is_again': is_again,
        }

        return feature_dict

    @staticmethod
    def _compute_thermo_bounds(reaction, temp_original, only_icsd=False):
        max_temp = min_temp = temp_original
        true_cascade = compute_cascade(reaction, [temp_original] * 10, only_icsd=only_icsd)

        def cascade_is_good():
            return len(cascade) >= 1 and \
                   len(cascade[-1]['current_vessel']) == 1 and \
                   list(cascade[-1]['current_vessel'])[0] == \
                   list(true_cascade[-1]['current_vessel'])[0]

        for min_temp in numpy.arange(temp_original, 500, -50):
            cascade = compute_cascade(reaction, [min_temp] * 10, only_icsd=only_icsd)
            if not cascade_is_good():
                break
        for max_temp in numpy.arange(temp_original, 3000, 50):
            cascade = compute_cascade(reaction, [max_temp] * 10, only_icsd=only_icsd)
            if not cascade_is_good():
                break

        return min_temp, max_temp

    @staticmethod
    def _get_beta(data_x, data_y):
        if len(data_x) <= 1:
            return float('nan')
        ones = numpy.ones_like(data_x)
        data = numpy.linalg.lstsq(numpy.stack([data_x, ones], axis=1), data_y, rcond=None)[0]
        return data[0]

    def _feature_thermo(
            self, reaction: ReactionEnergies, exp_t,
            calculate_thermo_bounds=True):
        cascades = {}
        for temp in self.thermo_temps:
            cascade = compute_cascade(reaction, [temp] * 10, only_icsd=False)
            cascades[temp] = cascade

        # print(cascades)

        if calculate_thermo_bounds:
            synthesis_temp_lb, synthesis_temp_ub = self._compute_thermo_bounds(
                reaction, exp_t + 273.15)
        else:
            synthesis_temp_lb, synthesis_temp_ub = float('nan'), float('nan')

        feature_dict = {
            'feature_thermo_lb_temp': synthesis_temp_lb,
            'feature_thermo_ub_temp': synthesis_temp_ub,
        }

        lstsq_data = {
            'total': {
                'x': [],
                'y': [],
            },
            'first': {
                'x': [],
                'y': [],
            },
            'last': {
                'x': [],
                'y': [],
            }
        }

        for temp in self.thermo_temps:
            cascade = cascades[temp]
            if len(cascade) == 0:
                total_df = last_df = first_df = float('nan')
            else:
                first_df = cascade[0]['driving_force']
                last_df = cascade[-1]['driving_force']
                total_df = sum(x['driving_force'] for x in cascade)

            if not numpy.isnan(first_df):
                lstsq_data['first']['x'].append(temp)
                lstsq_data['first']['y'].append(first_df)
            if not numpy.isnan(total_df):
                lstsq_data['total']['x'].append(temp)
                lstsq_data['total']['y'].append(total_df)
            if not numpy.isnan(last_df):
                lstsq_data['last']['x'].append(temp)
                lstsq_data['last']['y'].append(last_df)

            feature_dict['feature_thermo_%r_total_df' % temp] = total_df
            feature_dict['feature_thermo_%r_first_df' % temp] = first_df
            feature_dict['feature_thermo_%r_last_df' % temp] = last_df

            feature_dict['feature_thermo_%r_first_df_frac' % temp] = first_df / (total_df + 1e-3)
            feature_dict['feature_thermo_%r_last_df_frac' % temp] = last_df / (total_df + 1e-3)

        feature_dict.update({
            'feature_thermo_total_ddf': self._get_beta(lstsq_data['total']['x'],
                                                       lstsq_data['total']['y']),
            'feature_thermo_first_ddf': self._get_beta(lstsq_data['first']['x'],
                                                       lstsq_data['first']['y']),
            'feature_thermo_last_ddf': self._get_beta(lstsq_data['last']['x'],
                                                      lstsq_data['last']['y'])
        })

        return feature_dict

    @staticmethod
    def _feature_other(reaction: ReactionEnergies, precursors: Tuple[str]):
        n_target_mixture = len(reaction.species[0].thermo.compositions)
        precursors = tuple(sorted(precursors))

        melting_points = [wp_mp.get(C(x), wp_mp_median) for x in precursors]

        enthalpy_values = []
        hs_invalid = False
        gibbs_values = []
        dgf_invalid = False
        for precursor in precursors:
            try:
                enthalpy_values.append(exp_database[precursor].h(300))
            except KeyError:
                hs_invalid = True
            try:
                gibbs_values.append(exp_database[precursor].dgf(300))
            except (KeyError, ValueError):
                dgf_invalid = True

        feature_dict = {
            'feature_exp_min_mp': min(melting_points),
            'feature_exp_max_mp': max(melting_points),
            'feature_exp_mean_mp': numpy.mean(melting_points),
            'feature_exp_div_mp': max(melting_points) - min(melting_points),

            'feature_avg_p_h_300K': numpy.mean(enthalpy_values) if not hs_invalid else float('nan'),
            'feature_min_p_h_300K': numpy.min(enthalpy_values) if not hs_invalid else float('nan'),
            'feature_max_p_h_300K': numpy.max(enthalpy_values) if not hs_invalid else float('nan'),
            'feature_diff_p_h_300K': numpy.max(enthalpy_values) - numpy.min(
                enthalpy_values) if not hs_invalid else float(
                'nan'),
            'feature_avg_dgf_300K': numpy.mean(gibbs_values) if not dgf_invalid else float('nan'),
            'feature_min_dgf_300K': numpy.min(gibbs_values) if not dgf_invalid else float('nan'),
            'feature_max_dgf_300K': numpy.max(gibbs_values) if not dgf_invalid else float('nan'),
            'feature_diff_dgf_300K': numpy.max(gibbs_values) - numpy.min(
                gibbs_values) if not dgf_invalid else float(
                'nan'),

            'feature_syn_n_precursors': len(precursors),
            'feature_syn_n_carbonates': len([x for x in precursors if 'CO3' in x]),
            # 'feature_syn_n_cascades': len(cascade_data),
            'feature_syn_n_target_mixture': n_target_mixture,
        }

        return feature_dict

    def featurize(  # pylint: disable=too-many-arguments
            self,
            reaction: ReactionEnergies, exp_t: float, exp_time: float, k: str = '', i: int = 0,
            tms_data: Dict = None, calculate_thermo_bounds=True) -> Dict[str, float]:
        """
        Featurize a reaction using the 133 features.

        :param reaction: The reaction to compute features for.
        :param k: Key of the recipe.
        :param i: Index of the reaction.
        :param exp_t: Experimental temperature in Celsius.
        :param exp_time: Experimental time.
        :param tms_data: Text-mined dataset dictionary.
        :param calculate_thermo_bounds: Whether to calculate synthesis temperature bounds
            predicted by thermodynamic.
        :returns: Features for this reaction.
        """
        if tms_data is None:
            tms_data = defaultdict(lambda: '')

        precursors = [x.composition.reduced_formula
                      for x in reaction.species if
                      not x.is_target and not isinstance(x.thermo, GasMaterial)]
        precursors = tuple(sorted(precursors))

        entry = {
            'meta_doi': tms_data['doi'],
            'meta_k': k,
            'meta_i': i,
            'meta_precursors': precursors,
            'meta_text': ''.join(tms_data['ext_paragraph']),

            'y_speed': numpy.log10(1. / exp_time),
            'y_temperature': exp_t,
        }
        entry.update(self._feature_compositional(reaction=reaction))
        entry.update(self._feature_textual(tms_data=tms_data))
        entry.update(self._feature_thermo(
            reaction=reaction, exp_t=exp_t,
            calculate_thermo_bounds=calculate_thermo_bounds))
        entry.update(self._feature_other(reaction=reaction, precursors=precursors))
        return entry

    def __call__(self, args):
        return self.featurize(**args)
