import numpy as np
from rdkit.Chem import AllChem
from ast import literal_eval

class MorganDictVectorizer(object):
    '''
    Morgan fingerprint generator for organic molecules Source: https://programtalk.com/vs4/python/hcji/PyFingerprint/PyFingerprint/heteroencoder.py/
    '''
    def __init__(self, radius=2, augment=None):
        self.radius = radius
        self.augment = augment #Not used
        self.dims = None
        
    def fit(self, mols):
        """Analyses the molecules and creates the key index for the creation of the dense array"""
        keys=set()
        for mol in mols:
            fp = AllChem.GetMorganFingerprint(mol,self.radius)
            keys.update(fp.GetNonzeroElements().keys())
        keys = list(keys)
        keys.sort()
        self.keys= np.array(keys)
        self.dims = len(self.keys)
        
    def transform_mol(self, mol, misses=False):
        """ transforms the mol into a dense array using the fitted keys as index
        
            :parameter mol: the RDKit molecule to be transformed
            :parameter misses: wheter to return the number of key misses for the molecule
         """
        assert type(self.keys) is np.ndarray, "keys are not defined or is not an np.array, has the .fit(mols) function been used?"
        #Get fingerprint as a dictionary
        fp = AllChem.GetMorganFingerprint(mol,self.radius)
        fp_d = fp.GetNonzeroElements()
        
        #Prepare the array, and set the values
        #TODO is there a way to vectorize and speed up this?
        arr = np.zeros((self.dims,))
        _misses = 0
        for key, value in fp_d.items():
            if key in self.keys:
                arr[self.keys == key] = value
            else:
                _misses = _misses + 1
        
        if misses:
            return arr, _misses
        else:
            return arr
    
    def transform(self, mols, misses=False):
        """Transforms a list or array of RDKit molecules into a dense array using the key dictionary (see .fit())
        
        :parameter mols: list or array of RDKit molecules
        :parameter misses: Wheter to return the number of key misses for each molecule
        """
        arr = np.zeros((len(mols), self.dims))
        if misses:
            _misses = np.zeros((len(mols),1))
            for i, mol in enumerate(mols):
                arr[i,:], _misses[i] = self.transform_mol(mol, misses=misses)
            return arr, _misses
        else:
            for i, mol in enumerate(mols):
                arr[i,:] = self.transform_mol(mol, misses=False)
            return arr
        
def rename_disordered_interrupted(truncated_code):
    '''
    Maps 3-letter zeolite IZA code to codes with "*" and/or "-" (disordered/interrupted)

    Args:
        truncated_code (str): 3-letter zeolite IZA code without "*" and/or "-"
        
    Returns:
        str: 3-letter zeolite IZA code with "*" and/or "-"
    '''

    disordered_interrupted = {
                          'BEA': '*BEA',
                          'MRE': '*MRE',
                          'STO': '*STO',
                          'SVY': '*-SVY',
                          'CTH': '*CTH',
                          'ITV': '-ITV',
                          'CLO': '-CLO',
                          'ITN': '*-ITN',
                          'UOE': '*UOE',
                          'IRY': '-IRY',
                          'SVR': '-SVR',
                          'LIT': '-LIT',
                          'IFT': '-IFT',
                          'SFV': '*SFV',
                          'IFU': '-IFU',
                         }
    if truncated_code in disordered_interrupted.keys():
        return disordered_interrupted[truncated_code]
    else:
        return truncated_code
    
def adjacent_values(vals, q1, q3):
    '''
    Helper function for plotting boxplot whiskers are defined by the interquartile range (IQR: Q3-Q1) and the adjacent values.
    '''

    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value

def rmsd_matrix(A, B, squared=False, axis = 0):
    """
    Compute all pairwise distances between vectors in A and B.

    Parameters
    ----------
    A : np.array
        shape should be (M, K)
    B : np.array
        shape should be (N, K)

    Returns
    -------
    D : np.array
        A matrix D of shape (M, N).  Each entry in D i,j represnets the
        distance between row i in A and row j in B.

    See also
    --------
    A more generalized version of the distance matrix is available from
    scipy (https://www.scipy.org) using scipy.spatial.distance_matrix,
    which also gives a choice for p-norm.
    """
    M = A.shape[0]
    N = B.shape[0]

    assert A.shape[1] == B.shape[1], f"The number of components for vectors in A \
        {A.shape[1]} does not match that of B {B.shape[1]}!"

    A_dots = (A*A).sum(axis=1).reshape((M,1))*np.ones(shape=(1,N))
    B_dots = (B*B).sum(axis=1)*np.ones(shape=(M,1))
    if axis == 0:
        D_squared =  1/M*(A_dots + B_dots -2*A.dot(B.T))
    elif axis == 1:
        D_squared =  1/N*(A_dots + B_dots -2*A.dot(B.T))

    if squared == False:
        zero_mask = np.less(D_squared, 0.0)
        D_squared[zero_mask] = 0.0
        return np.sqrt(D_squared)

    return D_squared

def smoothsegment(seg, Nsmooth=100):
    '''
    Helper function for plotting circular dendrogram
    '''

    return np.concatenate([[seg[0]], np.linspace(seg[1], seg[2], Nsmooth), [seg[3]]])

def get_pore_type(pore_size):
    if (pore_size > 0.) & (pore_size <= 8.):
        return 'Small'
    elif (pore_size > 8.) and (pore_size <= 10.):
        return 'Medium'
    elif (pore_size > 10.) and (pore_size <= 12.):
        return 'Large'
    elif (pore_size > 12.):
        return 'Extra-large'
    else:
        return 'Failed'
    
def clean_cbus(cbu_str):
    # This fixes cbus column contains a str representation of list, instead of list itself
    if type(cbu_str) == str: # if str
        return literal_eval(cbu_str)
    else: # if not str, eg. NaN
        cbu_str

# Column names for t-SNE plots
tsne_cols = [
 'Si', 'Al','P','Na','K','Li','Sr','Rb','Cs','Ba','Ca','F','Ge','Ti','In','B','H2O','sda1','OH','cryst_time','cryst_temp',
 'osda_asphericity_mean_0','osda_axes_mean_0','osda_axes_mean_1','osda_bertz_ct_mean_0',
 'osda_eccentricity_mean_0','osda_free_sasa_mean_0',
 'osda_gyration_radius_mean_0','osda_inertial_shape_factor_mean_0', 'osda_npr1_mean_0', 'osda_npr2_mean_0',
 'osda_num_bonds_mean_0', 'osda_num_rot_bonds_mean_0', 'osda_pmi1_mean_0','osda_pmi2_mean_0','osda_pmi3_mean_0','osda_spherocity_index_mean_0','osda_volume_mean_0','osda_mol_weight','osda_formal_charge',
]

# Continuous zeolite descriptors
zeo_cols_cont = [
 'num_atoms',
 'a',
 'b',
 'c',
 'alpha',
 'beta',
 'gamma',
 'volume',
 'td',
 'rdls',
 'td_10',
 'ring_size_0',
 'ring_size_1',
 'ring_size_2',
 'ring_size_3',
 'ring_size_4',
 'ring_size_5',
 'ring_size_6',
 'framework_density',
 'largest_free_sphere',
 'accessible_volume_izc',
 'largest_free_sphere_a',
 'largest_free_sphere_b',
 'largest_free_sphere_c',
 'largest_free_sphere_izc',
 'largest_included_sphere',
 'largest_free_sphere_a_izc',
 'largest_free_sphere_b_izc',
 'largest_free_sphere_c_izc',
 'largest_included_sphere_a',
 'largest_included_sphere_b',
 'largest_included_sphere_c',
 'largest_included_sphere_fsp',
 'largest_included_sphere_izc',
 'num_atoms_per_vol',
'unitcell_vol',
 'density',
 'asa_a2',
 'asa_m2_cm3',
 'asa_m2_g',
 'nasa_a2',
 'nasa_m2_cm3',
 'nasa_m2_g',
 'chan_0_sa_a2',
 'chan_1_sa_a2',
 'chan_2_sa_a2',
 'chan_3_sa_a2',
 'chan_4_sa_a2',
 'chan_5_sa_a2',
 'chan_0_dim',
'sa_num_pockets',
 'av_a3',
 'av_frac',
 'av_cm3_g',
 'nav_a3',
 'nav_frac',
 'nav_cm3_g',
 'chan_0_vol_a3',
 'chan_1_vol_a3',
 'chan_2_vol_a3',
 'chan_3_vol_a3',
 'chan_4_vol_a3',
 'chan_5_vol_a3',
 'avol_num_pockets',
 'poav_a3',
 'poav_frac',
 'poav_cm3_g',
 'ponav_a3',
 'ponav_frac',
 'ponav_cm3_g',
 'probe_rad',
 'N_points',
 'probe_ctr_A_fract',
 'probe_ctr_NA_fract',
 'A_fract',
 'NA_fract',
 'narrow_fract',
 'ovlpvfract',
 'deriv_mean',
 'deriv_variance',
 'deriv_skewness',
 'deriv_kurtosis',
 'cum_mean',
 'cum_variance',
 'cum_skewness',
 'cum_kurtosis',
 'num_si',
]

# Discrete zeolite descriptors
zeo_cols_disc = [
 'chan_num_channels',
 'sa_num_channels',
 'avol_num_channels',
 'isdisordered',
 'isinterrupted',
]

cols_to_drop = [
    'doi',
    'normed',
    'seed',
    'aging_time',
    'aging_temp',
    'rotation',
    'Seed_type',
    'react_vol',
    'pH',
    'osda1',
    'osda2',
    'osda3',
    'product1',
    'product2',
    'product3',
    'precursors',
    'brands',
    'Si/Al',
    'yield',
    'percent cryst',
    'crystal size',
    'micropore volume',
    'micropore diameter',
    'bet area',
    'external surface area',
    'Notes',
    'title',
    'abstract_keywords',
    'recipe_keywords',
    'osda1 synonyms',
    'osda2 synonyms',
    'osda3 synonyms',
    'osda1 iupac',
    'osda2 iupac',
    'osda3 iupac',
    'osda1 smiles',
    'osda2 smiles',
    'osda3 smiles',
    'osda1 formula',
    'osda2 formula',
    'osda3 formula',
    'Code1',
    'Code2',
    'Code3',
    'year',
]

# OSDA features for classification model
osda_cols = [
 'asphericity',
 'axes',
 'bertz_ct',
 'binding',
 'eccentricity',
 'formal_charge', 
 'free_sasa',
 'getaway',
 'gyration_radius',
 'inertial_shape_factor',
 'mol_weight',
 'npr1',
 'npr2',
 'num_bonds', 
 'num_rot_bonds',
 'pmi1',
 'pmi2',
 'pmi3',
 'spherocity_index',
 'volume',
 'whim',
]

X_cols = {
          # Gel composition
         'Si': 'Si',
         'Al': 'Al',
         'P' : 'P',
         'Na': 'Na',
         'K' : 'K',
         'Li': 'Li',
         'Sr': 'Sr',
         'Rb': 'Rb',
         'Cs': 'Cs',
         'Ba': 'Ba',
         'Ca': 'Ca',
         'F' : 'F',
         'Ge': 'Ge',
         'Ti': 'Ti',
         'B' : 'B',
         'Mg': 'Mg',
         'Ga': 'Ga',
         'Zn': 'Zn',
         'Be': 'Be',
         'W' : 'W',
         'Cu': 'Cu',
         'Sn': 'Sn',
         'Zr': 'Zr',
         'V' : 'V',
         'H2O': 'H$_2$O',
         'sda1': 'OSDA amount',
         'OH': 'OH',
          
          # Reaction conditions
         'cryst_time': 'Crystallization time',
         'cryst_temp': 'Crystallization temp',
          
          # OSDA descriptors
         'osda1_asphericity_mean_0': 'OSDA asphericity',
         'osda1_axes_mean_0': 'OSDA axis 1',
         'osda1_axes_mean_1': 'OSDA axis 2',
         'osda1_formal_charge': 'OSDA charge',
         'osda1_free_sasa_mean_0': 'OSDA SASA',
         'osda1_mol_weight': 'OSDA molecular weight',
         'osda1_npr1_mean_0': 'OSDA NPR 1',
         'osda1_npr2_mean_0': 'OSDA NPR 2',
         'osda1_num_rot_bonds_mean_0': 'OSDA rotatable bonds',
         'osda1_pmi1_mean_0': 'OSDA PMI 1',
         'osda1_pmi2_mean_0': 'OSDA PMI 2',
         'osda1_pmi3_mean_0': 'OSDA PMI 3',
         'osda1_spherocity_index_mean_0': 'OSDA sphericity',
         'osda1_volume_mean_0': 'OSDA volume',
}

# Zeolites according to pore size
s_pore = ['CHA', 'LTA', 'AEI', 'LEV', 'ANA', 'SOD', 'AEN', 'GIS', 'NON', 'AST',
       'ERI', 'DDR', 'RTH', 'MTN', 'MRT', 'RHO', 'ITE', 'AFX', 'ITW', 'DOH',
       'AWO', 'RUT', 'MER', 'ABW', 'PHI', 'LOS', 'EDI', 'SAV', 'SAS', 'MTF',
       'KFI', 'ZON', 'ATN', 'AFN', 'PAU', 'NSI', 'SGT', 'JSW', 'AFV', 'AVL',
       'SWY', 'UFI', 'ETL', 'APC', 'RTE', 'JBW', 'APD', 'CDO', 'MWF', 'MSO',
       'SVV', 'AVE', 'IHW', 'UOZ', 'OWE', 'DFT', 'PWN', 'AWW', 'SFW', 'ATT',
       'LTJ', 'ESV', 'THO', 'LTN', 'EEI', 'POR', 'ACO', 'IRN', 'CAS', 'EAB',
       'BCT', 'ATV', 'SAT']
m_pore = ['MFI', 'MWW', 'TON', 'AEL', 'FER', 'MTT', 'MEL', 'STF', '*MRE', 'ITH',
       'EUO', 'STW', 'SZR', 'AFO', 'STT', 'IMF', 'CSV', 'TUN', 'CGS', 'MFS',
       'NES', 'STI', 'NAT', 'PON', 'LAU', 'HEU', '*UOE', 'UOS', 'PWW', 'ITR',
       '-SVR', '-LIT', 'SFF', 'IFW', 'JST', 'PTY', 'SBN', 'SFG', 'RSN', 'VSV',
       'PWO', 'EWS', 'RRO', 'JRY', 'MVY', 'CGF', 'ETV', 'AHT']
l_pore = ['*BEA', 'AFI', 'MTW', 'FAU', 'MOR', 'BEC', 'EMT', 'IWR', 'ATO', 'MAZ',
       '*STO', 'LTL', 'ISV', 'IFR', 'IWW', 'VET', 'ATS', 'OFF', 'MEI', 'GME',
       'CON', 'BPH', 'SFO', 'MSE', 'CAN', 'UOV', 'IWV', 'SFE', 'ITG', '*-ITN',
       'EZT', 'AFY', 'AFS', 'AFR', 'ASV', 'SFS', 'SAF', 'SSY', 'SBE', 'EON',
       'IWS', 'YFI', 'SAO', 'BOG', 'USI', 'DFO', '*SFV', 'BSV', 'CZP', 'RWY',
       'UWY', 'MOZ', 'GON', 'SOS', 'SSF', 'PUN', 'SOV', 'SOR', 'POS', 'LTF',
       'SOF', 'JSR']
xl_pore = ['UTL', 'IRR', '*-SVY', '*CTH', 'ITT', '-ITV', 'CFI', '-CLO', 'VFI',
       '-IRY', 'IFO', '-IFT', 'SFH', 'SFN', '-IFU', 'ETR', 'DON']

# Small cbus
small_cbus = [
            'lov', 'nat', 'vsv', 'mei', 'sti', 'bea', 'bre', 'jbw', 'mtt', 'afi', 'ats', 'bog', 'cas', 'lau', 'rth', 'bik', 'fer', 
#             'ifw', # not present in df_zeos
            'abw', 'bph', 'mel', 'mtw', 'non', 'ton', 'aww', 'ddr', 'rte', 'can', 'mso', 'gis', 'mtn', 'atn', 'gme', 'obw', 'phi', 'sod',
            'd3r', 'd4r', 'd6r', 'd8r']

# Large cbus
large_cbus = ['rut', 'lev', 'mwf', 'los', 'clo', 'pau', 'ast', 'ave', 'cha', 'doh', 'eab', 'pcr', 'boz', 'lio', 'aft', 'afy', 'lta', 'ltl']
