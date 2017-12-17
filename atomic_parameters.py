# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %(Emmanuel Haldoupis)s
"""

class atoms:
    atomic_number={
    'H':1,
    'D':2,
    'He':3,
    'Li':4,
    'Be':5,
    'B':6,
    'C':7,
    'N':8,
    'O':9,
    'F':10,
    'Ne':11,
    'Na':12,
    'Mg':13,
    'Al':14,
    'Si':15,
    'P':16,
    'S':17,
    'Cl':18,
    'Ar':19,
    'K':20,
    'Ca':21,
    'Sc':22,
    'Ti':23,
    'V':24,
    'Cr':25,
    'Mn':26,
    'Fe':27,
    'Co':28,
    'Ni':29,
    'Cu':30,
    'Zn':31,
    'Ga':32,
    'Ge':33,
    'As':34,
    'Se':35,
    'Br':36,
    'Kr':37,
    'Rb':38,
    'Sr':39,
    'Y':40,
    'Zr':41,
    'Nb':42,
    'Mo':43,
    'Tc':44,
    'Ru':45,
    'Rh':46,
    'Pd':47,
    'Ag':48,
    'Cd':49,
    'In':50,
    'Sn':51,
    'Sb':52,
    'Te':53,
    'I':54,
    'Xe':55,
    'Cs':56,
    'Ba':57,
    'La':58,
    'Ce':59,
    'Pr':60,
    'Nd':61,
    'Pm':62,
    'Sm':63,
    'Eu':64,
    'Gd':65,
    'Tb':66,
    'Dy':67,
    'Ho':68,
    'Er':69,
    'Tm':70,
    'Yb':71,
    'Lu':72,
    'Hf':73,
    'Ta':74,
    'W':75,
    'Re':76,
    'Os':77,
    'Ir':78,
    'Pt':79,
    'Au':80,
    'Hg':81,
    'Tl':82,
    'Pb':83,
    'Bi':84,
    'Po':85,
    'At':86,
    'Rn':87,
    'Fr':88,
    'Ra':89,
    'Ac':90,
    'Th':91,
    'Pa':92,
    'U':93,
    'Np':94,
    'Pu':95,
    'Am':96,
    'Cm':97
    }

    elements=['H', 'H', 'D', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
              'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc',
              'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge',
              'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc',
              'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
              'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb',
              'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os',
              'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr',
              'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm',
    ]
    # Covalent radii taken from DOI: Covalent radii revisited
    # Beatriz Cordero,   Verónica Gómez,   Ana E. Platero-Prats,
    # Marc Revés,   Jorge Echeverría,  Eduard Cremades, Flavia Barragána
    # and Santiago Alvarez Dalton Trans., 2008, 2832-2838
    # DOI: 10.1039/B801115J
    co = {
        'H': 0.31,
        'D': 0.31,
        'He': 0.28,
        'Li': 1.28,
        'Be': 0.96,
        'B': 0.84,
        'C': 0.73,
        'N': 0.71,
        'O': 0.66, # 0.8,
        'F': 0.57,
        'Ne': 0.58,
        'Na': 1.66,
        'Mg': 1.41,
        'Al': 1.21,
        'Si': 1.11,
        'P': 1.07,
        'S': 1.05,
        'Cl': 1.02,
        'Ar': 1.06,
        'K': 2.03,
        'Ca': 1.76,
        'Sc': 1.7,
        'Ti': 1.6,
        'V': 1.53,
        'Cr': 1.39,
        'Mn': 1.5,
        'Fe': 1.42,
        'Co': 1.38,
        'Ni': 1.24,
        'Cu': 1.32,
        'Zn': 1.22,
        'Ga': 1.22,
        'Ge': 1.2,
        'As': 1.19,
        'Se': 1.2,
        'Br': 1.2,
        'Kr': 1.16,
        'Rb': 2.2,
        'Sr': 1.95,
        'Y': 1.9,
        'Zr': 1.75,
        'Nb': 1.64,
        'Mo': 1.54,
        'Tc': 1.47,
        'Ru': 1.46,
        'Rh': 1.42,
        'Pd': 1.39,
        'Ag': 1.45,
        'Cd': 1.44,
        'In': 1.42,
        'Sn': 1.39,
        'Sb': 1.39,
        'Te': 1.38,
        'I': 1.39,
        'Xe': 1.4,
        'Cs': 2.44,
        'Ba': 2.15, #1.80
        'La': 2.07,
        'Ce': 2.04,
        'Pr': 2.03,
        'Nd': 2.01,
        'Pm': 1.99,
        'Sm': 1.98,
        'Eu': 1.98,
        'Gd': 1.96,
        'Tb': 1.94,
        'Dy': 1.92,
        'Ho': 1.92,
        'Er': 1.89,
        'Tm': 1.9,
        'Yb': 1.87,
        'Lu': 1.87,
        'Hf': 1.75,
        'Ta': 1.7,
        'W': 1.62,
        'Re': 1.51,
        'Os': 1.44,
        'Ir': 1.41,
        'Pt': 1.36,
        'Au': 1.36,
        'Hg': 1.32,
        'Tl': 1.45,
        'Pb': 1.46,
        'Bi': 1.48,
        'Po': 1.4,
        'At': 1.5,
        'Rn': 1.5,
        'Fr': 2.6,
        'Ra': 2.21,
        'Ac': 2.15,
        'Th': 2.06,
        'Pa': 2,
        'U': 1.96,
        'Np': 1.9,
        'Pu': 1.87,
        'Am': 1.8,
        'Cm': 1.69
    }

    def __init__(self):
        self.co = dict()
        self.atomic_number = dict()

        self.metals = list()
        self.define_metals()

    @classmethod
    def get_covelent_radius(cls, ele):
        return cls.co[ele]

    @classmethod
    def get_vdf_radius(cls, ele):
        #return 2.0 #get_covelent_radius(ele)+0.76
        return cls.get_covelent_radius(ele)+0.76

    @classmethod
    def check_if_metal(cls, ele):
        list_of_non_metals = ['H', 'D', 'B', 'C', 'N', 'O', 'F', 'Si', 'P', 'S',
                              'Cl', 'Br', 'I']
        if ele not in list_of_non_metals:
            return True
        return False

    @classmethod
    def check_if_first_row(cls, ele):
        list_of_first_row = ['Li', 'Na', 'K', 'Rb', 'Cs', 'Fr']
        if ele in list_of_first_row:
            return True
        return False


    @classmethod
    def is_heavy_metal(cls, ele): #\m/
        if cls.get_covelent_radius(ele) > 1.95:
             return True
        else:
            return False

    @classmethod
    def get_bond_tolerance(cls, ele1, ele2):
        if cls.check_if_heavy_metal_bond(ele1, ele2):
            return 0.2
        else:
            return 0.5  #0.4

    def define_metals(self):
        for e in self.elements:
            if self.check_if_metal(e):
                self.metals.append(e)

    @classmethod
    def is_lanthanide_or_actinide(cls, ele):
        if cls.is_lanthanide(ele) or cls.is_actinide(ele):
            return True
        else:
            return False

    @classmethod
    def is_lanthanide(cls, ele):
        if cls.atomic_number[ele] > 56 and cls.atomic_number[ele] < 72:
            return True
        else:
            return False

    @classmethod
    def is_actinide(cls, ele):
        if cls.atomic_number[ele] > 88 and cls.atomic_number[ele] <97:
            return True
        else:
            return False

    @classmethod
    def is_group1(cls, ele):
        list_of_non_metals = ['Li', 'Na', 'K', 'Rb', 'Cs', 'Fr']
        if ele in list_of_non_metals:
            return True
        else:
            return False

    @classmethod
    def is_d4_or_less(cls, ele):
        return cls.check_atomic_number(cls, ele, 54)

    @classmethod
    def check_atomic_number(cls, ele, max_atomic_number):
        if cls.atomic_number[ele] < max_atomic_number:
            return True
        else:
            return False

    @classmethod
    def check_if_heavy_metal_bond(cls, ele1, ele2):
        if cls.is_heavy_metal(ele1) or cls.is_heavy_metal(ele2):
            return True
        else:
            return False

    @classmethod
    def check_if_valid_bonds(cls, ligands, bond_tol, increase):
        if increase > 0.05:
            for l,fc in enumerate(ligands.frac_coords):
                if l == 0:
                    continue
                dist = ligands.lattice.get_all_distances(fc, ligands.frac_coords)
                count_bonds=0
                for i, dis in enumerate(dist[0]):
                    if i != l:
                        if not (cls.check_if_metal(str(ligands.species[l])) and
                                cls.check_if_metal(str(ligands.species[i]))):
                            if not (str(ligands.species[l]) == 'C' and
                                    str(ligands.species[i]) == 'C'):
                                if cls.bond_check(str(ligands.species[l]),
                                                  str(ligands.species[i]),
                                                  dis, bond_tol):
                                    count_bonds += 1
                if count_bonds > 0:
                    return False
        return True

    @classmethod
    def keep_valid_bonds(cls, ligands, center):
        if len(ligands) == 0:
            return ligands
        for l, fc in enumerate(ligands.frac_coords):
            if l == center:
                continue
            dist = ligands.lattice.get_all_distances(fc, ligands.frac_coords)
            for i, dis in enumerate(dist[0]):
                if i == l or i == center:
                    continue
                if not cls.is_valid(ligands, l, i, dis):
                    index_to_remove = cls.index_closer_to_center(ligands,
                                                                 l, i,
                                                                 center)
                    ligands.remove_sites([index_to_remove])
                    return cls.keep_valid_bonds(ligands, center)
        return ligands

    @classmethod
    def index_closer_to_center(cls, ligands, l, i, center):
        lat = ligands.lattice
        dist_l = lat.get_all_distances(ligands.frac_coords[l],
                                       ligands.frac_coords[center])[center]
        dist_i = lat.get_all_distances(ligands.frac_coords[i],
                                       ligands.frac_coords[center])[center]
        if dist_i > dist_l:
            return i
        if dist_i < dist_l:
            return l
        if dist_i == dist_l:
            return i

    @classmethod
    def is_valid(cls, ligands, l, i, dis):
        center = str(ligands.species[0])
        species_one = str(ligands.species[l])
        species_two = str(ligands.species[i])
        tol = cls.get_bond_tolerance(species_one, species_two)
        if cls.check_if_metal(species_one) and cls.check_if_metal(species_two):
            if species_one == center and species_two == center:
                return True
        if species_one == 'C' and species_two == 'C':
            return True
        if cls.bond_check(species_one, species_two, dis, tol):
            return False
        return True


    @classmethod
    def get_sum_of_cov_radii(cls, ele1, ele2):
        return cls.get_covelent_radius(ele1) + cls.get_covelent_radius(ele2)

    @classmethod
    def bond_check(cls, ele1, ele2, dist, bond_tol=None):
        #bond = dist #- get_sum_of_cov_radii(ele1,ele2)
        if bond_tol is None:
            bond_tol = cls.get_bond_tolerance(ele1, ele2)
        up_bound = cls.get_sum_of_cov_radii(ele1, ele2) + bond_tol
        low_bound = cls.get_sum_of_cov_radii(ele1, ele2) - bond_tol
        #if bond < bond_tol and abs(dist) > 0:
        if dist < up_bound : #and dist > low_bound: # and abs(dist) > 0:
            return True
        else:
            return False

    @classmethod
    def get_max_bond(cls):
        """Get the maximum possible distance
        Get the maximum possible covalent radius, multilpy by 1.1 to be safe
        and compute the maximum bond distance as twice that plus the maximum
        tolerance
        """
        max_co = max([cls.co[k] for k in cls.co])*1.1
        return max_co*2.0 + 0.4

