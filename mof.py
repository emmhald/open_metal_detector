
import json
import copy
from pymatgen import Structure
from atomic_parameters import atoms as ap
import numpy as np
import sys
import itertools
import os
import shutil


class MofStructure(Structure):

    def __init__(self, lattice, species, coords, charge=None, validate_proximity=False,
                 to_unit_cell=False, coords_are_cartesian=False,
                 site_properties=None):
        super(Structure, self).__init__(lattice, species, coords,
                                        charge=charge,
                                        validate_proximity=validate_proximity,
                                        to_unit_cell=to_unit_cell,
                                        coords_are_cartesian=coords_are_cartesian,
                                        site_properties=site_properties)
        self.metal = None
        self.organic = None
        self.metal_coord_spheres = None
        self.all_coord_spheres_indices = None
        self.name = "N/A"
        self.species_str = [str(s) for s in self.species]

        self.summary = dict()
        self.summary['material_name'] = self.name
        self.summary['problematic'] = False
        self.summary['max_surface_area_frac'] = 0.0
        self.summary['metal_sites'] = []
        self.summary['max_surface_area'] = 0.0
        self.summary['uc_volume'] = self.volume
        self.summary['open_metal_density'] = 0.0
        #We need to set to have only one of each metal atom, and then convert to
        #a list to store as JSON (JSON does not support sets)
        metal_set = set([s for s in self.species_str if ap.check_if_metal(s)])
        self.summary['metals'] = list(metal_set)

        self.tolerance = dict()
        self.tolerance['plane'] = 25  # 30 # 35
        self.tolerance['plane_5l'] = 25  # 30
        self.tolerance['tetrahedron'] = 10
        self.tolerance['plane_on_metal'] = 12.5

        self.metal_indexes = []
        self.all_distances = None
        self.max_bond = ap().get_max_bond()

    @classmethod
    def from_file(cls, filename, primitive=False, sort=False, merge_tol=0.0):
        s = Structure.from_file(filename, primitive=primitive, sort=sort,
                                merge_tol=merge_tol)
        s_mof = cls(s.lattice, s.species, s.frac_coords)
        s_mof.name = os.path.splitext(os.path.basename(filename))[0]
        s_mof.summary['material_name'] = s_mof.name
        return s_mof

    def set_name(self, name):
        self.summary['material_name'] = name

    def split_structure_to_organic_and_metal(self):
        self.metal = Structure(self.lattice, [], [])
        self.organic = Structure(self.lattice, [], [])
        i = 0
        for s, fc in zip(self.species, self.frac_coords):
            if ap.check_if_metal(str(s)):
                self.metal.append(s, fc)
                self.metal_indexes.append(i)
            else:
                self.organic.append(s, fc)
            i += 1

    def find_metal_coord_spheres(self):
        centers = self.metal
        coord_spheres = []
        for i, c in enumerate(centers):
            c_index = self.metal_indexes[i]
            coord_spheres.append(self.find_coord_sphere(c_index))
        self.metal_coord_spheres = coord_spheres

    def compute_distances(self):
        self.all_distances = self.lattice.get_all_distances(self.frac_coords,
                                                            self.frac_coords)

    def find_all_coord_sphere_indices(self):
        dist_all = self.all_distances
        self.all_coord_spheres_indices = []
        for a in range(0, len(self)):
            cs_i = self.find_coord_sphere_indices(a, dist_all[a])
            self.all_coord_spheres_indices.append(cs_i)

    def find_coord_sphere_indices(self, center, dist):
        if dist[center] > 0.0000001:
            sys.exit('The self distance appears to be non-zero')

        c_sphere_indices = [center]
        s_center = self.species_str[center]
        valid_dists = [(i, dis) for i, dis in enumerate(dist) if i != center
                       and dis <= self.max_bond]
        for i, dis in valid_dists:
            species_two = self.species_str[i]
            bond_tol = ap.get_bond_tolerance(s_center, species_two)
            if ap.bond_check(s_center, species_two, dis, bond_tol):
                c_sphere_indices.append(i)
        # TODO remove extra atoms from coordination sphere after keeping only valid bonds
        return c_sphere_indices

    def find_coord_sphere(self, center):
        dist = self.all_distances[center]
        if dist[center] > 0.0000001:
            sys.exit('The self distance appears to be non-zero')

        c_sphere = MetalSite(self.lattice, [self.species[center]],
                                           [self.frac_coords[center]],
                                           tolerance=self.tolerance)

        cs_i = self.find_coord_sphere_indices(center, dist)
        for i in cs_i[1:]:
            c_sphere.append(self.species_str[i], self.frac_coords[i])

        c_sphere = ap.keep_valid_bonds(c_sphere, 0)

        # ligands.insert(0, structure.species[center], structure.frac_coords[center])
        c_sphere = self.center_around_metal(c_sphere)
        return c_sphere

    def center_around_metal(self, coordination_sphere):
        gc = coordination_sphere.lattice.get_cartesian_coords
        center = coordination_sphere.frac_coords[0]
        center_cart_coords = gc(center)
        for i in range(1, coordination_sphere.num_sites):
            c_i = coordination_sphere.frac_coords[i]
            dist_vector = center - c_i
            dist_vector_r = []
            for j in range(0, 3):
                dist_vector_r.append(round(dist_vector[j]))
            dist_before = np.linalg.norm(center_cart_coords - gc(c_i))
            c_i_centered = c_i + dist_vector_r
            dist_after = np.linalg.norm(center_cart_coords - gc(c_i_centered))
            if dist_after > dist_before:
                for j in range(0, 3):
                    dist_vector_r[j] = np.rint(dist_vector[j])
                c_i_centered = c_i + dist_vector_r
                if dist_after > dist_before:
                    c_i_centered = c_i
            coordination_sphere.replace(i, coordination_sphere.species[i],
                                        c_i_centered)
        return coordination_sphere

    def analyze_metals(self, output_folder=None, output_level='normal'):
        if output_folder is not None:
            Helper.make_folder(output_folder)
            running_ = output_folder + "/analysis_running"
            open(running_, 'w').close()

        self.compute_distances()
        self.split_structure_to_organic_and_metal()
        self.find_metal_coord_spheres()
        self.find_all_coord_sphere_indices()
        oms_cs_list = []  # list of coordination sequences for each open metal found
        cms_cs_list = []  # list of coordination sequences for each closed metal found
        # cs_list = []
        for m, omc in enumerate(self.metal_coord_spheres):
            m_index = self.metal_indexes[m]

            omc.check_if_open()
            self.summary['metal_sites'].append(omc.site_summary)
            if omc.is_problematic:
                self.summary['problematic'] = True

            cs = self.find_coordination_sequence(m_index)

            if omc.is_open:
                cs_list = oms_cs_list
            else:
                cs_list = cms_cs_list

            if self.check_if_new_site(cs_list, cs):
                self.summary['metal_sites'][-1]['unique'] = True
                cs_list.append(cs)

        unique_sites = [1 for s in self.summary['metal_sites'] if s['unique']]
        open_sites = [s['is_open'] for s in self.summary['metal_sites']]

        self.summary['open_metal_density'] = sum(unique_sites) / self.volume
        self.summary['metal_sites_found'] = any(open_sites)
        if output_folder is not None:
            self.write_results(output_folder, output_level)
            os.remove(running_)
        return self.summary

    def check_if_new_site(self, cs_list, cs):
        """Check if a given site is unique based on its coordination sequence"""
        for i, cs_i in enumerate(cs_list):
            if Helper.compare_lists(cs_i, cs):
                return False  #i,
        return True  #len(cs_list),

    def find_coordination_sequence(self, center):
        """computes the coordination sequence up to the Nth coordination shell
        as input it takes the MOF as a pymatgen Structure and the index of the
        center metal in the Structure
        """

        shell_list = {(center, (0, 0, 0))}
        shell_list_prev = set([])
        all_shells = set(shell_list)
        n_shells = 6
        cs = []
        count_total = 0
        for n in range(0, n_shells):
            c_set = set([])
            for a_uc in shell_list:
                a = a_uc[0]
                lattice = a_uc[1]
                coord_sphere = self.all_coord_spheres_indices[a]
                count_total += 1
                coord_sphere_with_uc = []
                for c in coord_sphere:
                    diff = self.frac_coords[a] - self.frac_coords[c]
                    new_lat_i = [round(d, 0) for d in diff]
                    uc = tuple(l-nl for l, nl in zip(lattice, new_lat_i))
                    coord_sphere_with_uc.append((c, uc))
                coord_sphere_with_uc = tuple(coord_sphere_with_uc)
                c_set = c_set.union(set(coord_sphere_with_uc))
            for a in shell_list_prev:
                c_set.discard(a)
            for a in shell_list:
                c_set.discard(a)

            cs.append(len(c_set))
            all_shells = all_shells.union(c_set)
            shell_list_prev = shell_list
            shell_list = c_set

        return cs

    def write_results(self, output_folder, output_level):
        Helper.make_folder(output_folder)
        for index, mcs in enumerate(self.metal_coord_spheres):
            mcs.write_xyz_file(output_folder, index)

        output_fname = output_folder + '/' + self.name + '_metal.cif'
        if self.metal:
            self.metal.to(filename=output_fname)
        output_fname = output_folder + '/' + self.name + '_organic.cif'
        self.organic.to(filename=output_fname)

        json_file_out = output_folder + '/' + self.name + '.json'
        summary = copy.deepcopy(self.summary)
        if output_level == 'normal':
            for ms in summary["metal_sites"]:
                ms.pop('all_dihedrals', None)
                ms.pop('min_dihedral', None)
        with open(json_file_out, 'w') as outfile:
            json.dump(summary, outfile, indent=3)


class MetalSite(Structure):

    def __init__(self, lattice, species, coords, validate_proximity=False,
                 to_unit_cell=False, coords_are_cartesian=False,
                 site_properties=None, tolerance="default"):
        super(Structure, self).__init__(lattice, species, coords,
                                        validate_proximity=validate_proximity,
                                        to_unit_cell=to_unit_cell,
                                        coords_are_cartesian=coords_are_cartesian,
                                        site_properties=site_properties)
        self._sites = list(self._sites)

        self.site_summary = dict()
        self.site_summary["metal"] = str(self.species[0])
        self.site_summary["type"] = "closed"
        self.site_summary["is_open"] = False
        self.site_summary["unique"] = False
        self.site_summary["problematic"] = False
        if tolerance == 'default':
            self.tolerance = dict()
            self.tolerance['plane'] = 25  # 30 # 35
            self.tolerance['plane_5l'] = 25  # 30
            self.tolerance['tetrahedron'] = 10
            self.tolerance['plane_on_metal'] = 12.5
        else:
            self.tolerance = tolerance

        self.all_dihedrals = []
        self.all_indeces = []
        self.all_dihedrals_m = []
        self.all_indeces_m = []
        self.is_open = False
        self.is_problematic = False

    def check_if_open(self):
        self.site_summary["number_of_linkers"] = self.num_sites - 1
        self.get_t_factor()

        num_l = self.num_sites - 1
        min_coordination = 3
        if ap.is_lanthanide_or_actinide(str(self.species[0])):
            min_coordination = 5

        if num_l < min_coordination:
            self.site_summary["problematic"] = True

        if num_l <= 3:  # min_cordination:
            self.site_summary["is_open"] = True
            self.site_summary["type"] = '3_or_less'
            self.site_summary["min_dihedral"] = 0.0
            self.site_summary["all_dihedrals"] = 0.0
        else:
            self.check_non_metal_dihedrals()

        if self.site_summary['is_open']:
            self.is_open = True

        if self.site_summary['problematic']:
            self.is_problematic = True

    def get_t_factor(self):

        index_range = range(1, self.num_sites)
        all_angles = []
        for i in itertools.combinations(index_range, 2):
            angle = self.get_angle(i[0], 0, i[1])
            all_angles.append([angle, i[0], i[1]])

        all_angles.sort(key=lambda x: x[0])
        # beta is the largest angle and alpha is the second largest angle
        # in the coordination sphere; using the same convention as Yang et al.
        # DOI: 10.1039/b617136b
        if self.num_sites > 3:
            beta = all_angles[-1][0]
            alpha = all_angles[-2][0]

        if self.num_sites - 1 == 6:
            max_indeces_all = all_angles[-1][1:3]
            l3_l4_angles = [x for x in all_angles if
                            x[1] not in max_indeces_all and
                            x[2] not in max_indeces_all]
            max_indeces_all_3_4 = max(l3_l4_angles, key=lambda x: x[0])[1:3]
            l5_l6_angles = [x for x in l3_l4_angles
                            if x[1] not in max_indeces_all_3_4 and
                            x[2] not in max_indeces_all_3_4]
            gamma = max(l5_l6_angles, key=lambda x: x[0])[0]
            tau = self.get_t6_factor(gamma)
        elif self.num_sites - 1 == 5:
            tau = self.get_t5_factor(alpha, beta)
        elif self.num_sites - 1 == 4:
            tau = self.get_t4_factor(alpha, beta)
        else:
            tau = -1
        self.site_summary['t_factor'] = tau

    @staticmethod
    def get_t4_factor(a, b):
        return (360 - (a + b)) / 141.0

    @staticmethod
    def get_t5_factor(a, b):
        return (b - a) / 60.0

    @staticmethod
    def get_t6_factor(c):
        return c / 180

    def check_non_metal_dihedrals(self):
        #TODO revise this method. Simplify calculation and bookeeping
        num_l = self.num_sites - 1
        crit = dict()
        tol = dict()

        oms_test = dict()
        oms_test['plane'] = False
        oms_test['same_side'] = False
        oms_test['metal_plane'] = False
        oms_test['3_or_less'] = False
        oms_test['non_TD'] = False

        crit['plane'] = 180
        tol['plane'] = self.tolerance['plane']  # 35
        crit['plane_5l'] = 180
        tol['plane_5l'] = self.tolerance['plane_5l']  # 30
        crit['tetrahedron'] = 70.528779  # 70
        tol['tetrahedron'] = self.tolerance['tetrahedron']  # 10

        if num_l == 4:
            test_type = 'plane'  # 'tetrahedron'
            om_type = 'non_TD'
        elif num_l == 5:
            test_type = 'plane_5l'
            om_type = 'plane_5l'
        elif num_l > 5:
            test_type = 'plane'
            om_type = 'same_side'
        oms_type = ','.join([test_type, om_type])

        self.obtain_dihedrals()
        self.site_summary["min_dihedral"] = min(self.all_dihedrals)
        self.site_summary["all_dihedrals"] = self.all_dihedrals

        for dihedral, indeces in zip(self.all_dihedrals, self.all_indeces):
            [i, j, k, l] = indeces
            test1 = abs(dihedral - crit[test_type])
            test2 = abs(dihedral - crit[test_type] + 180)
            if test1 < tol[test_type] or test2 < tol[test_type]:
                if num_l <= 5:
                    self.site_summary["type"] = oms_type
                    self.site_summary["is_open"] = True
                elif num_l > 5:
                    if self.check_if_plane_on_metal(0, [i, j, k, l]):
                        other_i = self.find_other_indeces([0, i, j, k, l])
                        if self.check_if_atoms_on_same_side(other_i, j, k, l):
                            self.site_summary["type"] = oms_type
                            self.site_summary["is_open"] = True

        if (num_l >= 4) and not self.site_summary["is_open"]:
            self.check_metal_dihedrals()

    def obtain_dihedrals(self):
        indices_1 = range(1, self.num_sites)
        indices_2 = range(1, self.num_sites)
        for i, l in itertools.combinations(indices_1, 2):
            for j, k in itertools.combinations(indices_2, 2):
                if len({i, j, k, l}) == 4:
                    dihedral = abs(self.get_dihedral(i, j, k, l))
                    self.all_dihedrals.append(dihedral)
                    self.all_indeces.append([i, j, k, l])

    def check_if_plane_on_metal(self, m_i, indeces):
        crit = 180
        # Set to 12.5 so that ferocene type coordination spheres are detected
        # correctly. eg. BUCROH
        tol = self.tolerance['plane_on_metal']  # 12.5
        # tol = 25.0
        for i in range(1, len(indeces)):
            for j in range(1, len(indeces)):
                for k in range(1, len(indeces)):
                    if i == j or i == k or j == k:
                        pass
                    else:
                        d = self.get_dihedral(m_i, indeces[i], indeces[j],
                                              indeces[k])
                        dihedral = abs(d)
                        if abs(dihedral - crit) < tol or abs(
                                                dihedral - crit + 180) < tol:
                            return True
        return False

    def find_other_indeces(self, indeces):
        other_indeces = []
        for index in range(0, self.num_sites):
            if index not in indeces:
                other_indeces.append(index)
        return other_indeces

    def check_if_atoms_on_same_side(self, other_indeces, j, k, l):
        dihedrals_other = []
        for o_i in other_indeces:
            dihedrals_other.append(self.get_dihedral(j, k, l, o_i))
        dihedrals_p = [d > 0.0 for d in dihedrals_other]
        dihedrals_n = [d < 0.0 for d in dihedrals_other]
        if all(dihedrals_p) or all(dihedrals_n):
            return True
        return False

    def check_metal_dihedrals(self):
        crit = 180
        tol = self.tolerance['plane']  # 30
        oms_type = ','.join([str(self.species[0]), 'metal_plane'])
        self.obtain_metal_dihedrals()
        number_of_planes = 0
        for dihedral, indices in zip(self.all_dihedrals_m, self.all_indeces_m):
            [i, j, k, l] = indices
            if abs(dihedral - crit) < tol or abs(dihedral - crit + 180) < tol:
                number_of_planes += 1
                all_indices = self.find_other_indeces([0, j, k, l])
                if self.check_if_atoms_on_same_side(all_indices, j, k, l):
                    self.site_summary["type"] = oms_type
                    self.site_summary["is_open"] = True

    def obtain_metal_dihedrals(self):
        indices_1 = range(1, self.num_sites)
        i = 0
        for l in indices_1:
            for j, k in itertools.permutations(indices_1, 2):
                if len({i, j, k, l}) == 4:
                    dihedral = abs(self.get_dihedral(i, j, k, l))
                    self.all_dihedrals_m.append(dihedral)
                    self.all_indeces_m.append([i, j, k, l])

    def write_xyz_file(self, output_folder, index):
        Helper.make_folder(output_folder)
        output_fname = output_folder
        output_fname += '/first_coordination_sphere'+str(index)+'.cif'
        self.to(filename=output_fname)


class Helper:

    @classmethod
    def make_folder(cls, output_folder):
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

    @classmethod
    def compare_lists(cls, l1, l2):
        if len(l1) != len(l2):
            return False
        return all([i == j for i, j in zip(l1, l2)])

    @classmethod
    def copy_folder(cls, dest, src):
        if not os.path.exists(dest):
            os.makedirs(dest)
        s = src.split('/')[-1]
        d = os.path.join(dest, s)
        if not os.path.exists(d):
            shutil.copytree(src, d)
