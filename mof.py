import json
import copy
from pymatgen import Structure
from atomic_parameters import atoms as ap
import numpy as np
import sys
import itertools
import glob
import os


class MofCollection:

    def __init__(self, collection_folder=""):
        self.collection_folder = collection_folder
        self.mof_coll = dict()
        self.load_mofs()

    def load_mofs(self):
        for mof_file in glob.glob(self.collection_folder + "*.cif"):
            mof = MofStructure.from_file(mof_file, primitive=False)
            self.mof_coll[mof.name] = mof

    def analyse_mofs(self, output_folder="output", overwrite=False):
        Helper.make_folder(output_folder)

        for mofname in self.mof_coll:
            mof_folder = '/'.join([output_folder, mofname])
            results_exist = self.check_if_results_exist(mof_folder, mofname)
            if not overwrite and results_exist:
                print("Skiping {}. Results already exist "
                      "and overwrite is False.".format(mofname))
                continue

            mof = self.mof_coll[mofname]
            mof.analyze_metals(output_folder=output_folder,
                               output_level='debug')

    def check_if_results_exist(self, mof_folder, mofname):
        if os.path.isfile(mof_folder+'/'+mofname+'.json'):
            if not os.path.isfile(mof_folder + '/' + 'analysis_running'):
                return True
        return False


class MofStructure(Structure):

    def __init__(self, lattice, species, coords, validate_proximity=False,
                 to_unit_cell=False, coords_are_cartesian=False,
                 site_properties=None):
        super(Structure, self).__init__(lattice, species, coords,
                                        validate_proximity=validate_proximity,
                                        to_unit_cell=to_unit_cell,
                                        coords_are_cartesian=coords_are_cartesian,
                                        site_properties=site_properties)
        self.metal = None
        self.organic = None
        self.metal_coord_spheres = None
        self.all_coord_spheres = None
        self.mof_name = "N/A"

        self.summary = dict()
        self.summary['material_name'] = self.mof_name
        self.summary['problematic'] = False
        self.summary['max_surface_area_frac'] = 0.0
        self.summary['metal_sites'] = []
        self.summary['max_surface_area'] = 0.0
        self.summary['uc_volume'] = self.volume
        self.summary['open_metal_density'] = 0.0

        self.tolerance = dict()
        self.tolerance['plane'] = 25  # 30 # 35
        self.tolerance['plane_5l'] = 25  # 30
        self.tolerance['tetrahedron'] = 10
        self.tolerance['plane_on_metal'] = 12.5


    @classmethod
    def from_file(cls, filename, primitive=False, sort=False, merge_tol=0.0):
        s = super(Structure, cls).from_file(filename, primitive=primitive,
                                            sort=sort, merge_tol=merge_tol)
        s.mof_name = filename.split('/')[-1].split('.')[0]
        s.summary['material_name'] = s.mof_name
        return s

    def set_name(self, name):
        self.summary['material_name'] = name

    def split_structure_to_organic_and_metal(self):
        coords = self.frac_coords
        elements = self.species
        coords_metal = []
        coords_organic = []
        elements_metal = []
        elements_organic = []
        for element, coord in zip(elements, coords):
            if ap.check_if_metal(str(element)):
                elements_metal.append(element)
                coords_metal.append(coord)
            else:
                elements_organic.append(element)
                coords_organic.append(coord)
        self.metal = Structure(self.lattice, elements_metal, coords_metal)
        self.organic = Structure(self.lattice, elements_organic, coords_organic)

    def find_metal_coord_spheres(self):
        centers = self.metal
        coord_spheres = []
        for i, c in enumerate(centers):
            s_ = Structure(self.lattice, [centers.species[i]],
                           [centers.frac_coords[i]])
            c_index = self.match_index(s_)
            coord_spheres.append(self.find_coord_sphere(c_index)[1])
        self.metal_coord_spheres = coord_spheres
        return self.metal_coord_spheres

    def find_all_coord_sphere(self):

        dist_all = self.lattice.get_all_distances(self.frac_coords,
                                                  self.frac_coords)

        self.all_coord_spheres = []
        for a in range(0, len(self)):
            cs = self.find_coord_sphere_using_dist(a, dist_all[a])[0]
            self.all_coord_spheres.append(cs)

    def find_coord_sphere(self, center):
        dist = self.lattice.get_all_distances(self.frac_coords[center],
                                              self.frac_coords)
        res = self.find_coord_sphere_using_dist(center, dist[0])
        coord_sphere, coord_sphere_structure = res
        return coord_sphere, coord_sphere_structure

    def find_coord_sphere_using_dist(self, center, dist):
        if dist[center] > 0.0000001:
            sys.exit('The self distance appears to be non-negative')

        # ligands_structure_new = None
        ligands_structure_new = Structure(self.lattice,
                                          [self.species[center]],
                                          [self.frac_coords[center]])
        max_bond = ap.get_max_bond()
        coord_sphere = [center]
        for i, dis in enumerate(dist):
            if i == center:
                continue
            if dis > max_bond:
                continue
            species_one = str(self.species[center])
            species_two = str(self.species[i])
            tol = ap.get_bond_tolerance(species_one, species_two)
            bond_tol = tol
            if ap.bond_check(species_one, species_two, dis, bond_tol):
                coord_sphere.append(i)
                ligands_structure_new.append(species_two,
                                             self.frac_coords[i])

        ligands = ap.keep_valid_bonds(ligands_structure_new, 0)
        # TODO remove extra atoms from coordination sphere after keeping only valid bonds
        # coord_sphere_kept = []
        # for i, l in enumerate(ligands_structure_new):
        #     for j, lv in enumerate(ligands):
        #         if l == lv:
        #             coord_sphere_kept.append(coord_sphere[i])
        #             continue

        # ligands.insert(0, structure.species[center], structure.frac_coords[center])
        coord_sphere_structure = self.center_around_metal(ligands)
        return coord_sphere, coord_sphere_structure

    def center_around_metal(self, system):
        # return system
        center = system.frac_coords[0]
        tmp1 = []
        tmp2 = []
        tmp1.append(str(system.species[0]))
        tmp2.append([system.frac_coords[0][0], system.frac_coords[0][1],
                     system.frac_coords[0][2]])
        system_centered = MetalSite(system.lattice, tmp1, tmp2)
        for i in range(1, system.num_sites):
            c_i = system.frac_coords[i]
            dist_vector = center - c_i
            dist_vector_r = []
            for j in range(0, 3):
                dist_vector_r.append(round(dist_vector[j]))
            dist_before = np.linalg.norm(
                system.lattice.get_cartesian_coords(center)
                - system.lattice.get_cartesian_coords(c_i))
            c_i_centered = c_i + dist_vector_r
            dist_after = np.linalg.norm(
                system.lattice.get_cartesian_coords(center)
                - system.lattice.
                get_cartesian_coords(c_i_centered))
            if dist_after > dist_before:
                for j in range(0, 3):
                    dist_vector_r[j] = np.rint(dist_vector[j])
                c_i_centered = c_i + dist_vector_r
                if dist_after > dist_before:
                    c_i_centered = c_i
            system_centered.append(system.species[i], c_i_centered)
        return system_centered

    def match_index(self, metal):
        ele = str(metal.species[0])
        f_coords = metal.frac_coords[0]
        dist = self.lattice.get_all_distances(f_coords, self.frac_coords)
        for i, d in enumerate(dist[0]):
            if d < 0.001 and str(self.species[i]) == ele:
                return i

    def analyze_metals(self, output_folder=None, output_level='normal'):
        if output_folder is not None:
            Helper.make_folder(output_folder)
            running_ = output_folder + "/analysis_running"
            open(running_, 'w').close()

        self.split_structure_to_organic_and_metal()
        self.find_metal_coord_spheres()
        self.find_all_coord_sphere()
        oms_cs_list = []  # list of coordination sequences for each open metal found
        cms_cs_list = []  # list of coordination sequences for each closed metal found
        for m, omc in enumerate(self.metal_coord_spheres):
            omc.check_if_open()
            self.summary['metal_sites'].append(omc.site_summary)
            m_index = self.match_index(omc)
            cs = self.find_coordination_sequence(m_index)
            if omc.is_open:
                cs_list = oms_cs_list
            else:
                cs_list = cms_cs_list
            if omc.is_problematic:
                self.summary['problematic'] = True

            m_id, new_site = self.find_metal_id(cs_list, cs)
            if new_site:
                self.summary['metal_sites'][-1]['unique'] = True
                print('New site found')
                cs_list.append(cs)

        open_sites = [s['is_open'] for s in self.summary['metal_sites']]
        unique_sites = [1 for s in self.summary['metal_sites'] if s['unique']]

        self.summary['open_metal_density'] = sum(unique_sites) / self.volume
        self.summary['metal_sites_found'] = any(open_sites)
        if output_folder is not None:
            self.write_results(output_folder, output_level)
            os.remove(running_)
        return self.summary

    def find_metal_id(self, cs_list, cs):
        """Check if a given site is unique based on its coordination sequence"""
        for i, cs_i in enumerate(cs_list):
            if Helper.compare_lists(cs_i, cs):
                return i, False
        return len(cs_list), True

    def find_coordination_sequence(self, center):
        """computes the coordination sequence up to the Nth coordination shell
        as input it takes the MOF as a pymatgen Structure and the index of the
        center metal in the Structure
        """
        import time
        # dist_all = structure.lattice.get_all_distances(structure.frac_coords,
        #                                                structure.frac_coords)

        # all_coord_spheres = []
        # for a in range(0, len(structure)):
        #     all_coord_spheres.append(find_coord_sphere_using_dist(a, structure, dist_all[a])[0])
        # The shell_list is a set with the index of each atom and its unit
        # cell index realtive to a cetral unit cell
        shell_list = {(center, (0, 0, 0))}
        shell_list_prev = set([])
        all_shells = set(shell_list)
        n_shells = 6
        cs = []
        ele = [(str(self.species[center]))]
        coords = [[self.frac_coords[center][0],
                   self.frac_coords[center][1],
                   self.frac_coords[center][2]]]
       # coordination_structure = (Molecule(ele, coords))
        coord_sphere_time = 0.0
        count_total = 0
        for n in range(0, n_shells):
            c_set = set([])
            for a_uc in shell_list:
                a = a_uc[0]
                lattice = a_uc[1]
                t0 = time.time()
                #TODO make finding coordination sphere faster
                # coord_sphere = find_coord_sphere_using_dist(a, structure,
                #                                             dist_all[a])[0]
                # print(a)
                # print(coord_sphere)
                # print(find_coord_sphere_using_dist(a, structure,
                #                                             dist_all[a])[1])
                # input()
                coord_sphere = self.all_coord_spheres[a]
                count_total += 1
                t1 = time.time()
                coord_sphere_time += t1-t0
                coord_sphere_with_uc = []
                for c in coord_sphere:
                    # dist = structure.lattice.\
                    #     get_all_distance_and_image(structure.frac_coords[a],
                    #                                structure.frac_coords[c])
                    #
                    # uc = tuple(l-nl for l, nl in zip(lattice, min(dist)[1]))

                    diff = self.frac_coords[a] - self.frac_coords[c]
                    new_lat_i = [round(d, 0) for d in diff]
                    uc = tuple(l-nl for l, nl in zip(lattice, new_lat_i))

                    # check = [(d-md) < 1e-5 for d, md in zip(new_lat_i, min(dist)[1])]
                    # if not any(check):
                    #     print(min(dist)[1])
                    #     print(new_lat_i)
                    #     print(check)
                    #     input()

                    coord_sphere_with_uc.append((c, uc))
                coord_sphere_with_uc = tuple(coord_sphere_with_uc)
                c_set = c_set.union(set(coord_sphere_with_uc))
            for a in shell_list_prev:
                c_set.discard(a)
            for a in shell_list:
                c_set.discard(a)
            # for i_uc in c_set:
            #     i = i_uc[0]
            #     uc = i_uc[1]
            #     ele = ap().elements[n+3]
            #     coords = [structure.frac_coords[i][0] - uc[0],
            #               structure.frac_coords[i][1] - uc[1],
            #               structure.frac_coords[i][2] - uc[2]]
            #     coords = structure.lattice.get_cartesian_coords(coords)
            #     coordination_structure.append(ele, coords, validate_proximity=False)

            cs.append(len(c_set))
            all_shells = all_shells.union(c_set)
            shell_list_prev = shell_list
            shell_list = c_set
        # coordination_structure = center_around_metal(coordination_structure)
        # write_xyz_file('_cs.xyz', coordination_structure)
        # print('coordination sphere time:', coord_sphere_time, count_total, len(structure))
        # print('Coordination Sphere: ', cs)
        # input()
        return cs

    def write_results(self, output_folder, output_level):
        Helper.make_folder(output_folder)
        for index, mcs in enumerate(self.metal_coord_spheres):
            mcs.write_xyz_file(output_folder, index)

        output_fname = output_folder + '/' + self.mof_name + '_metal.cif'
        self.metal.to(filename=output_fname)
        output_fname = output_folder + '/' + self.mof_name + '_organic.cif'
        self.organic.to(filename=output_fname)

        json_file_out = output_folder + '/' + self.mof_name + '.json'
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
                 site_properties=None):
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

        self.tolerance = dict()
        self.tolerance['plane'] = 25  # 30 # 35
        self.tolerance['plane_5l'] = 25  # 30
        self.tolerance['tetrahedron'] = 10
        self.tolerance['plane_on_metal'] = 12.5

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
            # if num_l == -4:
            #     if not (abs(dihedral - crit[test_type]) < tol[test_type]
            #             or abs(dihedral - crit[test_type] + 180) < tol[
            #             test_type]):
            #         self.site_summary["type"] = type
            #         self.site_summary["is_open"] = True
            # else:
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
                        dihedral = abs(self.get_dihedral(m_i, indeces[i],
                                                         indeces[j],
                                                         indeces[k]))
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
        # if not (self.check_positive(dihedrals_other) and
        #         self.check_negative(dihedrals_other)):
            # return True
        return False

    # @staticmethod
    # def check_positive(n_list):
    #     return any([n > 0 for n in n_list])
        # for n in n_list:
        #     if n > 0:
        #         return True

    # @staticmethod
    # def check_negative(n_list):
    #     return any([n < 0 for n in n_list])
        # for n in n_list:
        #     if n < 0:
        #         return True

    def check_metal_dihedrals(self):
        crit = 180
        tol = self.tolerance['plane']  # 30
        oms_type = ','.join([str(self.species[0]), 'metal_plane'])
        oms_test = dict()
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
        import pymatgen.io.xyz as xyzio
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