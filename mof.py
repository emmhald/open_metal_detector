
from pymatgen import Structure
from atomic_parameters import atoms as ap
import numpy as np
import sys


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

    def find_all_coord_spheres(self):
        centers = self.metal
        coord_spheres = []
        for i, c in enumerate(centers):
            s_ = Structure(self.lattice, [centers.species[i]],
                           [centers.frac_coords[i]])
            c_index = self.match_index(s_)
            coord_spheres.append(self.find_coord_sphere(c_index)[1])
        return coord_spheres

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
        system_centered = Structure(system.lattice, tmp1, tmp2)
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
