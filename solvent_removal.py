"""
A python code to remove free and coordinated solvents.

* Features
- Remove free and coordinated solvents
- Keeps the fragments with largest number of atoms in the structure

* Usage:  $ python solvent_remove.py *input_file*

* Required Python packages:
- Numpy
- PyMatgen
- Atomic Simulation Environment (ASE)
- SciPy
- Python 2.7

* Code workflow

        input_file -> ASE_NeighborList ----> SparseMatrix ----> ClusterCalc -----> Output
                       |                                                      /
                       |                                                     /
                       -> CustomMatrix -> Modify Adj. Matrix -> ClusterCalc
"""
import sys, os, json, re, warnings, itertools, collections, timeit
warnings.filterwarnings('ignore')

# Numpy Library
import numpy as np

# PyMatgen Library
import pymatgen
from pymatgen.io.smartio import write_structure, read_structure
from pymatgen.io.cifio import CifWriter
from pymatgen.io.cifio import CifParser
from pymatgen.io.aseio import AseAtomsAdaptor
from pymatgen.transformations.standard_transformations import PrimitiveCellTransformation as pct

# Atomic Simulation Environment Library
from ase.io import read as ase_read
from ase.io import write as ase_write
from ase.calculators.neighborlist import NeighborList

# SciPy Library
from scipy.sparse import coo_matrix as sparse_matrix
from scipy.sparse.csgraph import connected_components

#Dictionary of atomic radii
cambridge_radii = {
'Ag': 1.45,
'Ac': 1.950,
'Al': 1.210,
'Am': 1.8,
'Sb': 1.390,
'Ar': 1.060,
'As': 1.190,
'At': 1.50,
'Ba': 2.150,
'Be': 0.960,
'Bi': 1.480,
'B': 0.84,
'Br': 1.20,
'Cd': 1.440,
'Cs': 2.440,'Ca': 1.760,'C': 0.690,'Ce': 2.040,
'Cl': 1.020,'Cr': 1.390,'Co': 1.260,'Cu': 1.320,'Cm': 1.690,
'Dy': 1.92,'Er': 1.890,
'Eu': 1.980,'F': 0.570,'Fr': 2.600,'Gd': 1.960,
'Ga': 1.220,'Ge': 1.20,'Au': 1.360,'Hf': 1.750,
'He': 0.280,'Ho': 1.920,'H': 0.310,'In': 1.420,'I': 1.390,
'Ir': 1.410,'Fe': 1.320,'Kr': 1.160,'La': 2.270,
'Pb': 1.460,'Li': 1.280,'Lu': 1.870,'Mg': 1.410,'Mn': 1.390,
'Hg': 1.320,'Mo': 1.540,'Nd': 2.010,
'Ne': 0.58,'Np': 1.900,'Ni': 1.240,'Nb': 1.640,'N': 0.710,
'Os': 1.440,'O': 0.660,'Pd': 1.390,'P': 1.070,
'Pt': 1.360,'Pu': 1.870,'Po': 1.400,'K': 2.030,'Pr': 2.030,
'Pm': 1.990,'Pa': 2.000,'Ra': 2.210,'Rn': 1.500,'Re': 1.510,
'Rh': 1.420,'Rb': 2.200,'Ru': 1.460,'Sm': 1.980,
'Sc': 1.700,'Se': 1.2,'Si': 1.110,
'Na': 1.660,'Sr': 1.950,'S': 1.050,'Ta': 1.700,'Tc': 1.470,
'Te': 1.380,'Tb': 1.940,'Tl': 1.450,'Th': 2.060,'Tm': 1.900,
'Sn': 1.390,'Ti': 1.600,'W': 1.620,'U': 1.960,'V': 1.530,
'Xe': 1.400,'Yb': 1.870,'Y': 1.900,'Zn': 1.220,'Zr': 1.750}

#List of metals which includes: 1A, 2A all Transition metals, 3A (excluding B), under the diagonal from Si to At, lanthanide, and actinide
metal_list = ['Li', 'Na', 'K', 'Rb', 'Cs', 'Fr', 'Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra', 'Sc',
'Y', 'La', 'Ac', 'Ti', 'Zr', 'Hf', 'V', 'Nb', 'Ta', 'Cr', 'Mo', 'W', 'Mn', 'Tc', 'Re', 'Fe',
'Ru', 'Os', 'Co', 'Rh', 'Ir', 'Ni', 'Pd', 'Pt', 'Cu', 'Ag', 'Au', 'Zn', 'Cd', 'Hg', 'Al',
'Ga', 'In', 'TL', 'Si', 'Ge', 'Sn', 'Pb', 'As', 'Sb', 'Bi', 'Te', 'Po', 'At', 'Ce', 'Pr',
'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Th', 'Pa', 'U',
'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']

def build_ASE_neighborlist(cif):
    """
    A function to build ASE neighbor list from a given cif.
    """
    radii = []
    for i in cif.get_chemical_symbols():
        #print i, cambridge_radii[i]
        radii.append(cambridge_radii[i])
#    ASE_neighborlist = NeighborList(radii,self_interaction=False, bothways=True)
    ASE_neighborlist = NeighborList(radii,self_interaction=False, bothways=True, skin=0.15)
    ASE_neighborlist.build(cif)
    return ASE_neighborlist

def find_clusters(adjacency_matrix):
    """
    A function to group the atoms together based on their connectivity.
    """
    clusters = []
    cluster_count, clusterIDs = connected_components(adjacency_matrix, directed=True)
    for n in range(cluster_count):
        clusters.append([i for i in range(atom_count) if clusterIDs[i] == n])
    return clusters

def find_metal_connected_atoms(struct, neighborlist):
    """
    A function to get indices of the atoms connected to the atoms belong to metal_list.
    """
    metal_connected_atoms = []
    metal_atoms = []
    for i, elem in enumerate(struct.get_chemical_symbols()):
        if elem in metal_list:
            neighbors, offsets = neighborlist.get_neighbors(i)
            #print "metal index: ", i, "neighbors: ", struct[neighbors].get_chemical_symbols()
            metal_connected_atoms.append(neighbors)
            metal_atoms.append(i)
    mca = metal_connected_atoms[0]
    #print struct[mca].get_chemical_symbols()
    #print struct[metal_connected_atoms[0]].get_chemical_symbols()
    return metal_connected_atoms, metal_atoms, struct

def CustomMatrix(neighborlist):
    """
    A (crude) function to create an adjacency matrix based on the neighbor list.
    """
    matrix = np.zeros((atom_count, atom_count), dtype=np.int)
    for i in range(atom_count):
        neighbors, offsets = neighborlist.get_neighbors(i)
        for j in neighbors:
            matrix[i][j] = 1
    return matrix

def mod_adjacency_matrix(adj_matrix, MetalConAtoms, MetalAtoms):
    """
    A function to only modify adjacency matrix elements if the total # of cluster changes.
    """
    clusters = find_clusters(adj_matrix)
    for i, element_1 in enumerate(MetalAtoms):
    	for j, element_2 in enumerate(MetalConAtoms[i]):
    		#print struct[element_2].symbol
    		if struct[element_2].symbol is "O":
    			#print "Examining bond between Metal - MetalCon: ", element_1, "-", element_2
    			tmp = len(find_clusters(adj_matrix))
    			adj_matrix[element_2][element_1] = 0
    			adj_matrix[element_1][element_2] = 0
    			new_clusters = find_clusters(adj_matrix)
    			#print "# of clusters: {} -> {}".format(tmp, len(find_clusters(adj_matrix)))
    			if tmp == len(find_clusters(adj_matrix)):
    				adj_matrix[element_2][element_1] = 1
    				adj_matrix[element_1][element_2] = 1
    			for ligand in new_clusters:
    				if ligand not in clusters:
    					tmp3 = struct[ligand].get_chemical_symbols()
    					#print len(tmp3)
    					if "O" and "H" in tmp3 and len(tmp3) == 2:
    						adj_matrix[element_2][element_1] = 1
    						adj_matrix[element_1][element_2] = 1
        #print "Loop #", i, "done"
    return adj_matrix

def makeP1CIF(fn):
    """
    A function to convert an input file to a formatted file from PyMatgen and return Atomic Simulation Environment Structure
    """
    tmp = read_structure(fn)
    pCell = pct().apply_transformation(tmp)
    dump = AseAtomsAdaptor.get_atoms(pCell)
    return dump

if __name__ == "__main__":

    # input
    fn = sys.argv[1]
    refcode = re.sub('[_fix.cif]','',fn)
    cif = makeP1CIF(fn)
    atom_count = len(cif.get_chemical_symbols())

    # Finds neighborlist based on atomic radii
    ASE_neighborlist = build_ASE_neighborlist(cif)
    # Build a custom adj. matrix
    a = CustomMatrix(ASE_neighborlist)
    b = find_clusters(a)
    b.sort(lambda x,y: cmp(len(x), len(y)))
    b.reverse()

    # a routine for catenated structures
    cluster_length=[]
    solvated_cluster = []
    for index, value in enumerate(b):
        #print cif[b[index]].get_chemical_formula(), "'s length: ", len(b[index])
        tmp = len(b[index])
        if len(cluster_length) > 0:
            if tmp > max(cluster_length):
                cluster_length = []
                solvated_cluster = []
                solvated_cluster.append(b[index])
                cluster_length.append(tmp)
            if tmp > 0.5 * max(cluster_length):
                #print "appending catenated!"
                solvated_cluster.append(b[index])
                cluster_length.append(tmp)
            else:
                continue
        else:
            solvated_cluster.append(b[index])
            cluster_length.append(tmp)

    # routine to merge solvated clusters
    solvated_merged = list(itertools.chain.from_iterable(solvated_cluster))
    atom_count = len(cif[solvated_merged].get_chemical_symbols())
    #print cif[solvated_merged].get_chemical_symbols()
    newASE_neighborlist = build_ASE_neighborlist(cif[solvated_merged])
    #sys.exit(0)
    MetalCon, MetalAtoms, struct = find_metal_connected_atoms(cif[solvated_merged], newASE_neighborlist)
    c = CustomMatrix(newASE_neighborlist)
    d = mod_adjacency_matrix(c, MetalCon, MetalAtoms)

    # finds clusters based on the modified adjacency matrix
    solvated_clusters2 = find_clusters(d)
    solvated_clusters2.sort(lambda x,y: cmp(len(x), len(y)))
    solvated_clusters2.reverse()

    # a routine for catenated structures
    cluster_length=[]
    final_clusters = []
    for index, value in enumerate(solvated_clusters2):
        #print struct[solvated_clusters2[index]].get_chemical_formula(), "'s cluster length: ", len(solvated_clusters2[index])
        tmp = len(solvated_clusters2[index])
        if len(cluster_length) > 0:
            if tmp > max(cluster_length):
                cluster_length = []
                final_clusters = []
                final_clusters.append(solvated_clusters2[index])
                cluster_length.append(tmp)
            if tmp > 0.5 * max(cluster_length):
                #print "appending catenated!"
                final_clusters.append(solvated_clusters2[index])
                cluster_length.append(tmp)
            else:
                continue
        else:
            final_clusters.append(solvated_clusters2[index])
            cluster_length.append(tmp)

    final_merged = list(itertools.chain.from_iterable(final_clusters))

    # look up stoichiometry in the structure
    bag_of_atoms = []
    tmp = struct[final_merged].get_chemical_symbols()
    tmp.sort()
    e = collections.Counter(tmp)

    for key, v in e.iteritems():
        if key == "C":
            num_of_c = float(v)
            #print key, v
        if key == "H":
            num_of_h = float(v)
            #print key, v
        else:
            continue
    print refcode, e.keys(), e.values()

    # output
    new_fn = refcode + '_clean.cif'
    # keep the largest cluster
    ase_write(new_fn, struct[final_merged])
