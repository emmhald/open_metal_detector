# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %(Emmanuel Haldoupis)s
"""
from __future__ import print_function
import sys,os
#sys.path.append('~/Dropbox/Work/python/modules')
path = os.path.expanduser('~/Dropbox/Work/python/modules')
#path = os.path.expanduser('~/dropbox/Work/python/modules')
if not path in sys.path:
    sys.path.insert(1, path)
del path
from xyz_read import xyz_file
from pymatgen import Lattice,Structure
from pymatgen.io.cifio import CifParser
from pymatgen.io.cifio import CifWriter
import numpy as np
import shlex
import shutil
import time
import math
import random
import pymatgen.io.smartio as sio
from atomic_parameters import atoms as ap
import json
import argparse
import pymatgen.io.xyzio as xyzio
import yappi

def main():

    parser = argparse.ArgumentParser(description = 'Split file into batches')
    parser.add_argument('-p','--parameter_file', nargs = '?',help='Name of parameters file', default = 'parameters.txt')
    parser.add_argument('-s','--summary_file', nargs = '?',help='Name of summary file', default = 'summary.out')
    parser.add_argument('-f','--folder_in', nargs = '?',help='Folder with structure files', default = 'CORE-MOF-DB-June2014')
    parser.add_argument('-o','--folder_out', nargs = '?',help='Folder with structure files', default = 'output')
    parser.add_argument('--continue_run', dest='continue_run', action='store_true')
    parser.add_argument('--attach_sorbate', dest='attach_sorbate', action='store_true')
    parser.add_argument('-m','--max_structures', nargs='?',default=10000000,const=10000000, type=int, help='The maximum number of structures to run')
    parser.set_defaults(continue_run = False)
    parser.set_defaults(attach_sorbate = False)

    args = parser.parse_args()
    cont = args.continue_run
    attach_ads = args.attach_sorbate
    params_filename = args.parameter_file
    source_folder = args.folder_in+'/'
    target_folder = args.folder_out+'/'
    number_of_structures = args.max_structures
    sfile = args.summary_file


    tolerance = dict()
    tolerance['plane'] = 35
    tolerance['plane_5l'] = 30
    tolerance['tetrahedron'] = 10
    tolerance['plane_on_metal'] = 12.5

    t0 = time.time()
    with open(params_filename,'r') as params_file:
        if not cont:
            clear_files(sfile,cont, target_folder)
        for i,struc in enumerate(params_file):
            t0s = time.time()
            uc_params=[]
            line_elements=shlex.split(struc)
            filename=line_elements[0]
            if filename.split('.')[0] == 'xyz':
                for j in range(1,7):
                    uc_params.append(float(line_elements[j]))
            analyze_structure(filename, uc_params, sfile, cont, source_folder, target_folder, attach_ads, tolerance)
            t1s = time.time()
            print('Time:',t1s-t0s)
            if i+1 >= number_of_structures:
                break
        params_file.close()
        t1 = time.time()
        print('Total Time',t1-t0)


def clear_files(sfile, cont, target_folder):
    make_folder('output')
    file_type = 'w'
    open_metal_mofs = open('output/open_metal_mofs.out',file_type)
    problematic_mofs = open('output/problematic.out',file_type)
    summary_mofs= open('output/'+sfile,file_type)
    open_metal_mofs.close()
    problematic_mofs.close()
    summary_mofs.close()

def make_folder(folder):
    if not os.path.exists(folder):
        os.makedirs(folder)

def delete_folder(folder_path):
    if os.path.exists(folder_path):
        for file_object in os.listdir(folder_path):
            file_object_path = os.path.join(folder_path, file_object)
            if os.path.isfile(file_object_path):
                os.unlink(file_object_path)
            else:
                shutil.rmtree(file_object_path)

def analyze_structure(filename, uc_params, sfile, cont, source_folder, target_folder, attach_ads, tolerance):

    mof_name = filename.split('.')[0]
    output_folder = target_folder+mof_name
    json_file_out = output_folder+'/'+mof_name+'.json'
    if cont and os.path.exists(json_file_out):
        print(mof_name,'has run already... skipping')
        return
    delete_folder(output_folder)
    make_folder(output_folder)
    make_folder(target_folder+'open_metal_mofs')
    make_folder(target_folder+'problematic_metal_mofs')

    open_metal_mofs = open(target_folder+'open_metal_mofs.out','a')
    problematic_mofs = open(target_folder+'problematic.out','a')
    summary_mofs = open(target_folder+sfile,'a')
    filetype = filename.split('.')[-1]
    if filetype=='xyz':
        xyz = xyz_file()
        xyz.filename_in = source_folder+filename #mof_name+'.xyz'
        if not os.path.isfile(xyz.filename_in):
            print('File not found',xyz.filename_in)
            return
        xyz.load_parameters(uc_params)
        xyz.open()
        lattice, system = make_system_from_xyz(xyz)
    elif filetype == 'cif':
        lattice, system = make_system_from_cif(source_folder+filename)
    else:
        sys.exit('Do not know this filetype')

    print("\n", filename)

    metal, organic = split_structure_to_organic_and_metal(system)
    if metal.num_sites == 0:
        print(mof_name+' not metal was found in structure', end="", file=summary_mofs)
        summary_mofs.close()
        return

    m_sa_frac,m_surface_area = 0.0,0.0
    #m_sa_frac,m_surface_areaget_metal_surface_areas(metal,system)

    #yappi.start()

    #first_coordination_structure, first_coordnation_structure_each_metal = find_first_coordination_sphere(metal, system)
    first_coordnation_structure_each_metal = find_all_coord_spheres(metal, system)
    output_json = get_output_dict(mof_name, m_surface_area, m_sa_frac, system.volume)

    cs_list = [] # list of coordination sequences for each open metal found
    for m,open_metal_candidate in enumerate(first_coordnation_structure_each_metal):
        site_dict = dict()
        op, pr, t, tf, min_dih, all_dih = check_if_open(open_metal_candidate, tolerance)
        site_dict = update_output_dict(site_dict, op, pr, t, tf, str(open_metal_candidate.species[0]), open_metal_candidate.num_sites-1, min_dih, all_dih)
        if op:
            if attach_ads:
                oms_index = match_index(str(open_metal_candidate.species[0]), open_metal_candidate.frac_coords[0], system)
                site_dict["oms_id"], cs_list = unique_site(oms_index, system, cs_list, output_folder, mof_name)
            if not output_json['metal_sites_found']:
                output_json['metal_sites_found'] = True
            if pr:
                output_json['problematic'] = True
        xyzio.XYZ(open_metal_candidate).write_file( output_folder+'/'+mof_name+'_first_coordination_sphere'+str(m)+'.xyz'  )
        output_json['metal_sites'].append(site_dict)

    print('Checking for OMSs done. Writing files')
    #yappi.get_func_stats().print_all()

    summary = write_summary(output_json)
    #if open_metal_site:
    if output_json['metal_sites_found']:
        print(summary, end="", file = open_metal_mofs)
        shutil.copyfile(source_folder+filename, target_folder+'open_metal_mofs/'+filename)
    print(summary, end="\n", file = summary_mofs)

    #if problematic_structure:
    if output_json['problematic'] :
        print(mof_name, end="", file = problematic_mofs)
        shutil.copyfile(source_folder+filename, target_folder+'problematic_metal_mofs/'+filename)

    write_xyz_file(output_folder+'/'+mof_name+'_metal.xyz',metal)
    write_xyz_file(output_folder+'/'+mof_name+'_organic.xyz',organic)
    write_xyz_file(output_folder+'/'+mof_name+'.xyz',system)
    system.make_supercell(2)
    write_xyz_file(output_folder+'/'+mof_name+'_super.xyz',system)

    open_metal_mofs.close()
    summary_mofs.close()
    problematic_mofs.close()
    with open(json_file_out, 'w') as outfile:
        json.dump(output_json, outfile,indent=3)

def get_output_dict(mof_name, m_surface_area, m_sa_frac, volume):
    output_json = dict()
    output_json['material_name'] = mof_name
    output_json['max_surface_area'] = m_surface_area
    output_json['max_surface_area_frac'] = m_sa_frac
    output_json['uc_volume'] = volume
    output_json['metal_sites_found'] = False
    output_json['problematic'] = False
    output_json['metal_sites'] = list()
    return output_json

def update_output_dict(site_dict, op, pr, t, tf, spec, open_metal_candidate_number, min_dih, all_dih):

    site_dict["is_open"] = op
    site_dict["t_factor"] = tf
    site_dict["metal"] = spec
    site_dict["type"] = 'closed'
    site_dict["number_of_linkers"] = open_metal_candidate_number
    site_dict["min_dihedral"] = min_dih
    site_dict["all_dihedrals"] = all_dih
    site_dict["unique"] = False
    site_dict["problematic"] = pr
    for ti in t:
        if t[ti]:
            site_dict["type"] = site_dict["type"]+','+str(ti)

    return site_dict

def unique_site(oms_index, system, cs_list, output_folder, mof_name):
    cs = find_coordination_sequence(oms_index, system)
    oms_id, new_site = find_oms_id(cs_list, cs)

    if new_site:
        print('New site found')
        cs_list.append(cs)

        end_to_end = 2.32
        eles = ['O', 'O', 'C']
        ads = add_co2_simple(system, oms_index, end_to_end, eles)
        mof_with_co2 = merge_structures(ads,system)
        cif = CifWriter(ads)
        cif.write_file(output_folder+'/'+mof_name+'_co2_'+str(oms_id)+'.cif')
        cif = CifWriter(mof_with_co2)
        cif.write_file(output_folder+'/'+mof_name+'_first_coordination_sphere_with_co2_'+str(oms_id)+'.cif')

        end_to_end = 1.1
        eles = ['N', 'N']
        ads = add_co2_simple(system, oms_index, end_to_end, eles)
        mof_with_co2 = merge_structures(ads,system)
        cif = CifWriter(ads)
        cif.write_file(output_folder+'/'+mof_name+'_n2_'+str(oms_id)+'.cif')
        cif = CifWriter(mof_with_co2)
        cif.write_file(output_folder+'/'+mof_name+'_first_coordination_sphere_with_n2_'+str(oms_id)+'.cif')
    return oms_id, cs_list

def find_oms_id(cs_list,cs):
    '''Check if a given site is unique based on its coordination sequence'''
    for i,cs_i in enumerate(cs_list):
        if compare_lists(cs_i,cs):
            return i,False
    return len(cs_list),True

def compare_lists(l1,l2):
    if len(l1) != len(l2):
        return False
    for i,j in zip(l1,l2):
        if i != j:
            return False
    return True

def write_summary(output_json):
    counter_om=0
    om_details=" "
    for d in output_json['metal_sites']:
        if d["unique"]:
            counter_om+=1
        if d["is_open"]:
            om_details=om_details+d['metal']
            om_details=om_details+','+str(d['number_of_linkers'])+'L'
            om_details=om_details+','+str(d['t_factor'])
            om_details=om_details+','+str(d['min_dihedral'])
    output_json['open_metal_density']=counter_om/output_json['uc_volume']
    open_metal_str='no'
    if output_json['metal_sites_found']:
        open_metal_str='yes'
    return output_json['material_name']+' '+str(counter_om)+' '+str(output_json['open_metal_density'])+' '+str(output_json['max_surface_area_frac'])+' '+str(output_json['max_surface_area'])+' '+open_metal_str+' '+om_details


def write_xyz_file(filename,system):
    xyzio.XYZ(system).write_file(filename+'.xyz')

def find_all_coord_spheres(centers, structure):
    coord_spheres = []
    for i,c in enumerate(centers):
        c_index = match_index(str(centers.species[i]), centers.frac_coords[i], structure)
        coord_spheres.append(find_coord_sphere(c_index, structure)[1])
    return coord_spheres

def find_coord_sphere(center, structure):
    dist = structure.lattice.get_all_distances(structure.frac_coords[center], structure.frac_coords)

    increase=1.0
    while True:
        bonds_found = 0
        first_coordnation_list_species = []
        first_coordnation_list_coords = []
        coord_sphere = []
        for i,dis in enumerate(dist[0]):
            bond_tol = ap.get_bond_tolerance(str(structure.species[center]),str(structure.species[i]))*increase
            if ap.bond_check(str(structure.species[center]), str(structure.species[i]), dis, bond_tol):
                first_coordnation_list_species.append(structure.species[i])
                first_coordnation_list_coords.append(structure.frac_coords[i])
                coord_sphere.append(i)
                bonds_found+=1
        check_structure = Structure(structure.lattice, first_coordnation_list_species, first_coordnation_list_coords)
        if ap.check_if_valid_bonds(check_structure, bond_tol, increase):
            break
        else:
            increase -= 0.1
            if increase < 0 :
                print('something went terribly wrong')
                raw_input()
            #print('Increasing bond_tolerance by ',increase)
    return coord_sphere, check_structure

def find_first_coordination_sphere(metal, structure_full):
    tol = 0.3
    first_coordination_structure = Structure(metal.lattice, metal.species, metal.frac_coords)
    first_coordnation_structure_each_metal = []
    for m, m_f_coor in enumerate(metal.frac_coords):
        structure = structure_full.copy()
        tmp1 = []
        tmp2 = []
        tmp1.append(str(metal.species[m]))
        tmp2.append([m_f_coor[0],m_f_coor[1],m_f_coor[2]])
        first_coordnation_structure_each_metal.append(Structure(metal.lattice,tmp1,tmp2))
        dist = metal.lattice.get_all_distances(m_f_coor,structure.frac_coords)
        for i, d in enumerate(dist[0]):
            if d <0.01 and str(structure.species[i]) == str(metal.species[m]):
                structure.__delitem__(i)
        dist = metal.lattice.get_all_distances(m_f_coor,structure.frac_coords)
        min_dist = min(dist[0])

        increase = 1.0
        while True:
            bonds_found=0
            first_coordnation_list_species=[]
            first_coordnation_list_coords=[]
            for i,dis in enumerate(dist[0]):
                bond_tol = ap.get_bond_tolerance(str(metal.species[m]),str(structure.species[i]))*increase
                if ap.bond_check(str(metal.species[m]),str(structure.species[i]),dis,bond_tol):
                    first_coordnation_list_species.append(structure.species[i])
                    first_coordnation_list_coords.append(structure.frac_coords[i])
                    bonds_found+=1
            #if check_if_enough_bonds(bonds_found,increase,str(metal.species[m])):
            check_structure = Structure(metal.lattice,first_coordnation_list_species,first_coordnation_list_coords)
            #print increase,bond_tol,bonds_found
            if ap.check_if_valid_bonds(check_structure,bond_tol,increase):
                break
            else:
                increase-=0.5
                if increase < 0 :
                    print('something went terribly wrong')
                    raw_input()

        print('Increased bond_tolerance by ',increase)
        for cls,clc in zip(first_coordnation_list_species,first_coordnation_list_coords):
            first_coordnation_structure_each_metal[m].append(cls,clc)
            first_coordination_structure.append(cls,clc)
        first_coordnation_structure_each_metal[m]=center_around_metal(first_coordnation_structure_each_metal[m])
    return first_coordination_structure, first_coordnation_structure_each_metal

#def check_if_enough_bonds(bonds_found,increase,metal):
#    if ap.is_lanthanide_or_actinide(metal):
#        min_bonds=5
#    else:
#        min_bonds=4
#    if bonds_found  < min_bonds and increase <= 1.5:
#        return False
#    else:
#        return True

def make_system_from_cif(ciffile):
    cif=CifParser(ciffile)
    system=cif.get_structures(primitive=False)
    print(ciffile)
    return system[0].lattice,system[0]

def make_latice_from_xyz(xyz):
    lattice=Lattice.from_parameters(xyz.uc_params[0], xyz.uc_params[1], xyz.uc_params[2], xyz.uc_params[3], xyz.uc_params[4], xyz.uc_params[5])
    return lattice

def make_system_from_xyz(xyz):
    lattice=make_latice_from_xyz(xyz)
    elements,coords=xyz.return_mp_structure_lists()
    structure=Structure(lattice,elements,coords,coords_are_cartesian=True)
    return lattice,structure

def split_structure_to_organic_and_metal(structure):
    coords=structure.frac_coords
    elements=structure.species
    coords_metal=[]
    coords_organic=[]
    elements_metal=[]
    elements_organic=[]
    for element, coord in zip(elements,coords):
        if ap.check_if_metal(str(element)):
            elements_metal.append(element)
            coords_metal.append(coord)
        else:
            elements_organic.append(element)
            coords_organic.append(coord)
    structure_metal = Structure(structure.lattice,elements_metal,coords_metal)
    structure_organic = Structure(structure.lattice,elements_organic,coords_organic)
    return structure_metal,structure_organic

def match_index(ele, f_coords,system):
    dist = system.lattice.get_all_distances(f_coords,system.frac_coords)
    for i, d in enumerate(dist[0]):
        if d < 0.001 and str(system.species[i]) == ele:
            return i

def check_if_6_or_more(system):
    if system.num_sites > 6:
        return True
    else:
        return False

def check_if_open(system, tolerance):
    tf = get_t_factor(system)

    problematic = False
    test = dict()
    test['plane'] = False
    test['same_side'] = False
    test['metal_plane'] = False
    test['4_or_less']=False
    test['non_TD'] = False

    open_metal_mof = False
    num = system.num_sites
    num_l = system.num_sites
    min_cordination = 3
    if ap.is_lanthanide_or_actinide(str(system.species[0])):
        min_cordination = 5

    if num_l <= min_cordination:
        open_metal_mof = True
        problematic = True
        test['4_or_less'] = True
        min_dihid = 0.0
        all_dihidrals = 0.0
    #    return open_metal_mof, problematic, test, tf, 0.0, 0.0
    else:
        open_metal_mof, test, min_dihid, all_dihidrals = check_non_metal_dihedrals(system, test, tolerance)
        if num_l == 5 and not open_metal_mof:
            open_metal_mof,test = check_metal_dihedrals(system, test, tolerance)

    return open_metal_mof, problematic, test, tf, min_dihid, all_dihidrals

def get_t_factor(system):
    angles=[]
    num=system.num_sites
    max_angle=0.0
    for i in range(1,num-1):
        for j in range(i+1,num):
            angles.append(system.get_angle(i,0,j))
            if max_angle < angles[-1]:
                max_index=i
                sec_index=j
                max_angle=angles[-1]

    angles.sort()
    if num-1>3 and num-1<6:
        beta=angles[-1]
        alpha=angles[-2]
        gamma=angles[0]
    elif num-1 ==6:
        angles_max=[]
        for j in range(1,num):
            if j != max_index:
                angles_max.append(system.get_angle(max_index,0,j))
        angles_max.sort()
        beta=angles_max[-2]
        alpha=angles_max[-1]
        gamma=angles_max[2]
#New Angles VOG
        max_angle=0.0
        for i in range(1,num-1):
            for j in range(i+1,num):
                if i!= max_index and i!= sec_index and j!= max_index and j!= sec_index:
                    angles.append(system.get_angle(i,0,j))
                    if max_angle < angles[-1]:
                         thr_index=i
                         fou_index=j
                         max_angle=angles[-1]
                         delta=angles[-1]

        for i in range(1,num-1):
            for j in range(i+1,num):
                if i!= max_index and i!= sec_index and i!= thr_index and i!= fou_index and j!= max_index and j!= sec_index and j!= thr_index and j!= fou_index:
                   epsilon=system.get_angle(i,0,j)

    if num-1==6:
        tau=get_t6_factor(alpha,beta,gamma,delta,epsilon)
    elif num-1==5:
        tau=get_t5_factor(alpha,beta)
    elif num-1==4:
        tau=get_t4_factor(alpha,beta)
    else:
        tau=-1
    return tau

def get_t4_factor(a,b):
    return (360-(a+b))/141.0

def get_t5_factor(a,b):
    return (b-a)/60.0

def get_t6_factor(a,b,c,d,e):
#    return 1.0-(a-c)/120
    #print 'abcde',a,b,c,d,e
#    print (abs(90.0-(a-b))/180.0)
#    print (abs(90.0-(a-c))/90.0)
#    print ((e-d)/180)
    return (e/180)
#    return (1-(180-e)/140 )
#    return ((540-a-d-e)/180)
#    return (1- (abs(90.0-(a-b))/180.0) - (abs(90.0-(a-c))/90.0) - ((e-d)/180))
#    return (abs(90.0-(a-b))/180.0)+(abs(90.0-(a-c))/90.0)


def check_metal_dihedrals(system, test, tolerance):
    num = system.num_sites
    num_l = system.num_sites - 1
    open_metal_mof = False
    crit=dict()
    tol=dict()

    crit['plane'] = 180
    tol['plane'] = tolerance['plane'] #30

    all_dihedrals, all_indeces = obtain_metal_dihedrals(num, system)
    min_dihedral = min(all_dihedrals)
    number_of_planes = 0
    for dihedral, indeces in zip(all_dihedrals, all_indeces):
        [i, j, k, l] = indeces
        if abs(dihedral - crit['plane']) < tol['plane'] or abs(dihedral - crit['plane']+180) < tol['plane']:
            number_of_planes += 1
            other_indeces = find_other_indeces([0, j, k, l], num)
            if len(other_indeces) > 2:
                print('Something went terribly wrong')
                input()
            #dihedral_other1 = system.get_dihedral(other_indeces[0],j,k,l)
            #dihedral_other2 = system.get_dihedral(other_indeces[1],j,k,l)
            dihedrals_other = []
            for o_i in other_indeces:
                dihedrals_other.append(system.get_dihedral(j,k,l,o_i))
            if not (check_positive(dihedrals_other) and check_negative(dihedrals_other)):
                open_metal_mof=True
                test['metal_plane']=True

    #i=0
    #for j in range(1,num):
    #    for k in range(1,num):
    #        for l in range(1,num):
    #    if (i==j or i==k or i==l or j==k or j==l or k==l):
    #       pass
    #    else:
    #       dihedral=abs(system.get_dihedral(i,j,k,l))
    #   if num_l == 5:

    if number_of_planes == 4:
        if open_metal_mof:
            print('conflicting criteria')
            input()
        else:
            open_metal_mof = False

    return open_metal_mof, test

def check_non_metal_dihedrals(system, test, tolerance):
    num = system.num_sites
    num_l = system.num_sites - 1
    crit = dict()
    tol = dict()
    crit['plane'] = 180
    tol['plane'] = tolerance['plane'] #35
    crit['plane_5l'] = 180
    tol['plane_5l'] = tolerance['plane_5l'] #30
    crit['tetrahedron'] = 70.528779  #70
    tol['tetrahedron'] = tolerance['tetrahedron'] #10
    open_metal_mof = False
    if num_l == 4:
        test_type='tetrahedron'
        om_type='non_TD'
    elif num_l == 5:
        test_type='plane_5l'
        om_type='plane_5l'
    elif num_l > 5:
        test_type='plane'
        om_type='same_side'

    all_dihedrals, all_indeces = obtain_dihedrals(num, system)
    min_dihedral = min(all_dihedrals)
    for dihedral, indeces in zip(all_dihedrals, all_indeces):
        [i, j, k, l] = indeces
        if num_l == 4:
            if not (abs(dihedral - crit[test_type]) < tol[test_type] or abs(dihedral - crit[test_type]+180) < tol[test_type]):
                test[om_type] = True
                test[test_type] = True
                open_metal_mof = True
            else :
                if abs(dihedral - crit[test_type]) < tol[test_type] or abs(dihedral - crit[test_type]+180) < tol[test_type]:
                    if num_l == 5:
                        test[om_type] = True
                        test[test_type] = True
                        open_metal_mof = True
                    elif num_l > 5:
                        if check_if_plane_on_metal(0, [i, j, k, l], system, tolerance):
                            other_indeces = find_other_indeces([0, i, j, k, l], num)
                            #check if other atoms are all in the same side of the plane
                            dihedrals_other = []
                            for o_i in other_indeces:
                                dihedrals_other.append(system.get_dihedral(j,k,l,o_i))
                            if not (check_positive(dihedrals_other) and check_negative(dihedrals_other)):
                                test[test_type] = True
                                test[om_type] = True
                                open_metal_mof = True

    return open_metal_mof, test, min_dihedral, all_dihedrals

def obtain_dihedrals(num, system):
    all_dihedrals = []
    indeces = []
    for i in range(1,num):
        for j in range(1,num):
            for k in range(1,num):
                for l in range(1,num):
                    if (i==j or i==k or i==l or j==k or j==l or k==l):
                       pass
                    else:
                        dihedral = abs(system.get_dihedral(i,j,k,l))
                        all_dihedrals.append(dihedral)
                        indeces.append([i, j, k , l])
    return all_dihedrals, indeces

def obtain_metal_dihedrals(num, system):
    all_dihedrals = []
    indeces = []
    i = 0
    for j in range(1,num):
        for k in range(1,num):
            for l in range(1,num):
                if (i==j or i==k or i==l or j==k or j==l or k==l):
                   pass
                else:
                    dihedral = abs(system.get_dihedral(i,j,k,l))
                    all_dihedrals.append(dihedral)
                    indeces.append([i, j, k , l])
    return all_dihedrals, indeces

def add_co2(v1,v2,v3,v4,system):
    ads_dist = 2.2
    #old method for adding co2.
    if 1==2:
        p1 = calc_plane(v1,v2,v3)
        p2 = calc_plane(v2,v3,v4)
        p_avg_up = [ads_dist*(p_1+p_2)/2 for p_1,p_2 in zip(p1,p2)]
        p_avg_down = [-p for p in p_avg_up]
        p_avg_up = p_avg_up+system.cart_coords[0]
        p_avg_down = p_avg_down+system.cart_coords[0]
        p_avg_up_f = system.lattice.get_fractional_coords(p_avg_up)
        p_avg_down_f = system.lattice.get_fractional_coords(p_avg_down)
        dist_up = min(system.lattice.get_all_distances(p_avg_up_f,system.frac_coords)[0])
        dist_down = min(system.lattice.get_all_distances(p_avg_down_f,system.frac_coords)[0])
        if dist_up < dist_down:
            p_avg_f = p_avg_down_f
            direction = -1
        else:
            p_avg_f = p_avg_up_f
            direction = 1
        co2_vector_C = [direction*(ads_dist+1.16)*(p_1+p_2)/2 for p_1,p_2 in zip(p1,p2)]+system.cart_coords[0]
        co2_vector_O = [direction*(ads_dist+2*1.16)*(p_1+p_2)/2 for p_1,p_2 in zip(p1,p2)]+system.cart_coords[0]
        co2_vector_C_f = system.lattice.get_fractional_coords(co2_vector_C)
        co2_vector_O_f = system.lattice.get_fractional_coords(co2_vector_O)

def add_co2_simple(structure, oms_index, end_to_end, eles):

    ads_dist = 2.2
    #end_to_end = 2.32
    bond = end_to_end/2
    #eles = ['O', 'O', 'C']
    adsorption_site = []
    adsorption_pos = []
    ads_vector = []
    ads_vector_f = []

    ads_vector.append(find_adsorption_site(structure,structure.cart_coords[oms_index], ads_dist))
    ads_vector.append(find_adsorption_site(structure, ads_vector[0],  end_to_end))
    if len(eles) > 2:
        r = [i-j for i,j in zip(ads_vector[1], ads_vector[0])]
        r_l = np.linalg.norm(np.array(r))
        r_unit = [i/r_l for i in r]
        pos = list(rr*bond for rr in r_unit)
        ads_vector.append([j-i for i,j in zip(pos, ads_vector[1])])

    ads_vector_f = [ structure.lattice.get_fractional_coords(ads_v) for ads_v in ads_vector]

    for e in eles:
        adsorption_site.append(e)
    for a in ads_vector_f:
        adsorption_pos.append(a)
    dists=[]
    for s,p in zip(adsorption_site,adsorption_pos):
        dists.append(min(structure.lattice.get_all_distances(p,structure.frac_coords)[0]))
    print('Min adsorbate distance from framework:',min(dists),max(dists))
    return Structure(structure.lattice, adsorption_site, adsorption_pos)

def find_adsorption_site(system,center,prob_dist):
    #find the adsorption site by maximizing the distance from all the atoms while
    #keep the distance fixed at some predefined distace
    tries = 10000
    sum_distance = []
    probe_positions = []
    for i in range(0,tries):
        #probe_pos=generate_random_position(system.cart_coords[0],center,(prob_dist)-atom.get_vdf_radius(center_e))
        probe_pos = generate_random_position(center,prob_dist)
        probe_pos_f=system.lattice.get_fractional_coords(probe_pos)
        sum_distance.append(sum([1.0/(r**12) for r in system.lattice.get_all_distances(probe_pos_f,system.frac_coords)[0]]))
        probe_positions.append(probe_pos)
    new_position = probe_positions[sum_distance.index(min(sum_distance))]
    return new_position

def calc_plane(x, y, z):
    v1 = [y[0] - x[0], y[1] - x[1], y[2] - x[2]]
    v2 = [z[0] - x[0], z[1] - x[1], z[2] - x[2]]
    cross_product = [v1[1]*v2[2]-v1[2]*v2[1], v1[2] * v2[0] - v1[0] * v2[2], v1[0] * v2[1] - v1[1] * v2[0]]
    d = cross_product[0] * x[0] - cross_product[1] * x[1] + cross_product[2] * x[2]
    a = cross_product[0]
    b = cross_product[1]
    c = cross_product[2]
    plane_v=[a,b,c]
    plane_v_norm = plane_v/np.linalg.norm(plane_v)
    return plane_v_norm

def check_if_plane_on_metal(m_i, indeces, system, tolerance):
    crit = 180
    tol = tolerance['plane_on_metal'] #12.5   #set to 12.5 so that ferocene type coordination spheres are dected correctly. eg. BUCROH
    #tol = 25.0
    for i in range(1,len(indeces)):
        for j in range(1,len(indeces)):
            for k in range(1,len(indeces)):
                if (i==j or i==k or j==k):
                    pass
                else:
                    dihedral=abs(system.get_dihedral(m_i,indeces[i],indeces[j],indeces[k]))
                    if abs(dihedral-crit)< tol or abs(dihedral-crit+180)< tol:
                        return True
    return False


def check_positive(N):
    for n in N:
        if n > 0:return True

def check_negative(N):
    for n in N:
        if n < 0:return True

def find_other_indeces(indeces,num):
    other_indeces=[]
    for index in range(0,num):
        if not index in indeces:
            other_indeces.append(index)
    return other_indeces

def center_around_metal(system):
    #return system
    center=system.frac_coords[0]
    tmp1=[]
    tmp2=[]
    tmp1.append(str(system.species[0]))
    tmp2.append([system.frac_coords[0][0],system.frac_coords[0][1],system.frac_coords[0][2]])
    system_centered=Structure(system.lattice,tmp1,tmp2)
    for i in range(1,system.num_sites):
        c_i=system.frac_coords[i]
        dist_vector=center-c_i
        dist_vector_r=[]
        for j in range(0,3):
            dist_vector_r.append(round(dist_vector[j]))
        dist_before=np.linalg.norm(system.lattice.get_cartesian_coords(center)-system.lattice.get_cartesian_coords(c_i))
        c_i_centered=c_i+dist_vector_r
        dist_after=np.linalg.norm(system.lattice.get_cartesian_coords(center)-system.lattice.get_cartesian_coords(c_i_centered))
        if dist_after>dist_before:
            for j in range(0,3):
                dist_vector_r[j]=np.rint(dist_vector[j])
            c_i_centered=c_i+dist_vector_r
            if dist_after>dist_before:
                c_i_centered=c_i
        system_centered.append(system.species[i],c_i_centered)
    return system_centered


def merge_structures(s1,s2):
    sites=[]
    posistions=[]
    for e,c in zip(s1.species,s1.frac_coords):
        sites.append(e)
        posistions.append(c)
    for e,c in zip(s2.species,s2.frac_coords):
        sites.append(e)
        posistions.append(c)
    if s1.lattice != s2.lattice:
        sys.exit('Trying to merger two structures with different lattices')
    return Structure(s1.lattice,sites,posistions)

def get_metal_surface_areas(metal,system):
    sa_list=[]
    for m,m_coor in enumerate(metal.frac_coords):
        sub_system=make_subsystem(m_coor,system,5.0) #make a structure of atoms withn 7.0 the metal
        sa_list.append(get_metal_surface_area(m_coor,str(metal.species[m]),sub_system))
    s_max=max(sa_list)
    return s_max

def make_subsystem(coord,system,dist_check):
    distances=system.lattice.get_all_distances(coord,system.frac_coords)
    coords=[]
    elements=[]
    for i,dist in enumerate(distances[0]):
        if dist<dist_check:
            if dist>0.1:#exclude the central metal
                elements.append(system.species[i])
                coords.append(system.frac_coords[i])
    return Structure(system.lattice,elements,coords)

def get_metal_surface_area(fcenter,metal_element,system):
    center=system.lattice.get_cartesian_coords(fcenter)

    vdw_probe=1.8405 #1.52 #0.25 #2.1 #1.52 #vdw radius for oxygen
    metal_full_surface_area=sphere_area(metal_element,vdw_probe) #use 0.0 for vdw_probe to get sphere of metal
    count=0
    mc_tries=5000
    params_file= open('test.txt','w')
    for i in range(0,mc_tries):
        dist=ap.get_vdf_radius(metal_element)+vdw_probe
        pos=generate_random_position(center,dist) #vdw_probe
        #pos=generate_random_position(center,metal_element,vdw_probe) #vdw_probe
        print >>params_file,'xx',pos[0],pos[1],pos[2]
        pos_f=system.lattice.get_fractional_coords(pos)
        if not check_for_overalp(center,pos_f,system,vdw_probe):
            count+=1
    sa_frac=float(count)/float(mc_tries) #metal_full_surface_area
    sa=metal_full_surface_area*sa_frac
    return sa_frac,sa

def check_for_overalp(center,pos,system,r_probe):
    #pos=[0.0,0.0,0.0]
    #print  'N 1.0',pos[0],pos[1],pos[2],'biso 1.0 N'
    distances=system.lattice.get_all_distances(pos,system.frac_coords)
    for i,dist in enumerate(distances[0]):
        #if not check_if_center(center,system.cart_coords[i],system):
        #print dist-r_probe+atom.get_vdf_radius(str(system.species[i]))
        #raw_input()
        if (dist - (r_probe+ap.get_vdf_radius(str(system.species[i]))) ) < -1e-4:
            return True
    return False

def check_if_center(center,test_coords,system):
    distances=system.lattice.get_all_distances(test_coords,center)
    if distances[0]< 0.1:
        return True
    else:
        return False

def sphere_area(metal_element,probe_r):
    r=ap.get_vdf_radius(metal_element)
    r=r+probe_r
    return 4*math.pi*r*r

def generate_random_position(center,dist):
    r=generate_random_vector()
    #dist=atom.get_vdf_radius(metal)+vdw_probe
    pos=list(rr*dist for rr in r)
    pos=[i+j for i,j in zip(pos, center)]
    return pos

def generate_random_vector():
    zeta_sq=2.0
    ran_vec=[]
    while zeta_sq>1.0:
        xi1=random.random()
        xi2=random.random()
        zeta1=1.0-2.0*xi1
        zeta2=1.0-2.0*xi2
        zeta_sq=(zeta1*zeta1+zeta2*zeta2)

    ranh=2.0*math.sqrt(1.0-zeta_sq)
    ran_vec.append(zeta1*ranh)
    ran_vec.append(zeta2*ranh)
    ran_vec.append(1.0-2.0*zeta_sq)
    return ran_vec

def find_coordination_sequence(center, structure):
    '''computes the coordination sequence up to the Nth coordination shell
    as input it takes the MOF as a pymatgen Structure and the index of the center metal in the Structure
    '''
    #The shell_list is a set with the index of each atom and its unit cell index realtive to a cetral unit cell
    shell_list = set([(center,(0,0,0))])
    shell_list_prev = set([])
    all_shells = set(shell_list)
    n_shells = 6
    cs = []
    ele = [(str(structure.species[center]))]
    coords = [ [structure.frac_coords[center][0], structure.frac_coords[center][1], structure.frac_coords[center][2]]]
    coordination_structure = (Structure(structure.lattice,ele,coords))
    for n in range(0, n_shells):
        c_set = set([])
        for a_uc in shell_list:
            a = a_uc[0]
            lattice = a_uc[1]
            coord_sphere = find_coord_sphere(a ,structure)[0]
            coord_sphere_with_uc = []
            for c in coord_sphere:
                dist = structure.lattice.get_all_distance_and_image(structure.frac_coords[a],structure.frac_coords[c])
                uc = lattice -min(dist)[1]
                coord_sphere_with_uc.append((c,tuple(uc)))
            coord_sphere_with_uc = tuple(coord_sphere_with_uc)
            c_set = c_set.union(set(coord_sphere_with_uc))
        for a in shell_list_prev:
            c_set.discard(a)
        for a in shell_list:
            c_set.discard(a)
        for i_uc in c_set:
            i = i_uc[0]
            ele = ap().elements[n+3]
            coords = [structure.frac_coords[i][0], structure.frac_coords[i][1], structure.frac_coords[i][2]]
            coordination_structure.append(ele,coords)

        cs.append(len(c_set))
        all_shells = all_shells.union(c_set)
        shell_list_prev = set(shell_list)
        shell_list = set(c_set)
    #coordination_structure = center_around_metal(coordination_structure)
    write_xyz_file('temp.xyz',coordination_structure)
    return cs


if __name__=='__main__':
    main()
