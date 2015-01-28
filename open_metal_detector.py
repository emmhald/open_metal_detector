# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %(Emmanuel Haldoupis)s
"""
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
from atomic_parameters import atoms
import json
import argparse

atom=atoms()
target_folder='CORE-MOF-DB-June2014/'
filetype='cif'
def main():
    global filetype

    parser = argparse.ArgumentParser(description='Split file into batches')
    parser.add_argument('-p','--parameter_file', nargs='?',help='Name of parameters file')
    parser.add_argument('-s','--summary_file', nargs='?',help='Name of summary file')
    parser.add_argument('--continue_run', dest='continue_run', action='store_true')
    parser.add_argument('-m','--max_structures', nargs='?',default=10000000,const=10000000, type=int, help='The maximum number of structures to run')
    parser.set_defaults(continue_run=False)

    args=parser.parse_args()
    cont=args.continue_run
    params_filename='parameters_'+filetype+'.txt'
    if args.parameter_file:
        params_filename=args.parameter_file

    number_of_structures=args.max_structures
    sfile='summary.out'
    if args.summary_file:
        sfile=args.summary_file
    t0 = time.time()
    with open(params_filename,'r') as params_file:
        if not cont:
            clear_files(sfile,cont)
        for i,struc in enumerate(params_file):
            t0s = time.time()
            uc_params=[]
            line_elements=shlex.split(struc)
            filename=line_elements[0]
            if filetype=='xyz':
                for j in range(1,7):
                    uc_params.append(float(line_elements[j]))
            analyze_structure(filename,uc_params,sfile,cont)
            t1s = time.time()
            print 'Time:',t1s-t0s
            if i+1 >= number_of_structures:
                break
        params_file.close()
        t1 = time.time()
        print 'Total Time',t1-t0


def clear_files(sfile,cont):
    make_folder('output')
    file_type='w'
    open_metal_mofs= open('output/open_metal_mofs.out',file_type)
    problematic_mofs= open('output/problematic.out',file_type)
    summary_mofs= open('output/'+sfile,file_type)
    open_metal_mofs.close()
    problematic_mofs.close()
    summary_mofs.close()

def make_folder(folder):
    if not os.path.exists(folder):
        os.makedirs(folder)

def analyze_structure(filename,uc_params,sfile,cont):
    global target_folder

    mof_name=filename.split('.')[0]
    output_folder='output/'+mof_name
    json_file_out=output_folder+'/'+mof_name+'.json'
    if cont and os.path.exists(json_file_out):
        print mof_name,'has run already... skipping'
        return
    make_folder(output_folder)
    make_folder('output/open_metal_mofs')
    make_folder('output/problematic_metal_mofs')

    open_metal_mofs= open('output/open_metal_mofs.out','a')
    problematic_mofs= open('output/problematic.out','a')
    summary_mofs= open('output/'+sfile,'a')

    if filetype=='xyz':
        print 'Reading Xyzs'
        xyz=xyz_file()
        xyz.filename_in=target_folder+filename #mof_name+'.xyz'
        if not os.path.isfile(xyz.filename_in):
            print 'File not found',xyz.filename_in
            return
        xyz.load_parameters(uc_params)
        xyz.open()
        lattice,system=make_system_from_xyz(xyz)
    elif filetype=='cif':
        print 'Reading Cifs'
        lattice,system=make_system_from_cif(target_folder+filename)
    else:
        sys.exit('Do not know this filetype')
    print filename

    metal,organic=split_structure_to_organic_and_metal(system)
    if metal.num_sites == 0:
        print>>summary_mofs,mof_name+' not metal was found in structure'
        summary_mofs.close()
        return
    first_coordination_structure,first_coordnation_structure_each_metal=find_first_coordination_sphere(metal,system)
    m_sa_frac,m_surface_area=0.0,0.0
    #m_sa_frac,m_surface_areaget_metal_surface_areas(metal,system)

    open_metal_site=False
    problematic_structure=False
    test=[] #dict()
    output_json=dict()
    output_json['material_name']=mof_name
    output_json['max_surface_area']=m_surface_area
    output_json['max_surface_area_frac']=m_sa_frac
    output_json['uc_volume']=system.volume
    output_json['metal_sites_found']=False
    output_json['metal_sites']=list()

    #summary=str(m_sa_frac)+' '+str(m_surface_area)+' '
    #summary=summary+' '+'no'
    count_omsites=0
    cs_list = [] # list of coordination sequences for each open metal found
    for m,open_metal_candidate in enumerate(first_coordnation_structure_each_metal):
        site_dict=dict()
        op,pr,t,ads,tf,min_dih,all_dih=check_if_open(open_metal_candidate)
        site_dict["is_open"]=op
        site_dict["t_factor"]=tf
        site_dict["metal"]=str(open_metal_candidate.species[0])
        site_dict["type"]='closed'
        site_dict["number_of_linkers"]=open_metal_candidate.num_sites-1
        site_dict["min_dihedral"]=min_dih
        site_dict["all_dihedrals"]=all_dih
        site_dict["unique"]=False

        if op:
            count_omsites+=1
            oms_index = match_index(str(open_metal_candidate.species[0]),open_metal_candidate.frac_coords[0],system)
            cs = find_coordination_sequence(oms_index, system)
            unique_site = check_if_unique(cs_list,cs)
            if len(ads) > 0 and unique_site:
                cs_list.append(cs)
                mof_with_co2=merge_structures(ads,system)
                cif=CifWriter(ads)
                cif.write_file(output_folder+'/'+mof_name+'_ads'+str(count_omsites)+'.cif')
                cif=CifWriter(mof_with_co2)
                cif.write_file(output_folder+'/'+mof_name+'_first_coordination_sphere_with_ads'+str(count_omsites)+'.cif')
            open_metal_site=op
            problematic_structure=pr
            test.append(t)
            if not output_json['metal_sites_found']:
                output_json['metal_sites_found']=True
            #if not 'yes' in summary:
            #    summary=summary.replace('no','yes') #=summary+' '+'yes'
            #string=' '+str(open_metal_candidate.species[0])+','+str(open_metal_candidate.num_sites-1)+'L'
            site_dict["type"]=''
            for ti in t:
                if t[ti]:
                    #string=string+','+str(ti)
                    site_dict["type"]=site_dict["type"]+','+str(ti)
            #if not string in summary and not str(tf) in summary:
            #    unique_site=True
            #    site_dict["unique"]=True
        #else:
            #string=' '+str(open_metal_candidate.species[0])+','+str(open_metal_candidate.num_sites-1)+'L'
            #unique_site=False
            #if str(tf) in summary: #not string in summary and not
            #    unique_site=True
            #    site_dict["unique"]=True
        #if unique_site:
        #    summary=summary+' '+string
        #    summary=summary+','+str(tf)
        #    summary=summary+','+str(min_dih)

        xyz_first_coord_sphere=xyz_file()
        xyz_first_coord_sphere.make_xyz(output_folder+'/'+mof_name+'_first_coordination_sphere'+str(m)+'.xyz',open_metal_candidate.cart_coords,open_metal_candidate.species)
        output_json['metal_sites'].append(site_dict)
    #summary=mof_name+' '+str(count_omsites)+' '+str(count_omsites/system.volume)+' '+summary
    summary=write_summary(output_json)

    if open_metal_site:
        print>>open_metal_mofs,summary
        shutil.copyfile(target_folder+filename, 'output/open_metal_mofs/'+filename)
#    else:
#        summary=summary+' '+'no'
    print>>summary_mofs,summary
    print 'max metal surface area',m_surface_area
    print 'open_metal_site',open_metal_site
    print 'problematic_structure',problematic_structure,"\n-----"

    if problematic_structure:
        print>>problematic_mofs,mof_name
        shutil.copyfile(target_folder+filename, 'output/problematic_metal_mofs/'+filename)

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


def check_if_unique(cs_list,cs):
    '''Check if a given site is unique based on its coordination sequence'''
    for cs_i in cs_list:
        if compare_lists(cs_i,cs):
            return False
    return True

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
    xyz=xyz_file()
    xyz.make_xyz(filename,system.cart_coords,system.species)


def find_coord_sphere(center,structure):
    dist=structure.lattice.get_all_distances(structure.frac_coords[center],structure.frac_coords)

    increase=1.3
    while True:
        bonds_found=0
        first_coordnation_list_species=[]
        first_coordnation_list_coords=[]
        coord_sphere=[]
        for i,dis in enumerate(dist[0]):
            bond_tol=get_bond_tolerance(str(structure.species[center]),str(structure.species[i]))*increase
            if bond_check(str(structure.species[center]),str(structure.species[i]),dis,bond_tol):
                first_coordnation_list_species.append(structure.species[i])
                first_coordnation_list_coords.append(structure.frac_coords[i])
                coord_sphere.append(i)
                bonds_found+=1
        check_structure=Structure(structure.lattice,first_coordnation_list_species,first_coordnation_list_coords)
        if check_if_valid_bonds(check_structure,bond_tol,increase):
            break
        else:
            increase-=0.1
            if increase < 0 :
                print 'something went terribly wrong'
                raw_input()
        #    print 'Increasing bond_tolerance by ',increase
    #print coord_sphere
    #i = center
    #print structure.species[i],structure.cart_coords[i][0],structure.cart_coords[i][1],structure.cart_coords[i][2]
    #for i in coord_sphere:
    #    print structure.species[i],structure.cart_coords[i][0],structure.cart_coords[i][1],structure.cart_coords[i][2]
    #raw_input()
    return coord_sphere


def find_first_coordination_sphere(metal,structure_full):
    tol=0.3
    first_coordination_structure=Structure(metal.lattice,metal.species,metal.frac_coords)
    first_coordnation_structure_each_metal=[]
    for m,m_f_coor in enumerate(metal.frac_coords):
        structure=structure_full.copy()
        tmp1=[]
        tmp2=[]
        tmp1.append(str(metal.species[m]))
        tmp2.append([m_f_coor[0],m_f_coor[1],m_f_coor[2]])
        first_coordnation_structure_each_metal.append(Structure(metal.lattice,tmp1,tmp2))
        dist=metal.lattice.get_all_distances(m_f_coor,structure.frac_coords)
        for i,d in enumerate(dist[0]):
            if d <0.01 and str(structure.species[i]) == str(metal.species[m]):
                structure.__delitem__(i)
        dist=metal.lattice.get_all_distances(m_f_coor,structure.frac_coords)
        min_dist=min(dist[0])

        increase=1.3
        while True:
            bonds_found=0
            first_coordnation_list_species=[]
            first_coordnation_list_coords=[]
            for i,dis in enumerate(dist[0]):
                bond_tol=get_bond_tolerance(str(metal.species[m]),str(structure.species[i]))*increase
                if bond_check(str(metal.species[m]),str(structure.species[i]),dis,bond_tol):
                    first_coordnation_list_species.append(structure.species[i])
                    first_coordnation_list_coords.append(structure.frac_coords[i])
                    bonds_found+=1
            #if check_if_enough_bonds(bonds_found,increase,str(metal.species[m])):
            check_structure=Structure(metal.lattice,first_coordnation_list_species,first_coordnation_list_coords)
            #print increase,bond_tol,bonds_found
            if check_if_valid_bonds(check_structure,bond_tol,increase):
                break
            else:
                increase-=0.1
                if increase < 0 :
                    print 'something went terribly wrong'
                    raw_input()
                print 'Increasing bond_tolerance by ',increase

        for cls,clc in zip(first_coordnation_list_species,first_coordnation_list_coords):
            first_coordnation_structure_each_metal[m].append(cls,clc)
            first_coordination_structure.append(cls,clc)
        first_coordnation_structure_each_metal[m]=center_around_metal(first_coordnation_structure_each_metal[m])
    return first_coordination_structure,first_coordnation_structure_each_metal

def check_if_enough_bonds(bonds_found,increase,metal):
    if atom.is_lanthanide_or_actinide(metal):
        min_bonds=5
    else:
        min_bonds=4
    if bonds_found  < min_bonds and increase <= 1.5:
        return False
    else:
        return True

def check_if_valid_bonds(ligands,bond_tol,increase):
    if increase >0.05:
        for l,fc in enumerate(ligands.frac_coords):
            dist=ligands.lattice.get_all_distances(fc,ligands.frac_coords)
            count_bonds=0
            for i,dis in enumerate(dist[0]):
                if i != l:
                    if bond_check(str(ligands.species[l]),str(ligands.species[i]),dis,bond_tol):
                        count_bonds+=1
            if count_bonds > 0:
                return False
    return True


def make_system_from_cif(ciffile):
    cif=CifParser(ciffile)
    system=cif.get_structures(primitive=False)
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
    for element,coord in zip(elements,coords):
        if atom.check_if_metal(str(element)):
            elements_metal.append(element)
            coords_metal.append(coord)
        else:
            elements_organic.append(element)
            coords_organic.append(coord)
    structure_metal=Structure(structure.lattice,elements_metal,coords_metal)
    structure_organic=Structure(structure.lattice,elements_organic,coords_organic)

    return structure_metal,structure_organic

def match_index(ele,f_coords,system):
    dist=system.lattice.get_all_distances(f_coords,system.frac_coords)
    for i,d in enumerate(dist[0]):
        if d <0.001 and str(system.species[i]) == ele:
            return i

def check_if_6_or_more(system):
    if system.num_sites > 6:
        return True
    else:
        return False


def check_if_open(system):
#def check_if_pyramidal_square(system):
    tf=get_t_factor(system)

    problematic=False
    test=dict()
    test['plane']=False
    test['same_side']=False
    test['metal_plane']=False
    test['4_or_less']=False
    test['non_TD']=False

    open_metal_mof=False
    dihedral_tolerance=5
    num=system.num_sites
    max_cordination=3

    #if num-1 ==4:
    #    print 'tfactor',tf
    if atom.is_lanthanide_or_actinide(str(system.species[0])):
        max_cordination=5
    if num-1 <= max_cordination:
        open_metal_mof=True
        problematic=True
        test['4_or_less']=True
        if num-1 > 2:
            v1=system.cart_coords[0]
            v2=system.cart_coords[1]
            v3=system.cart_coords[2]
            v4=system.cart_coords[2]
        else:
            v1=[0,0,0]
            v2=[0,0,0]
            v3=[0,0,0]
            v4=[0,0,0]
        ads=add_co2(v1,v2,v3,v4,system)
        return open_metal_mof,problematic,test,ads,tf,0.0,0.0
    ads=[]
    open_metal_mof,test,ads,min_dihid,all_dihidrals=check_non_metal_dihedrals(system,test)
    if num-1 == 5 and not open_metal_mof:
        open_metal_mof,test,ads=check_metal_dihedrals(system,test)

    print num-1,open_metal_mof,tf
    #raw_input()
    return open_metal_mof,problematic,test,ads,tf,min_dihid,all_dihidrals

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
#                print i,j,angles[-1]

    angles.sort()
    if num-1>3 and num-1<6:
        beta=angles[-1]
        alpha=angles[-2]
        gamma=angles[0]
    elif num-1 ==6:
        print angles
        angles_max=[]
        for j in range(1,num):
            if j != max_index:
                angles_max.append(system.get_angle(max_index,0,j))
        angles_max.sort()
#        print angles_max
#        print angles_max[0], angles_max[1]
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
    print 'abcde',a,b,c,d,e
#    print (abs(90.0-(a-b))/180.0)
#    print (abs(90.0-(a-c))/90.0)
#    print ((e-d)/180)
    return (e/180)
#    return (1-(180-e)/140 )
#    return ((540-a-d-e)/180)
#    return (1- (abs(90.0-(a-b))/180.0) - (abs(90.0-(a-c))/90.0) - ((e-d)/180))
#    return (abs(90.0-(a-b))/180.0)+(abs(90.0-(a-c))/90.0)


def check_metal_dihedrals(system,test):
    open_metal_mof=False
    crit=dict()
    tol=dict()

    crit['plane']=180
    tol['plane']=30

    crit['tetrahedron']=70
    tol['tetrahedron']=10

    num=system.num_sites

    number_of_planes=0
    i=0
    for j in range(1,num):
        for k in range(1,num):
            for l in range(1,num):
                if (i==j or i==k or i==l or j==k or j==l or k==l):
                   pass
                else:
                   dihedral=abs(system.get_dihedral(i,j,k,l))
                   if num-1==5:
                       if abs(dihedral-crit['plane'])< tol['plane'] or abs(dihedral-crit['plane']+180)< tol['plane']:
                         number_of_planes+=1
                         other_indeces=find_other_indeces([0,j,k,l],num)
                         if len(other_indeces) > 2:
                             print 'Something went terribly wrong'
                             raw_input()
                         #dihedral_other1=abs(system.get_dihedral(other_indeces[0],j,k,l))
                         #dihedral_other2=abs(system.get_dihedral(other_indeces[0],j,k,l))
                         dihedral_other1=system.get_dihedral(other_indeces[0],j,k,l)
                         dihedral_other2=system.get_dihedral(other_indeces[1],j,k,l)
                         if (dihedral_other1*dihedral_other2)  >= 0:
                             open_metal_mof=True
                             test['metal_plane']=True

    if number_of_planes == 4:
        if open_metal_mof:
            print 'conflicting criteria'
            raw_input()
        else:
            open_metal_mof=False
    metal_cluster_with_adsorbate=[]
    if open_metal_mof and test['metal_plane']:
        v1=system.cart_coords[0]
        v2=system.cart_coords[1]
        v3=system.cart_coords[2]
        v4=system.cart_coords[3]
        metal_cluster_with_adsorbate=add_co2(v1,v2,v3,v4,system)
    return open_metal_mof,test,metal_cluster_with_adsorbate


def check_non_metal_dihedrals(system,test):
    num=system.num_sites
    open_metal_mof=False
    crit=dict()
    tol=dict()
    crit['plane']=180
    tol['plane']=35
    crit['plane_5l']=180
    tol['plane_5l']=30
    crit['tetrahedron']=70
    tol['tetrahedron']=10
    if num-1==4:
        test_type='tetrahedron'
        om_type='non_TD'
    elif num-1 == 5:
        test_type='plane_5l'
        om_type='plane_5l'
    elif num-1 > 5:
        test_type='plane'
        om_type='same_side'

#    print num-1,om_type
#    raw_input()
    number_of_tetras=0
    pyramidal_dihedrals=0
    all_dihedrals=[]
    min_dihedral=1000
    plane_found=[0,0,0,0]
    for i in range(1,num):
        for j in range(1,num):
            for k in range(1,num):
                for l in range(1,num):
                    if (i==j or i==k or i==l or j==k or j==l or k==l):
                       pass
                    else:
                        #all_dihedrals+=1
                        dihedral=abs(system.get_dihedral(i,j,k,l))
                        min_dihedral=min(min_dihedral,abs(dihedral))
                        all_dihedrals.append(dihedral)
                        if abs(dihedral-crit[test_type])< tol[test_type] or abs(dihedral-crit[test_type]+180)< tol[test_type]:
                            if test_type == 'tetrahedron':
                                pass
                            else:
                                check=True
                                if num-1>5:
                                    test['plane']=True
                                    plane_found=[i,j,k,l]
                                    check=False
                                    test[test_type]=True
                                    if check_if_plane_on_metal(0,[i,j,k,l],system):
                                        other_indeces=find_other_indeces([0,i,j,k,l],num)
                                        #check if other atoms are all in the same side of the plane
                                        dihedrals_other=[]
                                        for o_i in other_indeces:
                                            dihedrals_other.append(system.get_dihedral(j,k,l,o_i))
                                        if check_positive(dihedrals_other) and check_negative(dihedrals_other):
                                            pass
                                        else:
                                            check=True
                                if check:
                                    plane_found=[i,j,k,l]
                                    test[om_type]=True
                                    open_metal_mof=True
                        else:
                            if test_type == 'tetrahedron':
                                test[om_type]=True
                                open_metal_mof=True
                                indeces=[i,j,k,l]
                                n=4
                                for ii in range(0,n-2):
                                    for jj in range(ii+1,n-1):
                                        for kk in range(jj+1,n):
                                            iii=indeces[ii]
                                            jjj=indeces[jj]
                                            kkk=indeces[kk]
                                            dihedral=abs(system.get_dihedral(0,iii,jjj,kkk))
                                            if abs(dihedral-180)< 10 or abs(dihedral) < tol[test_type]:
                                                plane_found=[0,iii,jjj,kkk]
    #if num-1 ==4:
    #    print 'min_dihedral:',min_dihedral,open_metal_mof
    metal_cluster_with_adsorbate=[]
    if open_metal_mof:
        #if test['plane']:
        v1=system.cart_coords[plane_found[0]]
        v2=system.cart_coords[plane_found[1]]
        v3=system.cart_coords[plane_found[2]]
        v4=system.cart_coords[plane_found[3]]
        #elif test['non_TD']:
        #    v1=system.cart_coords[1]
        #    v2=system.cart_coords[2]
        #    v3=system.cart_coords[3]
        #    v4=system.cart_coords[4]
        metal_cluster_with_adsorbate=add_co2(v1,v2,v3,v4,system)
    return open_metal_mof,test,metal_cluster_with_adsorbate,min_dihedral,all_dihedrals

def add_co2(v1,v2,v3,v4,system):
    ads_dist = 2.2
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


    adsorption_site = []
    adsorption_pos = []
    co2_vector_O=find_adsorption_site(system,system.cart_coords[0],2.2)
    for i in range (1,2):
        #+r = generate_random_vector()
        #r = [i-j for i,j in zip(co2_vector_O, system.cart_coords[0])]
        #r_l = numpy.linalg.norm(numpy.array(r))
        #r = [i/r_l for i in r]
        #pos = list(rr*ads_dist for rr in r)
        #co2_vector_O = [i+j for i,j in zip(pos, system.cart_coords[0])]
        #pos = list(rr*1.16 for rr in r)
        #co2_vector_C = [i+j for i,j in zip(pos, co2_vector_O)]
        co2_vector_C =find_adsorption_site(system,co2_vector_O,1.16)
        r = [i-j for i,j in zip(co2_vector_C, co2_vector_O)]
        r_l = np.linalg.norm(np.array(r))
        r = [i/r_l for i in r]
        pos = list(rr*1.16 for rr in r)
        co2_vector_O_2 = [i+j for i,j in zip(pos, co2_vector_C)]

        p_avg_f = system.lattice.get_fractional_coords(co2_vector_O)
        co2_vector_C_f = system.lattice.get_fractional_coords(co2_vector_C)
        co2_vector_O_f = system.lattice.get_fractional_coords(co2_vector_O_2)

        #for e,c in zip(system.species,system.frac_coords):
        #    adsorption_site.append(e)
        #    adsorption_pos.append(c)
        adsorption_site.append('O')
        adsorption_pos.append(p_avg_f)
        adsorption_site.append('C')
        adsorption_pos.append(co2_vector_C_f)
        adsorption_site.append('O')
        adsorption_pos.append(co2_vector_O_f)

    return Structure(system.lattice,adsorption_site,adsorption_pos)

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

def check_if_plane_on_metal(m_i,indeces,system):
    crit=180
    #tol=12.5
    tol=25.0
    for i in range(1,len(indeces)):
        for j in range(1,len(indeces)):
            for k in range(1,len(indeces)):
                if (i==j or i==k or j==k):
                    pass
                else:
                    dihedral=abs(system.get_dihedral(m_i,indeces[i],indeces[j],indeces[k]))
                    if abs(dihedral-crit)< tol or abs(dihedral-crit+180)< tol:
                        #print dihedral
                        return True
    return False


def check_positive(N):
    for n in N:
        if n>0:return True

def check_negative(N):
    for n in N:
        if n<0:return True


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
        #print dist_before
        #print dist_after
        if dist_after>dist_before:
            for j in range(0,3):
                dist_vector_r[j]=np.rint(dist_vector[j])
            c_i_centered=c_i+dist_vector_r
            if dist_after>dist_before:
                c_i_centered=c_i
        system_centered.append(system.species[i],c_i_centered)
    return system_centered

def bond_check(ele1,ele2,dist,bond_tol):
    #bond = dist #- get_sum_of_cov_radii(ele1,ele2)
    up_bound = get_sum_of_cov_radii(ele1,ele2) + bond_tol
    low_bound = get_sum_of_cov_radii(ele1,ele2) - bond_tol
    #if bond < bond_tol and abs(dist) > 0:
    if dist < up_bound and dist > low_bound: # and abs(dist) > 0:
        return True
    else:
        return False

def find_adsorption_site(system,center,prob_dist):
    #find the adsorption site by maximizing the distance from all the atoms while
    #keep the distance fixed at some predefined distace
    tries=2500
    #prob_dist=2.2
    #center = str(system.species[0])
    new_position=0.0
    max_min_dist=0.0
    for i in range(0,tries):
        #probe_pos=generate_random_position(system.cart_coords[0],center,(prob_dist)-atom.get_vdf_radius(center_e))
        probe_pos=generate_random_position(center,prob_dist)
        probe_pos_f=system.lattice.get_fractional_coords(probe_pos)
        min_distance = min(system.lattice.get_all_distances(probe_pos_f,system.frac_coords)[0])
        if min_distance > max_min_dist:
            new_position = probe_pos
            max_min_dist = min_distance
    print max_min_dist
    adsorption_site=[]
    adsorption_pos=[]
    for e,c in zip(system.species,system.frac_coords):
        adsorption_site.append(e)
        adsorption_pos.append(c)
    adsorption_site.append('C')
    adsorption_pos.append(probe_pos_f)
    metal_cluster_with_adsorbate=Structure(system.lattice,adsorption_site,adsorption_pos)
    #return metal_cluster_with_adsorbate
    return new_position


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

def get_sum_of_cov_radii(ele1,ele2):
    return atom.get_covelent_radius(ele1)+atom.get_covelent_radius(ele2)


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
        dist=atom.get_vdf_radius(metal_element)+vdw_probe
        pos=generate_random_position(center,dist) #vdw_probe
        #pos=generate_random_position(center,metal_element,vdw_probe) #vdw_probe
        print >>params_file,'xx',pos[0],pos[1],pos[2]
        pos_f=system.lattice.get_fractional_coords(pos)
        if not check_for_overalp(center,pos_f,system,vdw_probe):
            count+=1
    #print count
    #raw_input()
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
        if (dist - (r_probe+atom.get_vdf_radius(str(system.species[i]))) ) < -1e-4:
            return True
    return False

def check_if_center(center,test_coords,system):
    distances=system.lattice.get_all_distances(test_coords,center)
    if distances[0]< 0.1:
        return True
    else:
        return False

def sphere_area(metal_element,probe_r):
    r=atom.get_vdf_radius(metal_element)
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


def get_bond_tolerance(ele1,ele2):
    if check_if_heavy_metal_bond(ele1,ele2):
        return 0.2
    else:
        return 0.4

def check_if_heavy_metal_bond(ele1,ele2):
    if is_heavy_metal(ele1) or is_heavy_metal(ele2):
        return True
    else:
        return False
def is_heavy_metal(ele): #\m/
    if atom.get_covelent_radius(ele) > 1.95:
         return True
    else:
        return False


def find_coordination_sequence(center,structure):
    '''computes the coordination sequence up to the Nth coordination shell
    as input it takes the MOF as a pymatgen Structure and the index of the center metal in the Structure
    '''
    shell_list = set([(center,(0,0,0))])
    shell_list_prev = set([])
    all_shells = set(shell_list)
    n_shells = 4
    cs = []
    print structure.species[center],structure.cart_coords[center][0],structure.cart_coords[center][1],structure.cart_coords[center][2]
    ele = [(str(structure.species[center]))]
    coords = [ [structure.frac_coords[center][0], structure.frac_coords[center][1], structure.frac_coords[center][2]]]
    coordination_structure = (Structure(structure.lattice,ele,coords))
    for n in range(0,n_shells):
        c_set = set([])
        for a_uc in shell_list:
            a = a_uc[0]
            lattice = a_uc[1]
            #print lattice
            coord_sphere = find_coord_sphere(a,structure)
            coord_sphere_with_uc = []
            for c in coord_sphere:
                dist = structure.lattice.get_all_distance_and_image(structure.frac_coords[a],structure.frac_coords[c])
                uc = lattice -min(dist)[1]
                coord_sphere_with_uc.append((c,tuple(uc)))
            #print n,a,coord_sphere
            coord_sphere_with_uc = tuple(coord_sphere_with_uc)
            #print type(coord_sphere_with_uc)
            #raw_input()
            c_set = c_set.union(set(coord_sphere_with_uc))
        #print 'a',c_set,shell_list_prev
        for a in shell_list_prev:
            c_set.discard(a)
        #print 'b',c_set,shell_list
        for a in shell_list:
            c_set.discard(a)
        #print 'c',c_set
        #raw_input()

        for i_uc in c_set:
            i = i_uc[0]
            #ele = (str(structure.species[i]))
            ele = atom.elements[n+3]
            coords = [structure.frac_coords[i][0], structure.frac_coords[i][1], structure.frac_coords[i][2]]
            coordination_structure.append(ele,coords)
            #print atom.elements[n],structure.cart_coords[i][0],structure.cart_coords[i][1],structure.cart_coords[i][2]
            #print atom.elements[n],structure.cart_coords[i][0],structure.cart_coords[i][1],structure.cart_coords[i][2]

        #cs.append(len(all_shells.union(c_set)) - len(all_shells))
        #print len(shell_list.union(c_set)),len(shell_list)
        #cs.append(len(shell_list.union(c_set)) - len(shell_list))
        cs.append(len(c_set))
        all_shells = all_shells.union(c_set)
        shell_list_prev = set(shell_list)
        shell_list = set(c_set)
    #coordination_structure = center_around_metal(coordination_structure)
    write_xyz_file('temp.xyz',coordination_structure)
    print cs
    return cs


if __name__=='__main__':
    main()
