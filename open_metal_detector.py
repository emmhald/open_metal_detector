# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %(Emmanuel Haldoupis)s
"""
import sys,os
#sys.path.append('~/Dropbox/Work/python/modules')
path = os.path.expanduser('~/Dropbox/Work/python/modules')
if not path in sys.path:
    sys.path.insert(1, path)
del path
from xyz_read import xyz_file
from pymatgen import Lattice,Structure
from pymatgen.io.cifio import CifParser
import numpy as np
import shlex
import shutil
import time
import math
import random
from atomic_parameters import atoms

atom=atoms()
target_folder='CORE-MOF-DB-June2014/'
filetype='cif'
number_of_structures=5107
def main():
    t0 = time.time()
    params_file= open('parameters_'+filetype+'.txt','r')
    count=1
    clear_files()
    for struc in params_file:
        t0s = time.time()
        count+=1
        uc_params=[]
        line_elements=shlex.split(struc)
        filename=line_elements[0]
        if filetype=='xyz':
            for j in range(1,7):
                uc_params.append(float(line_elements[j]))
        analyze_structure(filename,uc_params)
        t1s = time.time()
        print 'Time:',t1s-t0s
        if count >number_of_structures:
            break
    params_file.close()
    t1 = time.time()
    print 'Total Time',t1-t0

def clear_files():
    make_folder('output')
    open_metal_mofs= open('output/open_metal_mofs.out','w')
    problematic_mofs= open('output/problematic.out','w')
    summary_mofs= open('output/summary.out','w')
    open_metal_mofs.close()
    problematic_mofs.close()
    summary_mofs.close()

def make_folder(folder):
    if not os.path.exists(folder):
        os.makedirs(folder)

def analyze_structure(filename,uc_params):

    mof_name=filename.split('.')[0]
    output_folder='output/'+mof_name

    make_folder(output_folder)
    make_folder('output/open_metal_mofs')
    make_folder('output/problematic_metal_mofs')

    open_metal_mofs= open('output/open_metal_mofs.out','a')
    problematic_mofs= open('output/problematic.out','a')
    summary_mofs= open('output/summary.out','a')

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
    m_sa_frac,m_surface_area=0.0,0.0 #get_metal_surface_areas(metal,system)

    open_metal_site=False
    problematic_structure=False
    test=[] #dict()
    summary=str(m_sa_frac)+' '+str(m_surface_area)+' '

    count_omsites=0
    for m,open_metal_candidate in enumerate(first_coordnation_structure_each_metal):
        op,pr,t=check_if_open(open_metal_candidate)
        if op:
            count_omsites+=1
            op_with_adsorbate=find_adsorption_site(open_metal_candidate)
            xyz_op_with_adsorbate=xyz_file()
            xyz_op_with_adsorbate.make_xyz(output_folder+'/'+mof_name+'_first_coordination_sphere_with_ads'+str(count_omsites)+'.xyz',op_with_adsorbate.cart_coords,op_with_adsorbate.species)
            open_metal_site=op
            problematic_structure=pr
            test.append(t)
            if not 'yes' in summary:
                summary=summary+' '+'yes'
            string=str(open_metal_candidate.species[0])+' '+str(open_metal_candidate.num_sites-1)+'L'
            for ti in t:
                if t[ti]:
                    string=string+' '+str(ti)
            if not string in summary:
                summary=summary+' '+string
        xyz_first_coord_sphere=xyz_file()
        xyz_first_coord_sphere.make_xyz(output_folder+'/'+mof_name+'_first_coordination_sphere'+str(m)+'.xyz',open_metal_candidate.cart_coords,open_metal_candidate.species)
        
    
    summary=mof_name+' '+str(count_omsites)+' '+str(count_omsites/system.volume)+' '+summary
    

    if open_metal_site:
        print>>open_metal_mofs,summary
        shutil.copyfile(target_folder+filename, 'output/open_metal_mofs/'+filename)
    else:
        summary=summary+' '+'no'
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


def write_xyz_file(filename,system):
    xyz=xyz_file()
    xyz.make_xyz(filename,system.cart_coords,system.species)

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

        increase=1.0
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
            if check_if_enough_bonds(bonds_found,increase,str(metal.species[m])):
                break
            else:
                increase+=0.1
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

def check_if_6_or_more(system):
    if system.num_sites > 6:
        return True
    else:
        return False


def check_if_open(system):
#def check_if_pyramidal_square(system):
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
    if atom.is_lanthanide_or_actinide(str(system.species[0])):
        max_cordination=5
    if num-1 <= max_cordination:
        open_metal_mof=True
        problematic=True
        test['4_or_less']=True
        return open_metal_mof,problematic,test

    open_metal_mof,test=check_non_metal_dihedrals(system,test)
    if num-1 == 4 and not open_metal_mof:
        open_metal_mof,test=check_metal_dihedrals(system,test)

    return open_metal_mof,problematic,test


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
                         dihedral_other1=abs(system.get_dihedral(other_index[0],j,k,l))
                         dihedral_other2=abs(system.get_dihedral(other_index[1],j,k,l))
                         if (dihedral_other1*dihedral_other2)  >= 0:
                             open_metal_mof=True
                             test['metal_plane']=True

    if number_of_planes == 4:
        if open_metal_mof:
            print 'conflicting criteria'
            raw_input()
        else:
            open_metal_mof=False
    return open_metal_mof,test

def check_non_metal_dihedrals(system,test):
    open_metal_mof=False
    crit=dict()
    tol=dict()
    crit['plane']=180
    tol['plane']=30
    crit['tetrahedron']=70
    tol['tetrahedron']=10

    number_of_tetras=0
    all_dihedrals=0
    pyramidal_dihedrals=0

    all_dihedrals=0
    num=system.num_sites
    for i in range(1,num):
        for j in range(1,num):
            for k in range(1,num):
                for l in range(1,num):
                    if (i==j or i==k or i==l or j==k or j==l or k==l):
                       pass
                    else:
                        #all_dihedrals+=1
                        dihedral=abs(system.get_dihedral(i,j,k,l))
                        if num-1==4:
                            if abs(dihedral-crit['tetrahedron'])< tol['tetrahedron'] or abs(dihedral-crit['tetrahedron']+180)< tol['tetrahedron']:
                                pass #pyramidal_dihedrals+=1
                            else:
                                test['non_TD']=True
                                open_metal_mof=True
                        if num-1==5:
                            if abs(dihedral-crit['plane'])< tol['plane'] or abs(dihedral-crit['plane']+180)< tol['plane']:
                                open_metal_mof=True
                                test['plane']=True
                        if num-1>5:
                            if abs(dihedral-crit['plane'])< tol['plane'] or abs(dihedral-crit['plane']+180)< tol['plane']:
                                test['plane']=True
                                if check_if_plane_on_metal(0,[i,j,k,l],system):
                                    other_indeces=find_other_indeces([0,i,j,k,l],num)
                                    #check if other atoms are all in the same side of the plane
                                    dihedrals_other=[]
                                    for o_i in other_indeces:
                                        dihedrals_other.append(system.get_dihedral(j,k,l,o_i))
                                    if check_positive(dihedrals_other) and check_negative(dihedrals_other):
                                        pass
                                    else:
                                        open_metal_mof=True
                                        test['same_side']=True

    return open_metal_mof,test


def check_if_plane_on_metal(m_i,indeces,system):
    crit=180
    tol=12.5
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
    if dist-get_sum_of_cov_radii(ele1,ele2) < bond_tol:
        return True
    else:
        return False

def find_adsorption_site(system):
    #find the adsorption site by maximizing the distance from all the atoms while
    #keep the distance fixed at some predefined distace
    tries=25000
    min_dist=0.0
    prob_dist=2.2
    metal=str(system.species[0])
    new_position=0.0
    for i in range(0,tries):
        probe_pos=generate_random_position(system.cart_coords[0],metal,(prob_dist)-atom.get_vdf_radius(metal)) 
        probe_pos_f=system.lattice.get_fractional_coords(probe_pos)
        sum_distance=0.0
        for e,c in zip(system.species,system.frac_coords):
#            print e,c
            sum_distance+=(system.lattice.get_all_distances(probe_pos_f,c))
        if sum_distance > min_dist:
            new_position=probe_pos
    print sum_distance
    adsorption_site=[]
    adsorption_pos=[]
    for e,c in zip(system.species,system.frac_coords):
        adsorption_site.append(e)
        adsorption_pos.append(c)
    adsorption_site.append('C')
    adsorption_pos.append(probe_pos_f)
    metal_cluster_with_adsorbate=Structure(system.lattice,adsorption_site,adsorption_pos)
    return metal_cluster_with_adsorbate



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
        pos=generate_random_position(center,metal_element,vdw_probe) #vdw_probe
        print >>params_file,'xx',pos[0],pos[1],pos[2]
        pos_f=system.lattice.get_fractional_coords(pos)
        if not check_for_overalp(center,pos_f,system,vdw_probe):
            count+=1
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

def generate_random_position(center,metal,vdw_probe):
    r=generate_random_vector()
    dist=atom.get_vdf_radius(metal)+vdw_probe
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



if __name__=='__main__':
    main()

