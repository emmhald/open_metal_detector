# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %(Emmanuel Haldoupis)s
"""

import sys,os
path = os.path.expanduser('~/Dropbox/Work/python/modules')
if not path in sys.path:
    sys.path.insert(1, path)
del path
from pymatgen import Lattice,Structure
from pymatgen.io.cifio import CifParser
import pymatgen.io.smartio as sio
import pymatgen.io.vaspio.vasp_input as vasp
import shutil
from atomic_parameters import atoms
import matplotlib.pyplot as plt
import numpy as np


def main():
    #analyse_results('Mg')
    for m in atoms().metals:
        analyse_results(m)
#        analyse_results('U')
    #    print m
#    
    #make_plot()
    #analyse_results_sa()

def make_plot():
    summary= open('summary.out','r')
    
    frequency1=[]
    elements1=[]
    frequency2=[]
    elements2=[]
    group=15
    for l in summary:
        #print l
        f=int(l.split()[0])
        ele=l.split()[3]
        if f < group and f>0:
            print f,ele
            frequency1.append(f)
            elements1.append(ele)
        if f > group and f>0:
            print f,ele
            frequency2.append(f)
            elements2.append(ele)
    barPlot(elements1,frequency1,1)
    barPlot(elements2,frequency2,2)
    plt.show()
    
def barPlot(elements,frequency,fig):
    ind = np.arange(len(frequency))
    width = 0.5
    fig = plt.figure(fig)
    ax = plt.subplot(111)
    ax.set_xticks(ind+width)
    ax.set_xticklabels( elements )
    ax.bar(ind , frequency, width=width)
    

def analyse_results_sa():
    nbin=50
    summary_file= open('../detect_open_metal_sites/July_2014_runs/output_all/summary.out','r')
    sa_sum= open('sa_summary.txt','w')
    sa_f=[]
    sa=[]
    sa_f_o=[]
    sa_o=[]
    for line in summary_file:
        if 'yes' in line:
            sa_f_o.append(float(line.split(' ')[3]))
            sa_o.append(float(line.split(' ')[4]))
            print>>sa_sum,float(line.split(' ')[4])
        else :
            try:
                sa_f.append(float(line.split(' ')[3]))
                sa.append(float(line.split(' ')[4]))
            except:
                print 'no metal found'
    bins=[]
    for i in xrange(-5,250,5):
        bins.append(i)
    nbin=len(bins)
    sa_hist=np.histogram(sa,bins,density=True)
    sa_o_hist=np.histogram(sa_o,bins) #,density=True
    fig = plt.figure(1)
    ind=range(0,nbin-1)
    for i,b in zip(bins,sa_o_hist[0]):
        print i,b
    raw_input()
    #print len(sa_o_hist[0]),len(ind)
    
    width = 0.5
    ind_shift=[l+width for l in ind]
    #plt.bar(ind_shift, sa_hist[0],color='blue',width=width)
    plt.bar(ind, sa_o_hist[0],color='red',width=width)
    plt.show()


def analyse_results(element):

    folder='selected_mofs/'+element
    if not os.path.exists(folder):
        os.makedirs(folder)
    cif_folder='../detect_open_metal_sites/July_2014_runs/output_all/open_metal_mofs/'
    
    filetype='cif'
    params_file= open('../detect_open_metal_sites/July_2014_runs/output_all/open_metal_mofs.out','r')
    summary_file= open('../detect_open_metal_sites/July_2014_runs/output_all/summary.out','r')
    parameters= open(folder+'/parameters_'+element+'.txt','w')
    summary= open('summary.out','a')
    summary_ele= open(folder+'/summary.out','w')
    
    
    count=0
    for line in summary_file:
        if 'yes' in line:
            if ' '+element+' ' in line:
                count+=1
                struc=line.split(' ')[0]
                print struc
                if not os.path.exists(folder+'/'+struc):
                    os.makedirs(folder+'/'+struc)
                cif_name=struc+'.cif'
                print>>summary_ele,line.rstrip()
                cif=CifParser(cif_folder+cif_name)
                system=cif.get_structures()[0]
                sio.write_mol(system,folder+'/'+struc+'/'+struc+'.xyz')
                folder+'/summary.out'
                shutil.copyfile(cif_folder+cif_name, folder+'/'+struc+'/'+cif_name)
                dynamics=[]
                for atom in range(0,system.num_sites):
                    dynamics.append([False,False,False])
                pos=vasp.Poscar(system,selective_dynamics =dynamics)
                pos.write_file(folder+'/'+struc+'/POSCAR_'+struc)
                ads_folder='../detect_open_metal_sites/July_2014_runs/output/'+struc
                om_count=0
                if 1==2:
                    print 'what'
                    while True:
                        om_count+=1
                        ads_file=struc+'_first_coordination_sphere_with_ads'+str(om_count)+'.cif'
                        ads_path=ads_folder+'/'+ads_file
                        if not os.path.isfile(ads_path):
                            break
                        cif=CifParser(ads_path)
                        try:
                            system_ads=cif.get_structures()[0]
                        except:
                            print 'Cannot read ads cif file'
                            continue
                        shutil.copyfile(ads_path, folder+'/'+struc+'/'+ads_file)
    
                        dynamics=[]
                        for atom in range(0,system_ads.num_sites):
                            dynamics.append([True,True,True])
                        i=-1
                        for e,c in zip(system_ads.species,system_ads.frac_coords):
                            i+=1
                            for es,cs in zip(system.species,system.frac_coords):
                                dist=system.lattice.get_all_distances(c,cs)
                                if es == e and dist < 0.5:
                                    dynamics[i]=([False,False,False])
                                    
                        #for atom1 in range(0,system_ads[0].num_sites):
                        #    for atom2 in range(0,system[0].num_sites):
                        #        system.lattice.get_all_distances(m_f_coor,structure.frac_coords)
                        pos=vasp.Poscar(system_ads,selective_dynamics =dynamics)
                        pos.write_file(folder+'/'+struc+'/POSCAR_'+ads_file)
                
    print count,' MOFs with ',element,' open metal sites were found'
    print>>summary, count,' MOFs with ',element,' open metal sites were found'
    
    
if __name__=='__main__':
    main()    
    
'''
for line in params_file:
    struc=line.split(' ')[0]
    cif_name=struc+'.cif'
    cif=CifParser(cif_folder+cif_name)
    system=cif.get_structures()
    dynamics=[]
    if element in str(system[0].species):
        print struc
        print 'Found ',element
            shutil.copyfile(cif_folder+cif_name, folder+'/'+cif_name)
        for atom in range(0,system[0].num_sites):
            dynamics.append([False,False,False])
        pos=vasp.Poscar(system[0],selective_dynamics =dynamics)
        print>>parameters,struc,system[0].lattice.abc[0],system[0].lattice.abc[1],system[0].lattice.abc[2],system[0].lattice.angles[0],system[0].lattice.angles[1],system[0].lattice.angles[2]
        print struc,system[0].lattice.abc[0],system[0].lattice.abc[1],system[0].lattice.abc[2],system[0].lattice.angles[0],system[0].lattice.angles[1],system[0].lattice.angles[2]
        #raw_input()
        #print pos
        pos.write_file(folder+'/POSCAR_'+struc)
        sio.write_mol(system[0],folder+'/'+struc+'.xyz')
'''
