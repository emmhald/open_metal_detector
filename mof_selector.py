# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %(Emmanuel Haldoupis)s
"""

import sys,os
path = os.path.expanduser('~/Dropbox/Work/python/modules')
path = os.path.expanduser('~/dropbox/Work/python/modules')
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
    #analyse_results('Fe')
    analyze_for_td()
#    for m in atoms().metals:
#        analyse_results(m)
#        analyse_results('U')
    #    print m
#
    #make_plot()
    #analyse_results_sa()

def analyze_for_td():
    td_analysis_folder='analysis/td_analysis_with_angles'
    td_analysis_folder='analysis/td_analysis_redo_all'
    td_analysis_folder='analysis/td_analysis_t5'


    summary= open('output/summary.out','r') #new_tol_first\
    summary= open('output/summary.out_first','r') ;filetype=1 #new_tol_first\
    summary= open('output/summary.out_5t','r') ;filetype=2 #new_tol_first\


    yes_4= open(td_analysis_folder+'/yes_4.out','w')
    no_4= open(td_analysis_folder+'/no_4.out','w')
    yes_5= open(td_analysis_folder+'/yes_5.out','w')
    no_5= open(td_analysis_folder+'/no_5.out','w')
    yes_6= open(td_analysis_folder+'/yes_6.out','w')
    no_6= open(td_analysis_folder+'/no_6.out','w')
    list_of_om_crit=['plane','tetrahedron','non_TD','same_side','4_or_less']
    for l in summary:
        struc=l.split(' ')[0]
        if 'yes' in l:
            yes_or_no='yes'
        else:
            yes_or_no='no'
        open_metals=l.split(yes_or_no)[1].rstrip('\n').rstrip().lstrip().split('  ')
        for om in open_metals:
            if ',' in om:
                L=om.lstrip().split(',')[1]
                L=L.split('L')[0]
                is_open_check=om.lstrip().split(',')[2]
               # print struc,is_open_check,is_open_check in list_of_om_crit
               # raw_input()
                if is_open_check in list_of_om_crit:
                    yes_or_no='yes'
                else:
                    yes_or_no='no'
        #        if struc == 'JUTCUW_clean':
        #            print yes_or_no
        #            raw_input()
                if int(L) > 3 and int(L) < 7:
                    outfile=yes_or_no+'_'+str(L)
                    print outfile,struc,om.lstrip().split(',')[-1]
                    if filetype==2:
                        print>>eval(outfile),struc,om.lstrip().split(',')[-1],om.lstrip().split(',')[-2]
                    else:
                        print>>eval(outfile),struc,om.lstrip().split(',')[-1]

    yes_no_list=['yes','no']
    for yn in yes_no_list:
        for i in range(4,7):
            outfile=yn+'_'+str(i)
            eval(outfile).close()
            datafile=open(td_analysis_folder+'/'+outfile+'.out','r')
            td_list=[]
            for l in datafile:
                td=l.split(' ')[filetype].rstrip('\n')
                td_list.append(float(td))
            eval(outfile).close()
            hist, edges = np.histogram(td_list,bins=40,range=(0,1),density=True)
            hist_file=open(td_analysis_folder+'/'+outfile+'_hist.out','w')
            w=(edges[1]-edges[0])/2
            for e,h in zip(edges,hist):
                print>>hist_file,e+w,h


#    outfile=yes_or_no+'_'+str(L)
#    data_file=open(outfile,r)


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

    folder='analysis/selected_mofs/'+element
    cif_folder_out='analysis/selected_mofs/'+element+'/cif_files'
    if not os.path.exists(folder):
        os.makedirs(folder)
    if not os.path.exists(cif_folder_out):
        os.makedirs(cif_folder_out)
    cif_folder='output_all/open_metal_mofs/'
    #cif_folder='../detect_open_metal_sites/July_2014_runs/output_all/open_metal_mofs/'

    filetype='cif'
    #params_file= open('../detect_open_metal_sites/July_2014_runs/output_all/open_metal_mofs.out','r')
    params_file= open('output_all/open_metal_mofs.out','r')
    #summary_file= open('../detect_open_metal_sites/July_2014_runs/output_all/summary.out','r')
    summary_file= open('output_all/summary.out','r')
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
                shutil.copyfile(cif_folder+cif_name, cif_folder_out+'/'+cif_name)
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
