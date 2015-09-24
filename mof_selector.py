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
#import matplotlib.pyplot as plt
import numpy as np
import json
import re
from pymatgen.io.smartio import read_structure, write_structure
import argparse
from atomic_parameters import atoms as ap


def main():
    parser = argparse.ArgumentParser(description='Split file into batches')
    parser.add_argument('-p','--parameter_file', nargs='?',help='Name of parameters file', default = 'parameters.txt')
    parser.add_argument('-o','--output_folder', nargs='?',help='Name of output folder', default = 'output')
    parser.add_argument('-a','--analysis_folder', nargs='?',help='Name of analysis folder', default = 'analysis')
    parser.add_argument('-e','--target_element', nargs='?',help='Specify element to collect structures containing it', default = None)

    args = parser.parse_args()
    output_folder = args.output_folder
    analysis_folder = args.analysis_folder
    parameters_file = args.parameter_file
    ele = args.target_element
    json_dicts = load_structures(output_folder, parameters_file)
    analyze_for_tfac_using_json(json_dicts, analysis_folder)
    collect_statistics(json_dicts, analysis_folder)
    if ele:
        analyse_results(json_dicts, ele, analysis_folder, output_folder)

#    analyze_for_tfac()
#    for m in atoms().metals:
#        analyse_results(m)
#        analyse_results('U')
#        print m
#
    #make_plot()
    #analyse_results_sa()

def read_json(filename):

    try:
        json_filename = filename.split('.')[0]+'.JSON'
        json_dict = json.load(open(json_filename))
    except:
        try:
            json_filename = filename.split('.')[0]+'.json'
            json_dict = json.load(open(json_filename))
        except:
            return False

    return json_dict


def load_structures(output_folder, parameter_file):
    print('Reading structures...',end=' ')
    json_dicts=[]
    with open(parameter_file, 'r') as parameters:
        for l in parameters:
            struc = l.split('.')[0]
            json_dict = read_json(output_folder+'/'+struc+'/'+struc)
            #if not json_dict:
            #    continue
            if json_dict:
                json_dict['source_name'] = output_folder+'/'+struc
                json_dicts.append(json_dict)
    return json_dicts


def analyze_for_tfac_using_json(json_dicts, analysis_folder):
    tfac_analysis_folder = analysis_folder+'/tfac_analysis'
    make_folder(analysis_folder)
    make_folder(tfac_analysis_folder)

#    tfac_analysis_folder='analysis/tfac_analysis_redo_all'
#    tfac_analysis_folder='analysis/tfac_analysis_t5'
    yes_4 = open(tfac_analysis_folder+'/yes_4.out','w')
    no_4 = open(tfac_analysis_folder+'/no_4.out','w')
    yes_5 = open(tfac_analysis_folder+'/yes_5.out','w')
    no_5 = open(tfac_analysis_folder+'/no_5.out','w')
    yes_6 = open(tfac_analysis_folder+'/yes_6.out','w')
    no_6 = open(tfac_analysis_folder+'/no_6.out','w')
    #r = re.compile("([a-zA-Z]+)(-?(?:\d+())?(?:\.\d*())?(?:e-?\d+())?(?:\©|\1\3))")

#    with open(output_folder+'/summary.out','r') as summary:
    yes_or_no = dict()
    yes_or_no[True] = 'yes'
    yes_or_no[False] = 'no'
    #with open('parameters_'+'cif'+'.txt','r') as parameters:
    #    for l in parameters:
    #        struc = l.split('.')[0]
    #        print struc
    #        json_dict = read_json(output_folder+'/'+struc+'/'+struc)
    #        if not json_dict:
    #            continue

    for json_dict in json_dicts:
        struc = json_dict['material_name']
        num_of_ligands = fetch_num_of_ligands(json_dict)
        om_type = fetch_if_open(json_dict)
        tfactors = fetch_t_factor(json_dict)
        for L, typ, tf in zip(num_of_ligands,om_type,tfactors):
            if int(L) > 3 and int(L) < 7:
                outfile = yes_or_no[typ]+'_'+str(L)
                # print(outfile,struc,yes_or_no[typ],tf)
                print(struc, yes_or_no[typ], tf, file=eval(outfile))

    print('okay')
    yes_no_list = ['yes','no']
    for yn in yes_no_list:
        for i in range(4,7):
            outfile = yn+'_'+str(i)
            eval(outfile).close()
            datafile = open(tfac_analysis_folder+'/'+outfile+'.out','r')
            tfac_list=[]
            for l in datafile:
                tfac = l.split(' ')[2].rstrip('\n')
                tfac_list.append(float(tfac))
            eval(outfile).close()
            print(yn,i)
            if len(tfac_list) > 0:
                hist, edges = np.histogram(tfac_list, bins=50, range=(0, 1),
                                           density=True)
                hist_file = open(tfac_analysis_folder+'/'+outfile+'_hist.out','w')
                w=(edges[1]-edges[0])/2
                for e,h in zip(edges,hist):
                    print(e+w,h,file=hist_file)


def collect_statistics(json_dicts, analysis_folder):
    count_open = 0.0
    stats = dict()
    for json_dict in json_dicts:
        metal_sites=json_dict['metal_sites']
        for ms in metal_sites:
            metal = ms['metal']
            stats[metal] = dict()
            stats[metal]['count_open_sites'] = 0
            stats[metal]['count_open'] = 0
            stats[metal]['count'] = 0
            stats[metal]['t0'] = 0
            stats[metal]['t1'] = 0
            stats[metal]['t2'] = 0
            stats[metal]['t3'] = 0
            stats[metal]['t4'] = 0
            stats[metal]['t5'] = 0
            stats[metal]['t6'] = 0
            stats[metal]['tover'] = 0

    for json_dict in json_dicts:
        struc = json_dict['material_name']
        if json_dict['metal_sites_found']:
            count_open+=1
        metal_sites = json_dict['metal_sites']
        metals = set()
        metals_open = set()
        metals = set()
        number_of_linkers = set()
        for ms in metal_sites:
            metal = ms['metal']
            metals.add(metal)
            if ms['is_open'] :
                num_of_linkers=ms['number_of_linkers']
                stats[metal]['count_open_sites']+=1
                metals_open.add(metal)
                number_of_linkers.add(num_of_linkers)
        #stats[metal]['count'] = len(metals)
        #stats[metal]['count_open'] = len(metals_open)
        for m in metals:
            stats[m]['count']+=1
        for m in metals_open:
            stats[m]['count_open']+=1
        for l in number_of_linkers:
            #if l < 6 and l > 2:
            if l == 0:
                print(struc)
            if l > 6:
                tfrac='tover'
            else:
                tfac = 't'+str(l)
            stats[metal][tfac]+=1

        if len(number_of_linkers) > 0:
            if min(number_of_linkers) == 0 :
                dest=analysis_folder+'/L_'+str(min(number_of_linkers))+'/'
                copy_folder(dest , json_dict['source_name'])

    # This returns a sorted tuple based on keyfunc, which uses the key count_open to reverse sort the stats dictionary
    # A more consise but less clear sollution would be
    #stats_sorted = sorted(stats.items(), key = lambda tup : (-tup[1]["count_open"]))
    stats_sorted = sorted(stats.items(), key = keyfunc)
    printouts=[]
    for stat in stats_sorted:
        #Since stats_sorted is a sorted tuple of the dictionary stats,
        #the first element corresponds to the keys from the stas dictionary
        #and the second element to the value, in this dictionary holding the stats for each metal
        metal = stat[0]
        s = stat[1]

        if ap.is_group1(metal):
            continue

        printout = []
        printout.append(metal)
        printout.append(s['count'])
        printout.append(s['count_open'])
        printout.append(s['count_open_sites'])
        #print metal
        percent = 100*float(s['count_open'])/float(s['count'])
        percent_s = "{0:.2f} %".format(percent)
        printout.append(percent_s)
        #printout.append(1.1)
        printout.append(s['t0'])
        printout.append(s['t1'])
        printout.append(s['t2'])
        printout.append(s['t3'])
        printout.append(s['t3'])
        printout.append(s['t4'])
        printout.append(s['t5'])
        printout.append(s['t6'])
        printout.append(s['tover'])
        #string.center(s, wid)
        printouts.append(printout)

    print("Total MOFs:     ",len(json_dicts))
    print("Open Metal MOFs: {0:} {1:.2f} %".format(int(count_open),100.0*count_open/len(json_dicts)))
    titles=['Metal','All Found','Open-Metal','Open-Metal-Site','Per.-Open','l0','l1','l2','l3','l4','l5','l6','l-over-6']
    #print "{0:6}{1:12}{2:18}{3:6}{4:6}{5:6}".format(*titles)
    print("{0:6}{1:^12}{2:^18}{3:^12}{4:^12}{5:^6}{6:^6}{7:^6}{8:^6}{9:^6}{10:^6}{11:^8}{12:^8}".format(*titles))
    for p in printouts:
        #print "{0:6}{1:11}{2:15}{3:6}{4:6}{5:6}".format(*p)
        #print "{0:3}{1:6}{2:8}{3:>8}{4:6}{5:6}{6:6}{7:6}{8:6}{9:6}{10:6}{11:6}".format(*p)
        print("{0:6}{1:^12}{2:^18}{3:^12}{4:^16}{5:^6}{6:^6}{7:^6}{8:^6}{9:^6}{10:^6}{11:^8}{12:^8}".format(*p))


def copy_folder(dest, src):
    if not os.path.exists(dest):
        os.makedirs(dest)
    s = src.split('/')[-1]
    d = os.path.join(dest, s)
    #print(src, d, not os.path.exists(d))
    #input()
    if not os.path.exists(d):
        shutil.copytree(src, d)



def keyfunc(tup):
        key, d = tup
        return -d["count_open"]

def fetch_num_of_ligands(json_dict):
    num_of_ligands = []
    for ms in json_dict["metal_sites"]:
        num_of_ligands.append(ms["number_of_linkers"])
    return num_of_ligands

def fetch_if_open(json_dict):
    om_type = []
    for ms in json_dict["metal_sites"]:
        om_type.append(ms["is_open"])
    return om_type

def  fetch_t_factor(json_dict):
    t_factor = []
    for ms in json_dict["metal_sites"]:
        t_factor.append(ms["t_factor"])
    return t_factor


#older methods. conisder deleting

def analyze_for_tfac():
    tfac_analysis_folder = 'analysis/tfac_analysis_with_angles'
    tfac_analysis_folder = 'analysis/tfac_analysis_redo_all'
    tfac_analysis_folder = 'analysis/tfac_analysis_t5'


    summary = open('output/summary.out','r') #new_tol_first\
    summary = open('output/summary.out_first','r') ;filetype=1 #new_tol_first\
    summary = open('output/summary.out_5t','r') ;filetype=2 #new_tol_first\


    yes_4 = open(tfac_analysis_folder+'/yes_4.out','w')
    no_4 = open(tfac_analysis_folder+'/no_4.out','w')
    yes_5 = open(tfac_analysis_folder+'/yes_5.out','w')
    no_5 = open(tfac_analysis_folder+'/no_5.out','w')
    yes_6 = open(tfac_analysis_folder+'/yes_6.out','w')
    no_6 = open(tfac_analysis_folder+'/no_6.out','w')
    list_of_om_crit = ['plane','tetrahedron','non_tD','same_side','4_or_less']
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
                    print(outfile,struc,om.lstrip().split(',')[-1])
                    if filetype==2:
                        print(struc,om.lstrip().split(',')[-1],om.lstrip().split(',')[-2],file=eval(outfile))
                    else:
                        print(struc,om.lstrip().split(',')[-1],file=eval(outfile))

    yes_no_list=['yes','no']
    for yn in yes_no_list:
        for i in range(4,7):
            outfile=yn+'_'+str(i)
            eval(outfile).close()
            datafile=open(tfac_analysis_folder+'/'+outfile+'.out','r')
            tfac_list=[]
            for l in datafile:
                tfac=l.split(' ')[filetype].rstrip('\n')
                tfac_list.append(float(tfac))
            eval(outfile).close()
            hist, edges = np.histogram(tfac_list,bins=80,range=(0,1),density=True)
            hist_file=open(tfac_analysis_folder+'/'+outfile+'_hist.out','w')
            w=(edges[1]-edges[0])/2
            for e,h in zip(edges,hist):
                print>>hist_file,e+w,h


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
            print(f,ele)
            frequency1.append(f)
            elements1.append(ele)
        if f > group and f>0:
            print(f,ele)
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
                print('no metal found')
    bins=[]
    for i in xrange(-5,250,5):
        bins.append(i)
    nbin=len(bins)
    sa_hist = np.histogram(sa,bins,density=True)
    sa_o_hist = np.histogram(sa_o,bins) #,density=True
    fig = plt.figure(1)
    ind=range(0,nbin-1)
    for i,b in zip(bins,sa_o_hist[0]):
        print(i,b)
    raw_input()
    #print len(sa_o_hist[0]),len(ind)

    width = 0.5
    ind_shift=[l+width for l in ind]
    #plt.bar(ind_shift, sa_hist[0],color='blue',width=width)
#    plt.bar(ind, sa_o_hist[0],color='red',width=width)
#    plt.show()


def analyse_results(json_dicts, element, analysis_folder, output_folder):

    folder = 'analysis/selected_mofs/'+element
    cif_folder_out = 'analysis/selected_mofs/'+element+'/cif_files'
    if not os.path.exists(folder):
        os.makedirs(folder)
    if not os.path.exists(cif_folder_out):
        os.makedirs(cif_folder_out)
    summary_ele = open(folder+'/summary.out','w')
    cif_folder = analysis_folder+'/open_metal_mofs/'
    summary = open('summary.out','a')

    print('MOF','#OMS','#OMS_types','OMS ids','#OMS_per_type')
    count=0
    for json_dict in json_dicts:
        oms_ids = []
        struc_contains_open_metal = False
        metal_sites = json_dict['metal_sites']
        oms_found = False
        contains_metal_but_closed = False
        for ms in metal_sites:
            if element in ms['metal'] and not ms['is_open'] :
                contains_metal_but_closed = True
            if ms['is_open'] and element in ms['metal']:
                struc_contains_open_metal = True
            if 'oms_id' in ms:
                oms_ids.append(ms['oms_id'])
        if contains_metal_but_closed:
            print(json_dict['material_name'],'closed')
        if struc_contains_open_metal:
            count += 1

            struc = json_dict['material_name']
            ads_folder = output_folder+'/'+struc
            copy_folder(folder+'/' , ads_folder)
            #copy_folder(folder+'/'+struc , ads_folder)
            if not os.path.exists(folder+'/'+struc):
                os.makedirs(folder+'/'+struc)
            if 1==2:

                cif_name = struc+'.cif'
                print(cif_name,'-',len(oms_ids),'-',len(set(oms_ids)),'-', oms_ids,'-', *[oms_ids.count(i) for i in set(oms_ids)] )
                print(struc+'.cif',file=summary_ele)
                cif = CifParser(cif_folder+cif_name)
                system = cif.get_structures()[0]
                system.to(fmt='xyz',filename = folder+'/'+struc+'/'+struc+'.xyz')
                shutil.copyfile(cif_folder+cif_name, folder+'/'+struc+'/'+cif_name)
                shutil.copyfile(cif_folder+cif_name, cif_folder_out+'/'+cif_name)
                dynamics=[]
                for atom in range(0,system.num_sites):
                    dynamics.append([False,False,False])
                pos=vasp.Poscar(system,selective_dynamics =dynamics)
                pos.write_file(folder+'/'+struc+'/POSCAR_'+struc)
                #print(folder+'/'+struc+'/POSCAR_'+struc)
                #ads_folder='../detect_open_metal_sites/July_2014_runs/output/'+struc


            if 1 == 2:
                for ads_type in ['co2', 'n2']:
                    om_count = -1
                    while True:
                        om_count+=1
                        ads_file = struc+'_first_coordination_sphere_with_'+ads_type+'_'+str(om_count)+'.cif'
                        ads_path = ads_folder+'/'+ads_file
                        if not os.path.isfile(ads_path):
                            break
                        cif=CifParser(ads_path)
                        try:
                            system_ads=cif.get_structures()[0]
                        except:
                            print('Cannot read ads cif file')
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
                        pos=vasp.Poscar(system_ads,selective_dynamics =dynamics)
                        pos.write_file(folder+'/'+struc+'/POSCAR_'+ads_file)

    print(count,' MOFs with ',element,' open metal sites were found')
    print(count,' MOFs with ',element,' open metal sites were found',file=summary)



def make_folder(folder):
    if not os.path.exists(folder):
        os.makedirs(folder)

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
