# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %(Emmanuel Haldoupis)s
"""

import sys
import os
path = os.path.expanduser('~/Dropbox/Work/python/modules')
path = os.path.expanduser('~/dropbox/Work/python/modules')
if path not in sys.path:
    sys.path.insert(1, path)
del path
from pymatgen import Lattice,Structure
from pymatgen.io.cif import CifParser
import pymatgen.io.smart as sio
import pymatgen.io.vasp.inputs as vasp
import shutil
from atomic_parameters import atoms
#import matplotlib.pyplot as plt
import numpy as np
import json
import re
from pymatgen.io.smart import read_structure, write_structure
import argparse
from atomic_parameters import atoms as ap


def main():
    parser = argparse.ArgumentParser(description='Split file into batches')
    parser.add_argument('-p', '--parameter_file', nargs='?',
                        help='Name of parameters file',
                        default='parameters.txt')
    parser.add_argument('-o', '--output_folder', nargs='?',
                        help='Name of output folder',
                        default='output')
    parser.add_argument('-a', '--analysis_folder',
                        nargs='?', help='Name of analysis folder',
                        default='analysis')
    parser.add_argument('-e', '--target_element',
                        nargs='?',
                        help='Specify element to collect '
                             'structures containing it',
                        default=None)

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


def read_json(filename):

    try:
        json_filename = filename+'.JSON'
        json_dict = json.load(open(json_filename))
    except:
        try:
            json_filename = filename+'.json'
            json_dict = json.load(open(json_filename))
        except:
            return False

    return json_dict


def load_structures(output_folder, parameter_file):
    print('Reading structures...', end=' ')
    json_dicts = []
    with open(parameter_file, 'r') as parameters:
        for l in parameters:
            filetype = l.split('.')[-1]
            struc = l.split('.'+filetype)[0]
            #struc = l.split('.')[0]
            json_file = output_folder+'/'+struc+'/'+struc
            json_dict = read_json(json_file)
            # if not json_dict:
            #    continue
            if isinstance(json_dict, dict):
                json_dict['source_name'] = output_folder+'/'+struc
                json_dicts.append(json_dict)
    return json_dicts


def analyze_for_tfac_using_json(json_dicts, analysis_folder):
    tfac_analysis_folder = analysis_folder+'/tfac_analysis'
    make_folder(analysis_folder)
    make_folder(tfac_analysis_folder)

#    tfac_analysis_folder='analysis/tfac_analysis_redo_all'
#    tfac_analysis_folder='analysis/tfac_analysis_t5'
    yes_4 = open(tfac_analysis_folder+'/yes_4.out', 'w')
    no_4 = open(tfac_analysis_folder+'/no_4.out', 'w')
    yes_5 = open(tfac_analysis_folder+'/yes_5.out', 'w')
    no_5 = open(tfac_analysis_folder+'/no_5.out', 'w')
    yes_6 = open(tfac_analysis_folder+'/yes_6.out', 'w')
    no_6 = open(tfac_analysis_folder+'/no_6.out', 'w')
    # r = re.compile("([a-zA-Z]+)(-?(?:\d+())?(?:\.\d*())?(?:e-?\d+())?(?:\©|\1\3))")

    # with open(output_folder+'/summary.out','r') as summary:
    yes_or_no = dict()
    yes_or_no[True] = 'yes'
    yes_or_no[False] = 'no'
    # with open('parameters_'+'cif'+'.txt','r') as parameters:
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
        for L, typ, tf in zip(num_of_ligands, om_type, tfactors):
            if 7 > int(L) > 3:  #
                outfile = yes_or_no[typ]+'_'+str(L)
                # print(outfile,struc,yes_or_no[typ],tf)
                print(struc, yes_or_no[typ], tf, file=eval(outfile))

    print('okay')
    yes_no_list = ['yes', 'no']
    for yn in yes_no_list:
        for i in range(4, 7):
            outfile = yn+'_'+str(i)
            eval(outfile).close()
            datafile = open(tfac_analysis_folder+'/'+outfile+'.out', 'r')
            tfac_list = []
            for l in datafile:
                tfac = l.split(' ')[2].rstrip('\n')
                tfac_list.append(float(tfac))
            eval(outfile).close()
            print(yn, i)
            if len(tfac_list) > 0:
                hist, edges = np.histogram(tfac_list, bins=50, range=(0, 1),
                                           density=True)
                hist_filename = tfac_analysis_folder+'/'+outfile+'_hist.out'
                with open(hist_filename, 'w') as hist_file:
                    w = (edges[1]-edges[0])/2
                    for e, h in zip(edges, hist):
                        print(e+w, h, file=hist_file)
                hist, edges = np.histogram(tfac_list, bins=50, range=(0, 1),
                                           density=False)
                folder_outfile = tfac_analysis_folder + '/' + outfile
                hist_filename = folder_outfile + '_hist_abs.out'
                with open(hist_filename, 'w') as hist_file:
                    w = (edges[1]-edges[0])/2
                    for e, h in zip(edges, hist):
                        print(e+w, h, file=hist_file)


def collect_statistics(json_dicts, analysis_folder):
    count_open = 0.0
    stats = dict()
    for json_dict in json_dicts:
        metal_sites = json_dict['metal_sites']
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
            count_open += 1
        metal_sites = json_dict['metal_sites']
        metals_open = set()
        metals = set()
        oms_set = set()
        number_of_linkers = set()
        for ms in metal_sites:
            metal = ms['metal']
            metals.add(metal)
            if ms['is_open']:
                num_of_linkers = ms['number_of_linkers']
                stats[metal]['count_open_sites'] += 1
                metals_open.add(metal)
                number_of_linkers.add(num_of_linkers)
                oms_set.add((metal, num_of_linkers))
        # stats[metal]['count'] = len(metals)
        # stats[metal]['count_open'] = len(metals_open)
        for m in metals:
            stats[m]['count'] += 1
        for om in oms_set:
            met, l = om
            if l > 6:
                tfac = 'tover'
            else:
                tfac = 't'+str(l)
            # stats[met]['count_open'] += 1
            stats[met][tfac] += 1
        for m in metals_open:
            stats[m]['count_open'] += 1
        #for l in number_of_linkers:
        #    # if l == 0:
        #    #     print(struc)
        #    if l > 6:
        #        tfac = 'tover'
        #    else:
        #        tfac = 't'+str(l)
        #    stats[metal][tfac] += 1

        if len(number_of_linkers) > 0:
            # if 6 in list(number_of_linkers) and 'Zr' in metals_open:
            if 5 in list(number_of_linkers) and len(metals_open) != len(metals):
                dest = analysis_folder+'/L_'+str(min(number_of_linkers))+'/'
                copy_folder(dest, json_dict['source_name'])

    # This returns a sorted tuple based on keyfunc,
    # which uses the key count_open to reverse sort the stats dictionary
    # A more consise but less clear sollution would be
    # stats_sorted = sorted(stats.items(), key = lambda tup :
    # (-tup[1]["count_open"]))
    stats_sorted = sorted(stats.items(), key=keyfunc)
    printouts = []
    for stat in stats_sorted:
        # Since stats_sorted is a sorted tuple of the dictionary stats,
        # the first element corresponds to the keys from the stas dictionary
        # and the second element to the value, in this dictionary holding
        # the stats for each metal
        metal = stat[0]
        s = stat[1]

        if ap.is_group1(metal):
            continue
        if not ap.is_d4_or_less(metal):
            continue
        percent = 100*float(s['count_open'])/float(s['count'])
        percent_s = "{0:.2f} %".format(percent)
        printout = [metal,
                    s['count'],
                    s['count_open'],
                    s['count_open_sites'],
                    percent_s,
                    s['t0'],
                    s['t1'],
                    s['t2'],
                    s['t3'],
                    s['t4'],
                    s['t5'],
                    s['t6'],
                    s['tover']]
        printouts.append(printout)

    print("Total MOFs:     ", len(json_dicts))
    if len(json_dicts) == 0:
        return
    print("Open Metal MOFs: {0:} {1:.2f} %"
          .format(int(count_open), 100.0 * count_open/len(json_dicts)))
    titles = ['Metal', 'All-Found', 'Open-Metal', 'Open-Metal-Site',
              'Per.-Open', 'l0', 'l1', 'l2', 'l3', 'l4', 'l5', 'l6',
              'l-over-6']
    # print "{0:6}{1:12}{2:18}{3:6}{4:6}{5:6}".format(*titles)
    print("{0:6}{1:^12}{2:^18}{3:^12}{4:^12}{5:^6}{6:^6}{7:^6}"
          "{8:^6}{9:^6}{10:^6}{11:^8}{12:^8}".format(*titles))
    for p in printouts:
        # print "{0:6}{1:11}{2:15}{3:6}{4:6}{5:6}".format(*p)
        # print "{0:3}{1:6}{2:8}{3:>8}{4:6}{5:6}{6:6}{7:6}{8:6}{9:6}{10:6}
        # {11:6}".format(*p)
        print("{0:6}{1:^12}{2:^18}{3:^12}{4:^16}{5:^6}{6:^6}{7:^6}"
              "{8:^6}{9:^6}{10:^6}{11:^8}{12:^8}".format(*p))


def copy_folder(dest, src):
    if not os.path.exists(dest):
        os.makedirs(dest)
    s = src.split('/')[-1]
    d = os.path.join(dest, s)
    # print(src, d, not os.path.exists(d))
    # input()
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


def fetch_t_factor(json_dict):
    t_factor = []
    for ms in json_dict["metal_sites"]:
        t_factor.append(ms["t_factor"])
    return t_factor

# older methods. conisder deleting


def make_plot():
    summary = open('summary.out', 'r')

    frequency1 = []
    elements1 = []
    frequency2 = []
    elements2 = []
    group = 15
    for l in summary:
        f = int(l.split()[0])
        ele = l.split()[3]
        if 0 < f < group:
            print(f, ele)
            frequency1.append(f)
            elements1.append(ele)
        if 0 < f > group:
            print(f, ele)
            frequency2.append(f)
            elements2.append(ele)
    barplot(elements1, frequency1, 1)
    barplot(elements2, frequency2, 2)
    plt.show()


def barplot(elements, frequency, figu):
    ind = np.arange(len(frequency))
    width = 0.5
    fig = plt.figure(figu)
    ax = plt.subplot(111)
    ax.set_xticks(ind+width)
    ax.set_xticklabels(elements)
    ax.bar(ind, frequency, width=width)


def analyse_results_sa():
    nbin = 50
    summary_file = open('../detect_open_metal_sites/July_2014_runs'
                        '/output_all/summary.out', 'r')
    sa_sum = open('sa_summary.txt', 'w')
    sa_f = []
    sa = []
    sa_f_o = []
    sa_o = []
    for line in summary_file:
        if 'yes' in line:
            sa_f_o.append(float(line.split(' ')[3]))
            sa_o.append(float(line.split(' ')[4]))
            print(sa_sum,float(line.split(' ')[4]))
        else :
            try:
                sa_f.append(float(line.split(' ')[3]))
                sa.append(float(line.split(' ')[4]))
            except:
                print('no metal found')
    bins = []
    for i in range(-5, 250, 5):
        bins.append(i)
    nbin = len(bins)
    sa_hist = np.histogram(sa, bins, density=True)
    sa_o_hist = np.histogram(sa_o, bins)  # ,density=True
    fig = plt.figure(1)
    ind = range(0, nbin-1)
    for i, b in zip(bins, sa_o_hist[0]):
        print(i, b)
    input()
    # print len(sa_o_hist[0]),len(ind)

    width = 0.5
    ind_shift = [l+width for l in ind]
    # plt.bar(ind_shift, sa_hist[0],color='blue',width=width)
    # plt.bar(ind, sa_o_hist[0],color='red',width=width)
    # plt.show()


def analyse_results(json_dicts, element, analysis_folder, output_folder):

    folder = analysis_folder+'/selected_mofs/'+element
    vaspfolder = analysis_folder+'/calcs/'
    cif_folder_out = analysis_folder+'/selected_mofs/'+element+'/cif_files'
    make_folder(folder)
    make_folder(vaspfolder)
    make_folder(cif_folder_out)

    open(folder+'/summary.out', 'w').close()
    cif_folder = output_folder+'/open_metal_mofs/'
    open('summary.out', 'w').close()

    print('MOF', '#OMS', '#OMS_types', 'OMS ids', '#OMS_per_type')
    count = 0
    for json_dict in json_dicts:
        oms_ids = []
        struc_contains_open_metal = False
        metal_sites = json_dict['metal_sites']
        oms_found = False
        contains_metal_but_closed = False
        for ms in metal_sites:
            if element in ms['metal'] and not ms['is_open']:
                contains_metal_but_closed = True
            if ms['is_open'] and element in ms['metal']:
                struc_contains_open_metal = True
            if 'oms_id' in ms:
                oms_ids.append(ms['oms_id'])
        if contains_metal_but_closed and not struc_contains_open_metal:
            print(json_dict['material_name'], 'closed')
        if struc_contains_open_metal:
            count += 1

            struc = json_dict['material_name']
            ads_folder = output_folder+'/'+struc
            copy_folder(folder+'/', ads_folder)
            # copy_folder(folder+'/'+struc , ads_folder)
            make_folder(folder+'/'+struc)

            cif_name = struc+'.cif'
            with open(folder+'/summary.out', 'a') as summary_ele:
                print(struc+'.cif', file=summary_ele)
            cif = CifParser(cif_folder+cif_name)
            system = cif.get_structures()[0]
            struc_xyz_ = folder + '/' + struc + '/' + struc + '.xyz'
            system.to(fmt='xyz', filename=struc_xyz_)

            name = folder + '/' + struc + '/' + cif_name
            shutil.copyfile(cif_folder+cif_name, name)
            out_cif_name = cif_folder_out + '/' + cif_name
            shutil.copyfile(cif_folder+cif_name, out_cif_name)

            dynamics = []
            for atom in range(0, system.num_sites):
                dynamics.append([False, False, False])
            pos = vasp.Poscar(system, selective_dynamics=dynamics)
            pos.write_file(folder+'/'+struc+'/POSCAR_'+struc)

            make_folder(vaspfolder+'/'+struc)
            make_folder(vaspfolder+'/'+struc+'/framework')
            pos.write_file(vaspfolder+'/'+struc+'/framework/POSCAR')

            potcar_file = vaspfolder+'/'+struc+'/framework/POTCAR'
            make_potcar(pos, potcar_file)
            copy_vasp_files('vasp_input_files',
                            vaspfolder+'/'+struc+'/framework',
                            'INCAR_opt', struc)
            # replace_in_file('name', struc, 'vasp_input_files/pbs.multi')


            ads_folder = output_folder+'/'+struc
            om_count = dict()
            adss = ['co2', 'n2']
            for ads_type in adss:
                om_count[ads_type] = -1
                while True:
                    with__ = struc+'_first_coordination_sphere_with_'+ads_type
                    ads_file = with__+'_'+str(om_count[ads_type]+1)+'.cif'
                    ads_path = ads_folder+'/'+ads_file
                    if not os.path.isfile(ads_path):
                        break
                    om_count[ads_type] += 1
                    cif = CifParser(ads_path)
                    try:
                        system_ads = cif.get_structures()[0]
                    except:
                        print('Cannot read ads cif file')
                        continue
                    file = folder + '/' + struc + '/' + ads_file
                    shutil.copyfile(ads_path, file)
                    dynamics = []
                    for atom in range(0, system_ads.num_sites):
                        dynamics.append([True, True, True])
                    i = -1
                    for e, c in zip(system_ads.species,
                                    system_ads.frac_coords):
                        i += 1
                        for es, cs in zip(system.species,
                                          system.frac_coords):
                            dist = system.lattice.get_all_distances(c, cs)
                            if es == e and dist < 0.5:
                                dynamics[i] = ([False, False, False])
                    pos = vasp.Poscar(system_ads, selective_dynamics=dynamics)
                    pos.write_file(folder+'/'+struc+'/POSCAR_'+ads_file)

                    site_folder = vaspfolder + '/' + struc + '/site_'
                    site_folder = site_folder + str(om_count[ads_type])

                    ads_vaspfolder = site_folder+'/'+ads_type
                    adsframe_vaspfolder = site_folder+'/framework_'+ads_type
                    make_folder(site_folder)
                    make_folder(ads_vaspfolder)
                    make_folder(adsframe_vaspfolder)

                    ads_potcar = 'vasp_input_files/potcar_files/POTCAR_'+\
                                 ads_type
                    shutil.copyfile(ads_potcar, ads_vaspfolder+'/POTCAR')

                    pos.write_file(adsframe_vaspfolder+'/POSCAR')
                    potcar_file = adsframe_vaspfolder+'/POTCAR'
                    make_potcar(pos, potcar_file)

                    copy_vasp_files('vasp_input_files',
                                    adsframe_vaspfolder, 'INCAR_opt',
                                    struc+'_framework_'+ads_type)
                    copy_vasp_files('vasp_input_files',
                                    ads_vaspfolder, 'INCAR_fixed',
                                    struc+'_'+ads_type)

            print(cif_name, '-', len(oms_ids), '-', len(set(oms_ids)),
                  '-', oms_ids, '-', [oms_ids.count(i) for i in set(oms_ids)])
            #print(om_count[adss[0]], om_count[adss[1]])
    print(count, ' MOFs with ', element, ' open metal sites were found')
    with open('summary.out', 'w') as summary:
        print(count, ' MOFs with ', element, ' open metal sites were found',
              file=summary)

def replace_in_file(str_match, str_replace, filename):
    lines = []
    with open(filename, 'r') as file_:
        for l in file_:
            new_l = re.sub(str_match, str_replace, l.rstrip())
            lines.append(new_l)
            # print(l.rstrip())
            # print(new_l)
            # input()
    with open(filename, 'w') as file_:
        for l in lines:
            print(l, file=file_)

def copy_vasp_files(vaspfiles, targetfolder, incarfile, struc):
    shutil.copyfile(vaspfiles+'/KPOINTS', targetfolder+'/KPOINTS')
    shutil.copyfile(vaspfiles+'/'+incarfile, targetfolder+'/INCAR')
    shutil.copyfile(vaspfiles+'/INCAR_fixed', targetfolder+'/INCAR')
    shutil.copyfile(vaspfiles+'/run_vasp', targetfolder+'/run_vasp')
    shutil.copyfile(vaspfiles+'/pbs.job_itasca', targetfolder+'/pbs.job_itasca')
    replace_in_file('name', struc, targetfolder+'/pbs.job_itasca')
    shutil.copyfile(vaspfiles+'/pbs.job_mesabi', targetfolder+'/pbs.job_mesabi')
    replace_in_file('name', struc, targetfolder+'/pbs.job_mesabi')

def make_potcar(pos, potcar_file):
    elements = []
    for s in pos.as_dict()['structure']['sites']:
        ele = s['species'][0]['element']
        if ele not in elements:
            elements.append(ele)

    with open(potcar_file, 'w') as outfile:
        for ele in elements:
            potcar_ele = 'vasp_input_files/potcar_files/POTCAR_'+ele
            with open(potcar_ele) as infile:
                for line in infile:
                    outfile.write(line)

def make_folder(folder):
    if not os.path.exists(folder):
        os.makedirs(folder)

if __name__ == '__main__':
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
        print>>parameters,struc,system[0].lattice.abc[0],
        system[0].lattice.abc[1],system[0].lattice.abc[2],
        system[0].lattice.angles[0],system[0].lattice.angles[1],
        system[0].lattice.angles[2]
        print struc,system[0].lattice.abc[0],system[0].lattice.abc[1],
        system[0].lattice.abc[2],system[0].lattice.angles[0],
        system[0].lattice.angles[1],system[0].lattice.angles[2]
        #raw_input()
        #print pos
        pos.write_file(folder+'/POSCAR_'+struc)
        sio.write_mol(system[0],folder+'/'+struc+'.xyz')
'''
