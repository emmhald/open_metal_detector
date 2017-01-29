# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %(Emmanuel Haldoupis)s
"""

import sys
import os
from pymatgen import Lattice,Structure
from pymatgen.io.cif import CifParser
# import pymatgen.io.smart as sio
import pymatgen.io.vasp.inputs as vasp
import shutil
from atomic_parameters import atoms
#import matplotlib.pyplot as plt
import numpy as np
import json
import re
# from pymatgen.io.smart import read_structure, write_structure
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
    parser.add_argument('-v', '--prepare_vasp_files',
                        nargs='?',
                        help='Prepare vasp input files for MOFs that contain '
                             'the target_element specified with the -e flag',
                        default=False)

    args = parser.parse_args()
    output_folder = args.output_folder
    analysis_folder = args.analysis_folder
    parameters_file = args.parameter_file
    ele = args.target_element
    prep_vasp = args.prepare_vasp_files
    json_dicts = load_structures(output_folder, parameters_file)
    analyze_for_tfac_using_json(json_dicts, analysis_folder)
    collect_statistics(json_dicts, analysis_folder)
    if ele:
        collect_files(json_dicts, ele, analysis_folder,
                      output_folder, prep_vasp)


def read_json(filename):
    """Load a json file to a dictionary and return this dictionary"""
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
    """Loop over all json files in output_folder and return
    a list of dictionaries for each json file"""
    print('Reading structures...', end=' ')
    json_dicts = []
    with open(parameter_file, 'r') as parameters:
        for l in parameters:
            filetype = l.split('.')[-1]
            struc = l.split('.'+filetype)[0]
            json_file = output_folder+'/'+struc+'/'+struc
            json_dict = read_json(json_file)
            if isinstance(json_dict, dict):
                json_dict['source_name'] = output_folder+'/'+struc
                json_dicts.append(json_dict)
    print('Done.')
    return json_dicts


def analyze_for_tfac_using_json(json_dicts, analysis_folder):
    """Read all the tfactors for all the open metal sites found in all the
    mofs in the json files (json_dicts) write them into files for each type
    yes_t4, no_t4, yes_t5 etc. where yes means the site has been found open.
    In addition, make histograms using this data."""
    tfac_analysis_folder = analysis_folder+'/tfac_analysis'
    make_folder(analysis_folder)
    make_folder(tfac_analysis_folder)

    # r = re.compile("([a-zA-Z]+)(-?(?:\d+())?(?:\.\d*())?(?:e-?\d+())?(?:\©|\1\3))")

    import itertools
    for yn, nl in itertools.product(['yes', 'no'], [4, 5, 6]):
        outfilename = yn + '_' + str(nl) + '.out'
        outpath = tfac_analysis_folder + '/' + outfilename
        open(outpath, 'w').close()
    yes_or_no = {True: 'yes', False: 'no'}

    for json_dict in json_dicts:
        struc = json_dict['material_name']
        metal_sites = json_dict['metal_sites']
        for ms in metal_sites:
            nl = ms["number_of_linkers"]
            if 7 > int(nl) > 3 and ms['unique']:  #
                outfilename = yes_or_no[ms["is_open"]]+'_'+str(nl) + '.out'
                outpath = tfac_analysis_folder + '/' + outfilename
                with open(outpath, 'a') as outfile:
                    print(struc, yes_or_no[ms["is_open"]], ms["t_factor"],
                          file=outfile)

    for yn in ['yes', 'no']:
        for i in range(4, 7):
            filename = yn+'_'+str(i)
            outpath = tfac_analysis_folder + '/' + filename + '.out'
            with open(outpath, 'r') as datafile:
                tfac_list = []
                for l in datafile:
                    tfac = l.split(' ')[2].rstrip('\n')
                    tfac_list.append(float(tfac))

            if len(tfac_list) > 0:
                hist, edges = np.histogram(tfac_list, bins=50, range=(0, 1),
                                           density=True)
                fhist = tfac_analysis_folder+'/'+filename+'_hist.out'
                with open(fhist, 'w') as hist_file:
                    w = (edges[1]-edges[0])/2
                    for e, h in zip(edges, hist):
                        print(e+w, h, file=hist_file)

                hist, edges = np.histogram(tfac_list, bins=50, range=(0, 1),
                                           density=False)

                fhist = tfac_analysis_folder+'/'+filename+'_hist_abs.out'
                with open(fhist, 'w') as hist_file:
                    w = (edges[1]-edges[0])/2
                    for e, h in zip(edges, hist):
                        print(e+w, h, file=hist_file)


def collect_statistics(json_dicts, analysis_folder):
    if len(json_dicts) == 0:
        return
    count_open = 0.0
    stats = dict()
    mof_open = dict()
    # check_count = dict()
    for json_dict in json_dicts:
        metal_sites = json_dict['metal_sites']
        for ms in metal_sites:
            metal = ms['metal']
            stats[metal] = dict()
            stats[metal].setdefault('count_mofs', 0)
            stats[metal].setdefault('count_sites', 0)
            stats[metal].setdefault('count_open_sites', 0)
            stats[metal].setdefault('count_open_mofs', 0)
            for t in range(0, 7):
                stats[metal].setdefault('l'+str(t), 0)
            stats[metal].setdefault('lover', 0)

    for json_dict in json_dicts:
        struc = json_dict['material_name']
        if json_dict['metal_sites_found']:
            count_open += 1
        metal_sites = json_dict['metal_sites']

        all_metals = [ms['metal'] for ms in metal_sites if ms['unique']]

        # oms_type = [(ms['metal'], ms['number_of_linkers']) for ms in metal_sites
        #             if ms['unique'] and ms['is_open']]

        oms_metals = [ms['metal'] for ms in metal_sites
                      if ms['unique'] and ms['is_open']]

        oms_num_linkers = [ms['number_of_linkers'] for ms in metal_sites
                           if ms['unique'] and ms['is_open']]

        for m in set(all_metals):
            stats[m]['count_mofs'] += 1

        for m in set(oms_metals):
            stats[m]['count_open_mofs'] += 1

        for m in all_metals:
            stats[m]['count_sites'] += 1

        for l, m in zip(oms_num_linkers, oms_metals):
            # m, l = om
            nl = 'l'+str(l)
            if l > 6:
                nl = 'lover'
            stats[m]['count_open_sites'] += 1
            stats[m][nl] += 1


    # This returns a sorted tuple based on keyfunc,
    # which uses the key count_open to reverse sort the stats dictionary
    # A more consise but less clear sollution would be
    # stats_sorted = sorted(stats.items(), key = lambda tup :
    # (-tup[1]["count_open"]))
    stats_sorted = sorted(stats.items(), key=keyfunc)
    printouts = []
    printouts_less_d4 = []
    for stat in stats_sorted:
        # Since stats_sorted is a sorted tuple of the dictionary stats,
        # the first element corresponds to the keys from the stats dictionary
        # and the second element to the value, in this dictionary holding
        # the stats for each metal
        metal = stat[0]
        s = stat[1]

        percent = 100*float(s['count_open_mofs'])/float(s['count_mofs'])
        percent_mof = "{0:.2f} %".format(percent)
        percent = 100*float(s['count_open_sites'])/float(s['count_sites'])
        percent_oms = "{0:.2f} %".format(percent)
        printout = [metal, s['count_mofs'], s['count_open_mofs'],
                    s['count_sites'], s['count_open_sites'],
                    percent_mof, percent_oms,
                    s['l0'], s['l1'], s['l2'], s['l3'], s['l4'], s['l5'],
                    s['l6'], s['lover']]
        printouts.append(printout)
        # if ap.is_group1(metal):
        #     continue
        if ap.is_d4_or_less(metal):
            printouts_less_d4.append(printout)


    tot_mofs = len(json_dicts)
    precent_open_mofs = 100.0 * count_open / tot_mofs
    tot_sites = sum([s[1]['count_sites'] for s in stats_sorted])
    tot_sites_open = sum([s[1]['count_open_sites'] for s in stats_sorted])

    try:
        precent_open_sites = 100.0 * tot_sites_open / tot_sites
    except ZeroDivisionError:
        "No metal sites were found"
        precent_open_sites = 0.0


    print("\nTotal MOFs: {0:}\nMOFs with Open Metal Site: {1:}\n"
          "Percentage of MOFs with OMS: {2:2.2f}%\n"
          "".format(tot_mofs, count_open, precent_open_mofs))
    print("Total Metal Sites: {0:}\nOpen Metal Sites: {1:}\n"
          "Percentage of Metal Sites which are open: {2:2.2f}% \n"
          "".format(tot_sites, tot_sites_open, precent_open_sites))

    titles = ['Metal', 'All-MOFs', 'All-Open-MOFs',
              'All-Sites', 'Open-Metal-Sites', 'Per.-Open-MOF',
              'Per.-Open-Site', 'l0', 'l1', 'l2', 'l3', 'l4', 'l5', 'l6',
              'l-over-6']

    print_stats(titles, printouts, sys.stdout)

    with open(analysis_folder+'/stats.out', 'w') as fstats:
        print_stats(titles, printouts, fstats)

    with open(analysis_folder+'/stats_less_d4.out', 'w') as fstats:
        print_stats(titles, printouts_less_d4, fstats)


def print_stats(titles, printouts, fstats):

    print("{0:6}{1:^12}{2:^12}{3:^12}{4:^18}{5:^14}{6:^14}{7:^6}{8:^6}{9:^6}"
          "{10:^6}{11:^6}{12:^6}{13:^8}{14:^8}".format(*titles), file=fstats)
    for p in printouts:
        print("{0:6}{1:^12}{2:^12}{3:^12}{4:^18}{5:^14}{6:^15}{7:^6}{8:^6}"
              "{9:^6}{10:^6}{11:^6}{12:^6}{13:^8}{14:^8}".format(*p),
              file=fstats)


def copy_folder(dest, src):
    if not os.path.exists(dest):
        os.makedirs(dest)
    s = src.split('/')[-1]
    d = os.path.join(dest, s)
    if not os.path.exists(d):
        shutil.copytree(src, d)


def keyfunc(tup):
    key, d = tup
    return -d["count_mofs"]


def fetch_list_of_n_omtype_tf(json_dict):
    num_of_ligands = []
    om_type = []
    t_factor = []
    for ms in json_dict["metal_sites"]:
        num_of_ligands.append(ms["number_of_linkers"])
        om_type.append(ms["is_open"])
        t_factor.append(ms["t_factor"])
    return num_of_ligands, om_type, t_factor


def collect_files(json_dicts, element, analysis_folder, output_folder,
                  prep_vasp):

    folder = analysis_folder+'/selected_mofs/'+element
    vaspfolder = analysis_folder+'/calcs/'
    cif_folder_out = analysis_folder+'/selected_mofs/'+element+'/cif_files'
    make_folder(folder)
    make_folder(vaspfolder)
    make_folder(cif_folder_out)

    open(folder+'/summary.out', 'w').close()
    cif_folder = output_folder+'/open_metal_mofs/'
    open('summary.out', 'w').close()

    # print('\nMOF', '#OMS', '#OMS_types', 'OMS ids', '#OMS_per_type')
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
        # if contains_metal_but_closed and not struc_contains_open_metal:
        #     print(json_dict['material_name'], 'closed')
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
            system = cif.get_structures(primitive=False)[0]
            struc_xyz_ = folder + '/' + struc + '/' + struc + '.xyz'
            system.to(fmt='xyz', filename=struc_xyz_)

            name = folder + '/' + struc + '/' + cif_name
            shutil.copyfile(cif_folder+cif_name, name)
            out_cif_name = cif_folder_out + '/' + cif_name
            shutil.copyfile(cif_folder+cif_name, out_cif_name)

            if not prep_vasp:
                continue

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
                        system_ads = cif.get_structures(primitive=False)[0]
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
                    try:
                        shutil.copyfile(ads_potcar, ads_vaspfolder+'/POTCAR')
                    except:
                        print(ads_potcar, 'File is missing')

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
            try:
                with open(potcar_ele) as infile:
                    for line in infile:
                        outfile.write(line)
            except:
                print(potcar_ele, 'File is missing.')


def make_folder(folder):
    if not os.path.exists(folder):
        os.makedirs(folder)

if __name__ == '__main__':
    main()
