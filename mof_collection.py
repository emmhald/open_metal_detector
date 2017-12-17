
import os
import glob
import json
import time
import pickle
import warnings
import itertools
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from multiprocessing import Process, cpu_count
from mof import Helper
from mof import MofStructure
from atomic_parameters import atoms as ap


class MofCollection:

    def __init__(self, collection_folder=".", output_folder=".",
                 structure_list=None, max_size=None):
        """
        An object to

        :param collection_folder: Full path to the folder of the cif files.
        :param output_folder: Full path to the folder where the results will
        be stored.
        :param structure_list:
        :param max_size:
        """
        self.collection_folder = collection_folder
        self.output_folder = output_folder
        self.analysis_folder = output_folder+'/analysis'
        self.load_balance_filename = output_folder+'/load_balance_info'
        self.structure_list = structure_list
        self.mof_coll = []
        self.batches = []
        self.json_dicts = []
        self.max_size = max_size
        self.load_balance_index = {}
        self.site_info = {}
        self.site_df = None

        self.load_mofs()
        Helper.make_folder(output_folder)
        Helper.make_folder(self.analysis_folder)

    def load_mofs(self):
        if self.structure_list:
            print('Using only in structures list: {}'.format(self.structure_list))
            d = self.collection_folder
            mof_list = [d+'/'+file.strip()
                        for file in open(self.structure_list, 'r')]
        else:
            print('Using all the structures in the collection folder')
            mof_list = glob.glob(self.collection_folder + "/*.cif")

        for mof_file in mof_list:
            mof_name = os.path.splitext(os.path.basename(mof_file))[0]
            mof_info = {"mof_name": mof_name, "mof_file": mof_file}
            self.mof_coll.append(mof_info)

    def analyse_mofs(self, overwrite=False, num_batches=1, redo_balance=False):
        t0 = time.time()

        self.make_batches(num_batches, overwrite, redo_balance)

        processes = []
        for i, batch in enumerate(self.batches):
            p = Process(target=self.run_batch, args=(i, batch, overwrite,))
            p.start()
            processes.append(p)

        for p in processes:
            p.join()

        t1 = time.time()
        print('\nAnalysis Finished. Time required:', t1-t0)

    def run_batch(self, b, batch, overwrite):
        lb = len(batch)/100.0
        for i, mi in enumerate(batch):
            print('Batch {} {:.2f}% : '.format(b+1, i/lb), end='')
            self.analyse(mi, overwrite)

        print('FINISHED Batch {}.'.format(b + 1))

    def analyse(self, mi, overwrite):
        mof_folder = '/'.join([self.output_folder, mi['mof_name']])
        results_exist = self.check_if_results_exist(mi['mof_name'])
        if not overwrite and results_exist:
            print("Skipping {}. Results already exist and overwrite is set "
                  "to False.".format(mi['mof_name']))
            return
        print("Analysing {}".format(mi['mof_name']))
        mof = MofStructure.from_file(mi['mof_file'], primitive=False)
        mof.analyze_metals(output_folder=mof_folder, output_level='debug')

    def check_if_results_exist(self, mof_name):
        mof_folder = '/'.join([self.output_folder, mof_name])
        if os.path.isfile(mof_folder+'/'+mof_name+'.json'):
            if not os.path.isfile(mof_folder + '/' + 'analysis_running'):
                return True
        return False

    def make_batches(self, num_batches=3, overwrite=False, redo_balance=False):
        if cpu_count() < num_batches:
            warnings.warn('You requested {} batches but there are only'
                          '{} CPUs available.'.format(num_batches, cpu_count()))

        self.validate_load_balance_info(redo_balance)
        # Remove any structures not in load balancing index.
        print('\nSkipping any structures with no load balancing info.')
        subset = [mc for mc in self.mof_coll if mc['mof_name'] in
                  self.load_balance_index]
        # Remove any structures already completed
        if not overwrite:
            print('Overwrite is set to False.')
            print('Skipping any MOFs for which results already exist...',
                  end='')
            subset = [mc for mc in subset if not
                      self.check_if_results_exist(mc['mof_name'])]
            print('Done')

        self.batches = [[] for b in range(0, num_batches)]
        # Sort mof list using the load balancing index
        subset.sort(key=lambda x: self.load_balance_index[x['mof_name']])

        sum_load_balance = sum(self.load_balance_index[mi["mof_name"]]
                               for mi in subset)
        lb_per_batch = sum_load_balance / num_batches

        # Select only up to max_size to work with
        if self.max_size and len(subset) > self.max_size:
            subset = subset[0:self.max_size]
        sum_lb = 0.0
        batch = 0
        for mi in subset:
            lb = self.load_balance_index[mi["mof_name"]]
            sum_lb += lb
            self.batches[batch].append(mi)
            batch = int(sum_lb / lb_per_batch)

        for i, batch in enumerate(self.batches):
            print("Batch {0} has {1} MOFs".format(i+1, len(batch)))
        print("\n")

    def validate_load_balance_info(self, redo_balance=False):
        if os.path.isfile(self.load_balance_filename) and not redo_balance:
            with open(self.load_balance_filename, 'rb') as load_balance_file:
                self.load_balance_index = pickle.load(load_balance_file)
        else:
            if redo_balance:
                print('\nRecomputing Load Balancing Info.')
            else:
                print('\nLoading Load Balancing Info not found.')
            print('This might a while but needs to be done once.')
            print('Computing now...')
            self.compute_load_balance_index()
            print('Done\n')

        print('Validating Load Balancing Info...')
        for mi in self.mof_coll:
            if mi['mof_name'] not in self.load_balance_index:
                print("Load Balance Info is missing computing "
                      "now for {:40}".format(mi['mof_name']), end="\r")
                self.compute_load_balance_info(mi)

        with open(self.load_balance_filename, 'wb') as load_balance_file:
            pickle.dump(self.load_balance_index, load_balance_file)
        print('Done\n')

    def compute_load_balance_index(self):
        to_remove = []
        mof_files_n_100 = len(self.mof_coll) / 100
        for i, mi in enumerate(self.mof_coll):
            if not self.compute_load_balance_info(mi):
                to_remove.append(i)
            print("{:5.2f}%".format(i/mof_files_n_100), end="\r")
        with open(self.output_folder+'/could_not_load.out', 'w') as fout:
            for i in to_remove:
                print(self.mof_coll[i]['mof_file'], file=fout)
        self.mof_coll = [mc for mc in self.mof_coll if mc not in to_remove]
        print("\n")

    def compute_load_balance_info(self, mi):
        try:
            mof = MofStructure.from_file(mi['mof_file'], primitive=False)
        except Exception as e:
            print('\nAn Exception occured: {}'.format(e))
            print('Cannot load {}\n'.format(mi['mof_name']))
            return False
        self.load_balance_index[mi['mof_name']] = len(mof)*len(mof)
        return True

    def check_analysis_status(self):
        done = 0
        not_done = []
        for mi in self.mof_coll:
            mf = '/'.join([self.output_folder, mi['mof_name']])
            if self.check_if_results_exist(mf):
                done += 1
            else:
                not_done.append(mf)
        if done > 0:
            print('\nAnalysis for {} out of {} structures have been completed.'
                  .format(done, len(self.mof_coll)))
        else:
            print('\nAnalysis for no structures has been completed.'
                  .format(done, len(self.mof_coll)))
        if done != len(self.mof_coll):
            print('The following structures are missing.')
            for nd in not_done:
                print(nd)

    def load_results(self):
        """Loop over all json files in output_folder and return
        a list of dictionaries for each json file"""
        print('Reading structures...', end=' ')
        self.json_dicts = []
        for mi in self.mof_coll:
            mof_name = mi["mof_name"]
            json_file = "{0}/{1}/{1}.json".format(self.output_folder, mof_name)
            # json_file = self.output_folder+'/'+mof_name+'/'+mof_name + '.json'
            if os.path.isfile(json_file):
                json_dict = json.load(open(json_file))
            if isinstance(json_dict, dict):
                json_dict['source_name'] = self.output_folder+'/'+mof_name
                self.json_dicts.append(json_dict)
        print('Done.')

    def summarize_tfactors(self):
        """Read all the tfactors for all the open metal sites found in all the
        mofs in the json files (json_dicts) write them into files for each type
        yes_t4, no_t4, yes_t5 etc. where yes means the site has been found open.
        In addition, make histograms using this data."""
        tfac_analysis_folder = self.analysis_folder+'/tfac_analysis'
        Helper.make_folder(self.analysis_folder)
        Helper.make_folder(tfac_analysis_folder)

        self.make_df()
        sites_u = self.site_df[self.site_df['unique']]

        for n in range(4, 7):
            self.write_t_factors(sites_u, n, tfac_analysis_folder)

    def write_t_factors(self, sites, n, target):
        s_n = sites.loc[sites['number_of_linkers'] == n].copy()
        s_n['is_open_yn'] = np.where(s_n['is_open'], 'yes', 'no')
        s_n = s_n[['mof_name', 'is_open_yn', 't_factor']]
        for flag in ['yes', 'no']:
            outpath = "{}/{}_{}.out".format(target, flag, str(n))
            s = s_n[s_n['is_open_yn'] == flag]
            s.to_csv(outpath, index=False)
            fout = "{}/{}_{}_hist.out".format(target, flag, n)
            self.write_histogram(s['t_factor'], True, fout)
            fout = "{}/{}_{}_hist_abs.out".format(target, flag, n)
            self.write_histogram(s['t_factor'], False, fout)

        fig = plt.figure(figsize=(10, 5))
        plt.title('t-{} factor'.format(n))
        s_yes = s_n[s_n['is_open_yn'] == 'yes']
        s_yes['t_factor'].hist(bins=50, range=(0, 1), normed=False)
        s_no = s_n[s_n['is_open_yn'] == 'no']
        s_no['t_factor'].hist(bins=50, range=(0, 1), normed=False)
        plt.show()

    def write_histogram(self, sites, dens, target):
        hist, edges = np.histogram(sites, bins=50, range=(0, 1), density=dens)
        with open(target, 'w') as hist_file:
            w = (edges[1] - edges[0]) / 2
            for e, h in zip(edges, hist):
                print(e + w, h, file=hist_file)

    def collect_statistics(self, max_atomic_number=54):
        self.site_info = {}
        for json_dict in self.json_dicts:
            metal_sites = json_dict['metal_sites']
            mof_name = json_dict['material_name']
            for i, ms in enumerate(metal_sites):
                key = mof_name + '_' + str(i)
                self.site_info[key] = ms
                if 'all_dihedrals' in ms:
                    del self.site_info[key]['all_dihedrals']
                if 'min_dihedral' in ms:
                    del self.site_info[key]['min_dihedral']
                self.site_info[key]['mof_name'] = mof_name
        self.make_df()
        site_df_u = self.site_df.loc[self.site_df['unique']]
        site_df_o = site_df_u.loc[site_df_u['is_open']]

        all_sites = self.group_and_summarize(site_df_u, ['All-MOFs',
                                                         'All-Sites'])
        open_sites = self.group_and_summarize(site_df_o, ['All-Open-MOFs',
                                                          'All-Open-Sites'])
        s_df = pd.concat([all_sites, open_sites], axis=1)
        s_df.fillna(0.0, inplace=True)
        s_df['Per.-Open-MOF'] = 100.0 * s_df['All-Open-MOFs']/s_df['All-MOFs']
        s_df['Per.-Open-Site'] = 100.0 * s_df['All-Open-Sites'] / s_df['All-Sites']
        cols = ['All-MOFs', 'All-Open-MOFs', 'All-Sites', 'All-Open-Sites',
                'Per.-Open-MOF', 'Per.-Open-Site']
        s_df = s_df[cols]

        s_df['Per.-Open-MOF'] = s_df['Per.-Open-MOF'].apply('{:.2f} %'.format)
        s_df['Per.-Open-Site'] = s_df['Per.-Open-Site'].apply('{:.2f} %'.format)
        s_df.sort_values("All-MOFs", inplace=True, ascending=False)
        s_df.to_csv(self.analysis_folder+'/stats.out', sep=' ')

        #Collect statistics for only up to max_atomic_number
        subset = pd.Series(s_df.index).apply(ap.check_atomic_number,
                                             args=(max_atomic_number,))
        site_df_subset = s_df.loc[subset.values]
        fname = "{0}/stats_less_{1}.out".format(self.analysis_folder,
                                              max_atomic_number)
        site_df_subset.to_csv(fname, sep=' ')
        print(site_df_subset)

    def group_and_summarize(self, df, names=None):
        rename = {"mof_name": names[0], "is_open": names[1]}
        agg_dict = {"mof_name": pd.Series.nunique, "is_open": "count"}
        return df.groupby('metal').agg(agg_dict).rename(columns=rename)

    def make_df(self):
        if self.site_df is not None:
            return
        self.site_info = {}
        for json_dict in self.json_dicts:
            metal_sites = json_dict['metal_sites']
            mof_name = json_dict['material_name']
            for i, ms in enumerate(metal_sites):
                key = mof_name + '_' + str(i)
                self.site_info[key] = ms
                if 'all_dihedrals' in ms:
                    del self.site_info[key]['all_dihedrals']
                if 'min_dihedral' in ms:
                    del self.site_info[key]['min_dihedral']
                self.site_info[key]['mof_name'] = mof_name
        self.site_df = pd.DataFrame.from_dict(self.site_info, orient='index')
