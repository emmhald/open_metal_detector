
import os
import glob
import json
import time
import pickle
import shutil
import warnings
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from multiprocessing import Process, cpu_count
from mof import Helper
from mof import MofStructure
from atomic_parameters import atoms as ap
pd.options.display.max_rows = 1000


class MofCollection:
    separator = "".join(['-'] * 50)

    def __init__(self, path_list, analysis_folder='analysis'):
        """
        An object to

        :param analysis_folder: Full path to the folder where the results will
        be stored.
        """
        self.analysis_folder = analysis_folder
        self.summary_folder = self.analysis_folder + '/summary'
        self.load_balance_filename = self.analysis_folder + '/load_balance_info'
        self.path_list = path_list
        self.mof_coll = []
        self.batches = []
        self._results = None
        self._results_df = None
        self.load_balance_index = {}
        self.analysis_limit = None

        self.load_mofs()
        Helper.make_folder(self.analysis_folder)
        Helper.make_folder(self.summary_folder)

    def __len__(self):
        return len(self.mof_coll)

    def __repr__(self):
        prnt_str = self.separator
        prnt_str += "\nThis collection holds information for " \
                    "{} MOFs.\n".format(len(self))
        if self.analysis_folder is None:
            prnt_str += "Analysis folder is not " \
                        "set.\n".format(self.analysis_folder)
        else:
            f = os.path.abspath(self.analysis_folder)
            prnt_str += "Analysis folder is: {}\n\n".format(f)
        prnt_str += "List of cif files in collection:\n\n"
        for mc in self.mof_coll:
            prnt_str += "{}\n".format(mc['mof_file'])
        prnt_str += self.separator
        return prnt_str

    @classmethod
    def from_folder(cls, collection_folder, analysis_folder='analysis',
                    name_list=None):

        if name_list:
            print(cls.separator)
            print('Using only Mofs in the name list.')
            print(cls.separator)
            d = collection_folder
            path_list = [d+'/'+name for name in name_list]
        else:
            path_list = glob.glob(collection_folder + "/*.cif")

        return cls(path_list, analysis_folder)

    def load_mofs(self):
        for mof_file in self.path_list:
            mof_name = os.path.splitext(os.path.basename(mof_file))[0]
            mof_info = {"mof_name": mof_name, "mof_file": mof_file}
            self.mof_coll.append(mof_info)

    def analyse_mofs(self, overwrite=False, num_batches=1, redo_balance=False,
                     analysis_folder=None, analysis_limit=None):
        self.analysis_limit = analysis_limit
        if analysis_folder is not None:
            self.analysis_folder = analysis_folder
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
        print('\nAnalysis Finished. Time required:{:.2f} sec'.format(t1 - t0))

    def run_batch(self, b, batch, overwrite):
        lb = (len(batch)-1)/100.0
        for i, mi in enumerate(batch):
            print('Batch {} {:.2f}% : '.format(b+1, i/lb), end='')
            self.analyse(mi, overwrite)

        print('FINISHED Batch {}.'.format(b + 1))

    def analyse(self, mi, overwrite):
        mof_folder = "{}/{}".format(self.analysis_folder, mi['mof_name'])
        results_exist = self.check_if_results_exist(mi['mof_name'])
        if not overwrite and results_exist:
            print("Skipping {}. Results already exist and overwrite is set "
                  "to False.".format(mi['mof_name']))
            return
        print("Analysing {}".format(mi['mof_name']))
        mof = self.create_mof_from_cifile(mi['mof_file'])
        # mof = MofStructure.from_file(mi['mof_file'], primitive=False)
        if mof:
            mof.analyze_metals(output_folder=mof_folder, output_level='debug')

    def check_if_results_exist(self, mof_name):
        mof_folder = "{}/{}".format(self.analysis_folder, mof_name)
        if os.path.isfile(mof_folder+'/'+mof_name+'.json'):
            if not os.path.isfile(mof_folder + '/' + 'analysis_running'):
                return True
        return False

    def make_batches(self, num_batches=1, overwrite=False, redo_balance=False):
        if cpu_count() < num_batches:
            warnings.warn('You requested {} batches but there are only {}'
                          ' CPUs available.'.format(num_batches, cpu_count()))
        b_s = {1: 'batch', 2: 'batches'}[min(num_batches, 2)]
        print('{} {} requested. '.format(num_batches, b_s))
        print('Overwrite is set to {}. '.format(overwrite))
        print(self.separator)
        self.validate_load_balance_info(redo_balance)
        # Remove any structures not in load balancing index.
        subset = [mc for mc in self.mof_coll if mc['mof_name'] in
                  self.load_balance_index]
        # If there is no balancing info for a MOF at this point it means
        # that it could not be read.
        if len(self.mof_coll) != len(subset):
            print('\nSkipping {} structures that could not be read.'
                  ' '.format(len(self.mof_coll)-len(subset)))

        # Remove any structures already completed
        if not overwrite:
            print('\nChecking if results for any of the MOFs exist...')
            all_ = len(subset)
            subset = [mc for mc in subset if not
                      self.check_if_results_exist(mc['mof_name'])]
            print('Skipping {} MOFs because results were found. '
                  ''.format(all_ - len(subset)))

        self.batches = [[] for b in range(0, num_batches)]
        # Sort mof list using the load balancing index
        subset.sort(key=lambda x: self.load_balance_index[x['mof_name']])

        sum_load_balance = sum(self.load_balance_index[mi["mof_name"]]
                               for mi in subset)
        lb_per_batch = sum_load_balance / num_batches

        # Select only up to analysis_limit to work with
        if self.analysis_limit and len(subset) > self.analysis_limit:
            subset = subset[0:self.analysis_limit]
        sum_lb = 0.0
        batch = 0
        for mi in subset:
            lb = self.load_balance_index[mi["mof_name"]]
            sum_lb += lb
            self.batches[batch].append(mi)
            batch = int(sum_lb / lb_per_batch)
        print(self.separator)
        for i, batch in enumerate(self.batches):
            print("Batch {0} has {1} MOFs".format(i+1, len(batch)))
        print(self.separator)

    def validate_load_balance_info(self, redo_balance=False):
        if os.path.isfile(self.load_balance_filename) and not redo_balance:
            with open(self.load_balance_filename, 'rb') as load_balance_file:
                self.load_balance_index = pickle.load(load_balance_file)
            print('Validating Load Balancing Info...')
            for mi in self.mof_coll:
                if mi['mof_name'] not in self.load_balance_index:
                    print("Load Balance Info is missing for {}. Computing"
                          " now".format(mi['mof_name']))
                    self.compute_load_balance_info(mi)
        else:
            if redo_balance:
                print('Recomputing Load Balancing Info.')
            else:
                print('Load Balancing Info not found.')
                print('This might take a while but needs to be done once.')
            print('\nComputing now...')
            self.compute_load_balance_index()

        with open(self.load_balance_filename, 'wb') as load_balance_file:
            pickle.dump(self.load_balance_index, load_balance_file)
        print('Done')

    def compute_load_balance_index(self):
        mof_files_n_100 = (len(self.mof_coll) - 1) / 100
        for i, mi in enumerate(self.mof_coll):
            print("{:5.2f}% {:100}".format(i / mof_files_n_100, mi['mof_name']),
                  end="\r")
            self.compute_load_balance_info(mi)
        print("")

    def compute_load_balance_info(self, mi):
        mof = self.create_mof_from_cifile(mi['mof_file'])
        if not mof:
            return False
        self.load_balance_index[mi['mof_name']] = len(mof)*len(mof)
        return True

    @staticmethod
    def create_mof_from_cifile(path_to_mof):
        try:
            mof = MofStructure.from_file(path_to_mof, primitive=False)
        except Exception as e:
            print('\nAn Exception occured: {}'.format(e))
            print('Cannot load {}\n'.format(path_to_mof))
            return False
        return mof

    def check_structures(self):
        print('\nChecking structures...')
        self.load_balance_index = {}
        self.compute_load_balance_index()
        read = [mc for mc in self.mof_coll if mc['mof_name'] in
                self.load_balance_index]
        not_read = [mc for mc in self.mof_coll if mc['mof_name'] not in
                    self.load_balance_index]
        print('Done')
        print('\nChecked {} structures.'.format(len(self.mof_coll)))
        print('{} were read.'.format(len(read)))
        print('{} were NOT read.'.format(len(not_read)))

        print('\nThe following structures could not be read:')
        for i, mi in enumerate(not_read):
            print("{}".format(mi['mof_name']))
        print('\nFinished checking structures.')

    def check_analysis_status(self):
        done = 0
        not_done = []
        for mi in self.mof_coll:
            mf = '/'.join([self.analysis_folder, mi['mof_name']])
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

    def filter_collection(self, using_filter=None,
                          new_collection_folder=None,
                          new_analysis_folder=None):
        """
        Filter a collection given a set of filters. Calling this method of a
         MofCollection applies the filter and creates a new collection for the
         Mofs that match the filter. The cif files that match the filter are
         copied to the new_collection_folder.
         For now the filter can only be a list of metal atoms.
        :param using_filter: Filter used to identify MOFs with certain
        characteristics
        :param new_collection_folder: Folder where the cif files matching the
        collection will be copied to.
        :param new_analysis_folder: Folder where the analysis for files matching
        the collection will be copied to.
        :return: A MofCollection with only the filtered MOFs pointing to the
        new_collection_folder
        """
        path_list = []
        mof_files_n_100 = (len(self.mof_coll) - 1) / 100
        for i, mi in enumerate(self.mof_coll):
            print("{:5.2f}% {:100}".format(i / mof_files_n_100, mi['mof_name']),
                  end="\r")
            mof = self.create_mof_from_cifile(mi['mof_file'])
            if not mof:
                continue
            for key in using_filter:
                value = using_filter[key]
                if all([v in mof.summary[key] for v in value]):
                    path_list.append(mi['mof_file'])

        print('\n')
        found_s = {0: "No", 1: len(path_list)}[min(1, len(path_list))]
        print('{} MOFs were matched using the provided filter.'.format(found_s))
        if len(path_list) == 0:
            print('No collection returned.')
            return None
        print('Returning a new collection using the matched MOFs.')

        sub_collection = self.make_subset(path_list)
        if new_collection_folder:
            print('The matched cif files for these MOFs will be copied to'
                  ' the specified folder:\n{}'.format(new_collection_folder))
            sub_collection.copy_cifs(new_collection_folder)
        else:
            print('The cif files will point to original locations.')
        print(self.separator)
        return sub_collection

    def make_subset(self, path_list):
        subset = MofCollection(path_list,
                               analysis_folder=self.analysis_folder)
        return subset

    def copy_cifs(self, target_folder):
        tf_abspath = os.path.abspath(target_folder)
        Helper.make_folder(tf_abspath)
        print('The cif files for this collection will be copied to'
              ' the specified folder:\n{}'.format(tf_abspath))
        print('The cif paths will be updated.')

        mofl_col_tmp = list(self.mof_coll)
        for i, mi in enumerate(mofl_col_tmp):
            destination_path = "{}/{}.cif".format(tf_abspath, mi['mof_name'])
            updated_mof_info = {"mof_name": mi['mof_name'],
                                "mof_file": destination_path}
            self.mof_coll[i] = updated_mof_info
            if not os.path.isfile(destination_path):
                shutil.copyfile(mi['mof_file'], destination_path)
        print(self.separator)

    def copy_results(self, target_folder):
        tf_abspath = os.path.abspath(target_folder)
        Helper.make_folder(tf_abspath)
        print('The result files for this collection will be copied to'
              ' the specified folder:\n{}'.format(tf_abspath))
        print('The analysis folder will be updated.')

        for i, mi in enumerate(self.mof_coll):
            mof_name = mi['mof_name']
            if self.check_if_results_exist(mof_name):
                destination_path = tf_abspath
                source_path = "{}/{}".format(self.analysis_folder, mof_name)
                Helper.copy_folder(destination_path, source_path)
        self.analysis_folder = target_folder
        print(self.separator)

    @property
    def results(self):
        """Loop over all json files in output_folder and return
        a list of dictionaries for each json file"""
        if self._results is not None:
            print('Results are already loaded.')
            return self._results

        print('Loading results from \n {}'.format(self.analysis_folder))
        self._results = []
        lb = (len(self.mof_coll) - 1)/100.0
        for i, mi in enumerate(self.mof_coll):
            mof_name = mi["mof_name"]
            json_file = "{0}/{1}/{1}.json".format(self.analysis_folder,
                                                  mof_name)
            print("{:5.2f}% {:100}".format(i / lb, mi['mof_name']), end="\r")
            json_dict = None
            if os.path.isfile(json_file):
                json_dict = json.load(open(json_file))
            if isinstance(json_dict, dict):
                json_dict['source_name'] = self.analysis_folder + '/' + mof_name
                self._results.append(json_dict)
            else:
                print('Problem reading results file for{}'.format(mof_name))
        print('Found results for {} MOFs.'.format(len(self._results)))
        return self._results

    @property
    def results_df(self):
        if self._results_df is not None:
            return self._results_df
        site_info = {}
        count_sites = 0
        count_mofs = 0
        mofs = []
        for json_dict in self.results:
            metal_sites = json_dict['metal_sites']
            mof_name = json_dict['material_name']
            mofs.append(mof_name)
            if len(metal_sites) == 0:
                print('No Metal Found in {}'.format(mof_name))
            for i, ms in enumerate(metal_sites):
                key = mof_name + '_' + str(i)
                site_info[key] = ms
                if 'all_dihedrals' in ms:
                    del site_info[key]['all_dihedrals']
                if 'min_dihedral' in ms:
                    del site_info[key]['min_dihedral']
                site_info[key]['mof_name'] = mof_name
                count_sites += 1
            count_mofs += 1
        self._results_df = pd.DataFrame.from_dict(site_info, orient='index')
        return self._results_df

    def summarize_tfactors(self):
        """Read all the tfactors for all the open metal sites found in all the
        mofs in the json files (json_dicts) write them into files for each type
        yes_t4, no_t4, yes_t5 etc. where yes means the site has been found open.
        In addition, make histograms using this data."""
        tfac_analysis_folder = self.summary_folder + '/tfac_analysis'
        Helper.make_folder(self.summary_folder)
        Helper.make_folder(tfac_analysis_folder)

        df = self.results_df.copy()
        sites_u = df[df['unique']]

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
        df = self.results_df.copy()
        site_df_u = df.loc[df['unique']]
        site_df_o = site_df_u.loc[site_df_u['is_open']]

        all_sites = self.group_and_summarize(site_df_u, ['MOFs',
                                                         'Sites'])
        open_sites = self.group_and_summarize(site_df_o, ['MOFs_with_OMS',
                                                          'OMS'])

        s_df = pd.concat([all_sites, open_sites], axis=1)
        s_df.fillna(0.0, inplace=True)
        s_df = s_df.astype(int)

        s_df['MOFs_with_OMS(%)'] = 100.0 * s_df['MOFs_with_OMS']/s_df['MOFs']
        s_df['OMS (%)'] = 100.0 * s_df['OMS'] / s_df['Sites']
        cols = ['MOFs', 'MOFs_with_OMS', 'Sites', 'OMS',
                'MOFs_with_OMS(%)', 'OMS (%)']
        s_df = s_df[cols]

        s_df['MOFs_with_OMS(%)'] = s_df['MOFs_with_OMS(%)'].apply('{:.2f} %'
                                                                  ''.format)
        s_df['OMS (%)'] = s_df['OMS (%)'].apply('{:.2f} %'.format)
        s_df.sort_values("MOFs", inplace=True, ascending=False)
        s_df.to_csv(self.summary_folder + '/stats.out', sep=' ')

        #Collect statistics for only up to max_atomic_number
        subset = pd.Series(s_df.index).apply(ap.check_atomic_number,
                                             args=(max_atomic_number,))
        site_df_subset = s_df.loc[subset.values]
        fname = "{0}/stats_less_{1}.out".format(self.summary_folder,
                                                max_atomic_number)
        site_df_subset.to_csv(fname, sep=' ')
        num_mofs = df['mof_name'].nunique()
        num_oms_mofs = df[df['is_open']]['mof_name'].nunique()
        num_sites = len(site_df_u)
        num_oms_sites = len(site_df_u[site_df_u['is_open']])
        print(self.separator)
        print('Number of total MOFs: {}'.format(num_mofs))
        print('Number of total MOFs with open metal sites: {}'
              ''.format(num_oms_mofs))
        print('Number of total unique sites: {}'.format(num_sites))
        print('Number of total unique open metal sites: {}'
              ''.format(num_oms_sites))
        print(self.separator)
        print("Summary Table\n")
        print(site_df_subset)

    def group_and_summarize(self, df, names=None):
        rename = {"mof_name": names[0], "is_open": names[1]}
        agg_dict = {"mof_name": pd.Series.nunique, "is_open": "count"}
        return df.groupby('metal').agg(agg_dict).rename(columns=rename)
