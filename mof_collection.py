
import pickle
from multiprocessing import Process,cpu_count
import time
import warnings
from mof import MofStructure
from mof import Helper
import glob
import os


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
        self.load_balance_filename = output_folder+'/load_balance_info'
        self.structure_list = structure_list
        self.mof_coll = []
        self.batches = []
        self.max_size = max_size
        self.load_balance_index = {}
        self.load_mofs()

        Helper.make_folder(output_folder)

    def load_mofs(self):
        if self.structure_list:
            print('Using only in structures list: {}'.format(self.structure_list))
            d = self.collection_folder
            mof_list = [d+'/'+file.strip()
                        for file in open(self.structure_list, 'r')]
        else:
            print('Using all the structures in the collection folder')
            mof_list = glob.glob(self.collection_folder + "*.cif")

        for mof_file in mof_list:
            mof_name = mof_file.split('/')[-1].split('.')[0]

            mof_info = {"mof_name": mof_name, "mof_file": mof_file}
            self.mof_coll.append(mof_info)

    def analyse_mofs(self, overwrite=False, num_batches=1, redo_balance=False):
        t0 = time.time()

        # if not self.batches:
        self.make_batches(num_batches, redo_balance)

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
        for i, mi in enumerate(batch):
            print('Batch {} MOF {} : '.format(b+1, i+1, mi['mof_name']),
                  end='')
            self.analyse(mi, overwrite)

        print('FINISHED Batch {}.'.format(b + 1))

    def analyse(self, mi, overwrite):
        mof_folder = '/'.join([self.output_folder, mi['mof_name']])
        results_exist = self.check_if_results_exist(mof_folder)
        if not overwrite and results_exist:
            print("Skiping {}. Results already exist and overwrite is set "
                  "to False.".format(mi['mof_name']))
            return
        print("Analysing {} ...".format(mi['mof_name']))
        mof = MofStructure.from_file(mi['mof_file'], primitive=False)
        mof.analyze_metals(output_folder=mof_folder, output_level='debug')

    def check_if_results_exist(self, mof_folder):
        mof_name = mof_folder.split('/')[-1]
        if os.path.isfile(mof_folder+'/'+mof_name+'.json'):
            if not os.path.isfile(mof_folder + '/' + 'analysis_running'):
                return True

        return False

    def make_batches(self, num_batches=3, redo_balance=False):
        if cpu_count() < num_batches:
            warnings.warn('You requested {} batches but there are only'
                          '{} CPUs available.'.format(num_batches, cpu_count()))

        self.validate_load_balance_info(redo_balance)

        self.batches = [[] for b in range(0, num_batches)]
        # Sort mof list using the load balancing index
        self.mof_coll.sort(key=lambda x: self.load_balance_index[x['mof_name']])
        # Select only up to max_size to work with
        self.mof_coll = self.mof_coll[0:self.max_size]
        sum_load_balance = sum(self.load_balance_index[mi["mof_name"]]
                               for mi in self.mof_coll)
        lb_per_batch = sum_load_balance / num_batches

        sum_lb = 0.0
        batch = 0
        for mi in self.mof_coll:
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
        mof_files_n_100 = len(self.mof_coll) / 100
        for i, mi in enumerate(self.mof_coll):
            self.compute_load_balance_info(mi)
            print("{:5.2f}%".format(i/mof_files_n_100), end="\r")
        print("\n")

    def compute_load_balance_info(self, mi):
        mof = MofStructure.from_file(mi['mof_file'], primitive=False)
        self.load_balance_index[mi['mof_name']] = len(mof)  #*len(mof)

    def check_analysis_status(self):
        # mof_folders = [f.path for f in os.scandir(self.output_folder) if f.is_dir()]
        done = 0
        not_done = []
        # for mf in mof_folders:
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
