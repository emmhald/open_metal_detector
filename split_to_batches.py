"""
Created on 12/8/2014
@author: %(Emmanuel Haldoupis)s
"""
import sys,os
import shlex
import argparse
from pymatgen.io.cifio import CifParser
from atomic_parameters import atoms as ap

target_folder='CORE-MOF-DB-June2014/'
filetype='cif'
num_of_batches=16

parser = argparse.ArgumentParser(description='Split file into batches')
parser.add_argument('-p','--parameter_file', nargs='?',help='Name of parameters file')
args=parser.parse_args()
params_file='parameters_'+filetype+'.txt'
if args.parameter_file:
    params_file=args.parameter_file
print(params_file)
print('Spliting files into batches')
with open(params_file,'r') as params_f:
    lines=params_f.readlines()

all_files = len(lines)
atoms = []
metals = []
for i,struc in enumerate(lines):
    line_elements=shlex.split(struc)
    filename=line_elements[0]
    ciffile = target_folder+'/'+filename
    cif = CifParser(ciffile)
    system = cif.get_structures(primitive=False)[0]
    # print(len(system.sites))
    count_metals = 0
    for s in system.species:
        s = str(s)
        if ap.check_if_metal(s):
            count_metals += 1
    metals.append(count_metals)
    atoms.append(len(system.sites))
    print(float(i)/float(all_files), end="\r")
    # print(metals[-1])
    # print(atoms[-1])
# sum_num = sum(atoms)
# print(sum_num)
# input()
sum_num = 0
for n,m in zip(atoms, metals):
    sum_num += n*m

# sum_num=0
# for i,struc in enumerate(lines):
#     line_elements=shlex.split(struc)
#     filename=line_elements[0]
#     num_of_atoms=int(line_elements[1])
#     sum_num+=num_of_atoms

atoms_per_batch = sum_num/num_of_batches
print(filename,sum_num,atoms_per_batch)

#params_prefix=params_file.split('.')[0]
params_batch_files = []
for i in range (1,num_of_batches+1):
    params_batch_files.append('parameters_batch'+str(i)+'.in')
    open(params_batch_files[-1],'w')

batch=0
sum_num=0
for i,struc in enumerate(lines):
    if sum_num >= atoms_per_batch*(batch):
        batch += 1
        print(batch)
        # with open('run_omsd'+str(batch)+'.sh','w') as run_f:
        #     print('#!/bin/bash', file=run_f)
        #     print('cd $1', file=run_f)
        #     print('python open_metal_detector.py -p '+params_batch_files[batch-1],' -s summary_'+str(batch)+'.out', file=run_f)
        # os.chmod('run_omsd'+str(batch)+'.sh', 760)

    line_elements=shlex.split(struc)
    filename = line_elements[0]
    num_of_atoms=int(line_elements[1])
    sum_num += atoms[i] * metals[i]
    with open(params_batch_files[batch-1],'a') as params_f:
        print(filename,num_of_atoms, file=params_f)



