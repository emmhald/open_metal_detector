"""
Created on 12/8/2014
@author: %(Emmanuel Haldoupis)s
"""
import sys,os
import shlex
import argparse

target_folder='CORE-MOF-DB-June2014/'
filetype='cif'
num_of_batches=16

parser = argparse.ArgumentParser(description='Split file into batches')
parser.add_argument('-p','--parameter_file', nargs='?',help='Name of parameters file')
args=parser.parse_args()
params_file='parameters_'+filetype+'.txt'
if args.parameter_file:
    params_file=args.parameter_file
print params_file
print 'Spliting files into batches'
with open(params_file,'r') as params_f:
    lines=params_f.readlines()

sum_num=0
for i,struc in enumerate(lines):
    line_elements=shlex.split(struc)
    filename=line_elements[0]
    num_of_atoms=int(line_elements[1])
    sum_num+=num_of_atoms
atoms_per_batch=sum_num/num_of_batches
print filename,sum_num,atoms_per_batch

#params_prefix=params_file.split('.')[0]
params_batch_files=[]
for i in range (1,num_of_batches+1):
    params_batch_files.append('parameters_batch'+str(i)+'.in')
    open(params_batch_files[-1],'w')

batch=0
sum_num=0
for i,struc in enumerate(lines):
    if sum_num >= atoms_per_batch*(batch):
        batch+=1
        print batch
        with open('run_omsd'+str(batch)+'.sh','w') as run_f:
            print>>run_f,'#!/bin/bash'
            print>>run_f,'cd $1'
            print>>run_f,'python open_metal_detector.py -p '+params_batch_files[batch-1],' -s summary_'+str(batch)+'.out'
        os.chmod('run_omsd'+str(batch)+'.sh',0760)

    line_elements=shlex.split(struc)
    filename=line_elements[0]
    num_of_atoms=int(line_elements[1])
    sum_num+=num_of_atoms
    with open(params_batch_files[batch-1],'a') as params_f:
        print>>params_f, filename,num_of_atoms



