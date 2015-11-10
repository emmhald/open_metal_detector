# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 19:42:43 2013

@author: emmhald
"""
from __future__ import print_function
import shlex
import sys

class xyz_file:
    """A class to read and write an xyz file"""
    filename_in=''
    filename_out=''
    filename_out_fractional=''
    comment_out=''
    coords=[]
    coords_frac=[]
    elements=[]
    charge=[]
    labels=[]
    uc_params=[]
    number_of_atoms=0

    def load_parameters(self,params):
        self.uc_params=params

    def open_xyz(self, filename_in, *args):
        self.filename_in = filename_in
        self.clear()
        xyz_file_obj= open(self.filename_in,'r')
        line = xyz_file_obj.readline().strip()
        self.number_of_atoms=int(float(line)) #(line_counter+1)-1
        line = xyz_file_obj.readline()
        line_elements=shlex.split(line)
        if len(args) > 0 and args[0]==1:
            if len(line_elements) ==7:
                self.uc_params.append(float(line_elements[1]))
                self.uc_params.append(float(line_elements[2]))
                self.uc_params.append(float(line_elements[3]))
                self.uc_params.append(float(line_elements[4]))
                self.uc_params.append(float(line_elements[5]))
                self.uc_params.append(float(line_elements[6]))
            else:
                print("\nDoes the second line of the xyz file contain the uc parameters?")
                print('The following format should be used:')
                print("Name a b c alpha beta gamma \n")
                sys.exit('Exiting')

        line_counter=0
        for line in xyz_file_obj:
            if not line.strip():
                break
            line_counter+=1
            #if line_counter>1:
            line_elements=shlex.split(line)
            self.elements.append(line_elements[0])
            for j in range(1,4):
                self.coords.append(float(line_elements[j]))
            if (len(line_elements) >5):
                self.charge.append(float(line_elements[4]))
            if (len(line_elements) >4):
                self.labels.append(line_elements[-1])
        if line_counter != self.number_of_atoms:
            sys.exit( 'The number of atoms does not match the number of lines in the xyz file. Stopping')
        self.check_charges()
        xyz_file_obj.close()
        #self.number_of_atoms=(line_counter+1)-1

    def get_atoms(self):
        atoms = []
        for atom in range(len(self.coords)//3):
            x = atom*3
            y = atom*3+1
            z = atom*3+2
            atoms.append([self.elements[atom], self.coords[x], self.coords[y],
                         self.coords[z]])
        return atoms


    def check_charges(self) :
        if len(self.charge) >0 :
            sum_of_charges=0
            sum_of_charges=sum(self.charge)
            print('*The total charge of the cluster is:',sum_of_charges)
            print(' ')

    def write(self):
        xyz_file_obj_out= open(self.filename_out,'w')
        print(len(self.coords)//3,end = '\n\n', file = xyz_file_obj_out)
        for atom in range(int(len(self.coords)/3)):
            x=atom*3;y=atom*3+1;z=atom*3+2
            if (len(self.labels) ==  len(self.coords)/3 and len(self.charge) != len(self.coords)/3):
                print('{:<3s} {:>9.5f} {:>9.5f} {:>9.5f} {:>s}'.format(self.elements[atom],self.coords[x],self.coords[y],self.coords[z],self.labels[atom]),file = xyz_file_obj_out) #,self.charge[atom]
            elif (len(self.charge) ==  len(self.coords)/3 and len(self.labels) != len(self.coords)/3):
                print('{:<3s} {:>9.5f} {:>9.5f} {:>9.5f} {:>9.5f}'.format(self.elements[atom],self.coords[x],self.coords[y],self.coords[z],self.charge[atom]),file = xyz_file_obj_out)
            elif (len(self.charge) ==  len(self.coords)/3 and len(self.labels) == len(self.coords)/3):
                print('{:<3s} {:>9.5f} {:>9.5f} {:>9.5f} {:>9.5f} {:>s} '.format(self.elements[atom],self.coords[x],self.coords[y],self.coords[z],self.charge[atom],self.labels[atom]),file = xyz_file_obj_out)
            else:
                print('{:<3s} {:>9.5f} {:>9.5f} {:>9.5f}'.format(self.elements[atom],self.coords[x],self.coords[y],self.coords[z]),file = xyz_file_obj_out)
        xyz_file_obj_out.close()


    def write_fractional(self):
        self.filename_out_fractional
        xyz_file_obj_out= open(self.filename_out_fractional,'w')
        print(len(self.coords_frac)//3,'\n\n',file = xyz_file_obj_out)
        for atom in range(len(self.coords_frac)/3):
            x=atom*3;y=atom*3+1;z=atom*3+2
            if (len(self.labels) ==  len(self.coords_frac)/3 and len(self.charge) != len(self.coords_frac)/3):
                print('{:<3s} {:>9.5f} {:>9.5f} {:>9.5f} {:>s}'.format(self.elements[atom],self.coords_frac[x],self.coords_frac[y],self.coords_frac[z],self.labels[atom]),file = xyz_file_obj_out) #,self.charge[atom]
            elif (len(self.charge) ==  len(self.coords)/3 and len(self.labels) != len(self.coords)/3):
                print('{:<3s} {:>9.5f} {:>9.5f} {:>9.5f} {:>9.5f}'.format(self.elements[atom],self.coords_frac[x],self.coords_frac[y],self.coords_frac[z],self.charge[atom]),file = xyz_file_obj_out)
            elif (len(self.charge) ==  len(self.coords)/3 and len(self.labels) == len(self.coords)/3):
                print('{:<3s} {:>9.5f} {:>9.5f} {:>9.5f}  {:>9.5f} {:>s}'.format(self.elements[atom],self.coords_frac[x],self.coords_frac[y],self.coords_frac[z],self.charge[atom],self.labels[atom]),file = xyz_file_obj_out)
            else:
                print>>xyz_file_obj_out, '{:<3s} {:>9.5f} {:>9.5f} {:>9.5f}'.format(self.elements[atom],self.coords_frac[x],self.coords_frac[y],self.coords_frac[z])
        xyz_file_obj_out.close()

    def write_fort77(self):
        filename=self.filename_out.split('.xyz')[0]
        fort77_file_obj_out= open(filename+'_fort.77','w')
        #print>>fort77_file_obj_out,len(self.coords)/3,'\n'
        box=str(1)
        for atom in range(len(self.coords)//3):
            x=atom*3;y=atom*3+1;z=atom*3+2
 #Nchain  Box  Label Type   X       Y         Z      Charge
            if (len(self.charge) ==  len(self.coords)/3 and len(self.labels) == len(self.coords)/3):
                Type=self.labels[atom].split('_')[0]
                Label=self.labels[atom].split('_')[1]
                Nchain=self.labels[atom].split('_')[2]
                print('{:<3s} {:<3s} {:<3s} {:<3s} {:>9.5f} {:>9.5f} {:>9.5f} {:>9.5f}'.format(Nchain,box,Label,self.elements[atom],self.coords[x],self.coords[y],self.coords[z],self.charge[atom]), file=fort77_file_obj_out)
            else:
                print(len(self.charge), len(self.coords)/3, len(self.labels), len(self.coords)/3)
                exit('Cannot print fort.77')

        fort77_file_obj_out.close()


    #def make_fractional_coordinates(self):


    def xyz_from_ref(self,refrence):
        self.number_of_atoms=refrence.number_of_ref_points*refrence.number_of_atoms_ads
        self.coords=[]
        self.elements=[]
        for j in range(refrence.number_of_ref_points):
            for k in range(refrence.number_of_atoms_ads):
              self.coords.append(refrence.coords[str(k)][j*3])
              self.coords.append(refrence.coords[str(k)][j*3+1])
              self.coords.append(refrence.coords[str(k)][j*3+2])
              self.elements.append(refrence.ads_label[k].split('_')[0])

    def xyz_append(self,xyz_name):

        xyz_obj=xyz_file()
        xyz_obj.filename_in=xyz_name
        xyz_obj.open_file()

        self.number_of_atoms=self.number_of_atoms+xyz_obj.number_of_atoms
        for coord in xyz_obj.coords:
            self.coords.append(coord)
        for element in xyz_obj.elements:
            self.elements.append(element)
    def make_xyz(self,xyz_name_out,coords,elements):
        self.coords=[]
        self.elements=[]
        for i in range(len(coords)):
            self.coords.append(coords[i][0])
            self.coords.append(coords[i][1])
            self.coords.append(coords[i][2])
        self.elements=elements
        self.filename_out=xyz_name_out
        self.write()
        #print 'an xyz file with the name: ',xyz_name_out,' was created.'
    def return_mp_structure_lists(self):
        coords=[]
        elements=[]
        for atom in range(len(self.coords)/3):
            x=atom*3;y=atom*3+1;z=atom*3+2
            coord=[self.coords[x],self.coords[y],self.coords[z]]
            coords.append(coord)
            elements.append(self.elements[atom])
        return elements,coords


    def clear(self):
        #self.filename_in=''
        #self.filename_out=''
        #self.comment_out=''
        self.coords=[]
        self.elements=[]
        self.charge=[]
        self.labels=[]
        self.number_of_atoms=0


class cssr_file:
    """A class to read and write an cssr file for use with th MCCCS code"""
    filename_in=''
    filename_out=''
    comment_out=''
    coords=[]
    elements=[]
    element_type=[]
    number_atoms_per_group=[]
    rexcl=[]
    unknown_variable_2=[]
    uc_params=[]
    super_cell=[1,1,1]
    number_atoms=0
    number_of_groups=0
    read_params= False
    def open_file(self):
        cssr_file_obj= open(self.filename_in,'r')
        #line = cssr_file_obj.readline()
        line_counter=-1
        for line in cssr_file_obj:
            if not line.strip():
                break
            line_counter+=1
            if line_counter==0:
                line_elements=shlex.split(line)
                self.uc_params.append(line_elements[0])
                self.uc_params.append(line_elements[1])
                self.uc_params.append(line_elements[2])
                self.uc_params.append(line_elements[3])
                self.uc_params.append(line_elements[4])
                self.uc_params.append(line_elements[5])
                self.super_cell.append(line_elements[6])
                self.super_cell.append(line_elements[7])
                self.super_cell.append(line_elements[8])


            if line_counter==1:
                line_elements=shlex.split(line)
                self.number_atoms=int(line_elements[0])
                self.number_of_groups=int(line_elements[1])
                self.read_params=True

            if self.read_params:
                if line_counter>1:
                    line_elements=shlex.split(line)
                    self.element_type.append(line_elements[0])
                    self.number_atoms_per_group.append(line_elements[1])
                    self.rexcl.append(line_elements[2])


            if self.read_params:
                test=self.number_of_groups+1;
                if line_counter>test:
                    line_elements=shlex.split(line)
                    self.elements.append(line_elements[1])
                    for j in range(2,5):
                        self.coords.append(float(line_elements[j]))
                    if (len(line_elements) >5):
                        self.unknown_variable_2.append(line_elements[5])

        cssr_file_obj.close()

    def write(self):
        cssr_file_obj_out= open(self.filename_out,'w')
        print>>cssr_file_obj_out,self.uc_params[0],self.uc_params[1],self.uc_params[2],\
            self.uc_params[3],self.uc_params[4],self.uc_params[5],self.super_cell[0],self.super_cell[1],self.super_cell[2]

        print>>cssr_file_obj_out,self.number_atoms,self.number_of_groups
        for group in range(self.number_of_groups):
             print>>cssr_file_obj_out, '{0:3s}'.format(self.element_type[group]), '{0:3s}'.format(self.number_atoms_per_group[group]),self.rexcl[group]

        for atom in range(self.number_atoms):
            x=atom*3;y=atom*3+1;z=atom*3+2
            #MCCCS reads it in Fortran like this
            #(i5,1x,a4,2x,3(f9.5,1x))
            print>>cssr_file_obj_out,'{0:5d}'.format(atom+1),\
            '{0:>4s}'.format(self.elements[atom])," ",\
            '{:9.5f}'.format(self.coords[x]),'{:9.5f}'.format(self.coords[y]),'{:9.5f}'.format(self.coords[z])
        cssr_file_obj_out.close()
