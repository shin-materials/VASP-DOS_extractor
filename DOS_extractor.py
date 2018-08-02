# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 22:12:08 2018

@author: Yongjin
0717: Error fixing + energy offset to Fermi energy
0718: Implement ISPIN 1 and ISPIN 2 + Printing band gap

System arguments: XML file name, output filename, Type of elements
"""
#from pymatgen.electronic_structure import dos
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.core.periodic_table import Element
from pymatgen.electronic_structure.core import Spin
import sys
#import os
#import numpy as np

###### System arguments: XML file name, output filename, Type of elements #################
xml_filename=sys.argv[1]
out_filename1=sys.argv[2]

list_element=[]
list_sites=[]
for i in range(3,len(sys.argv)):
    if sys.argv[i][-1].isdigit():
        list_sites.append(sys.argv[i])
    else:
        list_element.append(sys.argv[i])


dos_vrun=Vasprun(xml_filename)
total_dos = dos_vrun.complete_dos

#get structure from vasprun
#from pymatgen.electronic_structure.plotter import DosPlotter
struct = dos_vrun.structures[-1]

#site dos are labeled as Element + site position in poscar (i.e. BaTiO3 --> Ba1, Ti2, O1, O2, O3)
n_atom_count_dict=dict()
site_dos = dict()
for i in range(0,struct.num_sites):
    # Update label for each element
    if struct[i].specie in list(n_atom_count_dict.keys()):
        n_atom_count_dict.update({struct[i].specie:n_atom_count_dict[struct[i].specie]+1})
    else:
        n_atom_count_dict.update({struct[i].specie:1})
    #Obtain site dos and update site_dos_dict
    site_dos.update({'{0}{1}'.format(struct.species[i],n_atom_count_dict[struct[i].specie]): total_dos.get_site_dos(struct[i])})

#### Define ISPIN #####
ispin=len(total_dos.densities)
spin_list=[Spin.up, Spin.down]

####### Get DOS of each element ############        
element_dos=total_dos.get_element_dos()
# Usage: element_dos[Element["Fe"]].densities[Spin.up]

####### Print with p4vasp data format ###########

#### Print header, including band gap, and list of printed DOS ######
out_file1=open(out_filename1,'w')
#out_file1.write('# Total(Up, Down), Sr(Up,Down), Fe(Up,Down), O(Up,Down), Fe1(Up,Down)\n')
list_printed=list_element+list_sites
if ispin==2:
    out_file1.write('# BandGap: {0:.3f}eV, Label: Spin up/down for Total, '.format(total_dos.get_gap())+', '.join(list_printed)+'\n')
else:
    out_file1.write('# BandGap: {0:.3f}eV, Label: Spin up for Total, '.format(total_dos.get_gap())+', '.join(list_printed)+'\n')
                    
# Energy grid point is obtained from site_dos of the first atom
n_E_grid=len(total_dos.get_site_dos(struct[0]).energies)

#Printing Total_up
for i in range(n_E_grid):
    out_file1.write('{0:.4f}\t{1:.4f}\n'.format(total_dos.energies[i]-total_dos.efermi,total_dos.densities[Spin.up][i]))
#p4vasp format: Spacer between different DOS
out_file1.write('\n')

#Printing Total_down
if ispin == 2:    
    for i in range(n_E_grid):
        out_file1.write('{0:.4f}\t{1:.4f}\n'.format(total_dos.energies[i]-total_dos.efermi,-1.0*total_dos.densities[Spin.down][i]))
    #p4vasp format: Spacer between different DOS
    out_file1.write('\n')

#Printing each Element
Element_list=list_element
for k_element in Element_list:
    # Loop for each spin
    for j in range(ispin):
        #if ispin==2, this will go through Spin.down
        #j_spin is either Spin.up or Spin.down
        j_spin=spin_list[j]
        # Loop for each energy grid
        for i in range(n_E_grid):
            out_file1.write('{0:.4f}\t{1:.4f}\n'.format(total_dos.energies[i]-total_dos.efermi,float(int(j_spin))*element_dos[Element[k_element]].densities[j_spin][i]))
        #p4vasp format: Spacer between different array
        out_file1.write('\n')

# Printing inividual atom.

# Loop for each spin
for individual_atom in list_sites:
    for j in range(ispin):
        #if ispin==2, this will go through Spin.down
        #j_spin is either Spin.up or Spin.down
        j_spin=spin_list[j]
        # Loop for each energy grid
        for i in range(n_E_grid):
            #Spin down DOS is printed with negative number,by using float(int(j_spin))
            out_file1.write('{0:.4f}\t{1:.4f}\n'.format(total_dos.energies[i]-total_dos.efermi,float(int(j_spin))*site_dos[individual_atom].densities[j_spin][i]))
        #p4vasp_format: Spacer between different array
        out_file1.write('\n')


out_file1.close()
print ('Interpolated Band gap is: ')
if ispin==1:
    print('Only up-spin is printed for following components')
else:
    print('Spin up/down is printed for following components')
print ('--Elements: '+str(list_element))
print ('--Atoms: '+str(list_sites))

