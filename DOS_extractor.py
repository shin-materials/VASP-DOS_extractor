# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 2018
@author: Yongjin Shin
Research Group: Materials Design and Theory Group
Advisor: Prof. James Rondinelli
Affiliation: Northwestern University

Command:
 $ python DOS_extractor.py [xml_filename] [out_filename] [entries_or_options]
-[xml_filename]: name of vasprun.xml file.
-[out_filename]: name of DOS data file. By default it follows the printing format of p4vasp. For spin polarized system, spin up/down
-p4vasp format: This format is collection of data arrays. Each data array represents each entry, including total DOS, which are composed of two columns; column1 for energy and column2 for DOS of entry. line1 is header starting with #, which shows the Band gap, and list of data in the output file.
block format: This format is intuitive structured data with row of energy grid and column of each entry. The line1 contains labels.
-[entries]: Entry can be element or specific atom. Atom label follows the VESTA, for example atoms in BaTiO3 unit cell would be Ba1, Ti1, O1, O2, O3. The script also supports orbital projection by specifiying orbital after dash (-). Examples and nomenclatures are as follows.
Ex) Fe: DOS of all Fe atoms
Ex) Fe-d: d-orbitals of Fe
Ex) Fe-dxy: dxy-orbital of Fe
Ex) Fe1: DOS of Fe1 atom
Ex) Fe1-dxy: dxy-orbital of Fe1 atom
Ex) O-px: px-orbital of O
Note that these nomenclatures are basically the order in which the orbitals are reported in VASP and has no special meaning.
p-orbitals: px, py, pz
d-orbitals: dxy, dyz, dz2, dxz, dx2
f-orbitals: f_3, f_2, f_1, f0, f1, f2, f3
-[options]: option can be stated with --. At this moment there are three options.
--elements: Include all elements in the system to the entries.
--atoms: Include all individual atoms in the system to the entries.
--block: Change printing option to block data. Otherwise the printing option will follow p4vasp.

Maintenance update
180816: struct.sites is not recognized as pdos keys. So I change the struct=dos_vrun.structures[-1] to dos_vrun.final_structure
Other way to resolve this, I make list_sites=list(pdos.keys()) and use it to get pdos. 
To make sure if list_sites has the same order with struct.sites, I compare their ._fcoords and if they are different print warning message.
"""
#from pymatgen.electronic_structure import dos
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.core.periodic_table import Element
from pymatgen.electronic_structure.core import Spin, Orbital, OrbitalType
import numpy as np
import copy
import sys
import re

#########################################################################
###########   PART1: PREPARATION. LOAD FILE AND READ DOS ################
#########################################################################

### System arguments: XML file name, output filename ####
xml_filename=sys.argv[1]
out_filename=sys.argv[2]
##########################################################

### Read xmlfile and load dos object #####################
dos_vrun=Vasprun(xml_filename)
#get dos object
total_dos = dos_vrun.complete_dos
#get structure from vasprun
struct = dos_vrun.final_structure
#n_E
n_E=len(total_dos.energies)
#band_gap
band_gap=total_dos.get_gap()
##########################################################

### Make dictionary: from label to site index ############
# For easy use, Labe them in conventional manner. labeled as Element + site position in poscar (i.e. BaTiO3 --> Ba1, Ti2, O1, O2, O3)
# Plus, label2site dictionary will connect atomic label to site number
# Example: BaTiO3 --> Ba1:0 Ti1:1 O1:2 O2:3 O3:4
n_atom_count_dict=dict()
label2site_index=dict()
list_atoms=[]
for i in range(0,struct.num_sites):
    # Update label for each element
    if struct[i].specie in list(n_atom_count_dict.keys()):
        n_atom_count_dict.update({struct[i].specie:n_atom_count_dict[struct[i].specie]+1})
    else:
        n_atom_count_dict.update({struct[i].specie:1})
    label2site_index.update({'{0}{1}'.format(struct.species[i],n_atom_count_dict[struct[i].specie]):i})
    list_atoms.append('{0}{1}'.format(struct.species[i],n_atom_count_dict[struct[i].specie]))
##########################################################

### Define ISPIN #########################################
# ISPIN is based on the size of total_DOS
ispin=len(total_dos.densities)
# spin_list will be used when printing dos
spin_list=[Spin.up, Spin.down]
##########################################################

### Get DOS of each element ##############################
element_dos=total_dos.get_element_dos()
# Usage: element_dos[Element["Fe"]].densities[Spin.up]
##########################################################

### GET SITE DOS #########################################
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
##########################################################

### Function for getting orbital densities ###############
def orbital_dos(object_orbital):
    """
    input: 
        object_orbital - str, represented as either 'element-orbital' or 'site-orbital'.
        Ex) Fe-d, O-p, Fe1-d, Fe1-dx2, O-pz
    output:
        dos_dict - dictionary, key:Spin.up, Spin.down (if ispin=2), value: density_array with size (n_E,)
    external_variables:
        pymatgen Orbital, OrbitalType, struct, total_dos
        label2site_index, ispin, spin_list
    """
    object_str,orbital_str=re.split('-',object_orbital)
    #checking last letter of object will determine whether site or element
    if object_str[-1].isnumeric():
        #case: site
        # length of orbital_str can distinguish whether specific orbital or orbital type
        if len(orbital_str)==1:
            # case for orbital type
            site_spd_dos=total_dos.get_site_spd_dos(struct.sites[label2site_index[object_str]])
            dos_dict=site_spd_dos[OrbitalType[orbital_str]].densities        
        else:
            # orbital is one specific orbital. Ex) s, px, py, dx2, dxy, dxz
            site_orbital_dos=total_dos.get_site_orbital_dos(struct.sites[label2site_index[object_str]],Orbital[orbital_str])
            dos_dict=site_orbital_dos.densities
    else:
        #case: element
        if len(orbital_str)==1:
            # case for orbital type
            elemental_spd_dos=total_dos.get_element_spd_dos(Element[object_str])
            dos_dict=elemental_spd_dos[OrbitalType[orbital_str]].densities
        else: 
            # orbital is one specific orbital. Ex) s, px, py, dx2, dxy, dxz
            # Caution: there is no specific function in pymatgen, so I import all atoms and sum them up
            dos_dict=dict()
            dos_dict={Spin.up:np.zeros(shape=(len(total_dos.energies)))}
            # if ispin==2, entry for down spin is needed
            if ispin == 2:
                dos_dict.update({Spin.down:np.zeros(shape=(len(total_dos.energies)))})
            # now sum inidivdual dos of designate element
            for site_index in list(struct.indices_from_symbol(object_str)):
                site_orbital_dos=total_dos.get_site_orbital_dos(struct.sites[site_index],Orbital[orbital_str])
                for j in range(ispin):
                    spin=spin_list[j]
                    #update for each spin
                    temp_array=dos_dict[spin]+site_orbital_dos.densities[spin]
                    dos_dict.update({spin:temp_array})
    return dos_dict
# Usage example: dos.densities[Spin.up]
#############################################################

#########################################################################
###########   PART2: PRINTING FORMATS. BLOCK OR P4VASP ##################
#########################################################################

##### Block Print ###########################################
def dos_block_print(list_entries, total_dos, out_filename):
    """
    input:
        list_entries - list, elements are entries with object-orbital
            Ex) 'Sr', 'Fe', 'O', 'Fe-d', 'O1-px', 'Fe1-dx2', etc.
    print_output:
        data_array : numpy_array with shape (n_E,num_col)
             case: ispin1
             col(0): E_grid, col(1): Total(up), col(2-N+1): n_th entry
             case: ispin2
             col(0): E_grid, col(1): Total(up), col(2): Total(down), Col(3-2N+2): n_th entry (up or down)
             Total num_col is 1(E_grid) + ispin*(1(total)+number_of_entries)
    external_variables:
        pymatgen total_dos
        label2site_index, n_E, ispin, spin_list, list_printed
    """
    # Reason to copy total_dos is because numpy.reshape makes total_dos object somewhat dirty
    total_dos1=copy.deepcopy(total_dos)
    data_array=total_dos1.energies.reshape(n_E,1)
    data_array[:, 0] = data_array[:, 0] - total_dos1.efermi
    # Append total_dos first
    for j in range(ispin):
        spin = spin_list[j]
        temp_array = copy.copy(total_dos1.densities[spin])
        temp_array=float(int(spin))*temp_array
        data_array=np.concatenate((data_array,temp_array.reshape(n_E,1)),axis=1)
    # Loop for each entry                    
    for entry in list_entries:
        # case divided depending on the entry type
        if '-' not in entry and not entry[-1].isnumeric():
            #case: pure element, all orbital
            dos_dict=element_dos[Element[entry]].densities
        elif '-' not in entry and entry[-1].isnumeric():
            #case: one site, all orbital
            dos_dict=site_dos[entry].densities
        else:
            #case: orbital defined in entry
            dos_dict=orbital_dos(entry)
        # Printing for entries
        # Loop for available spins in entry
        for j in range(ispin):
            spin=spin_list[j]
            temp_array=copy.copy(dos_dict[spin])
            temp_array=float(int(spin))*temp_array
            data_array=np.concatenate((data_array,temp_array.reshape(n_E,1)),axis=1)
    ### Printing part #####
    out_file=open(out_filename,'w')
    ## Labels
    np.set_printoptions(precision=4)
    np.set_printoptions(formatter={'float': '{:0.4f}'.format})
    list_label=['E-E_f']
    spin_label=['(Up)','(Dn)']
    for j in range(ispin):
        list_label.append('Total'+spin_label[j])
    for entry in list_entries:
        for j in range(ispin):
            list_label.append(entry+spin_label[j])
    out_file.write(''.join('%-11s' % entry for entry in list_label)+'\n')
    np.savetxt(out_filename, data_array,fmt='% -11.4f',delimiter='')
    out_file.close()
    return None
#############################################################
   
### Print with p4vasp data format ###########################
def dos_p4v_print(list_entries, total_dos, out_filename1):
    """
    input:
        list_entries - list, elements are entries with object-orbital
            Ex) 'Sr', 'Fe', 'O', 'Fe-d', 'O1-px', 'Fe1-dx2', etc.
    print_output: p4v format for each entry, Energy_grid and DOS_value arrays
    external_variables:
        pymatgen total_dos
        label2site_index, n_E, ispin, spin_list, list_printed
    """
    total_dos1=copy.deepcopy(total_dos)
    # Print header, including band gap, and list of printed DOS ######
    out_file1=open(out_filename1,'w')
    #out_file1.write('# Total(Up, Down), Sr(Up,Down), Fe(Up,Down), O(Up,Down), Fe1(Up,Down)\n')
    list_label=[]
    spin_label=['(Up)','(Dn)']
    for j in range(ispin):
        list_label.append('Total'+spin_label[j])
    for entry in list_entries:
        for j in range(ispin):
            list_label.append(entry+spin_label[j])
    out_file1.write('# BandGap: {0:.3f}eV, Label: '.format(band_gap)+', '.join(list_label)+'\n')
    # print total dos first
    for j in range(ispin):
        # spin is Spin.up or Spin.down
        spin=spin_list[j]
        for i in range(n_E):
            out_file1.write('{0:.4f}\t{1:.4f}\n'.format(total_dos1.energies[i]-total_dos1.efermi,\
                            float(int(spin))*total_dos1.densities[spin][i]))
        #p4vasp_format: Spacer between different array
        out_file1.write('\n')
    # Loop for each entry                    
    for entry in list_entries:
        # case divided depending on the entry type
        if '-' not in entry and not entry[-1].isnumeric():
            #case: pure element, all orbital
            dos_dict=element_dos[Element[entry]].densities
        elif '-' not in entry and entry[-1].isnumeric():
            #case: one site, all orbital
            dos_dict=site_dos[entry].densities
        else:
            #case: orbital defined in entry
            dos_dict=orbital_dos(entry)
        # Printing for entries
        # Loop for available spins in entry
        for j in range(ispin):
            spin=spin_list[j]
            for i in range(n_E):
                out_file1.write('{0: .4f}\t{1:.4f}\n'.format(total_dos1.energies[i]-total_dos1.efermi,\
                                float(int(spin))*dos_dict[spin][i]))
            #p4vasp_format: Spacer between different array
            out_file1.write('\n')
    out_file1.close()
    return None
#############################################################

#########################################################################
###########   PART3: DEFINE ENTRIES AND PRINT OUT #######################
#########################################################################

 
### Construct entries and options ############################
list_entries=[]
#default for p4v format
print_function=dos_p4v_print
# Make list of entries
for i in range(3,len(sys.argv)):
    #option arguments
    if sys.argv[i][:2] == '--':        
        if sys.argv[i] == '--elements':
            list_entries=list_entries+list(struct.symbol_set)
        elif sys.argv[i] == '--atoms':
            list_entries=list_entries+list_atoms
        
        if sys.argv[i] == '--block':
            print_function=dos_block_print
    # normal entry arguments
    else:
        list_entries.append(sys.argv[i])
##############################################################

### Print on file ############################################

print_function(list_entries,total_dos,out_filename)
##############################################################

### Print out the summary #####################################
print ('Interpolated Band gap: {0:.4f} eV'.format(band_gap))
if ispin==1:
    print('Only up-spin is printed for following components')
else:
    print('Spin up/down is printed for following components')
print ('--Entries: '+', '.join(list_entries))
##############################################################

