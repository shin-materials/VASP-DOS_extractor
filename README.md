# VASP_DOS_extractor
------------------
### Simple description:

Pymatgen-based python script to extract density of states (DOS) and projected DOS from vasprun.xml file
  
### What does this do:
This script replaces one of the most common use of p4vasp: extracting DOS data from VASP (http://cms.mpi.univie.ac.at/vasp/) output, which is a density functional theory (DFT) calculation program. Though p4vasp is with intuitive GUI, extracting DOS data can easily take time especially when the number of atom, band, and energy grid are larger. Note that you should install pymatgen (http://pymatgen.org/) before using this script.
  
### The process with p4vasp usually happens as follows:
  - Download xml file to local computer ***(often over 100 MB)***
  - Open p4vasp and open xml file
  - Plot element or individual atom of interest to plot the local projection of DOS (PDOS). ***(This process takes most of the time)***
  - Export data to designated filename ***(This file is usually less than 1MB)***

### With DOS_extractor.py, these four process can be done with one command line:
```
  $ python DOS_extractor.py [xml_filename] [out_filename] [entries_or_options]
 ```
  </br>
***[xml_filename]***: name of vasprun.xml file.</br>
***[out_filename]***: name of DOS data file. By default it follows the printing format of p4vasp, except for the line1 containing header. Header starts with # showing the Band gap, and list of data in the output file.</br>
***[entries]***: entry can be element or specific atom. The script also supports orbital projection.</br>
Example) Fe-d: d-orbitals of Fe, Fe1: DOS of Fe1 atom, Fe-dxy: dxy-orbital of Fe, Fe1-dxy: dxy-orbital of Fe 1 atom, O-px: p-orbital of O. </br>
***[options]***: option can be stated with --. At this moment there are three options,</br>
--elements: </br>
'--atoms', '--block'.


### Example case
For example, if you want to extract PDOS of all elements in Sr2Fe2O5, and Fe1 atom, you might do the command as follows.</br>
  **$ python DOS_extractor.py vasprun.xml Sr2Fe2O5.dat Sr Fe O Fe1</br>**
  
  
### The output file, Sr2Fe2O5.dat will be look like as follows (basically same with p4vasp format):

'# BandGap: 0.348, Label: Spin up/down for Total, Sr, Fe, O, Fe1<br/>
-59.5136	0.0000  <-- This is start of the first block<br/>
-59.4556	0.0000<br/>
.<br/>
.<br/>
.<br/>
9.9503	-0.0000<br/>
10.0082	-0.0000  <-- This is end of the first block<br/>
<br/>
-59.5136	0.0000  <-- This is start of the second block<br/>
-59.4556	0.0000<br/>
.<br/>
.<br/>
<br/>
Left column of each block is energy grid, and the right column is the DOS.
In this example, there will be 10 blocks with Total DOS (up-spin/down-spin), Sr (up-spin/down-spin), Fe (up-spin/down-spin), O (up-spin/down-spin), and Fe1 (up-spin/down-spin). If the system was non-magnetic calculation, there will be 5 blocks by not including down-spin DOS.
