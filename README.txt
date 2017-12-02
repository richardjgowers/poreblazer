================================================
poreblazer_v3.0.2: Computational structure characterization tools
Copyright (C) 2012 Lev Sarkisov

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
================================================

Welcome to poreblazer_v3.0.2, 07-03-2013.

Lev.Sarkisov@ed.ac.uk
================================================


1. About the program:

This is a new, integrated version of the computational structure 
characterization tools originally presented in 

------------------------------------------------
Sarkisov and Harrison, Computational structure characterisation tools in 
application to ordered and disordered porous materials, Molecular 
Simulation, Vol 37, Issue 15, pages 1248-1257, 2011.
------------------------------------------------

Please, cite this article if you use poreblazer_v3.0.2. 

2. Poreblazer version 3.x.x vs 1.2.x
 
In the previous versions of poreblazer, a separate code/tool existed for the 
calculation of the accessible surface area, pore size distribution and 
other properties. It required separate runs, separate inputs and 
separate compilation for each property.

The new generation of poreblazer is different:

One code, one run. All you need to specify is the name of the file 
containing the xyz coordinates of the structure and the dimensions 
(length and angles) of the unit cell.

For a given structure poreblazer_v3.0 will calculate:

- Accessible surface area (with Nitrogen atom as a probe)
- Geometric pore size distribution
- Structure density
- Helium pore volume
- Maximum pore size
- Pore limiting diameter

3. Poreblazer v3.0.2 vs v3.0

In the updated version of poreblazer several modifications are introduced
to increase the efficiency and simplicity of the program:

- Hoshen-Kopelman algorithm is properly implemented to reduce unnecessary
  redundant cycles (substantial modifications of  percolation.f90)

- for non-orthorhombic cells, calculations are now done on non-orthorhombic 
  lattice: this allows to use correct dimensions of the non-orthorhombic 
  unit cell and simplifies a number of calculations involving 
  fundamental cell subroutine 

- positions of the structure atoms in the slanted representation are pre-calculated 
  which now saves a bit of time

- calculation of the surface area using potential surface is hard-coded now as
  it has become the standard. No hard sphere surface option! The default parameter
  file is slightly different now to reflect this change
 
- other optimizations throughout


4. Contents of the distribution:

The distribution contains the original source code files written in Fortran 90, 
compilation scripts and pre-compiled static (library independent) version of 
the program (poreblazer.exe).

4.1 Executable poreblazer.exe was compiled on Intel Westmere E5620 processor, 
using gfortran compiler GNU Fortran (GCC) 4.4.6 20110731 (Red Hat 4.4.6-3). 
This executable should be compatible with most Intel platforms. If you wish to compile your 
own version go to the next step.


The Makefile is also provided. To change compiler from gfortran to intel, in the line

FORTRAN_COMPILER= gfortran

change gfortran to the intel compiler executable name on your system. 

4.2 To compile the code you can simply issue a "make" command, which 
will follow the instructions in the Makefile to assemble the code into the
poreblazer.exe executable file.

4.3 Other compilers. The code has been tested with Intel fortran compiler and 
gfortran. Compilation and testing using other compilers is at the 
discretion of the users. The mutual file dependency is summarized at the 
end of the Makefile for custom compilations.

4.4 Windows. The code has been also compiled for Windows using gfortran under
MinGW environment. The static executable file is located in WinXP. Also, in the same
directory an example of using the program under Windows is provided. Change you location
to HKUST1, and double click on run.bat script. This will run the program (located
one directory up from your test case location). The results will be saved  in
results.txt.

The distribution also contains the following case studies:

HKUST1
IRMOF1
MOF180
MIL47V

with reference results and performance.

5. How to run

5.1 Basic mode

In basic mode the program uses default parameter setting. For this, in the 
location of the run, you need to put files defaults.dat and UFF.atoms. 
They can be copied from the case studies provided with the 
distribution and are the same for all runs. All that a user needs to 
specify is the name of the file containing the xyz coordinates of the 
structure (which uses the standard names for the elements, H, C, 
O etc) and dimensions of the unit cell (side lengths and yz, xz, xy angles. 
For example, for HKUST1, it is specified in the input.dat file:

HKUST1.xyz
26.3430000000000        26.3430000000000        26.3430000000000
90             90             90

To run simply issue:

./poreblazer.exe < input.dat

Examples of the output are summarized in the provided case studies. 

5.2 Advanced mode

In the advanced mode, a user controls various parameters of the simulation, 
such as the interaction parameters of Nitrogen and Helium atoms, cut-off 
radius, grid size etc. All the parameters are specified in the defaults.dat 
file and can be changed there. This file also specifies a source of the 
force field parameters for the atoms of the structure. By default it is the 
UFF forcefield, provided in UFF.atoms file (note that not all possible 
atoms from the UFF are there, but the most common ones). Other force fields 
and files can be used, with the customized convention for the atom names 
defined by the user. 

This mode assumes the user knows what s/he does. Here is an example of the 
defaults.dat file:

UFF.atoms                  
2.58, 10.22, 3.314, 298    
12.8, 500           
0.2
20.0, 0.25
21908391 

These values are (in order of their appearance):

- Default forcefield: UFF
- Helium atom sigma (A), helium atom epsilon (K), nitrogen atom sigma (A), temperature (K)
Cutoff distance (A),  number of trials for surface area 
calculation
- 0.2: Cubelet size (A)
- Largest anticipated pore diameter (A), size of the bin for PSD (A)
- Random number seed

 
6. Understanding and interpreting the results

Results are printed on the screen, and can be redirected to a text file. 
PSD is saved in the psd.txt file.

The results for the references= cases are included in the distribution. 
The general convention for the names of the reference result files is as follows:

results_intel: results from intel compiler
results_gfort: results from gfortran compiler
psd.txt_intel/gfort: pore size distribution from intel/gfortran compiler
time_cpu_intel/gfortran: cpu time for intel/gfortran compiler


These files are located in each of the respective material folder (HKUST1,
IRMOF1, MOF180, MIL47)

Below we summarize some key results for these case studies as a function 
of the cubelet size. The data is for hard-sphere surface area, crystal density,
Helium pore volume, maximum pore size, pore limiting diameter, and CPU (Linux)
time of the calculation. 

Note that only the maximum pore size, pore limiting diameter, CPU time and pore
size distribution depend on the cubelet size.



                   SA (m2/g)     rho (g/cm3)  Vpore (He, cm3/g)  Dmax (A)  PLD    CPU(gfortran) v.3.0.2 (v3.0)  
--------------------------------------------------------------------------------------------------------------
HKUST1             1910.59       0.879        0.853              12.74     6.37  1m32.725s (3m3.888s)   
--------------------------------------------------------------------------------------------------------------
IRMOF1             3378.17       0.604        1.353              14.85     7.65  1m1.104s  (1m47.379s)        
--------------------------------------------------------------------------------------------------------------
MOF180             5829.66       0.271        3.378              15.22     12.80 31m11.043s (-)
--------------------------------------------------------------------------------------------------------------
MIL47              1324.04       1.000        0.629              7.56      7.01  3m2.679s  (4m26.634s)    



Comparison with v1.2

IRMOF1

                   SA (m2/g)     rho (g/cm3)  Vpore (He, cm3/g)  Dmax (A)  PLD    CPU(gfortran)   
--------------------------------------------------------------------------------------------------------------
v1.2               3312.75       0.593        1.363              -         7.25   hours   
--------------------------------------------------------------------------------------------------------------
v3.02              3378.17       0.604        1.353              14.85     7.65   1m1.104s  (1m47.379s)        


MOF180

                   SA (m2/g)     rho (g/cm3)  Vpore (He, cm3/g)  Dmax (A)  PLD    CPU(gfortran)  
--------------------------------------------------------------------------------------------------------------
v1.2               5940.20       0.267        3.2941              -        12.25   hours   
--------------------------------------------------------------------------------------------------------------
v3.02              5829.66       0.271        3.378              14.85     12.80   1m1.104s  (1m47.379s)        


MOF180 is an important case study. The results provided in the 
results files show that the calculated density is higher than the 
literature value and as a consequence the reported surface 
area is lower than the literature value. 

The reason for this that the MOF180.xyz structure contains additional 
unresolved duplicates of atoms, making the structure denser. 
This provides an illustration that files obtained directly from various 
crystallographic packages must be used carefully to ensure that 
the xyz file generated from CIF file contains an appropriate 
number of atoms and corresponds to a correct crystal density.

Using HKUST1 as an example, it seems:

- ATEN gives an XYZ file with the correct number of atoms (624)
- VESTA and CrystalMaker give XYZ files with 648 atoms. This is 
because the 8 Cu atoms on each unit cell face are double-counted.
- Mercury gives an XYZ file with way too many atoms, because it extends 
the displayed range of the structure well beyond the unit cell boundaries

(We would like to thank Dr. Michael Fischer, UCL, for carrying out these
important tests).

5.  Limitations/other

- The current upper pore diameter is set to 20A. 
Change it in the defaults.dat if you anticipate larger pores 
in the system

- The accuracy of the PSD, pore limiting diameter etc is governed 
by the default lattice size currently set at 0.2A. 

- Density of the system is calculated from the xyz file of the structure. 
Make sure that the file and the specified dimensions of the cell are correct.



-----------------------------------------




As of March 2013, poreblazer_v3.0.2 is in the testing stage. 

Please, send your comments and corrections to

Lev.Sarkisov@ed.ac.uk

