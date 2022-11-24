# Protein-protein-interaction

This set of python scripts is mainly prepared to work with a PDB file with a single structure containing
two main chains.It has been firstly developed to analyze the protein-to-protein interaction between SPIKE
Protein and ACE2.
For more information about the PDB structure being used in this project:
https://www.rcsb.org/structure/6M0J


## PdbChecker.py
Python script used to remove some atoms/molecules of a given PDB structure.

## RamachandranPlot.py
Python script used to create a Ramachandran Plot of a given PDB structure.

## SurfaceInteractions.py
Python script used to detect which residues and/or atoms conform the 
interaction surface between two protein chains of a given PDB structure.

## VDWParameters.py
Python script used to process aaLib and vdwprm and add the required extra parameters to all atoms within the structure.

## Energies.py
Python script used to compute the interaction energy of two protein chains of 
a given PDB structure.

## P2PAnalysis.py
Python Script that imports an executes all required code of the python scripts above.



