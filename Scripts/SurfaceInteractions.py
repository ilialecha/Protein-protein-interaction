#----------------------------------------------------------------------------
# Created By    : Ilia Lecha, Marc Romera and Marçal Vazquez.
# Contributions : Jose Luis Gelpi, Irene Acero, Alberto Meseguer
# Date          : 04/11/2022
# version       = '1.0'
# --------------------------------------------------------------------------- 
from Bio.PDB import NeighborSearch, Selection
from os import system

class SurfaceInteractions():
    def __init__(self,chain_A,chain_E,distance):
        self.chain_A   = chain_A
        self.chain_E   = chain_E
        self.distance   = distance

    def get_atoms_in_contact(self,surface_neighbors):
        '''
            This function return the atoms at a distance lower or equal to a given number of amstrongs.
        '''
        atoms_in_contact = {}
        for k in surface_neighbors.keys():
            E_residues = [self.chain_E[resi] for resi in surface_neighbors.get(k)]
            atoms_E_list = Selection.unfold_entities(E_residues, 'A')
            for atom in self.chain_A[k].get_atoms():
                for atom_E in atoms_E_list:
                    if (atom - atom_E) <= self.distance:
                        if not atom.get_name() in atoms_in_contact:
                            atoms_in_contact[atom.get_name()] = [atom_E.get_name(),atom - atom_E]
                        else:
                            atoms_in_contact[atom.get_name()] += [atom_E.get_name(),atom - atom_E]
        return atoms_in_contact

    def get_neighbors(self):
        '''
            This function returns a dictionary where the keys are the residues from the chain_A and 
            the items are its neighbors given a distance in amstrongs.
        '''
        chain_A_residues  = Selection.unfold_entities(self.chain_A, 'R')# Obtaining el residues on chain A
        chain_E_atoms     = Selection.unfold_entities(self.chain_E, 'A')# Obtaining all atoms in chain E
        neighbor_search   = NeighborSearch(chain_E_atoms)               # Setting up all atoms on chain E as a reference for the search.  
        surface_neighbors = {}

        for res in chain_A_residues:
            for atom in res.get_atoms():
                residue_list = Selection.unfold_entities(neighbor_search.search(atom.get_coord(), self.distance), 'R') # R for residues
                if len(residue_list)>0:
                    surface_neighbors[res.get_id()[1]] = set(res.get_id()[1] for res in residue_list)
        return surface_neighbors  
