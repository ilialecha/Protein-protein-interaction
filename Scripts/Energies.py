#----------------------------------------------------------------------------
# Created By    : Ilia Lecha, Marc Romera and Marçal Vazquez.
# Contributions : Jose Luis Gelpi, Irene Acero, Alberto Meseguer
# Date          : 04/11/2022
# version       = '1.0'
# --------------------------------------------------------------------------- 
import math
from Bio.PDB import Selection

#Possible Atom names that correspond to Ala atoms"
ala_atoms = {'N', 'H', 'CA', 'HA', 'C', 'O', 'CB', 'HB', 'HB1', 'HB2', 'HB3', 'HA1', 'HA2', 'HA3'}

class Energies():
    def residue_id(self,res):
        '''Returns readable residue id'''
        return '{} {}{}'.format(res.get_resname(), res.get_parent().id, res.id[1])

    def atom_id(self,at):
        '''Returns readable atom id'''
        return '{}.{}'.format(self.residue_id(at.get_parent()), at.id)

    def calc_int_energies(st, res):
        '''Returns interaction energies (residue against other chains)
            for all atoms and for Ala atoms
        '''
        elec  = 0.
        vdw   = 0.
        solv  = 0.
        for at1 in res.get_atoms():
            for at2 in st.get_atoms():
            # skip same chain atom pairs
                if at2.get_parent().get_parent() != res.get_parent():
                    r = at1 - at2
                    elec += Energies.elec_int(at1, at2, r)
                    vdw += Energies.vdw_int(at1, at2, r)          
        return elec, vdw, Energies.calc_solvation(res)


    def calc_int_energies_ala(st, res):
        '''Returns interaction energies (residue against other chains)
            for all atoms and for Ala atoms
        '''
        elec_ala = 0.
        vdw_ala = 0.
        solv_ala = 0.

        for at1 in res.get_atoms():
            for at2 in st.get_atoms():
            # skip same chain atom pairs
                if at2.get_parent().get_parent() != res.get_parent():
                    r = at1 - at2
                    if at1.id in ala_atoms:
                        elec_ala += Energies.elec_int(at1, at2, r) 
                        vdw_ala  += Energies.vdw_int(at1, at2, r) 
                        
        return elec_ala, vdw_ala

    def MH_diel(r):
        '''Mehler-Solmajer dielectric'''
        return 86.9525 / (1 - 7.7839 * math.exp(-0.3153 * r)) - 8.5525

    def elec_int(at1, at2, r):
        '''Electrostatic interaction energy between two atoms at r distance'''
        return 332.16 * at1.xtra['charge'] * at2.xtra['charge'] / Energies.MH_diel(r) / r

    def vdw_int(at1, at2, r):
        '''Vdw interaction energy between two atoms'''
        eps12 = math.sqrt(at1.xtra['vdw'].eps * at2.xtra['vdw'].eps)
        sig12_2 = at1.xtra['vdw'].sig * at2.xtra['vdw'].sig
        return 4 * eps12 * (sig12_2**6/r**12 - sig12_2**3/r**6)
    
    def calc_solvation(res):
        '''Solvation energy based on ASA'''
        '''Return solv.Energy of a given residue, have to execute this for all residues.'''
        solv = 0.
        for at in res.get_atoms():
            if 'EXP_NACCESS' not in at.xtra: continue
            else: solv += float(at.xtra['EXP_NACCESS'])* at.xtra['vdw'].fsrf
        return solv
    
    def calc_solvation_ala(res):
        '''Solvation energy based on ASA'''
        '''Return solv.Energy of a given residue, have to execute this for all residues.'''
        solv_ala = 0.
        for at in res.get_atoms():
            if 'EXP_NACCESS' not in at.xtra: continue
            else:
                if at.id in ala_atoms: solv_ala += float(at.xtra['EXP_NACCESS'])* at.xtra['vdw'].fsrf     
        return solv_ala
