# ----------------------------------------------------------------------------
# Created By    : Ilia Lecha, Marc Romera and Marçal Vazquez.
# Contributions : Jose Luis Gelpi, Irene Acero, Alberto Meseguer
# Date          : 04/11/2022
# version       = '1.0'
# ---------------------------------------------------------------------------
import math
import tqdm
from Bio.PDB import Selection

# Possible Atom names that correspond to Ala atoms"
ala_atoms = {'N', 'H', 'CA', 'HA', 'C', 'O', 'CB',
             'HB', 'HB1', 'HB2', 'HB3', 'HA1', 'HA2', 'HA3'}


class Energies():
    def residue_id(self, res):
        '''Returns readable residue id'''
        return '{} {}{}'.format(res.get_resname(), res.get_parent().id, res.id[1])

    def atom_id(self, at):
        '''Returns readable atom id'''
        return '{}.{}'.format(self.residue_id(at.get_parent()), at.id)


    '''
        This function computes electrostatic, vdw and solvation energies 
        of the interaction surface.
    '''
    def calc_surface_energy(st, chain_A, srf_A, chain_E, srf_E):
        elec = 0.
        vdw = 0.
        solv = 0
        # For all resA € Residues in chain A
        for resA in Selection.unfold_entities(chain_A, 'R'):
            # Residue is on the interaction surface.
            if resA.get_id()[1] in srf_A:
                # For all resE € Residues in chain E
                for resE in Selection.unfold_entities(chain_E, 'R'):
                    # Residue is on the interaction surface.
                    if resE.get_id()[1] in srf_E:
                        # Compute energy between all atoms in resA and resE.
                        solv += Energies.calc_solvation2(st, resA)[0]
                        solv += Energies.calc_solvation2(st, resE)[0]

                        for atmA in Selection.unfold_entities(chain_A[resA.get_id()[1]], 'A'):
                            for atmE in Selection.unfold_entities(chain_E[resE.get_id()[1]], 'A'):
                                r = atmA-atmE
                                e = Energies.elec_int(atmA, atmE, r)
                                elec += e
                                e = Energies.vdw_int(atmA, atmE, r)
                                vdw += e          
        return elec, vdw, solv

    '''
        This function computes electrostatic and vdw energies of a given
        residue against ALL atoms in the PDB structure.
        Given a residue it returns electrostatic and vdw energies for all
        residues but alanine and all electrostatic and vdw energies for
        alanine residues only.
    '''
    def calc_int_energies2(st, res):
        elec = 0.
        elec_ala = 0.
        vdw = 0.
        vdw_ala = 0.

        for at1 in res.get_atoms():
            for at2 in st.get_atoms():
                # skip same chain atom pairs
                if at2.get_parent().get_parent() != res.get_parent():
                    r = at1 - at2
                    e = Energies.elec_int(at1, at2, r)
                    elec += e
                    if at1.id in ala_atoms:  # GLY are included implicitly
                        elec_ala += e
                    e = Energies.vdw_int(at1, at2, r)
                    vdw += e
                    if at1.id in ala_atoms:  # GLY are included implicitly
                        vdw_ala += e
        return elec, elec_ala, vdw, vdw_ala

    '''
        This function computes solvation energy based on ASA.
        Given a residue it returns both solvation for any residue but alanine
        and solvation for alanine residues only.
    '''
    def calc_solvation2(st, res):
        solv = 0.
        solv_ala = 0.
        for at in res.get_atoms():
            if 'EXP_NACCESS' not in at.xtra:
                continue
            s = float(at.xtra['EXP_NACCESS']) * at.xtra['vdw'].fsrf
            solv += s
            if at.id in ala_atoms:
                solv_ala += s
        return solv, solv_ala

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
