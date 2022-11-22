# ----------------------------------------------------------------------------
# Created By    : Ilia Lecha, Marc Romera and Marçal Vazquez.
# Contributions : Jose Luis Gelpi, Irene Acero, Alberto Meseguer
# Date          : 04/11/2022
# version       = '1.0'
# ---------------------------------------------------------------------------
import os
import getopt
from Energies import *
from VDWParameters import *
from SurfaceInteractions import *
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NACCESS import NACCESS_atomic
from Bio.PDB.PDBIO import PDBIO
from tqdm import tqdm


def man():
    print('main.py -i <inputdir> -d <distance>')
    print('Inputdir must contain the following files:')
    print('aaLib.lib:\tAtoms library')
    print('Pdb_file:\tContaining a model with two chains.')
    print('vdwprm:\t\tParameters for Van der Waals forces.')


def main(argv):
    # Setting up of the required parameters.
    INPUTDIR = ''
    DISTANCE = 0.0
    try:
        opts, args = getopt.getopt(argv, "hi:d:", ["idir=", "dist="])
    except getopt.GetoptError:
        man()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            man()
            sys.exit()
        elif opt in ("-i", "--idir"):
            INPUTDIR = arg
        elif opt in ("-d", "--dist"):
            DISTANCE = float(arg)

    if 0.0 == DISTANCE:
        man()
        sys.exit(2)
    if "/" != INPUTDIR[-1]:
        INPUTDIR += "/"

    # Necessary for adding required parameters to atoms.
    NACCESS_BINARY = INPUTDIR+"soft/NACCESS/naccess"

    print("(#) Adding Van der Waals parameters to PDB Structure.\n")

    # loading residue library from files/aaLib.lib
    res_lib = ResiduesDataLib(INPUTDIR+'aaLib.lib')

    # loading VdW parameters
    ff_params = VdwParamset(INPUTDIR+'vdwprm')

    # set the pdb_path and load the structure
    pdb_path = INPUTDIR+"6m0j_fixed.pdb"

    # Setting the Bio.PDB.Parser object
    parser = PDBParser(PERMISSIVE=1)

    # Loading structure
    st = parser.get_structure('st', pdb_path)

    # Adding VDW parameters to atoms.
    add_atom_parameters(st, res_lib, ff_params)

    srf = NACCESS_atomic(st[0], naccess_binary=NACCESS_BINARY)

    io = PDBIO()

    st_chains = {}
    # Using BioIO trick (see tutorial) to select chains

    for ch in st[0]:
        io.set_structure(st)
        io.save('tmp.pdb', SelectChain(ch.id))
        st_chains[ch.id] = parser.get_structure('stA', 'tmp.pdb')
        srfA = NACCESS_atomic(
            st_chains[ch.id][0], naccess_binary=NACCESS_BINARY)
    os.remove('tmp.pdb')

    f = open(INPUTDIR+"energies.tsv", "w")
    print("Type;Res;Electrostatic_AE;vdw_AE;solv_AE;total", file=f)

    print("(#) Identifying interaction surface residues.\n")

    chain_A = st[0]["A"]  # Obtaining chain A from model 0 --> ACE2
    chain_E = st[0]["E"]  # Obtaining chain E from model 0 --> Spike
    surfaceInt = SurfaceInteractions(chain_A, chain_E, DISTANCE)
    surfaceRes = surfaceInt.get_neighbors()
    surface_chain_A = []
    surface_chain_E = []


    # Storing into a list all residues of chain A within the surface interaction.
    for k, v in surfaceRes.items():
        surface_chain_A.append(int(k))
        [surface_chain_E.append(res) for res in v]

    print("(#) Computing energy of the interaction surface.\n")
    
    srf_energy = Energies.calc_surface_energy(st, chain_A, surface_chain_A, chain_E, surface_chain_E)
    print(f"(#) Int. Surf. energy = (Electrostatic) {srf_energy[0]} + (VDW) {srf_energy[1]} + (Solvation) {srf_energy[2]} = {srf_energy[0]+srf_energy[1]+srf_energy[2]}" ) 

    print("(#) Computing all energies: Electrostatic, Van der waals and Solvation.\n")

    # Solv_AE = Solv_A + Solv_E so we can omit some computation.
    solv_E = sum([Energies.calc_solvation2(st, res)[0]
                 for res in Selection.unfold_entities(st[0]['E'], 'R')])
    

    # Computing energies for all atoms against all atoms. More information related in Energies.py
    AAG_ = dict()
    for res in tqdm(chain_A):
        tmp = Energies.calc_int_energies2(st, res)
        tmp2 = Energies.calc_solvation2(st, res)
        #                 (elec,elec_ala,vdw,vdw_ala,solv,solv_ala)
        AAG_[res.id[1]] = (tmp[0], tmp[1], tmp[2], tmp[3], tmp2[0], tmp2[1])

    int_elec = 0.
    int_elec_ala = 0.
    int_vdw = 0.
    int_vdw_ala = 0.
    int_solv = 0.
    int_solv_ala = 0.

    for res in AAG_:
        int_elec += AAG_[res][0]
        int_elec_ala += AAG_[res][1]
        int_vdw += AAG_[res][2]
        int_vdw_ala += AAG_[res][3]
        int_solv += AAG_[res][4]
        int_solv_ala += AAG_[res][5]
    
    print(f"(#) Total energy = (Electrostatic) {int_elec} + (VDW) {int_vdw} + (Solvation A) {int_solv} + (Solvation E) {solv_E}= {int_elec+int_vdw+int_solv+solv_E}" ) 

    print(f"(#) Writting into output file.\n")
    print("Normal;", "All", ";", int_elec, ";", int_vdw, ";", int_solv +
          solv_E, ";", int_elec+int_vdw+int_solv+solv_E, file=f)

    print("(#) Running Ala-Scanning:\n")
    print(
        f"(#) For each execution, one of the following residues: \n(#) {surface_chain_A} is treated as an Alanine.\n")

    for surf_res in tqdm(surface_chain_A):
        int_elec_ala = 0.
        int_vdw_ala = 0.
        int_solv_ala = 0.
        for res in AAG_.keys():
            if surf_res == res:  # Treated as ALA
                int_elec_ala += AAG_[res][1]  # Elect_ala
                int_vdw_ala += AAG_[res][3]  # vdw_ala
                int_solv_ala += AAG_[res][5]  # solv_ala
            else:
                int_elec_ala += AAG_[res][0]  # Elect
                int_vdw_ala += AAG_[res][2]  # vdw
                int_solv_ala += AAG_[res][4]  # solv
        print("AlaScan;", chain_A[surf_res].get_resname(), int_elec_ala, ";", int_vdw_ala,
              ";", int_solv_ala+solv_E, ";", int_elec_ala+int_vdw_ala+int_solv_ala+solv_E, file=f)

    f.close()


if __name__ == "__main__":
    main(sys.argv[1:])