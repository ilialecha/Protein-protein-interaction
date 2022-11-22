#! /usr/bin/python3
# ----------------------------------------------------------------------------
# Created By    : Ilia Lecha, Marc Romera and Marçal Vazquez.
# Contributions : Jose Luis Gelpi, Irene Acero, Alberto Meseguer
# Date          : 04/11/2022
# version       : '1.0'
# ---------------------------------------------------------------------------

"""
    Initial setup a structure for Energy evaluation
"""
import argparse
import os
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NACCESS import NACCESS_atomic
from Bio.PDB.PDBIO import PDBIO, Select
from VDWParameters import *
from SurfaceInteractions import * 
import Energies as en
from tqdm import tqdm


def interaction_energies(args, st, st_chains, distance):
    
    # Interface residues
    if distance > 0.:
        print(f"(#) Computing interaction surface given a cutoff of {distance} angstroms.")
        interface = SurfaceInteractions.get_interface(st, args.cutoff_dist)
        
        f = open(args.reslib_file[0:args.reslib_file.rfind("/")+1]+"interaction_energies.tsv", "w")
        print("RES_ID","ELECT.", "VDW.", "SOLV_AB", "SOLV_A","ELECT.ALA", "VDW.ALA", "SOLV_AB.ALA", "SOLV_A.ALA",sep="\t",file=f)

        f2 = open(args.reslib_file[0:args.reslib_file.rfind("/")+1]+"ala_interaction_energies.tsv", "w")
        print("RES.ID","CHAIN.ID","RES.Num","E+E.Ala","V+V.Ala","sAB+sAB.Ala","sA+sA.Ala","Total.AAG.Change",sep="\t",file=f2)

    elif distance == -1:
        interface = SurfaceInteractions.get_interface(st, 100000)
        print("(#) Computing energy of the whole structure.")
    

    # Initiatlize Energy aggregates
    elec = {} ; elec_ala = {} 
    vdw = {} ; vdw_ala = {}
    solvAB = {} ; solvAB_ala = {}
    solvA = {} ; solvA_ala = {}
    totalIntElec = 0. ;totalIntVdw = 0.
    totalSolv = 0.  ;totalSolvMon = {}

    # We get the chain ids,not always they are A and B
    print(f"(#) Obtaining chain identifiers. ")
    chids = []
    for ch in tqdm(st[0]):
        chids.append(ch.id)
        totalSolvMon[ch.id] = 0

    total = 0.

    print("(#) Computing AAG of the interaction surface.")
    for ch in st[0]:
        for res in tqdm(ch.get_residues()):
            # Checking if residue is located on the interaction surface.
            if args.cutoff_dist > 0 and res not in interface[ch.id]:
                continue

            # Obtaining electrostatic and vdw energies for each residue. Also computing the supposed case
            # that the residue is an ALANINE.
            elec[res], elec_ala[res], vdw[res], vdw_ala[res] = en.calc_int_energies(
                st[0], res)
            # Obtaining the total solvation energy.
            solvAB[res], solvAB_ala[res] = en.calc_solvation(st[0], res)
            # Obtaining the solvation of chain A.
            solvA[res], solvA_ala[res] = en.calc_solvation(
                st_chains[ch.id],
                st_chains[ch.id][0][ch.id][res.id[1]]
            )
            totalIntElec        += elec[res]
            totalIntVdw         += vdw[res]
            totalSolv           += solvAB[res]
            totalSolvMon[ch.id] += solvA[res]
            total               += elec[res] + vdw[res] + solvAB[res] - solvA[res]

            if distance > 0.:
                print(
                    en.residue_id(res),
                        elec[res], vdw[res], solvAB[res], solvA[res],
                        elec_ala[res], vdw_ala[res], solvAB_ala[res], solvA_ala[res],
                        sep="\t",file=f
                    )
            
    print("(#) Summary of interaction surface energies: ")
    print('{:20}: {:11.4f}'.format('Total Elec Int.', totalIntElec))
    print('{:20}: {:11.4f}'.format('Total Vdw Int.', totalIntVdw))
    print('{:20}: {:11.4f}'.format('Total Solv AB', totalSolv))
    print('{:19}{}: {:11.4f}'.format('Total Solv ', chids[0], totalSolvMon[chids[0]]))
    print('{:19}{}: {:11.4f}'.format('Total Solv ', chids[1], totalSolvMon[chids[1]]))
    print('{:20}: {:11.4f}'.format('DGintAB-A-B', total))


    print(                   
        "Total   ",
        totalIntElec, totalIntVdw, totalSolv, totalSolvMon[chids[0]],
        totalSolvMon[chids[1]], total, "-", "-",sep="\t",file=f
    )

    #------------------------------------------------------------------------------------

    if distance > 0: 
        print("(#) Ala Scanning: DDGs for X->Ala mutations on interface residues")

        for ch in st[0]:
            for res in tqdm(ch.get_residues()):
                if args.cutoff_dist > 0 and res not in interface[ch.id]:
                    continue
                tmp =en.residue_id(res).split()
                
                print(  tmp[0],
                        tmp[1][0:1],
                        tmp[1][1:],
                        - elec[res] + elec_ala[res],
                        - vdw[res] + vdw_ala[res],
                        - solvAB[res] + solvAB_ala[res],
                        - solvA[res] + solvA_ala[res],
                        - elec[res] + elec_ala[res] - vdw[res] + vdw_ala[res] - solvAB[res] +
                        solvAB_ala[res] - solvA[res] + solvA_ala[res],
                        sep="\t",file=f2)
    f.close()
    f2.close()

def main():

    class SelectChain(Select):
        def __init__(self, chid):
            self.id = chid
        def accept_chain(self, chain):
            if chain.id == self.id:
                return 1
            else:
                return 0

    # Setting up required parameters. 
    parse_cmd = argparse.ArgumentParser(
        prog='structure_setup',
        description='basic structure setup'
    )
    # Residue Library
    parse_cmd.add_argument(
        '--rlib',
        action='store',
        dest='reslib_file',
        default='data/aaLib.lib',
        help='Residue Library'
    )
    # Van der Waals parameters.
    parse_cmd.add_argument(
        '--vdw',
        action='store',
        dest='vdwprm_file',
        default='data/vdwprm',
        help='Vdw parameters'
    )
    # Distance cutoff.
    parse_cmd.add_argument(
        '--dist',
        action='store',
        dest='cutoff_dist',
        default=8.0,
        type=float,
        help='Cutoff distance for determining the interface (0: use all residues):'
    )
    # Pdb file.
    parse_cmd.add_argument('pdb_file', help='Input PDB', type=open)

    args = parse_cmd.parse_args()


    # Setting up location of naccess executable.
    #NACCESS_BINARY = '/home/ilia/Escritorio/Bioinformatics/1_Second_year/BioPhysics/Seminars/Protein-protein-interaction/Scripts/files/soft/NACCESS/naccess'
    NACCESS_BINARY = args.reslib_file[0:args.reslib_file.rindex("/")+1] + "soft/NACCESS/naccess"


    # Loading Libraries
    # loading residue library from files/aaLib.lib
    residue_library = ResiduesDataLib(args.reslib_file)

    # loading VdW parameters
    ff_params = VdwParamset(args.vdwprm_file)

    parser = PDBParser(PERMISSIVE=1)

    # load structure from PDB file of PDB ifle handler
    st = parser.get_structure('STR', args.pdb_file.name)

    # assign data types, and charges from libraries
    # We will use the xtra attribute in Bio.PDB.Atom to hold the new data
    # Possible errors on N-term and C-Term atoms
    # Possible errors on HIS alternative forms
    print("(#) Adding vdw parameters.")
    add_atom_parameters(st, residue_library, ff_params)

    # Calculating surfaces
    # The specific PATH to naccess script (in soft) is needed
    # ASA goes to .xtra field directly

    srf = NACCESS_atomic(st[0], naccess_binary=NACCESS_BINARY)

    # Prepare surfaces for the separate chains
    # Alternatively the twp PDB files can be prepared outside and parsed here

    io = PDBIO()
    st_chains = {}
    # Using BioIO trick (see tutorial) to select chains
    for ch in st[0]:
        io.set_structure(st)
        io.save('tmp.pdb', SelectChain(ch.id))
        st_chains[ch.id] = parser.get_structure('stA', 'tmp.pdb')
        add_atom_parameters(st_chains[ch.id], residue_library, ff_params)
        srfA = NACCESS_atomic(
            st_chains[ch.id][0], naccess_binary=NACCESS_BINARY)
    os.remove('tmp.pdb')
    

  #--------------------------------------------------------------------------------
  # COMPUTING INTERACTION SURFACE ENERGIES OF TWO PROTEINS GIVEN A CUTOFF DISTANCE
  # AND RUNNING AN ALA-SCANNING IN ORDER TO DETERMINE THE RELATIVE IMPORTANCE OF
  # EACH RESIDUE LOCATED ON THE INTERACTION SURFACE WHEN CHANGED BY AN ALANINE.
    interaction_energies(args, st, st_chains,args.cutoff_dist)
  #--------------------------------------------------------------------------------

if __name__ == "__main__":
    main()

'''
/bin/python3 /home/ilia/Escritorio/Bioinformatics/1_Second_year/BioPhysics/Seminars/Protein-protein-interaction/Scripts/int_analysis.py 
--rlib /home/ilia/Escritorio/Bioinformatics/1_Second_year/BioPhysics/Seminars/Protein-protein-interaction/Scripts/files/aaLib.lib 
--vdw /home/ilia/Escritorio/Bioinformatics/1_Second_year/BioPhysics/Seminars/Protein-protein-interaction/Scripts/files/vdwprm 
--dist 4.0 
/home/ilia/Escritorio/Bioinformatics/1_Second_year/BioPhysics/Seminars/Protein-protein-interaction/Scripts/files/6m0j_fixed.pdb
'''
