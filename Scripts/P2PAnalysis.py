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

    if distance > 0. : print("(#) Computing AAG of the interaction surface.")
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

            #Only if we are analyzing surface interactions we want to write to the output file.
            if distance > 0.:
                print(
                    en.residue_id(res),
                        elec[res], vdw[res], solvAB[res], solvA[res],
                        elec_ala[res], vdw_ala[res], solvAB_ala[res], solvA_ala[res],
                        sep="\t",file=f
                    )
            
    if distance > 0. : print("(#) Summary of interaction surface energies: ")
    else: print("(#) Summary of whole structure energies: ")

    print('{:20}: {:11.4f}'.format('Total Elec Int.', totalIntElec))
    print('{:20}: {:11.4f}'.format('Total Vdw Int.', totalIntVdw))
    print('{:20}: {:11.4f}'.format('Total Solv AB', totalSolv))
    print('{:19}{}: {:11.4f}'.format('Total Solv ', chids[0], totalSolvMon[chids[0]]))
    print('{:19}{}: {:11.4f}'.format('Total Solv ', chids[1], totalSolvMon[chids[1]]))
    print('{:20}: {:11.4f}'.format('DGintAB-A-B', total))

    # Only if we are analyzing surface interactions we want to write to the output file.
    if distance > 0.: 
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
    if distance > 0. : f.close() ; f2.close()

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
  # COMPUTING WHOLE STRUCTURE ENERGIES OF TWO PROTEINS.
    interaction_energies(args, st, st_chains,-1)
  #--------------------------------------------------------------------------------  

  #--------------------------------------------------------------------------------
  # COMPUTING INTERACTION SURFACE ENERGIES OF TWO PROTEINS GIVEN A CUTOFF DISTANCE
  # AND RUNNING AN ALA-SCANNING IN ORDER TO DETERMINE THE RELATIVE IMPORTANCE OF
  # EACH RESIDUE LOCATED ON THE INTERACTION SURFACE WHEN CHANGED BY AN ALANINE.
    interaction_energies(args, st, st_chains,args.cutoff_dist)
  # YOU SHALL FIND 2 FILES IN THE SAME DIRECTORY THAN THE GIVEN LIBRARY AND VDW 
  # PARAMETERS FILES CONTAINING THE RESULTS OF THE ANALYSYS FOR NEXT STEPS.   
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
'''
[20,21,23,24,25,26,27,28,29,30,31,32,33,34,35,37,38,39,41,42,43,45,48,49,75,76,78,79,80,81,82,83,84,324,325,326,327,329,330,351,352,353,354,355,356,357,386,387,388,393]

select (resid 20 and chain A ) , (resid 21 and chain A ) , (resid 23 and chain A ) , (resid 24 and chain A ) , (resid 25 and chain A ) , (resid 26 and chain A ) , (resid 27 and chain A ) , (resid 28 and chain A ) , (resid 29 and chain A ) , (resid 30 and chain A ) , (resid 31 and chain A ) , (resid 32 and chain A ) , (resid 33 and chain A ) , (resid 34 and chain A ) 
select (resid 35 and chain A ) , (resid 37 and chain A ) , (resid 38 and chain A ) , (resid 39 and chain A ) , (resid 41 and chain A ) , (resid 42 and chain A ) , (resid 43 and chain A ) , (resid 45 and chain A ) , (resid 48 and chain A ) , (resid 49 and chain A ) , (resid 75 and chain A ) , (resid 76 and chain A ) , (resid 78 and chain A ) , (resid 79 and chain A )
select (resid 80 and chain A ) , (resid 81 and chain A ) , (resid 82 and chain A ) , (resid 83 and chain A ) , (resid 84 and chain A ) , (resid 324 and chain A ) , (resid 325 and chain A ) , (resid 326 and chain A ) , (resid 327 and chain A ) , (resid 329 and chain A ) , (resid 330 and chain A ) , (resid 351 and chain A ) , (resid 352 and chain A ) , (resid 353 and chain A ) 
select (resid 354 and chain A ) , (resid 355 and chain A ) , (resid 356 and chain A ) , (resid 357 and chain A ) , (resid 386 and chain A ) , (resid 387 and chain A ) , (resid 388 and chain A ) , (resid 393 and chain A )
select (resid 20 and chain A ) or (resid 21 and chain A ) or (resid 23 and chain A ) or (resid 24 and chain A ) or (resid 25 and chain A ) or (resid 26 and chain A ) or (resid 27 and chain A ) or (resid 28 and chain A ) or (resid 29 and chain A ) or (resid 30 and chain A ) or (resid 31 and chain A ) or (resid 32 and chain A ) or (resid 33 and chain A ) or (resid 34 and chain A ) 
or (resid 35 and chain A ) or (resid 37 and chain A ) or (resid 38 and chain A ) or (resid 39 and chain A ) or (resid 41 and chain A ) or (resid 42 and chain A ) or (resid 43 and chain A ) or (resid 45 and chain A ) or (resid 48 and chain A ) or (resid 49 and chain A ) or (resid 75 and chain A ) or (resid 76 and chain A ) or (resid 78 and chain A ) or (resid 79 and chain A ) or 
(resid 80 and chain A ) or (resid 81 and chain A ) or (resid 82 and chain A ) or (resid 83 and chain A ) or (resid 84 and chain A ) or (resid 324 and chain A ) or (resid 325 and chain A ) or (resid 326 and chain A ) or (resid 327 and chain A ) or (resid 329 and chain A ) or (resid 330 and chain A ) or (resid 351 and chain A ) or (resid 352 and chain A ) or (resid 353 and chain A ) 
or (resid 354 and chain A ) or (resid 355 and chain A ) or (resid 356 and chain A ) or (resid 357 and chain A ) or (resid 386 and chain A ) or (resid 387 and chain A ) or (resid 388 and chain A ) or (resid 393 and chain A ) or (resid 403 and chain E ) or (resid 405 and chain E ) or (resid 406 and chain E ) or (resid 417 and chain E ) or (resid 421 and chain E ) or (resid 439 and chain E ) 
or (resid 444 and chain E ) or (resid 445 and chain E ) or (resid 446 and chain E ) or (resid 447 and chain E ) or (resid 448 and chain E ) or (resid 449 and chain E ) or (resid 453 and chain E ) or (resid 454 and chain E ) or (resid 455 and chain E )
403,405,406,417,421,439,444,445,446,447,448,449,453,454,455
for res in [403,405,406,417,421,439,444,445,446,447,448,449,453,454,455]:
    print(f"(resid {res} and chain E ) or ",end="")
'''

