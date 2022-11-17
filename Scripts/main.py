#----------------------------------------------------------------------------
# Created By    : Ilia Lecha, Marc Romera and Marçal Vazquez.
# Contributions : Jose Luis Gelpi, Irene Acero, Alberto Meseguer
# Date          : 04/11/2022
# version       = '1.0'
# --------------------------------------------------------------------------- 
import  os
import  getopt
from    Energies            import *
from    VDWParameters       import *
from    SurfaceInteractions import *
from    Bio.PDB.PDBParser   import PDBParser
from    Bio.PDB.NACCESS     import NACCESS_atomic
from    Bio.PDB.PDBIO       import PDBIO
from    tqdm                import tqdm

def man():
    print ('main.py -i <inputdir> -d <distance>')
    print ('Inputdir must contain the following files:')
    print ('aaLib.lib:\tAtoms library')
    print ('Pdb_file:\tContaining a model with two chains.')
    print ('vdwprm:\t\tParameters for Van der Waals forces.')

def main(argv):
    # Setting up of the required parameters. 
    INPUTDIR   = ''
    DISTANCE    = 0.0
    
    try:
        opts, args = getopt.getopt(argv,"hi:d:",["idir=","dist="])
    except getopt.GetoptError:
        man() ; sys.exit(2)
    
    for opt, arg in opts:
        if opt == '-h':
            man() ;sys.exit()
        elif opt in ("-i", "--idir"):
            INPUTDIR = arg
        elif opt in ("-d", "--dist"):
            DISTANCE = float(arg)

    if 0.0==DISTANCE:
        man() ; sys.exit(2)
    if "/" != INPUTDIR[-1]:
        INPUTDIR+="/"

    # Necessary for adding required parameters to atoms.
    NACCESS_BINARY =  INPUTDIR+"soft/NACCESS/naccess"

    print("(#) Adding Van der Waals parameters to PDB Structure.")
    
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
        srfA = NACCESS_atomic(st_chains[ch.id][0], naccess_binary=NACCESS_BINARY)
    os.remove('tmp.pdb')


    f = open(INPUTDIR+"energies.tsv", "w")
    print("Type;Res;Electrostatic_AE;vdw_AE;solv_AE;total",file=f)

    print("(#) Identifying interaction surface residues.")
    
    chain_A         = st[0]["A"] # Obtaining chain A from model 0 --> ACE2
    chain_E         = st[0]["E"] # Obtaining chain E from model 0 --> Spike
    surfaceInt      = SurfaceInteractions(chain_A,chain_E,DISTANCE)
    surfaceRes      = surfaceInt.get_neighbors()
    surface_chain_A = set()
    surface_chain_E = set()

    for k,v in surfaceRes.items():
        surface_chain_A.add(int(k))
        [surface_chain_E.add(res) for res in v]

    print("(#) Computing Electrostatic and Van der Waals interaction energies.")
    int_elect_A = 0.
    int_vdw_A   = 0.
    solv_A      = 0.

    # For all surface residues in chain A we compute elec. and vdw energies
    # against all atoms in the structure. 

    AAG = set()

    for res in tqdm(chain_A):
        tmp = Energies.calc_int_energies(st, res)
        int_elect_A       += tmp[0]
        int_vdw_A         += tmp[1]
        solv_A            += tmp[2]
        # Adding (elec,vdw) to AAG =  {(elec,vdw),(elec,vdw),(elec,vdw),...}
        AAG.add(tmp)

    print(f"\tElectrostatic {int_elect_A}")
    print(f"\tVDW {int_vdw_A}")

    print("(#) Computing solvation energies.")

    # Solv_AE = Solv_A + Solv_E so we can omit some computation.
    solv_AE = sum([Energies.calc_solvation(res) for res in Selection.unfold_entities(st[0], 'R')])
    solv_E = sum([Energies.calc_solvation(res) for res in Selection.unfold_entities(st[0]['E'], 'R')])
    
    print(f"\tSolv AE {solv_AE}")
    print(f"{solv_A + solv_E}")

    print(f"(#) Writting into energies.tsv")
    print("Normal;","All",";",int_elect_A,";",int_vdw_A,";",solv_AE,";", int_vdw_A+int_elect_A+solv_AE,file=f)

    print("(#) Running Ala-Scanning:")
    print(f"(#) For each execution, one of these residues({surface_chain_A}) shall be treated as an Alanine.")

    for surf_A in surface_chain_A:
        
        int_elect_ala_A   = 0.
        int_vdw_ala_A     = 0.
        int_solv_ala      = 0.
        solv_AE_ala       = 0.
        for res in tqdm(chain_A):
            if res.get_id()[1] == surf_A:
                # Found residue to be treated as an Alanine.
                tmp = Energies.calc_int_energies_ala(st, res)
                int_elect_ala_A   += tmp[0]
                int_vdw_ala_A     += tmp[1]
            else:
                tmp = Energies.calc_int_energies(st, res)
                int_elect_ala_A   += tmp[0]
                int_vdw_ala_A     += tmp[1]

        for res in Selection.unfold_entities(st[0],"R"):
            if res.get_id()[1] == surf_A:
                # Found residue to be treated as an Alanine.
                solv_AE_ala += Energies.calc_solvation_ala(res)
            else:
                solv_AE_ala += Energies.calc_solvation(res)

        print(f"\tModified res {chain_A[surf_A].get_resname()}")
        print(f"\tElectrostatic {int_elect_ala_A}")
        print(f"\tVDW {int_vdw_ala_A}")
        print(f"\tSolvation {solv_AE_ala}")
        print(f"(#) Writting into energies.tsv")
        print("AlaScan;",chain_A[surf_A].get_resname(),int_elect_ala_A,";",int_vdw_ala_A,";",solv_AE_ala,";", int_vdw_ala_A+int_elect_ala_A+solv_AE_ala,file=f)
    f.close()

    '''
    Results before Ala-Scanning.
    Electrostatic -2.3990267908564156
    VDW -86.33619619342456
    Solv AE -517.432398000001
    '''
    '''
    Results after Ala-Scanning.
    Electrostatic -3.1085559070555813
    VDW -80.4134868382032
    Solvation -517.1682650000009
    '''

if __name__ == "__main__":
   main(sys.argv[1:])
