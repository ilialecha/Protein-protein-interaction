import  time
import  os
from    Energies            import *
from    VDWParameters       import *
from    SurfaceInteractions import *
from    Bio.PDB.PDBParser   import PDBParser
from    Bio.PDB.NACCESS     import NACCESS_atomic
from    Bio.PDB.PDBIO       import PDBIO

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
    
    try:  opts, args = getopt.getopt(argv,"hi:d:",["idir=","dist="])
    
    except getopt.GetoptError: man() ; sys.exit(2)
    
    for opt, arg in opts:
        if opt == '-h': man() ;sys.exit()
        elif opt in ("-i", "--idir"): INPUTDIR = arg
        elif opt in ("-d", "--dist"): DISTANCE = float(arg)

    if 0.0==DISTANCE: man() ; sys.exit(2)
    if "/" != INPUTDIR[-1]: INPUTDIR+="/"
    # Necessary for adding required parameters to atoms.
    NACCESS_BINARY =  INPUTDIR+"soft/NACCESS/naccess"

    print("(#) Adding Van der Waals parameters to PDB Structure.")
    time.sleep(1)
    
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


    print("(#) Identifying interaction surface residues.")
    time.sleep(1)
    
    chain_A         = st[0]["A"] # Obtaining chain A from model 0 --> ACE2
    chain_E         = st[0]["E"] # Obtaining chain E from model 0 --> Spike
    chain_A_at      = Selection.unfold_entities(chain_A,'A')
    chain_E_at      = Selection.unfold_entities(chain_E,'A')
    surfaceInt      = SurfaceInteractions(chain_A,chain_E,DISTANCE)
    surfaceRes      = surfaceInt.get_neighbors()
    surfaceAtm      = surfaceInt.get_atoms_in_contact(surfaceRes)
    surface_chain_A = []
    surface_chain_E = []

    '''
     This structure correlates the hypothetical interaction between key and item residues.
     Afterwards, instead of requiring Theta(n*m) where n and m are the lengths of the chains,
     the number of comparations is narrowed down to Theta(k*i) where k is a residue and i the
     number of neighbors of k.
    '''
    for k,v in surfaceRes.items():
        surface_chain_A.append(int(k));[surface_chain_E.append(res) for res in v]


    print("(#) Computing Electrostatic and Van der Waals interaction energies.")
    int_elect_A       = 0.
    int_elect_ala_A   = 0.
    int_vdw_A         = 0.
    int_vdw_ala_A     = 0.

    # For all surface residues in chain A we compute elec. and vdw energies
    # against all atoms in the structure. 

    for res in chain_A:
        tmp = Energies.calc_int_energies(st, res)
        
        int_elect_A       += tmp[0]
        int_elect_ala_A   += tmp[1]
        int_vdw_A         += tmp[2]
        int_vdw_ala_A     += tmp[3]

    
    print(f"\tElectrostatic {int_elect_A}")
    print(f"\tElectrostatic Ala {int_elect_ala_A}")
    print(f"\tVDW {int_vdw_A}")
    print(f"\tVDW Ala {int_vdw_ala_A}")


    print("(#) Computing solvation energies.")
    solv_AE = sum([Energies.calc_solvation(res) for res in Selection.unfold_entities(st[0], 'R')])
    solv_A  = sum([Energies.calc_solvation(res) for res in Selection.unfold_entities(chain_A, 'R')])
    solv_E  = sum([Energies.calc_solvation(res) for res in Selection.unfold_entities(chain_E, 'R')])
    print(f"\tSolv AE {solv_AE}")
    print(f"\tSolv A  {solv_A} ")
    print(f"\tSolv E  {solv_E} ")


if __name__ == "__main__":
   main(sys.argv[1:])
