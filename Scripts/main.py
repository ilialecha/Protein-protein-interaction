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


    f = open(INPUTDIR+"energies.tsv", "w"); print("Type;Res;Electrostatic_AE;vdw_AE;solv_AE;total",file=f)

    print("(#) Identifying interaction surface residues.")
    
    chain_A         = st[0]["A"] # Obtaining chain A from model 0 --> ACE2
    chain_E         = st[0]["E"] # Obtaining chain E from model 0 --> Spike
    surfaceInt      = SurfaceInteractions(chain_A,chain_E,DISTANCE)
    surfaceRes      = surfaceInt.get_neighbors()
    surface_chain_A = []

    for k,v in surfaceRes.items():
        surface_chain_A.append(int(k))
        #[surface_chain_E.add(res) for res in v]

    print("(#) Computing all energies: Electrostatic, Van der waals and Solvation.")

    # Solv_AE = Solv_A + Solv_E so we can omit some computation.
    solv_E = sum([Energies.calc_solvation2(st,res)[0] for res in Selection.unfold_entities(st[0]['E'], 'R')])

    AAG_ = dict()
    # Computing energies for all atoms against all atoms. More information related in Energies.py
    for res in tqdm(chain_A):
        tmp  = Energies.calc_int_energies2(st, res)
        tmp2 = Energies.calc_solvation2(st,res) 
        #                 (elec,elec_ala,vdw,vdw_ala,solv,solv_ala)
        AAG_[res.id[1]] = (tmp[0],tmp[1],tmp[2],tmp[3],tmp2[0],tmp2[1])

    int_elec = 0.
    int_elec_ala = 0.
    int_vdw = 0.
    int_vdw_ala = 0.
    int_solv = 0.
    int_solv_ala = 0.

    for res in AAG_:
        int_elec        += AAG_[res][0]
        int_elec_ala    += AAG_[res][1]
        int_vdw         += AAG_[res][2]
        int_vdw_ala     += AAG_[res][3]
        int_solv        += AAG_[res][4]
        int_solv_ala    += AAG_[res][5]

    print(f"(#) Writting into output file.")
    print("Normal;","All",";",int_elec,";",int_vdw,";",int_solv+solv_E,";", int_elec+int_vdw+int_solv+solv_E,file=f)

    print("(#) Running Ala-Scanning:")
    print(f"(#) For each execution, one of the following residues: \n(#) {surface_chain_A} is treated as an Alanine.")


    for surf_res in tqdm(surface_chain_A):
        int_elec_ala = 0.
        int_vdw_ala = 0.
        int_solv_ala = 0.
        for res in AAG_.keys():
            if surf_res == res: #Treated as ALA
                int_elec_ala += AAG_[res][1] # Elect_ala
                int_vdw_ala  += AAG_[res][3] # vdw_ala
                int_solv_ala += AAG_[res][5] # solv_ala    
            else:
                int_elec_ala += AAG_[res][0] # Elect
                int_vdw_ala  += AAG_[res][2] # vdw
                int_solv_ala += AAG_[res][4] # solv
        print("AlaScan;",chain_A[surf_res].get_resname(),int_elec_ala,";",int_vdw_ala,";",int_solv_ala+solv_E,";", int_elec_ala+int_vdw_ala+int_solv_ala+solv_E,file=f)

    f.close()
if __name__ == "__main__":
   main(sys.argv[1:])


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

'''
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
            int_solv_ala      += Energies.calc_solvation_ala(res)
        else:
            tmp = Energies.calc_int_energies(st, res)
            int_elect_ala_A   += tmp[0]
            int_vdw_ala_A     += tmp[1]
            int_solv_ala      += Energies.calc_solvation(res)'''


'''for res in Selection.unfold_entities(st[0],"R"):
    if res.get_id()[1] == surf_A:
        # Found residue to be treated as an Alanine.
        solv_AE_ala += Energies.calc_solvation_ala(res)
    else:
        solv_AE_ala += Energies.calc_solvation(res)'''

'''print(f"\tModified res {chain_A[surf_A].get_resname()}")
print(f"\tElectrostatic {int_elect_ala_A}")
print(f"\tVDW {int_vdw_ala_A}")
print(f"\tSolvation {int_solv_ala+solv_E}")
print(f"(#) Writting into energies.tsv")
print("AlaScan;",chain_A[surf_A].get_resname(),int_elect_ala_A,";",int_vdw_ala_A,";",int_solv_ala,";", int_vdw_ala_A+int_elect_ala_A+int_solv_ala,file=f)'''