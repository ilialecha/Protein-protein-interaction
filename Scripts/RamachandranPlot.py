#----------------------------------------------------------------------------
# Created By    : Marçal Vázquez, Marc Romera and Ilia Lecha.
# Contributions : https://pypi.org/project/RamachanDraw/
# Created Date  : 04/11/2022
# version       = '1.0'
# --------------------------------------------------------------------------- 
import math
from Bio.PDB import *
from RamachanDraw import fetch, phi_psi, plot
from sys import argv

def degrees(rad_angle) :
    """Converts any angle in radians to degrees.

    If the input is None, the it returns None.
    For numerical input, the output is mapped to [-180,180]
    """
    if rad_angle is None :
        return None
    angle = rad_angle * 180 / math.pi
    while angle > 180 :
        angle = angle - 360
    while angle < -180 :
        angle = angle + 360
    return angle

def ramachandran_type(residue, next_residue) :
    """Expects Bio.PDB residues, returns ramachandran 'type'

    If this is the last residue in a polypeptide, use None
    for next_residue.

    Return value is a string: "General", "Glycine", "Proline"
    or "Pre-Pro".
    """
    if residue.resname.upper()=="GLY" :
        return "Glycine"
    elif residue.resname.upper()=="PRO" :
        return "Proline"
    elif next_residue is not None \
    and next_residue.resname.upper()=="PRO" :
        #exlcudes those that are Pro or Gly
        return "Pre-Pro"
    else :
        return "General"    

pdb_code = "6m0j"

try:
    structure = PDBParser().get_structure("6m0j","6m0j.pdb")

    plot(fetch(pdb_code))

    phi_psi_dict, ignored_res = phi_psi(fetch(pdb_code), return_ignored=True)

    print ("About to save angles to file...")
    output_file = open("%s_biopython.tsv" % pdb_code,"w") # 
    output_file.write("ID\tPHI\tPSI\TYPE\n")
    for model in structure :
        for chain in model:
            polypeptides = CaPPBuilder().build_peptides(chain)
            for poly_index, poly in enumerate(polypeptides) :
                phi_psi = poly.get_phi_psi_list()
                for res_index, residue in enumerate(poly) :
                    phi, psi = phi_psi[res_index]
                    if phi and psi :
                        #Don't write output when missing an angle
                        output_file.write("%s:Chain%s:%s%i\t%f\t%f\t%s\n" \
                            % (pdb_code, str(chain.id), residue.resname,
                            residue.id[1], degrees(phi), degrees(psi),
                            ramachandran_type(residue, poly[res_index+1])))
    output_file.close()
except Exception as e:
    print(e)


