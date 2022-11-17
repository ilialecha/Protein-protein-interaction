#----------------------------------------------------------------------------
# Created By    : Ilia Lecha, Marc Romera and Marçal Vazquez.
# Contributions : Jose Luis Gelpi, Irene Acero and Alberto Meseguer
# Date          : 04/11/2022
# version       = '1.0'
# --------------------------------------------------------------------------- 
import biobb_structure_checking
import biobb_structure_checking.constants as cts
from biobb_structure_checking.structure_checking import StructureChecking

class PdbChecker():
    def check(PATH):
        base_dir_path=biobb_structure_checking.__path__[0]
        args = cts.set_defaults(base_dir_path,{'notebook':True})
        with open(args['commands_help_path']) as help_file: print(help_file.read())
        
        args['input_structure_path'] = PATH + '6m0j.pdb'
        args['output_format'] = 'pdb'
        args['output_structure_path'] = PATH + '6m0j_fixed.pdb'
        args['output_structure_path_charges'] = PATH + '6m0j_fixed.pdbqt'
        args['keep_canonical'] = False
        args['debug'] = False
        args['verbose'] = False

        #Initializing checking engine, loading structure and showing statistics
        st_c = StructureChecking(base_dir_path, args)

        # Show models
        st_c.models()

        # Show chains
        st_c.chains()

        # Residues  with alternative locations
        st_c.altloc()

        # Choosing an alternative forms for each residue.
        st_c.altloc('occupancy')

        # Checking if there are no more Alt. Loc.
        st_c.altloc()

        # Checking for the presence of metals.
        st_c.metals()

        # Checking for ligands
        st_c.ligands()

        # Deleting ligands.
        st_c.ligands('All')

        # Removing possible hydrogens
        st_c.rem_hydrogen("yes")

        # Removing possible Water molecules.
        st_c.water("yes")

        # Removing all possible Amides.
        st_c.amide('all')

        # Looking for res. with incorrect side-chain chirality found
        st_c.chiral()

        # Looking for and fixing several problems with the backbone.
        st_c.backbone()
        st_c.backbone('--fix_atoms All --fix_chain none --add_caps none')

        # Detecting and fixing possible missing protein side chains.
        st_c.fixside()

        # Detecting possible -S-S bonds 
        st_c.getss()
        st_c.getss('all')

        # Adding hydrogens
        # st_c.add_hydrogen()
        # st_c.add_hydrogen('auto')

        # Detecting steric clashes base on distance
        st_c.clashes()

        # Complete check process in a single method
        # st_c.checkall()

        #Saving structure
        st_c._save_structure(args['output_structure_path'])
