# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 14:47:14 2023

@author: danist
"""

from biopandas.pdb import PandasPdb

def convert_to_polyglycine(input_pdb_file, output_pdb_file):
    # Load PDB file using BioPandas
    ppdb = PandasPdb().read_pdb(input_pdb_file)

    # Change all residues to alanine
    ppdb.df['ATOM']['residue_name'] = 'GLY'


 # Filter only the N, O, CA, and CB atoms
    valid_atoms = ['N','C', 'O', 'CA']
    ppdb.df['ATOM'] = ppdb.df['ATOM'][ppdb.df['ATOM']['atom_name'].isin(valid_atoms)]

    # Update atom serial numbers
    ppdb.df['ATOM']['atom_number'] = range(1, len(ppdb.df['ATOM']) + 1)


    # Update insertion codes
    ppdb.df['ATOM']['insertion'] = ''
    

    # Save the new structure to a PDB file using BioPandas
    ppdb.to_pdb(output_pdb_file, records=['ATOM'])
    
input_path= "C:/Users/danist/OneDrive - post.bgu.ac.il/Ph.D/Bfr/relion/post/284/relaxed-structure/" # Change this to your input path
file_name= "Bfr12-UNK-lwrf_relax284.pdb" # Change this to your input PDB file name
input_pdb1= input_path+file_name

output_file_name= 'Bfr12_polygly_test.pdb' # Change this to the desired output PDB file name
output_pdb= input_path+output_file_name


if __name__ == "__main__":
    input_pdb_file = input_pdb1  
    output_pdb_file = output_pdb 

    convert_to_polyglycine(input_pdb_file, output_pdb_file)

