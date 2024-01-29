# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 10:46:07 2023

@author: danist
"""
from biopandas.pdb import PandasPdb
input_path= "C:/Users/danist/OneDrive - post.bgu.ac.il/Ph.D/Bfr/relion/post/284/relaxed-structure/finals_with_Raz/"
file_name_Bfr1= "Bfr12-UNK-lwrf_relax284_real_space_refined_031.pdb"


input_pdb1= input_path+file_name_Bfr1

read_PDB1 = PandasPdb().read_pdb(input_pdb1)
PDB1_all_ATOM_df = read_PDB1.df['ATOM']#all atoms from bfr1
PDB1_all_ATOM_df.occupancy = 1 #change occupancy of all atoms in PDB1


# Update the 'ATOM' DataFrame in PDB_mix_read
read_PDB1.df['ATOM'] = PDB1_all_ATOM_df

output_file_name= 'Bfr12_rsr031_occ1.pdb'
output_pdb= input_path+output_file_name

read_PDB1.to_pdb(path=output_pdb, #save pdb
            records=None, 
            gz=False, 
            append_newline=True)