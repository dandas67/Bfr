# -*- coding: utf-8 -*-
"""
Created on Sun Jun 18 13:42:23 2023

@author: danda
"""
#seqs with 1st position changed from M to the aa of the vector(same as pdb)
Bfr1= "PRGSPKVISVLNGLLTGELTAADQYFVHARMLENWGFKVLYERIEHERHDELDHAGLLINRILFLEGVPDVASRAALNIGSDVPKMMANDLAYELQVVDELKAAIALCESERDYDTRRILVHLLEETEQDHVRWLEVQVGLIDKLGLKNYLQSAAGEIA"
Bfr2= "DKANRTVLAALNDVLRHQLTAINQYFLHARMMKNWGFNALGKHEYKESIEEMKAADKLIERILLLEGLPNLQDLGKLLIGENVPEMLKNDFAMEKDAHADLVKTIALCEKQADYVSRDLLSEFLEECEERMDFYETQLELVKKMGEQNYLQSAVGALED"

diff_positions_list= [] #list for the seq positions where Bfr1 and Bfr2 are different

for i in range(len(Bfr1)-1):
    if i+1 <= 129: # to match index with aa position
        Bfr1_i = i # up to aa 129, same index
    else:
        Bfr1_i = i+1 # for 130 and upwards, we want to incrase A index by 1      
    if Bfr1[Bfr1_i] != Bfr2[i]:
        diff_positions_list.append(i+1) #correction of the indexing for the pdb.
 
 
#change the pdb based on the position from the former part:
#import pandas as pd
from biopandas.pdb import PandasPdb

input_path= "C:/Users/danist/OneDrive - post.bgu.ac.il/Ph.D/Bfr/relion/post/284/rerefined/"
file_name= "284_Bfr1_real_space_refined_005.pdb"
input_pdb= input_path+file_name

read_PDB1 = PandasPdb().read_pdb(input_pdb)
Bfr1_all_ATOM_df = read_PDB1.df['ATOM']#all atoms from bfr1
  
#the chain I want to edit
chain="a"
#a data frame of only the chain I want to edit
chain_df=Bfr1_all_ATOM_df[Bfr1_all_ATOM_df.chain_id ==chain]
#iterate only on the lines of a specific chain
for index, line in chain_df.iterrows():
    if line.residue_number in diff_positions_list: #if the position is the the list
        #print(line.residue_number)
        print(index)
        Bfr1_all_ATOM_df.at[index, "residue_name"] = "UNK"
        #Bfr1_all_ATOM_df.iloc[index, chain_df.columns.get_loc("residue_name")]= "UNK" #dosnt work. #change the residue name to A or UNK
        #delete all extra atoms that are not in alanine:
        if line.atom_name != 'N' and line.atom_name != 'CA' and line.atom_name != 'C' and line.atom_name != 'O' and line.atom_name != 'CB':
            Bfr1_all_ATOM_df=Bfr1_all_ATOM_df.drop(index) 
     
#update the "atom" df.            
read_PDB1.df['ATOM'] = Bfr1_all_ATOM_df 
read_PDB1.to_pdb(path='C:/Users/danist/OneDrive - post.bgu.ac.il/Ph.D/Bfr/relion/post/284/relaxed-structure/284_Bfr1_005_polyALA.pdb', #save pdb
            records=None, 
            gz=False, 
            append_newline=True)

#print the residues and positions of the differences.
for i in diff_positions_list:
    print(Bfr1[i-1], i)
