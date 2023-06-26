# -*- coding: utf-8 -*-
"""
Created on Sun Jun 25 12:19:25 2023

@author: danist
"""

from biopandas.pdb import PandasPdb
from Bio import SeqIO

def compare_sequences(input_pdb, fasta_file):
    # Read PDB file
    ppdb = PandasPdb()
    ppdb.read_pdb(input_pdb)
    # Find the first occurrence of the specified chain ID
    chain_df = ppdb.df['ATOM']
    chain_df = chain_df[chain_df['chain_id'] == chain_id].reset_index(drop=True)

    # Exclude non-standard residues (e.g., HETATM)
    standard_residues = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
        'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
        'TYR': 'Y', 'VAL': 'V'
    }

    chain_df = chain_df[chain_df['residue_name'].isin(standard_residues.keys())]

    # Get unique residue sequence numbers
    residue_sequence_numbers = chain_df['residue_number'].unique()

    # Sort the residue sequence numbers
    sorted_residue_sequence_numbers = sorted(residue_sequence_numbers)

    # Extract amino acid sequence in one-letter code
    aa_sequence_pdb = ''
    for residue_number in sorted_residue_sequence_numbers:
        residue_df = chain_df[chain_df['residue_number'] == residue_number]
        residue_name = residue_df['residue_name'].values[0]
        amino_acid_code = standard_residues[residue_name]
        aa_sequence_pdb += amino_acid_code
    
    
    # Read FASTA file
    fasta_sequences = SeqIO.parse(fasta_file, 'fasta')
    for fasta in fasta_sequences:
        # Extract sequence from FASTA
        seq_id = fasta.id
        seq = str(fasta.seq)

        # Compare sequences
        if seq == aa_sequence_pdb:
            print(f"The sequence from PDB chain {chain_id} matches the sequence in {fasta_file} ({seq_id}).")
        else:
            # Find differing positions
            differing_positions = [i + 1 for i, (pdb_aa, fasta_aa) in enumerate(zip(aa_sequence_pdb, seq)) if pdb_aa != fasta_aa]
            print(f"The sequence from PDB chain {chain_id} and the sequence in {fasta_file} ({seq_id}) differ at positions:")
            print(differing_positions)

    print("Sequence comparison completed.")
    

# Extract amino acid sequence from PDB chain "a"
chain_id = "a"
    
   


# Provide the paths to your PDB and FASTA files
input_path= "C:/Users/danist/OneDrive - post.bgu.ac.il/Ph.D/Bfr/relion/post/284/rerefined/"
file_name= "284_Bfr1_real_space_refined_005.pdb"
input_pdb= input_path+file_name

fasta_file= "C:/Users/danist/OneDrive - post.bgu.ac.il/Ph.D/Bfr/Bfr1_seq.fasta"

compare_sequences(input_pdb, fasta_file)