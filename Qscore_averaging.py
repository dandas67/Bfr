"""
Created on Tue Jun 20 15:26:40 2023

@author: danist
"""

import pandas as pd
#from biopandas.pdb import PandasPdb
excel_file = 'C:/Users/danist/OneDrive - post.bgu.ac.il/Ph.D/Bfr/relion/post/284/Qscore/Bfr2_scores.xlsx'
xl_df = pd.read_excel(excel_file)
xl_out= pd.DataFrame()



# Group the dataframe by 'residue_number' and 'chain_id' and calculate the average Qscore
averaged_df = xl_df.groupby(['residue_number', 'chain_id'])['Qscore'].mean().reset_index()

# Print the averaged dataframe
print(averaged_df)
# Save the averaged dataframe to an Excel file
averaged_df.to_excel('C:/Users/danist/OneDrive - post.bgu.ac.il/Ph.D/Bfr/relion/post/284/Qscore/Bfr2_averaged_scores.xlsx', index=False)

# Print a message to confirm the file has been saved
print("Averaged scores saved to 'averaged_scores.xlsx'")
