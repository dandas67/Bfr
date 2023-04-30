import chimera
import os
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand


location= "C:/Users/danist/OneDrive - post.bgu.ac.il/Ph.D/Bfr/relion/post/284/arrangemnts/"

sel_Bfr1_chains = "sel #1:.f:.r:.x:.m:.c:.v:.t:.j:.p:.e:.o:.b:.i:.s"
sel_Bfr2_chains = "sel #2:.K:.G:.Q:.A:.W:.L:.U:.D:.N:.H"
rc(sel_Bfr1_chains)
write_sel_Bfr1 = "write selected #1 " + location + "sel_Bfr1.pdb"
rc(write_sel_Bfr1)
rc("~select")
rc(sel_Bfr2_chains)
write_sel_Bfr2 = "write selected #2 " + location + "sel_Bfr2.pdb"
rc(write_sel_Bfr2)
rc("~select")
open_Bfr1 = "open #31 " + location + "sel_Bfr1.pdb"
rc(open_Bfr1)
open_Bfr2 = "open #32 " + location + "sel_Bfr2.pdb"
rc(open_Bfr2)

z=1
y=33
basename = "Structure7-5"
new_name = basename + "_" + str(z) 
model= "#" + str(y)

combine = "combine #31-#32 modelId " + model + " name " + new_name
rc(combine)
write = "write " + model +" " + location + new_name + ".pdb"
rc(write)
rc("close #31")
rc("close #32")
rc("~select")
#close = "close " + model
#rc(close)    
