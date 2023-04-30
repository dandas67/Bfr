# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 15:09:26 2022

@author: danist
"""

### "arrangements7-5.txt"
import chimera
import os

from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
rc("open #0 C:/Users/danist/Desktop/not_drive/Bfr_246_chimera_bruno-structures/postprocess246.mrc")
rc("open #1 C:/Users/danist/Desktop/not_drive/Bfr_246_chimera_bruno-structures/Bfr1_189-coot-RAZ.pdb")
rc("open #2 C:/Users/danist/Desktop/not_drive/Bfr_246_chimera_bruno-structures/Bfr2_189-coot-all2_heme.pdb")

with open("arrangements7-5.txt") as f:
    lines = f.readlines()
Bfr1_lines = lines[1:144:3]
Bfr2_lines = lines[2:145:3]

z=1
y=33
for line_Bfr1,line_Bfr2  in zip(Bfr1_lines,Bfr2_lines):
    sel_Bfr1_chains= str(line_Bfr1)
    sel_Bfr2_chains= str(line_Bfr2)
    rc(sel_Bfr1_chains)
    rc("write selected #1 C:/Users/danist/Desktop/not_drive/Bfr_246_chimera_bruno-structures/sel_Bfr1.pdb")
    rc("~select")
    rc(sel_Bfr2_chains)
    rc("write selected #2 C:/Users/danist/Desktop/not_drive/Bfr_246_chimera_bruno-structures/sel_Bfr2.pdb")
    rc("~select")
    rc("open #31 C:/Users/danist/Desktop/not_drive/Bfr_246_chimera_bruno-structures/sel_Bfr1.pdb")
    rc("open #32 C:/Users/danist/Desktop/not_drive/Bfr_246_chimera_bruno-structures/sel_Bfr2.pdb")
    basename= "Structure7-5"
    new_name = basename + "_" + str(z) 
    model= "#" + str(y)
    combine = "combine #31-#32 modelId " + model + " name " + new_name
    rc(combine)
    write = "write " + model + " C:/Users/danist/Desktop/not_drive/Bfr_246_chimera_bruno-structures/" + new_name + ".pdb"
    rc(write)
    rc("close #31")
    rc("close #32")
    close= "close " + model
    rc(close)    
    rc("~select")
    z=z+1
    y=y+1

z=1
y=33
for z in range(1:49):
    model= "#" + str(y)
    #structure_id = "#" + str(z)
    new_name = basename + "_" + str(z) 
    open_structure = "open " + model + " C:/Users/danist/Desktop/not_drive/Bfr_246_chimera_bruno-structures/" + new_name + ".pdb"
    rc(open_structure)
    color_blue = "color medium blue ,a ,r " + model + ":.a:.b:.c:.d:.e:.f:.g:.h:.i:.j:.k:.l:.m:.n:.o:.p:.q:.r:.s:.t:.u:.v:.w:.x"
    color_orange = "color #ffff6db60000 ,a ,r " + model + ":.A:.B:.C:.D:.E:.F:.G:.H:.I:.J:.K:.L:.M:.N:.O:.P:.Q:.R:.S:.T:.U:.V:.W:.X"
    rc(color_blue)
    rc(color_orange)
    y=y+1
