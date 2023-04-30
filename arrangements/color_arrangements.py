"""
Created on Mon Aug 22 13:18:20 2022

@author: danist
"""
import chimera
import os
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand

#model= "#" + str(4)
#structure_id = "#" + str(z)
#new_name = basename + "_" + str(z) 
#open_structure = "open " + model + " C:/Users/danist/Desktop/not_drive/Bfr_246_chimera_bruno-structures/" + new_name + ".pdb"
#rc(open_structure)
color_blue = "color medium blue ,a ,r " + model + ":.a:.b:.c:.d:.e:.f:.g:.h:.i:.j:.k:.l:.m:.n:.o:.p:.q:.r:.s:.t:.u:.v:.w:.x"
color_orange = "color #ffff6db60000 ,a ,r " + model + ":.A:.B:.C:.D:.E:.F:.G:.H:.I:.J:.K:.L:.M:.N:.O:.P:.Q:.R:.S:.T:.U:.V:.W:.X"
rc(color_blue)
rc(color_orange)
