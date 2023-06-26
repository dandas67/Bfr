import chimera
import os

from chimera import runCommand as rc # use 'rc' as shorthand for runCommand

##open map, Bfr1, Bfr2, and Bfr1polyA-UNK  in that order
#write the path to the output folder:
#location= "C:/Users/danist/Desktop/not_drive/relax_arrangement/"
location= "C:/Users/danist/OneDrive - post.bgu.ac.il/Ph.D/Bfr/relion/post/284/relaxed-structure/"

#select chains of Bfr1:
rc("sel #1:.c:.v:.m:.x:.s:.i:.t:.j:.e:.p:.b:.o")
write_sel_Bfr1 = "write selected #1 " + location + "sel_Bfr1.pdb"
rc(write_sel_Bfr1)

#select chains of Bfr2 Including the hemes:
rc("sel #2:.A:.Q:.D:.U:.G:.K:.N:.H:101._:16._:36._:56._")
write_sel_Bfr2 = "write selected #2 " + location + "sel_Bfr2.pdb"
rc(write_sel_Bfr2)


#select chains of Bfr1polyA-UNK:
rc("sel #3:.l:.w:.f:.r")
write_sel_Bfr1_UNK = "write selected #3 " + location + "sel_Bfr1-UNK.pdb"
rc(write_sel_Bfr1_UNK)

#open the 3 files with the chains
open_Bfr1 = "open #31 " + location + "sel_Bfr1.pdb"
rc(open_Bfr1)
open_Bfr2 = "open #32 " + location + "sel_Bfr2.pdb"
rc(open_Bfr2)
open_Bfr1_UNK = "open #33 " + location + "sel_Bfr1-UNK.pdb"
rc(open_Bfr1_UNK)

#combine all the chains into one structure
rc("combine #31-#33 modelId #34 name  Bfr12-UNK-lwrf_relax284.pdb")

#save the combined structure to a pdb file
save_combined_structure = "write #34 " + location + "Bfr12-UNK-lwrf_relax284_test.pdb"
rc(save_combined_structure)

#color by chains
rc("color medium blue ,a ,r #34:.a:.b:.c:.d:.e:.g:.h:.i:.j:.k:.m:.n:.o:.p:.q:.s:.t:.u:.v:.x")
rc("color #ffff6db60000 ,a ,r #34:.A:.B:.C:.D:.E:.F:.G:.H:.I:.J:.K:.L:.M:.N:.O:.P:.Q:.R:.S:.T:.U:.V:.W:.X")
