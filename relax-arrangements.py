# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 16:39:33 2022

@author: danist
"""

import chimera
import os

from chimera import runCommand as rc # use 'rc' as shorthand for runCommand

##open map, Bfr1 and Bfr2(after heme edit) in that order
# rc("open #0 C:/Users/danist/Desktop/not_drive/Bfr_246_chimera_bruno-structures/postprocess246.mrc")
# rc("open #1 C:/Users/danist/Desktop/not_drive/Bfr_246_chimera_bruno-structures/Bfr1_189-coot-RAZ.pdb")
# rc("open #2 C:/Users/danist/Desktop/not_drive/Bfr_246_chimera_bruno-structures/Bfr2_189-coot-RAZ.pdb")

# rc("sel #1:.c:.v:.m:.x:.s:.i:.t:.j:.e:.p:.b:.o")
# rc("write selected #1 C:/Users/danist/Desktop/not_drive/relax_arrangement/sel_Bfr1.pdb")
# rc("sel #2:.A:.Q:.D:.U:.G:.K:.N:.H")
# rc("write selected #2 C:/Users/danist/Desktop/not_drive/relax_arrangement/sel_Bfr2.pdb")
# rc("open #31 C:/Users/danist/Desktop/not_drive/relax_arrangement/sel_Bfr1.pdb")
# rc("open #32 C:/Users/danist/Desktop/not_drive/relax_arrangement/sel_Bfr2.pdb")
# rc("combine #31-#32 modelId #33 name  Bfr12_relax284.pdb")
# rc("write #33 C:/Users/danist/Desktop/not_drive/relax_arrangement/Bfr12_relax284.pdb")

# rc("close #31")
# rc("close #32")
# rc("color medium blue ,a ,r #33:.a:.b:.c:.d:.e:.f:.g:.h:.i:.j:.k:.l:.m:.n:.o:.p:.q:.r:.s:.t:.u:.v:.w:.x")
# rc("color #ffff6db60000 ,a ,r #33:.A:.B:.C:.D:.E:.F:.G:.H:.I:.J:.K:.L:.M:.N:.O:.P:.Q:.R:.S:.T:.U:.V:.W:.X")


rc("sel #1:.c:.v:.m:.x:.s:.i:.t:.j:.e:.p:.b:.o")
rc("write selected #1 C:/Users/danist/Desktop/not_drive/relax_arrangement/sel_Bfr1.pdb")
rc("sel #2:.A:.Q:.D:.U:.G:.K:.N:.H:101._:16._:36._:56._")
rc("write selected #2 C:/Users/danist/Desktop/not_drive/relax_arrangement/sel_Bfr2.pdb")
rc("sel #3:.l:.w:.f:.r")##polyA-UNK
rc("write selected #3 C:/Users/danist/Desktop/not_drive/relax_arrangement/sel_Bfr1-UNK.pdb")
rc("open #31 C:/Users/danist/Desktop/not_drive/relax_arrangement/sel_Bfr1.pdb")
rc("open #32 C:/Users/danist/Desktop/not_drive/relax_arrangement/sel_Bfr2.pdb")
rc("open #33 C:/Users/danist/Desktop/not_drive/relax_arrangement/sel_Bfr1-UNK.pdb")
rc("combine #31-#33 modelId #34 name  Bfr12-UNK-lwrf_relax284.pdb")
rc("write #34 C:/Users/danist/Desktop/not_drive/relax_arrangement/Bfr12-UNK-lwrf_relax284.pdb")

sel #2:.B:11._
sel #2:.O:11._

sel #2:.N:101._
sel #2:.H:101._

sel #2:.C:81._
sel #2:.V:81._

sel #2:.G:16._
sel #2:.K:16._

sel #2:.M:61._
sel #2:.X:61._

sel #2:.F:31._
sel #2:.R:31._

sel #2:.S:1._
sel #2:.I:1._

sel #2:.D:36._
sel #2:.U:36._

sel #2:.T:21._
sel #2:.J:21._

sel #2:.L:41._
sel #2:.W:41._

sel #2:.E:51._
sel #2:.P:51._

sel #2:.Q:56._
sel #2:.A:56._