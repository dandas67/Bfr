
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 12:53:25 2022

@author: danist
"""


import chimera
import os

from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
#rc("open #0 C:/Users/danist/Desktop/not_drive/Hotspots_284/284postprocess.mrc")
#rc("open #1 C:/Users/danist/Desktop/not_drive/Hotspots_284/Bfr1_284_raz.pdb")
#rc("open #2 C:/Users/danist/Desktop/not_drive/Hotspots_284/Bfr2_284_raz_heme201.pdb")

#C:/Users/danist/Desktop/not_drive/Hotspots_284/refined_Hotspots_matlab/rerefined_outputs
#rc("sel #1-2")
#rc("molmap sel 3 onGrid #0 modelId #50")
#rc("volume #50 save C:/Users/danist/Desktop/not_drive/Hotspots_284/Hotspots_24aug/Bfr1_and_Bfr1-molmap_res3")
#sel_hotspots = "sel #1:54:26:134:30:122:41:45:98 #2:201:54:26:133:30:122:41:45:98"
Bfr1_hotspots_list = (41,64,134,43,118,123,26,54,98,45,56,74,72,30,91,102,93,135,112)
Bfr2_hotspots_list = (41,64,133,43,118,123,26,54,98,45,56,74,72,30,91,102,93,134,112)

#Bfr1_hotspots_list = (41,64,134,43,118,123,54,98,45,56,74,72,30,91,102,93,135)
#Bfr2_hotspots_list = (41,64,133,43,118,123,54,98,45,56,74,72,30,91,102,93,134)

#Bfr1_hotspots_list = (41,56)
#Bfr2_hotspots_list = (41,56)

for hot_res1 in Bfr1_hotspots_list:
    sel_hotspots1 = "sel #1: " + str(hot_res1)
    rc(sel_hotspots1)
    #rc("~sel @CA,N,O,C")
    if hot_res1 ==30:
        rc("~sel @NH2")
    rc("molmap sel 5 onGrid #0 modelId #42")
    ##threshold?
    #file_name_Bfr1 = "Bfr1_molmap_" + str(hot_res1) +".mrc"
    file_name_Bfr1 = "Bfr1_molmap_" + str(hot_res1) +"bb_res5.mrc"
#    save_Bfr1 = "volume #42 save C:/Users/danist/Desktop/not_drive/Hotspots_284/phenix_refined_Hotspots/" + file_name_Bfr1
    save_Bfr1 = "volume #42 save C:/Users/danist/Desktop/not_drive/hotspots_292/" + file_name_Bfr1
    rc(save_Bfr1)
    rc("close #42")
    
for hot_res2 in Bfr2_hotspots_list:
    sel_hotspots2 = "sel #2: " + str(hot_res2)
    rc(sel_hotspots2)
    rc("~sel #2:26._")
    rc("~sel #2:41._")
    rc("~sel #2:56._")    
    #rc("~sel @CA,N,O,C")
    if hot_res2 ==56:
        rc("~sel @OD1")
    rc("molmap sel 5 onGrid #0 modelId #42")
    ##threshold?
    #file_name_Bfr2 = "Bfr2_molmap_" + str(hot_res2) +".mrc"
    file_name_Bfr2 = "Bfr2_molmap_" + str(hot_res2) +"bb_res5.mrc"
#    save_Bfr2 = "volume #42 save C:/Users/danist/Desktop/not_drive/Hotspots_284/phenix_refined_Hotspots/" + file_name_Bfr2
    save_Bfr2 = "volume #42 save C:/Users/danist/Desktop/not_drive/hotspots_292/" + file_name_Bfr2
    rc(save_Bfr2)
    rc("close #42")

