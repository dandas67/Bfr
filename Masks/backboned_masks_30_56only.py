# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 11:57:52 2022

@author: danist
"""
import chimera
import os

from chimera import runCommand as rc # use 'rc' as shorthand for runCommand


##For Arg30 of Bfr1:
Bfr1_hotspots_list = (41,64,134,43,118,123,26,54,98,45,56,74,72,30)
for hot_res1 in Bfr1_hotspots_list:
    sel_hotspots1 = "sel #1: " + str(hot_res1)
    rc(sel_hotspots1)
    if hot_res1 ==30:
        rc("~sel @NH2")
    rc("molmap sel 5 onGrid #0 modelId #42")
    ##threshold?
    #file_name_Bfr1 = "Bfr1_molmap_" + str(hot_res1) +".mrc"
    file_name_Bfr1 = "Bfr1_molmap_" + str(hot_res1) +"bb_res5.mrc"
    save_Bfr1 = "volume #42 save C:/Users/danist/Desktop/not_drive/Hotspots_284/refined_Hotspots_matlab/backboned_masks/" + file_name_Bfr1
    rc(save_Bfr1)
    rc("close #42")
##For Asp56 of Bfr2:
Bfr2_hotspots_list = (56)
for hot_res2 in Bfr2_hotspots_list:
    sel_hotspots2 = "sel #2: " + str(hot_res2)
    rc(sel_hotspots2)
    rc("~sel #2:26._")
    rc("~sel @OD1")
    rc("molmap sel 5 onGrid #0 modelId #42")
    ##threshold?
    #file_name_Bfr2 = "Bfr2_molmap_" + str(hot_res2) +".mrc"
    file_name_Bfr2 = "Bfr2_molmap_" + str(hot_res2) +"bb_res5.mrc"
    save_Bfr2 = "volume #42 save C:/Users/danist/Desktop/not_drive/Hotspots_284/refined_Hotspots_matlab/backboned_masks/" + file_name_Bfr2
    rc(save_Bfr2)
    rc("close #42")
    
    
##For Arg30 of Bfr1:
Bfr1_hotspots_list = (30)
for hot_res1 in Bfr1_hotspots_list:
    sel_hotspots1 = "sel #1: " + str(hot_res1)
    rc(sel_hotspots1)
    rc("~sel @NH2")
    rc("molmap sel 5 onGrid #0 modelId #42")
    ##threshold?
    #file_name_Bfr1 = "Bfr1_molmap_" + str(hot_res1) +".mrc"
    file_name_Bfr1 = "Bfr1_molmap_" + str(hot_res1) +"bb_res5.mrc"
    save_Bfr1 = "volume #42 save C:/Users/danist/Desktop/not_drive/Hotspots_284/refined_Hotspots_matlab/backboned_masks/" + file_name_Bfr1
    rc(save_Bfr1)
    rc("close #42")
##For Asp56 of Bfr2:
Bfr2_hotspots_list = (56)
for hot_res2 in Bfr2_hotspots_list:
    sel_hotspots2 = "sel #2: " + str(hot_res2)
    rc(sel_hotspots2)
    rc("~sel #2:26._")
    rc("~sel @OD1")
    rc("molmap sel 5 onGrid #0 modelId #42")
    ##threshold?
    #file_name_Bfr2 = "Bfr2_molmap_" + str(hot_res2) +".mrc"
    file_name_Bfr2 = "Bfr2_molmap_" + str(hot_res2) +"bb_res5.mrc"
    save_Bfr2 = "volume #42 save C:/Users/danist/Desktop/not_drive/Hotspots_284/refined_Hotspots_matlab/backboned_masks/" + file_name_Bfr2
    rc(save_Bfr2)
    rc("close #42")