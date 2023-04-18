# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 16:29:52 2022

@author: danist
"""


import chimera
import os

from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
#rc("open #0 C:/Users/danist/Desktop/not_drive/Hotspots_284/284postprocess.mrc")
#rc("open #1 C:/Users/danist/Desktop/not_drive/Hotspots_284/Bfr1_284_raz.pdb")
#rc("open #2 C:/Users/danist/Desktop/not_drive/Hotspots_284/Bfr2_284_raz_heme201.pdb")


import string
#Bfr1_chains_list=(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x)
low_alphabet = list(string.ascii_lowercase)
Bfr1_chains_list = low_alphabet[0:24]

#Bfr2_chains_list=(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X)
up_alphabet = list(string.ascii_uppercase)
Bfr2_chains_list = up_alphabet[0:24]

for Bfr1_chain,Bfr2_chain in zip(Bfr1_chains_list,Bfr2_chains_list):
    select_chains = "sel #1:." + str(Bfr1_chain) + " #2:." + str(Bfr2_chain)
    rc(select_chains)
    rc("molmap sel 5 onGrid #0 modelId #42")
    file_name_chains = "Bfr_molmap_chain" + str(Bfr2_chain) +"_res5.mrc"
    save_chain = "volume #42 save C:/Users/danist/Desktop/not_drive/Hotspots_284/phenix_refined_Hotspots/" + file_name_chains
    rc(save_chain)
    rc("close #42")

#C:/Users/danist/Desktop/not_drive/Hotspots_284/refined_Hotspots_matlab/rerefined_outputs
