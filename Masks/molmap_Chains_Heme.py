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
    save_chain = "volume #42 save C:/Users/danist/Desktop/not_drive/Hotspots_284/refined_Hotspots_matlab/rerefined_outputs/" + file_name_chains
    rc(save_chain)
    rc("close #42")

#C:/Users/danist/Desktop/not_drive/Hotspots_284/refined_Hotspots_matlab/rerefined_outputs

chains = list(B,O,N,H,C,V,G,K,M,X,F,R,S,I,D,U,T,J,L,W,E,P,Q,A)
test= "blabla"
heme_ids= list(11,11,101,101,81,81,16,16,61,61,31,31,1,1,36,36,21,21,41,41,51,51,56,56)
for chain,heme_id in zip(chains,heme_ids):
    command = "sel #2:." + chain + ":" + heme_id + "._"
    print(command)
    

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