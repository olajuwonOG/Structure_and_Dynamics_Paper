#!/usr/bin/env python

import sys
sys.path.append('/blue/jsampath/tejuosho.taofeek/PhD_Project_Polymer_Dispersity/AIChE_2023/Schulz_Zimm_Replica_2/D_1.4/140') #/Users/brownj/Research/Scripts
import numpy as np
import pppmd
import fileinput

skip = 0
file = 'e2e_Dataskip_NEW%d.csv' % skip

r, ir, timestep, boxbds, id2type, id2mol, mol2ids = pppmd.read_lammpstrj('production.dump', skip_beginning=0, skip_between=0)
file
e2e_mols = []
for mol in range(len(mol2ids)):
	if len(mol2ids[mol]) > 1:
		e2e_mols.append(mol2ids[mol])

e2e_autocorr = pppmd.end2end_autocorr(r, ir, boxbds, e2e_mols)
OUT = open(file, 'w')
OUT.write("timestep, E2E\n")

print("t Reet_dot_Ree0")
for t in range(len(timestep)):
    OUT.write("%7i, %8.4f\n" % ((timestep[t]-timestep[0]), e2e_autocorr[t]))
	#print(timestep[t]-timestep[0], e2e_autocorr[t])
OUT.close()
