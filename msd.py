#!/usr/bin/python

# Script:  msd.py
# Purpose: calculate mean sq displacement of various types
#g1(t) - for middle monomers in selected selected_chains
#g3(t) - center of mass msd for selected chains
# Syntax:  msdions.py < filename
# Example: msdions.py < test.dump (dump file with scaled coordinates)
# Author:  Lisa Hall from Mark's g(r) code
# Modified by: Taofeek Tejuosho

# derived from fortran code
# -------------------------------------------------------

import sys,string
from numpy import *
from math import *
import fileinput
import numpy
import warnings
# INPUT PARAMETERs
nconf =  2015  #number of configurations to take data for
npoly =    2000 #need this for com calc
warnings.filterwarnings("ignore", category=RuntimeWarning)

def msdcalc():
    global msdvec1, msdvec2, msdvec3, msdvec4, msdvec5, msdvec6, msdvec7, timestep, msd, xx, yy, zz, x0, y0, z0, x0com, y0com, z0com, xcom, ycom, zcom, step, num1, num2, num3, num4, num5, num6, npoly, selected_mol_ids
    global natoms, typea,msdcall, num_counts, atom_id_mol_id_dict, match_list, middle_monomer_count, num_counts, msdvec8, overall_monomer_middle, msdvec9
    msd[1]=0
    msd[2]=0
    msd[3]=0
    msd[4]=0
    msd[5]=0
    msd[6]=0
    msd[7]=0
    msd[8]=0
    msd[9]=0
    #overall monomer MSD - middle; g1(t)
    overall_monomer_middle=0
    for ii in range(1,natoms+1): #msd loops for monomer/bead
        chain_length1 = num_counts[atom_id_mol_id_dict[ii]]
        middle_index1 = chain_length1//2
        if ii % chain_length1 == middle_index1:
            overall_monomer_middle+=1
            xr = xx[ii] - x0[ii]
            yr = yy[ii] - y0[ii]
            zr = zz[ii] - z0[ii]
            r2 = xr*xr + yr*yr + zr*zr
        #msd[typea[ii]] = msd[typea[ii]]+r2
            msd[8] = msd[8]+r2
    #specific monomer MSD - middle; g1(t) for specific chains
    middle_monomer_count = 0
    for ii in range(1, natoms + 1):
        if ii in match_list:  # Check if the atom belongs to a selected chain
            chain_length = num_counts[atom_id_mol_id_dict[ii]]
            middle_index = chain_length // 2
            if ii % chain_length == middle_index:
                middle_monomer_count +=1
                xr = xx[ii] - x0[ii]
                yr = yy[ii] - y0[ii]
                zr = zz[ii] - z0[ii]
                r2 = xr * xr + yr * yr + zr * zr
                msd[typea[ii]] = msd[typea[ii]] + r2

    #center of mass MSD - Overall
    for ii in range(1,npoly+1): #msd loop for chains to calculate their center of mass
        xr = xcom[ii] - x0com[ii]
        yr = ycom[ii] - y0com[ii]
        zr = zcom[ii] - z0com[ii]
        r2 = xr*xr + yr*yr + zr*zr
        #print(ii, "for", r2)
        msd[9] = msd[9]+r2
    #center of mass MSD - Specific chains
    for ii in range(1,npoly+1): #msd loop for chains to calculate their center of mass
        if ii in selected_mol_ids:
            xr = xcom[ii] - x0com[ii]
            yr = ycom[ii] - y0com[ii]
            zr = zcom[ii] - z0com[ii]
            r2 = xr*xr + yr*yr + zr*zr
            #print(ii, "for", r2)
            msd[7] = msd[7]+r2
    msd[1]=msd[1]/(middle_monomer_count)
    msd[2]=msd[2]/num2
    msd[3]=msd[3]/num3
    msd[4]=msd[4]/num4
    msd[5]=msd[5]/num5
    msd[6]=msd[6]/num6
    msd[7]=msd[7]/len(selected_mol_ids)
    msd[8]=msd[8]/overall_monomer_middle
    msd[9]=msd[9]/npoly
    msdvec1[msdcall]=msdvec1[msdcall]+msd[1]
    msdvec2[msdcall]=msdvec2[msdcall]+msd[2]
    msdvec3[msdcall]=msdvec3[msdcall]+msd[3]
    msdvec4[msdcall]=msdvec4[msdcall]+msd[4]
    msdvec5[msdcall]=msdvec5[msdcall]+msd[5]
    msdvec6[msdcall]=msdvec6[msdcall]+msd[6]
    msdvec7[msdcall]=msdvec7[msdcall]+msd[7]
    msdvec8[msdcall]=msdvec8[msdcall]+msd[8]
    msdvec9[msdcall]=msdvec9[msdcall]+msd[9]
    timestep[msdcall]=step

    #OUT = open(file, 'a')
   # OUT.write("%7i %8.4f %8.4f %8.4f\n" % (step,msd[1],msd[2],msd[3]))
   # for ii in range(1,natoms+1):
    #    OUT.write("%7i %7i %7i %7i %7i %8.4f %8.4f %8.4f\n" % (mol[ii],typea[ii],xi[ii],yi[ii],zi[ii],xc[ii],yc[ii],zc[ii]))
   # OUT.close()
# end def msdcalc
warnings.resetwarnings()
def wholecalculation():
    global msdvec1, msdvec2, msdvec3, msdvec4, msdvec5, msdvec6,msdvec7, timestep, msd, xx, yy, zz, x0, y0, z0, x0com, y0com, z0com, xcom, ycom, zcom, step, num1, num2, num3, num4, num5, num6, npoly
    global skip
    global natoms, typea,msdcall, atom_id_mol_id_dict, selected_mol_ids, match_list, num_counts, middle_monomer_count, msdvec8, msdvec9

    infiles = ['production.dump'] #PSsorted.lammpstrj
    file = 'msdcm_R1skip%dg1_N=30.csv'% skip
    dummy=0

    msdcall=0
    msdvec1=zeros(nconf,float)
    msdvec2=zeros(nconf,float)
    msdvec3=zeros(nconf,float)
    msdvec4=zeros(nconf,float)
    msdvec5=zeros(nconf,float)
    msdvec6=zeros(nconf,float)
    msdvec7=zeros(nconf,float)
    msdvec8=zeros(nconf,float)
    msdvec9 = zeros(nconf,float)
    timestep=zeros(nconf,int)

    IN = fileinput.input(infiles)
    natoms = 0
    num_counts = {}


    for loopnum in range(0,1): #read starting timestep

        IN.readline()
        line = IN.readline()      # time step
        fields = line.split()
        step = int(fields[0])
        print(("read step %i" % step))
        IN.readline()
        line = IN.readline()      # number of atoms
        fields = line.split()
        natoms = int(fields[0])  #reads number of atoms from first configuration; don't support changing num atoms
        dim=natoms+1
        xc=zeros(dim,float32)
        yc=zeros(dim,float32)
        zc=zeros(dim,float32)
        xx=zeros(dim,float32)
        yy=zeros(dim,float32)
        zz=zeros(dim,float32)
        x0=zeros(dim,float32)
        y0=zeros(dim,float32)
        z0=zeros(dim,float32)
        x0com=zeros(npoly+1,float32)
        y0com=zeros(npoly+1,float32)
        z0com=zeros(npoly+1,float32)
        nbead=zeros(npoly+1)
        xi=zeros(dim)
        yi=zeros(dim)
        zi=zeros(dim)
        typea=[0]*dim
        mol=[0]*dim
        msd=zeros(10,float32)
        IN.readline()
        line = IN.readline()
        [xm,xp] = list(map(float,line.split()))
        line = IN.readline()
        [ym,yp] = list(map(float,line.split()))
        line = IN.readline()
        [zm,zp] = list(map(float,line.split()))
        line = IN.readline()
        xbox = xp - xm
        ybox = yp - ym
        zbox = zp - zm
        vol = xbox*ybox*zbox
        xbox2 = xbox/2.0
        xcom=zeros(npoly+1,float32) #initialize these before every config read
        ycom=zeros(npoly+1,float32)
        zcom=zeros(npoly+1,float32)
        xcm=ycm=zcm=0.0
        num1 = 0
        num2 = 0
        num3 = 0
        num4 = 0
        num5 = 0
        num6 = 0
        selected_atoms = []
        atom_id_mol_id_dict = {}
        for j in range(1,dim):
            line = IN.readline()
            #todo: make a list called "dumpstyle" so this can be changed easily in one place?
            [ii,molj,typej,q,x1,x2,x3,n1,n2,n3] = line.split()
            k=int(ii)
            atom_id_mol_id_dict[k] = int(molj)
            #key_value = [k, int(molj)]
            #if key_value[1] in selected_mol_ids:
                #selected_atoms.append(k)
            chain_id = int(molj)
            if chain_id in num_counts:
                num_counts[chain_id] +=1
            elif chain_id not in num_counts:
                num_counts[chain_id] = 1
            typea[k] = int(typej)
            if typea[k] == 1:
                num1=num1+1
            elif typea[k] == 2:
                num2=num2+1
            elif typea[k] == 3:
                num3=num3+1
            elif typea[k] == 4:
                num4=num4+1
            elif typea[k] == 5:
                num5=num5+1
            elif typea[k] == 6:
                num6=num6+1
            mol[k] = int(molj)
            xc[k] = xbox*(float(x1)-0.5) #scaled coords go from 0 to 1; want =\0[mm     //to go from -xbox/2 to xbox/2
            yc[k] = ybox*(float(x2)-0.5)
            zc[k] = zbox*(float(x3)-0.5)
            xi[k] = int(round(float(n1)))
            yi[k] = int(round(float(n2)))
            zi[k] = int(round(float(n3)))
            xx[k] = xc[k] + xbox*xi[k]
            yy[k] = yc[k] + ybox*yi[k]
            zz[k] = zc[k] + zbox*zi[k]
        selected_mol_ids = []
        a = 25
        b = 35
        for key, value in num_counts.items():
            #if  value == 140:
            if a <= value <= b:
            #if  10 < value < 215:
                selected_mol_ids.append(key)
            #else:
            #    print("No chain with this length in the melt")
        print("list of IDs of chains with 60 beads is",selected_mol_ids)
        match_list = [key for key, value in atom_id_mol_id_dict.items() if value in selected_mol_ids]

        #print("num_counts", (num_counts))
        print("\n")
        print("\n")
        #print("num_counts", len(num_counts))
        print("\n")
        #print("total atoms and chains are", atom_id_mol_id_dict)
        print("\n")
        #print("length of the dictionary", len(atom_id_mol_id_dict))
        print("\n")
        #print("selected atoms IDs", match_list)
        print("\n")
        for j in range(1,dim): #position of the center of mass of the beads at time, t = 0
            xx[j] = xx[j]
            yy[j] = yy[j]
            zz[j] = zz[j]
            #xcm = xcm + xx[j]
            #ycm = ycm + yy[j]
            #zcm = zcm + zz[j]
        #xcm=xcm/natoms
        #ycm=ycm/natoms
        #zcm=zcm/natoms
        #print(("xcm is ", xcm))
        #for j in range(1,dim): #position of beads relative to their center of mass for t = 0
        #    xx[j] = xx[j] - xcm
        #    yy[j] = yy[j] - ycm
        #    zz[j] = zz[j] - zcm
        #selected_chains = []
        #xcom1 = np.zeros(len(selected_chains), float32)
        #ycom1 = np.zeros(len(selected_chains), float32)
        #zcom1 = np.zeros(len(selected_chains), float32)

        for j in range(1,dim): #save initial config
            x0[j] = xx[j]
            y0[j] = yy[j]
            z0[j] = zz[j]
            if mol[j] < npoly+1: #implies that polymers come first in molecule number; positions of the chains from their center of mass
                xcom[mol[j]] = xcom[mol[j]]+xx[j]
                ycom[mol[j]] = ycom[mol[j]]+yy[j]
                zcom[mol[j]] = zcom[mol[j]]+zz[j]
                nbead[mol[j]] = nbead[mol[j]]+1
        for m in range(1,npoly+1):
            #if m in selected_mol_ids:
            xcom[m] = xcom[m]/float(nbead[m])
            ycom[m] = ycom[m]/float(nbead[m])
            zcom[m] = zcom[m]/float(nbead[m])
        #save initial config
            x0com[m] = xcom[m]
            y0com[m] = ycom[m]
            z0com[m] = zcom[m]

        istep=step

        msdcalc()
        msdcall+=1

        print("Reading config file....")

        istart = loopnum + 1
        print("istart is ",istart)
        # Read configuration from zconfig
        for kconf in range(istart,nconf):
          try:
            IN.readline()
            line = IN.readline()      # time step
            fields = string.split(line)
            step = int(fields[0])
            IN.readline()
            line = IN.readline()      # number of atoms
            fields = string.split(line)
            num = fields[0]
            IN.readline()
            line = IN.readline()
            [xm,xp] = list(map(float,line.split()))
            line = IN.readline()
            [ym,yp] = list(map(float,line.split()))
            line = IN.readline()
            [zm,zp] = list(map(float,line.split()))
            line = IN.readline()
            xbox = xp - xm
            ybox = yp - ym
            zbox = zp - zm
            vol = xbox*ybox*zbox
            xbox2 = xbox/2.0
            xcom=zeros(npoly+1,float32) #initialize these before every config read
            ycom=zeros(npoly+1,float32)
            zcom=zeros(npoly+1,float32)
            xcm=ycm=zcm=0.0 #initialize center of mass of system
            nbead=zeros(npoly+1) #initialize center of mass each chain
            for j in range(1,dim):
                line = IN.readline()
                #todo: make a list called "dumpstyle" so this can be changed easily in one place?
                [ii,molj,typej,q,x1,x2,x3,n1,n2,n3] = string.split(line)
                k=int(ii)
                typea[k] = int(typej)
                mol[k] = int(molj)
                #print("coordinate x1 is ",x1)
                xc[k] = xbox*(float(x1)-0.5) #scaled coords go from 0 to 1; want to go from -xbox/2 to xbox/2
                yc[k] = ybox*(float(x2)-0.5)
                zc[k] = zbox*(float(x3)-0.5)
                xi[k] = int(round(float(n1)))
                yi[k] = int(round(float(n2)))
                zi[k] = int(round(float(n3)))
                xx[k] = xc[k] + xbox*xi[k]
                yy[k] = yc[k] + ybox*yi[k]
                zz[k] = zc[k] + zbox*zi[k]
            """for j in range(1,dim): #loop for calculating the center of mass system; at time, t > initial time
                xcm = xcm + xx[j]
                ycm = ycm + yy[j]
                zcm = zcm + zz[j]
            xcm=xcm/natoms
            ycm=ycm/natoms
            zcm=zcm/natoms"""
            for j in range(1,dim): #position of each bead relative to the center of mass at time, t > initial time ----Change here if you want to calc g2(t)
                #xx[j] = xx[j] - xcm #xx[j] = xx[j] - x0[j] for g2(t)
                #yy[j] = yy[j] - ycm
                #zz[j] = zz[j] - zcm
                xx[j] = xx[j]#xx[j] = xx[j] - x0[j] for g2(t)
                yy[j] = yy[j]
                zz[j] = zz[j]
            for j in range(1,dim):
                if mol[j] < npoly+1: #implies that polymers come first in molecule number
                    #rint("mol is",mol[j])
                    xcom[mol[j]] = xcom[mol[j]]+xx[j]
                    ycom[mol[j]] = ycom[mol[j]]+yy[j]
                    zcom[mol[j]] = zcom[mol[j]]+zz[j]
                    nbead[mol[j]] = nbead[mol[j]]+1
            for m in range(1,npoly+1):
                xcom[m] = xcom[m]/float(nbead[m])
                ycom[m] = ycom[m]/float(nbead[m])
                zcom[m] = zcom[m]/float(nbead[m])
            # calculate msd for this config
            msdcalc()
            msdcall+=1
            print(("%7i %7i %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n" % (istep,step,msd[1],msd[2],msd[3],msd[4],msd[5],msd[6],msd[7], msd[8], msd[9])))
          except:break
        fileinput.close()

    # end of loop
    OUT = open(file, 'w')
    #OUT.write("msd\n")
    OUT.write("step, msd_Monomer_Specific, msd_type2, msd_type3, msd_type4, msd_type5, msd_type6, msd_COM_specific, msdMonomerOverall, msdCOMoverall\n")
    for i in range(0,msdcall):
    ##    msdvec1[i]=msdvec1[i]/(nconf-i)
    ##    msdvec2[i]=msdvec2[i]/(nconf-i)
    ##    msdvec3[i]=msdvec3[i]/(nconf-i)
    ##    msdvec4[i]=msdvec4[i]/(nconf-i)
    ##    msdvec5[i]=msdvec5[i]/(nconf-i)
        OUT.write("%7i, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f\n" % (timestep[i],msdvec1[i],msdvec2[i],msdvec3[i],msdvec4[i],msdvec5[i],msdvec6[i],msdvec7[i], msdvec8[i], msdvec9[i]))
    OUT.close()

#main:
skip=0
wholecalculation()
