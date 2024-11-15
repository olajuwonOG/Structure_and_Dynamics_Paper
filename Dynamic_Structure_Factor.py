# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.


#Author: Taofeek Tejuosho

import math
import numpy as np
import csv
import time

# Function to read timesteps from a file
def readTimesteps(inputFile, timestepsCount=None):
    timesteps = {}
    with open(inputFile) as f:
        line = f.readline()
        while line != '' and (not timestepsCount or len(timesteps) < timestepsCount):
            if line.split()[0] == "ITEM:":
                timestep = int(f.readline().split()[0])
                timesteps[timestep] = {}
                f.readline()
                atoms = int(f.readline().split()[0])
                timesteps[timestep]["atoms"] = atoms
                f.readline()
                x = f.readline()
                x_pos = float(x.split()[1])
                x_neg = float(x.split()[0])
                timesteps[timestep]["-x"] = x_neg
                timesteps[timestep]["x"] = x_pos
                y = f.readline()
                y_pos = float(y.split()[1])
                y_neg = float(y.split()[0])
                timesteps[timestep]["-y"] = y_neg
                timesteps[timestep]["y"] = y_pos
                z = f.readline()
                z_pos = float(z.split()[1])
                z_neg = float(z.split()[0])
                timesteps[timestep]["-z"] = z_neg
                timesteps[timestep]["z"] = z_pos
                headers = f.readline()
                timesteps[timestep]["headers"] = headers
                timesteps[timestep]["mols"] = {}
                for i in range(atoms):
                    line = f.readline()
                    id = int(line.split()[0])
                    mol = int(line.split()[1])
                    type = int(line.split()[2])
                    q = int(line.split()[3])
                    xs = float(line.split()[4]) * (timesteps[timestep]['x'] - timesteps[timestep]['-x'])
                    ys = float(line.split()[5]) * (timesteps[timestep]['y'] - timesteps[timestep]['-y'])
                    zs = float(line.split()[6]) * (timesteps[timestep]['z'] - timesteps[timestep]['-z'])
                    ix = int(line.split()[7])
                    iy = int(line.split()[8])
                    iz = int(line.split()[9])
                    xbox = (timesteps[timestep]['x'] - timesteps[timestep]['-x'])
                    ybox = (timesteps[timestep]['y'] - timesteps[timestep]['-y'])
                    zbox = (timesteps[timestep]['z'] - timesteps[timestep]['-z'])
                    if xs > xbox:
                        xs = xs - xbox
                    elif xs < 0:
                        xs = xs + xbox
                    if ys > ybox:
                        ys = ys - ybox
                    elif ys < 0:
                        ys = ys + ybox
                    if zs > zbox:
                        zs = zs - zbox
                    elif zs < 0:
                        zs = zs + zbox
                    xs = xbox * float(line.split()[4]) + ix * xbox
                    ys = ybox * float(line.split()[5]) + iy * ybox
                    zs = zbox * float(line.split()[6]) + iz * zbox
                    if mol in timesteps[timestep]["mols"]:
                        timesteps[timestep]["mols"][mol]["atoms"].append((id, type, q, xs, ys, zs, ix, iy, iz))
                    else:
                        timesteps[timestep]["mols"][mol] = {}
                        timesteps[timestep]["mols"][mol]["atoms"] = []
                        timesteps[timestep]["mols"][mol]["atoms"].append((id, type, q, xs, ys, zs, ix, iy, iz))
                for mol in timesteps[timestep]["mols"]:
                    x = 0
                    y = 0
                    z = 0
                    length = len(timesteps[timestep]["mols"][mol]["atoms"])
                    timesteps[timestep]["mols"][mol]["length"] = length
                    for i in range(length):
                        x += timesteps[timestep]["mols"][mol]["atoms"][i][3]
                        y += timesteps[timestep]["mols"][mol]["atoms"][i][4]
                        z += timesteps[timestep]["mols"][mol]["atoms"][i][5]
                    timesteps[timestep]["mols"][mol]["xcm"] = round(x / length, 4)
                    timesteps[timestep]["mols"][mol]["ycm"] = round(y / length, 4)
                    timesteps[timestep]["mols"][mol]["zcm"] = round(z / length, 4)
            line = f.readline()
    return timesteps

# Function to calculate Rg^2
def Rg2Calculator(timesteps):
    for timestep in timesteps:
        for mol in timesteps[timestep]["mols"]:
            length = timesteps[timestep]["mols"][mol]["length"]
            rg2 = 0
            xcm = timesteps[timestep]["mols"][mol]["xcm"]
            ycm = timesteps[timestep]["mols"][mol]["ycm"]
            zcm = timesteps[timestep]["mols"][mol]["zcm"]
            for i in range(length):
                x = timesteps[timestep]["mols"][mol]["atoms"][i][3] - xcm
                y = timesteps[timestep]["mols"][mol]["atoms"][i][4] - ycm
                z = timesteps[timestep]["mols"][mol]["atoms"][i][5] - zcm
                rg2 += (x**2 + y**2 + z**2)
            timesteps[timestep]["mols"][mol]["rg2"] = round(rg2 / length, 4)
    return timesteps

# Function to write Rg^2 to a file
def writeRg2ToFile(timesteps, outFile):
    with open(outFile, 'w') as f:
        f.write("rg2, molNumber, timestep\n")
        for timestep in timesteps:
            for mol in timesteps[timestep]["mols"]:
                rg2 = str(timesteps[timestep]["mols"][mol]["rg2"])
                molNumber = str(mol)
                timestepWrite = str(timestep)
                f.write(rg2 + ", " + molNumber + ", " + timestepWrite + "\n")

# Function to write average Rg^2 to a file
def writeRg2AverageToFile(timesteps, outFile, frames=1):
    with open(outFile, 'w') as f:
        f.write("rg2\n")
        mols = {}
        frameCount = 0
        for timestep in timesteps:
            if frameCount < frames:
                for mol in timesteps[timestep]["mols"]:
                    if mol not in mols:
                        mols[mol] = timesteps[timestep]["mols"][mol]["rg2"]
                    else:
                        mols[mol] += timesteps[timestep]["mols"][mol]["rg2"]
                frameCount += 1
            else:
                continue
        sum = 0
        for mol in mols:
            sum += mols[mol] / frames
        f.write('# ' + str(frames) + '\n')
        f.write('polymer number, rg squared \n')
        for mol in mols:
            f.write(str(mol) + ', ' + str(round(mols[mol] / frames, 4)) + '\n')
        avg_sum = round(sum / len(mols), 4)
        rg_square_root = round(math.sqrt(avg_sum), 4)
        f.write('------------------AVERAGE Square radius of gyration----------------\n')
        f.write(str(avg_sum))
        f.write('\n--------------Rg------------\n')
        f.write(str(rg_square_root))

# Function to calculate the coherent dynamic structure factor
def calculate_Sqt_single_q(timesteps, q_values, chain_length_range):
    Sqt = {q: [] for q in q_values}

    initial_timestep = min(timesteps.keys())
    initial_positions = {}
    for mol in timesteps[initial_timestep]["mols"]:
        initial_positions[mol] = {atom[0]: np.array(atom[3:6]) for atom in timesteps[initial_timestep]["mols"][mol]["atoms"]}

    for q in q_values:
        for timestep in sorted(timesteps.keys()):
            Sq_timestep = 0
            N = 0  # Initialize the number of atoms counter

            for mol in timesteps[timestep]["mols"]:
                mol_atoms = timesteps[timestep]["mols"][mol]["atoms"]
                mol_length = timesteps[timestep]["mols"][mol]["length"]
                
                # Filter based on chain length range
                if chain_length_range[0] <= mol_length <= chain_length_range[1]:
                    N += len(mol_atoms)

                    # Current positions of atoms in this molecule at current timestep
                    r_current = np.array([atom[3:6] for atom in mol_atoms])
                    # Initial positions of atoms in this molecule
                    r_initial = np.array([initial_positions[mol][atom[0]] for atom in mol_atoms])

                    # Compute all pairwise differences r_current - r_initial in a vectorized way
                    r_diff = r_current[:, np.newaxis, :] - r_initial[np.newaxis, :, :]

                    # Compute dot product q · (r_current - r_initial) for all pairs
                    dot_product = np.dot(r_diff, q)

                    # Compute sum of exponential terms in a vectorized way
                    sum_exp = np.sum(np.exp(1j * dot_product))

                    Sq_timestep += sum_exp.real

            if N > 0:  # Avoid division by zero if no chains are in the length range
                Sq_timestep /= N  # Normalize by the total number of atoms
            else:
                Sq_timestep = 0

            Sqt[q].append(round(Sq_timestep, 4))

    return Sqt

def write_Sqt_to_file(Sqt, q_values, timesteps, outFile):
    with open(outFile, 'w') as f:
        # Write header
        f.write("timestep," + ",".join(f"S(q,t) q={q}" for q in q_values) + "\n")
        
        # Write data for each timestep
        for i, timestep in enumerate(sorted(timesteps.keys())):
            f.write(f"{timestep}," + ",".join(f"{Sqt[q][i]}" for q in q_values) + "\n")

# Main execution
start_time = time.time()
result = readTimesteps("production.dump", 500)
q_values = [0.1, 0.2, 0.3, 0.4, 0.6, 1.0]
chain_length_range = (330, 395) #(350, 370)  # Specify the range of chain lengths to include
Sqt = calculate_Sqt_single_q(result, q_values, chain_length_range)
write_Sqt_to_file(Sqt, q_values, result, "Sqt_results_single_q_1.4_360.csv")
end_time = time.time()
elapsed_time_seconds = end_time - start_time
elapsed_time_minutes = elapsed_time_seconds / 60
print(f"Execution time: {elapsed_time_minutes:.2f} minutes")
print("DONE")
