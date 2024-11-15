# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.


#This script includes contributions of multiple authors; an original version was written by Mark Stevens, 
#and it was modified in the research group of Lisa Hall and later by the Taofeek Tejuosho from the Janani Sampath group


import math
import numpy as np
import csv
import time 

start_time = time.time()

def readTimesteps(inputFile, timestepsCount = None):
	if timestepsCount:
		counter = 0
	else:
		timestepsCount = 5 #999999999999999999999999
	timesteps = {} #created an empty dictionary
	with open(inputFile) as f:
		line = f.readline()
		#print(line)
		while line != '':
			if line.split()[0] == "ITEM:":
				timestep = int(f.readline().split()[0]) # particular Timestep number
				timesteps[timestep] = {}
				filler = f.readline()
				atoms = int(f.readline().split()[0]) #number of atoms
				timesteps[timestep]["atoms"] = atoms #
				filler = f.readline()
				x = f.readline() #box dimension in x direction
				x_pos = float(x.split()[1]) #positive x direction
				#print(x_pos)
				x_neg = float(x.split()[0])
				#print("x negative is", x_neg)
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
				#print(timesteps[timestep]['z'])
				headers = f.readline()
				timesteps[timestep]["headers"] = headers
				timesteps[timestep]["mols"] = {}
				#print("The timestep is", timesteps, "\n", "\n", "\n  ")
				#print("This is the dictionary that stores our data ", timesteps[timestep])
				#print("\n")
				#print("\n")
				for i in range(atoms):
					line = f.readline()
					id = int(line.split()[0])
					#print(id)
					mol = int(line.split()[1])
					type = int(line.split()[2])
					q = int(line.split()[3])
					xs = float(line.split()[4]) * (timesteps[timestep]['x']-timesteps[timestep]['-x'])
					ys = float(line.split()[5]) * (timesteps[timestep]['y']-timesteps[timestep]['-y'])
					zs = float(line.split()[6]) * (timesteps[timestep]['z']-timesteps[timestep]['-z'])
					ix = int(line.split()[7])
					iy = int(line.split()[8])
					iz = int(line.split()[9])
					xbox = (timesteps[timestep]['x']-timesteps[timestep]['-x'])
					ybox = (timesteps[timestep]['y']-timesteps[timestep]['-y'])
					zbox = (timesteps[timestep]['z']-timesteps[timestep]['-z'])
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
					xs = xbox*float(line.split()[4]) + ix*xbox
					ys = ybox*float(line.split()[5]) + iy*ybox
					zs = zbox*float(line.split()[6]) + iz*zbox
					#print("this is where the timestep", timesteps[timestep])
					if mol in timesteps[timestep]["mols"]: #arranging the trajectory into dictionary
						#print("this is the mol", timesteps[timestep]["mols"])
						#print("Before",timesteps[timestep]["mols"][mol]["atoms"])
						timesteps[timestep]["mols"][mol]["atoms"].append((id, type, q, xs, ys, zs, ix, iy, iz))
						#print("After",timesteps[timestep]["mols"][mol]["atoms"])
					else:
						timesteps[timestep]["mols"][mol] = {}
						#print("from the else",timesteps[timestep]["mols"][mol] )
						timesteps[timestep]["mols"][mol]["atoms"] = []
						#print("from the else",timesteps[timestep]["mols"][mol]["atoms"] )
						timesteps[timestep]["mols"][mol]["atoms"].append((id, type, q, xs, ys, zs, ix, iy, iz))
				for mol in timesteps[timestep]["mols"]:
					x = 0
					y = 0
					z = 0
					length = len(timesteps[timestep]["mols"][mol]["atoms"]) #Accessing the dictionary; dic:key:value:key:value
					#print("this is the length", len(timesteps[timestep]["mols"][mol]["atoms"]))
					timesteps[timestep]["mols"][mol]["length"] = length
					#print("the length of a particular chain", timesteps[timestep]["mols"][mol]["length"])
					#print(timesteps[timestep]["mols"])
					for i in range(length):
						x = x + timesteps[timestep]["mols"][mol]["atoms"][i][3]/length
						y = y + timesteps[timestep]["mols"][mol]["atoms"][i][4]/length
						z = z + timesteps[timestep]["mols"][mol]["atoms"][i][5]/length
					timesteps[timestep]["mols"][mol]["xcm"] = x
					timesteps[timestep]["mols"][mol]["ycm"] = y
					timesteps[timestep]["mols"][mol]["zcm"] = z
					#print("the length of this", timesteps[timestep]["mols"])
			line = f.readline()
			#print("This is the line ",line)
			counter = counter + 1
			if counter == timestepsCount:
				return timesteps
			#print("read timestep: ", timestep)
			#print("final timestep", timesteps)
	return timesteps

#print(readTimesteps("production.dump", 1))
# Calculate Rg2 (radius of gyration) of whatever system inputed.
# Timesteps input MUST be object like dictionary from readTimesteps(....)
def Rg2Calculator(timesteps):
	for timestep in timesteps:
		for mol in timesteps[timestep]["mols"]:
			length = timesteps[timestep]["mols"][mol]["length"]
			rg2 = 0
			for i in range(length):
				x = timesteps[timestep]["mols"][mol]["atoms"][i][3]
				y = timesteps[timestep]["mols"][mol]["atoms"][i][4]
				z = timesteps[timestep]["mols"][mol]["atoms"][i][5]
				x = x - timesteps[timestep]["mols"][mol]["xcm"]
				y = y - timesteps[timestep]["mols"][mol]["ycm"]
				z = z - timesteps[timestep]["mols"][mol]["zcm"]
				x2 = x ** 2
				y2 = y ** 2
				z2 = z ** 2
				rg2 = ((x2 + y2 + z2)/length)  + rg2
				squarerootrg = math.sqrt(rg2)
				#print(squarerootrg)

			timesteps[timestep]["mols"][mol]["rg2"] = rg2
			#print("The final timesteps is ", timesteps)
	return timesteps

def writeRg2ToFile(timesteps, outFile):
	with open(outFile, 'w') as f:
		f.write("rg2, molNumber, timestep")
		f.write("\n")
		for timestep in timesteps:
			for mol in timesteps[timestep]["mols"]:
				rg2 = str(timesteps[timestep]["mols"][mol]["rg2"])
				molNumber = str(mol)
				timestepWrite = str(timestep)
				f.write(rg2 + ", " + molNumber + ", " + timestepWrite)
				f.write("\n")
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
            f.write(str(mol) + ', ' + str(mols[mol] / frames) + '\n')
        sum = sum / len(mols)
        rg_square_root = math.sqrt(sum)
        f.write('------------------AVERAGE Square radius of gyration----------------\n')
        f.write(str(sum))
        f.write('\n--------------Rg------------\n')
        f.write(str(rg_square_root))

def calculate_structure_factor(timesteps, q_values, min_length=None, max_length=None):
    structure_factors = {}
    for timestep in timesteps:
        for mol in timesteps[timestep]["mols"]:
            N = timesteps[timestep]["mols"][mol]["length"]
            if (min_length is not None and N < min_length) or (max_length is not None and N > max_length):
                continue
            positions = np.array([(atom[3], atom[4], atom[5]) for atom in timesteps[timestep]["mols"][mol]["atoms"]])
            S_q = np.zeros_like(q_values, dtype=np.float64)
            for i, q in enumerate(q_values):
                sum_exp = np.sum(np.exp(1j * q * np.linalg.norm(positions, axis=1)))
                S_q[i] = (1 / N) * np.abs(sum_exp) ** 2
            if (timestep, mol) not in structure_factors:
                structure_factors[(timestep, mol)] = S_q
            else:
                structure_factors[(timestep, mol)] += S_q
    averaged_S_q = np.mean(np.array(list(structure_factors.values())), axis=0)
    return q_values, averaged_S_q

def write_structure_factor_to_csv(q_values, averaged_S_q, outFile):
    with open(outFile, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["q", "S(q)"])
        for q, S_q_value in zip(q_values, averaged_S_q):
            writer.writerow([q, S_q_value])

# Main execution
input_file = "production.dump"
output_file = "sofk_1.4_360_SPECIFIC.csv"
q_values = np.logspace(-2, 2, 100)  # Define a range of q values

result = readTimesteps(input_file, 100)
#result = Rg2Calculator(result)

# Define the chain length range for static structure factor calculation
min_length = 350  # Example minimum chain length
max_length = 370  # Example maximum chain length

q_values, averaged_S_q = calculate_structure_factor(result, q_values, min_length, max_length)
write_structure_factor_to_csv(q_values, averaged_S_q, output_file)

end_time = time.time()
elapsed_time_seconds = end_time - start_time
elapsed_time_minutes = elapsed_time_seconds/60
print(f"Execution time: {elapsed_time_minutes:.2f} minutes")

print("DONE calculation")
