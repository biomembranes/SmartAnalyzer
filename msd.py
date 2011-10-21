import cfg
import math
import numpy
import time
import sys
import os
import re
from numpy import ndarray

def computeCOMxy(molarray, molstart, molend):
	x = y = 0.0
	for mol in range(molstart, molend+1):
		x+= float(molarray[mol][0])
		y+= float(molarray[mol][1])
	return (x/float(molend - molstart), y/float(molend - molstart))

def computeCOMdiff(molstart, molend, molarray, COMarray):
	diffdict = {}
	for mol in range(molstart, molend+1):
		diffdict[mol] = [float(molarray[mol][0]) - float(COMarray[0]), float(molarray[mol][1]) - float(COMarray[1])]
	return diffdict

def summationMSD(molarray, diffarr, diffarrinit, molstart, molend):
	sum = 0.0
	for mol in range(molstart, molend+1):
		xdelta = float(diffarr[mol][0]) - float(diffarrinit[mol][0])
		ydelta = float(diffarr[mol][1]) - float(diffarrinit[mol][1])
		sum+= xdelta**2 + ydelta**2
	return (sum/float(molend - molstart))

def modulus(atom, atomspermol):
	op = int(atom-1) % (int(atomspermol))
	return op+1

def computeMSD(filename, opfilename, molstart, molend, atomselect, atomspermol, startframe):
	f = open(filename, "r")
	w = open(opfilename, "w")
	molecules = molend - molstart + 1
	msddict = {}; moldict = {}
	collecttime = collectdata = gotnumtime = False
	for line in f:
		if "TIMESTEP" in line:
			collecttime = True
		if collecttime and not collectdata and not gotnumtime and re.search('\d+', line):
			time = int(re.search('\d+', line).group())
			gotnumtime = True
		if collecttime and "ITEM: ATOMS" in line:
			collectdata = True
		if collectdata and re.search('\d+', line):
			if (int(line.split()[2]) <= int(molecules)):
				atom_data = line.split()
				atomno = modulus(int(atom_data[0]), atomspermol)
				vectpos = atom_data[3], atom_data[4], atom_data[5]
				mol = int(atom_data[2])
				if (mol >= molstart) and (mol <= molend) and (atomno == atomselect):
					moldict[mol] = vectpos
				if int(line.split()[2]) == int(molecules) and (atomno == atomspermol):
					collecttime = collectdata = gotnumtime = False
					print "Finished frame: ", time
					#raw_input("Press ENTER to exit")
					if (int(time) == int(startframe)):
						initCOM = computeCOMxy(moldict, molstart, molend)
						initdiffdict = computeCOMdiff(molstart, molend, moldict, initCOM)
						diffdict = computeCOMdiff(molstart, molend, moldict, initCOM)
						
					else:
						currCOM = computeCOMxy(moldict, molstart, molend)
						diffdict = computeCOMdiff(molstart, molend, moldict, currCOM)
					if (int(time) > int(startframe)):
						MSD = summationMSD(moldict, diffdict, initdiffdict, molstart, molend)
						msddict[time] = MSD
						line = "%f %f\n" % (time, MSD)
						w.write(line)
					
						

	#print msddict
	w.close()
	f.close()

def main():
	computeMSD(cfg.inputtraj, cfg.outputMSD, cfg.molstart, cfg.molend, cfg.headatom, cfg.atomspermol, cfg.startframe)

if __name__ == "__main__":
	main()



	
"""def Print_Pos_arr(arr_pos):
	for mol in range(0,int_diffusion_num_molecules):
		print arr_pos[0][mol], arr_pos[1][mol] 
	raw_input("Press ENTER to exit")
	return ""


def DoMSD():
	## Now loop through all the times and generate MSD for all these times and OP

	MSD_array = []
	FILE_frame = open(cfg.frame_file, "r")
	line = FILE_frame.readline().strip()
	count = 0
	while (line != "</OpenMD>"):    
		#print line    
		if (line == "<FrameData>"): 
			count+=1
	   		time = float(FILE_frame.readline().split()[1])
			if (time >= float(cfg.int_start_frame)):







				#print "Time: ", time
				line = FILE_frame.readline().strip()
				while (line != "<StuntDoubles>"):
					line = FILE_frame.readline().strip()
				#print "Reading data..."
				mol = 0
				for my_bead in range(0,cfg.int_beads_per_mol*cfg.int_number_molecules_per_leaflet):
					line = FILE_frame.readline().strip()
					my_split = line.split()
					bead_number = int(my_split[0]) # bead number
					if ((bead_number % cfg.int_beads_per_mol) == 0):
						#print "Bead number: ", bead_number
						
						if (time == 0 or time == float(cfg.int_start_frame)):
							arr_initial_pos[0][mol] = float(my_split[2])
							arr_initial_pos[1][mol] = float(my_split[3])
							##arr_pos[0][mol] = float(my_split[2])
							##arr_pos[1][mol] = float(my_split[3])
							#print "",mol, arr_initial_pos[0][mol], arr_initial_pos[1][mol]
							## Get initial positions...
							#raw_input("Press ENTER to exit")
						else:
							arr_pos[0][mol] = float(my_split[2])
							arr_pos[1][mol] = float(my_split[3])
							#arr_pos[2][mol] = float(my_split[4])
					
						mol+=1
				if (time == 0 or time == float(cfg.int_start_frame)):
					initial_COM = Calculate_COM(arr_initial_pos)
					#print "Initial time=0 COM: ", initial_COM
					arr_initial_diff = Calculate_COM_diff(arr_initial_pos,initial_COM)
					arr_diff = Calculate_COM_diff(arr_initial_pos,initial_COM)
			
					#raw_input("Press ENTER to exit")
				else:
					arr_COM = Calculate_COM(arr_pos)
					#print "COM: ", time, arr_COM
					arr_diff = Calculate_COM_diff(arr_pos,arr_COM)
					#Print_Pos_arr(arr_diff)
					#print "t=: ",time, " COM:= ",arr_COM, " ", arr_diff
					#print "", arr_pos[0][mol], arr_pos[1][mol]

				MSD = Summation_MSD(arr_initial_diff, arr_diff)
				time_MSD_pair = [time,MSD]
				MSD_array.append(time_MSD_pair) 
				#print "MSD: ", time, MSD
				#raw_input("Press ENTER to exit")
		line = FILE_frame.readline().strip()
	FILE_frame.close()

	#print MSD_array

	FILE_MSD = open(cfg.output_MSD_file, "w")
	for each in MSD_array:
		print each
		FILE_MSD.write(str(each[0]) + "   "+ str(each[1]) + str("\n"))
	FILE_MSD.close()"""
