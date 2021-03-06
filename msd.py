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
	collecttime = collectdata = gotnumtime = initialized = False
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
					if (int(time) >= int(startframe)) and not initialized:
						print "Getting initial configuration..."
						initCOM = computeCOMxy(moldict, molstart, molend)
						initdiffdict = computeCOMdiff(molstart, molend, moldict, initCOM)
						diffdict = computeCOMdiff(molstart, molend, moldict, initCOM)
						initialized = True
					else:
						currCOM = computeCOMxy(moldict, molstart, molend)
						diffdict = computeCOMdiff(molstart, molend, moldict, currCOM)
					if (int(time) > int(startframe)):
						print time, startframe
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

