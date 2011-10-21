#######################################################################################
## :: Compute Tail Ordering of Lipid Systems in LAMMPS ::                            ##
##                                                                                   ##
##  Reference for algorithm: 		                                             ##
##  Molecular dynamics simulation of a smectic liquid crystal with atomic detail     ##
##                                                            Brett Donovan 2011     ##
#######################################################################################

import cfg
import math
import numpy
import time
import sys
import os
import re
from numpy import ndarray

### Maths and Vector ###

def kroneckerdelta(alpha,beta):
	op=0.0
	if (alpha==beta):
		op=1.0;
	return op
	
def crossproduct(u,v):
	opv = [0, 0, 0]
	opv[0] =   u[1]*v[2] - u[2]*v[1]
	opv[1] = -(u[0]*v[2] - u[2]*v[0])
	opv[2] =   u[0]*v[1] - u[1]*v[0]
	return opv

def dotproduct(u,v):
	opv = u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
	return opv

def normalizevector(u):
	mag = math.sqrt(dotproduct(u,u))
	opv = [0, 0, 0]
	for i in range(3):
		opv[i] = u[i]/mag
	return opv

def vectordiff(u,v):
	opv = [0, 0, 0]
	for i in range(3):
		opv[i] = v[i] - u[i]
	return opv

def anglebetweenvect(u,v):
	u = normalizevector(u)
	v = normalizevector(v)
	mag = dotproduct(u,v)
	theta = math.acos(mag)
	return theta

def computecomponents(atom_n_minus_1, atom_n, atom_n_plus_1):
	#z: vector from Cn _ 1 to Cn + 1
	#y: vector, perpendicular to z, and in the plane through Cn - 1 , Cn, and Cn + 1;
	#x: vector, perpendicular to z and y
	z = normalizevector(vectordiff(atom_n_minus_1, atom_n_plus_1))
	tempA = vectordiff(atom_n_minus_1, atom_n)
	tempB = vectordiff(atom_n, atom_n_plus_1)
	tempA = normalizevector(tempA)
	tempB = normalizevector(tempB)
	x = crossproduct(tempA,tempB)
	x = normalizevector(x)
	y = crossproduct(z,x)
	y = normalizevector(y)
	return x, y, z

def getCalphas(molarray, atomindex):
	C_n_m1 = molarray[atomindex-1][0], molarray[atomindex-1][1], molarray[atomindex-1][2]
	C_n    = molarray[atomindex][0], molarray[atomindex][1], molarray[atomindex][2]
	C_n_p1 = molarray[atomindex+1][0], molarray[atomindex+1][1], molarray[atomindex+1][2]
	return (C_n_m1, C_n, C_n_p1)

def computeSij(bilayernormalv, i, j, atom_n_minus_1, atom_n, atom_n_plus_1):
	vectarr = computecomponents(atom_n_minus_1, atom_n, atom_n_plus_1)
	anglei = math.cos(anglebetweenvect(vectarr[i], bilayernormalv))
	anglej = math.cos(anglebetweenvect(vectarr[j], bilayernormalv))
	op = (1.0/2.0)*(3*anglei*anglej - kroneckerdelta(i,j))
	return op

def selectatomsforanalysis(bilayernormalv, intstartn, intendn, molarray, molecules, Sijarr):
	count = 0
	for mol in range(1, molecules+1):
		for atom in range(intstartn, intendn+1):
			C_n_m1, C_n, C_n_p1 = getCalphas(molarray[mol], int(atom))
			for i in range(3):
				for j in range(3):
					Sijarr[atom][i][j]+= computeSij(bilayernormalv, i, j, C_n_m1, C_n, C_n_p1);

			count+=1
	return count
		
def modulus(atom, atomspermol):
	op = int(atom-1) % (int(atomspermol))
	return op+1

def computetailordering(bilayernormalv, trajectoryfile, molecules, atomspermol):
	### Process the current frame of data...
	### Stubbed out time range
	vdimmol = molecules+1, atomspermol+1, 3
	vdimSij = atomspermol+1, 3, 3
	molarray = ndarray(vdimmol)
	Sijarr = ndarray (vdimSij)	
	collecttime = collectdata = gotnumtime = False
	maincount = 0
	f = open(trajectoryfile)
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
				molarray[mol][atomno] = vectpos
				if int(line.split()[2]) == int(molecules) and (atomno == atomspermol):
					collecttime = collectdata = gotnumtime = False
					print "Finished frame: ", time
					if (int(time) > 0):
						maincount += selectatomsforanalysis(bilayernormalv, (cfg.start), (cfg.end), molarray, molecules, Sijarr)
	Sijarr = Sijarr/maincount
	f.close()
	return Sijarr
	
		
def printScd(Sijarr, atom):
	S_cd = (2.0/3.0)*Sijarr[atom][0][0] + (1.0/3.0)*Sijarr[atom][1][1]
	return S_cd

def outputcomponents(Sijarr, prefactor, i, j, startint, endint, filename):
	f = open(filename, "w")
	for n in range(startint, endint+1):
		line = "%d %f \n" % (n, prefactor*Sijarr[n][i][j])	
		f.write(line)
	f.close()
	return 0

def outputcombined(prefactor, Sijarr, startint, endint, filename):
	f = open(filename, "w")
	for n in range(startint, endint+1):
		line = "%d %f \n" % (n, prefactor*printScd(Sijarr, n))	
		f.write(line)
	f.close()
	return 0

def main():
	Sijarr = computetailordering((cfg.bilayernormal), cfg.inputtraj, 128, 14)
	outputcomponents(Sijarr, -2.0, 0, 0, (cfg.start), (cfg.end), (cfg.output2Sxx))
	outputcomponents(Sijarr,  1.0, 2, 2, (cfg.start), (cfg.end), (cfg.outputSzz))
	outputcombined(1.0, Sijarr, (cfg.start), (cfg.end), (cfg.outputScd))
	outputcombined(-1.0, Sijarr, (cfg.start), (cfg.end), (cfg.outputmScd))

if __name__ == "__main__":
	main()


 
