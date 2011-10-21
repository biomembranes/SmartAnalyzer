import math

##

inputtraj = "dumplong.trj"
nummols = 128
atomspermol = 14


### SCD ###

### Output data files for Scd Tail ordering analysis

outputScd  = "Scd.dat"
outputmScd = "mScd.dat"
output2Sxx = "2Sxx.dat"
outputSzz  = "Szz.dat"

### Start and end atoms
### Remember atoms n-1 and n+1 will be pulled in

bilayernormal = [0, 0, -1]
start = 2
end = 7


### MSD ###

molstart = 1
molend = 128
startsteps = 0 ### When to start the MSD from.
headatom = 1 ### This is what we take as our reference position.
startframe = 1000

outputMSD = "MSD.xy"
 
### Data Formats

ROWDATA  = 7 ## This is the number of lumps of data we take from LAMMPS data file
ROWBOND  = 4 ## Bond data from the LAMMPS data file
ROWANGLE = 5 ## Angle data from the LAMMPS data file
