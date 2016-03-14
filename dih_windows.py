#!/Users/martinmccullagh/anaconda/bin/python
# NOTE: will have to point to version of python and have MDAnalysis library

# USAGE: dih_windows.py [config file]
# CONFIG FILE FORMAT:
# psffile = [psf, prmtop, gro, or pdb file]
# dcdfile = [traj file in format trr, dcd, etc]
# atom_sel_1 = [CHARMM style atom selection for first group of atoms]
# atom_sel_2 = [CHARMM style atom selection for second group of atoms]
# atom_sel_3 = [CHARMM style atom selection for third group of atoms]
# atom_sel_4 = [CHARMM style atom selection for fourth group of atoms]
# dih_max = [maximum dihedral in degrees]
# dih_min = [minimum dihedral in degrees]
# dih_delta = [dihedral bin size in degrees]

# load libraries
import scipy
import sys
import os
import numpy
import math
import MDAnalysis
import numpy.linalg

#################################################################################################################
##############################################     SUBROUTINEs     ##############################################
#################################################################################################################

# read the configuration file and populate the global variables
def ParseConfigFile(cfg_file):
	global inp_psf_file,inp_dcd_file,user_sel_1,user_sel_2,user_sel_3,user_sel_4,dih_min,dih_max,dih_delta
	f = open(cfg_file)
	for line in f:
		# first remove comments
		if '#' in line:
			line, comment = line.split('#',1)
		if '=' in line:
			option, value = line.split('=',1)
			option = option.strip()
			value = value.strip()
			print "Option:", option, " Value:", value
			# check value
			if option.lower()=='psffile':
				inp_psf_file = value
			elif option.lower()=='dcdfile':
				inp_dcd_file = value
			elif option.lower()=='atom_sel_1':
				user_sel_1 = value
			elif option.lower()=='atom_sel_2':
				user_sel_2 = value
			elif option.lower()=='atom_sel_3':
				user_sel_3 = value
			elif option.lower()=='atom_sel_4':
				user_sel_4 = value
			elif option.lower()=='dih_min':
				dih_min = float(value)
			elif option.lower()=='dih_max':
				dih_max = float(value)
			elif option.lower()=='dih_delta':
				dih_delta = float(value)
			else :
				print "Option:", option, " is not recognized"
	
# compute the distance between two position vectors taking into account PBC
def computePbcDist(r1,r2,box):
	dist = 0

	for j in range(0,3):
		temp = r1[j]-r2[j]
		if temp < -box[j]/2.0:
			temp += box[j]
		elif temp > box[j]/2.0:
			temp -= box[j]
		dist += temp*temp

	dist = sqrt(dist)
	return dist;

#compute the diheral angle in degrees between four ordered position vectors
def computeDih(r1,r2,r3,r4):

	# define the distance vectors
	b1 = r1 - r2
	b2 = r2 - r3
	b3 = r3 - r4
	# compute the cross product vectors
	A = numpy.cross(b1,b2)
	A = A/math.sqrt(numpy.dot(A,A))
	B = numpy.cross(b2,b3)
	B = B/math.sqrt(numpy.dot(B,B))
	C = numpy.cross(b2,A)
	C = C/math.sqrt(numpy.dot(C,C))
	# now compute the dihedral
	dih = -numpy.arctan2(numpy.dot(C,B),numpy.dot(A,B))

	# convert to degrees
	dih *= 180.0/3.1415926535
	return dih

#################################################################################################################
##############################################    MAIN PROGRAM     ##############################################
#################################################################################################################

# read in command line argument
cfg_file = sys.argv[1]

# read cfg file
ParseConfigFile(cfg_file)

print "PSF file:", inp_psf_file
print "Coord DCD file:", inp_dcd_file

# start MDAnalysis with a universal
#coordinate universe
coord = MDAnalysis.Universe(inp_psf_file,inp_dcd_file)

# print some general log info
print "Numer of time steps in coordinate trajectory:", len(coord.trajectory)

# define the first selection
sel1 = coord.selectAtoms(user_sel_1)
print "First atom selection:", user_sel_1
print sel1
# define the second atom selection
sel2 = coord.selectAtoms(user_sel_2)
print "Second atom selection:", user_sel_2
print sel2
# define the third atom selection
sel3 = coord.selectAtoms(user_sel_3)
print "Third atom selection:", user_sel_3
print sel3
# detine the fourth atom selection
sel4 = coord.selectAtoms(user_sel_4)
print "Fourth atom selection:", user_sel_4
print sel4

r1 = numpy.empty(3,dtype=float)
r2 = numpy.empty(3,dtype=float)
r3 = numpy.empty(3,dtype=float)
r4 = numpy.empty(3,dtype=float)

all_atoms = coord.selectAtoms("all")

# determine number of dihedral bins
n_dih_bins = int((dih_max-dih_min)/dih_delta)
dih_pop = numpy.empty(n_dih_bins,dtype=int)
for i in range(0,n_dih_bins):
	dih_pop[i]=0

# loop through all time steps in trajectory
for ts in coord.trajectory:
	# get box dimensions
	box = coord.trajectory.ts.dimensions[:3]
	# compute 4 COM position vectors from atom selections 
	r1 = sel1.centerOfMass()
	r2 = sel2.centerOfMass()
	r3 = sel3.centerOfMass()
	r4 = sel4.centerOfMass()
	# compute distance
	dih = computeDih(r1,r2,r3,r4)
	# determine dihedral bin
	dih_bin = int((dih-dih_min)/dih_delta)
	if dih_bin >=0 and dih_bin < n_dih_bins and dih_pop[dih_bin]==0:
		dih_pop[dih_bin]=1
		window = dih_bin+1
		window_pdb_file = "window"+str(window)+".pdb"
		print "Writing pdb file:", window_pdb_file, "with a dihedral of:", dih
		all_atoms.write(window_pdb_file)


	

