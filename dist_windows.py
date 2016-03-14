# NOTE: will have to point to version of python and have MDAnalysis library

# USAGE: dist_windows.py [config file]
# CONFIG FILE FORMAT:
# psffile = [psf, prmtop, gro, or pdb file]
# dcdfile = [traj file in format trr, dcd, etc]
# atom_sel_1 = [CHARMM style atom selection for first group of atoms]
# atom_sel_2 = [CHARMM style atom selection for second group of atoms]
# dist_max = [maximum distance]
# dist_min = [minimum distance]
# dist_delta = [distance bin size]


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


# define subroutines
def ParseConfigFile(cfg_file):
	global inp_psf_file,inp_dcd_file,user_sel_1,user_sel_2,user_sel_3,user_sel_4,dist_min,dist_max,dist_delta
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
			elif option.lower()=='dist_min':
				dist_min = float(value)
			elif option.lower()=='dist_max':
				dist_max = float(value)
			elif option.lower()=='dist_delta':
				dist_delta = float(value)
			else :
				print "Option:", option, " is not recognized"
	

def computePbcDist(r1,r2,box):
	dist = 0

	for j in range(0,3):
		temp = r1[j]-r2[j]
		if temp < -box[j]/2.0:
			temp += box[j]
		elif temp > box[j]/2.0:
			temp -= box[j]
		dist += temp*temp

	dist = math.sqrt(dist)
	return dist;

def computeDist(r1,r2,box):
	dist = 0

	for j in range(0,3):
		temp = r1[j]-r2[j]
		#if temp < -box[j]/2.0:
		#	temp += box[j]
		#elif temp > box[j]/2.0:
			#temp -= box[j]
		dist += temp*temp

	dist = math.sqrt(dist)
	return dist;

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

def wrapPositions(old_positions,n_atoms,box,target):
	new_positions = numpy.empty( (n_atoms,3), dtype=float)
	for i in range(0,n_atoms):
		for j in range(0,3):
			if (old_positions[i,j]-target[j]) < -box[j]/2.0:
				#print "Changing position of atom", i+1
				new_positions[i,j] = old_positions[i,j] + box[j]
			elif (old_positions[i,j]-r1[j]) > box[j]/2.0:
				#print "Changing position of atom", i+1
				new_positions[i,j] = old_positions[i,j] - box[j]
			else:
				new_positions[i,j] = old_positions[i,j]
	# update the MDAnalysis positions
	return new_positions

#################################################################################################################
##############################################     MAINE PROGRAM   ##############################################
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

# allocate COM position vectors
r1 = numpy.empty(3,dtype=float)
r2 = numpy.empty(3,dtype=float)

# define ATOMGROUP that contains all atoms for printing purposes
all_atoms = coord.selectAtoms("all")

# determine number of distance bins
n_dist_bins = int((dist_max-dist_min)/dist_delta)
dist_pop = numpy.empty(n_dist_bins,dtype=int)
for i in range(0,n_dist_bins):
	dist_pop[i]=0

for ts in coord.trajectory:
	# get box dimensions
	box = coord.trajectory.ts.dimensions[:3]
	# wrap coordinates
	r1 = sel1.centerOfMass()
	#sel1.positions = wrapPositions(sel1.positions,len(sel1),box,r1)
	#sel2.positions = wrapPositions(sel2.positions,len(sel2),box,r2)
	# compute new centers of mass 
	r1 = sel1.centerOfMass()
	r2 = sel2.centerOfMass()
	# compute distance
	dist = computeDist(r1,r2,box)
	print ts.frame, dist
	# determine distance bin
	dist_bin = int((dist-dist_min)/dist_delta)
	# check to see if it falls in range and if we have already populated that window
	if dist_bin >=0 and dist_bin < n_dist_bins and dist_pop[dist_bin]==0:
		#all_atoms.positions = wrapPositions(all_atoms.positions,len(all_atoms),box,r1)
		dist_pop[dist_bin]=1
		window = dist_bin+1
		window_pdb_file = "window"+str(window)+".pdb"
		print "Writing pdb file:", window_pdb_file, "with a distance of:", dist
		all_atoms.write(window_pdb_file)


	

