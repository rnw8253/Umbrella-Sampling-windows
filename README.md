# Umbrella-Sampling-windows

These scripts produce windows from a trajectory for Umbrella Sampling. 

dih_windows.py and dih_sample.cfg choose windows based off of dihedral angles. 

dist_windows.py and dist_sample.cfg choose windows based off of center of mass distance between atoms/molecules.

FHF2_cpptraj.sh and cpptraj.in will take the windows produced and converts them into AMBER restart files. 