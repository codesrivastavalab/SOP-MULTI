#	Code to generate pair list using lammps data file
#	Krishnakanth B, 
#	Theoretical Biophysics Laboratory, Molecular Biophysics Unit,
#	Indian Institute of Science, Bangalore - 560012
#
#	Last Modified: 21 March 2023

# Script usage syntax
# python3 SOP-data2pair_list.py data.hnrnpa1 lists.txt 

#  Load the required packages
import numpy as np
import sys
from essential_functions import *

data_fname = sys.argv[1]
list_fname = sys.argv[2]
data_file = open(data_fname,'r').readlines()
list_file  = open(list_fname,"w")

epsilon_BB = 0.45 # kcal/mol
epsilon_BS = 0.45 # kcal/mol
epsilon_SS = 0.18 # kcal/mol
d_cut = 20.0 #distance cutoff 
r_cut = 2 #residue_cutoff 



atom_locs = []
for i in range(len(data_file)):
	if('atoms' in data_file[i]):
		atom_locs.append(i)
	if('Atoms' in data_file[i]):
		atom_locs.append(i)

natoms = int(data_file[atom_locs[0]].split(" ")[0])

atom_ids = np.loadtxt(data_fname, skiprows=atom_locs[1]+1, max_rows = natoms, dtype= np.int64, usecols = [0])
domain_id = np.loadtxt(data_fname, skiprows=atom_locs[1]+1, max_rows = natoms, dtype= np.int64, usecols = [2])
residue_id = np.loadtxt(data_fname, skiprows=atom_locs[1]+1, max_rows = natoms, dtype= np.int64, usecols = [3])
atom_type = np.loadtxt(data_fname, skiprows=atom_locs[1]+1, max_rows = natoms, dtype= np.int64, usecols = [4])
x = np.loadtxt(data_fname, skiprows=atom_locs[1]+1, max_rows = natoms, dtype= np.float64, usecols = [6,7,8])
atom_BS = np.repeat("CB", len(atom_type))
atom_BS[np.where(atom_type==1)] = "CA"
atom_BS[np.where(atom_type==2)] = "CA"


bead_types = np.loadtxt("sop_model_library/sop_bead_types.dat",delimiter="\t",dtype = np.int64, usecols = (0),skiprows =1 )
bead_sigma = np.loadtxt("sop_model_library/sop_bead_types.dat",delimiter="\t",dtype = np.float64, usecols = (3),skiprows =1 )
bead_BS = np.loadtxt("sop_model_library/sop_bead_types.dat",delimiter="\t",dtype = str, usecols = (1),skiprows =1 )


np.stack((bead_types, bead_sigma),axis =1)

def sigma_reader(a, b):
	if a[1] == "CA":
		a[0] ="CA"
	if b[1] == "CA":
		b[0] ="CA"
	s1 = float(bead_data[np.where(bead_data[:,1]==a[0])[0],2][0])
	s2 = float(bead_data[np.where(bead_data[:,1]==b[0])[0],2][0])	# what a mistake
	return str(round(s1+s2,3))

# check this
#/media/krishnakanth/krishnakanth-disk-3/others/clone/tia
#/media/krishnakanth/krishnakanth-disk-2/lammps/test/lammps-sop-idp/SOP-HYBRID
for i in range(len(atom_ids)):
	for j in range(len(atom_ids)):
		if(i<j):
			atom_i = atom_ids[i]
			atom_j = atom_ids[j]
			domain_i = domain_id[i]
			domain_j = domain_id[j]
#			atom_type_i = atom_type[i]
#			atom_type_j = atom_type[j]
			coord_i = x[i,:]
			coord_j = x[j,:]
			i_BS = atom_BS[i]
			j_BS = atom_BS[j]
			residue_i = residue_id[i]
			residue_j = residue_id[j]
#			sigma = float(bead_sigma[np.where(bead_types==atom_type_i)]+bead_sigma[np.where(bead_types==atom_type_j)])
#			sigma = round(sigma,2)
			e_dist = euclidean_distance(coord_i, coord_j)
			c_1 = (domain_i == domain_j) and (domain_i!=0) and (domain_j!=0) # condition 1
			c_2 = (e_dist < d_cut) # condition 2
			c_3 = (abs(residue_i-residue_j)>r_cut) # condition 3

			if c_1 and c_2 and c_3 :
				if (i_BS == "CA") and  (j_BS =="CA" ):
					print(str(atom_i)+"\t"+str(atom_j)+"\t"+"ljsop\t"+str(epsilon_BB)+"\t"+str(e_dist),file =list_file)
				if (i_BS == "CB") and  (j_BS =="CB" ):
					print(str(atom_i)+"\t"+str(atom_j)+"\t"+"ljsop\t"+str(epsilon_SS)+"\t"+str(e_dist),file =list_file)
				if(i_BS != j_BS):
					print(str(atom_i)+"\t"+str(atom_j)+"\t"+"ljsop\t"+str(epsilon_SS)+"\t"+str(e_dist),file =list_file)

list_file.close()
