#	Code to convert all-atom pdb file to SOP CG model pdb, considering COM of side chain atoms for CB atoms
#	Krishnakanth B, 
#	Theoretical Biophysics Laboratory, Molecular Biophysics Unit,
#	Indian Institute of Science, Bangalore - 560012
#
#	Last Modified:  10 April 2023
#
# Usage example: python3 SOP-AA2CG-builder.py hnrnpa1-full_with_H.pdb hnrnpa1-full_CG.pdb



#  Load the required packages
import sys
import numpy as np
from essential_functions import *

aa_filename = sys.argv[1]
cg_pdb = open(sys.argv[2],"w")
mass_dict={'H':1.00800,'C':12.01100,'N':14.00700,'O':15.99900,'S':32.06000}
side_chain_dict = {}
side_chain_table = np.loadtxt("sop_model_library/com_bead_without_H.dat",delimiter="\t", dtype=str)
for i in side_chain_table:
	side_chain_dict[i[0]+"_"+i[1]] = i[2].split()
	

x,y,z, atom_id, atom_type, res_type, res_no, chain_id = pdb_reader(aa_filename)

res_diff = 0
# adjust res numbers if required
if (int(res_no[0])-1 != 0):
	res_diff = int(res_no[0])-1

res_no = np.array([str(int(i)-res_diff) for i in res_no])
x = np.array(x)
y = np.array(y)
z = np.array(z)
atom_id = np.array(atom_id)
atom_type = np.array(atom_type)
res_type = np.array(res_type)
clean_H(atom_type)   # using referencing concept in python


atom_counter = 0
res_counter = 1

add_on ="C"
for k in range(1, int(res_no[-1])+1):
	res_atoms = np.where(res_no==str(k))[0]
	r_dict = {}
	r_rtype = np.unique(res_type[res_atoms])
	for i in res_atoms:
		r_dict[atom_type[i]] = [x[i],y[i], z[i]]

	print(r_rtype)
	if (r_rtype[0]=="GLY"):
		atom_counter = atom_counter+1
		atom_info = str.rjust(str(atom_counter),5)+" "+str.rjust("CA",3)
		res_info = str.rjust((r_rtype[0]),4) + " "+"A"+str.rjust(str(k),4)
		coord_info = str.rjust('%.3f'%(r_dict['CA'][0]),8)+str.rjust('%.3f'%(r_dict['CA'][1]),8)+str.rjust('%.3f'%(r_dict['CA'][2]),8)
		print("ATOM  "+atom_info+" "+res_info+"    "+coord_info+"  1.00  0.00           "+add_on,file = cg_pdb)
	
	else:
		atom_counter = atom_counter+1
		atom_info = str.rjust(str(atom_counter),5)+" "+str.rjust("CA",3)
		res_info = str.rjust(r_rtype[0],4) + " "+"A"+str.rjust(str(k),4)
		coord_info = str.rjust('%.3f'%(r_dict['CA'][0]),8)+str.rjust('%.3f'%(r_dict['CA'][1]),8)+str.rjust('%.3f'%(r_dict['CA'][2]),8)
		print("ATOM  "+atom_info+" "+res_info+"    "+coord_info+"  1.00  0.00           "+add_on,file = cg_pdb)

		atom_counter = atom_counter+1
		com_x, com_y, com_z = 0, 0, 0
		mass = 0
		for j in side_chain_dict[r_rtype[0]+"_CB"]:
			com_x = com_x + mass_dict[j[0]]*r_dict[j][0]
			com_y = com_y + mass_dict[j[0]]*r_dict[j][1]
			com_z = com_z + mass_dict[j[0]]*r_dict[j][2]
			mass = mass + mass_dict[j[0]]
		com_x, com_y, com_z = com_x/mass, com_y/mass, com_z/mass

		atom_info = str.rjust(str(atom_counter),5)+" "+str.rjust("CB",3)
		res_info = str.rjust(r_rtype[0],4) + " "+"A"+str.rjust(str(k),4)
		coord_info = str.rjust('%.3f'%(com_x),8)+str.rjust('%.3f'%(com_y),8)+str.rjust('%.3f'%(com_z),8)
		print("ATOM  "+atom_info+" "+res_info+"    "+coord_info+"  1.00  0.00           "+add_on,file = cg_pdb)			

#ATOM      1   CA  ASP    1       2.901   5.646  -1.753  1.00  0.00      C
#ATOM      1  N   MET A   1     -18.993   2.728  20.883  1.00 69.10           N  
cg_pdb.close()
