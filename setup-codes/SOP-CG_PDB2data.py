#	Code to convert multi-molecule SOP-CG pdb file to LAMMPS data file
#	Krishnakanth B, 
#	Theoretical Biophysics Laboratory, Molecular Biophysics Unit,
#	Indian Institute of Science, Bangalore - 560012
#
#	Last Modified: 30 March 2024

#	Usage format example: 
#	python3 SOP-PDB2data.py hnrnpa1_CG.pdb data.hnrnpa1

#  Load the required packages
import os
import sys
import numpy as np
from essential_functions import *

# Read python code arguments
cg_filename = sys.argv[1]
nmol = 1
data_file = open(sys.argv[2],"w")


box_dimensions = [500.0, 500.0, 500.0] # box dimensions in angstroms

R0 = 2.0 # Tolerance of fene bond in ang
spring_constant_idp = 20.0  # kcal mol^{-1} ang^{2}
spring_constant_folded = 50.0  # kcal mol^{-1} ang^{2}


x,y,z, atom_id, atom_type, res_type, res_no, chain_id =pdb_reader(cg_filename)


natoms = int(len(x)/nmol)

if(natoms*nmol != len(x)):
	print("Something's wrong : Mismatch in atoms and number of molecules in the files")
	sys.exit(1)
	

# generate molecule id information
molecule_id = []
for i in range(1,nmol+1):
	molecule_id = molecule_id + [str(i)]*natoms

# generate charge information
charge = [0.0]*len(x)
for i in range(len(res_type)):
	if(atom_type[i]=="CB"):
		if ((res_type[i]=="ASP") or (res_type[i]=="GLU")):
			charge[i] = -1.0
		if ((res_type[i]=="ARG") or (res_type[i]=="LYS")):
			charge[i] = 1.0
	
# generate domain id information
if (os.path.exists("domain_data.dat")):
	domain_info  = np.loadtxt("domain_data.dat",dtype=int)
	domain_id = []
	for i in res_no:
		for j in range(np.shape(domain_info)[0]):
			if(int(i)<=domain_info[j,1] and int(i)>=domain_info[j,0]):
				domain_id.append(domain_info[j,2])
else:
	domain_id = [0]*len(x)

		

# check  size mismatches
if(len(domain_id)*nmol != len(x)):
	print("Something's wrong : Missing residue information in domain data file")
	sys.exit(1)
else:
	domain_id = domain_id*nmol

# Get the Bead info from bead_types.dat
f=open("sop_model_library/sop_bead_types.dat").readlines()
nbeads=int(f[0].strip().split()[0])
print("Detected "+str(nbeads))
bead_type_dict={}
for i in range(1,nbeads+1):
	temp=f[i].strip().split()
	bead_type_dict[temp[2]+"_"+temp[1]]=temp[0]



mols = unique_list(molecule_id)
multi_molecule_dict = {} 
print("Detected "+str(len(mols))+" molecules:")
print("Extracting topology of protein molecules..")
for mol in mols:
	print(mol+"\t",end="")
	molecule_dict={}
	bead_ids=which(mol, molecule_id)
	init_res=bead_ids[0]
	term_res=bead_ids[-1]+1
	residues=unique_list(res_no[init_res:term_res])
	for res in residues:
		temp=which(res,res_no)
		init_bead=temp[0]
		term_bead=temp[-1]
		res_dict={}
		for bead in temp:
			if(atom_type[bead]=="CA"):
				if res_type[bead]=="GLY":
					res_dict[res_type[bead]+"_"+atom_type[bead]]=[atom_id[bead],domain_id[bead],res_no[bead],bead_type_dict[res_type[bead]+"_"+atom_type[bead]],charge[bead],x[bead],y[bead],z[bead]]	
				else:
					res_dict["____"+atom_type[bead]]=[atom_id[bead],domain_id[bead],res_no[bead],bead_type_dict["____"+atom_type[bead]],charge[bead], x[bead],y[bead],z[bead]]
			else:
				res_dict[res_type[bead]+"_"+atom_type[bead]]=[atom_id[bead],domain_id[bead],res_no[bead],bead_type_dict[res_type[bead]+"_"+atom_type[bead]],charge[bead],x[bead],y[bead],z[bead]]
		molecule_dict[res]=res_dict
	multi_molecule_dict[mol] = molecule_dict	


#----  Bead detection -----------------------------------#
#chains -> bead type -> nucleotide -> capture
lammps_bead_data=[]
bead_counter=0
for i in multi_molecule_dict.keys():
	for j in multi_molecule_dict[i].keys():
		for k in multi_molecule_dict[i][j].keys():
			bead_counter=bead_counter+1
			lammps_bead_data.append([bead_counter,i]+multi_molecule_dict[i][j][k][1:8])
			
			

#----------------------------------------------------------------------------------------#

#---- Bond detection -----------------------------------#
lammps_bond_data=[]
bond_counter=0
for i in mols:
	res_list = list(multi_molecule_dict[i].keys())
	for j in range(len(res_list)):
		c_bead = list(multi_molecule_dict[i][res_list[j]].keys())
		c_bead_ca, c_bead_cb,n_bead_ca = [], [], []
		c_bead_ca_xyz, c_bead_cb_xyz, n_bead_ca_xyz = [], [], []
		for k in c_bead:
			if "CA" in k:
				c_bead_ca = multi_molecule_dict[i][res_list[j]][k][0]
				c_bead_ca_xyz = multi_molecule_dict[i][res_list[j]][k][5:8]
				

			if "CB" in k:
				c_bead_cb = multi_molecule_dict[i][res_list[j]][k][0]
				c_bead_cb_xyz = multi_molecule_dict[i][res_list[j]][k][5:8]

		if(j!=(len(res_list)-1)):
			n_bead = list(multi_molecule_dict[i][res_list[j+1]].keys())	
			for k in n_bead:
				if "CA" in k:
					n_bead_ca =  multi_molecule_dict[i][res_list[j+1]][k][0]
					n_bead_ca_xyz = multi_molecule_dict[i][res_list[j+1]][k][5:8]
					
			if(len(c_bead_cb)!=0):
				bond_counter=bond_counter+1
				dist = euclidean_distance(c_bead_ca_xyz,c_bead_cb_xyz)
				lammps_bond_data.append([bond_counter,bond_counter,c_bead_ca, c_bead_cb, dist])
		

			bond_counter=bond_counter+1
			dist = euclidean_distance(c_bead_ca_xyz,n_bead_ca_xyz)
			lammps_bond_data.append([bond_counter, bond_counter,c_bead_ca, n_bead_ca, dist])
		else:
			if(len(c_bead_cb)!=0):
				bond_counter=bond_counter+1
				dist = euclidean_distance(c_bead_ca_xyz,c_bead_cb_xyz)
				lammps_bond_data.append([bond_counter, bond_counter,c_bead_ca, c_bead_cb, dist])				
	
#-------------------------------------Writing the topology file-------------------------#
print("#LAMMPS datafile",file=data_file)
print(str(bead_counter)+" atoms",file=data_file)
print(str(bond_counter)+" bonds",file=data_file)
print(file=data_file)
print(str(nbeads)+" atom types",file=data_file)
print(str(bond_counter)+" bond types",file=data_file)
print(file=data_file)
box_dim_str = str(-round(box_dimensions[0]/2.0,1))+" "+str(round(box_dimensions[0]/2.0,1))
print(box_dim_str+" xlo xhi",file=data_file)
print(box_dim_str+" ylo yhi",file=data_file)
print(box_dim_str+" zlo zhi",file=data_file)
print("\nAtoms\n",file=data_file)
for i in lammps_bead_data:
    print(i[0]," ",i[1]," ", i[2], " ", i[3], " ",i[4], " ", i[5], " ",i[6], " ",i[7]," ",i[8],file=data_file)

print("\nBonds\n",file=data_file)
for i in lammps_bond_data:
    print(i[0],i[1],i[2],i[3],sep="\t",file=data_file)


print("\nBond Coeffs\n",file=data_file)
for i in lammps_bond_data:
    b1 = int(i[2])-1
    b2 = int(i[3])-1
    if( (domain_id[b1] != 0) and (domain_id[b2] != 0) and (domain_id[b1] == domain_id[b2]) ):
        if( charge[b1]+charge[b2] != 0 ):
            spring_constant = spring_constant_folded
        else:
            spring_constant = spring_constant_idp
    else:            
        spring_constant = spring_constant_idp
    
    print(i[0],spring_constant,R0,i[4],file=data_file)
data_file.close()

