#	Essential supporting functions required in other codes
#	Krishnakanth B, 
#	Theoretical Biophysics Laboratory, Molecular Biophysics Unit,
#	Indian Institute of Science, Bangalore - 560012
#
#	Last Modified:  28 April 2024

import numpy as np
from scipy.spatial import distance_matrix

# Euclidean distance rounded to 4 digits after decimal
def euclidean_distance(a, b):
	return round(((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)**0.5,4)



#Function to read a .pdb file, obeys 1996 nomenclature of .pdb file
def pdb_reader(pdb_filename):
	print("Reading pdb file: "+pdb_filename)
	x, y, z,res_id,atom_id ,atom_type,res_type,res_no,chain_id = [],[],[],[],[],[],[],[],[]
	f=open(pdb_filename,"r").readlines()
	for i in range(len(f)):
		if(f[i][0:4]=="ATOM"):
			x.append(float(f[i][30:38]))
			y.append(float(f[i][38:46]))
			z.append(float(f[i][46:54]))
			atom_id.append(str.strip(f[i][6:11]))
			atom_type.append(str.strip(f[i][12:16]))
			res_type.append(str.strip(f[i][17:20]))
			res_no.append(str.strip(f[i][22:26]))
			chain_id.append(f[i][21])
	return (x,y,z, atom_id, atom_type, res_type, res_no, chain_id)

# Get the unique elements of a list : Used to get the total number of chains
def unique_list(my_list):
	uniq_list=[]
	for i in  my_list:
		if not (i  in uniq_list):
			uniq_list.append(i)
	return uniq_list

# which like function (R language) in python
def which(x,my_list):
	req=[]
	for i in range(len(my_list)):
		if(x==my_list[i]):
			req.append(i)
	return req
	
# This function takes in an np array of atom types and cleans the H atom types starting with numeric value
def clean_H(y):
	for i in range(len(y)):
		if(str.isnumeric(y[i][0])):
			y[i] = y[i][1:]+y[i][0]
#	return(y)
	
# Function to read a large file in chunks
def rows(f, chunksize=1024, sep='|'):
    """
    Read a file where the row separator is '|' lazily.

    Usage:

    >>> with open('big.csv') as f:
    >>>     for r in rows(f):
    >>>         process(r)
    """
    row = ''
    while (chunk := f.read(chunksize)) != '':   # End of file
        while (i := chunk.find(sep)) != -1:     # No separator found
            yield row + chunk[:i]
            chunk = chunk[i+1:]
            row = ''
        row += chunk
    yield row
	
	
def rg_ree_calc(x,y,z, atom_type, mass_dict):
	rg = 0.0
	cmx = 0.0
	cmy = 0.0
	cmz = 0.0
	t_mass = 0.0
	Ixx, Ixy, Ixz, Iyy, Iyz, Izz = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	for i in range(len(x)):
		cmx = cmx + x[i]*mass_dict[atom_type[i]]
		cmy = cmy + y[i]*mass_dict[atom_type[i]]
		cmz = cmz + z[i]*mass_dict[atom_type[i]]
		t_mass = t_mass + mass_dict[atom_type[i]]
	cmx = cmx/t_mass 
	cmy = cmy/t_mass
	cmz = cmz/t_mass
	for i in range(len(x)):
		Ixx = Ixx + (x[i]-cmx)*(x[i]-cmx)
		Ixy = Ixy + (x[i]-cmx)*(y[i]-cmy)
		Ixz = Ixz + (x[i]-cmx)*(z[i]-cmz)
		Iyy = Iyy + (y[i]-cmy)*(y[i]-cmy)
		Iyz = Iyz + (y[i]-cmy)*(z[i]-cmz)
		Izz = Izz + (z[i]-cmz)*(z[i]-cmz)
		rg = rg + mass_dict[atom_type[i]]*((x[i]-cmx)**2 + (y[i]-cmy)**2 + (z[i]-cmz)**2)
	Ixx = Ixx/len(x)
	Ixy = Ixy/len(x)
	Ixz = Ixz/len(x)
	Iyy = Iyy/len(x)
	Iyz = Iyz/len(x)
	Izz = Izz/len(x)
	gyration_tensor  = np.array([[Ixx,Ixy,Ixz],[Ixy,Iyy,Iyz],[Ixz,Iyz,Izz]])
	evals, evecs = np.linalg.eig(gyration_tensor)
	L1 = evals[0]
	L2 = evals[1]
	L3 = evals[2]
	asphr = ((L1 - L2)**2 + (L2 - L3)**2 + (L1 - L3)**2)/(2.0*((L1 + L2 + L3)**2))
	L_mean = (L1+L2+L3)/3
	shape = ((L1-L_mean)*(L2-L_mean)*(L3-L_mean))/(L_mean**3)
	rg = np.sqrt(rg/t_mass)
	ree = np.sqrt(((x[0]-x[-1])**2) + ((y[0]-y[-1])**2) + ((z[0]-z[-1])**2))
	req = np.transpose(np.vstack([x,y,z]))
	req = distance_matrix(req, req)
	req = 1/req[np.triu_indices(len(x), k=1)]
	rh = 1.0/(np.mean(req))
	return([rg/10.0, ree/10.0,rh/10.0, asphr,shape])




def flory_exponent_calc(x,y,z, atom_type):
	sep_vals = [0]
	CA_indices = np.where((atom_type=='1') | (atom_type=='2') )
	x_mod = x[CA_indices]
	y_mod = y[CA_indices]
	z_mod = z[CA_indices]

	req = np.transpose(np.vstack([x_mod,y_mod,z_mod]))
	req = distance_matrix(req, req)
	seq_len = len(x_mod)
	for i in range(1,seq_len):
		m = []
		for j in range(seq_len-i):
			m.append(req[j,j+i])
		sep_vals.append(np.mean(m)/10.0)
#	print(sep_vals)
	return(np.array(sep_vals))

def pdb_dumper(x,y,z, atom_type, aa_dict, pdb_num, stride, pdb_start):
	if((pdb_num > pdb_start)  and (pdb_num%stride == 0)):
		pdb_file = open("./pdb_files/"+str(pdb_num)+".pdb","w")
		res_no = 0
		atom_counter = 1
		for i in range(len(x)):
			if(atom_type[i]=="1"):
				res_no = res_no + 1
				res_type = "CA"
				res_name = aa_dict[atom_type[i+1]]
			else:
				res_name = aa_dict[atom_type[i]]
				res_type = "CB"
				#res_no = res_no + 1
			
			if(atom_type[i]=="2"):
				res_type = "CA"
				res_no = res_no + 1
			
			ps = "ATOM  "+"{0:>5}".format(atom_counter)+"{0:>4}".format(res_type)+ \
				"  "+res_name+" A"+"{0:>4}".format(res_no)+ \
				"    "+"{0:>8}".format(np.round(x[i],3))+"{0:>8}".format(np.round(y[i],3))+ \
				"{0:>8}".format(np.round(z[i],3))
			atom_counter = atom_counter+1
			print(ps, file =pdb_file)
		pdb_file.close()
	else:
		return(0)
		
		
	
	
def rg_calc(x,y,z, atom_type, mass_dict):
	rg = 0.0
	cmx = 0.0
	cmy = 0.0
	cmz = 0.0
	t_mass = 0.0
	for i in range(len(x)):
		cmx = cmx + x[i]*mass_dict[atom_type[i]]
		cmy = cmy + y[i]*mass_dict[atom_type[i]]
		cmz = cmz + z[i]*mass_dict[atom_type[i]]
		t_mass = t_mass + mass_dict[atom_type[i]]
	cmx = cmx/t_mass 
	cmy = cmy/t_mass
	cmz = cmz/t_mass
	for i in range(len(x)):
		rg = rg + mass_dict[atom_type[i]]*((x[i]-cmx)**2 + (y[i]-cmy)**2 + (z[i]-cmz)**2)

	rg = np.sqrt(rg/t_mass)
	return(rg/10.0)
	

def contact_map(x,y,z, atom_type):
	ca_indices = np.where((atom_type=="1") | (atom_type=="2"))[0]
	x = x[ca_indices]
	y = y[ca_indices]
	z = z[ca_indices]
	req = np.transpose(np.vstack([x,y,z]))
	req = distance_matrix(req, req)
	return(req)
		
		


		
#def local_rg(pdb_fname,mass_dict window_size):
#	tmp = []
#	x,y,z, atom_id, atom_type, res_type, res_no, chain_id =pdb_reader(cg_filename)
#	n_residues = 
#	for i in range(window_size - 1, n_residues):
#		m_begin = i - (window_size-1)
#		m_end = i+1
#	tmp.append(rg_calc(x[m_begin:m_end],y[m_begin:m_end],z[m_begin:m_end],atom_type[m_begin:m_end],mass_dict))
#	print(tmp)
#	return(tmp)
	
	
def rg_ree_calc_quick(x,y,z, atom_type, mass_dict):
	rg = 0.0
	cmx = 0.0
	cmy = 0.0
	cmz = 0.0
	t_mass = 0.0
	for i in range(len(x)):
		cmx = cmx + x[i]*mass_dict[atom_type[i]]
		cmy = cmy + y[i]*mass_dict[atom_type[i]]
		cmz = cmz + z[i]*mass_dict[atom_type[i]]
		t_mass = t_mass + mass_dict[atom_type[i]]
	cmx = cmx/t_mass 
	cmy = cmy/t_mass
	cmz = cmz/t_mass
	for i in range(len(x)):
		rg = rg + mass_dict[atom_type[i]]*((x[i]-cmx)**2 + (y[i]-cmy)**2 + (z[i]-cmz)**2)
	rg = np.sqrt(rg/t_mass)
	ree = np.sqrt(((x[0]-x[-1])**2) + ((y[0]-y[-1])**2) + ((z[0]-z[-1])**2))
	print(str(rg/10.0)+","+str(ree/10.0))
