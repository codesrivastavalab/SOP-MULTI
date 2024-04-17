#	Force-field writer for SOP-Multi
#	Krishnakanth B, 
#	Theoretical Biophysics Laboratory, Molecular Biophysics Unit,
#	Indian Institute of Science, Bangalore - 560012
#
#	Usage: python3 ff_writer.py 25 150

#  Load the required packages
import numpy as np
import matplotlib.pylab as plt
import pandas
import sys

temperature = float(sys.argv[1]) # in celsius
sim_temp = temperature+273.15 # to kelvin conversion
print("Use temperature : "+ str(sim_temp))
salt_conc = float(sys.argv[2]) # in mM
lj_cutoff = 25.0 # angstroms
electro_cutoff = 35.0  # angstroms
excluded_volume_res_cutoff = 2 # integer
excluded_volume_epsilon = 1.0 # kcal per mole
epsilon_coul = 1.0
epsilon_rep = 1.0

# KELVIN_TO_KT = unit.AVOGADRO_CONSTANT_NA * unit.BOLTZMANN_CONSTANT_kB / unit.kilocalorie_per_mole
# = 6.02214076e23 * 1.380649e-23 / 4.184e3
# = 0.0019872042586408316
kelv_to_kt = 0.0019872042586408316
temp_unitless = kelv_to_kt*sim_temp

sim_eps = 296.0736276 - 619.2813716 * temp_unitless + 531.2826741 * temp_unitless**2 - 180.0369914 * temp_unitless**3
sim_bjerrum =  332.0637/sim_eps
kappa = np.sqrt (4*3.14159 * sim_bjerrum * 2*salt_conc*6.022e-7 / (temp_unitless ))
dielectric_idp = sim_eps


ff_file = open("forcefield.sop_multi","w")
#print("kappa in angstroms : "+str(kappa))

#pair_style hybrid/overlay lj/sop 25.0 27.0 1.0 2 debye/sop 78.0 20.0 0.127 40.0
#pair_coeff * * debye/sop 332.0637 0.127 40.0

print("T_unitless  ", temp_unitless)
#simu.epsilon = 296.0736276 - 619.2813716 *temp_unitless+ 531.2826741 * temp_unitless**2 - 180.0369914 * temp_unitless**3;
print("epsilon  ", sim_eps)
#simu.l_Bjerrum = 332.0637*unit.angstroms / simu.epsilon
print("Bjerrum length  ", sim_bjerrum)
#simu.kappa = unit.sqrt (4*3.14159 * simu.l_Bjerrum * 2*simu.Kconc*6.022e-7 / (T_unitless * unit.angstrom**3))
print("kappa   ", kappa)

hyb_str = "pair_style hybrid/overlay lj/sop "+str(lj_cutoff)
hyb_str = hyb_str + " debye/sop "+str(np.round(electro_cutoff,4))
hyb_str = hyb_str + " list/sop folded_pair_list.txt 20.0"

print("dielectric "+str(np.round(dielectric_idp,4)), file = ff_file)
print(hyb_str, file = ff_file)
print("pair_coeff * * list/sop", file = ff_file)
#print("pair_coeff * * debye/sop "+str(np.round(epsilon_coul,4))+" "+str(np.round(kappa,4)) + " " + str(np.round(electro_cutoff,4)), file=ff_file)
#print("pair_coeff * * debye/sop 332.0637 "+str(np.round(kappa,4)) + " " + str(np.round(electro_cutoff,4)), file=ff_file)
# epsilon and scaling matrices adapted from Debyan et al. 2023
epsilon = [
      # G     A      R      N      D      C      E      Q      H      I      L      K      M      F      P      S      T      W      Y      V     B   (internal unit - kJ/mol)
      #GLY   ALA    ARG    ASN    ASP    CYS    GLU    GLN    HIS    ILE    LEU    LYS    MET    PHE    PRO    SER    THR    TRP    TYR    VAL    BBone
     -0.50, -0.07,  0.35,  0.25,  0.42, -0.22,  1.19,  0.50,  0.57,  0.52,  0.35,  0.30,  0.20,  0.27, -0.03,  0.25,  0.00, -0.59, -0.10,  0.10,  0.74,  # G  GLY
     -0.07, -0.50,  0.67,  0.59,  0.74, -0.64,  1.07,  0.52,  0.52, -0.87, -0.92,  0.50, -0.57, -0.82,  0.17,  0.37,  0.00, -0.99, -0.37, -0.94,  0.74,  # A  ALA
      0.35,  0.67,  0.32,  0.05, -1.76,  0.79, -1.86, -0.30,  0.10,  0.45,  0.22,  1.24,  0.42,  0.20, -0.05,  0.30,  0.00, -1.02, -0.92,  0.42,  0.74,  # R  ARG
      0.25,  0.59,  0.05, -0.10, -0.30,  0.69, -0.02, -0.12,  0.25,  1.36,  0.89, -0.35,  0.79,  0.72,  0.32,  0.35,  0.00, -0.22,  0.02,  0.97,  0.74,  # N  ASN
      0.42,  0.74, -1.76, -0.30,  0.67,  0.94,  0.99,  0.30, -0.55,  1.34,  1.54, -1.71,  1.54,  1.19,  0.62,  0.02,  0.00,  0.15, -0.17,  1.64,  0.74,  # D  ASP
     -0.22, -0.64,  0.79,  0.69,  0.94, -3.32,  1.14,  0.10, -0.47, -1.19, -1.24,  0.87, -1.21, -1.31, -0.45,  0.22,  0.00, -1.83, -0.40, -1.26,  0.74,  # C  CYS
      1.19,  1.07, -1.86, -0.02,  0.99,  1.14,  1.12,  0.25, -0.27,  0.94,  0.92, -2.16,  0.59,  0.84,  0.64,  0.25,  0.00, -0.37, -0.40,  1.02,  0.74,  # E  GLU
      0.50,  0.52, -0.30, -0.12,  0.30,  0.10,  0.25,  0.35,  0.55,  0.35,  0.20, -0.50, -0.02, -0.10, -0.12,  0.62,  0.00, -0.27, -0.45,  0.42,  0.74,  # Q  GLN
      0.57,  0.52,  0.10,  0.25, -0.55, -0.47, -0.27,  0.55, -0.82,  0.47,  0.25,  0.64, -0.42, -0.47, -0.12,  0.37,  0.00, -1.14, -0.52,  0.45,  0.74,  # H  HIS
      0.52, -0.87,  0.45,  1.36,  1.34, -1.19,  0.94,  0.35,  0.47, -1.49, -1.96,  0.52, -1.49, -1.61,  0.12,  0.87,  0.00, -1.61, -0.82, -1.69,  0.74,  # I  ILE
      0.35, -0.92,  0.22,  0.89,  1.54, -1.24,  0.92,  0.20,  0.25, -1.96, -2.01,  0.40, -1.69, -1.93, -0.20,  0.64,  0.00, -1.74, -1.09, -1.98,  0.74,  # L  LEU
      0.30,  0.50,  1.24, -0.35, -1.71,  0.87, -2.16, -0.50,  0.64,  0.52,  0.40,  0.94,  0.55,  0.27,  0.30,  0.25,  0.00, -0.69, -0.99,  0.40,  0.74,  # K  LYS
      0.20, -0.57,  0.42,  0.79,  1.54, -1.21,  0.59, -0.02, -0.42, -1.49, -1.69,  0.55, -1.39, -2.21, -0.40,  0.79,  0.00, -2.33, -1.26, -1.17,  0.74,  # M  MET
      0.27, -0.82,  0.20,  0.72,  1.19, -1.31,  0.84, -0.10, -0.47, -1.61, -1.93,  0.27, -2.21, -2.03, -0.47,  0.25,  0.00, -1.93, -1.21, -1.66,  0.74,  # F  PHE
     -0.03,  0.17, -0.05,  0.32,  0.62, -0.45,  0.64, -0.12, -0.12,  0.12, -0.20,  0.30, -0.40, -0.47, -0.17,  0.42,  0.00, -1.81, -0.99, -0.20,  0.74,  # P  PRO
      0.25,  0.37,  0.30,  0.35,  0.02,  0.22,  0.25,  0.62,  0.37,  0.87,  0.64,  0.25,  0.79,  0.25,  0.42,  0.32,  0.00,  0.17,  0.17,  0.62,  0.74,  # S  SER
      0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.74,  # T  THR
     -0.59, -0.99, -1.02, -0.22,  0.15, -1.83, -0.37, -0.27, -1.14, -1.61, -1.74, -0.69, -2.33, -1.93, -1.81,  0.17,  0.00, -1.83, -1.36, -1.54,  0.74,  # W  TRP
     -0.10, -0.37, -0.92,  0.02, -0.17, -0.40, -0.40, -0.45, -0.52, -0.82, -1.09, -0.99, -1.26, -1.21, -0.99,  0.17,  0.00, -1.36, -0.67, -0.67,  0.74,  # Y  TYR
      0.10, -0.94,  0.42,  0.97,  1.64, -1.26,  1.02,  0.42,  0.45, -1.69, -1.98,  0.40, -1.17, -1.66, -0.20,  0.62,  0.00, -1.54, -0.67, -1.78,  0.74,  # V  VAL
      0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74]  # B  BBone
      
scaling = [0.4352, 0.5731,
               0.5731, 0.3035]
epsilon = np.array(epsilon)               
epsilon = epsilon.reshape((21,21))

mapping_vector = np.array([2,3,4,11,7,13,8,12,6,15,16,5,17,18,14,9,10,19,20,21,1])

prefac = [0.4352,0.5731,0.3035]
# Read bead details
bead_details = pandas.read_csv("sop_bead_details.csv")

#print(bead_details)

for i in range(len(bead_details)):
	for j in range(len(bead_details)):
		if(i<=j):
			c1 = bead_details["charge_e"][i]
			c2 = bead_details["charge_e"][j]
			if ((c1!="0") and (c2!="0")):
				b1 = bead_details["LAMMPS_bead_id"][i]
				b2 = bead_details["LAMMPS_bead_id"][j]
				print("pair_coeff "+str(b1)+" "+str(b2)+" debye/sop "
				+str(np.round(epsilon_coul,4))
				+" "+str(np.round(kappa,4)) + " " + str(np.round(electro_cutoff,4)), file=ff_file)

			
			





for i in bead_details["LAMMPS_bead_id"]:
	for j in bead_details["LAMMPS_bead_id"]:
		if(i<=j):
			p_i = np.where(mapping_vector==i)[0][0]
			p_j = np.where(mapping_vector==j)[0][0]
			if(i == 1 and j == 1):
				epsilon_lj = prefac[0]*np.absolute(epsilon[p_i,p_j]-1.74)/4.184 # SOP-IDP
			
			if(i == 1 and j != 1):
				epsilon_lj = prefac[1]*np.absolute(epsilon[p_i,p_j]-1.74)/4.184 # SOP-IDP

			if(i != 1 and j != 1):
				epsilon_lj = prefac[2]*np.absolute(epsilon[p_i,p_j]-1.74)/4.184 # SOP-IDP
			
			v1 = np.float64(bead_details["vdW_ang"][bead_details["LAMMPS_bead_id"]==i])
			v2 = np.float64(bead_details["vdW_ang"][bead_details["LAMMPS_bead_id"]==j])
			
			A = np.round(epsilon_rep,4)
			B = np.round(epsilon_lj,5)
			C = np.round(v1+v2,5)
			print("pair_coeff "+str(i)+" "+str(j)+" lj/sop "+str(A)+" "+str(B)+" "+str(C), file =ff_file)

		
ff_file.close()
