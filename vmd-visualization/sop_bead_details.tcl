# This .tcl script allows the user to visualize SOP CG beads with their real charges, massses and VdW radii
# Author : Krishnakanth B, Theoretical Biophysics Laboratory, 
# Molecular Biophysics Unit, Indian Institute of Science, Bangalore - 560012
# Please feel free to contact me at krishnakanth.baratam@gmail.com for any issues
# Usage
# Open VMD's TK Console and source this .tcl script with full path
# In VMD's TK console
# source /home/krishnakanth/vmd-visualization/sop_bead_details.tcl


#1,C_alpha ,1.9,0,56.04298
set sel [atomselect top "type 1"]
$sel set radius 1.9
$sel set charge 0.0
$sel set mass 56.04298


#2,Gly,2.25,0,57.051
set sel [atomselect top "type 2"]
$sel set radius 2.25
$sel set charge 0.0
$sel set mass 57.051

#3,Ala,2.52,0,15.03452
set sel [atomselect top "type 3"]
$sel set radius 2.52
$sel set charge 0.0
$sel set mass 15.03452

#4,Arg,3.28,1,100.1423
set sel [atomselect top "type 4"]
$sel set radius 3.28
$sel set charge 1.0
$sel set mass 100.1423

#5,Lys,3.18,1,72.1289
set sel [atomselect top "type 5"]
$sel set radius 3.18
$sel set charge 1.0
$sel set mass 72.1289

#6,His,3.04,0/1,82.10384
set sel [atomselect top "type 6"] 
$sel set radius 3.18
$sel set charge 0.0
$sel set mass 82.10384

#7,Asp,2.79,-1,59.04322
set sel [atomselect top "type 7"] 
$sel set radius 2.79
$sel set charge -1.0
$sel set mass 59.04322


#8,Glu,2.96,-1,73.0698
set sel [atomselect top "type 8"] 
$sel set radius 2.96
$sel set charge 0.0
$sel set mass 73.0698

#9,Ser,2.59,0,31.03352
set sel [atomselect top "type 9"] 
$sel set radius 2.59
$sel set charge 0.0
$sel set mass 31.03352


#10,Thr,2.81,0,45.0601
set sel [atomselect top "type 10"] 
$sel set radius 2.81
$sel set charge 0.0
$sel set mass 45.0601

#11,Asn,2.84,0,58.05886
set sel [atomselect top "type 11"] 
$sel set radius 2.84
$sel set charge 0.0
$sel set mass 58.05886

#12,Gln,3.01,0,72.08544
set sel [atomselect top "type 12"] 
$sel set radius 3.01
$sel set charge 0.0
$sel set mass 72.08544

#13,Cys,2.74,0,47.09952
set sel [atomselect top "type 13"] 
$sel set radius 2.74
$sel set charge 0.0
$sel set mass 47.09952

#14,Pro,2.78,0,41.0718
set sel [atomselect top "type 14"] 
$sel set radius 2.78
$sel set charge 0.0
$sel set mass 41.0718

#15,Ile,3.09,0,57.11426
set sel [atomselect top "type 15"] 
$sel set radius 3.09
$sel set charge 0.0
$sel set mass 57.11426

#16,Leu,3.09,0,57.11426
set sel [atomselect top "type 16"] 
$sel set radius 3.09
$sel set charge 0.0
$sel set mass 57.11426

#17,Met,3.09,0,75.15268
set sel [atomselect top "type 17"] 
$sel set radius 3.09
$sel set charge 0.0
$sel set mass 75.15268

#18,Phe,3.18,0,91.13048
set sel [atomselect top "type 18"] 
$sel set radius 3.18
$sel set charge 0.0
$sel set mass 91.13048

#19,Trp,3.39,0,130.16652
set sel [atomselect top "type 19"] 
$sel set radius 3.39
$sel set charge 0.0
$sel set mass 130.16652

#20,Tyr,3.23,0,107.12948
set sel [atomselect top "type 20"] 
$sel set radius 3.23
$sel set charge 0.0
$sel set mass 107.12948

#21,Val,2.93,0,43.08768
set sel [atomselect top "type 21"] 
$sel set radius 2.93
$sel set charge 0.0
$sel set mass 43.08768
