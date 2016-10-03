'''
Script that get the clouds of a docking that are near any defined site (macro)
'''

import numpy as np
import fileinput
from pylbtc.vmd import *
from macros import macro_ncp, macro_docking_ncp
from reps import rep_ncp

# Main header
DOCK_IN = "centers.pdb"
RECEPTOR = "/home/vfernandes/Simulations/chromatin/analysis/docking/ncp/with_DNA.dcd_22-27/frames/av.pdb"
COMPLEX = "complex.pdb"
CUTOFF = 3
OUTPUT = "check_sites.dat"

#Selections header
SEL_RECEPTOR = 'ncp'
SEL_LIGAND = 'resname CHL'
SEL_SITES = ['site2', 'site2s', 'site3', 'site4', 'site5s', 'site3s', 'site3s_6sn']

#Representation header
REP = True
LIG_COLOR = 2 # ID from VMD, can't use White (8)

# Joining receptor and dock clouds results in a single pdb file
with open('{}'.format(COMPLEX), 'w') as fout, fileinput.input([DOCK_IN, RECEPTOR]) as fin:
    for line in fin:
        fout.write(line)

mol = System(COMPLEX)
mol.load(COMPLEX)

fd = open(OUTPUT, "w")
fd.write("ANALYZING IF CLOUDS FROM {} ARE AT {} DISTANCE FROM SITES {}\n\n".format(SEL_LIGAND, CUTOFF, SEL_SITES))

if REP:
    mol.clearReps()
    rep_ncp.sysRep(mol)
    tmp = ") and not ("
    sel_none = ""

for site in SEL_SITES:
    fd.write("################# SITE: {} #################\n".format(site))
    sel = mol.selectAtoms("{} and same residue as within {} of {}".format(SEL_LIGAND, CUTOFF, site))
    if len(sel) == 0:
        print("\n\nWARNING!!! NO CLOUDS WERE FOUND IN THE SITE {} AT **{}A** CUTOFF DISTANCE\n\n".format(site.upper(), CUTOFF))
        fd.write("NO CLOUDS WERE FOUND IN THE SITE {} AT **{}A** CUTOFF DISTANCE\n\n".format(site.upper(), CUTOFF))
    else:
        dock_seg = np.unique(sel["segname"]).tolist() # Getting the cloud segname
        for d_seg in dock_seg:
            fd.write("Docking cloud segname: {}\n".format(d_seg))
            d_num = mol.selectAtoms('segname {}'.format(d_seg))["occupancy"][0]
            fd.write("Docking cloud number: {}\n".format(d_num))
            d_energy = mol.selectAtoms('segname {}'.format(d_seg))["beta"][0]
            fd.write("Docking cloud energy: {}\n\n".format(d_energy))

        if REP:
            rep = VDWRep(color="colorId {}".format(LIG_COLOR), selection=sel)
            mol.addRep(rep)
            sel_none += rep.selection+tmp

if REP:
    rep_none = VDWRep(color="colorId 8", selection="{} and not({} water)".format(SEL_LIGAND, sel_none))
    mol.addRep(rep_none)



if not REP:
    quit()





