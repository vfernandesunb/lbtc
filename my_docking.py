'''
Script that get the clouds of a docking that are near any defined site (macro)
'''

import numpy as np
import fileinput
from pylbtc.vmd import *

#TODO CHANGE THIS
from reps import rep_ncp

def clouds_in_sites(DOCK_IN, RECEPTOR, CUTOFF, CENTER, OUTPUT, SEL_RECEPTOR, SEL_LIGAND, SEL_SITES, REP, REP_RECEPTOR, LIG_COLOR, EXCLUDE_SEL=None):
# Joining receptor and dock clouds results in a single pdb file
    with open('complex.pdb', 'w') as fout, fileinput.input([DOCK_IN, RECEPTOR]) as fin:
        for line in fin:
            fout.write(line)

    mol = System("complex")
    mol.load("complex.pdb")

    fd = open(OUTPUT, "w")
    fd.write("ANALYZING IF CLOUDS FROM {} ARE AT {} DISTANCE FROM SITES {}\n\n".format(SEL_LIGAND, CUTOFF, SEL_SITES))

    # some controlers for atoms selections
    tmp = ") and not ("
    sel_none = ""
    if REP:
        mol.clearReps()

        if REP_RECEPTOR.upper() == 'NCP':
            rep_ncp.sysRep(mol)
        else:
            print("\n\n SPECIFIC REPRESENTATION FOR RECEPTOR NOT FOUND. USING STANDARD REPRESENTATION\n\n") 
            receptor = NewCartoonRep(color="colorId 2", selection="{}".format(REP_RECEPTOR))
            mol.addRep(receptor)


    for site in SEL_SITES:
        fd.write("################# SITE: {} #################\n".format(site))
        if CENTER:
            coords = mol.selectAtoms("{}".format(site)).center()
            #coords = mol.selectAtoms("{}").centerOfMass
            sel = mol.selectAtoms("{} and same residue as (sqr(x{:+.03f}) + sqr(y{:+.03f}) + sqr(z{:+.03f}) < sqr({}))".format(SEL_LIGAND, coords[0]*(-1), coords[1]*(-1), coords[2]*(-1), CUTOFF))
        else:
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

            # TODO 
            rep = VDWRep(color="colorId {}".format(LIG_COLOR), selection=sel)
            sel_none += rep.selection+tmp
            if REP:
                mol.addRep(rep)

    if REP:
        rep_none = VDWRep(color="colorId 3", selection="{} and not ({} water)".format(SEL_LIGAND, sel_none))
        mol.addRep(rep_none)

    # Getting total number of clouds
    lig = mol.selectAtoms("{}".format(SEL_LIGAND))
    lig_num = len(np.unique(lig["segname"]))
    fd.write("\n\n{} CLOUDS WERE FOUND. \n".format(lig_num))

    if EXCLUDE_SEL:
        sel_exc = mol.selectAtoms("{} and (same residue as within {} of {}) and not ({} water)".format(SEL_LIGAND, EXCLUDE_SEL[1], EXCLUDE_SEL[0], sel_none))
        exc_num = len(np.unique(sel_exc["segname"]))
        fd.write("\n{} CLOUDS WERE FOUND AT {} OF THE EXCLUDE SELECTION: {}.\n\n".format(exc_num, EXCLUDE_SEL[1], EXCLUDE_SEL[0]))
        if REP:
            rep_exc = VDWRep(color="colorId 1", selection=sel_exc)
            mol.addRep(rep_exc)

    if not REP:
        quit()





