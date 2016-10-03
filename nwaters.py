"""
Script to calculate number of water around a com 
"""
from __future__ import division, print_function
import math
import numpy as np
from  pylbtc.vmd import *
from pylbtc.vmd.Wrap3 import wrap
import path_sys
from macros import macro_docking_ncp, macro_pep

def calc_wat(NAME, FIRST, LAST, SKIP, SELECTIONS, DISTANCES, CENTER_SEL, OUT):

    path = path_sys.path(NAME)
# Loading molecules
    PSF = NAME+".0.psf"
    PDB = NAME+".0.pdb"
    DCD = NAME+".{}.dcd"

# Loading Molecule
    mol = molecule.load("psf", path+PSF)
    for dcd in range(FIRST, LAST+1):
        molecule.read(mol, "dcd", path+DCD.format(dcd), waitfor=-1, skip=SKIP)
    nframes = molecule.numframes(mol)

# Wrapping
    #evaltcl('package require pbctools; pbc wrap -centersel "{0}" -center com -compound fragment -all'.format(CENTER_SEL))
    wrap(CENTER_SEL, mol)

# Centralizing frames
    center_sel = atomsel(CENTER_SEL, mol)
    for frame in range(nframes):
        center_sel.frame = frame
        vmdnumpy.timestep(mol, frame)[:] += np.array(center_sel.center())*-1

    for selection in SELECTIONS:
            # Moving dummy atom into place
            sel = atomsel(selection, mol)
            #dummyAtom = atomsel(DUMMY, mol)
            # To get all water in distance delete ("+rad+" >= {3} and ) and ( DISTANCES[i-1]**2,)
            waters = np.zeros((nframes, len(DISTANCES)-1), dtype="float")
            rad = "((x-{0})*(x-{0}) + (y-{1})*(y-{1}) + (z-{2})*(z-{2}))"
            sel_str = "name OH2 and ("+rad+"< {4})"
            for frame in range(nframes):
                    sel.frame = frame
                    com = sel.center(sel.get("mass"))
                    for i in range(1, len(DISTANCES)):
                            waters[frame, i-1] = len(atomsel(sel_str.format(com[0], com[1], com[2], DISTANCES[i-1]**2, DISTANCES[i]**2), mol, frame))

            data = np.zeros((len(DISTANCES)-1, 3), dtype="float")
            #data[:,0] = np.array(DISTANCES[1:], dtype="float") - 0.5
            data[:,0] = DISTANCES[1:]
            data[:,1] = np.mean(waters, axis=0)
            data[:,2] = np.std(waters, axis=0)
            np.savetxt("{}{}.dcd{}-{}_skip{}.nwaters.dat".format(OUT, selection, FIRST, LAST, SKIP), waters)
            np.savetxt("{}{}.dcd{}-{}_skip{}.av_nwaters.dat".format(OUT, selection, FIRST, LAST, SKIP), data)

