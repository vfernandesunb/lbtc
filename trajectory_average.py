"""
"""

from __future__ import print_function, division
import numpy as np
from pylbtc.vmd import *
from pylbtc.vmd.Wrap3 import wrap
import path_sys
from macros import macro_docking_ncp, macro_pep


def gen_av(NAME, FIRST, LAST, SKIP=1, FIT_SEL='backbone'):

    path = path_sys.path(NAME)
# Loading molecules
    PSF = NAME+".0.psf"
    PDB = NAME+".0.pdb"
    DCD = NAME+".{}.dcd"

# Calculating the average 
            # Loading structures
    mol = molecule.load("psf", path+PSF)
    av = molecule.load("psf", path+PSF)
    molecule.read(av, "dcd", path+DCD.format(FIRST), waitfor=-1, skip=100)

# Initializations
    all_sel = atomsel("all", mol)
    ref_sel = atomsel(FIT_SEL, mol, frame=0)
    ref_sel = atomsel(FIT_SEL, av, frame=0)
    fit_sel = atomsel(FIT_SEL, mol)

# Only centralize first frame, from the reference, the others will be fitted
    vmdnumpy.timestep(av, 0)[:] += np.array(ref_sel.center()) * -1

    av_coords = np.zeros([molecule.numatoms(mol), 3])

    n_frames = 0

    for dcd in range(FIRST, LAST+1):
        molecule.read(mol, "dcd", path+DCD.format(dcd), waitfor=-1, skip=SKIP)
        wrap(FIT_SEL, mol)
        for frame in range(molecule.numframes(mol)):
            fit_sel.frame = all_sel.frame = frame
            all_sel.move(fit_sel.fit(ref_sel))
            av_coords += vmdnumpy.timestep(mol, frame)
            n_frames += 1
        molecule.delframe(mol)       

# Duplicate the frame to avoid losing any frame data    
    molecule.dupframe(mol, -1)
#Making average
    vmdnumpy.timestep(mol, -1)[:] = av_coords / n_frames
    all_sel.write("psf", "av.{}.psf".format(NAME))
    all_sel.write("pdb", "av.{}.pdb".format(NAME))

