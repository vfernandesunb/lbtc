"""
This script calculates RMSD between given selections along a simulation.
The RMSD for the first selection in SELECTIONS will be written in the
first file name in OUTPUTS. Folder paths must end with a slash (/).
DT is the timestep in fentoseconds and NSAVED is the number of timesteps
between each frame. CONVERSION defines the time conversion from
fentoseconds to the desired time scale (1000000 = 1 nanosecond).
"""

from __future__ import print_function, division
from collections import defaultdict
import numpy as np
from pylbtc.vmd import *
import path_sys
from macros import macro_docking_ncp
from macros import macro_pep

__all__ = ("rmsd",)

def rmsd(name, first, last, sels, outs, fits, path_out='./', which_run='all', sub_folders=0, skip=1, nsaved=5000, dt=2, conversion=1000000): 

    # Initializations
    path = path_sys.path(name)
    psf = name+".0.psf"
    pdb = name+".0.pdb"
    dcds = name+".{}.dcd"

    d_compute = {}
    if which_run == "all": 
        d_compute[""] = (1, "", path_out+'all/' if sub_folders != 0 else path_out) 
        print("\n\n######WARNING: It was't choosed for which atoms will run. So calculating only for all Atoms######\n") 
    else:
        w_run = [i.lower() for i in which_run] # Let case insensitive
        #s_folders = [j.lower() for j in sub_folders] # If want to especify sub folders name
        d_compute[""] = (1, "", path_out+'all/' if sub_folders != 0 else path_out) if "all" in w_run else 0
        d_compute["_b"] = (1, "and backbone", path_out+'backbone/' if sub_folders != 0 else path_out) if "backbone" in w_run else 0
        d_compute["_noh"] = (1, "and noh", path_out+'noh/' if sub_folders != 0 else path_out) if "heavy" in w_run else 0

    # Some verifications
    if len(sels) != len(outs):
        raise IndexError("SELECTIONS and OUTPUTS must have the same length.")
    if not d_compute:
        raise IndexError("Only accepts python list with one or more of the following options: All, Backbone or Heavy!")
    #if sub_folders != 0 and len(sub_folders) != 3:
        #raise IndexError("You must specify the outputs path for all selections at which runs")
    if sub_folders == 0 and type(which_run) is list:
        print("\n\nSaving all outputs (for all which_run) in the actual folder\n\n")

    # Loading molecules
    mol = System()
    mol_ref = System()
    mol.load(path+psf)
    mol_ref.load(path+psf)
    mol_ref.load(path+pdb)

    #ref_sels = []
    #fit_sels = []
    d_selections = defaultdict(list)
    d_fit = defaultdict(list)

    # Checking which selections will run
    for key in d_compute:
        if d_compute[key]:
            for sel, out, fit in zip(sels, outs, fits):
                # Changing output name and opening files
                ext = out.split('.')[-1]
                fname = out.split('.')[0:-1]
                out = '.'.join(fname)+'{}.'+ext
                tmp_file = open(d_compute[key][2]+out.format(key), "w")

                # Making selections for fit 
                d_fit[tmp_file].append(mol.selectAtoms("({}) ".format(fit)))
                d_fit[tmp_file].append(mol_ref.selectAtoms("({})".format(fit)))

                # Making selections for rmsd 
                d_selections[tmp_file].append(mol.selectAtoms("({} {})".format(sel, d_compute[key][1])))
                d_selections[tmp_file].append(mol_ref.selectAtoms("({} {})".format(sel, d_compute[key][1])))

    # Setting time             
    step = (skip*nsaved*dt)/conversion
    time = 0

    # Main loop
    for dcd in range(first, last+1):
        mol.load(path+dcds.format(dcd), "dcd", step=skip, waitfor=-1)
        for frame in mol.trajectory:
            time += step
            for s_key, f_key in zip(d_selections, d_fit):
                rmsd = []
                fit_sel = d_fit[f_key][0]
                fit_ref = d_fit[f_key][1]                
                rmsd_sel = d_selections[s_key][0]
                rmsd_ref = d_selections[s_key][1]
                # Fitting frames
                mol.fit(fit_sel, fit_ref)
                # Computing RMSD
                rmsd = (rmsd_sel.rmsd(rmsd_ref))
                s_key.write("{} {} \n".format(time, rmsd))
        mol.delFrame()

