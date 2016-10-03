import numpy as np
from pylbtc.vmd import System
from pylbtc.vmd.Timeseries import BaseTimeseries

class RmsdTimeseries(BaseTimeseries):
    # TODO
    """
    documentation
    """

    def __init__(self, fit_sel, fit_sel_ref, system_fit, system_fit_ref, rmsd_sel=None, rmsd_sel_ref=None, system_rmsd=None, system_rmsd_ref=None):
        super().__init__((0, 2))    
        self.__system_fit = system_fit
        self.__system_fit_ref = system_fit_ref
        self.__fit_sel = self.__system_fit.selectAtoms(fit_sel)
        self.__fit_sel_ref = self.__system_fit_ref.selectAtoms(fit_sel_ref)
        self.__nextLine = 0

        if system_rmsd is None:
            self.__system_rmsd = system_fit
        else:
            self.__system_rmsd = system_rmsd

        if system_rmsd_ref is None:
            self.__system_rmsd_ref = system_fit_ref
        else:
            self.__system_rmsd_ref = system_rmsd_ref

        if rmsd_sel is None:
            self.__rmsd_sel = self.__system_rmsd.selectAtoms(fit_sel)
        else:
            self.__rmsd_sel = self.__system_rmsd.selectAtoms(rmsd_sel)

        if rmsd_sel_ref is None:
            self.__rmsd_sel_ref = self.__system_rmsd_ref.selectAtoms(fit_sel_ref)
        else:
            self.__rmsd_sel_ref = self.__system_rmsd_ref.selectAtoms(rmsd_sel_ref)

    def initialize(self, nframes):
        self._data.extend(nframes, 0)

    def run(self, frame):
        # Fitting
        self.__system_fit.fit(self.__fit_sel, self.__fit_sel_ref)

        # Calculate RMSD
        self.data[self.__nextLine, 0] = self.__nextLine + 1
        self.data[self.__nextLine, 1:] = self.__rmsd_sel.rmsd(self.__rmsd_sel_ref)
        self.__nextLine += 1

    def finalize(self):
        pass

    def clear(self):
        super()._clear((0, 2))
        self.__nextLine = 0

