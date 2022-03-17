from gurobipy import GRB, quicksum
from top_callbk import dfjcut_callback
from cycflw import Cycflw
from ximblnce import Ximblnce
from utils import sol2stage, make_mheur
import csv
import os

class CycflwRun(Cycflw):
    """
    Class to solve the Cycflw base class model without adding any extras.
    """ 
    def __init__(self, args, RESLIM=None, OPTCR=0, MIPSTART=False, MIPHEUR = None):    
        
        mipstart = MIPSTART
        mip_heur = MIPHEUR        
        Cycflw.__init__(self, args, RESLIM, OPTCR, MIPSTART=mipstart, MIPHEUR = mip_heur)
        if self.model._MipHeur['enabled']:
            self.model = make_mheur(self.model)
            self.model._mheur._ztruu = True
            self.model._mheur._kvlfl = True

    def solve(self):
        self = sol2stage(self, dfjcut_callback)
        ##self.model.Params.MIPFocus = 1
        self = Cycflw.solve(self, dfjcut_callback)
        self.totCUT = 0
        for i,j in self.uarcs:
            y_ij = self.model._y[i,j].x
            self.totCUT += self.flow0[(i,j)]*round(y_ij)
        print("TOTAL RUNTIME IS: " + str(self.model.Runtime))
        return self    