from gurobipy import GRB, quicksum
from top_callbk import dfjcut_callback
from cycflw import Cycflw

class CycflwRun(Cycflw):
    """
    Class to solve the Cycflw base class model without adding any extras.
    """ 
    def __init__(self, args, RESLIM=None, OPTCR=0, MIPSTART=False):    
        
        mipstart = MIPSTART
        Cycflw.__init__(self, args, RESLIM, OPTCR, MIPSTART=mipstart)

    def solve(self):
        self = Cycflw.solve(self, dfjcut_callback)
        self.totCUT = 0
        for i,j in self.uarcs:
            y_ij = self.model._y[i,j].x
            self.totCUT += self.flow0[(i,j)]*round(y_ij)
        print("TOTAL RUNTIME IS: " + str(self.model.Runtime))
        return self