from gurobipy import GRB, quicksum, tupledict
from cycflw import Cycflw
import csv
import os

def root_arcs(islsw,nodes):
    """
    Return fictitious k arcs between root nodes of generator groups (k is the
    total number of groups, k=len(islsw)) and the fictitious "superroot" sr.
    """
    sr = max(nodes)+1;
    roots = islsw.copy()
    uarcs = [(sr,root) for root in roots]
    darcs = uarcs + [(root,sr) for root in roots]
    return sr, uarcs, darcs

class CycflwCohen(Cycflw):
    """
    Class to solve the Cycflw base class model by adding graph connectivity  
    as formulated by Natann Cohen. The optimization objective is in cycflw.py.
    
    If no cycle subdivision is desired, simply modify the class to inherit from 
    Cycflw (change three class calls and one import).
    """ 
    def __init__(self, args, RESLIM=None, OPTCR=0, MIPSTART=False):
        
        mipstart = MIPSTART
        Cycflw.__init__(self, args, RESLIM, OPTCR, MIPSTART=mipstart)
        
        # Create the graph partitioning model (initial solve)
        modOPF = self.model
        modOPF.update()

        # Edge flow constraints
        islsw = modOPF._islsw.copy()
        sr, fuarcs, fdarcs = root_arcs(islsw,self.nodes)    
        f_arcs = fdarcs + list(self.arcs)
        rcprcl = 1/(len(self.nodes)+1);  # all nodes plus the "superroot" sr
        modOPF._f = modOPF.addVars( f_arcs, name="f", vtype=GRB.CONTINUOUS, lb=0, ub=(1-rcprcl))
        modOPF.addConstr(quicksum(modOPF._f[i,j] for i,j in f_arcs )==len(self.nodes), "sptreeF1")     
        EgFLW = tupledict()        
        for i,j in fuarcs:
            EgFLW[i,j] = modOPF.addConstr(modOPF._f[i,j] + modOPF._f[j,i] == 1, "EgFLOW[%s,%s]" % (i,j))
            modOPF._f[i,j].lb = rcprcl
            modOPF._f[j,i].lb = rcprcl
        for i,j in self.uarcs:
            EgFLW[i,j] = modOPF.addConstr(modOPF._f[i,j] + modOPF._f[j,i] == modOPF._z[i,j] + modOPF._z[j,i], "EgFLOW[%s,%s]" % (i,j))            
            modOPF.addConstr(modOPF._f[i,j] >= rcprcl*(modOPF._z[i,j] + modOPF._z[j,i]), "FLOWlo1[%s,%s]" % (i,j))            
            modOPF.addConstr(modOPF._f[j,i] >= rcprcl*(modOPF._z[i,j] + modOPF._z[j,i]), "FLOWlo1[%s,%s]" % (j,i))                       
            modOPF.addConstr(modOPF._f[i,j] <= (1-rcprcl)*(modOPF._z[i,j] + modOPF._z[j,i]), "FLOWHi1[%s,%s]" % (i,j)) 
            modOPF.addConstr(modOPF._f[j,i] <= (1-rcprcl)*(modOPF._z[i,j] + modOPF._z[j,i]), "FLOWHi1[%s,%s]" % (j,i)) 
        VxFLW = tupledict() 
        Nv = islsw  # neighbors of sr
        VxFLW[sr] = modOPF.addConstr(quicksum(modOPF._f[j,sr] for j in Nv) == 1-rcprcl, "VxFLOW[%s]" % (sr))               
        nonroot = set(self.nodes) - set(islsw)
        for v in islsw:
            Nv = set(modOPF._graph.neighbors(v))
            Nv.add(sr)
            VxFLW[v] = modOPF.addConstr(quicksum(modOPF._f[j,v] for j in Nv) == 1-rcprcl, "VxFLOW[%s]" % (v))
        for v in nonroot:
            Nv = set(modOPF._graph.neighbors(v))
            VxFLW[v] = modOPF.addConstr(quicksum(modOPF._f[j,v] for j in Nv) == 1-rcprcl, "VxFLOW[%s]" % (v))       
        
        # The objective is in cycflw.py. This function only modifies the power 
        # flow constraints. 
        modOPF._ztruu = False    # callback flags
        modOPF._kvlfl = True     # callback flags

        
    def solve(self, callbackfcn=None):
        self = Cycflw.solve(self, None) 
        y = list()     #branch status
        self.totCUT = 0
        for i,j in self.uarcs:
            y_ij = self.model._y[i,j].x
            self.totCUT += self.flow0[(i,j)]*round(y_ij)
            if y_ij<1e-2: y.append([self.model._i2n[i], self.model._i2n[j], 1])
            else: y.append([self.model._i2n[i], self.model._i2n[j], 0])  
        print("TOTAL RUNTIME IS: " + str(self.model.Runtime))
        ici_pth = os.path.abspath(os.path.join(os.path.dirname( __file__ ),'ici'))         
        with open(os.path.join(ici_pth,"y.csv"),"w",newline='') as f:
            csv_out=csv.writer(f)
            csv_out.writerows(y)
        return self  

