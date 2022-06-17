"""
This file contains second-level callback functions to be called from 
within the other callback functions (e.g., for the optimization of 
auxiliary MIP models used in MIP heuristics).  
"""

from gurobipy import GRB
from cut_callbk import zijcut
from cut_callbk import kvlcut
import os                    #for dbg
import sys                   #for dbg
rootpth = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..\\..\\'))
sys.path.append(rootpth)     #for dbg
import igraphmod             #for dbg

#@profile
def dfjhrs_callback(model, where):
    fl_zij = model._ztruu
    fl_kvl = model._kvlfl
    if where == GRB.Callback.MIPSOL:
        y_val = model.cbGetSolution(model._y)
        g_yyy = model._ntwrk.copy()
        e2i   = model._e2i
        a_del = [e2i[ij] for ij,w in y_val.items() if abs(w)>0.98]
        #ijdel = [ij for ij,w in y_val.items() if abs(w)>0.98]  #(dbg!)
        g_yyy.delete_edges(a_del)        
        ### FOR VIZ ONLY!
        ##g_yyy.es['label'] = ['%.2f'%w   for w in g_yyy.es['weight']]
        ##igraphmod.igraph2graphviz(g_yyy,'pyconncmp-dfj')
        ccs = g_yyy.clusters(mode="weak")
        islsw = model._islsw
        n_grp = len(islsw)        
        if len(ccs)<n_grp:
            raise RuntimeError('Infeasible incumbent with too few connected components')
        if len(ccs)==n_grp:
            if fl_kvl:
                kvlcut(model,a_del)                
                return
            else:
                return                

        if not fl_zij:
            raise RuntimeError('Without z_ij variables (DFJ constraints), incumbent should not have too many connected components...')
        else:
            z_val = model.cbGetSolution(model._z)
            zijcut(model, z_val, ccs, a_del)
