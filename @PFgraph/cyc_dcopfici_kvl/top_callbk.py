"""
This file contains first-level callback functions (with model and where args) 
for various main and auxiliary (e.g., MIP heuristics) functions.  
"""

from gurobipy import GRB, quicksum, tupledict
from timeit import default_timer as timer
from cut_callbk import zijcut
from cut_callbk import kvlcut
from hrs_callbk import dfjhrs_callback
import warnings
import math
import os                    #dbg
import sys                   #dbg
rootpth = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..\\..\\'))
sys.path.append(rootpth)     #dbg
import igraphmod             #dbg

#@profile
def dfjcut_callback(model, where):
    fl_zij = model._ztruu
    fl_kvl = model._kvlfl
    if where==GRB.Callback.MIPNODE:
        # if model.cbGet(GRB.Callback.MIPNODE_NODCNT)==0:
        #     start = timer()
        #     # Just add all possible connectivity cuts at the root node as long 
        #     # as the root node is being processed (it's usually not too long...)
        #     z_val = model.cbGetNodeRel(model._z)
        #     n_grp = len(model._islsw)
        #     subtours = get_dfjcuts(model._graph,z_val,n_grp,model._a2i,'mipnod')
        #     for subtour in subtours:
        #         arcs = subtour[0]
        #         nodes = subtour[1]
        #         model.cbLazy(quicksum(model._z[arc[0],arc[1]] for arc in arcs) <= len(nodes)-1)
        #         model.cbLazy(quicksum(model._z[arc[1],arc[0]] for arc in arcs) <= len(nodes)-1) 
        #     cutsets = get_cutsets(model._graph,z_val,model._a2i,model._rtd)
        #     for cutset in cutsets:
        #         model.cbLazy(quicksum(model._z[i,j] for i,j in cutset) >= 1)                    
        #     model._troot += timer() - start
        
        # Next proceed with heuristics regardless of the MIP node index		
        if model._MipHeur['enabled']:
            if model._Theur<0:
                return
            start = timer()
            # ----------------------------------------------------------------
            obj_bnd = model.cbGet(GRB.Callback.MIPNODE_OBJBND)                #(dbg!)
            obj_bst = model.cbGet(GRB.Callback.MIPNODE_OBJBST)                #(dbg!)
            mip_gap = abs(obj_bst - obj_bnd)/max(abs(obj_bst),abs(obj_bnd))   #(dbg!)
            # Push the heuristic until feasibility, then become more gentle
            gaplast = model._xlast  
            soln_cnt = model.cbGet(GRB.Callback.MIPNODE_SOLCNT)
            sol_time = model.cbGet(GRB.Callback.RUNTIME)
            if soln_cnt>0 or sol_time<7: 
                return 
            # ----------------------------------------------------------------
            nodes = model._i2n
            n_nod = len(nodes)
            x_lvl = model.cbGetNodeRel(model._x)
            islid = model._k2trm.keys()             
            a2i   = model._a2i
            thrsh = 0.99
            a_onn = [a2i[i,j]  for k in islid  for i,j in a2i.keys()  if x_lvl[i,k]>thrsh and x_lvl[j,k]>thrsh]
            g_xxx = model._graph.copy()
            a_off = set([a.index for a in g_xxx.es]) - set(a_onn)        
            g_xxx.delete_edges(a_off) 
            ccs = g_xxx.clusters(mode="weak")                
            roots = set(model._islsw)
            trm2k = model._trm2k
            k2ccr = dict()   # k to root connected components
            k2src = dict()
            for cc in ccs:
                R = roots.intersection(set(cc))            
                if len(R)==1:
                    R = R.pop()
                elif len(R)==0:
                    continue
                else:
                    print('Two root nodes in one cc - relaxation not nice..')
                    model._Theur = model._Theur - (timer() - start)
                    return                    
                k = trm2k[R]
                k2ccr[k] = cc 
                k2src[k] = R
            # ----------------------
            x_xxx = tupledict()
            y_xxx = tupledict()
            n_grp = len(islid)            
            x_fix = [0] * n_nod
            y_fix = 0
            e2i = model._e2i
            for k in islid:
                cc = k2ccr[k]
                for i,j in e2i.keys():
                    if x_lvl[i,k]>thrsh and x_lvl[j,k]>thrsh and i in cc and j in cc:
                        y_xxx[i,j] = 0
                        x_xxx[i,k] = 1
                        x_xxx[j,k] = 1
                        x_fix[i] = 1; x_fix[j] = 1;
                        y_fix += 1
            x_fix = sum(x_fix)
            m_heur = model._mheur
            if x_fix/n_nod<0.8:
                print("[ROUNDG] TOTAL TIME LEFT FOR HEURISTICS IS: " + str(model._Theur))
                model._Theur = model._Theur - (timer() - start)
                return   # too few nodes can be confidently fixed to reduce the problem
            for i,k in x_xxx.keys():
                m_heur._x[i,k].lb = x_xxx[i,k]
                m_heur._x[i,k].ub = x_xxx[i,k]
            for i,j in y_xxx.keys():
                m_heur._y[i,j].lb = y_xxx[i,j]
                m_heur._y[i,j].ub = y_xxx[i,j]
            if fl_zij:
                z_xxx = tupledict()   # commpute some z[i,j] for further speedups
                for k in islid:
                    R = k2src[k]
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        esIDs = g_xxx.get_shortest_paths(R, to=None, weights=g_xxx.es['weight'], mode='out', output='epath')
                    esIDs = {esID for path in esIDs for esID in path if len(path)>0}
                    for esID in esIDs: z_xxx[g_xxx.es[esID].tuple] = 1				
                for i,j in z_xxx.keys():
                    m_heur._z[i,j].lb = z_xxx[i,j]
                    m_heur._z[i,j].ub = z_xxx[i,j]
            if fl_zij:
                m_heur.optimize(dfjhrs_callback)
            else:
                m_heur.optimize()
            # Return if infeasible or not enough improvement            
            if m_heur.status==GRB.Status.INF_OR_UNBD or m_heur.status==GRB.Status.INFEASIBLE or m_heur.status==GRB.Status.UNBOUNDED:
                ## m_heur.computeIIS()                #(dbg!)
                ## m_heur.write('iis_heur.ilp')       #(dbg!)            
                # Restore the variable bounds
                for i,k in m_heur._x.keys():
                    m_heur._x[i,k].lb = m_heur._x_lbnd[i,k]
                    m_heur._x[i,k].ub = m_heur._x_ubnd[i,k]
                for i,j in m_heur._y.keys():
                    m_heur._y[i,j].lb = m_heur._y_lbnd[i,j]
                    m_heur._y[i,j].ub = m_heur._y_ubnd[i,j]
                if fl_zij:
                    for i,j in m_heur._z.keys():
                        m_heur._z[i,j].lb = m_heur._z_lbnd[i,j]
                        m_heur._z[i,j].ub = m_heur._z_ubnd[i,j]
                print("[INFEAS] TOTAL TIME LEFT FOR HEURISTICS IS: " + str(model._Theur))
                model._Theur = model._Theur - (timer() - start)                
                return
            else:            
                model._xlast = mip_gap

            x_xxx = tupledict()
            y_xxx = tupledict()
            p_xxx = tupledict()   
            t_xxx = tupledict()
            pGS_x = tupledict()
            pLS_x = tupledict()
            if fl_zij:
                z_xxx = tupledict()
                for i,j in model._z.keys():
                    z_xxx[i,j] = m_heur._z[i,j].x
            for i,k in model._x.keys():
                x_xxx[i,k] = m_heur._x[i,k].x
            for i,j in model._y.keys():
                y_xxx[i,j] = m_heur._y[i,j].x
            for i,j in model._p.keys():
                p_xxx[i,j] = m_heur._p[i,j].x
            for i in model._pGS.keys():
                pGS_x[i] = m_heur._pGS[i].x
            for i in model._pLS.keys():
                pLS_x[i] = m_heur._pLS[i].x
            if m_heur._T is not None:
                for i,j in model._T.keys():
                    t_xxx[i,j] = m_heur._T[i,j].x
            #z_onn = [(i+1,j+1) for i,j in z_xxx.keys() if z_xxx[i,j]>thrsh]    #(dbg!)
            model.cbSetSolution(model._x, x_xxx)
            model.cbSetSolution(model._y, y_xxx)
            if fl_zij:
                model.cbSetSolution(model._z, z_xxx)
            model.cbSetSolution(model._p, p_xxx)
            if m_heur._T is not None:
                model.cbSetSolution(model._T, t_xxx)
            model.cbSetSolution(model._pGS, pGS_x)
            model.cbSetSolution(model._pLS, pLS_x)                
            # Restore the variable bounds
            for i,k in m_heur._x.keys():
                m_heur._x[i,k].lb = m_heur._x_lbnd[i,k]
                m_heur._x[i,k].ub = m_heur._x_ubnd[i,k]
            for i,j in m_heur._y.keys():
                m_heur._y[i,j].lb = m_heur._y_lbnd[i,j]
                m_heur._y[i,j].ub = m_heur._y_ubnd[i,j]
            if fl_zij:
                for i,j in m_heur._z.keys():
                    m_heur._z[i,j].lb = m_heur._z_lbnd[i,j]
                    m_heur._z[i,j].ub = m_heur._z_ubnd[i,j]
            ##objval = model.cbUseSolution()      #(dbg!)
            foo = 5;         #(dbg!)
            print("[NEWSOL=" + str(m_heur.ObjVal) + "] TOTAL TIME LEFT FOR HEURISTICS IS: " + str(model._Theur))
            model._Theur = model._Theur - (timer() - start)
            
    #assert(model._graph.is_directed(), "The network graph must be directed")
    if where == GRB.Callback.MIPSOL:
        start = timer()
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
        if len(ccs)==n_grp and fl_kvl:
            kvlcut(model,a_del)
            model._ttree += timer() - start
            return
        elif len(ccs)<n_grp:
            raise RuntimeError('Infeasible incumbent with too few connected components')
        
        if not fl_zij:
            model._ttree += timer() - start
            return            
        else:
            zijcut(model, y_val, ccs, a_del)

        model._ttree += timer() - start



