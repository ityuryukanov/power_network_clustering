"""
This file contains callback functions for various main/first-level optimization
functions.
There is a option to change the solution strategy at a certain percentage of the 
full optimization time limit.
"""

from gurobipy import GRB, tupledict, quicksum
from timeit import default_timer as timer
from cut_callbk import zijcut
from cut_callbk import kvlcut
from cut_callbk import get_cutsets, get_dfjcuts
import warnings
import math
import os                    #for dbg
import sys                   #for dbg
rootpth = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..\\..\\'))
sys.path.append(rootpth)     #for dbg
import igraphmod             #for dbg

#@profile
def dfjcut_callback(model, where):
    if where==GRB.Callback.MIPNODE:
        # Change solution strategy when the runtime exceeds a certain limit
        if not model._changeParam and model.cbGet(GRB.Callback.RUNTIME)>model.Params.TimeLimit*0.01:
            model._changeParam = True
            model.terminate()

        # Add cuts at the root node
        if model.cbGet(GRB.Callback.MIPNODE_NODCNT)==0:
            fl_zij = model._ztruu
            if fl_zij:
                start = timer()
                z_val = model.cbGetNodeRel(model._z)
                # Just add all possible connectivity cuts at the root node as long
                # as the root node is being processed (it's usually not too long...)
                n_grp = len(model._islsw)
                subtours = get_dfjcuts(model._graph,z_val,n_grp,model._a2i,'mipnod',siz_max=7)
                for subtour in subtours:
                    arcs = subtour[0]
                    nodes = subtour[1]
                    model.cbLazy(quicksum(model._z[arc[0],arc[1]] for arc in arcs) <= len(nodes)-1)
                    model.cbLazy(quicksum(model._z[arc[1],arc[0]] for arc in arcs) <= len(nodes)-1) 
                cutsets = get_cutsets(model._graph,z_val,model._a2i,model._rtd,siz_max=10)
                for cutset in cutsets:
                    model.cbLazy(quicksum(model._z[i,j] for i,j in cutset) >= 1)                    
                model._troot += timer() - start        
        
        # Next proceed with heuristics regardless of the MIP node index
        if model._MipHeur['enabled']:
            if model._Theur<0:
                return
            start = timer()
            # --------------------------------------------------------------------
            # Push the heuristic until feasibility, then become more gentle 
            soln_cnt = model.cbGet(GRB.Callback.MIPNODE_SOLCNT)            
            if soln_cnt>0:
                return 
            if model._infeascycle>=8: 
                return                  # too many times MIP heuristic infeasibilities
            infeasskips=11 
            if model._infeascount>0 and model._infeascount<infeasskips:
                model._infeascount+=1   # wait one more call 
                return
            elif model._infeascount>=infeasskips:
                model._infeascount=0             
                model._infeascycle+=1    
            # --------------------------------------------------------------------
            nodes = model._i2n.keys()
            n_nod = len(nodes)
            x_lvl = model.cbGetNodeRel(model._x)
            islid = model._k2trm.keys()
            a2i   = model._a2i
            m_heur= model._mheur
            thrsh = m_heur._thrsh
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
            k_fix = [0] * n_nod
            e2i = model._e2i
            for k in islid:
                cc = k2ccr[k]
                for i,j in e2i.keys():
                    if x_lvl[i,k]>thrsh and x_lvl[j,k]>thrsh and i in cc and j in cc:
                        y_xxx[i,j] = 0
                        x_xxx[i,k] = 1
                        x_xxx[j,k] = 1
                        x_fix[i] = 1; x_fix[j] = 1;
                        k_fix[i] = k; k_fix[j] = k;
                        y_fix += 1                        
            # Roots are fixed regardless of the state of their neighbours
            for r in roots:
                x_xxx[r,trm2k[r]] = 1
                x_fix[r] = 1
                k_fix[r] = trm2k[r]             
            k_num = dict()
            for k in islid:
                k_num[k] = k_fix.count(k)    
            # #------------------------------Visualize------------------------------------  
            # COLRS = ['#AAAAAA','#FF0000','#00FF00','#0000FF','#FF00FF','#00FFFF','#00FFFF','#AAFFAA','#AAAAFF']
            # g_all = model._ntwrk.copy()                  #(dbg!)
            # g_all.vs['color'] = ['#AAAAAA']*n_nod        #(dbg!)
            # for i in nodes:
            #     g_all.vs[i]['color'] = COLRS[k_fix[i]]   #(dbg!)
            # igraphmod.igraph2graphviz(g_all,'heurist')   #(dbg!)
            # #------------------------------------------------------------------------------
            s_fix = sum(x_fix)             
            if s_fix/n_nod<m_heur._roudg:
                print("[ROUNDG] TOTAL TIME LEFT FOR HEURISTICS IS: " + str(model._Theur))
                model._infeascount+=1   # skip some relaxations
                # Reduce the threshold to attempt to solve the problem  
                if m_heur._roudg>0.4:
                    m_heur._roudg = m_heur._roudg - 0.1  
                model._Theur = model._Theur - (timer() - start)
                return   # too few nodes can be confidently fixed to reduce the problem
            # Fix the auxiliary model...
            for i,k in x_xxx.keys():
                m_heur._x[i,k].lb = x_xxx[i,k]
                m_heur._x[i,k].ub = x_xxx[i,k]
            for i,j in y_xxx.keys():
                m_heur._y[i,j].lb = y_xxx[i,j]
                m_heur._y[i,j].ub = y_xxx[i,j]    
            # Adjust the connectivity flow constraints to speed up the solution 
            # (changes each time equally for all branches)
            M_nod = n_nod-s_fix-1
            for i,j in m_heur._y.keys():
                m_heur._T[i,j].lb = -M_nod
                m_heur._T[i,j].ub = +M_nod
            # Switching & Commodity flow limitation
            # (is set each time equal for all branches)            
            for i,j in m_heur._y.keys():
                m_heur._NFanHI[i,j].RHS = +M_nod
                m_heur.chgCoeff(m_heur._NFanHI[i,j],m_heur._y[i,j],+M_nod) 
                m_heur._NFanLO[i,j].RHS = -M_nod
                m_heur.chgCoeff(m_heur._NFanLO[i,j],m_heur._y[i,j],-M_nod)                 
            # Flow balance
            for r in roots:
                m_heur._NFanCo[r].RHS = -k_num[trm2k[r]]
            for i in m_heur._nonsw:
                if x_fix[i]>0.2:
                    m_heur._NFanCo[i].RHS = 0
            # Optimize the auxiliary model....
            m_heur.optimize()
            # Return if infeasible or not enough improvement            
            if m_heur.status==GRB.Status.INF_OR_UNBD or m_heur.status==GRB.Status.INFEASIBLE or m_heur.status==GRB.Status.UNBOUNDED or m_heur.status==GRB.Status.TIME_LIMIT or m_heur.status==GRB.Status.NODE_LIMIT or m_heur.status==GRB.Status.ITERATION_LIMIT:
                # m_heur.computeIIS()                     #(dbg!)
                # m_heur.write('iis_heur.ilp')            #(dbg!)            
                model._infeascount+=1
                # Restore the variable bounds
                for i,k in m_heur._x.keys():
                    m_heur._x[i,k].lb = m_heur._x_lbnd[i,k]
                    m_heur._x[i,k].ub = m_heur._x_ubnd[i,k]
                for i,j in m_heur._y.keys():
                    m_heur._y[i,j].lb = m_heur._y_lbnd[i,j]
                    m_heur._y[i,j].ub = m_heur._y_ubnd[i,j]
                for i in nodes:
                    m_heur._NFanCo[i].RHS = -1       
                if m_heur.status==GRB.Status.INF_OR_UNBD or m_heur.status==GRB.Status.INFEASIBLE or m_heur.status==GRB.Status.UNBOUNDED:
                    # Increase the threshold to round up more carefully 
                    # and thus lower the chance of infeasibility.
                    if m_heur._thrsh<0.9995:
                        m_heur._thrsh = m_heur._thrsh + 10**(-m_heur._thrshacc-1)*9
                        m_heur._thrshacc+=1
                    # Also increase the bounds for voltage angles. This
                    # will only be effective if the main model itself   
                    # has no volt. angle variables (e.g., a KVL approach).
                    for i in nodes:
                        m_heur._ph[i].lb = 2*m_heur._ph[i].lb
                        m_heur._ph[i].ub = 2*m_heur._ph[i].ub
                    for i,j in m_heur._y.keys():
                        U_ij = m_heur.getCoeff(m_heur._AngLimU[i,j], m_heur._y[i,j])
                        L_ij = m_heur.getCoeff(m_heur._AngLimL[i,j], m_heur._y[i,j])
                        m_heur.chgCoeff(m_heur._AngLimU[i,j], m_heur._y[i,j], 2*U_ij)
                        m_heur.chgCoeff(m_heur._AngLimL[i,j], m_heur._y[i,j], 2*L_ij)
                    print("[INFEAS] TOTAL TIME LEFT FOR HEURISTICS IS: " + str(model._Theur))
                if m_heur.status==GRB.Status.TIME_LIMIT or m_heur.status==GRB.Status.NODE_LIMIT or m_heur.status==GRB.Status.ITERATION_LIMIT:                    
                    print("[TLIMIT] TOTAL TIME LEFT FOR HEURISTICS IS: " + str(model._Theur))                
                model._Theur = model._Theur - (timer() - start)               
                return            
            ## Restore the variable bounds
            for i,k in m_heur._x.keys():
                m_heur._x[i,k].lb = m_heur._x_lbnd[i,k]
                m_heur._x[i,k].ub = m_heur._x_ubnd[i,k]
            for i,j in m_heur._y.keys():
                m_heur._y[i,j].lb = m_heur._y_lbnd[i,j]
                m_heur._y[i,j].ub = m_heur._y_ubnd[i,j]
            for i in nodes:
                m_heur._NFanCo[i].RHS = -1                                   
        #--------------------------------------------------------------------------------------------                     
            x_xxx = tupledict()
            y_xxx = tupledict()
            p_xxx = tupledict()               
            pGS_x = tupledict()
            pL_x  = tupledict()
            for i,k in model._x.keys():
                x_xxx[i,k] = m_heur._x[i,k].x
            for i,j in model._y.keys():
                y_xxx[i,j] = m_heur._y[i,j].x
            for i,j in model._p.keys():
                p_xxx[i,j] = m_heur._p[i,j].x
            for i in model._pGS.keys():
                pGS_x[i] = m_heur._pGS[i].x
            for i in model._pLS.keys():
                pL_x[i] = m_heur._pLS[i].x
            model.cbSetSolution(model._x, x_xxx)
            model.cbSetSolution(model._y, y_xxx)
            model.cbSetSolution(model._p, p_xxx)
            model.cbSetSolution(model._pGS, pGS_x)
            model.cbSetSolution(model._pLS, pL_x)  
            if hasattr(model,"_ph"):
                phi_x = tupledict()
                for i in model._ph.keys():
                    phi_x[i] = m_heur._ph[i].x
                model.cbSetSolution(model._ph, phi_x)
            if hasattr(model,"_z") or model._ztruu:
                z_xxx = tupledict()
                g_sol = model._graph.copy()
                a_off = list()
                for i,j in y_xxx:
                    z_xxx[i,j]=0
                    z_xxx[j,i]=0
                    if y_xxx[i,j]>0.95:
                        a_off.append(a2i[i,j])
                        a_off.append(a2i[j,i])
                g_sol.delete_edges(a_off)
                for k in islid:
                    R = k2src[k]
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        esIDs = g_sol.get_shortest_paths(R, to=None, weights=g_sol.es['weight'], mode='out', output='epath')
                    esIDs = {esID for path in esIDs for esID in path if len(path)>0}
                    for esID in esIDs: z_xxx[g_sol.es[esID].tuple] = 1     
                model.cbSetSolution(model._z, z_xxx)               
            ##objval = model.cbUseSolution()      #(dbg!)
            print("[NEWSOL=" + str(m_heur.ObjVal) + "] TOTAL TIME LEFT FOR HEURISTICS IS: " + str(model._Theur))
            model._Theur = model._Theur - (timer() - start)
    
    if where == GRB.Callback.MIPSOL:
        start = timer()
        fl_zij = model._ztruu
        fl_kvl = model._kvlfl
        if not fl_zij and not fl_kvl:
            if model._Tfeas<0:
                sol_time = model.cbGet(GRB.Callback.RUNTIME)
                model._Tfeas = sol_time
            return
        y_val = model.cbGetSolution(model._y)
        g_yyy = model._ntwrk.copy()
        e2i   = model._e2i
        a_del = [e2i[ij] for ij,w in y_val.items() if abs(w)>0.98]
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
            model._ttree += timer() - start
            if fl_kvl:
                dc_pf_okk = kvlcut(model,a_del)
                if dc_pf_okk and model._Tfeas<0:
                    sol_time = model.cbGet(GRB.Callback.RUNTIME)
                    model._Tfeas = sol_time
                return
            else:
                if model._Tfeas<0:
                    sol_time = model.cbGet(GRB.Callback.RUNTIME)
                    model._Tfeas = sol_time
                return

        if not fl_zij:
            raise RuntimeError('Without z_ij variables (DFJ constraints), incumbent should not have too many connected components...')
        else:
            z_val = model.cbGetSolution(model._z)
            connected = zijcut(model, z_val, ccs, a_del)
            if connected and model._Tfeas<0:
                sol_time = model.cbGet(GRB.Callback.RUNTIME)
                model._Tfeas = sol_time                  
        
        model._ttree += timer() - start


