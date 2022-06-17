"""
This module contains functions for generation of all needed cuts, to be used 
in MIP callback routines of various types (e.g., in callbacks of heuristic model
solves and callbacks of main model solves).
"""

from gurobipy import quicksum, GRB
from utils import mst_cb, check_KVL
from termcolor import colored
import warnings
import math
import os                    #for dbg
import sys                   #for dbg
rootpth = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..\\..\\'))
sys.path.append(rootpth)     #for dbg
import igraphmod             #for dbg

#@profile
def get_cutsets(g_inp, z_lvl, a2i, rtd, siz_max=1e9):
    """
    Although the returned cuts can be non-obvious (a minimum cut value, but 
    not a cut closest to the source or target),  they are the same for both 
    inverted and uninverted edges.
    """
    g_all = g_inp.copy()      
    for ij,w in z_lvl.items():
        g_all.es[a2i[(ij[1],ij[0])]]['weight'] = abs(w)             # revert edges for back cuts 
    g_all.es['weight'] = [w+0.00001  for w in g_all.es['weight']]   # for creep flow add a small extra weight!
    ### FOR VIZ ONLY!     
    ##g_all.es['label'] = ['%.2f'%w   for w in g_all.es['weight']]
    ##igraphmod.igraph2graphviz(g_all,'pydicuts-all')
    cuts = list() 
    for R in rtd.keys():                    
        for td in rtd[R]:
            t = td[0]   #d=td[1]             
            cut = g_all.mincut(t,R,capacity=g_all.es['weight'])       
            if cut.value<0.98:
                tst = [z_lvl[(g_all.es[arcID].target,g_all.es[arcID].source)] for arcID in cut.cut]        #to test arc existence                
                z_c = frozenset({(g_all.es[arcID].target,g_all.es[arcID].source) for arcID in cut.cut})    #revert the back-cut for the constraint      
                if z_c not in cuts and len(z_c)<=siz_max:
                    cuts.append(z_c)
            ##if len(cuts)>20: break  
    return cuts
    
#@profile
def cycles_dfj(graph, modd='mipsol', siz_max=1e9):    
    esIDs = graph.spanning_tree(return_tree=False, weights=graph.es['weight'])
    esALL = [a.index for a in graph.es]    
    esCYC = set(esALL) - set(esIDs)
    cycles = list()
    if len(esCYC)==0:
        return cycles
    for esID in esCYC:
        f = graph.es[esID].source
        t = graph.es[esID].target
        # the t->f path direction is both valid for arborescences formed by the 
        # arcs "flowing" out of the root node and into the root node        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")        
            escyc = graph.get_shortest_paths(t,to=f,mode="out",output="epath")[0]  #mode=OUT !            
            cycle = [(graph.es[esix].source, graph.es[esix].target) for esix in escyc] + [(graph.es[esID].source,graph.es[esID].target)]
        try: 
            cycOK = [cycle[idx][1]==cycle[idx+1][0] for idx,arc in enumerate(cycle) if idx<=len(cycle)-2]            
            if not all(cycOK) or not cycle[-1][1]==cycle[0][0]: 
                continue
        except: 
            continue
        if len(cycle)>1 and modd=='mipsol' and len(cycle)<=siz_max:
            cycles.append(set(cycle))
        elif len(cycle)>1:            
            z_sum = sum([abs(graph.es[esix]['weight']) for esix in escyc])+abs(graph.es[esID]['weight'])
            dviol = z_sum - len(cycle) + 1
            if dviol>0 and len(cycle)<=siz_max:
                cycles.append((dviol,set(cycle)))
    # As any non-root vertex has at most one "inflowing" or "outflowing" arc 
    # (depending on the "inflowing" or "outflowing" formulation type), nested
    # cycles are impossible for an integer solution, and each cycle certainly 
    # contains a unique set of nodes.                
    if modd=='mipsol':            
        return cycles
    # Return only most strongly violated set packing constraints when dealing  
    # with fractional solutions (i.e., modd!='mipsol').
    cycles = sorted(cycles, key=lambda x: (-x[0],len(x[1])), reverse=False)
    vxseen = set()
    cycdel = list()
    cycles = [cycle[1] for cycle in cycles]
    for cycle in cycles:
        if len(cycle.intersection(vxseen))>0:
            cycdel.append(cycle)
        else:
            [vxseen.add(vx) for vx in cycle]
    for cycle in cycdel:
        cycles.remove(cycle)
    return cycles

        
#@profile
def get_dfjcuts(g_inp, z_lvl, n_grp, a2i, modd, siz_max=1e9):   
    """
    Finds and returns DFJ subtours.
    """
    subtours = list()
    g_dfj = g_inp.copy()
    for ij,w in z_lvl.items():
        g_dfj.es[a2i[ij]]['weight'] = abs(w)
    a_del = [a.index for a in g_dfj.es if abs(a['weight'])<1e-5]
    g_dfj.delete_edges(a_del)        
    ### FOR VIZ ONLY!     
    ##g_dfj.es['label'] = ['%.2f'%w   for w in g_dfj.es['weight']]
    ##igraphmod.igraph2graphviz(g_dfj,'pydicuts-dfj')
    cc = g_dfj.clusters(mode="weak")  #mode=WEAK
    if modd!='mipsol' or len(cc)>n_grp:
        cycles = cycles_dfj(g_dfj, modd, siz_max)
        for arcs in cycles:
            subtour = list()
            subtour.append(arcs)
            cycle = {v for arc in arcs for v in arc}
            subtour.append(cycle)          
            subtours.append(tuple(subtour))                                          
        return subtours    
    elif len(cc)==n_grp:
        return subtours
    else:
        raise RuntimeError('Number of connected components is lower than requested number of groups of terminals')


#@profile
def zijcut(model, z_val, ccs, a_del):    
    n_grp = len(model._islsw)
    subtours = get_dfjcuts(model._graph,z_val,n_grp,model._a2i,'mipsol')
    for subtour in subtours:
        arcs = subtour[0]
        nodes = subtour[1]
        model.cbLazy(quicksum(model._z[arc[0],arc[1]] for arc in arcs) <= len(nodes)-1)
        model.cbLazy(quicksum(model._z[arc[1],arc[0]] for arc in arcs) <= len(nodes)-1)            
   
    # These s-t cuts do help to speed up certain cases a lot. A reasonable
    # cost for a possible slight average runtime increase.  Comment out if 
    # they are inserted a level above.
    cutsets = get_cutsets(model._graph,z_val,model._a2i,model._rtd)
    for cutset in cutsets:
        model.cbLazy(quicksum(model._z[i,j] for i,j in cutset) >= 1)
    
    # Output flag
    connected = len(subtours)==0

    # For each connected component without a root node, at least one input 
    # arc must be enabled.
    vx2ic = dict()
    for icc,cc in enumerate(ccs):
        for vx in cc:
            vx2ic[vx] = icc
    srcic = [icc for ir,r in enumerate(model._islsw) for icc,cc in enumerate(ccs) if r in cc]
    mnric = set(range(len(ccs))) - set(srcic)    #ccs without root nodes
    kr2ij = dict((ic,list()) for ic in mnric)
    i2e   = model._i2e
    for ieg in a_del:
        i,j = i2e[ieg]
        if vx2ic[i] in mnric:
            kr2ij[vx2ic[i]].append((j,i))
        if vx2ic[j] in mnric:    
            kr2ij[vx2ic[j]].append((i,j))
    for icc in kr2ij.keys():
        if len(kr2ij[icc])>0:
            model.cbLazy(quicksum(model._z[i,j] for i,j in kr2ij[icc])>=1)
    return connected


def kvlcut(model,a_del,modd='mipsol'):
    """
    Check KVL and insert violated KVL cycles as lazy cuts.
    """
    n_grp = len(model._islsw)
    grf_sp = model._grfsp.copy()
    grf_sp.delete_edges(a_del)
    esIDs = grf_sp.spanning_tree(return_tree=False, weights=grf_sp.es['weight'])        
    mstcb = mst_cb(grf_sp, esIDs)
    assert(len(mstcb)==len(grf_sp.es)-len(grf_sp.vs)+n_grp)
    p_val = model.cbGetSolution(model._p)
    b_val = model._B    # for cuts insertion in callbacks
    M_val = model._M    # for cuts insertion in callbacks
    fails = check_KVL(mstcb,p_val,b_val)
    cycfail = fails[0]
    csmfail = fails[1]
    for cyc in cycfail:
        th_ij = list()
        for i,j,pm in cyc:
            th_ij.append(abs(M_val[(i,j)])/abs(b_val[(i,j)]))   # abs()/abs() is the "worst case"
        th_ij.sort(reverse=True)
        th_ij.pop(); th_ij.pop();
        M_cyc = sum(th_ij)
        model.cbLazy(quicksum(model._p[(i,j)]/b_val[(i,j)]*pm for i,j,pm in cyc) - 0.5*M_cyc*quicksum(model._y[i,j] for i,j,pm in cyc) <= 0)
        model.cbLazy(quicksum(model._p[(i,j)]/b_val[(i,j)]*pm for i,j,pm in cyc) + 0.5*M_cyc*quicksum(model._y[i,j] for i,j,pm in cyc) >= 0)
        model._cbkcyc+=1
        #csm = sum([p_val[(i,j)]/b_val[(i,j)]*pm for i,j,pm in cyc]);  #(dbg!)
        #ysm = sum([y_val[(i,j)] for i,j,pm in cyc])   #(dbg!)                    
        #print('')
        #print(colored(str(cyc),'red'),colored(str(csm),'magenta'),colored(str(ysm),'cyan'))  #(dbg!)
        #print('')
    # Output flag
    dc_pf_okk = len(cycfail)==0
    return dc_pf_okk


