import os
import re
import warnings  
from timeit import default_timer as timer
from gurobipy import GRB, quicksum, tupledict
from collections import Counter

def order_checking(cyc_ord):    
    arc_ord = [(cyc_ord[0][0],cyc_ord[0][1])]
    for idx in range(1,len(cyc_ord)):
        if arc_ord[idx-1][1]==cyc_ord[idx][0]:
            arc = (cyc_ord[idx][0],cyc_ord[idx][1])
        else:
            arc = (cyc_ord[idx][1],cyc_ord[idx][0])  
        arc_ord.append(arc)   
    assert(all([arc_ord[idx][1]==arc_ord[idx+1][0] for idx in range(len(arc_ord)-1)]))
    assert(arc_ord[0][0]==arc_ord[-1][1])    
    return arc_ord

def mst_cb(netw, esIDs):
    graph = netw.copy()
    ##for esID in esIDs:
    ##    print((graph.es[esID].source+1,graph.es[esID].target+1))        
    esALL = [a.index for a in graph.es]
    esCYC = set(esALL) - set(esIDs)
    cycles = list()
    for es in graph.es:
        es['weight'] = 1
    big_C = len(graph.es)+5
    if len(esCYC)==0:
        return cycles
    for esID in esCYC:
        f = graph.es[esID].source
        t = graph.es[esID].target
        # the t->f path direction is both valid for arborescences formed by the 
        # arcs "flowing" out of the root node and into the root node        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            graph.es[esID]['weight']=big_C  ## also excludes that cycle for the future
            escyc = graph.get_shortest_paths(t,to=f,mode="all",output="epath",weights=graph.es['weight'])[0]              
            cycle = [(graph.es[esix].source, graph.es[esix].target) for esix in escyc] + [(graph.es[esID].source,graph.es[esID].target)]
        try: 
            cycOK = all([val==2 for val in Counter(list(sum(cycle, ()))).values()])            
            if not cycOK: 
                continue
        except: 
            continue
        if len(cycle)>1:
            cycles.append(tuple(sorted(cycle)))
    # As any non-root vertex has at most one "inflowing" or "outflowing" arc 
    # (depending on the "inflowing" or "outflowing" formulation type), nested
    # cycles are impossible for an integer solution, and each cycle certainly 
    # contains a unique set of nodes.
    return cycles

def check_KVL(mstcb,p_sol,b_val):
    """
    Cycles should be filtered to exclude those with switched-off branches.
    Otherwise the result can be invalid.
    """
    mcb_lin = list()
    for cyc0 in mstcb:
        cyc_sgn = cycle_signing(cyc0)
        mcb_lin.append(cyc_sgn)    
    cycfail = list()
    csmfail = list()
    for cyc in mcb_lin:       
        csm = sum([p_sol[(i,j)]/b_val[(i,j)]*pm for i,j,pm in cyc]);
        if abs(csm)/len(cyc)>1e-4:
            cycfail.append(cyc)
            csmfail.append(csm)
    return (cycfail,csmfail)

def cycle_ordering(cyc0):
    # cyc0 must be a cycle for this function
    # to succeed.
    cyc_lin = list(cyc0)            
    lin0 = cyc_lin.pop(0)
    jnxt = lin0[1]
    cyc_ord = list()  # ordered cyc_lin
    cyc_ord.append(lin0)
    while cyc_lin:
        lin_nxt = [lin for lin in cyc_lin if lin[0]==jnxt or lin[1]==jnxt]
        assert(len(lin_nxt)==1)
        lin_nxt = lin_nxt[0]
        cyc_ord.append((lin_nxt[0],lin_nxt[1]))
        if lin_nxt[0]==jnxt:
            jnxt = lin_nxt[1]
        elif lin_nxt[1]==jnxt:
            jnxt = lin_nxt[0]
        cyc_lin.remove(lin_nxt)
    order_checking(cyc_ord)
    return cyc_ord

def cycle_signing(cyc0):    
    cyc_ord = list(cyc0)
    cyc_sgn = list()
    lin0 = cyc_ord.pop(0)
    cyc_sgn.append((lin0[0],lin0[1],+1))
    while cyc_ord:
        lin0 = cyc_sgn[-1]
        lin_nxt = [lin for lin in cyc_ord if (lin[0] in lin0[0:2]) or (lin[1] in lin0[0:2])]
        if len(lin_nxt)==0:
            lin0 = cyc_sgn[0]
            lin_nxt = [lin for lin in cyc_ord if (lin[0] in lin0[0:2]) or (lin[1] in lin0[0:2])]
        assert(len(lin_nxt)>=1 and len(lin_nxt)<=2)
        lin_nxt = lin_nxt[0]
        if lin0[1]==lin_nxt[0] or lin_nxt[1]==lin0[0]:
            cyc_sgn.append((lin_nxt[0],lin_nxt[1],+lin0[2]))
        else:
            cyc_sgn.append((lin_nxt[0],lin_nxt[1],-lin0[2]))
        cyc_ord.remove(lin_nxt)
    # Test the result to always be sure
    cyc_ord = [(i,j)  if pm==1 else (j,i) for i,j,pm in cyc_sgn]
    vxcyc = Counter([vx for lin in cyc_ord for vx in lin])
    chord =  [vx for vx in vxcyc if vxcyc[vx]==1]
    if len(chord)==0:
        cyc_ord = cycle_ordering(cyc_ord)   #cyc_ord needs to be re-ordered for order_checking()
        order_checking(cyc_ord)   #it should have been done in cycle_ordering(), but still
    else:
        arc_ord = [lin for lin in cyc_ord if lin[0] in chord]
        cyc_ord.remove(arc_ord[0])
        while cyc_ord:
            arc = [(i,j) for i,j in cyc_ord if i==arc_ord[-1][1]][0]
            arc_ord.append(arc)
            cyc_ord.remove(arc)
        assert(all([arc_ord[idx][1]==arc_ord[idx+1][0] for idx in range(len(arc_ord)-1)]))
        assert(set((arc_ord[0][0],arc_ord[-1][1]))==set(chord))
    return cyc_sgn

def write_sol(model, varname, isnodid, folder):
    var = getattr(model,'_'+varname)                 
    keys = var.keys()
    key0 = keys[0]
    if hasattr(var[key0],'x'): isvar=True        
    else: isvar=False      
    if isvar:
        fl_path = os.path.join(folder, var.values()[0].VarName.split('[')[0]+".csv")
    else:
        fl_path = os.path.join(folder, varname[0:-1]+".csv")
    with open(fl_path,"w") as f:
        for key in var:
            if hasattr(key,"__len__"): 
                keys = key
            else:
                keys = [key]
            ln = ''
            for idx,ind in enumerate(keys):
                if isnodid[idx]:
                    ln += str(model._i2n[ind])+','  #refuse to convert to node labels at all 
                else:
                    ln += str(ind)+','
            if isvar:        
                ln = ln+str(var[key].x)
            else:
                ln = ln+str(var[key])
            f.write(ln+'\n')

def sol2stage(self,callbackfcn):
    self.model.Params.MIPFocus = 1
    self.model.Params.LazyConstraints = 1
    # DEVIATING FROM DEFALT SETTINGS IS NOT WORTH IT!
    #self.model.Params.FlowCoverCuts = 2
    #self.model.Params.FlowPathCuts = 2
    #self.model.Params.NetworkCuts = 2
    #for i in self.nodes:
    #    for k in self.islid:
    #        self.model._x[i,k].BranchPriority = 1
    # for i,j in self.uarcs:
    #     self.model._y[i,j].VarHintVal = 0
    # DEVIATING FROM DEFALT SETTINGS IS NOT WORTH IT!
                 
    # Try to obtain an initial feasible solution satisfying the KVL 
    self.model.update()
    MIPGap0 = self.model.Params.MIPGap
    TimeLim = self.model.Params.TimeLimit
    self.model.Params.MIPGap = 0.999999999999   # accept first feasible solution
    start = timer()
    while True:
        if self.model._rtd is not None: 
            self.model.optimize(callbackfcn)
        else:
            self.model.optimize()
        if self.model.status==GRB.Status.INF_OR_UNBD or self.model.status==GRB.Status.INFEASIBLE or self.model.status==GRB.Status.UNBOUNDED: 
            self.model.computeIIS()
            self.model.write("IIS_GP.ilp")
            return None
        y_sol = dict()
        p_sol = dict()
        for i,j in self.uarcs:
            y_ij = self.model._y[i,j].x
            p_ij = self.model._p[i,j].x
            y_sol[(i,j)] = y_ij;
            p_sol[(i,j)] = p_ij;
        # Check cycles KVL
        grf_sp = self.model._ntwrk.copy()  # assign sign to each cycle edge
        e2i    = self.model._e2i
        nc     = len(self.islid)
        a_del  = [e2i[ij] for ij,w in y_sol.items() if abs(w)>0.98]
        grf_sp.delete_edges(a_del)
        esIDs  = grf_sp.spanning_tree(return_tree=False, weights=None)
        mstcb  = mst_cb(grf_sp, esIDs)
        assert(len(mstcb)==len(grf_sp.es)-len(grf_sp.vs)+nc)
        cycfail = check_KVL(mstcb,p_sol,self.bdcpf)[0]
        if len(cycfail)>0:
            for cyc in cycfail:
                th_ij = list()
                for i,j,pm in cyc:
                    th_ij.append(abs(self.brlim[(i,j)])/abs(self.bdcpf[(i,j)]))   # abs()/abs() is the "worst case"
                th_ij.sort(reverse=True)
                th_ij.pop(); th_ij.pop()
                M_cyc = sum(th_ij)
                self.model.addConstr(quicksum(self.model._p[i,j]/self.bdcpf[(i,j)]*pm for i,j,pm in cyc)<=+0.5*M_cyc*quicksum(self.model._y[i,j] for i,j,pm in cyc))
                self.model.addConstr(quicksum(self.model._p[i,j]/self.bdcpf[(i,j)]*pm for i,j,pm in cyc)>=-0.5*M_cyc*quicksum(self.model._y[i,j] for i,j,pm in cyc))
            theur = timer() - start
        else:
            theur = timer() - start
            break
        if theur>TimeLim:
            return None
    self.model.Params.TimeLimit = TimeLim - theur
    self.model.Params.StartNumber = 1
    for i,k in self.model._x:
        self.model._x[i,k].start = self.model._x[(i,k)].x
    for i,j in self.model._z:
        self.model._z[i,j].start = self.model._z[(i,j)].x
    for i,j in self.model._y:
        self.model._y[i,j].start = self.model._y[(i,j)].x
    self.model.Params.MIPGap = MIPGap0
    return self         


def make_mheur(model):
    # Create MIP model for heuristics
    MIPHEUR = model._MipHeur
    if MIPHEUR is not None and MIPHEUR['enabled']:
        model.update()
        m_heur      = model.copy()
        # # Remove extra cycle constraints (more useful for m_heur than for main)
        # # (DON'T REMOVE, SHORT DFJ CYCLES ARE OK!)
        # cnstrs = model.getConstrs()
        # for cnstr in cnstrs:
        #     if 'CBC_MCB_FWD[' in cnstr.ConstrName:
        #         model.remove(cnstr)
        #     if 'CBC_MCB_BCK[' in cnstr.ConstrName:
        #         model.remove(cnstr)
        #     if 'CBC_CYC_FWD[' in cnstr.ConstrName:
        #         model.remove(cnstr)
        #     if 'CBC_CYC_BCK[' in cnstr.ConstrName:
        #         model.remove(cnstr)
        # model.update()        
        # --------------------------------------------------------
        m_heur._graph  = model._graph.copy() 
        m_heur._ntwrk  = model._ntwrk.copy()         
        if model._rtd is not None:                                
            m_heur._rtd = model._rtd.copy()
        else:
            m_heur._rtd = None
        m_heur._i2n    = model._i2n.copy()
        m_heur._n2i    = model._n2i.copy()        
        m_heur._i2a    = model._i2a.copy()
        m_heur._a2i    = model._a2i.copy()
        m_heur._i2e    = model._i2e.copy()
        m_heur._e2i    = model._e2i.copy()
        m_heur._trm2k  = model._trm2k.copy()
        m_heur._k2trm  = model._k2trm.copy() 
        m_heur._altrm  = model._altrm.copy()
        m_heur._troot  = model._troot
        m_heur._ttree  = model._ttree
        m_heur._relobj = model._relobj
        m_heur._heurst = False
        m_heur._Theur  = -1              
        m_heur._islsw  = model._islsw.copy() 
        m_heur.Params.PreCrush = model.Params.PreCrush
        m_heur.Params.LazyConstraints = model.Params.LazyConstraints
        #---------------------------------------------------------
        varbls         = m_heur.getVars()
        m_heur._x      = tupledict()
        m_heur._y      = tupledict()
        m_heur._z      = tupledict()
        m_heur._p      = tupledict()         
        m_heur._T      = tupledict()
        m_heur._pGS    = tupledict()
        m_heur._pLS    = tupledict()
        m_heur._PwB    = tupledict()
        m_heur._IMBLOD = tupledict()
        for var in varbls:
            if 'x[' in var.VarName:
                idx = re.findall(r'\d+', var.VarName)
                idx = tuple([int(i) for i in idx])
                m_heur._x[idx] = var
            elif 'y[' in var.VarName:
                idx = re.findall(r'\d+', var.VarName)
                idx = tuple([int(i) for i in idx])
                m_heur._y[idx] = var
            elif 'z[' in var.VarName:
                idx = re.findall(r'\d+', var.VarName)
                idx = tuple([int(i) for i in idx])
                m_heur._z[idx] = var
            elif 'pGS[' in var.VarName:
                idx = re.findall(r'\d+', var.VarName)
                idx = int(idx[0])
                m_heur._pGS[idx] = var
            elif 'IMBLOD[' in var.VarName:
                idx = re.findall(r'\d+', var.VarName)
                idx = int(idx[0])
                m_heur._IMBLOD[idx] = var                    
            elif 'pLS[' in var.VarName:
                idx = re.findall(r'\d+', var.VarName)
                idx = int(idx[0])
                m_heur._pLS[idx] = var
            elif 'p[' in var.VarName:
                idx = re.findall(r'\d+', var.VarName)
                idx = tuple([int(i) for i in idx])
                m_heur._p[idx] = var                            
            elif 'T[' in var.VarName:
                idx = re.findall(r'\d+', var.VarName)
                idx = tuple([int(i) for i in idx])
                m_heur._T[idx] = var                     
        #     elif 'z[' in var.VarName:
        #         m_heur.remove(var)            
        # cnstrs = m_heur.getConstrs()
        # for cnstr in cnstrs:
        #     if 'sptree' in cnstr.ConstrName:
        #         m_heur.remove(cnstr)
        #     if 'rootOUT' in cnstr.ConstrName:
        #         m_heur.remove(cnstr)   
        #     if 'nodeIN' in cnstr.ConstrName:
        #         m_heur.remove(cnstr)
        #     if 'validIJ' in cnstr.ConstrName:
        #         m_heur.remove(cnstr)        
        if len(m_heur._T)==0:
            m_heur._T = None
        x_lbnd = dict()
        for i,k in m_heur._x:
            x_lbnd[i,k] = m_heur._x[i,k].lb
        x_ubnd = dict()
        for i,k in m_heur._x:
            x_ubnd[i,k] = m_heur._x[i,k].ub            
        y_lbnd = dict()
        for i,j in m_heur._y:
            y_lbnd[i,j] = m_heur._y[i,j].lb
        y_ubnd = dict()
        for i,j in m_heur._y:
            y_ubnd[i,j] = m_heur._y[i,j].ub
        z_lbnd = dict()
        for i,j in m_heur._z:
            z_lbnd[i,j] = m_heur._z[i,j].lb
        z_ubnd = dict()
        for i,j in m_heur._z:
            z_ubnd[i,j] = m_heur._z[i,j].ub                      
        # To reset the original bounds after each run of the heuristic 
        m_heur._x_lbnd = x_lbnd
        m_heur._x_ubnd = x_ubnd        
        m_heur._y_lbnd = y_lbnd
        m_heur._y_ubnd = y_ubnd
        m_heur._z_lbnd = z_lbnd
        m_heur._z_ubnd = z_ubnd          
        m_heur.Params.LogToConsole = 0       
        m_heur.Params.TimeLimit = MIPHEUR['T_lim']   # max time per heuristic run
        model._mheur = m_heur  
        model._Theur = MIPHEUR['Theur']     # max time available for MIP heuristics
        model._heurst = MIPHEUR['enabled']
        model._xlast = None           
    else:
        model._heurst = False
    return model


def xtrct_filam(net,root_mlt,root_sgl):
    nods = list()
    grf = net.copy()
    for vi in grf.vs:
        # Identify a graph leaf node
        if vi.indegree()==1:
            nod = list()
            nod.append(vi.index);
            vj = grf.neighbors(vi)[0]
            nod.append(vj); 
            vj = grf.vs[vj]
            while vj.indegree()==2:
                v0 = vj.index
                vx = set([v for v in grf.neighbors(vj)]) - set(nod)
                assert(len(vx)==1)
                vx = vx.pop()
                nod.append(vx)
                vj = grf.vs[vx]
            nods.append(nod)
   
    # Remove complex and highly unlikely root combinations
    i_del = list()
    roots = root_mlt + root_sgl
    for i_nod,nod in enumerate(nods):
        vx_root = [vx  for vx in nod  if vx in roots]
        if len(vx_root)>1: 
            i_del.append(i_nod)        
    for i in sorted(i_del, reverse=True): 
        del nods[i]
    
    # Classify the remaining filaments
    flms_in = list()
    flmsout = list()         
    for i_nod,nod in enumerate(nods):
        ix_root_sgl = [ix  for ix,vx in enumerate(nod)  if vx in root_sgl]
        ix_root_mlt = [ix  for ix,vx in enumerate(nod)  if vx in root_mlt]
        assert(len(ix_root_sgl)<=1 and len(ix_root_mlt)<=1)
        if len(ix_root_sgl)==0 and len(ix_root_mlt)==0:
            flms_in.append(nod) 
        if len(ix_root_sgl)==1:
            ix_root_sgl = ix_root_sgl[0]
            # Don't add a single node with no "before" arcs
            if ix_root_sgl>0: flms_in.append(nod[0:ix_root_sgl+1])
        if len(ix_root_mlt)==1:
            ix_root_mlt = ix_root_mlt[0]
            # Don't add a single node with no "before" arcs            
            if ix_root_mlt>0: flms_in.append(nod[0:ix_root_mlt+1])
            flmsout.append(nod[ix_root_mlt:])
    
    # Output arcs
    arcs = list()        
    for filam in flms_in:
        arc = [(filam[ix+1],filam[ix]) for ix in range(len(filam)-1)]
        arcs.append(arc)
    for filam in flmsout:
        arc = [(filam[ix],filam[ix+1]) for ix in range(len(filam)-1)]
        arcs.append(arc)    
    
    return arcs, flms_in, flmsout, nods

            