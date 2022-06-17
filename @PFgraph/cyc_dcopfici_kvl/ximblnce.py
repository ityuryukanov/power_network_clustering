import os
import re
from cycopf import Cycopf
from utils import write_sol, mst_cb
from gurobipy import GRB, quicksum, tupledict 

class Ximblnce(Cycopf):
    
    def __init__(self, args, RESLIM=None, OPTCR=0, MIPSTART=False):        
        
        Cycopf.__init__(self, args, RESLIM, OPTCR)                   
        x_ini = args['x_ini']
        trm2k = args['trm2k']
        islid = self.islid
        n2i = self.model._n2i
        i2n = self.model._i2n
        net = self.model._graph

        # Convert all data from labels to indices        
        x_ini = dict(((n2i[i],int(k)), int(v)) for i,k,v in x_ini)
        trm2k = dict((n2i[n],int(k)) for (n,k) in trm2k.items()) 
        
        # Defines x to be a SOS1 variable
        x = self.model.addVars(self.nodes, islid, name="x", vtype=GRB.BINARY, lb=0, ub=1)
        for i in self.nodes:
            self.model.addSOS(GRB.SOS_TYPE1, [x[i,k] for k in islid])
        # Defines x[i,k] to have exactly single "1" per node
        for i in self.nodes:
            self.model.addConstr(quicksum([x[i,k] for k in islid])== 1, "DefSOSx[%s]" % (i2n[i]))
        # Set initial values of x[i,k] to x_ini
        if MIPSTART:
            self.model.NumStart = 2; self.model.Params.StartNumber = 0;                     
            for i,k in x_ini:
                x[i,k].start = x_ini[(i,k)]
        
        # Form a dictionary from group ids to corresponding terminal nodes
        k2trm = dict()
        for i in trm2k:
            if trm2k[i] not in k2trm:
                k2trm[trm2k[i]] = []
                k2trm[trm2k[i]].append(i)
            else:
                k2trm[trm2k[i]].append(i)
                
        # Form the list of all terminals to save time
        altrm = set([vx for k in k2trm.keys() for vx in k2trm[k]])        
        
        # Find farthest terminals in each group and set them as tree roots         
        grf = self.model._ntwrk
        grf_sp = grf.copy()
        es_wgt = len(grf_sp.es)*[1]                
        alllen = grf_sp.shortest_paths(source=None, target=None, weights=es_wgt, mode='out')
        rtd = dict()    # root sets        
        for k in k2trm:
            trm_curr = set(k2trm[k])
            trm_rest = altrm - trm_curr
            trm_curr = list(trm_curr)
            trm_rest = list(trm_rest)
            avg_dist = list()
            for f in trm_curr:
                dst_rest = [alllen[f][t] for t in trm_rest]
                dst_curr = [alllen[f][t] for t in trm_curr]
                avg_dist.append(min(dst_rest)+0.1*sum(dst_rest)/len(trm_rest) + 0.02*sum(dst_curr)/len(dst_curr))
            max_dist = max(avg_dist)
            idx_maxx = avg_dist.index(max_dist)
            rootcurr = trm_curr[idx_maxx] 
            rtd[rootcurr] = list()
            for t in trm_curr:
                if abs(alllen[rootcurr][t])>1e-6:
                    rtd[rootcurr].append((t,alllen[rootcurr][t]))
            
        # Find basic MST cycle bases for each generator group
        esIDs = grf_sp.spanning_tree(return_tree=False, weights=None)        
        mstcb = mst_cb(grf_sp, esIDs)
        assert(len(mstcb)==len(grf_sp.es)-len(grf_sp.vs)+1)
        mstcb = [tuple([(i,j) if (i,j) in self.uarcs else (j,i) for i,j in cyc]) for cyc in mstcb] 
        for es in grf_sp.es:
            if es.index in esIDs:
                es['weight'] = 1e-6;
            else:
                es['weight'] = 1;

        # Minimize worst power imbalance (which must be positive: at least 1 load-rich island)
        ub_IMB = max(+sum(self.pG0[i] for i in self.nodes), sum(self.pL0[i] for i in self.nodes))
        IMBLOD = self.model.addVars(self.islid, name="IMBLOD", vtype=GRB.CONTINUOUS, lb=0, ub=ub_IMB)
        for k in self.islid:
            gen_k = k2trm[k]
            pwr_k = [self.pG0[i] for i in gen_k]
            assert(all([pwr>=-1e-12 for pwr in pwr_k]))
            self.model.addConstr(IMBLOD[k]>= quicksum((self.pL0[i]-self.pG0[i])*x[i,k] for i in self.nodes), "IMBPOS[%s]" % k) 
            self.model.addConstr(IMBLOD[k]>=-quicksum((self.pL0[i]-self.pG0[i])*x[i,k] for i in self.nodes), "IMBNEG[%s]" % k) 
        
        self.model.update()      
        self.model._x = x
        self.model._IMBLOD = IMBLOD
        self.model._relobj = list()
        self.model._imprv = 0
        self.model._troot = 0
        self.model._ttree = 0
        self.model.Params.PreCrush = 1
        self.model.Params.LazyConstraints = 1
        self.model._trm2k = trm2k
        self.model._k2trm = k2trm
        self.model._altrm = altrm
        self.model._rtd = rtd
        self.model._islsw = list(rtd.keys())
        self.model._grfsp = grf_sp   # graph of MST cycle basis
        self.mstcb = mstcb
        
    
    def conf_modgp(self):
        """
        Configure the auxiliary modgp model (pure graph partitioning without 
        power flow constraints). It might be useful for finding initial feasible
        solutions, but it is not used in the final versions of the model.  
        """        
        modOPF = self.model
        modOPF.update()    
        mod_gp = modOPF.copy()    
        mod_gp._graph  = modOPF._graph.copy() 
        mod_gp._ntwrk  = modOPF._ntwrk.copy()         
        if modOPF._rtd is not None:                                
            mod_gp._rtd = modOPF._rtd.copy()
        else:
            mod_gp._rtd = None
        mod_gp._i2n    = modOPF._i2n.copy()
        mod_gp._n2i    = modOPF._n2i.copy()        
        mod_gp._i2a    = modOPF._i2a.copy()
        mod_gp._a2i    = modOPF._a2i.copy()
        mod_gp._i2e    = modOPF._i2e.copy()
        mod_gp._e2i    = modOPF._e2i.copy()
        mod_gp._trm2k  = modOPF._trm2k.copy()
        mod_gp._k2trm  = modOPF._k2trm.copy() 
        mod_gp._altrm  = modOPF._altrm.copy()
        mod_gp._troot  = modOPF._troot
        mod_gp._ttree  = modOPF._ttree
        mod_gp._relobj = modOPF._relobj
        mod_gp._heurst = modOPF._heurst      # OK if auxiliary model copies the original       
        mod_gp._Theur  = modOPF._Theur       # max time available for MIP heuristics
        mod_gp._MipHeur = modOPF._MipHeur.copy()               
        mod_gp._islsw  = modOPF._islsw.copy() 
        mod_gp.Params.PreCrush = modOPF.Params.PreCrush
        mod_gp.Params.LazyConstraints = modOPF.Params.LazyConstraints                   
        mod_gp._x      = tupledict() 
        mod_gp._y      = tupledict() 
        mod_gp._z      = tupledict() 
        mod_gp._l      = tupledict()
        mod_gp._p      = tupledict()
        mod_gp._pLS    = tupledict()
        mod_gp._pGS    = tupledict()
        mod_gp._IMBLOD = tupledict()
        varbls = mod_gp.getVars()      
        for var in varbls:
            if 'x[' in var.VarName:
                idx = re.findall(r'\d+', var.VarName)
                idx = tuple([int(i) for i in idx])
                mod_gp._x[idx] = var            
            elif 'y[' in var.VarName: 
                idx = re.findall(r'\d+', var.VarName)
                idx = tuple([int(i) for i in idx])
                mod_gp._y[idx] = var
            elif 'z[' in var.VarName: 
                idx = re.findall(r'\d+', var.VarName)
                idx = tuple([int(i) for i in idx])
                mod_gp._z[idx] = var                
            elif 'p[' in var.VarName: 
                idx = re.findall(r'\d+', var.VarName)
                idx = tuple([int(i) for i in idx])
                mod_gp._p[idx] = var                  
            elif 'l[' in var.VarName: 
                nme = var.VarName
                idx = nme[1:]
                mod_gp._l[idx] = var                                
            elif 'pLS[' in var.VarName:
                idx = re.findall(r'\d+', var.VarName)
                idx = int(idx[0])
                mod_gp._pLS[idx] = var           
            elif 'pGS[' in var.VarName:
                idx = re.findall(r'\d+', var.VarName)
                idx = int(idx[0])
                mod_gp._pGS[idx] = var
            elif 'IMBLOD' in var.VarName:
                idx = re.findall(r'\d+', var.VarName)
                idx = int(idx[0])
                mod_gp._IMBLOD[idx] = var       
        return mod_gp
        
    def solve(self, callbackfcn=None):                         
        self = Cycopf.solve(self, callbackfcn)
        ici_pth = os.path.abspath(os.path.join(os.path.dirname( __file__ ), 'ici'))
        write_sol(self.model, 'x', [True, False], ici_pth)        
        return self 
