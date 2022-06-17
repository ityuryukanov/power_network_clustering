from bnc_dcopf import Bnc_Dcopf
from utils import write_sol, mst_cb
from gurobipy import GRB, quicksum 
import os

class Bnc_Xbase(Bnc_Dcopf):
    
    def __init__(self, args, RESLIM=None, OPTCR=0, MIPSTART=False):        
        
        Bnc_Dcopf.__init__(self, args, RESLIM, OPTCR)                   
        x_ini = args['x_ini']        
        trm2k = args['trm2k']        
        islid = self.islid
        n2i = self.model._n2i
        i2n = self.model._i2n
        net = self.model._graph
        grf = self.model._ntwrk

        # Convert all data from labels to indices        
        x_ini = dict(((n2i[i],int(k)), int(v)) for i,k,v in x_ini)
        trm2k = dict((n2i[n],int(k)) for (n,k) in trm2k.items()) 
        
        # Introduce x as binary SOS1 variables        
        x = self.model.addVars(self.nodes, islid, name="x", vtype=GRB.BINARY, lb=0, ub=1)
        for i in trm2k.keys():
            x[i,trm2k[i]].lb = 1
            x[i,trm2k[i]].ub = 1
        # Defines x to be a SOS1 variable
        for i in self.nodes:
            self.model.addSOS(GRB.SOS_TYPE1, [x[i,k] for k in islid])
        # Defines x[i,k] to have exactly single "1" per node
        for i in self.nodes:
            self.model.addConstr(quicksum([x[i,k] for k in islid])== 1, "DefSOSx[%s]" % (i2n[i]))
        # Set initial values of x[i,k] to x_ini
        if MIPSTART:
            self.model.Params.StartNumber = 0         
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
        self.model._trm2k = trm2k
        self.model._k2trm = k2trm
        self.model._altrm = altrm      
        self.model._rtd = rtd
        self.model._islsw = list(rtd.keys())
        self.model._grfsp = grf_sp   # graph of MST cycle basis
        self.mstcb = mstcb                  
        
                        
    def solve(self, callbackfcn):                         
        self = Bnc_Dcopf.solve(self, callbackfcn)
        ici_pth = os.path.abspath(os.path.join(os.path.dirname( __file__ ), 'ici'))
        write_sol(self.model, 'x', [True, False], ici_pth)
        return self              


