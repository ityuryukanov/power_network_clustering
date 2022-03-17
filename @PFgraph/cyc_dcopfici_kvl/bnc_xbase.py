from bnc_dcopf import Bnc_Dcopf
from utils import write_sol
from gurobipy import GRB, quicksum 
import os

class Bnc_Xbase(Bnc_Dcopf):
    
    def __init__(self, args, RESLIM=None, OPTCR=0, DIGRAPH=False, MIPSTART=False, MIPHEUR = None):        
        
        mipheurt = MIPHEUR
        Bnc_Dcopf.__init__(self, args, RESLIM, OPTCR, DIGRAPH, MIPHEUR = mipheurt)                   
        x_ini = args['x_ini']        
        trm2k = args['trm2k']        
        islid = self.islid
        n2i = self.model._n2i
        i2n = self.model._i2n
        net = self.model._graph

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
        ### Fixes x[i,k] to x_ini (for debug only)
        ##for i,k in x_ini:
        ##    self.model.addConstr(x[i,k] == x_ini[(i,k)], "FixXINI[%s]" % (i2n[i]))
        
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
                
        # # For terminals in non-singleton groups, at least one "credible" node neighbour should be in the group
        # for k in k2trm:
        #     grp = set(k2trm[k])
        #     if len(grp)<2: continue
        #     for t in grp:
        #         Nt = set(net.neighbors(t))
        #         if len(Nt.intersection(grp))>0: 
        #             continue   # terminal is already connected; forcing extra connections is wrong due to RHS==1 in TermNbr[%s,%s]
        #         At = Nt.copy()
        #         Ct = Nt.copy()
        #         At.add(t)                
        #         for v in Nt:
        #             Av = set(net.neighbors(v))
        #             if Av.issubset(At):  
        #                 Ct.remove(v)       
        #         self.model.addConstr(quicksum([x[i,k] for i in Ct])>=1, "TermNbr[%s,%s]" % (i2n[t],k))    
        
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
        
                        
    def solve(self, callbackfcn):                         
        self = Bnc_Dcopf.solve(self, callbackfcn)
        ici_pth = os.path.abspath(os.path.join(os.path.dirname( __file__ ), 'ici'))
        write_sol(self.model, 'x', [True, False], ici_pth)
        return self              


