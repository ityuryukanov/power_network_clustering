from gurobipy   import GRB, quicksum, tupledict
from bnc_xbase  import Bnc_Xbase
from top_callbk import dfjcut_callback
from utils import check_KVL, mst_cb, cycle_ordering, order_checking, xtrct_filam, write_sol
import pickle
import math
import os 

class Pnm_Cflow(Bnc_Xbase):
    # Implements the model based on artificial commodity flows.
    
    def __init__(self, args, RESLIM=None, OPTCR=0, MIPSTART=False, CONNTYP='flow'):
        
        mipstart = MIPSTART
        Bnc_Xbase.__init__(self, args, RESLIM, OPTCR, MIPSTART=mipstart)                
        
        # Create the graph partitioning model (initial solve)
        mod_gp = self.model
        num_vx = len(self.nodes)
        trm2k  = args['trm2k']
        n2i = self.model._n2i
        trm2k = dict((n2i[n],int(k)) for (n,k) in trm2k.items())
        mod_gp._kvlfl = False     # callback flags
        mod_gp._ztruu = False     # callback flags
        
        # Extra parameters/variables
        mod_gp._y = mod_gp.addVars(list(self.uarcs), name="y", vtype=GRB.BINARY, lb=0, ub=1)
        # Set initial values of z[i,j] to z_ini
        if MIPSTART:
            y_ini = args['y_ini']
            y_ini = dict(((mod_gp._n2i[i],mod_gp._n2i[j]), int(w)) for i,j,w in y_ini)            
            mod_gp.NumStart = 2; mod_gp.Params.StartNumber = 0;
            for i,j in y_ini: 
                mod_gp._y[i,j].start = y_ini[(i,j)] 
#---------------------Cycles (to check KVL at least)---------------------------
        grf_sp = mod_gp._grfsp
        assert(not grf_sp.is_directed())
        mstcb = self.mstcb
        cyc_pth = os.path.join(os.path.dirname(os.path.dirname( __file__ )), "cycledata", 'cyc_'+self.icase+'.pickle')
        fle_cyc = open(cyc_pth,'rb')        
        cyc = pickle.load(fle_cyc)         # cyc is already converted to igraph node indices
        self.cyc = set()                   # short cycles for initial DFJ constraints        
        self.cycmax = 7                    # maximum cycle length for cycle breaking constraints and other cycle constraints besides KVL
        self.srt = set()                   # extra "short" cycles for KVL
        self.cycmin = 3
        for loop in cyc:
            for i,j in loop: assert((i,j) in self.uarcs or (j,i) in self.uarcs)
            loop = [(i,j) if (i,j) in self.uarcs else (j,i) for i,j in loop]            
            if len(loop)<=self.cycmax:     # don't add too long cycles
                self.cyc.add(tuple(sorted(loop)))      
            if len(loop)<=self.cycmin:     # 
                self.srt.add(tuple(sorted(loop)))
        assert(len(mstcb)==len(grf_sp.es)-len(grf_sp.vs)+1)
#-----------------------------CONNECTIVITY-------------------------------------
        if CONNTYP=='flow':
            # Connectivity constraints through artificial commodity flows
            mod_gp._T = mod_gp.addVars(list(self.uarcs), name="T", vtype=GRB.CONTINUOUS, lb=-(num_vx-1), ub=(num_vx-1))
            mod_gp._NFanCo = tupledict()
            for r in mod_gp._islsw:
                Nr = set(mod_gp._graph.neighbors(r))
                k = trm2k[r]
                mod_gp._NFanCo[r] = mod_gp.addConstr(quicksum(mod_gp._T[r,j] for j in Nr if (r,j) in self.uarcs)-quicksum(mod_gp._T[j,r] for j in Nr if (j,r) in self.uarcs) == quicksum(mod_gp._x[i,k] for i in self.nodes)-1, "NFanCo1[%s]"%(r))
            nonsw = set(self.nodes) - set(mod_gp._islsw)
            for i in nonsw:        
                Ni = set(mod_gp._graph.neighbors(i))
                mod_gp._NFanCo[i] = mod_gp.addConstr(quicksum(mod_gp._T[i,j] for j in Ni if (i,j) in self.uarcs)-quicksum(mod_gp._T[j,i] for j in Ni if (j,i) in self.uarcs) == -1, "NFanCo2[%s]" % (i))
            # Switching & Commodity flow limitation
            mod_gp._NFanHI = tupledict()
            mod_gp._NFanLO = tupledict()
            for i,j in self.uarcs:
                mod_gp._NFanHI[i,j] = mod_gp.addConstr(mod_gp._T[i,j]<= (num_vx-1)*(1-mod_gp._y[i,j]), "NFanCo3[%s,%s]" % (i,j))
                mod_gp._NFanLO[i,j] = mod_gp.addConstr(mod_gp._T[i,j]>=-(num_vx-1)*(1-mod_gp._y[i,j]), "NFanCo4[%s,%s]" % (i,j))        
        else:
            mod_gp._ztruu = True         # callback flags
            # Connectivity constraints through spanning forests
            self.arcs = [(mod_gp._n2i[i],mod_gp._n2i[j]) for i,j in self.arcs]
            mod_gp._z = mod_gp.addVars(list(self.arcs), name="z", vtype=GRB.BINARY, lb=0, ub=1)
            arccyc = {(i,j) if (i,j) in self.uarcs else (j,i) for cyc in mstcb for i,j in cyc}
            arcarb = set(self.uarcs) - arccyc     
            assert(len(arccyc)+len(arcarb)==len(self.uarcs))
            # Set initial values of z[i,j] to z_ini
            if MIPSTART:
                z_ini = args['z_ini']
                z_ini = dict(((mod_gp._n2i[i],mod_gp._n2i[j]), int(w)) for i,j,w in z_ini)            
                for i,j in z_ini:
                    mod_gp._z[i,j].start = z_ini[(i,j)]        
            # Spanning tree cardinality constraint
            mod_gp.addConstr(quicksum([mod_gp._z[i,j] for i,j in self.arcs])==len(self.nodes)-len(self.islid), "sptree")
            # Root nodes have no incoming arcs
            for r in mod_gp._islsw:
                Nr = set(mod_gp._graph.neighbors(r))            
                for j in Nr:
                    mod_gp._z[j,r].lb = 0   #rootIN[%s]
                    mod_gp._z[j,r].ub = 0   
            # Root nodes of MULTITERMINAL groups have at least one outgoing arc
            assert(mod_gp._rtd is not None)
            root_mlt = [r for r in mod_gp._rtd if len(mod_gp._rtd[r])>0]
            root_sgl = list(set(mod_gp._islsw) - set(root_mlt))
            for r in root_mlt:
                Nr = set(mod_gp._graph.neighbors(r))
                mod_gp.addConstr(quicksum([mod_gp._z[r,j] for j in Nr])>=1, "rootOUT[%s]" % (r))
            # Non-root nodes have strictly one incoming arc
            nonsw = set(self.nodes) - set(mod_gp._islsw)
            for v in nonsw:
                Nv = set(mod_gp._graph.neighbors(v))
                mod_gp.addConstr(quicksum([mod_gp._z[j,v] for j in Nv])==1, "nodeIN[%s]" % (v))
            # Fix z_ij for network filaments:
            arcsfil, fils_in, filsout, nods = xtrct_filam(mod_gp._graph,root_mlt,root_sgl)
            arcsfix = [arc  for fil in arcsfil for arc in fil]
            for i,j in arcsfix:
                assert(((i,j) in arcarb) or ((j,i) in arcarb))
                mod_gp._z[i,j].lb = 1
                mod_gp._z[i,j].ub = 1
                mod_gp._z[j,i].lb = 0
                mod_gp._z[j,i].ub = 0
            # Relationship between y and z
            for i,j in arccyc:
                mod_gp.addConstr(mod_gp._z[i,j]+mod_gp._z[j,i] <= 1-mod_gp._y[i,j], "validIJ[%s,%s]" % (i,j))
            for i,j in arcarb:
                mod_gp.addConstr(mod_gp._z[i,j]+mod_gp._z[j,i] == 1-mod_gp._y[i,j], "validIJ[%s,%s]" % (i,j))
            # Cycle breaking constraints for available cycles (KEEP IT, AS SHORT 
            # CYCLES HELP THE HEURISTICS AND THE SOLUTION ITSELF)
            len_lim = self.cycmax  
            for icyc,cyc in enumerate(self.cyc):
                if len(cyc)<=len_lim:
                    cyc_ord = cycle_ordering(cyc)
                    arc_ord = order_checking(cyc_ord)
                    cycle = {v for arc in arc_ord for v in arc}
                    mod_gp.addConstr(quicksum(mod_gp._z[arc[0],arc[1]] for arc in arc_ord)<=len(cycle)-1, "CBC_CYC_FWD[%s]" % (icyc))
                    mod_gp.addConstr(quicksum(mod_gp._z[arc[1],arc[0]] for arc in arc_ord)<=len(cycle)-1, "CBC_CYC_BCK[%s]" % (icyc))        
            # Definition of SOS1 variables
            for i,j in self.uarcs:
                mod_gp.addSOS(GRB.SOS_TYPE1, [mod_gp._z[i,j],mod_gp._z[j,i]])                
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
        # Switching constraints
        for i,j in self.uarcs:
            for k in self.islid:
                mod_gp.addConstr(mod_gp._x[i,k]-mod_gp._x[j,k] <= mod_gp._y[i,j], "lblXij1[%s,%s,%s]" % (i,j,k))
                mod_gp.addConstr(mod_gp._x[j,k]-mod_gp._x[i,k] <= mod_gp._y[i,j], "lblXij2[%s,%s,%s]" % (j,i,k))
                mod_gp.addConstr(mod_gp._x[i,k]+mod_gp._x[j,k]-1+mod_gp._y[i,j] <= 1, "SWITCHG[%s,%s,%s]" % (i,j,k))   #line switching inside islands is OK                
        
        # Switching & Power flow limitation
        for i,j in self.uarcs:
            if (i,j) in self.brlim and self.brlim[(i,j)]>0:
                mod_gp.addConstr(mod_gp._p[i,j]<= self.brlim[(i,j)]*(1-mod_gp._y[i,j]), "LinLimU[%s,%s]" % (i,j))
                mod_gp.addConstr(mod_gp._p[i,j]>=-self.brlim[(i,j)]*(1-mod_gp._y[i,j]), "LinLimL[%s,%s]" % (i,j))        
        # Power flows & Phase differences
        mod_gp._AngLimU = tupledict()
        mod_gp._AngLimL = tupledict()
        for i,j in self.uarcs:
            mod_gp._AngLimU[i,j] = mod_gp.addConstr(mod_gp._p[i,j]-self.bdcpf[(i,j)]*(mod_gp._ph[i]-mod_gp._ph[j])<= 2*math.pi*self.bdcpf[(i,j)]*mod_gp._y[i,j], "AngLimU[%s,%s]" % (i,j))                
            mod_gp._AngLimL[i,j] = mod_gp.addConstr(mod_gp._p[i,j]-self.bdcpf[(i,j)]*(mod_gp._ph[i]-mod_gp._ph[j])>=-2*math.pi*self.bdcpf[(i,j)]*mod_gp._y[i,j], "AngLimL[%s,%s]" % (i,j))
        
        # Objective
        alpha = 0.00
        betta = 1.00
        gamma = 0.01
        muuuu = self.mu
        mod_gp.setObjective(
          + alpha*quicksum(mod_gp._IMBLOD[k] for k in self.islid)\
          + betta*quicksum(mod_gp._pLS[l] for l in self.loads) \
          + gamma*quicksum(mod_gp._pGS[g] for g in self.genss) \
          + muuuu*quicksum(self.flow0[(i,j)]*mod_gp._y[i,j] for i,j in self.uarcs), GRB.MINIMIZE)
        self.model = mod_gp
        self.cntyp = CONNTYP
        mod_gp._nonsw = nonsw
        mod_gp._Tfeas = -1     # time to feasibility
    
    def solve(self):
        self.model.update()      
        if self.cntyp=='flow':
            self = Bnc_Xbase.solve(self, dfjcut_callback)  # callback is for MIP heuristics, recording Tfeas, changing solution strategy...
        elif self.cntyp=='tree':
            self = Bnc_Xbase.solve(self, dfjcut_callback)
        else:
            raise RuntimeError        
        self.totPLS = sum([self.model._pLS[i].x  for i in self.loads])
        print("TOTAL LOAD SHEDDING IS: " + str(self.totPLS))
        self.totPGS = 0
        for i in self.genss:
            self.totPGS += self.model._pGS[i].x
        print("TOTAL GENERATION SHEDDING IS: " + str(self.totPGS))
        self.totIMB = 0
        for k in self.islid:
            self.totIMB += self.model._IMBLOD[k].x
        print("TOTAL LOAD GENERATION IMBALANCE IS: " + str(self.totIMB))
        self.islIMB = dict((k,0) for k in self.islid)
        for i in self.nodes:
            for k in self.islid:
                if self.model._x[i,k].x>0.99:
                    self.islIMB[k]+=self.pL0[i]-self.pG0[i]
        print("TOTAL POWER IMBALANCE IS: " + str(sum([abs(self.islIMB[k]) for k in self.islid])))
        self.totCUT = 0
        y_sol = dict()
        p_sol = dict()
        for i,j in self.uarcs:
            y_ij = self.model._y[i,j].x
            p_ij = self.model._p[i,j].x
            y_sol[(i,j)] = y_ij
            p_sol[(i,j)] = p_ij
            self.totCUT += self.flow0[(i,j)]*round(y_ij)
        print("TOTAL POWER FLOW CUT IS: " + str(self.totCUT))

        # Check cycles KVL
        grf_sp = self.model._grfsp
        assert(not grf_sp.is_directed())
        e2i = self.model._e2i
        nc = len(self.islid)
        a_del = [e2i[ij] for ij,w in y_sol.items() if abs(w)>0.98]
        grf_sp.delete_edges(a_del)
        esIDs = grf_sp.spanning_tree(return_tree=False, weights=None)
        mstcb = mst_cb(grf_sp, esIDs)
        mstcb = [tuple([(i,j) if (i,j) in self.uarcs else (j,i) for i,j in cyc]) for cyc in mstcb]        
        assert(len(mstcb)==len(grf_sp.es)-len(grf_sp.vs)+nc)
        cycls = set(mstcb).union(self.cyc).union(self.srt)
        valid_cycls = set()
        for cyc in cycls:
            n_sw = [1 for i,j in cyc if y_sol[i,j]<0.02]
            if len(n_sw)==0:
                valid_cycls.add(tuple(sorted(cyc)))
        cycfail = check_KVL(valid_cycls,p_sol,self.bdcpf)[0]
        if len(cycfail)>0:
            print("KVL is not OK!")
            return self
        else:
            print("KVL is OKAY!")
            
        print("TOTAL RUNTIME IS: " + str(self.model.Runtime))
        if hasattr(self.model,'_Theur'):
            T_heur = self.model._MipHeur['Theur'] - self.model._Theur
            print("TOTAL TIME SPENT IN HEURISTIC IS: " + str(T_heur))
        else:
            print("TOTAL TIME SPENT IN HEURISTIC IS: -1")
        # Extract the solution     
        self.GAP    = abs(self.model.ObjVal-self.model.ObjBound)/abs(self.model.ObjVal)
        self.islIMB = sum([abs(self.islIMB[k]) for k in self.islid])            
        ici_pth = os.path.abspath(os.path.join(os.path.dirname( __file__ ), 'ici'))
        write_sol(self.model, 'p',   [True, True], ici_pth)
        write_sol(self.model, 'x',   [True, False],ici_pth)
        write_sol(self.model, 'y',   [True, True], ici_pth)
        write_sol(self.model, 'pGS', [True], ici_pth)
        write_sol(self.model, 'pLS', [True], ici_pth)
        with open(os.path.join(ici_pth, "runtime.csv"), "w", newline='') as f:
            f.write(str(self.runtime))
        with open(os.path.join(ici_pth, "UB.csv"), "w", newline='') as f:
            f.write(str(self.model.ObjVal))
        with open(os.path.join(ici_pth, "GAP.csv"), "w", newline='') as f:
            f.write(str(self.GAP))
        with open(os.path.join(ici_pth, "islIMB.csv"), "w", newline='') as f:
            f.write(str(self.islIMB))
        with open(os.path.join(ici_pth, "totIMB.csv"), "w", newline='') as f:
            f.write(str(self.totIMB))
        with open(os.path.join(ici_pth, "totPLS.csv"), "w", newline='') as f:
            f.write(str(self.totPLS))
        with open(os.path.join(ici_pth, "totPGS.csv"), "w", newline='') as f:
            f.write(str(self.totPGS))
        with open(os.path.join(ici_pth, "totCUT.csv"), "w", newline='') as f:
            f.write(str(self.totCUT))
        with open(os.path.join(ici_pth, "feastime.csv"), "w", newline='') as f:
            f.write(str(self.model._Tfeas))            
        return self


