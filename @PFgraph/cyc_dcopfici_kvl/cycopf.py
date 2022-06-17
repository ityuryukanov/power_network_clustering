from gurobipy import Model, GRB, quicksum, tupledict
from utils import write_sol, check_KVL, mst_cb, opf_lbl2idx_grf
import pickle
import igraph
import os 
import re
import sys
sys.path.append(os.path.dirname(os.path.dirname( __file__ )))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..\\..\\')))  #for dbg
import igraphmod    #for dbg
from grb_read_graph import check_graph

class Cycopf(object):
    """
    Define a dc network flow model (no voltage angles).
    The parameter epsL adds an extra tiny load at each load node. This may improve 
    convergence by reducing the number of required lazy connectivity DFJ onstraints.
    However, a non-zero epsL may result in infeasibilities that were not present in 
    the original problem (but this typically only happens with peculiar study cases  
    such as the IEEE 24 bus reliability test system)...    
    """
    
    def __init__(self, args, RESLIM=None, OPTCR=0):
		
        uarcs, arcs = check_graph(args['bdcpf'],args['nodes'],args['genss'],args['loads'])
        islid = [int(k) for k in args['islid']]
        islid.sort(reverse=False)
        self.bdcpf = args['bdcpf']
        self.pL0   = args['pL0']
        self.pG0   = args['pG0']
        self.flow0 = args['flow0']
        self.brlim = args['brlim']
        self.nodes = args['nodes']
        self.genss = args['genss']
        self.loads = args['loads']
        self.costL = args['costL']
        self.pGmin = args['pGmin']
        self.pGmax = args['pGmax']
        self.pLmin = args['pLmin']
        self.pLmax = args['pLmax']
        self.icase = args['icase']
        self.mu    = args['mu']
        self.nu    = args['nu']
        model = Model('bnc_dcopf')
        model.Params.MIPGap = OPTCR
        if RESLIM: model.Params.TimeLimit = RESLIM
        
        # Prepare the common DC OPF and graph partitioning parameters:
        self,net,grf,i2n,n2i,i2a,a2i,e2i,i2e,arcs,uarcs,lo_LS,up_LS,lo_GS,up_GS,brLlo,brLup = opf_lbl2idx_grf(self, arcs, uarcs)

        # Introduce p as line power flows in each uarc        
        p = model.addVars(uarcs, name="p", vtype=GRB.CONTINUOUS, lb=brLlo, ub=brLup)   
        # Introduce pLS as load shedding at load nodes
        pLS = model.addVars(self.loads, name="pLS", vtype=GRB.CONTINUOUS, lb=lo_LS, ub=up_LS)
        # Introduce pGS as generator shedding at generator nodes
        pGS = model.addVars(self.genss, name="pGS", vtype=GRB.CONTINUOUS, lb=lo_GS, ub=up_GS)

        # Power balance        
        epsL = 0;   # 5e-3
        for i in self.nodes:
            if i in self.genss and i in self.loads:
                model.addConstr(p.sum(i,'*')-p.sum('*',i)==-epsL-self.pL0[i]+pLS[i]+self.pG0[i]-pGS[i], "PwBal[%s]" % (i))
            elif i in self.genss:    
                model.addConstr(p.sum(i,'*')-p.sum('*',i)==self.pG0[i]-pGS[i], "PwBal[%s]" % (i))
            elif i in self.loads:
                model.addConstr(p.sum(i,'*')-p.sum('*',i)==-epsL-self.pL0[i]+pLS[i], "PwBal[%s]" % (i))
            else:
                model.addConstr(p.sum(i,'*')-p.sum('*',i)==-epsL, "PwBal[%s]" % (i))
                
        # # Introduce ph as bus voltage angle difference over each branch at each bus
        # ph = model.addVars(self.nodes, name="ph", vtype=GRB.CONTINUOUS, lb=-math.pi, ub=math.pi)
        # for i in islsw:     
        #     ph[i].lb = 0
        #     ph[i].ub = 0                
        
        # Adjust "initial power flows weights"    
        for i,j in uarcs:
            if (i,j) not in self.flow0:
                self.flow0[(i,j)]=0  
				
        # Load precomputed fundamental cycles
        mcb_pth = os.path.join(os.path.dirname(os.path.dirname( __file__ )), "cycledata", 'mcb_'+self.icase+'.pickle')
        try:
            fle_mcb = open(mcb_pth,'rb')
            mcb = pickle.load(fle_mcb)   # mcb is already converted to igraph node indices
            self.mcb = set()
            for loop in mcb:
                for i,j in loop:
                    assert((i,j) in uarcs or (j,i) in uarcs)
                loop = [(i,j) if (i,j) in uarcs else (j,i) for i,j in loop]
                self.mcb.add(tuple(sorted(loop)))
            lll = [len(loop) for loop in self.mcb]   # loop lengths
        except FileNotFoundError:
            self.mcb = set()
        cyc_pth = os.path.join(os.path.dirname(os.path.dirname( __file__ )), "cycledata", 'cyc_'+self.icase+'.pickle')
        fle_cyc = open(cyc_pth,'rb')        
        cyc = pickle.load(fle_cyc)   # cyc is already converted to igraph node indices
        self.cyc = set()   # short cycles for initial DFJ constraints        
        self.cycmax = 7    # maximum cycle length for cycle breaking constraints and other cycle constraints besides KVL
        self.srt = set()   # extra "short" cycles for KVL
        self.cycmin = 3
        for loop in cyc:
            for i,j in loop: assert((i,j) in uarcs or (j,i) in uarcs)
            loop = [(i,j) if (i,j) in uarcs else (j,i) for i,j in loop]            
            if len(loop)<=self.cycmax:     # don't add too long cycles
                self.cyc.add(tuple(sorted(loop)))      
            if len(loop)<=self.cycmin:     # 
                self.srt.add(tuple(sorted(loop)))
        
        # Store the variables
        model._p = p
        model._pLS = pLS
        model._pGS = pGS
        #model._ph  = ph        
        # BnC information             
        model._graph = net
        model._ntwrk = grf
        model._i2n = i2n
        model._n2i = n2i
        model._i2a = i2a
        model._a2i = a2i
        model._e2i = e2i       
        model._i2e = i2e
        model.NumStart = 2 		
        model._ztruu = False  # default values of callback flags
        model._kvlfl = False  # default values of callback flags		
        # Finalize self
        self.arcs  = arcs
        self.uarcs = uarcs
        self.model = model            
        self.islid = islid                         
                        
    def solve(self, callbackfcn=None):
        model = self.model
        model.Params.MIPFocus = 3
        model._changeParam = False  # solution strategy change inside of MIP callback
        if model.isMIP == 0:
            print('Model is not a MIP')
            exit(0)
        try:
            if callbackfcn is not None:
                model.optimize(callbackfcn)  #callback can be used to change the solution strategy or for MIP heuristics
            else:
                model.Params.MIPFocus = 0
                model.optimize()
        except GurobiError:
            print(GurobiError.message)
        t_st1 = None
        if model._changeParam:
            model.Params.MIPFocus = 0
            model.Params.TimeLimit -= model.Runtime
            t_st1 = model.Runtime   # solution time of the first stage
            try:
                if callbackfcn is not None:
                    model.optimize(callbackfcn)
                else:
                    raise RuntimeError
            except GurobiError:
                print(GurobiError.message)
        
        if model.status==GRB.Status.INF_OR_UNBD or model.status==GRB.Status.INFEASIBLE or model.status==GRB.Status.UNBOUNDED:
            model.computeIIS()
            model.write("test.ilp")
            return None 
        if t_st1 is None:
            t_st1 = 0
        
        self.optimal    = model.Status == GRB.OPTIMAL
        self.runtime    = model.Runtime + t_st1
        self.node_count = model.nodecount
        self.mip_gap    = model.mipgap
        self.objective  = model.ObjVal
        self.totLS = sum([self.model._pLS[i].x  for i in self.loads])
        self.cstLS = self.nu*sum([self.costL[l]*self.model._pLS[l].x for l in self.loads])
        MIPHEUR = model._MipHeur
        if MIPHEUR is not None and MIPHEUR['enabled']:
            T_heur = MIPHEUR['Theur'] - model._Theur
            print("TOTAL TIME SPENT IN HEURISTICS IS: " + str(T_heur))
        else:
            print("TOTAL TIME SPENT IN HEURISTICS IS: -1")
        print("TOTAL TIME SPENT IN MIPSOL CALLBACKS IS: " + str(model._ttree))
        print("TOTAL LOAD SHEDDING IS: " + str(self.totLS))
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
            y_sol[(i,j)] = y_ij;
            p_sol[(i,j)] = p_ij;
            self.totCUT += self.flow0[(i,j)]*round(y_ij)
        print("TOTAL POWER FLOW CUT IS: " + str(self.totCUT)) 
        # Check cycles KVL
        grf_sp = model._ntwrk.copy()  # assign sign to each cycle edge
        e2i = model._e2i
        nc = len(self.islid)
        a_del = [e2i[ij] for ij,w in y_sol.items() if abs(w)>0.98]
        grf_sp.delete_edges(a_del)
        esIDs = grf_sp.spanning_tree(return_tree=False, weights=None)
        mstcb = mst_cb(grf_sp, esIDs)
        assert(len(mstcb)==len(grf_sp.es)-len(grf_sp.vs)+nc)  
        cycls = set(mstcb).union(self.cyc).union(self.srt).union(self.mcb) 
        valid_cycls = set()
        for cyc in cycls:
            n_sw = [1 for i,j in cyc if y_sol[i,j]>0.98]
            if len(n_sw)==0:
                valid_cycls.add(tuple(sorted(cyc)))
        cycfail = check_KVL(valid_cycls,p_sol,self.bdcpf)[0]
        if len(cycfail)>0:
            print("KVL is not OK!")
            return self
        else:
            print("KVL is OKAY!")
        print("The total number of ìnserted cycles: "+str(model._cbkcyc))
       
        # Extract the solution  
        self.GAP = abs(self.model.ObjVal-self.model.ObjBound)/abs(self.model.ObjVal)
        self.islIMB = sum([abs(self.islIMB[k]) for k in self.islid])
        ici_pth = os.path.abspath(os.path.join(os.path.dirname( __file__ ), 'ici'))
        write_sol(model, 'y', [True, True], ici_pth)
        write_sol(model, 'p', [True, True], ici_pth)
        write_sol(model, 'pGS', [True], ici_pth)
        write_sol(model, 'pLS', [True], ici_pth)
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
            f.write(str(self.totLS))
        with open(os.path.join(ici_pth, "totPGS.csv"), "w", newline='') as f:
            f.write(str(self.totPGS))
        with open(os.path.join(ici_pth, "totCUT.csv"), "w", newline='') as f:
            f.write(str(self.totCUT))
        with open(os.path.join(ici_pth, "feastime.csv"), "w", newline='') as f:
            f.write(str(self.model._Tfeas))
            
        return self
        ##model.write("out.sol")

