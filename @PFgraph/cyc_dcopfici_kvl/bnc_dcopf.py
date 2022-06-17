from gurobipy import Model, GRB, tupledict 
from utils import write_sol, opf_lbl2idx_grf
import igraph  
import os 
import sys
sys.path.append(os.path.dirname(os.path.dirname( __file__ )))
from grb_read_graph import check_graph
import math

class Bnc_Dcopf(object):
    """
    Define a dc network flow model (classic dc opf with  voltage angles).
    """
    def __init__(self, args, RESLIM=None, OPTCR=0):
        uarcs, arcs = check_graph(args['bdcpf'],args['nodes'],args['genss'],args['loads'])
        islid = [int(k) for k in args['islid']]
        islid.sort(reverse=False)				
        islsw = args['islsw']
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
        # Introduce pGS as generator shedding at generator nodes.
        pGS = model.addVars(self.genss, name="pGS", vtype=GRB.CONTINUOUS, lb=lo_GS, ub=up_GS)

        # Power balance
        for i in self.nodes:
            if i in self.genss and i in self.loads:
                model.addConstr(p.sum(i,'*')-p.sum('*',i)==self.pG0[i]-pGS[i]-self.pL0[i]+pLS[i], "PowrBal[%s]" % (i))
            elif i in self.genss:
                model.addConstr(p.sum(i,'*')-p.sum('*',i)==self.pG0[i]-pGS[i], "PowrBal[%s]" % (i))
            elif i in self.loads:
                model.addConstr(p.sum(i,'*')-p.sum('*',i)==0-self.pL0[i]+pLS[i], "PowrBal[%s]" % (i))
            else:
                model.addConstr(p.sum(i,'*')-p.sum('*',i)==0, "PowrBal[%s]" % (i))
                
        # Introduce ph as bus voltage angle difference over each branch at each bus
        islsw = [n2i[i] for i in islsw]
        ph = model.addVars(self.nodes, name="ph", vtype=GRB.CONTINUOUS, lb=-math.pi, ub=math.pi)
        for i in islsw:
            ph[i].lb = 0
            ph[i].ub = 0
        
        # Adjust "initial power flows weights"    
        for i,j in uarcs:
            if (i,j) not in self.flow0:
                self.flow0[(i,j)]=0  
        
        # Store the variables
        model._p = p
        model._pLS = pLS
        model._pGS = pGS
        model._ph = ph
        # BnC information
        model._relobj = None
        model._imprv = 0
        model._graph = net
        model._ntwrk = grf
        model._i2n = i2n
        model._n2i = n2i
        model._i2a = i2a
        model._a2i = a2i
        model._e2i = e2i
        model._i2e = i2e        
        # runtime information
        model._troot = 0
        model._ttree = 0
        model.Params.PreCrush = 1
        model.Params.LazyConstraints = 1
        model._islsw = islsw
        model._ztruu = False   # default values of callback flags
        model._kvlfl = False   # default values of callback flags
        # Finalize self
        model.NumStart = 4
        self.arcs  = arcs        
        self.uarcs = uarcs 
        self.model = model            
        self.islid = islid


    def solve(self, callbackfcn):                    
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
        self.cstLS = self.nu*sum([self.costL[l]*self.model._pLS[l].x  for l in self.loads])
            
        # Extract the solution     
        ici_pth = os.path.abspath(os.path.join(os.path.dirname( __file__ ), 'ici'))              
        write_sol(model, 'p', [True, True], ici_pth)     
        write_sol(model, 'pGS', [True], ici_pth)  
        write_sol(model, 'pLS', [True], ici_pth)  
        with open(os.path.join(ici_pth, "runtime.csv"), "w", newline='') as f:
            f.write(str(self.runtime))        
        return self
        ##model.write("out.sol")
