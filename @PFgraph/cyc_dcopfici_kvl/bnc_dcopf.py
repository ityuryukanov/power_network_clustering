from gurobipy import Model, GRB, quicksum 
from utils import write_sol
import igraph  
import os 
import sys
sys.path.append(os.path.dirname(os.path.dirname( __file__ )))
from grb_read_graph import check_graph
import math

class Bnc_Dcopf(object):
    
    def __init__(self, args, RESLIM=None, OPTCR=0, DIGRAPH=False, MIPHEUR = None):        
		
        uarcs, arcs = check_graph(args['bdcpf'],args['nodes'],args['genss'],args['loads'])
        islid = args['islid']
        islsw = args['islsw']
        islid = [int(k) for k in islid]
        islid.sort(reverse=False)
        dstTR = args['dstTR']
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
        
        # Create igraph to obtain a numeric label for each node
        net = igraph.Graph(directed=True)    # for DFJ cut callbacks       
        grf = igraph.Graph(directed=False)   # for quick connectivity checks ONLY                 
        for n in self.nodes:
            net.add_vertex(n)
            grf.add_vertex(n)
        i2n = dict(); n2i = dict();
        for vx in net.vs: 
            i2n[vx.index]=vx['name']
            n2i[vx['name']]=vx.index
        for nd in grf.vs: 
            assert(n2i[nd['name']]==nd.index)
            assert(i2n[nd.index]==nd['name'])
        net.add_edges(list(arcs))           
        grf.add_edges(list(uarcs))      
        ccs = net.clusters().sizes()     
        if len(ccs)>1 or ccs[0]<len(self.nodes):
            raise RuntimeError('The input graph is disconnected')        
        i2a = dict(); a2i = dict();   #for net                
        for arc in arcs:
            es = net.es.find(_source=n2i[arc[0]],_target=n2i[arc[1]])
            es['weight'] = 1    # to speedup callbacks (graph is connected)
            assert(es.tuple==(n2i[arc[0]],n2i[arc[1]]))
            i2a[es.index] = (n2i[arc[0]],n2i[arc[1]])
            a2i[(n2i[arc[0]],n2i[arc[1]])] = es.index
        e2i = dict(); i2e = dict();  #for grf
        for uarc in uarcs:
            es = grf.es.find(_source=n2i[uarc[0]],_target=n2i[uarc[1]])
            assert(len(es)==1)   # grf is undirected so find() works in both directions
            es['weight'] = 1    # to speedup callbacks (graph is connected)
            e2i[(n2i[uarc[0]],n2i[uarc[1]])] = es.index
            i2e[es.index] = (n2i[uarc[0]],n2i[uarc[1]])
            
        # # FOR VIZ ONLY!
        # grf.es['label'] = [None  for w in grf.es['weight']]
        # igraphmod.igraph2graphviz(grf,'grfviz')            
        
        # Convert all data from labels to indices        
        self.nodes = [n2i[i] for i in self.nodes]    # convert all sets of nodes from labels to indices
        self.genss = [n2i[i] for i in self.genss]
        self.loads = [n2i[i] for i in self.loads]
        islsw = [n2i[i] for i in islsw] 
        uarcs = [(n2i[i],n2i[j]) for i,j in uarcs]   # convert uarcs and all branch data from labels to indices       
                
        brLup = dict(((n2i[i],n2i[j]),+abs(float(w))) for i,j,w in self.brlim)       # line limits are symmetric and positive    
        brLlo = dict(((n2i[i],n2i[j]),-abs(float(w))) for i,j,w in self.brlim)       # line limits are symmetric and positive 
        self.bdcpf = dict(((n2i[i],n2i[j]),float(w)) for i,j,w in self.bdcpf)        
        self.brlim = dict(((n2i[i],n2i[j]),abs(float(w))) for i,j,w in self.brlim)   # line limits are symmetric and positive         
        self.flow0 = dict(((n2i[i],n2i[j]),float(w)) for i,j,w in self.flow0)        
        self.pL0 = dict((n2i[n],float(w)) for (n,w) in self.pL0.items())   # convert all nodal data from labels to indices
        self.pG0 = dict((n2i[n],float(w)) for (n,w) in self.pG0.items())
        self.pGmin = dict((n2i[n],float(w)) for (n,w) in self.pGmin.items())   
        self.pGmax = dict((n2i[n],float(w)) for (n,w) in self.pGmax.items())
        self.pLmin = dict((n2i[n],float(w)) for (n,w) in self.pLmin.items())
        self.pLmax = dict((n2i[n],float(w)) for (n,w) in self.pLmax.items())
        self.costL = dict((n2i[n],float(w)) for (n,w) in self.costL.items())                  
        if dstTR is not None: dstTR = dict(((n2i[i],n2i[j]),float(w)) for i,j,w in dstTR)
         
        # Define Hi and Lo bounds for generator and load shedding
        # self.pL0[i] and self.pG0[i] as values for maximum load and generator 
        # shedding are included to make some study cases feasible. They can be 
        # excluded if respecting generator and load power limits is crucial. 
        lo_LS = dict((i, min([0.000000000,self.pL0[i]-self.pLmax[i],self.pL0[i]-self.pLmin[i]])) for i in self.loads)
        up_LS = dict((i, max([self.pL0[i],self.pL0[i]-self.pLmax[i],self.pL0[i]-self.pLmin[i]])) for i in self.loads)
        lo_GS = dict((i, min([0.000000000,self.pG0[i]-self.pGmax[i],self.pG0[i]-self.pGmin[i]])) for i in self.genss)
        up_GS = dict((i, max([self.pG0[i],self.pG0[i]-self.pGmax[i],self.pG0[i]-self.pGmin[i]])) for i in self.genss)
            
        # Introduce p as line power flows in each uarc        
        p = model.addVars(uarcs, name="p", vtype=GRB.CONTINUOUS, lb=brLlo, ub=brLup)   
        # Introduce pLS as load shedding at load nodes
        pLS = model.addVars(self.loads, name="pLS", vtype=GRB.CONTINUOUS, lb=lo_LS, ub=up_LS, obj=1.0)
        # Introduce pGS as generator shedding at generator nodes
        pGS = model.addVars(self.genss, name="pGS", vtype=GRB.CONTINUOUS, lb=lo_GS, ub=up_GS, obj=0.0)

        # Power balance
        for i in self.nodes:
            if i in self.genss and i in self.loads:
                model.addConstr(p.sum(i,'*')-p.sum('*',i)==self.pG0[i]-pGS[i]-self.pL0[i]+pLS[i], "PowrBal[%s]" % (i2n[i]))
            elif i in self.genss:
                model.addConstr(p.sum(i,'*')-p.sum('*',i)==self.pG0[i]-pGS[i], "PowrBal[%s]" % (i2n[i]))
            elif i in self.loads:
                model.addConstr(p.sum(i,'*')-p.sum('*',i)==0-self.pL0[i]+pLS[i], "PowrBal[%s]" % (i2n[i]))
            else:
                model.addConstr(p.sum(i,'*')-p.sum('*',i)==0, "PowrBal[%s]" % (i2n[i]))
                
        # Introduce ph as bus voltage angle difference over each branch at each bus
        ph = model.addVars(self.nodes, name="ph", vtype=GRB.CONTINUOUS, lb=-math.pi, ub=math.pi)
        for i in islsw:
            ph[i].lb = 0
            ph[i].ub = 0
        
        # Adjust "initial power flows weights"    
        for i,j in uarcs:
            if (i,j) not in self.flow0:
                self.flow0[(i,j)]=0  

        # Sort terminals in dstTR by the distance to their roots
        if dstTR is not None:
            rtd = dict()         # root sets
            for sr in dstTR.keys(): 
                s = sr[0]  # "sink"
                R = sr[1]  # "root"
                d = dstTR[(s,R)]
                if R in rtd:
                    rtd[R].append((s,d))
                else:
                    rtd[R] = [(s,d)]
            for R in rtd.keys():        
                rtd[R] = sorted(rtd[R], key=lambda x: x[1], reverse=False)
        else: rtd = None
        
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
        model._rtd = rtd
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
        model._MipHeur = MIPHEUR
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
        if model.isMIP == 0:
            print('Model is not a MIP')
            exit(0)
        try:
            if model._rtd is not None:
                model.optimize(callbackfcn)
            else:
                model.optimize()
        except GurobiError:
            print(GurobiError.message)
        
        if model.status==GRB.Status.INF_OR_UNBD or model.status==GRB.Status.INFEASIBLE or model.status==GRB.Status.UNBOUNDED: 
            model.computeIIS()
            model.write("test.ilp")   
            return None                     
            
        self.optimal    = model.Status == GRB.OPTIMAL
        self.runtime    = model.Runtime
        self.node_count = model.nodecount
        self.mip_gap    = model.mipgap
        self.objective  = model.ObjVal        
        self.totLS      = sum([model._pLS[l].x for l in self.loads])
        self.cstLS      = self.nu*sum([self.costL[l]*model._pLS[l].x for l in self.loads])     
        ##self.powLOD = sum([model._pLS[l].x-self.pL0[l] for l in self.loads])
        ##self.powGEN = sum([self.pG0[g]-model._pGS[i].x for g in self.genss])
        ##self.iniLOD = sum([-self.pL0[l] for l in self.loads])
        ##self.iniGEN = sum([ self.pG0[g] for g in self.genss])
            
        # Extract the solution     
        ici_pth = os.path.abspath(os.path.join(os.path.dirname( __file__ ), 'ici'))              
        write_sol(model, 'p', [True, True], ici_pth)     
        write_sol(model, 'pGS', [True], ici_pth)  
        write_sol(model, 'pLS', [True], ici_pth)  
        with open(os.path.join(ici_pth, "runtime.csv"), "w", newline='') as f:
            f.write(str(self.runtime))        
        return self
        ##model.write("out.sol")
