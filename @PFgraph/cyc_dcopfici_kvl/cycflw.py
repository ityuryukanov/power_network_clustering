from gurobipy import GRB, quicksum, tupledict
from ximblnce import Ximblnce
from utils import cycle_ordering, cycle_signing, order_checking, xtrct_filam

class Cycflw(Ximblnce):
    """
    Here the KVL cycle flow constraints and the spanning tree constraints with 
    directed arc variables z[i,j] are introduced.
    """ 
    def __init__(self, args, RESLIM=None, OPTCR=0, MIPSTART=False):
        
        mipstart = MIPSTART
        Ximblnce.__init__(self, args, RESLIM, OPTCR, MIPSTART=mipstart)
        
        # Create the graph partitioning model (initial solve)
        modOPF = self.model
        modOPF.update()

        # Extra parameters/variables
        z_ini = args['z_ini']  
        z_ini = dict(((modOPF._n2i[i],modOPF._n2i[j]), int(w)) for i,j,w in z_ini)
        y_ini = args['y_ini']  
        y_ini = dict(((modOPF._n2i[i],modOPF._n2i[j]), int(w)) for i,j,w in y_ini)
        self.arcs = [(modOPF._n2i[i],modOPF._n2i[j]) for i,j in self.arcs]
        modOPF._ztruu = True
        modOPF._z = modOPF.addVars(list(self.arcs), name="z", vtype=GRB.BINARY, lb=0, ub=1)
        modOPF._y = modOPF.addVars(list(self.uarcs), name="y", vtype=GRB.BINARY, lb=0, ub=1)
        arccyc = {(i,j) if (i,j) in self.uarcs else (j,i) for cyc in self.mstcb for i,j in cyc}
        arcarb = set(self.uarcs) - arccyc 
        assert(len(arccyc)+len(arcarb)==len(self.uarcs))
        # Set initial values of z[i,j] to z_ini
        if MIPSTART:
            modOPF.NumStart = 2; modOPF.Params.StartNumber = 0;
            for i,j in z_ini:
                modOPF._z[i,j].start = z_ini[(i,j)]
            for i,j in y_ini:
                modOPF._y[i,j].start = y_ini[(i,j)]
                
        # Fix terminal constraints
        alltrm = modOPF._altrm
        trm2k  = modOPF._trm2k
        k2trm  = modOPF._k2trm
        for i in trm2k.keys():
            modOPF._x[i,trm2k[i]].lb = 1
            modOPF._x[i,trm2k[i]].ub = 1
        # Spanning tree cardinality constraint
        modOPF.addConstr(quicksum([modOPF._z[i,j] for i,j in self.arcs])==len(self.nodes)-len(self.islid), "sptree")
        # Root nodes have no incoming arcs
        for r in modOPF._islsw:
            Nr = set(modOPF._graph.neighbors(r))            
            for j in Nr:
                modOPF._z[j,r].lb = 0   #rootIN[%s]
                modOPF._z[j,r].ub = 0                   
        # Root nodes of MULTITERMINAL groups have at least one outgoing arc
        assert(modOPF._rtd is not None)
        root_mlt = [r for r in modOPF._rtd if len(modOPF._rtd[r])>0]
        root_sgl = [r for r in modOPF._rtd if len(modOPF._rtd[r])==0]
        for r in root_mlt:
            Nr = set(modOPF._graph.neighbors(r))
            modOPF.addConstr(quicksum([modOPF._z[r,j] for j in Nr])>=1, "rootOUT[%s]" % (r))
        # Non-root nodes have strictly one incoming arc
        nonsw = set(self.nodes) - set(modOPF._islsw)
        for v in nonsw:
            Nv = set(modOPF._graph.neighbors(v))
            modOPF.addConstr(quicksum([modOPF._z[j,v] for j in Nv])==1, "nodeIN[%s]" % (v))
        # Fix z_ij for network filaments        
        arcsfil, fils_in, filsout, nods = xtrct_filam(modOPF._graph,root_mlt,root_sgl)
        arcsfix = [arc  for fil in arcsfil for arc in fil]
        for i,j in arcsfix:
            assert(((i,j) in arcarb) or ((j,i) in arcarb))
            modOPF._z[i,j].lb = 1
            modOPF._z[i,j].ub = 1
            modOPF._z[j,i].lb = 0
            modOPF._z[j,i].ub = 0
                
        # Relationship between y and z
        for i,j in arccyc:
            modOPF.addConstr(modOPF._z[i,j]+modOPF._z[j,i] <= 1-modOPF._y[i,j], "validIJ[%s,%s]" % (i,j))
        for i,j in arcarb:
            modOPF.addConstr(modOPF._z[i,j]+modOPF._z[j,i] == 1-modOPF._y[i,j], "validIJ[%s,%s]" % (i,j))
		
        # Cycle breaking constraints for available cycles (KEEP IT, AS SHORT CYCLES 
        # HELP THE HEURISTICS AND THE SOLUTION ITSELF)
        len_lim = self.cycmax 
        for icyc,cyc in enumerate(self.cyc):
            if len(cyc)<=len_lim:
                cyc_ord = cycle_ordering(cyc)
                arc_ord = order_checking(cyc_ord)
                cycle = {v for arc in arc_ord for v in arc}
                modOPF.addConstr(quicksum(modOPF._z[arc[0],arc[1]] for arc in arc_ord)<=len(cycle)-1, "CBC_CYC_FWD[%s]" % (icyc))
                modOPF.addConstr(quicksum(modOPF._z[arc[1],arc[0]] for arc in arc_ord)<=len(cycle)-1, "CBC_CYC_BCK[%s]" % (icyc))

        # Definition of SOS1 variables
        for i,j in self.uarcs:
            modOPF.addSOS(GRB.SOS_TYPE1, [modOPF._z[i,j],modOPF._z[j,i]])
        # Switching constraints
        for i,j in self.uarcs:
            for k in self.islid:
                modOPF.addConstr(modOPF._x[i,k]-modOPF._x[j,k] <= modOPF._y[i,j], "lblXij1[%s,%s,%s]" % (i,j,k))
                modOPF.addConstr(modOPF._x[j,k]-modOPF._x[i,k] <= modOPF._y[i,j], "lblXij2[%s,%s,%s]" % (j,i,k))
                modOPF.addConstr(modOPF._x[i,k]+modOPF._x[j,k]-1+modOPF._y[i,j] <= 1, "SWITCHG[%s,%s,%s]" % (i,j,k))   #if removed, line switching inside islands is OK
    
        # Switching & power flow limitation
        for i,j in self.uarcs:           
            if (i,j) in self.brlim and self.brlim[(i,j)]>0:
                modOPF.addConstr(modOPF._p[i,j]<= self.brlim[(i,j)]*(1-modOPF._y[i,j]), "LinLim1[%s,%s]" % (i,j))
                modOPF.addConstr(modOPF._p[i,j]>=-self.brlim[(i,j)]*(1-modOPF._y[i,j]), "LinLim2[%s,%s]" % (i,j))

        # Assign sign to each cycle edge
        mcb_lin = list() 
        srt     = {tuple(sorted(cyc)) for cyc in self.srt}    #100% avoid repeats!
        mstmcb  = {tuple(sorted(cyc)) for cyc in self.mstcb}
        cyc_all = mstmcb.union(srt)    # always sync with kvl_all
        for cyc0 in cyc_all:
            cyc_sgn = cycle_signing(cyc0)
            mcb_lin.append(cyc_sgn)
        
        # Switching constraints for cycles  (cycle length is encoded into the 
        # name of each constraint)           
        modOPF._B = self.bdcpf   # for cuts insertion in callbacks
        modOPF._M = self.brlim   # for cuts insertion in callbacks
        for cyc in mcb_lin:
            th_ij = list()
            for i,j,pm in cyc:
                th_ij.append(abs(self.brlim[(i,j)])/abs(self.bdcpf[(i,j)]))   # abs()/abs() is the "worst case"
            th_ij.sort(reverse=True)
            th_ij.pop(); th_ij.pop();
            M_cyc = sum(th_ij)
            pthlb = tuple(sorted([(i,j) for i,j,pm in cyc]))            
            lname = "["+','.join(['('+str(i)+','+str(j)+')' for i,j in pthlb])+"]"
            if len(lname)>240:
                lname = lname[0:240]
            modOPF.addConstr(quicksum(modOPF._p[i,j]/self.bdcpf[(i,j)]*pm for i,j,pm in cyc)<=+0.5*M_cyc*quicksum(modOPF._y[i,j] for i,j,pm in cyc), "KVL_Y_UP[%s||%s]" % (len(pthlb),lname))
            modOPF.addConstr(quicksum(modOPF._p[i,j]/self.bdcpf[(i,j)]*pm for i,j,pm in cyc)>=-0.5*M_cyc*quicksum(modOPF._y[i,j] for i,j,pm in cyc), "KVL_Y_LO[%s||%s]" % (len(pthlb),lname))
            continue

        # Cycle constraints to ensure that 0 or more than 1 line is switched  
        # for each cycle
        len_lim = 4
        n_yyy = 0
        cyc_all = cyc_all.union(self.cyc)  # always sync with kvl_all
        for cyc0 in cyc_all:
            if len(cyc0)>len_lim:
                continue
            cyc = set(cyc0)
            cycid = set([i for i,j in cyc] + [j for i,j in cyc])
            cycid = tuple(sorted(cycid))
            lname = "["+','.join([str(i) for i in cycid])+"]"
            for i,j in cyc0:
                cyc_rst = cyc-{(i,j)}
                modOPF.addConstr(modOPF._y[i,j] - quicksum(modOPF._y[p,q] for p,q in cyc_rst)<=0, "cycsw[%s,(%s,%s)]" % (lname,i,j))
                n_yyy+=1
                continue

        # Switching constraints with cycle variables l never work because they 
        # are trivial! They only cause a large slowdown.
        alpha = 0.00
        betta = 1.00
        gamma = 0.01
        muuuu = self.mu
        modOPF.setObjective(
        + alpha*quicksum(modOPF._IMBLOD[k] for k in self.islid)\
        + betta*quicksum(modOPF._pLS[l] for l in self.loads)    \
        + gamma*quicksum(modOPF._pGS[g] for g in self.genss)   \
        + muuuu*quicksum(self.flow0[(i,j)]*modOPF._y[i,j] for i,j in self.uarcs), GRB.MINIMIZE)
            
        # mod_gp based on modOPF for solution restarts
        modOPF.update()
        modOPF._ztruu = True     # callback flags
        modOPF._kvlfl = True     # callback flags
        modOPF._kvlall = cyc_all
        modOPF._cbkcyc = 0
        modOPF._Tfeas = -1       # time to feasibility
        ##modOPF.write("test.lp")        
        self.model = modOPF

    def solve(self, callbackfcn=None):
        self = Ximblnce.solve(self, callbackfcn)
        return self

