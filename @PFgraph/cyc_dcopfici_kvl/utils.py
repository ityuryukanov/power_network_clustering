import os
import re
import warnings  
import igraph
from timeit import default_timer as timer
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

def opf_lbl2idx_grf(obj, arcs, uarcs):
    """
    Implements a common segment in cycopf.py and bnc_dcopf.py that is 
    meant to create the auxiliary directed and undirected graphs (using
    Python igraph), to create mappings between these graphs and decision
    variables in the model, and transform text label based variable 
    indices to integer based variable indices.
    """
    # Create igraph to obtain a numeric label for each node
    net = igraph.Graph(directed=True)    # for DFJ cut callbacks       
    grf = igraph.Graph(directed=False)   # for quick connectivity checks ONLY                 
    for n in obj.nodes:
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
    if len(ccs)>1 or ccs[0]<len(obj.nodes):
        raise RuntimeError('The input graph is disconnected')        
    i2a = dict(); a2i = dict();
    for arc in arcs:
        es = net.es.find(_source=n2i[arc[0]],_target=n2i[arc[1]])
        es['weight'] = 1     # to speedup callbacks (graph is connected)
        assert(es.tuple==(n2i[arc[0]],n2i[arc[1]]))
        i2a[es.index] = (n2i[arc[0]],n2i[arc[1]])
        a2i[(n2i[arc[0]],n2i[arc[1]])] = es.index        
    uarcs = [(i2n[es.source],i2n[es.target]) for es in grf.es]  # possibly redefine uarcs to be conform with grf for convenience later
    e2i = dict(); i2e = dict();  #for grf
    for uarc in uarcs:
        es = grf.es.find(_source=n2i[uarc[0]],_target=n2i[uarc[1]])
        assert(len(es)==1)   # grf is undirected so find() works in both directions
        es['weight'] = 1     # to speedup callbacks (graph is connected)
        e2i[(n2i[uarc[0]],n2i[uarc[1]])] = es.index
        i2e[es.index] = (n2i[uarc[0]],n2i[uarc[1]])
    # # FOR VIZ ONLY!
    # grf.es['label'] = [None  for w in grf.es['weight']]
    # igraphmod.igraph2graphviz(grf,'grfviz')
    
    # Convert all data from labels to indices        
    obj.nodes = [n2i[i] for i in obj.nodes]      # convert all sets of nodes from labels to indices
    obj.genss = [n2i[i] for i in obj.genss]
    obj.loads = [n2i[i] for i in obj.loads]  		
    uarcs = [(n2i[i],n2i[j]) for i,j in uarcs]   # convert uarcs and all branch data from labels to indices		
    brLup = dict()
    brLlo = dict()
    for i,j,w in obj.brlim:
        if (n2i[i],n2i[j]) in uarcs:
            brLup[(n2i[i],n2i[j])] = +abs(float(w))   # line limits are symmetric and positive
            brLlo[(n2i[i],n2i[j])] = -abs(float(w)) 
        else:
            brLup[(n2i[j],n2i[i])] = +abs(float(w))   # line limits are symmetric and positive
            brLlo[(n2i[j],n2i[i])] = -abs(float(w))
    bdcpf = dict()
    for i,j,w in obj.bdcpf:
        if (n2i[i],n2i[j]) in uarcs:
            bdcpf[(n2i[i],n2i[j])] = float(w)
        elif (n2i[j],n2i[i]) in uarcs:
            bdcpf[(n2i[j],n2i[i])] = float(w)
        else:
            raise RuntimeError('Edge does not exist!')                 
    obj.bdcpf = bdcpf;
    brlim = dict()            
    for i,j,w in obj.brlim:
        if (n2i[i],n2i[j]) in uarcs:
            brlim[(n2i[i],n2i[j])] = abs(float(w))
        elif (n2i[j],n2i[i]) in uarcs:
            brlim[(n2i[j],n2i[i])] = abs(float(w))
        else:
            raise RuntimeError('Edge does not exist!')
    obj.brlim = brlim;
    flow0 = dict() 
    for i,j,w in obj.flow0:
        if (n2i[i],n2i[j]) in uarcs:
            flow0[(n2i[i],n2i[j])] = float(w)
        elif (n2i[j],n2i[i]) in uarcs:
            flow0[(n2i[j],n2i[i])] = float(w)
        else:
            raise RuntimeError('Edge does not exist!')                
    obj.flow0 = flow0;
    obj.pL0   = dict((n2i[n],float(w)) for (n,w) in obj.pL0.items())       # convert all nodal data from labels to indices
    obj.pG0   = dict((n2i[n],float(w)) for (n,w) in obj.pG0.items())
    obj.pGmin = dict((n2i[n],float(w)) for (n,w) in obj.pGmin.items())   
    obj.pGmax = dict((n2i[n],float(w)) for (n,w) in obj.pGmax.items())
    obj.pLmin = dict((n2i[n],float(w)) for (n,w) in obj.pLmin.items())
    obj.pLmax = dict((n2i[n],float(w)) for (n,w) in obj.pLmax.items())
    obj.costL = dict((n2i[n],float(w)) for (n,w) in obj.costL.items())
        
    # Define Hi and Lo bounds for generator shedding and load values obj.pL0[i] and
    # obj.pG0[i] as values for maximum load and maximum generator shedding are included 
    # to make some study cases feasible. They can be removed if needed.
    assert all([obj.pL0[i]>=0 for i in obj.loads])   # no negative loads, which are, in fact, generators
    lo_GS = dict((i, min([0.000000000,obj.pG0[i]-obj.pGmax[i],obj.pG0[i]-obj.pGmin[i]])) for i in obj.genss)
    up_GS = dict((i, max([obj.pG0[i],obj.pG0[i]-obj.pGmax[i],obj.pG0[i]-obj.pGmin[i]])) for i in obj.genss)
    lo_LS = dict((i, min([0.000000000,obj.pL0[i]-obj.pLmax[i],obj.pL0[i]-obj.pLmin[i]])) for i in obj.loads)
    up_LS = dict((i, max([obj.pL0[i],obj.pL0[i]-obj.pLmax[i],obj.pL0[i]-obj.pLmin[i]])) for i in obj.loads)    
    return (obj,net,grf,i2n,n2i,i2a,a2i,e2i,i2e,arcs,uarcs,lo_LS,up_LS,lo_GS,up_GS,brLlo,brLup)


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

            