import os 
import sys
import pickle
import igraph
import warnings
import networkx as nx
sys.path.append(os.path.join(os.path.dirname(os.path.dirname( __file__ )),'@PFgraph'))
import grb_read_graph
from grb_read_graph import *

def findPaths(G,vx,dpt,strt,visited):
    if dpt==0 or len(set(G.neighbors(vx))-visited)==0 or (strt in set(G.neighbors(vx)) and len(visited)>2):
        return [[vx]]
    visnxt = visited.copy()
    visnxt.add(vx)    
    paths = []
    for neighbor in set(G.neighbors(vx))-visnxt:
        for path in findPaths(G,neighbor,dpt-1,strt,visnxt):
            if vx not in path:
                paths.append([vx]+path)
    return paths

#visited.add(vx)
#visited.add(strt)

def find_cycles(G,vx,dpt):
    paths = findPaths(G,vx,dpt,vx,set())
    bckvx = set(G.neighbors(vx))    
    paths = [list(path) for path in paths if path[-1] in bckvx and len(set(path))==len(path) and len(path)>2]
    return paths

def mcbs(ici_pth):
	"""
	Normally called from mcb_generator.m, which also generates the input data files. 
	"""

    # Case name
    f = open(os.path.join(ici_pth, 'caseid.txt'), 'r')
    caseid = f.read()
    f.close()
    
    # Graph parameters
    bdcpf, nodes = read_file(os.path.join(ici_pth, 'graph.csv'))
    flow0, _ = read_file(os.path.join(ici_pth, 'flow0.csv'))
    brlim, _ = read_file(os.path.join(ici_pth, 'brlim.csv'))
    buses = read_file(os.path.join(ici_pth, 'nodes.csv'))
    
    # Unidirect input graphs
    bdcpf = unidirect(bdcpf)
    flow0 = unidirect(flow0)
    brlim = unidirect(brlim)
    
    # Create optimization model
    assert set(nodes)==buses, "Set of nodes implied from admittance graph should be equal to set of nodes read from nodes.csv"
    args = {'nodes':nodes,
    		  'bdcpf':bdcpf,
    		  'flow0':flow0,
    		  'brlim':brlim}
    uarcs, arcs = check_graph(args['bdcpf'],args['nodes'])
    
    # Create igraph reference model to be copied into networkX
    net = igraph.Graph(directed=False);                    
    for n in args['nodes']:
        net.add_vertex(n)
    net.add_edges(list(uarcs))
    ccs = net.clusters().sizes()     
    if len(ccs)>1 or ccs[0]<len(args['nodes']):
        raise RuntimeError('The input graph is disconnected')   
    i2n = dict()
    n2i = dict()
    for i in range(len(net.vs)):
        i2n[i]=net.vs[i]['name'] 
        n2i[net.vs[i]['name']]=i        
    uarcs = [(n2i[i],n2i[j]) for i,j in uarcs]   # convert uarcs and all branch data from labels to indices    
	
    # Create networkx based on igraph to obtain minimal cycle basis
    grf = nx.Graph()
    for i in range(len(net.es)):
        grf.add_edge(net.es[i].source, net.es[i].target, weight = 1.0)
        net.es[i]['weight'] = 1
        
    # Prepare savefiles
    savedir = os.path.join(os.path.dirname(os.path.dirname( __file__ )), "@PFgraph", "cycledata");
    mcbfile = os.path.join(savedir, 'mcb_'+caseid+'.pickle')	
    cycfile = os.path.join(savedir, 'cyc_'+caseid+'.pickle')
    mstfile = os.path.join(savedir, 'mst_'+caseid+'.pickle')
    
    # # Save the minimal cycle basis		
    # print("Start computing MCB")
    # mcb_lst = nx.minimum_cycle_basis(grf)
    # print("Done computing MCB")
    # mcb = set()
    # for cyc0 in mcb_lst:
    #     n_cyc = len(cyc0)
    #     loop_edg = set()
    #     for i in cyc0:
    #         lst2 = cyc0.copy()
    #         lst1 = [i]*(n_cyc-1)
    #         lst2.remove(i)
    #         cnd_edg = zip(lst1,lst2) 
    #         for edg in cnd_edg:
    #             if edg in uarcs:
    #                 loop_edg.add(edg)
    #         cnd_edg = zip(lst2,lst1) 
    #         for edg in cnd_edg:
    #             if edg in uarcs:
    #                 loop_edg.add(edg)   
    #     mcb.add(tuple(sorted(loop_edg)))
    # mcb = list(mcb)        
    # with open(mcbfile, 'wb') as f:
    #     pickle.dump(mcb, f)
    # print("Done writing MCB")        
    
    # Save other cycles shorter than 9
    print("Start computing other cycles")
    cyc = set();
    for v in grf.nodes():
        vx = find_cycles(grf,v,7)
        for loop in vx:
            loop_edg = set()
            cnd_edg = zip(loop,loop[1:]+[loop[0]])        
            for edg in cnd_edg:
                if (edg[0],edg[1]) in uarcs:
                    loop_edg.add((edg[0],edg[1]))
                elif (edg[1],edg[0]) in uarcs:     
                    loop_edg.add((edg[1],edg[0]))
                else:
                    raise RuntimeError('No such edge in the graph')
            cyc.add(tuple(sorted(loop_edg)))
    cyc = list(cyc)
    with open(cycfile, 'wb') as f:
        pickle.dump(cyc, f)
    print("Done writing other cycles") 		
		
	# #Save the minimal spanning tree cycle basis
    # mst  = set()
    # esIDs = net.spanning_tree(return_tree=False, weights=net.es['weight'])
    # es_tr = [(net.es[esID].source,net.es[esID].target) for esID in esIDs]
    # esALL = [a.index for a in net.es]    
    # esCYC = set(esALL) - set(esIDs)    # closing edges of each mst cycle
    # upwgt = len(esALL)+1  
    # for esID in esCYC:
        # net.es[esID]['weight'] = upwgt # one non-mst edge per cycle
    # for esID in esCYC:
        # f = net.es[esID].source
        # t = net.es[esID].target            
        # with warnings.catch_warnings():
            # warnings.simplefilter("ignore")        
            # escyc = net.get_shortest_paths(f,to=t,mode="all",weights=net.es['weight'],output="epath")[0]              
            # cycle = [(net.es[esix].source, net.es[esix].target) for esix in escyc] + [(net.es[esID].source,net.es[esID].target)]
        # try: 
            # if cycle[-1][0]!=f or cycle[-1][1]!=t: 
                # continue               # cycle has not been properly closed 
        # except: 
            # continue
        # mst.add(tuple(sorted(cycle)))
    # mst = list(mst)
    # with open(mstfile, 'wb') as f:
        # pickle.dump(mst, f)        
    return cyc	

##mcbs('E:\\proj\\Data\\program_codes\\MATLAB\\Splitting_and_Clustering\\,code4papers\\Benders_2018\\ici')







