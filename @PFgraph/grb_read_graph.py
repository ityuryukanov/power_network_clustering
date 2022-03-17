def read_file(fname):
    nodes = set()
    edges = set()
    nodct = dict()
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if not line: continue
            xx = line.strip().split(',')
            nr = len(xx)
            if nr==1:
                edges.add(line)
            elif nr==2:
                nodct[xx[0]] = float(xx[1])
            elif nr>2:
                edges.add((xx[0],xx[1],xx[2]))
                nodes.add(xx[0])
                nodes.add(xx[1])
    if len(nodes)==0 and len(edges)==0 and len(nodct)==0:
        return None, None  # the file was empty
    if len(nodct)>0:
        return nodct
    # try to sort nodes to clearly map node labels with indices in igraph
    if len(nodes)>0:
        nodes = list(nodes)
        if all([isinstance(n,str) for n in nodes]) and all([n.isdigit() for n in nodes]):
            nodes.sort(key=lambda z:int(z), reverse=False)
        elif all([isinstance(n,(int,float,complex)) for n in nodes]):
            nodes.sort(reverse=False)
        return edges, nodes
    else:
        return edges

def unidirect(edgeset):
    """
    Creates unidirected edges ensuring f>t for every (weighted) edge.
    """
    temp = edgeset.copy()
    for u, v, w in edgeset:
        temp.add((v, u, w))
    edges = set()
    # Really try to ensure u>v if u and v are numeric strings or numbers
    for u, v, w in temp:
        if isinstance(u,str) and u.isdigit() and isinstance(v,str) and v.isdigit():
            if int(u)>int(v): 
                edges.add((u, v, w))
        elif u>v:
             edges.add((u, v, w))
    return edges
		
def check_graph(edges, nodes, genss=set(), loads=set()):
	vx  = set([i for (i, _, _) in edges])
	vx |= set([i for (_, i, _) in edges])
	assert(vx == set(nodes))        
	assert(genss <= set(nodes))
	assert(loads <= set(nodes))
	uarcs = set([(u,v)  for u,v,w in edges])
	arcs = uarcs.copy()  #define bi-directed arcs
	for i,j in uarcs:
		arcs.add((j,i))
	return uarcs, arcs	

