classdef GraphUtils
  %GraphUtils A static class containing some useful functions on graphs.
  %  These methods are too general-purpose to restrict them to work on an
  %  instance of a class (as an object should be always first created)
  
  properties
    
  end
  
  methods (Static)
    [S, C] = alecconncomp(G) 
    allSP = allShortestPaths(adj, mode) 
    [graphD, dm, nnInd, dst_srtd] = buildDistanceGraph(Y, adj, nm) 
    L = createLaplacian(adj, nm, vw) 
    ok = cmpIncAdj(inc, adj) 
    T = hclust(pts, dm, link, nc, varargin) 
    nSP = nShortestPaths(adj, m_last, mode) 
    adj_out = simadj2distadj(adj, mode) 
    inc = adj2inc( adj ) 
    vxlbl = indvec2vxlbl(busind) 
    graphviz(path, gv_command) 
    [adj_sep, C, ixs, eXp, pcut] = cc_info(adj, inc, vw, cutind) 
    [eXp, ixs, pcut] = lbl_info(adj, vw, T) 
    [cutind] = cutset(inc_matr, busind) 
    [idx_out, val] = edges2adj( adj, inc, linidx, NDEBUG) 
    [adj_out, vw_out, merge_map, ixsc] = merge_nodes(adj_in, ixs, vw_in) 
    [Y, v] = spClust(L, lo_up, nc) 
    [T] = sp_postpr(Y, postpr, pts_repr, nc, adj) 
    [Vr_out, Qt_out, t_run] = nc_ev_sym(V, rm, dm, batch_siz, k_min, debug) 
    [Z_out, J_out, t_run] = rotVslow(V, mode, dm, batch_siz, k_min) 
    pfgraph2hmetis(adj, inc, vw, coh) 
    [T] = min_cut_rb(adj, vs, vw, debug) 
    ixs = aggregate_mls( bran_ml, coh_ml )
    [T, i_best, contig] = best_part(cntg, all_bus, weXp, outl_gen, dif)
	mcb = mcb_generator(caseid, GRBDIR)
  end
  
end

