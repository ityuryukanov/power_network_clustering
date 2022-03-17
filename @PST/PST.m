classdef PST
  % 
  %PST A static class containing functions based on PST (Power System 
  %  Toolbox by J.H. Chow). These are the PST functions from PSTOLD
  %  related to power system computation tasks (Ybus, loadflow etc).
  %  
  % This class was created to host several MODIFIED PST functions with 
  % IMPROVED functionality (YYbus, LLoadflow) (mainly better  handling
  % of phase shifters). Then it was augmented with unchanged copies of  
  % other used PST functions to keep things together.
  %
  % The class copies/reproduces some key functionalities of PSTOLD (the
  % first version of PST written by J.H. Chow) and is focused on generator
  % coherency and aggregation. The later augmentations by G. Rogers (pstv3,
  % pstdat, statespace) are not in this class. 
  % 
  % However, savepst() and loadpst() and some others can be used with 
  % pstv3 as well. 
  % 
   
  properties
  
  end
  
  methods (Static)
	pst_var
	loadpst(pststate)
  pststate = savepst()	
  [Y,nSW,nPV,nPQ,SB] = YYbus(bus,line)
	[Y11,Y12,Y21,Y22,rec_V1,rec_V2,bus_order] = RED_YYbus(bus_sol,line)
  [mac_ord, tim_sim, mac_ang, mac_spd] = SIM_PST(pst_str, mac_ref, bus_flt, bus_adj, t_flt, t_clr, t_end)           
	[n_mac] = eqgen(mac,mac_list,basemva,bus_num,mac_num)
	[MK, K, M] = svm_em(pst)
	[f] = mac_em(i,k,bus,flag)
  [T] = sc_xform(area,nmach_a,mac_con)
	[V_s, lambda, MK, K, M] = pstVslow(pst, nc)	
	[bus_sol,line_flow] = LLoadflow(bus,line,tol,iter_max,vmin,vmax,acc,display,flag)
	[S1,S2] = LLine_PQ(V1,V2,R,X,B,tap,phi)
  red_mod = zhukovpp(pst,grp,nm,trmbusagg,modefit)	
	[n_bus,n_line,nmac_con] = mi_agg(bus,line,area,nmach_a,basemva)	       
	[Jac11,Jac12,Jac21,Jac22] = form_jac(V,ang,Y,ang_red,volt_red)
	[delP,delQ,P,Q,conv_flag] = calc(nbus,bus_type,V,ang,Y,Pg,Qg,Pl,Ql,tol)                 
	[area,nmach_a,L, mach_ord] = L_group(V_s)
	[area,nmach_a,L,mach_ord,LdLgN] = LRgroup(V_s)
	[area,nmach_a,arear,nmach_ar] = ex_group(grouping)
	[T, C] = coh_loose(Vs, tol)
  [mac_con,mac0ref] = agg_mac(mac_con,mac0ref,basmva)
	[grouping,c_map] = coh_map(V_s,tol)
	[area,nmach_a,areal,nmach_al,list,S] = weaklink(MK,tol)                                          
	pst = coh_PST(pst, ns, meth, tol)
  end
  
end

