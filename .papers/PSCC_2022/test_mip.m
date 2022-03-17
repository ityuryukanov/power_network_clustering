function res = test_mip(PHI_OR_KVL, nc, caseid)
clc;
solver = 'gurobi';
caseid0 = 'case2869pegase';
nc0 = 7;
if nargin<1, nc=nc0; PHI_OR_KVL='KVL'; caseid=caseid0; end
if nargin<2, nc=nc0; caseid=caseid0; end
if nargin<3, caseid=caseid0; end

% Amend MATLAB path
old_warning = warning;
cleanup2 = onCleanup(@() warning(old_warning));
warning on;
old_matlabpath = path;
cleanup1 = onCleanup(@() path(old_matlabpath));
addpath( fullfile(fileparts(fileparts(fileparts(mfilename('fullpath'))))),...
  fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'third_party'),...
  fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'third_party', 'matlab_bgl') ); 
warning off;
% Initialise the studycase from MATPOWER
if strcmp(caseid,'datanp48')
  eval(caseid);
  % Prepare the data of PST to be consistent with PFgraph
  [bus(:,1),I] = sort(bus(:,1),'ascend');  % sort with respect to busIDs to ensure that...
  bus(:,2:end) = bus(I,2:end);
  [mac_con(:,2),I] = sort(mac_con(:,2),'ascend');  % busIDs in bus[] and gen[] are in the same order
  ncol_g = size(mac_con,2); mac_con(:,[1,3:ncol_g]) = mac_con(I,[1,3:ncol_g]);
  % Aggregate machines sitting on the same bus (the initial approach)
  samebus = tabulate(mac_con(:,2));
  samebus = samebus(samebus(:,2)>1,1);
  for i = 1:1:numel(samebus)
    b = samebus(i);
    mac_lst = find(mac_con(:,2)==b);
    agg_mac = PST.eqgen(mac_con,mac_lst,100,b,mac_con(mac_lst(1),1));
    agg_mac(1,22:23) = 1;
    mac_con(mac_lst,:) = [];
    mac_con = [mac_con;agg_mac];
  end
  [~,ord] = sort(mac_con(:,2),'ascend');  %a new meaninful machine order
  mac_con = mac_con(ord,:);
  mac_con(:,1) = (1:1:size(mac_con,1))';
  tol = 2e-13; iter_max = 30; vmin = 0.5;  vmax = 1.5; acc = 1.0;
  [DUMMY, bus, line_flw] = evalc('PST.LLoadflow(bus,lin,tol,iter_max,vmin,vmax,acc,''n'',2);');
  flo_fr = 1:2:size(line_flw,1)-1;
  flo_to = 2:2:size(line_flw,1);
  lin(:,12) = 0.5*(abs(line_flw(flo_fr,4))+abs(line_flw(flo_to,4)));   % averaged AC line flow  
  bus(:,13) = 1;   % dummy nominal voltage
  mac_con(:,[4:6,8:15,18,20:21]) = 0;  %ensure classical model  
  pst = struct_pst(bus, lin, mac_con, 100, [], 60);   % 100 MVA base in PST toolbox coh_np48.m
  mpw = pst2mp(bus, lin, mac_con, 100, 60);
else
  obj = MatpowerIn('caseid',caseid,'n_pf',1);
  mpw = mp2bgl(obj,'acdc','ac');    %mpw is used below - same 'acdc' as for pst
  pst = mp2pst(mpw,'acdc','ac');
end
if size(pst.gen,2)<20
  pst.gen(:,[4:6,8:15,18]) = 0;   %ensure classical model
else
  pst.gen(:,[4:6,8:15,18,20:21]) = 0;  
end

% Convert to graph
g = pst2graph(pst,'mode','dcpf');
m = size(g.adj, 1); 

% Obtain coherency. Normalise columns by machine inertias
[V_s, lambda, MK, K, M] = PST.pstVslow(pst, nc);
%assert(max(lambda)<0.0001);
ng = size(V_s, 1);
mac_con = pst.gen;
V_s = Utils.normalize_rows(V_s')';
V_s = V_s*abs(diag(1./diag(sqrt(V_s'*M*V_s)))); 
V_s'*M*V_s;

% Get the aligned eigenspace
[Y, JJ, ~] = GraphUtils.rotVslow(V_s, 'YuShi', [1 0 0], 2, 2);

% Get coherency similarity graph
m000 = 1:ng;
Kful = K(:,:,1); 
Kful(1:ng+1:ng*ng) = 0;
Kfsm = abs((Kful + Kful')/2);  %discard potential negative entries
Gfsm = PFgraph('adj', abs(Kfsm), 'vw', full(diag(M)));

% Get coherency groups and set them in the graph
spcopt.card_min = 1; spcopt.improv = false;  
spcopt.thrs = sqrt(2)/2+0.0001;  % the "normal" value
V_z = Y{nc-1};
cores = spec_cores3(V_z, Gfsm.vw(:), Gfsm.adj, spcopt);
T1 = greedy_ncut(Gfsm.adj, Gfsm.vw, cores, setdiff(m000,vertcat(cores{:})),true);
[T2, EXP2] = ncut_refine(Gfsm.adj, Gfsm.vw, T1, true);
[~, ia, ~] = intersect(g.bus, mac_con(:,2)');  % re-express into graph internal indexing
g.coh(1,:) = ia(:)';
g.coh(2,:) = T2';
%}

% Set negative loads as generators and negative generators as loads to
% ensure positive load and generator shedding to ensure that the lowest  
% minimum of the LS minimization problem is zero. "Negative loads" are 
% indeed actually generators (e.g., solar or wind) and their reduction 
% indeed qualifies as generator shedding (GS) => "free" LS with no LS costs
pow_gen = g.powgen(:,1);  % actual generation from lossy ACPF
pow_lod = g.powlod(:,1);
cL = ones(1,m);
cL(pow_lod>0) = 1;  % only penalize LS on the original loads
gen2lod = pow_gen<0;
lod2gen = pow_lod<0;
pow_lod(gen2lod) = pow_lod(gen2lod)-pow_gen(gen2lod);
pow_gen(gen2lod) = 0;
g.gen(gen2lod) = false;   % from here, g.gen and g.coh don't share same buses!
pow_gen(lod2gen) = pow_gen(lod2gen)-pow_lod(lod2gen);
pow_lod(lod2gen) = 0;
g.powgen(:,1) = pow_gen;
g.powlod(:,1) = pow_lod;
vw(:,1) = min(0.2*pow_gen(:,1), 1.0*pow_gen(:,1));   % min gen output assuming only online generators are available and cannot ramp up a lot
vw(:,2) = max(0.2*pow_gen(:,1), 1.0*pow_gen(:,1));   % max gen output assuming only online generators are available and cannot ramp up a lot
vw(:,3) = min(pow_lod(:,1), 0.0*pow_lod(:,1));   % min load demand
vw(:,4) = max(pow_lod(:,1), 0.0*pow_lod(:,1));   % max load demand
assert(all(vw(:,2)>=vw(:,1)) && all(vw(:,4)>=vw(:,3)));
assert(all(pow_gen>=vw(:,1)) && all(pow_gen<=vw(:,2)));
%g.powlod(setdiff(1:m,[g.coh(1,:),find(pow_lod)'])) = mean(abs(g.powlod))/1e5;  % all potential loads should demand some power flow (for graph connectedness)
g.vw = vw; 

% Create the power flow graph
g_mw = pst2graph(pst, 'mode', 'apfg');
g_mw.coh = g.coh;
adj_mw = g_mw.adj;
[i_mw,j_mw,w_mw] = find(tril(adj_mw));
p_mw = [i_mw,j_mw,w_mw];
m = size(adj_mw, 1);

% Power flow limits for each line (large uniform limits below are harder to solve)
if ~exist('p_mx', 'var')
  p_mx = [i_mw,j_mw,max(abs(w_mw))*ones(size(w_mw))];  %(this is changed in coMILPdcOPF.m anyway)
end

g_mw.coh = g.coh;
g_mw = reduce_pml(g_mw);
g_dst = PFgraph('adj', double(logical(g_mw.adj)), 'inc', g_mw.inc, 'gen',...
  g_mw.gen, 'scal', g_mw.scal, 'bus', g_mw.bus, 'coh', g_mw.coh, 'ml',...
  g_mw.ml);  %, 'vw', g.vw    
trees = hclust_constrained_path(g_dst);
T = g_mw.coGrReduClust(trees, 'refine', true);
T = T(g_mw.merge_map);  
[cut, ~, T] = final_cutset(g, T);
[~, pcut, ~, ~, num_viols, viols, tot_outl_gen] = g.ici_info(cut);
%assert(tot_outl_gen==0);
x_ini = sparse(1:m, T, 1, m, nc);
y_ini = zeros(size(p_mw));
d_ini = zeros(size(p_mw,1)*nc, 4);
for p = 1:1:size(p_mw,1)
  i = p_mw(p,1);
  j = p_mw(p,2);
  y_ini(p,1:2) = [i,j];
  d_ini((p-1)*nc+1:p*nc,1:2) = [i*ones(nc,1),j*ones(nc,1)];
  for k = 1:1:nc
    d_ini((p-1)*nc+k,3) = k;
    d_ini((p-1)*nc+k,4) = x_ini(i,k).*x_ini(j,k);
  end
  y_ini(p,3) = sum( x_ini(i,:).*x_ini(j,:) );
end

model = 'NFan.gms';  %-dicut
[T, sol_stat, res] = coMILPdcOPF(g, 'model', model, 'p_mw', p_mw, 'p_mx', p_mx,...
  'mu', 0.01, 'nu', 1.0, 'x_ini', x_ini, 'y_ini', y_ini, 'd_ini', d_ini,...  %0.033
  'mpw', mpw, 'solver', [solver,'|',PHI_OR_KVL], 'cL', cL, 'caseid', caseid);
%[eXp, pcut, ixs, ~, num_viols, viols, tot_outl_gen] = ici_info( g, cutind );
end




