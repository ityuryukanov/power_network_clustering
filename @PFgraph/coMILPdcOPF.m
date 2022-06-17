function [T, sol_stat, res] = coMILPdcOPF(obj, varargin)
% 
% Syntax: [T, sol_stat, res] = coMILPdcOPF(obj, varargin)
% 
% Purpose: Implements DC OPF based controlled islanding using the GAMS
%   interface with MATLAB. It provides data for the GAMS model and fetches 
%   the results back to MATLAB. Noteworthy, transformers aren't forced to
%   stay switched on.
% 
% Input: 
% obj: a (sensible) PFgraph object (with adjacency, incidence matrices and 
%   constraints). The vw entry of obj is very important here: it should 
%   contain the following data (in an n_buses x 4 matrix, where n_buses is
%   the number of nodes in the graph):
%     1) Minimal generated active power at the bus
%     2) Maximal generated active power at the bus
%     3) Minimal load active power at the bus
%     4) Maximal load active power at the bus
% varargin: key-value pairs for this function (a struct or an isomorphic 
%   cell array). The keys and their possible values are listed below.
% 
%   model: The gms model file name representing the flavour of DF OPF based
%     islanding that needs to be calculated.
%   p_mx: Active power flow limitation for branches, in units of obj.vw
%     (active and reactive loads at each bus). The format is a 3 column 
%     matrix with the rows in the following format: [from to p_mx]
%     where from stands for the "from" terminal of the branch to stands for
%     the "to" terminal of the branch p_mx stands for the power limitation 
%     of the "from"-"to" branch.
%   mpw: initial MATPOWER network data to check the islanding solution for
%     the feasibility, validity and against the constraints.
%   solver: GAMS or Gurobi.
% 
% Output:
%   T: partition vertices assigned to each island
%   sol_stat: solution statistics
%   res: Gurobi results
% 
% MAIN REFERENCES:
%   Many islanding papers....
% 
% Author: Ilya Tyuryukanov
% Date of first version: 11 September 2016
% 

% Initialisation
GRBDIR = 'cyc_dcopfici_kvl'; 
adj = obj.adj;        % abs()
assert(issymmetric(adj));
inc = obj.inc;  
m = size(adj, 1);     % number of buses
n = size(inc, 2);     % number of branches/edges
[fr, t0, b] = find(tril(adj));
assert( size(obj.vw,1)==m && size(obj.vw,2)>=4,...
  [mfilename,':WrongPublicInput'],...
  ['[%s] The vw field of the input PFgraph object obj should be at least ',...
  '%d x 4, where %d is the number of buses. See the function help for the ',...
  'format description.'], mfilename, m, m);

% Setup the basic vertex labeling based on coherency
K = unique(obj.coh(2,:));
nc = numel(K);
x_fx0 = sparse(obj.coh(1,:), obj.coh(2,:), 1, m, nc);
x_fx0 = full(x_fx0);
y_fx0 = [fr, t0, zeros(length(b), 1)];
d_fx0 = [fr, t0, ones(length(fr), 1)];
for k = 2:1:nc
  d_fx0 = [d_fx0; [fr, t0, k*ones(length(fr), 1)]];
end
d_fx0 = [d_fx0, zeros(size(d_fx0,1),1)];

% Input processing
p = inputParser;
p = Utils.inputParserSetup(p);
p.addParameter('x_ini', x_fx0, @(x) isnumeric(x) && ismatrix(x) && size(x,2)==nc && size(x,1)==m);
p.addParameter('d_ini', d_fx0, @(x) isnumeric(x) && ismatrix(x) && size(x,2)==4);
p.addParameter('y_ini', y_fx0, @(x) isnumeric(x) && ismatrix(x) && size(x,2)==3);
p.addParameter('p_mx', [fr, t0, zeros(length(b), 1)], @(x) isnumeric(x) && ismatrix(x) && size(x,2)==3);  %zero means no limitation (this will keep things simple)
p.addParameter('model', 'NFan.gms', @(x) ischar(x) && strcmp(x(end-3:end), '.gms'));
p.addParameter('mpw', struct(), @isstruct);
p.addParameter('cL', ones(1,m), @(x) isnumeric(x) && isvector(x));
p.addParameter('mu', 0, @(x) isscalar(x) && x>=0 && x<=1);
p.addParameter('nu', 1, @(x) isscalar(x) && x>=0 && x<=1);
p.addParameter('p_mw', [fr, t0, zeros(length(b), 1)], @(x) isnumeric(x) && ismatrix(x) && size(x,2)==3);
p.addParameter('solver', 'gams', @(x)ischar(x)&&(strcmp(x(1:4),'gams')||strcmp(x(1:6),'gurobi')));
p.addParameter('caseid', NaN, @(x)ischar(x));
p.parse(varargin{:});
varinp = p.Results;
trash = p.Unmatched;
assert(isempty(fieldnames(trash)), [mfilename,':WrongKeyValueInput'],...
  ['[%s] Some unexpected key value pairs are detected. Please check your ',...
   'spelling. The allowed keys are: p_mx, model.'], mfilename);
p_mx = varinp.p_mx;
gmsfile = varinp.model;
mpw = varinp.mpw;
p_mw = varinp.p_mw;
mu = varinp.mu;
nu = varinp.nu;
cL = varinp.cL;
x_ini = varinp.x_ini;
y_ini = varinp.y_ini;
d_ini = varinp.d_ini;
solver = varinp.solver;
caseid = varinp.caseid;
[solver,launcher] = strtok(solver, '|');
launcher = erase(launcher, '|');
if nargout>2 && ~strcmp(solver,'gurobi')
  error(['[%s] Unacceptable combination of input and output parameters: ',...
    'three outputs are only programmed for the gurobi solver.'], mfilename);
end
[v_fx, k_fx, x_ini] = find(x_ini);

% Make sure that edge UELS are consistent for all edge parameters
if all(p_mx(:,3)==0),
  p_mx = [fr, t0, b*sin(pi/4)];   %limit maximum line angle to 45 degrees?
end
[fr_lim, t0_lim, p_mx, adj_mx] = get_edg_attr(p_mx, m);
[fr_mw, t0_mw, p_mw, adj_mw] = get_edg_attr(p_mw, m);
[fr_ini, t0_ini, y_ini] = get_edg_attr(y_ini, m);
d_ini = sortrows(d_ini, 3);
f_dini = []; t_dini = []; d_dini = []; k_dini = [];
for k = 1:1:nc
  rows = d_ini(:,3)==k;
  [from, to, val] = get_edg_attr(d_ini(rows,[1,2,4]), m);
  f_dini = [f_dini; from(:)];
  t_dini = [t_dini; to(:)];
  d_dini = [d_dini; val(:)];
  k_dini = [k_dini; k*ones(length(val),1)];
end
d_ini = [f_dini, t_dini, k_dini, d_dini];

% Compute initial rooted arborescence z_ini (only if there is a feasible solution to coherency)
v = obj.bus;
i2e = v(:);
[~,ref_idx] = unique(obj.coh(2,:));
[v_in,ord] = sort(v_fx,'ascend');
if ( all(diff(v_in)==1) && v_in(1)==1 && v_in(m)==m && isequal(v(:),i2e(v_in(:))) )
  T = k_fx(ord);
  [cutind] = final_cutset(obj, T);
  [eXp, pcut, ixs, adj_sep, num_viols, viols, tot_outl_gen] = obj.ici_info( cutind );  
else
  tot_outl_gen = 100;  %if no initial solution was provided
  eXp = ones(2,nc+1);
end
if size(eXp,2)==nc && tot_outl_gen==0
  adj_sep = logical(adj_sep).*adj_mw;   %use power flow weights
  %edg_dir = sparse(fr_ini, t0_ini, y_ini, m, m);  %trees only spanning terminal nodes
  edg_dir = sparse(fr_ini, t0_ini, 0, m, m);  %trees spanning each node
  for i = 1:1:nc
    r = obj.coh(1,ref_idx(i));
    %gr = obj.coh(2,ref_idx(i));
    %t = setdiff(obj.coh(1, obj.coh(2,:)==gr), r);
    k = T(r);
    adj_k = adj_sep;
    adj_k(T~=k,T~=k) = 0;
    mst_k = kruskal_mst(-adj_k);
    [~,~,pred] = bfs(abs(mst_k),r);
    pred(r) = 0;
    [ff, tt, ~] = find(mst_k);
    nk = length(ff);
    %trees spanning each node
    for j = 1:1:nk
      f = ff(j);
      t = tt(j);
      if f==pred(t)
        edg_dir(f,t) = 1;
        edg_dir(t,f) = 0;        
      end
    end        
    % % trees only spanning terminal nodes
    %for j = 1:1:length(t)
    %  q = t(j);
    %  while pred(q)>0
    %    edg_dir(pred(q),q) = 1;
    %    edg_dir(q,pred(q)) = 0;
    %    q = pred(q);
    %  end
    %end
  end
  [f_zini, t_zini, z_ini] = find(edg_dir);
else
  f_zini = fr_ini;
  t_zini = t0_ini;
  z_ini = y_ini;
end

% Get current directory
old_pwd = pwd;
reset_pwd = onCleanup(@() cd(old_pwd));
file_dir = fileparts(mfilename('fullpath'));

% Amend MATLAB path
old_matlabpath = path;
reset_path = onCleanup(@() path(old_matlabpath));
addpath(fullfile(fileparts(fileparts(mfilename('fullpath')))),...
  fullfile(fileparts(fileparts(mfilename('fullpath'))),'third_party'),...
  'D:\GAMS\win64\26.1\');

% Remove the old solution (for exit if GAMS doesn't generate a new solution):
delete('soln.gdx'); 

% Extract parameters for the GAMS model
obj.coh(1,:) = v(obj.coh(1,:));   % initially coh uses internal graph indexing...
pow_gen = obj.powgen;
pow_lod = obj.powlod;
vw = obj.vw;
g = v(vw(:,1)~=0 | vw(:,2)~=0);   % buses with generation possible inside of [vw(:,1), vw(:,2)] limits
l = v(pow_lod(:,1)~=0 | pow_lod(:,2)~=0);   % pow_lod(:,1) is equivalent to vw(:,4) // MATPOWER provides no explicit loads limits
[~, ord] = sort(obj.coh(2,ref_idx), 'ascend');
ref = obj.coh(1,ref_idx);
ref = ref(ord);   % reference buses for coherent groups in ascending order

% Visualize the problem using graphviz and python igraph
%{
obj_sep = PFgraph('adj',adj_sep,'bus',obj.bus); viz = ones(1,m); viz(ismember(obj.bus,obj.coh(1,obj.coh(2,:)==1))) = 2; viz(ismember(obj.bus,obj.coh(1,obj.coh(2,:)==2))) = 3; viz(ismember(obj.bus,obj.coh(1,obj.coh(2,:)==3))) = 4; 
part_viz(obj_sep,viz);
part_viz(obj,viz);
%}

% Find the distance between group roots and other terminals:
e2i = NaN(numel(i2e),1);
e2i(i2e) = 1:numel(i2e);
adj_X = sparse([fr;t0], [t0;fr], abs([1./b;1./b]), m, m);
spopt.algname = 'dijkstra';
spopt.inf = Inf;
R = obj.coh(1,ref_idx);
dst_R = [];
for k = 1:1:numel(R)
  grp = obj.coh(2,ref_idx(k));
  T = setdiff(obj.coh(1,obj.coh(2,:)==grp), R(k));
  [d] = shortest_paths(adj_X, e2i(R(k)), spopt);
  D = d(e2i(T));
  S = R(k)*ones(numel(T),1);
  I = k*ones(numel(T),1);
  dst_R = [dst_R; [T(:), S(:), I(:), D(:)]];  % T-S swapped for directed "back-cuts" (if undirected, this plays no role)
end
D = dst_R;

% Write graph data to GAMS\ici.gdx
if strcmp(solver,'gams')
  delete(fullfile(file_dir,'GAMS\ici.gdx'));
  in1.name = 'v'; in1.type = 'set'; in1.uels = cellfun(@num2str, num2cell(v), 'UniformOutput', false);
  in2.name = 'g'; in2.type = 'set'; in2.uels = cellfun(@num2str, num2cell(g), 'UniformOutput', false);
  in3.name = 'l'; in3.type = 'set'; in3.uels = cellfun(@num2str, num2cell(l), 'UniformOutput', false);
  in4.name = 'k'; in4.type = 'set'; in4.uels = cellfun(@num2str, num2cell(K), 'UniformOutput', false);
  in5.name = 'b'; in5.type = 'parameter'; in5.form = 'sparse'; in5.val  = [i2e(fr), i2e(t0), b];
  in6.name = 'coh'; in6.type = 'set'; in6.form = 'sparse'; in6.val  = obj.coh';
  in7.name = 'e'; in7.type = 'set'; in7.val  = [i2e(fr), i2e(t0)];
  in8.name = 'pG0'; in8.type = 'parameter'; in8.form = 'full'; in8.uels = in1.uels; in8.val  = pow_gen(:,1)';
  in9.name = 'pL0'; in9.type = 'parameter'; in9.form = 'full'; in9.uels = in1.uels; in9.val  = pow_lod(:,1)';
  in10.name = 'A'; in10.type = 'set'; in10.val  = [[i2e(fr);i2e(t0)], [i2e(t0);i2e(fr)]];
  in11.name = 'DT'; in11.type = 'parameter'; in11.form = 'sparse'; in11.val  = [D(:,1), D(:,2), D(:,3), D(:,4)];
  in12.name = 'pGmin'; in12.type = 'parameter'; in12.form = 'full'; in12.uels = in1.uels; in12.val  = vw(:,1)';
  in13.name = 'pGmax'; in13.type = 'parameter'; in13.form = 'full'; in13.uels = in1.uels; in13.val  = vw(:,2)';
  in14.name = 'pLmin'; in14.type = 'parameter'; in14.form = 'full'; in14.uels = in1.uels; in14.val  = vw(:,3)';
  in15.name = 'pLmax'; in15.type = 'parameter'; in15.form = 'full'; in15.uels = in1.uels; in15.val  = vw(:,4)';
  in20.name = 'TLIM'; in20.type = 'parameter'; in20.form = 'sparse'; in20.val  = [i2e(fr_lim), i2e(t0_lim), p_mx];
  in21.name = 'r'; in21.type = 'set'; in21.uels = cellfun(@num2str, num2cell(ref), 'UniformOutput', false);
  in22.name = 'mu'; in22.type = 'parameter'; in22.form = 'full'; in22.val  = mu;
  in23.name = 'p0'; in23.type = 'parameter'; in23.form = 'sparse'; in23.val  = [i2e(fr_mw), i2e(t0_mw), abs(p_mw)];
  in24.name = 'x_ini'; in24.type = 'parameter'; in24.form = 'sparse'; in24.val  = [i2e(v_fx), k_fx, x_ini]; %in24.uels = {in1.uels, in4.uels};
  in25.name = 'y_ini'; in25.type = 'parameter'; in25.form = 'sparse'; in25.val  = [i2e(fr_ini), i2e(t0_ini), ~logical(y_ini)];    % ON lines are at zero 
  in26.name = 'd_ini'; in26.type = 'parameter'; in26.form = 'sparse'; in26.val  = [i2e(d_ini(:,1)), i2e(d_ini(:,2)), d_ini(:,3), d_ini(:,4)];
  in27.name = 'z_ini'; in27.type = 'parameter'; in27.form = 'sparse'; in27.val  = [i2e(f_zini), i2e(t_zini), z_ini];
  in28.name = 'costL'; in28.type = 'parameter'; in28.form = 'full'; in28.val  = cL; in28.uels = in1.uels;
  wgdx(fullfile(file_dir, 'GAMS','ici'), in1, in2, in3, in4, in5, in6, in7,...
    in8, in9, in10, in11, in12, in13, in14, in15, in20, in21, in22, in23,...
    in24, in25, in26, in27, in28);
else  % write data to gurobi\ici (csv files for gurobi)
  if exist(fullfile(file_dir,GRBDIR,'ici'), 'dir')==7
    rmdir(fullfile(file_dir,GRBDIR,'ici'), 's');
  end
  mkdir(fullfile(file_dir,GRBDIR,'ici'));
  fid = fopen(fullfile(file_dir,GRBDIR,'ici','caseid.txt'),'wt');
  fprintf(fid, caseid); fclose(fid);  
  csvwrite(fullfile(file_dir,GRBDIR,'ici','nodes.csv'), v(:));
  csvwrite(fullfile(file_dir,GRBDIR,'ici','genss.csv'), g(:));
  csvwrite(fullfile(file_dir,GRBDIR,'ici','loads.csv'), l(:));
  csvwrite(fullfile(file_dir,GRBDIR,'ici','islid.csv'), K(:));
  csvwrite(fullfile(file_dir,GRBDIR,'ici','cohgr.csv'), obj.coh');
  csvwrite(fullfile(file_dir,GRBDIR,'ici','islsw.csv'), ref(:));
  csvwrite(fullfile(file_dir,GRBDIR,'ici','dstTR.csv'), [D(:,1), D(:,2), D(:,4)]);
  csvwrite(fullfile(file_dir,GRBDIR,'ici','x_ini.csv'), [i2e(v_fx), k_fx, x_ini]);
  csvwrite(fullfile(file_dir,GRBDIR,'ici','z_ini.csv'), [i2e(f_zini), i2e(t_zini), z_ini]);
  csvwrite(fullfile(file_dir,GRBDIR,'ici','y_ini.csv'), [i2e(fr_ini), i2e(t0_ini), ~logical(y_ini)]);  % ON lines are at zero  
  csvwrite(fullfile(file_dir,GRBDIR,'ici','graph.csv'), [i2e(fr),i2e(t0),b]);
  csvwrite(fullfile(file_dir,GRBDIR,'ici','brlim.csv'), [i2e(fr_lim),i2e(t0_lim),p_mx]);
  csvwrite(fullfile(file_dir,GRBDIR,'ici','flow0.csv'), [i2e(fr_mw),i2e(t0_mw),abs(p_mw)]);
  csvwrite(fullfile(file_dir,GRBDIR,'ici','pG0.csv'),   [v(:),pow_gen(:,1)]);
  csvwrite(fullfile(file_dir,GRBDIR,'ici','pL0.csv'),   [v(:),pow_lod(:,1)]);
  csvwrite(fullfile(file_dir,GRBDIR,'ici','pGmin.csv'), [v(:),vw(:,1)]);
  csvwrite(fullfile(file_dir,GRBDIR,'ici','pGmax.csv'), [v(:),vw(:,2)]);
  csvwrite(fullfile(file_dir,GRBDIR,'ici','pLmin.csv'), [v(:),vw(:,3)]);
  csvwrite(fullfile(file_dir,GRBDIR,'ici','pLmax.csv'), [v(:),vw(:,4)]);
  csvwrite(fullfile(file_dir,GRBDIR,'ici','costL.csv'), [v(:),cL(:)]);  
  csvwrite(fullfile(file_dir,GRBDIR,'ici','mu.csv'),    mu(:));
  csvwrite(fullfile(file_dir,GRBDIR,'ici','nu.csv'),    nu(:));
  %csvwrite(fullfile(file_dir,GRBDIR,'ici','H.csv'),    H(:));
end

% Try to solve with GAMS
isplanar = test_planar_graph(obj.adj);   %planarity test from matlabBGL
if strcmp(solver,'gams')
  %gdxInfo('ici');
  system(['gams ', fullfile(file_dir, 'GAMS', gmsfile), ' gdx=soln.gdx ', ' lo=3 ']);
else
  if strcmp(launcher,'PHI')
    system(['python ', fullfile(file_dir, GRBDIR, 'start_PHI.py')]);
    fprintf('\n');
    fprintf('-----------------------------------------------------------------------------');
    fprintf('\n');
  elseif strcmp(launcher,'KVL')
    system(['python ', fullfile(file_dir, GRBDIR, 'start_KVL.py')]);
    fprintf('\n');
    fprintf('-----------------------------------------------------------------------------');
    fprintf('\n');
  elseif strcmp(launcher,'ZIJ')
    system(['python ', fullfile(file_dir, GRBDIR, 'start_TMP.py')]);
    fprintf('\n');
    fprintf('-----------------------------------------------------------------------------');
    fprintf('\n');    
  else
    error('');
  end
end

% Read the optimization result
if strcmp(solver,'gams')
  fmt.form = 'sparse';
  fmt.compress = true;
  fmt.name = 'x';
  bussol = rgdx(fullfile(file_dir, 'GAMS\soln.gdx'), fmt);   %y, T
  busids = cellfun(@str2num, bussol.uels{1}(bussol.val(:,1)), 'UniformOutput', true);
  prtlbl = cellfun(@str2num, bussol.uels{2}(bussol.val(:,2)), 'UniformOutput', true);
else
  x = csvread(fullfile(file_dir, GRBDIR, 'ici', 'x.csv'));
  prtidx = round(x(:,3))>=0.995;  
  busids = x(prtidx,1);
  prtlbl = x(prtidx,2);    
end
busids = busids(:);
[bus_cmn,ia,ib] = intersect(obj.bus(:), busids, 'stable');
assert( isequal(sort(bus_cmn(:)),sort(obj.bus(:))) );  % ensures all buses are there
assert( isequal(obj.bus(:),busids(ib)) );
prtlbl = prtlbl(ib);  % reorder busids same as obj.bus

if strcmp(solver,'gams')
  fmt.form = 'full';
  fmt.name = 'stats';
  stats = rgdx(fullfile(file_dir, 'GAMS\soln.gdx'), fmt);
  switch length(stats.val)
    case 2
      stats.val(3)=0;
      stats.val(4)=0;
    case 3
      stats.val(4)=0;
  end
  sol_stat.mod_time = stats.val(3);
  sol_stat.sol_time = stats.val(4);
else
  ttt = csvread(fullfile(file_dir,GRBDIR,'ici','runtime.csv'));
  sol_stat.mod_time = ttt;
  sol_stat.sol_time = ttt;    
end

% Now reset the coherency to graph node indices from graph node labels:
gen_idx = e2i(obj.coh(1,:));
obj.coh = [gen_idx(:)'; obj.coh(2,:)];   

% Check the solution and the constraints!
[cutind, ~, T] = final_cutset(obj, prtlbl);
[eXp, pcut, ixs, adj_sep, num_viols, viols, tot_outl_gen] = obj.ici_info( cutind );
assert(isequal(T, prtlbl));
%assert(tot_outl_gen==0);
%assert(size(eXp,2)==nc);
T = T(:);

% Assert that cut via obj.final_cutset(prtlbl) is equal to cut returned by GAMS in y:
[fm, to] = find(triu(adj));
edgs = [fm, to];
cut1 = cutind; cut1(cut1~=1) = 0; cut1 = any(cut1,1);  
cut1 = obj.edges2adj(cut1);  
cut1 = sort(cut1, 2, 'ascend');
cut1 = sortrows(cut1, [1,2]);
if strcmp(solver,'gams')
  fmt.name = 'y';
  fmt.form = 'sparse';
  y = rgdx(fullfile(file_dir, 'GAMS\soln.gdx'), fmt);
  idx1 = y.val(:,1);
  idx2 = y.val(:,2);
  fm = cellfun(@str2num, y.uels{1}(idx1), 'UniformOutput', true);
  to = cellfun(@str2num, y.uels{2}(idx2), 'UniformOutput', true);
else
  y = csvread(fullfile(file_dir,GRBDIR,'ici','y.csv'));
  prtidx = round(y(:,3))==1;  
  fm = y(prtidx,1);
  to = y(prtidx,2);
end
ix_fm = arrayfun(@(x)find(obj.bus==x),fm);  % node_lbl 2 node_idx
ix_to = arrayfun(@(x)find(obj.bus==x),to);
cut2 = [ix_fm(:), ix_to(:)];
cut2 = sort(cut2, 2, 'ascend');
cut2 = sortrows(cut2,1);
assert(isequal(sortrows(cut2,[1,2]), cut1), ['The ',...
  'returned cuts from two different calculations are not the same!']);

if strcmp(solver,'gurobi') && nargout>2
  UB = csvread(fullfile(file_dir,GRBDIR,'ici','UB.csv'));
  GAP = csvread(fullfile(file_dir,GRBDIR,'ici','GAP.csv'));
  islIMB = csvread(fullfile(file_dir,GRBDIR,'ici','islIMB.csv'));
  totIMB = csvread(fullfile(file_dir,GRBDIR,'ici','totIMB.csv'));
  totPGS = csvread(fullfile(file_dir,GRBDIR,'ici','totPGS.csv'));
  totPLS = csvread(fullfile(file_dir,GRBDIR,'ici','totPLS.csv'));
  totCUT = csvread(fullfile(file_dir,GRBDIR,'ici','totCUT.csv'));
  runTme = csvread(fullfile(file_dir,GRBDIR,'ici','runtime.csv'));
  t_feas = csvread(fullfile(file_dir,GRBDIR,'ici','feastime.csv'));
  res.UB = UB;
  res.GAP    = GAP;
  res.islIMB = islIMB;
  res.totIMB = totIMB;
  res.totPGS = totPGS;
  res.totPLS = totPLS;
  res.totCUT = totCUT;
  res.runTme = runTme;  
  res.t_feas = t_feas; 
end

params.GRBDIR = GRBDIR;
params.solver = solver;
params.file_dir = file_dir;
end

function [fr_cnstr, t0_cnstr, p_cnstr, adj_cnstr] = get_edg_attr(p_cnstr, m)
adj_cnstr = sparse(p_cnstr(:,1), p_cnstr(:,2), p_cnstr(:,3), m, m);
if ~any(nonzeros(tril(adj_cnstr)))||~any(nonzeros(triu(adj_cnstr)))
  adj_cnstr = adj_cnstr + adj_cnstr';
end
assert(norm((adj_cnstr - adj_cnstr'),'fro')/(nnz(adj_cnstr)+eps)<1e-4,...
  [mfilename,':BadPublicInput'],...
  ['[%s] The provided p_mx input does not always account for the both ',...
  'directions of power flow (i->j and j->i). Setting the missing entries ',...
  'to their counterparts from the available data.'], mfilename);
adj_cnstr = 0.5*(adj_cnstr + adj_cnstr');  % if nearly symmetric, ensure 100% symmetric
[fr_cnstr, t0_cnstr, p_cnstr] = find(tril(adj_cnstr));  % same (sub)set of edges as for b
end

function fl = ici_sol(mpw, buslbl, ici_lin, ref)
  % Sets up a split network in MATPOWER (e.g., to test its feasibility).
  % Estimate if input splitting solution is feasible...
  mpw.branch(:,1:2) = sort(mpw.branch(:,1:2), 2, 'ascend');
  linlbl1 = buslbl(ici_lin(:,1));
  linlbl2 = buslbl(ici_lin(:,2));
  lin_lbl = [linlbl1(:),linlbl2(:)];
  lin_lbl = sort(lin_lbl, 2, 'ascend');
  branch_on = ismember(mpw.branch(:,1:2), lin_lbl, 'rows');
  mpw.branch = mpw.branch(branch_on,:);
  old_slacks = mpw.bus(:,2)==3;
  mpw.bus(old_slacks,2) = 2;
  bus_slacks = ismember(mpw.bus(:,1), ref);
  mpw.bus(bus_slacks,2) = 3;
  warning('off','MATLAB:subscripting:noSubscriptsSpecified');
  fl = rundcpf(mpw);
  fl = rundcopf(mpw);
end