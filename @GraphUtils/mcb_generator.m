function mcbj = mcb_generator(caseid, GRBDIR)
%
% Prepares the graph and call mcb_generator.py to find various types of cycles 
% in that graph. 
%
old_warning = warning;
cleanup2 = onCleanup(@() warning(old_warning));
old_matlabpath = path;
cleanup1 = onCleanup(@() path(old_matlabpath));
addpath( fullfile(fileparts(fileparts(fileparts(mfilename('fullpath'))))),...
  fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'third_party') );
old_pwd = pwd;
cleanup3 = onCleanup(@() cd(old_pwd));
file_dir = fileparts(which(mfilename));

% ----
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
  bus(:,13) = 1;   % dummy line flow, dummy nominal voltage
  mac_con(:,[4:6,8:15,18,20:21]) = 0;  %ensure classical model
  pst = struct_pst(bus, lin, mac_con, 100, [], 60);   % 100 MVA base in PST toolbox coh_np48.m
else
  obj = MatpowerIn('caseid',caseid,'n_pf',1);
  mpw = mp2bgl(obj,'acdc','ac');   %mpw is used below - same 'acdc' as for pst
  pst = mp2pst(mpw,'acdc','ac');
end
g = pst2graph(pst,'mode','dcpf');
adj = g.adj;  %abs() (dbg!)
assert(issymmetric(adj));
[fr, t0, b] = find(tril(adj));
v = g.bus;
i2e = v(:);
g_mw = pst2graph(pst, 'mode', 'apfg');
adj_mw = g_mw.adj;
[i_mw,j_mw,w_mw] = find(tril(adj_mw));
p_mw = [i_mw,j_mw,w_mw];
m = size(adj_mw, 1);
if ~exist('p_mx', 'var')
  p_mx = 2*abs(max(w_mw))*ones(length(w_mw),1);  % a dummy maximal power limitation
  p_mx = [i_mw,j_mw,p_mx];
end
[fr_lim, t0_lim, p_mx, adj_mx] = get_edg_attr(p_mx, m);
[fr_mw, t0_mw, p_mw, adj_mw] = get_edg_attr(p_mw, m);

% Write data to gurobi\ici (csv files for gurobi)
% if exist(fullfile(GRBDIR,'ici'), 'dir')==7
%   rmdir(fullfile(GRBDIR,'ici'), 's');
% end
mkdir(fullfile(GRBDIR,'ici'));
fid = fopen(fullfile(GRBDIR,'ici','caseid.txt'),'wt');
fprintf(fid, caseid); fclose(fid);  
csvwrite(fullfile(GRBDIR,'ici','nodes.csv'), v(:));
csvwrite(fullfile(GRBDIR,'ici','graph.csv'), [i2e(fr),i2e(t0),b]);
csvwrite(fullfile(GRBDIR,'ici','brlim.csv'), [i2e(fr_lim),i2e(t0_lim),p_mx]);
csvwrite(fullfile(GRBDIR,'ici','flow0.csv'), [i2e(fr_mw),i2e(t0_mw),abs(p_mw)]);
addpath(file_dir);
file_dir = strrep(file_dir,'\','\\');
if count(py.sys.path,'') == 0
  insert(py.sys.path,int32(0),'');
end
if count(py.sys.path,file_dir) == 0
  insert(py.sys.path,int32(0),file_dir);
end

mcbj = py.mcb_generator.mcbs(fullfile(GRBDIR,'ici'));
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


