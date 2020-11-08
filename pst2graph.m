function g = pst2graph( pst, varargin )
% Syntax: 
% g = pst2graph( pst, varargin )
% 
% Purpose:
% create load flow graph (represented by adjacency and incidence matrices)
% from load flow data CONTAINED IN the pst input, and output the PFgraph 
% object g.
% 
% Input:
% pst: a pst-formatted struct
%   It contains bus, line, gen data for the input network in PST format.
%   It also includes the solution of the load flow and generator buses 
%   assigned to their respective coherent generator groups. 
% varargin: key-value pairs for this function (a struct or an isomorphic 
%   cell array). The keys and their possible values are listed below.
%                    
%   mode: what kind of graph should be created. This only affects the
%   weights of the adjacency matrix.
%   'apfg' will create an active power flow graph (default, also the 
%    legacy behaviour); Power flows are in p.u. relative to g.basmva
%   'dcpf' will create a classical DC power flow admittance matrix with
%    weights beeing line reactances in per unit on the system base.
%   'dqdv' will create a full graph with edge weights equal to sensitivities 
%    of bus reactive powers to bus voltages.
% 
% Output: 
% g: a PFgraph object
%   load flow represented as graph, with initial bus labels, coherent nodes 
%   labeling, must-link branch labeling. P and Q injections/consumptions at
%   each bus in [MVA] are stored in g.vw
%   
% Remarks:           
% A. LINE DATA FORMAT 
% 1.from_bus, 2.to_bus, 3.resistance(pu), 4.reactance(pu), 
% 5.line_charging(pu), 6.tap_ratio, 7.tap_phase, 8.tapmax, 9.tapmin,
% 10.tapstep, 11.avg_active_linflow (NEW!), 12.MUST-LINK (NEW!)
% 
% Author: Ilya Tyuryukanov
% Date: 24 July 2015
% Last revision: 12 January 2015

% Input processing
p = inputParser;
p = Utils.inputParserSetup(p);
p.addParameter('mode', 'apfg', @(x) ischar(x) && length(x)==4);
p.parse(varargin{:});
varinp = p.Results;
trash = p.Unmatched;
assert(isempty(fieldnames(trash)), [mfilename,':WrongKeyValueInput'],...
  ['[%s] Some unexpected key value pairs are detected. Please check your ',...
   'spelling. The allowed keys are: mode.'], mfilename);
mode = varinp.mode;

lin = pst.lin;
gen = pst.gen;
bus = pst.bus;
coh = pst.coh;
bus = sortrows(bus,1,'ascend');
gen = sortrows(gen,2,'ascend');

% SANITY CHECKS
assert(issorted(gen(:,2)) && issorted(bus(:,1)),...
  [ mfilename,':InputsMustBeSorted'],...
  '[%s] data.gen(:,2) and data.bus(:,1) must be sorted (for acceleration)');
assert(all(lin(:,12)>=0), [ mfilename, ':InvalidProperty'],...
  '[%s] line active power flows in data.lin(:,12) should be non-negative', mfilename);
assert(all(all(lin(:,[6])>=0)), [ mfilename, ':InvalidProperty'],...
  '[%s] tap ratios in data.lin(:,[3,4,6]) should be non-negative', mfilename);  % -> line admittance may be negative, and line resistances too (take IEEE 145 bus...)

% ELEMENTARY PREPARATIONS
% Set branches from lower to higher busID to help unique(.., 'rows')
[lin(:,1:2),idx] = sort(lin(:,1:2),2,'ascend');
% Readjust transformer ratios (fr/to <-> to/fr)
row_swp = find(idx(:,1)==2 & idx(:,2)==1);
for i = 1:1:numel(row_swp)
  if abs(lin(row_swp(i),6))>3e-16
    row_old = lin(row_swp(i),:);
    lin(row_swp(i),6) = 1/row_old(1,6);
    lin(row_swp(i),7) = -row_old(1,7);
    lin(row_swp(i),3) = row_old(1,3)*row_old(1,6)^2;
    lin(row_swp(i),4) = row_old(1,4)*row_old(1,6)^2;    
    if size(lin,2)>7
      if abs(row_old(1,9))>3e-16
        lin(row_swp(i),8) = 1/row_old(1,9);
      end
      if abs(row_old(1,8))>3e-16
        lin(row_swp(i),9) = 1/row_old(1,8);
      end
    end    
  end
end
[frto,idx] = sortrows(lin(:,1:2), [1,2]); 
lin = lin(idx,:);   % standartize the order of lin

% Remove isolated buses in bus[]
actBus = lin(:,1:2);
actBus = unique(actBus(:), 'sorted');   % all interconnected buses
idx = ismember(bus(:,1), actBus);
bus = bus(idx,:);   % remove all disconnected buses
idx = ismember(gen(:,2), actBus);
gen = gen(idx,:);

% INITIALIZATIONS
m = size(bus, 1);  % number of buses/(graph vertices)
busmap = bus(:,1).';  % busmap(i) gives original number (from the PST/MATPOWER model) of ith graph node
genmap = ismember(busmap, gen(:,2).');  % tells whether node is a generator (logical one) or a load/not-a-generator (logical zero)

switch lower(mode)
  case 'apfg'
    col_wgt = 12;  % weights are active power flows in p.u.    
    [lin] = prep_adj(lin, col_wgt);
  case 'dcpf'
    mpw = pst2mp(bus, lin, gen, pst.basmva, pst.f0, 'ac', true);
    [i2e, bus_mpw, gen_mpw, lin_mpw] = ext2int(mpw.bus, mpw.gen, mpw.branch);
    B_BUS = makeBdc(mpw.baseMVA, bus_mpw, lin_mpw);
    assert(issymmetric(B_BUS));
    B_BUS(1:m+1:m*m) = 0;
    assert(isequal(i2e(:),busmap(:)));
    % Prepare lin(:,[1,2,4]) from B_BUS:
    [i,j,w] = find(triu(B_BUS));
    frto = sort([i2e(i), i2e(j)], 2, 'ascend');
    [frto,idx] = sortrows(frto, [1,2]); 
    w = w(idx);    
    % Prepare lin(:,[1,2,4]) from PST input:
    col_wgt = 4;
    [lin] = prep_adj(lin, col_wgt);
    lin(:,1:2) = sort(lin(:,1:2), 2, 'ascend');
    [lin(:,1:2),idx] = sortrows(lin(:,1:2), [1,2]);
    lin = lin(idx,:);    
    assert(isequal(frto,lin(:,1:2)));
    lin(:,col_wgt) = -w;  % after all checks, set weights as admittances from Bdc
  otherwise
    error([mfilename,':WrongKeyValueInput'],...
     ['[%s] The optional key-value input ''mode'' can only take values ',...
      '''apfg'' and ''dcpf''.'], mfilename);
end
n = size(lin, 1);  % final number of branches/(graph edges)

% Precalculate and store adjacency matrix indices of line buses
[~, idx] = ismember(lin(:,1), busmap);
i_adj = idx(:).';  % row indices of sparce adjacency matrix; row (bus1) indices of sparce incidence matrix
[~, idx] = ismember(lin(:,2), busmap);
j_adj = idx(:).';  % col indices of sparce adjacency matrix; row (bus2) indices of sparce incidence matrix

% CONSTRUCT ADJACENCY & INCIDENCE MATRIX
wgt = lin(:,col_wgt);  % weights vector of sparse adjacency matrix 
adj = sparse([i_adj, j_adj], [j_adj, i_adj], [wgt, wgt], m, m);
[S, ~] = GraphUtils.alecconncomp(adj);  % check connectivity
assert(S == 1, [ mfilename, ':NetIsNotConnected'],...
  ['[%s] The supplied power network in the PST format is not fully ',...
  'connected!'], mfilename);
lin_idx = [i_adj; j_adj];
i_inc = lin_idx(:);
j_inc = repmat(1:n, [2, 1]);
j_inc = j_inc(:);
v_inc = ones(1, 2*n);
inc = sparse(i_inc, j_inc, v_inc, m, n);

% Convert coherent groups information from PST-output format to more compact
% two-row representation (Now pst.coh should INITIALLY be in this format)
%{
cohmap = zeros(1, m);
for i = 1:size(coh, 2)
  [~, ~, ib] = intersect(coh(:,i), bus(:,1));
  cohmap(:, ib) = i;  % this was version 1
end
[~, j, v] = find(cohmap);
n_gen = nnz(coh);
coh = zeros(2, n_gen);
coh(1,:) = j;  % and this is version 2
coh(2,:) = v;
%}

% Set default coherency (all to one group) in case coh was not provided
if isempty(coh)
  coh = sort(find(genmap));
  coh(2,:) = 1*ones(1, nnz(genmap));
else
  % convert bus ids into consequitive node numbers
  assert(issorted(coh(1,:))); % bus ids are supposed to be sorted
  coh(1,:) = find( ismember(bus(:,1),coh(1,:)) );  
end

% Extract generator and load P and Q at the buses
pg = bus(:,4);
pL = bus(:,6);
qg = bus(:,5);
qL = bus(:,7);

% Extract computed voltage angles and convert to radians
ph0 = bus(:,3)/180*pi;

% Extract actual and nominal voltages at the buses
V = bus(:,2).*bus(:,13);
Vn = bus(:,13);

% Sanity check
assert(isequal(sort(find(genmap)), sort(coh(1,:))),...
  [mfilename,':InconsistentResults'],...
  ['[%s] The computed coherency bus indices do not agree with the actual ',...
  'generator bus indices, return empty.'], mfilename);

% WRITE THE OUTPUT
lin_idx = lin_idx';
lin_ml = lin_idx(logical(lin(:,11)), :);
g = PFgraph('adj', adj, 'inc', inc, 'bus', busmap, 'gen', genmap,...
  'coh', coh, 'ml', lin_ml, 'scal', 1, 'powgen', [pg, qg], 'powlod',... 
  [pL, qL], 'vmag', V, 'vnom', Vn, 'basmva', pst.basmva, 'vang', ph0);
end


function [lin] = prep_adj(lin, col_wgt)
% 
% Given parallel lines in lin(:,1:2), this function merges such parallel
% lines into single lines with their weights equal to the sum of the
% parallel lines' weights.  The line weights should be in the colums of lin
% provided in col_wgt, which can be a scalar or a vector.
%
% The line data in lin can be e.g. from Matpower or Power System Toolbox
% (Chow-Cheung-Rogers).   
% 
% Also reassigns branches with small weights to ensure network connectivity
% Preserving the original connectivity structure for all graph types is 
% CRITICAL for many functions that expect any graph returned by PST2GRAPH()
% to reflect the actual network interconnection structure.
% 

idx = abs(lin(:,col_wgt))<100*eps('double');
if nnz(idx)>0
  warning([mfilename,':IllCondInput'],...
    ['[%s] Some of the original graph weights are lower than %f. ',...
    'Setting them to %f...'], mfilename, 100*eps('double'),...
    100*eps('double'));
end
lin(idx,col_wgt) = 100*eps('double');

% Sum up similarity weights in multicircut lines
col_all = 1:size(lin,2); 
col_dif = setdiff(col_all,[1,2,col_wgt]);
lin(:,col_dif) = 0;
lin(:,1:2) = sort(lin(:,1:2), 2, 'ascend');
[~, ia, ic] = unique(lin(:,1:2), 'rows', 'stable');
acm_wgt = accumarray(ic, lin(:,col_wgt));
lin = lin(ia,:);
lin(:,col_wgt) = acm_wgt;
end