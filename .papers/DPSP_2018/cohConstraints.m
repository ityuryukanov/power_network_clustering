function [pst, Q, labels] = cohConstraints(pst, nc)

% 
% Purpose:
% Create the constraint matrix needed to perform constrained Spectral
% Clustering using the method developed by Xiang Wang.
% 
% Input:
% pst: a PST struct
%   It contains (among other things) bus, line and gen data for the input
%   network in the PST format (with the solution for the load flow), base
%   frequency and possibly the generator coherency groups.
% nc: number of groups to be created. 
%
% Output:
% Q: The consraint matrix (not normalized) [needed in Flexible Spectral Clustering]
% labels: label matrix of size nVertices x k. Labels(i, j)=1 if the vertex
% i should belong to the class j [needed in transductiveBalancedKCut and
% in constrained hmetis].
% 

[V_s, lambda] = PST.pstVslow(pst, nc);
assert(max(lambda)<0.0001);
mac_con = pst.gen;

% Normalize the eigenvector columns (feature normalisation)
for i = 1:nc-1
  V_s(:,i) = V_s(:,i)/norm(V_s(:,i));
end
V_s = V_s(:, 2:nc);  % drop the first eigenvector
c = kmedoids(V_s, nc);
c = clusterdata(V_s, 'maxclust', nc, 'linkage', 'average');   

n_mach = size(V_s, 1);
coh = zeros(n_mach, nc);
for i = 1:1:n_mach
  coh(i,c(i)) = mac_con(i,2);
end

% (robust implementation) Delete all possible and "impossible" non-unique buses
% "impossible": if generators at one bus go to different coherent groups -> columnwise unique() isn't enough
tmp = coh(:);
[~, i_uniq, ~] = unique(tmp);
i_rep = setdiff(1:1:numel(tmp), i_uniq);
tmp(i_rep) = 0;
coh = reshape(tmp, size(coh));

% Delete the possible groups with no generators
empty_groups = all(coh==0,1);
coh(:, empty_groups) = [];

% Fill in pst.coh according to the standard format
n_grp = size(coh,2);
pst.coh(1,:) = sort(nonzeros(coh),'ascend');
for j = 1:1:n_grp
  i_grp = ismember(pst.coh(1,:),nonzeros(coh(:,j)));
  pst.coh(2,i_grp) = j;
end

% Create the contraints matrix for Flexible Constrained Spectral Clustering
bus = pst.bus(:,1:10);
Q = zeros(size(bus,1),size(bus,1));
all_gen = nonzeros(coh);
for i=1:1:length(all_gen)
  all_gen(i) = find( bus(:,1)==all_gen(i) );
end
Q(all_gen,all_gen) = -1;
for i = 1:1:size(coh,2)
  gen_ind = ismember( bus(:,1), nonzeros(coh(:,i)) );
  Q(gen_ind,gen_ind) = 1;
end
Q(1:size(bus,1)+1:size(bus,1)*size(bus,1)) = 0;   % make diagonal zero
Q = sparse(Q);

% Create the contraints for Balanced k-Cut and hMetis
labels = zeros(size(bus,1),nc);
for i=1:1:size(coh,2)
  x = find(coh(:,i));
  z = zeros(1,length(x));
  for j=1:1:length(x)
    z(1,j) = find(pst.bus(:,1)==coh(x(j),i));
  end
  labels(z,i) = 1;
end
labels = sparse(labels);
end

