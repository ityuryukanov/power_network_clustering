% 
% "Get Minimal Cycle Basis" (in fact, the script computes multiple sets of 
%  graph cycles).
% 

% Amend MATLAB path
old_matlabpath = path;
cleanup1 = onCleanup(@() path(old_matlabpath));
addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath'))))));

% -----------
cases = {'case39'};
n_css = numel(cases);
file_dir = fileparts(mfilename('fullpath'));
for j = 1:1:n_css,
  mcb = GraphUtils.mcb_generator(cases{j}, file_dir);
end

