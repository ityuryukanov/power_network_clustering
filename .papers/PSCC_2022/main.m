%   
% A very basic wrapper around test_mip(). Many study parameters are still
% changed inside of test_mip() and the optimization routines that it calls
% (e.g., in start_KVL.py and start_PHI.py).
% This one is only needed to run multiple cases/group numbers at once.
%   
clc;
PHI_OR_KVL = {'PHI','ZIJ','KVL'};
cases = {'case89pegase', 'case300', 'case1354pegase', 'case1888rte', 'case2868rte'};  %, 'datanp48'
k_small = 5;
k_large = 8;
for p = 1:1:numel(PHI_OR_KVL) 
  PXX = PHI_OR_KVL{p};
  eval(['res_',PXX,'=containers.Map();'])
  for i = 1:1:numel(cases)
    caseid = cases{i};
    try
      % First try to load a Matpower case
      obj = MatpowerIn('caseid',caseid,'n_pf',1);
      mpw = mp2bgl(obj,'acdc','ac');
      n   = size(mpw.bus,1);
      m   = size(mpw.gen,1);
    catch
      % Try to load a PST data file, if it fails then accept error
      eval(caseid);
      n   = size(bus,1);
      m   = size(mac_con,1);
    end
    if m<100
      K = k_small;
    else
      K = k_large;
    end
    res_caseid = NaN(K-1,9);
    for k = 2:1:K
      try
        res = test_mip(PXX, k, caseid);
        res_caseid(k-1,:) = [k,res.UB,res.runTme,res.GAP*100,res.totPLS,res.islIMB,res.totPGS,res.totCUT,res.t_feas];
      catch
        res_caseid(k-1,:) = [k,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
      end      
    end
    eval(['res_',PXX,'(caseid) = res_caseid;']);
  end
end

foo=5;