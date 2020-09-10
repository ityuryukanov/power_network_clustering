

% 
% Reference: I. Tyuryukanov, M. Popov, M. A. M. M. van der Meijden, and V. 
% Terzija. "Slow Coherency Identification and Power System Dynamic Model 
% Reduction by using Orthogonal Structure of Electromechanical Eigenvectors". 
% In: IEEE Trans. Power Syst. 33.6 (2020), pp. ???????????.
% 

%% Adjust path
addpath(fileparts(fileparts(fileparts(which(mfilename)))));
addpath(fullfile(fileparts(fileparts(fileparts(which(mfilename)))),'third_party'));
addpath('netw_data');

%% Extract the electromechanical model of NPCC48 test system
format compact;   
clear;            
close all;        
clearvars -global;
warning off;
rng('default');
x_width = 241; y_width = 161;
allsym  = false;   % use the original inv(M)K matrix or inv(M)Ksym, where Ksym is symmetrized K
dynsim  = true;    % simulate or not the aggregated system in time-domain?
PST.pst_var;       % initialize the global PST system state  

%% Load the network
basemva = 100;  %system base for all networks
%datane; f0 = 60; mac0ref = 10; bus_f = 21; bus_t = 22; nc = Inf; kmax = 9; studyarea = 1:4;
load('datanp48.mat'); f0 = 60; nc = 16; kmax = 16; studyarea = 1:9; mac0ref = 48; genID = 4; bus_f = 7; bus_t = 6;
%load('datanp48.mat'); f0 = 60; nc =  3; kmax = 18; studyarea = 1:9; mac0ref = 43; genID = 4; bus_f = 7; bus_t = 6;
%d2adce; f0 = 60; nc = 2; kmax = 4; bus(1,10) = 2; bus(6,10) = 1; studyarea = 3:4; mac0ref = 4; genID = 3; bus_f = 101; bus_t = 3; 
%data3m9b; f0 = 60; nc = 3; kmax = 3;
%load('data16m.mat'); f0 = 60; nc = 5; kmax = 12;
%load('data50m.mat'); f0 = 60; kmax = 18;
[nrow,ncol] = size(mac_con);
if ncol <= 21,
  mac_con = [mac_con, ones(nrow,23-ncol)];
end
mac_con(:,[4:6,8:15,18,20:21]) = 0;  %ensure classical model
mac_con(:,19) = mac_con(:,1);

%% Optionally aggregate machines sitting on the same bus (the initial approach)
%{
samebus = tabulate(mac_con(:,2));
samebus = samebus(samebus(:,2)>1,1);
bus_sar = mac_con(ismember(mac_con(:,1),studyarea(:)),2);
bus_ref = mac_con(mac_con(:,1)==mac0ref,2);
bus_tst = mac_con(mac_con(:,1)==genID,2);
for i = 1:1:numel(samebus)
  b = samebus(i);
  mac_lst = find(mac_con(:,2)==b);
  agg_mac = PST.eqgen(mac_con,mac_lst,basemva,b,mac_con(mac_lst(1),1));  
  agg_mac(1,22:23) = 1;
  mac_con(mac_lst,:) = [];
  mac_con = [mac_con;agg_mac];
end
[~,ord] = sort(mac_con(:,2),'ascend');  %a new meaninful machine order
mac_con = mac_con(ord,:);
mac_con(:,1) = (1:1:size(mac_con,1))';
studyarea = mac_con(ismember(mac_con(:,2),bus_sar),1);
mac0ref = mac_con(mac_con(:,2)==bus_ref,1);
genID = mac_con(mac_con(:,2)==bus_tst,1);
%}

%% Re-order the data and run the loadflow
[mac_con(:,1),I] = sort(mac_con(:,1),'ascend');   %to easier map to eigenvector rows
mac_con(:,2:end) = mac_con(I,2:end);
tol = 2e-13; iter_max = 30;
vmin = 0.5;  vmax = 1.5; acc = 1.0; 
[bus(:,1),idx] = sort(bus(:,1),'ascend');   
bus(:,2:end) = bus(idx,2:end);
bus0 = bus; line0 = line; mac_con0 = mac_con;
[DUMMY, bus, ~] = evalc('PST.LLoadflow(bus0,line0,tol,iter_max,vmin,vmax,acc,''n'',2);');

%% Form the electromechanical model and compute its right eigenvectors
line(:,12) = eps;      %dummy line flow 
bus(:,13) = 1;         %dummy nominal voltage
m = size(mac_con,1);
pst = struct_pst(bus, line, mac_con, basemva, [], f0);
[V, lambda, MK, K, M] = PST.pstVslow(pst, m);
lambda = sort(sqrt(abs(lambda(:))));
V = V(:,1:kmax);
Ksym = K(:,:,2);

% Symmetric K(:,:,2); V and lambda from inv(M)K(:,:,2); eXp and Ncut from K(:,:,2)
if allsym
  MK = M\Ksym;       
  [V, lam] = eig(MK);
  lam = diag(lam); [~,idx] = sort(real(lam),'descend'); lam = lam(idx);
  V = V(:,idx); lam = sort(sqrt(abs(lam(:)))); V = V(:,1:kmax);
  K(:,:,1) = Ksym; 
  lambda = lam; 
end
V = Utils.normalize_rows(V')';   %Normalize V to length 1
Z = V*abs(diag(1./diag(sqrt(V'*M*V))));   %Normalize V with inertias

%% Defining various synchronizing torque graphs
Ksym(1:m+1:m*m) = 0;
Kful = K(:,:,1);
Kful(1:m+1:m*m) = 0;
Kfsm = (Kful + Kful')/2;
assert(mean(mean(abs(Kfsm-Ksym)))<1e-6);   %demonstrate Kfsm==Ksym for NCut
Gsyn = PFgraph('adj', abs(Ksym), 'vw', full(diag(M)));
Gfsm = PFgraph('adj', abs(Kfsm), 'vw', full(diag(M)));

%% Various rotated eigenvectors
[Y, JJ, ~] = GraphUtils.rotVslow(Z, 'YuShi', [1 0 0], 2, 2); 
[Y_ZM, J_ZM, t_ZM] = GraphUtils.rotVslow(Z, 'ZM', [1 0 0], 2, 2); 
[X, YS, ~] = GraphUtils.nc_ev_sym(Z, [1 1], [0 0 0 1], [], 2, true); Q = ev_qt(X);
spcopt.card_min = 1; spcopt.improv = false;  spcopt.thrs = sqrt(2)/2+0.0001;

%% Various ways to identify slow eigensubspace
eigIVN0 = zeros(1,kmax-1); ncutVN0 = zeros(1,kmax-1);  % Classic slow coherency algorithm
eigIVN2 = zeros(1,kmax-1);                             % Repeated classic slow coherency algorithm
eigIZN1 = zeros(1,kmax-1); ncutZN1 = zeros(1,kmax-1);  % Simple rounding of aligned rotated electromechanical eigenvectors + refinement
eigIZN2 = zeros(1,kmax-1); ncutZN2 = zeros(1,kmax-1);  % New alignment + greedy assignment 
eigIZN3 = zeros(1,kmax-1); ncutZN3 = zeros(1,kmax-1);  % New alignment + greedy assignment + refinement
eigIZN4 = zeros(1,kmax-1); ncutZN4 = zeros(1,kmax-1);  % "Discovering Clusters in Power Networks From Orthogonal Structure of Spectral Embedding"
eigIZN6 = zeros(1,kmax-1); ncutZN6 = zeros(1,kmax-1);  % Closest distance to cores + refinement
eigIZN7 = zeros(1,kmax-1); ncutZN7 = zeros(1,kmax-1);  % New alignment + random assignment + refinement (10 times?)
%eigIZN8 = zeros(1,kmax-1); ncutZN8 = zeros(1,kmax-1); % Exact ncut
%% Start the loop
mo000 = 1:m;
for k = 2:1:kmax
  %% Run various slow coherency identification methods
  % Classic slow coherency algorithm
  Vs = V(:,1:k);
  [aaVN0, naVN0, LdVN0, moVN0] = PST.L_group(Vs);
  T_VN0 = U2T( logical([eye(k);discretize(LdVN0)]), moVN0 );
  [EXP_VN0, ixs_VN0] = GraphUtils.lbl_info(Gfsm.adj, Gfsm.vw(:), T_VN0);
  % Repeated classic slow coherency algorithm
  [aaVN2, naVN2, LdVN2, moVN2, LdLgN] = PST.LRgroup(Vs);
  % "Greedy Min Ncut using orthogonal eigenvectors"
  spcopt.thrs = 0.707; spcopt.improv = false;
  cores_ZN2 = spec_cores3(Y{k-1}, Gfsm.vw(:), Gfsm.adj, spcopt);  %, 1:m
  clear('spec_cores3');  % clears persistent variables used for plotting
  [T_ZN2, EXP_ZN2] = greedy_ncut(Gfsm.adj, Gfsm.vw, cores_ZN2, setdiff(mo000,vertcat(cores_ZN2{:})), true); [aaZN2, naZN2] = T2A(T_ZN2);
  [EXP_ZN2, ixs_ZN2] = GraphUtils.lbl_info(Gfsm.adj, Gfsm.vw(:), T_ZN2);
  % "Greedy Min Ncut using orthogonal eigenvectors + refinement"
  [T_ZN3, EXP_ZN3] = ncut_refine(Gfsm.adj, Gfsm.vw, T_ZN2, true); [aaZN3, naZN3] = T2A(T_ZN3);
  [EXP_ZN3, ixs_ZN3] = GraphUtils.lbl_info(Gfsm.adj, Gfsm.vw(:), T_ZN3);
  % "Old algorithm using row-normalized aligned eigenvectors"
  spcopt.thrs = 0.707; spcopt.improv = true;
  cores_ZN4 = spec_cores2(X{k-1}, Gfsm.vw(:), Gfsm.adj, spcopt);  %, 1:m
  clear('spec_cores2');  % clears persistent variables used for plotting
  T_ZN4 = Gfsm.uckwaycut(cores_ZN4, 2, 1); [aaZN4, naZN4] = T2A(T_ZN4);
  [EXP_ZN4, ixs_ZN4] = GraphUtils.lbl_info(Gfsm.adj, Gfsm.vw(:), T_ZN4);
  % "Simple rounding of rotated mode shapes + refinement" (ANOTHER BAD ALTERNATIVE)
  [~, T_T] = max(Y{k-1},[],2);
  [T_ZN1, EXP_ZN1] = ncut_refine(Gfsm.adj, Gfsm.vw, T_T, true); [aaZN1, naZN1] = T2A(T_ZN1);
  [EXP_ZN1, ixs_ZN1] = GraphUtils.lbl_info(Gfsm.adj, Gfsm.vw(:), T_ZN1);
  % "Closest distance to cores + refinement" (ANOTHER BAD ALTERNATIVE)
  T_T = dist2cores(Y{k-1}, cores_ZN2);
  [T_ZN6, EXP_ZN6] = ncut_refine(Gfsm.adj, Gfsm.vw, T_T, true); [aaZN6, naZN6] = T2A(T_ZN6);
  [EXP_ZN6, ixs_ZN6] = GraphUtils.lbl_info(Gfsm.adj, Gfsm.vw(:), T_ZN6); 
  % Min Ncut using orthogonal eigenvectors + random assignment + refinement (10 times?)
  KZN7=17; ncut_ZN7=Inf(1,KZN7); T_ZN7=zeros(m,KZN7);
  rst=setdiff(mo000,vertcat(cores_ZN2{:})); nr=numel(rst);
  for i = 1:1:k, T_ZN7(cores_ZN2{i},:) = i; end
  T_ZN7(:,1) = T_ZN2; T_ZN7(:,2) = T_ZN1;
  for j = 3:1:KZN7, T_ZN7(rst,j) = randi(k,nr,1); end
  for i = 1:1:KZN7,
    [T_ZN7(:,i), EXP] = ncut_refine(Gfsm.adj, Gfsm.vw, T_ZN7(:,i), true); ncut_ZN7(i) = mean(EXP(1,:));
  end
  [best, best_ZN7] = min(ncut_ZN7);
  T_ZN7 = T_ZN7(:,best_ZN7);
  [EXP_ZN7, ixs_ZN7] = GraphUtils.lbl_info(Gfsm.adj, Gfsm.vw(:), T_ZN7); [aaZN7, naZN7] = T2A(T_ZN7);
  %[X_ZN8, ncutZN8(k-1)] = GraphUtils.mipncut(Kfsm, M, k, T_ZN7, 60*60); X_ZN8(abs(X_ZN8)<1.01e-3) = 0;
  %[I,T_ZN8,~] = find(X_ZN8); [~,I] = sort(I,'ascend'); T_ZN8 = T_ZN8(I); [aaZN8, naZN8] = T2A(T_ZN8);
  %--------
  ncutVN0(k-1) = mean(EXP_VN0(1,:)); ncutZN1(k-1) = mean(EXP_ZN1(1,:));
  ncutZN2(k-1) = mean(EXP_ZN2(1,:)); ncutZN3(k-1) = mean(EXP_ZN3(1,:));
  ncutZN4(k-1) = mean(EXP_ZN4(1,:)); ncutZN6(k-1) = mean(EXP_ZN6(1,:));
  ncutZN7(k-1) = mean(EXP_ZN7(1,:));
  ncut_currr = [ncutVN0(k-1),ncutZN1(k-1),ncutZN2(k-1),ncutZN3(k-1),ncutZN4(k-1),ncutZN6(k-1),ncutZN7(k-1)];
  Wexp_currr = [max(EXP_VN0(1,:)),max(EXP_ZN1(1,:)),max(EXP_ZN2(1,:)),max(EXP_ZN3(1,:)),max(EXP_ZN4(1,:)),max(EXP_ZN6(1,:)),max(EXP_ZN7(1,:))];
  %assert((ncutZN7(k-1)+1e-6)>ncutZN8(k-1))
  % Test the transformation (for theory only...)
  testS(mac_con,MK,M,LdVN0,moVN0);
  
  TR_VN0 = PST.sc_xform(aaVN0,naVN0,mac_con);
  A_VN0 = TR_VN0*MK/TR_VN0; A11_VN0 = A_VN0(1:k,1:k); A22_VN0 = A_VN0(k+1:end,k+1:end);
  TR_VN2 = PST.sc_xform(aaVN2,naVN2,mac_con);
  A_VN2 = TR_VN2*MK/TR_VN2; A11_VN2 = A_VN2(1:k,1:k); A22_VN2 = A_VN2(k+1:end,k+1:end);
  TR_ZN1 = PST.sc_xform(aaZN1,naZN1,mac_con);
  A_ZN1 = TR_ZN1*MK/TR_ZN1; A11_ZN1 = A_ZN1(1:k,1:k); A22_ZN1 = A_ZN1(k+1:end,k+1:end);
  TR_ZN2 = PST.sc_xform(aaZN2,naZN2,mac_con);
  A_ZN2 = TR_ZN2*MK/TR_ZN2; A11_ZN2 = A_ZN2(1:k,1:k); A22_ZN2 = A_ZN2(k+1:end,k+1:end);
  TR_ZN3 = PST.sc_xform(aaZN3,naZN3,mac_con);
  A_ZN3 = TR_ZN3*MK/TR_ZN3; A11_ZN3 = A_ZN3(1:k,1:k); A22_ZN3 = A_ZN3(k+1:end,k+1:end);
  TR_ZN4 = PST.sc_xform(aaZN4,naZN4,mac_con);
  A_ZN4 = TR_ZN4*MK/TR_ZN4; A11_ZN4 = A_ZN4(1:k,1:k); A22_ZN4 = A_ZN4(k+1:end,k+1:end);
  TR_ZN6 = PST.sc_xform(aaZN6,naZN6,mac_con);
  A_ZN6 = TR_ZN6*MK/TR_ZN6; A11_ZN6 = A_ZN6(1:k,1:k); A22_ZN6 = A_ZN6(k+1:end,k+1:end);
  TR_ZN7 = PST.sc_xform(aaZN7,naZN7,mac_con);
  A_ZN7 = TR_ZN7*MK/TR_ZN7; A11_ZN7 = A_ZN7(1:k,1:k); A22_ZN7 = A_ZN7(k+1:end,k+1:end);
  %TR_ZN8 = PST.sc_xform(aaZN8,naZN8,mac_con);
  %A_ZN8 = TR_ZN8*MK/TR_ZN8; A11_ZN8 = A_ZN8(1:k,1:k); A22_ZN8 = A_ZN8(k+1:end,k+1:end);
  [eigIVN0(k-1),errsAVN0] = eigErrr(A11_VN0,lambda(1:k));
  [eigIVN2(k-1),errsAVN2] = eigErrr(A11_VN2,lambda(1:k));
  [eigIZN1(k-1),errsAZN1] = eigErrr(A11_ZN1,lambda(1:k));
  [eigIZN2(k-1),errsAZN2] = eigErrr(A11_ZN2,lambda(1:k));
  [eigIZN3(k-1),errsAZN3] = eigErrr(A11_ZN3,lambda(1:k));
  [eigIZN4(k-1),errsAZN4] = eigErrr(A11_ZN4,lambda(1:k));
  [eigIZN6(k-1),errsAZN6] = eigErrr(A11_ZN6,lambda(1:k));
  [eigIZN7(k-1),errsAZN7] = eigErrr(A11_ZN7,lambda(1:k));
  %[eigIZN8(k-1),errsAZN8] = eigErrr(A11_ZN8,lambda(1:k));
  
  %% Dynamic simulations of various topologies for a given number of clusters
  if k==nc && dynsim
    % Visualize the original network
    %{
    pst0 = struct_pst(bus0, line0, mac_con0, basemva, [], f0);
    g0 = pst2graph(pst0,'mode','dcpf'); T_VZ = zeros(1, size(bus0,1));
    for j = 1:1:size(aaZN7,1),
      T_VZ(mac_con0(ismember(mac_con0(:,1),nonzeros(aaZN7(j,:))),2)) = j+1;
    end
    part_viz(g0, T_VZ);
    %}

    %Compare nonlinear inertial aggregate with the linear one (i.e., A_ZN7)
    [red_IAG.bus, red_IAG.lin, red_IAG.gen] = PST.mi_agg(bus0, line0, aaZN7, sum(logical(aaZN7),2), basemva);
    [DUMMY,red_IAG.bus,~] = evalc('PST.LLoadflow(red_IAG.bus,red_IAG.lin,tol,iter_max,vmin,vmax,acc,''n'',2);');
    pst_tst = struct_pst(red_IAG.bus, red_IAG.lin, red_IAG.gen, basemva, [], f0); 
    [~, lam_red, MK_RED, K_RED, M_RED] = PST.pstVslow(pst_tst, k);
    [err,errs] = eigErrr(MK_RED,lambda(1:nc));

    % Remove the study area
    areas = aaZN7; %
    [~,rowSTD] = ismember(areas,studyarea); rowSTD = any(rowSTD,2);
    arstd = areas(rowSTD,:);
    areas(rowSTD,:) = []; 
    [~,rowDEL] = ismember(aaVN0,studyarea); rowDEL = any(rowDEL,2);
    aaVN0(rowDEL,:) = [];
    
    % Simulate the full system model
    t_flt = 0.10; t_clr = t_flt + 1/60*6; t_end = 15.0;
    [ord_all, t_all, ang_all, spd_all] = PST.SIM_PST(pst, mac0ref, bus_f, bus_t, t_flt, t_clr, t_end);    
    
    % Plot coherent swings in each area
    %{
    figure; h1 = gca;
    set(gcf, 'pos', [100 100 x_width y_width]);
    figure; h2 = gca;
    areas = [areas;arstd];
    delar = sum(logical(areas),2)<=1;
    aaRED = areas; aaRED(delar,:) = [];
    dtcoi = zeros(size(aaRED,1),numel(t_all));
    spcoi = zeros(size(aaRED,1),numel(t_all));
    for i = 1:1:size(aaRED,1)      
      ixG_cur = ismember(mac_con0(:,1),nonzeros(aaRED(i,:)));
      H = mac_con0(ixG_cur,16).*mac_con0(ixG_cur,3)/basemva; Ha = sum(H);
      ang_agg = bsxfun(@times,ang_all(ixG_cur,:),H);
      ang_agg = sum(ang_agg,1)/Ha;
      spd_agg = bsxfun(@times,spd_all(ixG_cur,:),H);
      spd_agg = sum(spd_agg,1)/Ha;
      dtcoi(i,:) = ang_agg;
      spcoi(i,:) = spd_agg;
      leganga = plot(h1,t_all,ang_all(ixG_cur,:)*180/pi,'Color',rgb('Aqua'),'LineWidth',0.5); hold(h1,'on');
      legangc = plot(h1,t_all,ang_agg*180/pi,'k-','LineWidth',1.5); hold(h1,'off');
      xlabel(h1,'Time, [s]'); ylabel(h1,'\delta, [deg]'); xlim([0,t_end]);
      set(h1,'XTick',0:1:15); %set(h1,'YTick',-80:20:80);
      legend([leganga(1),legangc],'\delta_i', '\delta_{COI}', 'Box', 'off');      
      set(h1,'ygrid','on')      
      plot(h2,t_all,spd_all(ixG_cur,:),'b-','LineWidth',1.5); hold(h2,'on');
      plot(h2,t_all,spd_agg,'r-','LineWidth',1.5); hold(h2,'off');
      xlabel(h2,'Time, [s]'); ylabel(h2,'\omega, [rad/s]'); xlim([0,t_end]);
      set(h2,'XTick',0:1:15); %set(h2,'YTick',-80:20:80);
      legend(h2,'\omega_i', '\omega_{COI}');
    end
    for i = 1:1:size(aaRED,1)      
      plot(h1,t_all,dtcoi(i,:),'b-','LineWidth',1.5); hold(h1,'on');
      xlabel(h1,'Time, [s]'); ylabel(h1,'\delta, [deg]'); xlim([0,t_end]);
      set(h1,'XTick',0:1:15); %set(h1,'YTick',-80:20:80);
      plot(h2,t_all,spcoi(i,:),'b-','LineWidth',1.5); hold(h2,'on');
      xlabel(h2,'Time, [s]'); ylabel(h2,'\omega, [rad/s]'); xlim([0,t_end]);
      set(h2,'XTick',0:1:15); %set(h2,'YTick',-80:20:80);
      %%legend('all angles', 'agg. angle');
    end     
    %}
    
    % Produce reduced inertial aggregates
    [red_IAG.bus, red_IAG.lin, red_IAG.gen] = PST.mi_agg(bus0, line0, areas, sum(logical(areas),2), basemva);
    red_IAG.basmva = basemva; red_IAG.f0 = f0;
    red_ZHU = PST.zhukovpp(pst, areas, k, false, true);
    
    [red_IVN.bus, red_IVN.lin, red_IVN.gen] = PST.mi_agg(bus0, line0, aaVN0, sum(logical(aaVN0),2), basemva); %dbg
    red_IVN.basmva = basemva; red_IVN.f0 = f0;
    red_ZVN = PST.zhukovpp(pst, aaVN0, k, false, true); %dbg    
    %}
            
    ixG_all = mac_con0(:,1)==genID; ang_all = ang_all/pi*180;
    mac_ref = find(red_ZHU.gen(:,1)==mac0ref);
    [ord_VN0, t_VN0, ang_VN0, spd_VN0] = PST.SIM_PST(red_ZHU, mac_ref, bus_f, bus_t, t_flt, t_clr, t_end);
    mac_ref = find(red_IAG.gen(:,1)==mac0ref);
    [ord_IAG, t_IAG, ang_IAG, spd_IAG] = PST.SIM_PST(red_IAG, mac_ref, bus_f, bus_t, t_flt, t_clr, t_end);
    ixG_VN0 = red_ZHU.gen(:,1)==genID; ixG_IAG = red_IAG.gen(:,1)==genID;
    ang_VN0 = ang_VN0/pi*180; ang_IAG = ang_IAG/pi*180;
    
    mac_ref = find(red_ZVN.gen(:,1)==mac0ref); %dbg
    [ord_ZVN, t_ZVN, ang_ZVN, spd_ZVN] = PST.SIM_PST(red_ZVN, mac_ref, bus_f, bus_t, t_flt, t_clr, t_end); %dbg
    mac_ref = find(red_IVN.gen(:,1)==mac0ref); %dbg
    [ord_IVN, t_IVN, ang_IVN, spd_IVN] = PST.SIM_PST(red_IVN, mac_ref, bus_f, bus_t, t_flt, t_clr, t_end); %dbg
    ixG_ZVN = red_ZVN.gen(:,1)==genID; ixG_IVN = red_IVN.gen(:,1)==genID;   %dbg
    ang_ZVN = ang_ZVN/pi*180; ang_IVN = ang_IVN/pi*180;   %dbg
    %}
        
    figure;
    subplot(2,1,1);
    plot(t_all, ang_all(ixG_all,:), 'b-', 'LineWidth', 1.5); hold on;
    plot(t_ZVN, ang_ZVN(ixG_ZVN,:), '-.', 'Color', rgb('ForestGreen'), 'LineWidth', 1.5); hold on; %dbg
    plot(t_VN0, ang_VN0(ixG_VN0,:), 'r--', 'LineWidth', 1.5); hold off;    
    xlabel('Time, [s]'); ylabel('\delta, [deg]'); xlim([0,t_end]);    
    set(gca,'XTick',0:1:15); set(gca,'YTick',-80:20:80);    
    legend('Full Model', 'Classic Grouping & Alg. 2', 'Our Grouping & Alg. 2');
    subplot(2,1,2);
    plot(t_all, ang_all(ixG_all,:), 'b-', 'LineWidth', 1.5); hold on;
    plot(t_IVN, ang_IVN(ixG_IVN,:), '-.', 'Color', rgb('ForestGreen'), 'LineWidth', 1.5); hold on; %dbg
    plot(t_IAG, ang_IAG(ixG_IAG,:), 'r--', 'LineWidth', 1.5); hold off;   
    xlabel('Time, [s]'); ylabel('\delta, [deg]'); xlim([0,t_end]);    
    set(gca,'XTick',0:1:15); set(gca,'YTick',-80:20:80);    
    legend('Full Model', 'Classic Grouping & IAGG', 'Our Grouping & IAGG');
    %}   
           
    % Visualize the reduced network
    %{
    pstZHU = struct_pst(red_ZHU.bus, red_ZHU.lin, red_ZHU.gen, 100, [], 60);
    gZHU = pst2graph(pstZHU); gZHU_py = gZHU.pfgraph2igraph(); py.igraphmod.igraph2graphviz(gZHU_py, 'agrg_netw');
    GraphUtils.graphviz(pwd);
    %}    
  end
end
%% Various eigenvalue/spectralclustering plots...
%{
figure;
plot(JJ, 'o-', 'Color', 'k', 'MarkerSize', 3, 'LineWidth', 1.2);
xlabel('Nr. of groups, [-]', 'Interpreter', 'Latex', 'FontSize', 11);
ylabel('$J^*,~[-]$', 'Interpreter', 'Latex', 'FontSize', 12);
xlim([1 kmax-1]);
set(gca, 'XTick', 1:2:kmax-1);
set(gca, 'XTickLabel', cellfun(@num2str, num2cell(2:2:kmax), 'unif', false));
ax = gca; ax.GridColor = [0, 0, 0]; set(ax, 'xgrid', 'on'); set(gca, 'xminorgrid', 'off'); ax.GridAlpha = 0.25;
set(gcf, 'pos', [100 100 x_width y_width]);
close;
%}
hf = figure;
h1 = plot(eigIVN0, 'o-.', 'Color', rgb('Green'), 'MarkerSize', 4, 'LineWidth', 1.2); hold on;
h3 = plot(eigIZN3, '.-', 'Color', 'r', 'MarkerSize', 12, 'LineWidth', 1.2);
h2 = plot(eigIZN7, '.-', 'Color', 'k', 'MarkerSize', 12, 'LineWidth', 1.2);
%h4 = plot(eigIZN4, 'o--', 'Color', 'b', 'MarkerSize', 4, 'LineWidth', 1.2);
xlabel('Nr. of groups, [-]', 'Interpreter', 'Latex', 'FontSize', 13);
ylabel('$\overline{\delta\lambda}$, [\%]', 'Interpreter', 'Latex', 'FontSize', 14);
xlim([1 kmax-1]);
set(gca, 'XTick', 1:1:kmax-1);
set(gca, 'XTickLabel', cellfun(@num2str, num2cell(2:1:kmax), 'unif', false));
%set(gca, 'yscale', 'log');
%set(gca, 'YTick', round((10*ones(1,8)).^(linspace(0.7,2,8))));
legend([h1,h2,h3], {'Classic','Refined','Greedy'}, 'Box', 'off', 'Location', 'Best', 'FontSize', 8);  % h4--'[17]'
ax = gca; ax.GridColor = [0, 0, 0]; set(ax, 'ygrid', 'on'); set(gca, 'yminorgrid', 'off'); ax.GridAlpha = 0.25;
set(gcf, 'pos', [100 100 x_width y_width]);
close;
%ncutVN0(2) = 4.0556*1.35; ncutZN4(1) = 1.1455*1.4;
h1 = bar(ncutVN0./ncutZN7, 1.0, 'FaceColor', rgb('PaleGreen'), 'EdgeColor', 'k'); hold on;
h2 = bar(ncutZN3./ncutZN7, 0.67, 'FaceColor', rgb('LightSalmon'), 'EdgeColor', 'k');
%h4 = bar(ncutZN4./ncutZN7, 0.33, 'FaceColor', 'm', 'EdgeColor', 'k');
ylim([0.95 inf]);
plot([0 kmax], [1 1], 'k-', 'LineWidth', 1.5);
%set(gca, 'YTick', [0.95:0.05:1.3,1.35,1.4]);
%set(gca, 'YTickLabel', cellfun(@num2str, num2cell([0.95:0.05:1.3,2.3,4.15]), 'unif', false));
xlabel('Nr. of groups, [-]', 'Interpreter', 'Latex', 'FontSize', 13);
ylabel('$\mathrm{Ncut_M} / \mathrm{Ncut_M}^{\dagger}$, [-]', 'Interpreter', 'Latex', 'FontSize', 12);
xlim([1-0.5 kmax-1+0.5]);
set(gca, 'XTick', 1:1:kmax-1);
set(gca, 'XTickLabel', cellfun(@num2str, num2cell(2:1:kmax), 'unif', false));
legend([h1,h2], {'Classic','Our'}, 'Box', 'off', 'Location', 'Best');  % h4--'[17]'
ax = gca; ax.GridColor = [0, 0, 0]; set(ax, 'ygrid', 'on'); set(gca, 'yminorgrid', 'off'); ax.GridAlpha = 0.25;
set(gcf, 'pos', [100 100 x_width y_width]);  
close;


%% Internal functions
function Y = discretize(V)
  [n,k]=size(V);
  [Maximum,J]=max(V,[],2);
  Y = zeros(n,k);
  idx = sub2ind([n,k], 1:n, J');
  Y(idx) = 1;
end

function T = dist2cores(Y, cores)
m = size(Y,1);
nc = numel(cores);
T = zeros(m,1);
nodes_rst = setdiff(1:m, vertcat(cores{:}));
for k = 1:1:nc
  T(cores{k}) = k;
end
nr = numel(nodes_rst);
for i = 1:1:nr
  d = Inf(nc,1);
  for k = 1:1:nc
    dY = bsxfun(@minus, Y(cores{k},:), Y(nodes_rst(i),:));
    dD = sqrt(sum(dY.*dY,2));
    d(k) = mean(dD);
  end
  [xx,T(nodes_rst(i))] = min(d);
end
end

function [err,errs] = eigErrr(Ared,lambda)
  lambda = lambda(:);
  idx_del = lambda<1e-6;  %full system evals could be more robust than Ared's
  lambda(idx_del) = [];
  eigA = sort(sqrt(abs(eig(Ared))),'ascend');
  eigA(idx_del) = [];
  lambda = lambda/2/pi;
  eigA = eigA/2/pi;
  errs = abs(eigA-lambda)./min([eigA,lambda],[],2)*100;
  err = mean( errs );
end

function testS(mac_con,MK,M,Ld,mo)
  % To show that PST.sc_xform computes exactly the slow-fast transformation from Winkelman 1981/Avramovic 1980
  % Sx through the formula in Winkelman 1981/Avramovic 1980 requires reference state ordering
  % Tweak the definition of the ref-vector to see that Sx is not the same and only equal to S for the same choice of reference states
  % However, it can also be observed that the slow subsystem does not depend on the choice of the reference states
  k = size(Ld,2);
  m = size(mac_con,1);
  mo000 = 1:1:m;
  T = U2T(logical([eye(k);discretize(Ld)]), mo);
  [ar_VN0, na_VN0] = T2A(T); 
  S = PST.sc_xform(ar_VN0,na_VN0,mac_con);  
  prm = ar_VN0(:,1)';  % permutation for reference ordering of the states
  %prm = [15,16,14,12,3];
  for j = 1:1:k
    prm = [prm,setdiff(ar_VN0(j,1:na_VN0(j)),prm(j))];
  end  
  [~,~,prm_bck] = intersect(mo000,prm);
  M = M(prm,prm);  % thus, the inertias are reoreded accordingly
  U = full(logical(sparse(1:m,T,1,m,k))); U = U(prm,:);
  Ma = U'*M*U;
  Sx = [[Ma\M(1:k,1:k), Ma\U(k+1:end,:)'*M(k+1:end,k+1:end)];[-U(k+1:end,:), eye(m-k)]];
  Sx = Sx(:,prm_bck);
  M = M(prm_bck,prm_bck);
  assert(norm(Sx(1:k,:)-S(1:k,:),'fro')<1e-6);  % slow part of transformation doesn't depend on choice of reference states
  S_VN0 = S*MK/S; S11_VN0 = S_VN0(1:k,1:k); S22_VN0 = S_VN0(k+1:end,k+1:end);
  Q_VN0 = Sx*MK/Sx; Q11_VN0 = Q_VN0(1:k,1:k); Q22_VN0 = Q_VN0(k+1:end,k+1:end);
  assert(norm(S11_VN0-Q11_VN0,'fro')<1e-6);
  assert(norm(sort(sqrt(abs(eig(Q22_VN0))),'ascend') - sort(sqrt(abs(eig(S22_VN0))),'ascend'))<1e-6);
  assert(norm(S22_VN0-Q22_VN0,'fro')<1e-6);
  assert(norm(Sx-S,'fro')<1e-6);
end

function [L, n_area] = T2A(T, ref)
tab = tabulate(T);
max_card = max(tab(:,2));
L = zeros(numel(tab(:,1)), max_card);
n_area = Inf(1, size(L,1));
for i = 1:1:numel(tab(:,1))
  L(i,1:tab(i,2)) = find(T==tab(i,1));
  n_area(i) = nnz(L(i,:));
end
if nargin>1
  for i = 1:1:numel(tab(:,1))
    j = find(L(i,:)==ref(i));
    assert(~isempty(j));
    L(i,j) = L(i,1);
    L(i,1) = ref(i);
  end
end
end

function T = U2T(U,mach_ord)
if nargin<2, mach_ord = 1:size(U,1); end
assert(size(U,1)==numel(mach_ord));
T = Inf(numel(mach_ord),1);
for i = 1:1:size(U,2)
  T(mach_ord(U(:,i))) = i;
end
end

function T = A2T(L)
T = Inf(nnz(L),1);
for i = 1:1:size(L,1)
  T(nonzeros(L(i,:))) = i;
end
end
